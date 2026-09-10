# Open defects of EXHALE, as the source stands on 2026-09-03

**Second revision, 2026-09-03 (late).** Section 145 landed in the tree and eight
further sections (146-154) were drafted, measured and staged outside it. Every
Status and Disposition below now says which. The short version:

* **Landed in the tree (`Update_EXHALE_stage1.md` section 145):** the flux gate reads
  the Riemann face mass flux at a re-derived threshold `2.0e-5`; the residual
  norm is the cellwise maximum and `resid_vol` is gone; every marching stop says
  it is a marching stop; the update-map diagnostic exists **and it settles
  item 3's open question**; profile files carry provenance.
* **Closed:** item 6 (base energy budget, in the tree's measurement documents),
  item 8 (landed), item 9's sign (bounded), item 12.7's diagnosis.
* **Staged for block I, measured, not applied:** sections 146 (AB ghost order),
  147 (seed interface), 148 (PLM/WENO continuation, default off), 149 (element
  assertions and tests), 150 (LW geometry bound and an FUV beam defect), 151 (H2
  photochannels and a measured cross section), 152 (characteristic base
  boundary), 153/154 (trust region and acceptance merit, default off).
* **Two new items:** **D6**, an acoustic reflection at the lower boundary that
  no static reservoir pair can remove, and the **H2 front creep**, which is now
  measured and whose base-cell explanation is excluded.
* **One finding of ours retracted by measurement (section 155):** the base cells
  are reproducible and `calc_rho` does not move `rho`. In its place is a real
  and different defect -- the accepted residual is evaluated at the *previous
  iterate's* composition and understates the written state's own by 4 to 56
  times. Acceptance today is unaffected (these roots pass on the flux gate,
  which is composition-independent); what changes is what `||R||` means
  everywhere it is quoted.
* **The front creep is now partly decided (section 10):** it is motion per unit
  *physical time*, not a step bias -- a factor 4 in `dt` puts the front at the
  same radius at matched times. That measurement turned up **a further new
  defect**: the carrier marching step does not converge to the carrier steady
  residual as `dt -> 0`, so the marching and the steady solve do not solve the
  same problem in that row. Three candidates are excluded bit-identically; the
  base inflow Dirichlet was the remaining candidate -- and was NOT the cause:
  section 158 (2026-09-04) found the two paths already one operator, the flat
  `|G+R|` a measurement made in the wrong variable, and the real defect in the
  advected variable itself (a particle mixing ratio against a non-conserved
  `n_tot`). Landed.
* **Briefing corrections recorded in place:** the front creep does not reach the
  end of the domain (1.0485 -> 1.0586 R_p at `f = 0.5` over 200k steps), and the
  `dt` scan is a CFL scan at matched physical time rather than a `dt`/`dr`
  discriminator -- no `dr` leg was run at matched time.

**Revised after review, 2026-09-03.** The first edition of this catalog was
reviewed against the working tree in `docs/open_defects_20260903_review.md`
(2026-09-03, 18:11), and two of its items were settled by measurements taken
after it was written. This revision carries those corrections. Every section now
ends with a `Disposition` line giving the review's recommended disposition
(review section 5) beside our judgment of it -- **agree**, **corrected**, or
**settled by later measurement**. What changed substantively:

* **Item 6** no longer carries the `Heating_breakdown.txt` diagnostic: it is
  repaired in the tree, and the identifier collision on `P52` is resolved below.
* **Item 4**'s reasoning is corrected. Two reservoir conditions are the right
  count for subsonic inflow; the defects are component-wise imposition, a
  ghost-center radius used in a face state, and the `massflux` branch
  prescribing all three primitive variables. One table heading was wrong.
* **Item 3**'s claim that the two routes share one operator is weakened to what
  the code supports, and the discriminator that would settle it is named.
* **Item 8 is settled**, and not the way it was posed: the dilution is real, it
  hides the base cells on the hot Uranus and fifty wind cells on WASP-121b, and
  it does **not** hide the layer. The `3.4e-3` the layer was accused of is a
  functional mismatch, not a conservation failure.
* **Item 7**'s disposition is changed: these are three concepts with three
  roles, to be named and reported separately rather than merged.
* **Item 9**'s error sign and size are reported as unknown.
* **Item 11** moves to a migration and compatibility section.
* A **Program** section carries review section 6, phases A to F, with the status
  of the work already under way.

Every "Where" below was read in the working tree on 2026-09-03 and carries the
line numbers of that tree; where a line number is quoted, the expression beside
it is copied from the file. Measured numbers are attributed to the document and
section that measured them, and none of them was re-measured for this catalog
except where the text says so.

## Summary

The items below fall into two groups, and the distinction matters more than the
severity ordering.

**Items that make a computed number wrong.** The molecular layer is solved with
a closure -- H2 in local chemical equilibrium, no transport term -- that is
violated on its own solution by a factor 44 to 117 across the H2 front, and the
operator that would repair it exists but has never produced a steady state, so
it is default off. Everything the molecular layer reports below about 1.1 R_p
(the H2 front radius, the H3+ column, the layer temperature) rests on that
closure. Separately, the lower boundary is imposed component by component rather
than characteristically, a ghost-cell center radius is used as a face state
(a first-order-in-`dr` error, and the one part of the boundary account the review
confirms), the `massflux` branch prescribes all three primitive variables where
subsonic inflow admits two, and the base layer above it is not steady: a state
the steady solver accepts leaves it as soon as it is marched. A slab-geometry
fluorescent-trapping factor is applied to a spherical, outward-open layer, with
an error of unknown sign and size.

**Items that are questions of definition or of judgment, not of arithmetic.**
"Converged" has two definitions in this code (a zero of the steady residual, and
a fixed point of the marching loop) and they are demonstrably different states --
though, as the review points out and the first edition did not, the two routes
are not even the same discrete map, so the discrepancy has a second possible
source. The marching stop, the solver handoff and the acceptance gate are three
concepts sharing two windows and no names, and a run stopped by the first prints
"converged". These do not make any individual number wrong; they decide which
states the code is willing to call answers.

**One item that was a question and is now settled.** The residual norm's
weighting was under investigation in the first edition. It is answered: the
weighting does dilute, but it hides the base cells on the hot Uranus and fifty
wind cells on WASP-121b, **not** the layer -- and the `3.4e-3` the layer stood
accused of is the difference between a cell-centered product and the face flux
the scheme conserves, which above 1.10 R_p is one unique double-precision value
over three hundred cells. The acceptance norm becomes the maximum over cells of
each cell's own scaled residual, and the flux gate moves to the faces.

## Severity table

| # | Item | Calculation affected | Status | Review disposition (sec. 5) and our judgment |
|---|---|---|---|---|
| 1 | H2 solved in local equilibrium with no transport term (default) | every molecular-layer quantity: H2 front, H3+ column, layer T, and the P23 comparison against Koskinen and Frelikh | **open.** D1 (decouple from oxygen) not started. | Confirmed. Decouple from oxygen, add a closure-validity diagnostic now, make transport mandatory once the coupled solve is qualified. **Agree.** |
| 2 | The coupled carrier/wind steady solve does not converge | blocks the repair of item 1 | **open, and the obstruction has moved off the solver.** Trust region built (153/154), model exact, still stops at `\|\|R\|\| = 1.004` -- the BASE CELL's row. Default off, staged. | Confirmed. Add a scaled trust region and verify the model before changing defaults; do not loosen the carrier gate. **Agree.** |
| 3 | The base layer does not hold still: the steady root is not a fixed point of the marching loop (the excursion is a damped oscillation, period ~2100 steps) | anything read below ~1.1 R_p on a "converged" state; the meaning of convergence | **open on the hot Uranus.** The update map (145.6, landed) excludes an operator mismatch: `G_dt -> -R` at first order. Four candidates excluded, none is the cause. Section 152 makes the departure *worse*. | Reported measurement, explanation incomplete. Measure the complete update-map defect before attributing the motion to an acoustic mode. **Corrected** -- our "one shared operator" claim was too strong. |
| 4 | The lower boundary is imposed component by component, a ghost-center radius is used in a face state, and the `massflux` branch prescribes all three primitive variables | base-face state, base T, any diagnostic reading v(1) | **prescription found and measured, staged (152).** Sawtooth gone (`v(1)` `-251.9 -> +14.4` cm/s), reproducibility repaired; marching departure worsens 2.09x; new item **D6**. | Location mismatch confirmed, interpretation partly rejected: two reservoir conditions are appropriate for subsonic inflow. Replace component-wise face overwrites with a characteristic face solve. **Corrected.** |
| 5 | Layer quantities do not converge under grid refinement | the H2 front, the rate, the layer T minimum on the gate rung | open; not re-run. Section 152's three-grid test is the next measurement and was **not run**. | Reported measurement; cannot yet be assigned uniquely to item 3. Repeat only after items 3 and 4 are repaired. **Agree.** |
| 6 | The base layer's energy budget did not close on the accepted state | any energy statement below 1.08 r_base | **CLOSED.** `gamma_eff` in the advection term plus a face-flux-finished state: 0.99 / 1.01 / 1.02. Provisional against the landed norm. | Physical requirement confirmed, cause open; the cited heating-output defect is already fixed. Add local and shell-integrated energy gates. **Corrected** -- the diagnostic defect is removed from this item. |
| 7 | The marching stop, the solver handoff and the acceptance gate are three concepts measured through two windows | where an unconverged run is allowed to stop, and what may be called a solution | **partly landed (145.5).** Every marching stop now says it is one and prints all three gates; unification formally withdrawn. Three distinct *names* still not given. | Different definitions confirmed but not necessarily one defect: they serve termination, handoff and acceptance roles. Do not force one threshold onto all three; rename and document by role. **Corrected** -- our unification proposal is withdrawn. |
| 8 | What the residual norm hides; where the layer's `3.4e-3` comes from; and at which composition the accepted residual is read | which cells acceptance can average away; the meaning of the flux gate's number | **LANDED (145)** for the norm and the gate. **One finding RETRACTED (155):** the base cells reproduce (`drho/rho <= 1.7e-15`; `calc_rho` moves nothing). **A real defect in its place, staged:** the accepted `\|\|R\|\|` is read at the PREVIOUS iterate's composition and understates the written state's own by **4 to 56x** (Phase C cuts WASP 56 -> 5.7). | Open with a verified formulation concern: neither norm branch is a cellwise maximum of relative residuals; add dual gates and contribution diagnostics. **Settled by later measurement** -- and the review's formulation concern is exactly right. |
| 9 | Fluorescent trapping carries a plane-parallel geometry | H2 photodissociation rate at depth | **sign and size now BOUNDED (150).** The slab over-traps: 0.2% at the base, 7-10% through the H2 layer, 20% at top, 6.5% rate-weighted. Table deliberately **not** rebuilt. | Model mismatch confirmed, sign and size unverified. Add a geometry provenance flag; do a two-geometry comparison before replacing the table. **Corrected** -- our claimed direction is withdrawn. |
| 10 | The steady residual fills the base ghost before its ionization sweep | the residual during relaxation, not at the fixed point | **repaired, staged (146 + 147).** Ordering fixed by 146; `F(Y)` reproducibility needed 147's seed interface (1406/1500 -> 0/1500, bitwise). 146 needs a golden refresh; 147 is byte-identical. | Ordering confirmed, effect unmeasured. Make residual evaluation a deterministic composition fixed point; test repeated `F(Y)`. **Agree**, and both halves are now done. |
| 11 | Restart state transport | fixed 2026-09-03; legacy files and two diagnostic columns remain | fixed; `restart_schema=3` in the provenance header, and the round-trip fixture is rebuilt (149 / item 12.4), staged. | Main repairs confirmed in code. Add a restart schema version and an executable round-trip test. **Agree** -- reclassified as compatibility, not an open numerical defect. |
| 12 | Smaller open items (H2 photoionization branches, stale regression cases) | as stated in each | 12.1-12.3 **superseded by 151**. 12.7 reframed and addressed by 148. **12.8 NOT REPRODUCED** and still open. | Per item in section 12. **Agree**, with 12.8 raised: elemental double-counting blocks any coupled result. |
| D6 | The lower boundary reflects sound completely, and no static reservoir pair can stop it | the base-layer ringing of item (AA); whether a steady root and a marching run can share a boundary | **new, open (152).** `\|R\| = 0.956` against the legacy `0.953`, amplitude- and grid-independent | Raised by the Phase C implementation, not by the review. |
| FC | The H2 front creeps outward and never stops | every transported-H2 result; whether a transported steady state exists at all | **new, open** (`p23_transport_on_state.md` sections 8-10). Motion in PHYSICAL time (CFL 0.6/0.3/0.15 agree at matched times); base cells excluded; and a **separate new defect** -- the carrier marching step does not converge to `-R_carrier` as `dt -> 0` | Not a review item. Step-bias reading now excluded; physics vs numerics still **not** decided. |

---

## 1. H2 is solved in local chemical equilibrium with no transport term, and that closure is violated on its own solution

### Symptom

The molecular layer's H2 fraction falls four decades over a narrow radial
interval and is essentially gone above it, which no published profile of a
comparable planet shows. On the 9 microbar hot Uranus, `f(H2)` goes from 0.999
at the base to 1.42e-2 at 1.15 r_base and 1.5e-6 at 3.0 r_base; Koskinen et al.
(2022) Figure 8 reads 0.963 and 0.403 at the same two radii, and Frelikh's
super-Earth Figure 21 reads 0.269 and 0.013 (`docs/p23_published_profiles.md`
section 6). n(H3+) follows H2 down: our maximum is 3.09e5 cm^-3 at 1.023 r_base
falling to 1.7e-5 at 2.00, against a Koskinen maximum of 1.2e4 cm^-3 that is
flat over 1.23-1.54 and still 7.2e3 at 2.00.

### Where it arises

* `src/modules/radiation/ionization_equilibrium.f90:1145-1177`. The H2 partition
  is a root of the cell's own balance rows unless something overrides it. The
  default is stated at line 1145, `ieq_cell%x_h2_fixed = .false.`, and the
  override at lines 1163-1166 is reached only when the carriers are transported:

  ```fortran
  if (thereis_mol .and. carrier_transport                   &
      .and. (bg_ready .or. do_load_IC)) then
     ieq_cell%x_h2_fixed = .true.
     ieq_cell%x_h2_fix  = 2.0d0*nmol_eq(j,1)/nh(j)
  ```

* `src/modules/nonlinear_system_solver/System_HeH_mol.f90:553` and
  `System_HeH_mol_metals.f90:278` are the two places the flag acts. Where it is
  set the H2 row becomes an identity, `fvec(4) = x(4) - ieq_cell%x_h2_fix`;
  where it is not, that row is the local chemical balance and nothing in it
  contains `d(r^2 v n_H2)/dr`.

* `src/modules/files_IO/input_read.f90:1662` decides the default:

  ```fortran
  if (.not. carrier_transport_stated) carrier_transport = thereis_oxychem
  ```

  and `src/modules/init/parameters.f90:292` declares
  `logical :: carrier_transport = .false.`. So a molecular run **without** the
  oxygen cycle solves H2 in local equilibrium unless `Molecular carrier
  transport: True` is written by hand; a run with the oxygen cycle transports
  it. That is one physics closure selected by an unrelated key, and it is
  deliberate -- the comment at `input_read.f90:1655-1661` says so and points at
  `docs/supersonic_molecular_base.md` section 13.6.

* The transport operator itself is
  `src/modules/lower_atmosphere/diffusive_photochemistry.f90`
  (`photochemical_transport_step` at line 534, `carrier_steady_residual` at
  line 777), and the coupled Newton branch that makes n(H2) a fourth unknown is
  `Coupled carrier solve` (`input_read.f90:516-528`;
  `parameters.f90:298` declares `carrier_in_newton = .false.`).

### Why it is wrong

The continuity equation of a transported species is
`r^-2 d(r^2 v n_i)/dr = P_i - L_i`. Dropping the divergence is legitimate only
where it is small against the chemical terms, i.e. where the chemical time is
short against the flow time. In the region where `f(H2)` falls its four decades
that condition fails by two orders of magnitude, in the direction that keeps
more H2 than the closure allows.

### Measured size

From `docs/supersonic_molecular_base.md` section 13.6 (lines 1867-1890, the
cycle-resolved budget), the omitted divergence against the net chemical rate:

| r [R_p] | 1.05 | 1.15 | 1.25 |
|---|---|---|---|
| `f(H2)` | 0.966 | 1.43e-2 | 1.47e-4 |
| net H2 loss [cm^-3 s^-1] | 6.75e3 | 1.33e2 | 8.51e1 |
| `\|r^-2 d(r^2 v n_H2)/dr\|` | 2.96e5 | 1.55e4 | 4.74e1 |
| ratio | **44** | **117** | 0.6 |

`supersonic_molecular_base.md:1928-1934` states the direction and the limit of
the measurement: "The equilibrium profile is therefore not a steady state of the
full continuity equation anywhere below ~1.25 R_p ... The measurement is
one-sided and should be read that way: it shows the closure is violated on its
own solution, not how much H2 a transported solution would keep."

The one comparison that is like for like is the H2 replacement time, in
`docs/p23_published_profiles.md` section 7.6 (lines 677-683):

| | 1.05 r_base | 1.15 r_base | 1.25 r_base |
|---|---|---|---|
| Koskinen, n(H2)/\|T\| | 5.0e8 s | 1.2e7 s | 9.7e5 s |
| EXHALE, n(H2)/\|advection\| | 2.0e7 s | 1.1e5 s | 1.5e4 s |
| EXHALE, n(H2)/net chemical | 8.9e8 s | 1.3e7 s | 8.2e3 s |

The chemistry that consumes H2 is helium ions, not thermal dissociation:
across the front at 1.15 R_p, R20 (`He+ + H2 -> HeH+ + H`) and R17
(`He+ + H2 -> He + H + H+`) carry 95 percent of the net between them, and at
1.25 R_p R17 alone carries 99 percent
(`supersonic_molecular_base.md:1908-1913`). At the first cell of the He/H = 30
rung the helium ions destroy 99.9 percent of the H2 and R17 alone 93.5 percent,
while the residence time in that cell is 19 times shorter than the H2 chemical
lifetime (`docs/p23_thermal_budget.md:445-476`). The published thermal
dissociation barrier is not what truncates our H2: our front sits at 765 K
against the 48,350 K barrier of R12 that both codes take from Baulch et al.
(1992) (`supersonic_molecular_base.md:1789-1792`).

Turning the transport on moves the front outward and does not stop it. Over
20,000 marching steps (`supersonic_molecular_base.md:2062-2071`), `f(H2)=1e-4`
moves from 1.2698 to 1.3072 R_p, and after 32,000 further steps it is at
1.3375 R_p and still moving. With the wind JFNK-converged and the carriers
relaxed once, `f(H2)` at 1.20 R_p is 212 times its untransported value,
n(H3+) there goes from 0.16 to 157 cm^-3, and the cold dip at the front is gone
(T 768 -> 1677 K). A real feedback drives it: where H2 arrives, the
helium-ionizing continuum is absorbed and n(He II) falls by 50-100 at
1.15-1.20 R_p (`Update_EXHALE_stage1.md` section 128.6, lines 18087-18095).

### Remedy options

1. **Physics first: make the carrier transport the default and finish it.**
   This is the physically correct closure -- the term is not negligible, it is
   the larger one -- and the operator exists, is second-order accurate
   (measured orders 1.938 to 1.985 over five grids, `Update_EXHALE_stage1.md`
   section 134, lines 19452-19460) and is grid-converged in the front position
   to 0.2 percent at 150,000 steps. What blocks it is item 2 below: no steady
   state has been produced with it on. Cost: none in code; the whole cost is in
   item 2. Verification: `info = 0` from a steady solve started at a marched
   state whose front has stopped, and a golden refresh of the six `mol_*`
   cases.
2. **State the closure's validity at the code site and gate on it.** A cheap
   partial measure that does not wait on item 2: compute the ratio measured
   above at run time and refuse (or warn on) a molecular run in which the
   omitted term exceeds the chemical one anywhere the H2 fraction matters.
   Cost: small. Side effect: it does not make any number right, it stops the
   code from quoting a wrong one silently.
3. **Do not couple the closure to the oxygen key.** `carrier_transport =
   thereis_oxychem` makes the H2 closure depend on whether the oxygen cycle is
   on, which is a different physical question. Whatever the default becomes, it
   should be one rule. Cost: one line plus a golden refresh of the oxygen
   examples.

### Status

Open. The operator is in the tree and is default off for a molecular-only run;
`Update_EXHALE_stage1.md` section 134 lines 19537-19540 states the condition for
flipping it ("a marching run stopped the front and a steady solve from that
state returned `info = 0`; the first half is nearly done and the second is
blocked"). Every molecular-layer number in `docs/supersonic_molecular_base.md`
section 14, in `docs/p23_*`, and in the six `mol_*` regression goldens was
produced with the closure on.

**Disposition.** Review section 5: "Critical physics defect for molecular
results. Decouple the option from oxygen, add a closure-validity diagnostic
immediately, and make transport mandatory after the coupled solve is qualified."
**Our judgment: agree**, and this is the review's executive judgment as well --
"the inability of the coupled steady solver to converge does not make the
local-equilibrium closure acceptable." Note the review's ordering, which differs
from a naive reading of this catalog: the closure-validity diagnostic and the
oxygen decoupling are immediate, while making transport the default waits on
item 2 and on Phases A to C.

---

## 2. The coupled carrier and wind steady solve does not converge

### Symptom

With `Molecular carrier transport: True` and `Coupled carrier solve: True`, the
JFNK solve stops at `info = 2, STAGNATED`, holding a carrier row two orders
above its gate. The alternating (Picard) route it replaced does not converge
either, and moves the H2 front the wrong way.

### Where it arises

* `src/modules/time_step/steady_newton.f90`, `solve_steady_jfnk` (line 1445) and
  `solve_steady_ptc` (line 1135). The carrier unknown is switched in by
  `set_carrier_unknown` (line 240); `nvar_jac` goes 3 -> 4 and the band
  half-width 8 -> 11, with 23 Jacobian colors against 17
  (`Update_EXHALE_stage1.md:20515-20516`).
* The pseudo-transient term. `Update_EXHALE_stage1.md` section 142 identifies it as
  what turns the step into an ascent: at `dtau0 = 1` the model's predicted slope
  along the step turns positive at outer iteration 5 and reaches +0.80 by 12,
  while the same Jacobian, preconditioner and frozen model at `dtau0 = 1e6`
  give -36 to -14 there (lines 21545-21551).
* The step-length controls: `levenberg_marquardt_descent` (line 884), the
  Grippo non-monotone line search, and `n_stall_best = 20`
  (`steady_newton.f90:1513`) which is what makes the solve stop and say
  STAGNATED instead of burning 500 iterations.
* The carrier gate is separate from `||R||` by design:
  `steady_residual.f90:577-633`, `steady_gates_met`, third gate
  `carrier_relnorm < carrier_resid_th`, with
  `carrier_resid_th = 1.0d-3` at `parameters.f90:310`. The reason is written at
  `steady_newton.f90:1019-1027`: the hydrodynamic rows are divided by a bound on
  each row's largest term and converged states sit at 1e-6 there, while the
  carrier row is divided by the terms themselves and a converged row sits near
  1e-3, so one number cannot hold both.

### Why it is wrong

Numerics, not physics. A Newton step is a descent direction of the merit only
while the model that produced it is valid over the length of the step. The
measurements say the model is valid over a much shorter arc than the step the
solver takes, and that the pseudo-time term is what pushes the predicted slope
positive. A rule that raises `dtau` when the model says the step ascends was
built and measured and does not fix it, because after the raise the line search
still cuts to `lam = 9.5e-7`.

### Measured size

Coupled solve, entry from an accepted wind (`Update_EXHALE_stage1.md:20602-20610`):

| | entry | exit |
|---|---|---|
| `\|\|R\|\|` with the carrier row folded in | 1.221e-2 | 4.073e-3 |
| carrier row, worst cell / its own terms | 6.04e-1 | **1.36e-1** |
| carrier row, volume weighted | 7.72e-3 | 4.07e-3 |
| `info` | -- | **2, STAGNATED at 49 iterations** |

Against a gate of 1e-3, the exit is 136 times out on the worst cell. Three
independent Krylov settings all stop at a carrier row of 1.02e-2 to 1.03e-2
(`Update_EXHALE_stage1.md:21528-21538`). On the 1 microbar gate configuration
(`mol_carrier`) the solve stagnates at `||R|| = 2.30e-3` with a **flux spread of
0.244 against a gate of 5e-3** (`Update_EXHALE_stage1.md:20694-20703`). The
alternating route it replaced walks the front outward one cell per pass
indefinitely (0.082 R_p per pass of 0.1 flow times, about 2050 cm/s against a
local gas speed of 3000 cm/s) and gets the sign of the front's motion wrong
relative to the coupled solve, which moves `f(H2)=0.5` inward from 1.1421 to
1.0618 R_p (`Update_EXHALE_stage1.md:20384-20391`, 20631-20636).

Cost of the coupled route: 23 Jacobian colors and 778 residual evaluations at
14.51 s against 17 and 64 at 1.00 s for the three-unknown solve; 11 residual
evaluations per outer iteration against 39 (`Update_EXHALE_stage1.md:20655-20666`).

### Remedy options

1. **A trust region sized by the arc over which the model is valid.** This is
   what `TO_BE_DONE.md` (P51 continuation) names and it is the one candidate
   the measurements point at rather than away from. The ray scan of
   `docs/p51_coupled_stagnation.md` section 3 measures that arc directly, and
   `pgmres` can return `A_z x` from its own Arnoldi relation, so the
   directional derivative costs no extra residual evaluation. Cost: moderate,
   confined to `steady_newton.f90`. Side effect: it changes the step of every
   JFNK solve unless gated on the carrier branch, so it needs the byte-identity
   check on the atomic matrix. Verification: `info = 0` with the carrier row
   below 1e-3 on the hot Uranus, and 10/10 byte-identical `make check`.
2. **Raise `dtau` on a predicted ascent (built, measured, rejected).** It does
   exactly what it is specified to do and converges neither solve: it improves
   the section 139 hand-off by 17 percent and turns the P23 marching stop's
   500-iteration crawl at `||R|| = 1.582e-5` into a 17-iteration abort at
   3.502e-5. Kept as `p51b_dtau_rule.patch` in the P51 scratch directory, not
   applied.
3. **Stop coupling and integrate the carriers in time to a stopped front,**
   then solve the wind alone at that front. This is not a steady state of the
   coupled system and should not be quoted as one, but it would give a defensible
   H2 profile where the current default gives one that is known to be wrong. Cost:
   compute time (150,000 steps and the front is still drifting 0.0034 R_p per
   20,000 steps). Side effect: the answer would depend on when the run stopped.

### Status

**Updated 2026-09-03, late: the trust region was built, it works as a trust
region, and the coupled solve still does not converge -- because the obstruction
is not in the solver.** Sections 153 and 154, staged, default off.

**The model is not the problem, and that is the result.** Over 47 outer
iterations the reduction ratio has a median of 0.933, 42 of 47 lie in [0.5, 2],
and on the first four iterations it is 1.000, 1.000, 0.999, 0.998 -- the model is
exact to the printed digits. **The ray test refused zero models**, in either run.
Element-budget refusals fall from 18-19 under the line search to **0** under the
trust region.

**And the gates are still missed.** On the current tree binary, from the same two
hand-off states, the trust region and the line search both return
`||R|| = 1.004` and `1.005`. The D2 success list is not met: `||R||` by four
decades, the flux spread by a factor 8 to 24, the carrier row by a factor 9.

**Where the residual lives.** The worst cell is the **base cell** -- `j = 1`,
`r = 1.0002` -- in the energy row in 20 of 28 reported iterations and the mass
row in 6; below the escape radius the split reads mass `4.13e-1` and energy
`3.14e-1` against momentum `7.24e-3`. Under section 153's measure the census was
the base cell's momentum row in 29 of 46 iterations and the carrier row at the
front in 15. `p51_s154_acceptance_merit.md` puts it plainly: "what is left is a
single cell at the lower boundary whose row no step can balance because the state
there is imposed rather than solved. **That is the base-ghost item, not a Newton
problem.**"

A cellwise maximum of `1.004` means one cell's row is 100 percent out on its own
terms. **So item 2 is now downstream of item 4**, and the natural next
measurement is the coupled solve on the section 152 boundary, which has not been
run.

**And the state the coupled solve starts from is produced by an operator that
does not agree with the residual it is then asked to zero.** Item FC(2) measures
the marching's carrier step against `-R_carrier` over four decades of `dt` and
finds `|G + R|` `dt`-independent to four significant figures -- so in the carrier
row the marching and the steady solve are not solving the same problem, and every
coupled hand-off in this item begins at a state the marching produced under the
other one. The inconsistency is worst at the base cells (1.9e-2 of the row's own
terms), which is where this item's residual floor also sits. Whether the two are
the same defect is not established; the repair in progress is section 157.

### Status

Open. `Molecular carrier transport` and `Coupled carrier solve` are both default
off and `Update_EXHALE_stage1.md:21503-21504` states that nothing measured with them on
is quoted as a rate. (P51) is closed for the diagnosis; (P51 continuation) is now
answered -- the trust region exists (`Coupled carrier solve: True trust_region`,
`carrier_newton_trust_region` defaulting to false, an unrecognized third word
refused rather than ignored) and it is measured and **not adopted**.

**Disposition.** Review section 5: "Confirmed implementation and open convergence
failure ... Add a scaled trust region and verify the Jacobian/model before
changing defaults. Do not loosen the carrier gate to obtain a nominal success."
**Our judgment: agree, and done -- with the opposite outcome to the one the
recommendation anticipated.** The trust region was added, the model was verified
and is exact, the carrier gate was not loosened, and the solve still fails. The
review's own instinct -- item 2 "blocks the physically correct resolution of
item 1" -- survives; what changes is that item 2 is itself blocked by item 4.

**Disposition.** Review section 5: "Confirmed implementation and open
convergence failure; historical outcome not rerun. Add a scaled trust region and
verify the Jacobian/model before changing defaults. Do not loosen the carrier
gate to obtain a nominal success." **Our judgment: agree.** The review's Phase D2
specifies the trust region in more detail than this item did -- dogleg or
truncated Krylov step, predicted reduction from the same scaled model GMRES
uses, radius updated by the actual/predicted ratio, carrier headroom and
thermodynamic admissibility enforced inside the step, and a model rejected when
its directional derivative disagrees with a finite-difference ray test -- and
adds a success criterion that is not `info = 0`: all gates, elemental closure,
fixed point of the complete production update, and a front that no longer drifts
under further transport relaxation.

---

## 3. The base layer does not hold still: a state the steady solver accepts leaves it as soon as it is marched

### Symptom

Reload an `info = 0` state and march it with nothing able to stop the run, and
the mass flux in the base layer moves immediately. The flux at `r = 1.03`
executes an excursion of order the wind flux itself over a few thousand steps
and settles onto a different, still non-flat profile that both acceptance gates
would now reject.

### Where it arises

**Corrected after review.** The first edition said the spatial operator is
shared by construction, so this cannot be two different discretizations. The
first half is true and the conclusion does not follow. `assemble_residual`
(`src/modules/time_step/steady_residual.f90:69`) does call the same `Reconstruct`
and `RK_rhs` the marching loop calls -- but that is the Euler part only. One
production step also runs an operator-split carrier advance, an ionization
sweep, a semi-implicit energy update, repeated boundary fills, an optional
viscous and conductive stage, and an optional Shapiro filter. The steady
residual eliminates some of that physics locally and represents other parts as
source terms. **Sharing `Reconstruct` and `RK_rhs` does not by itself make the
two maps share a fixed point.**

* `src/EXHALE_main.f90:697-1406`, the marching loop (`do while` at 697, `enddo` at
1406); within one
step: the ionization sweep at 910, the
  semi-implicit energy step, `Apply_BC(u)` at line 940, the optional viscous
  and conduction stage, and the Shapiro filter at lines 960-965.
* `src/modules/states/Apply_BC.f90:259-276`, `shapiro_filter`, the 1-2-1 low-pass
  pass. It is opt-in (`shapiro_eps <= 0` returns immediately) and on this
  configuration it makes the departure worse, not better.
* `src/modules/time_step/RK_rhs.f90:29-114`, which is where the momentum row's
  near-cancellation lives: in a quasi-hydrostatic layer the flux divergence and
  the gravitational source are each of order 1.1e5 times the wind's mass flux
  and their difference is what drives the layer.

### Why it is wrong

The two statements the code makes about a solution -- "the steady residual is
zero" and "the marching loop does not move it" -- are not the same statement in
the base layer, and every published number rests on the first.

**What the measurements support, and what they do not.** They support that the
departure is a residual being integrated rather than a mode growing: it starts
at step 1 at a finite size, grows linearly in step number, is predicted by the
accepted state's own momentum-row residual to nine percent, and is insensitive
to the time step. That is evidence against a discrete instability and against a
physical instability. **It is not evidence that the two maps are the same
operator**, which is what the first edition inferred, and item 8's later
measurement removes the third candidate (the row scale and the norm) without
supplying a cause.

**The discriminator that would settle it**, from review section 3.3, is the
defect of the complete update map,

```
G_dt(q) = [Phi_dt(q) - q] / dt,
```

evaluated at the accepted steady state for decreasing `dt` and compared with
`assemble_residual` row by row and cell by cell. Agreement would establish that
the two routes solve the same discrete equations and would move the question
back onto the state, i.e. onto how accurately the root is solved. Disagreement
would locate a splitting, boundary-ordering, or state-refresh defect, and the
comparison should be broken down after each stage -- Euler, chemistry, energy,
carrier, diffusion, boundary, filter. A run with the Shapiro filter enabled
cannot be expected to share a fixed point with an unfiltered residual at all,
so either the filter's steady operator belongs in the residual or the filter
must be prohibited in steady-root validation.

### Measured size

`docs/p55_base_mode.md` sections 3-6 and `Update_EXHALE_stage1.md` section 144.

The departure is a residual being integrated, not a mode growing. The flux at
`r = 1.03` moves linearly in step number for the first twenty steps at
4.71e-04 of `F_0` per step, and the state's own momentum-row residual there
predicts 4.31e-04 -- nine percent agreement with no free parameter
(`p55_base_mode.md:213-223`). The accepted state was flat to 1.043e-04 over
`r >= 1.01`, so **one step undoes 4.5 times the flatness the solve achieved**.

Judged by the largest term each row contains, the accepted hot-Uranus state's
energy row is out by 1.2e-3 at `r = 1.03` and 6.0e-3 at 1.05, its momentum row
by 2.5e-5, and cell 1 by 1.4e-1 and 9.9e-1 (`p55_base_mode.md:11-16`).

The excursion itself, from `docs/p54_base_layer_mass_flux.md:239-252`
(face mass flux in units of `F_0`, over 20,000 steps):

| r [R_p] | min | max | peak to peak | period |
|---|---|---|---|---|
| 1.005 | 0.221 | 2.251 | 2.030 | 2050 steps (3446 s) |
| 1.03 | 0.702 | 1.667 | 0.964 | 2100 steps (3530 s) |
| 1.10 | 0.976 | 1.230 | 0.254 | 2150 steps (3614 s) |
| 1.20 | 0.998 | 1.060 | 0.062 | 2300 steps (3866 s) |

The oscillation as measured is **damped**, not sustained: the autocorrelation at
the first repeat is 0.37 and the layer settles by about step 12,000 onto a
different non-flat profile with the sign of the deviation reversed
(`p54_base_layer_mass_flux.md:262-276`). At 20,000 steps that state reads
`du = 6.22e-03`, `||R|| = 3.08e-05` volume weighted, and a gate flux spread of
2.69e-02 -- so both gates would reject it, by a factor 11 on the flux gate and
2 to 3 on the residual gate. The period is *consistent with* a base-to-1.10 R_p
acoustic round trip (measured 3288 s against 3446-3614 s) but the document is
explicit that this is one measured number against another and not an
identification of a mode (`p54:254-260`, `p55:422-425`).

What does and does not move it (`p55_base_mode.md:246-296`, 3000 steps from the
same accepted state, peak-to-peak at `r = 1.03` in units of `F_0`, compared at
equal model time):

| variant | hot Uranus | vs baseline | WASP-121b | vs baseline |
|---|---|---|---|---|
| A baseline, WENO3 | 0.534 | 1.00 | 2.148 | 1.00 |
| B Shapiro filter 0.5, every 4 steps | 22.55 | **42.3** | 1.976 | 0.92 |
| R restart from the marching's own 3000-step state | 0.013 | 0.025 | 0.675 | 0.31 |

Quartering the time step changes the excursion by 3 percent at equal model time.
First-order reconstruction at the base leaves WASP-121b unchanged and makes the
hot Uranus 12 times worse. Freezing the first 20 cells does nothing on
WASP-121b (`Update_EXHALE_stage1.md:22009-22014`).

**WASP-121b turned out to be a different problem and is closed.** Its departure
was the restart running physics its root was not converged under (section 11
below); with the coupling header the peak-to-peak at `r = 1.03` falls from 4.296
to 0.0174 and `du` ends at 5.187e-05 against the 5.03e-05 the accepted state
carries (`Update_EXHALE_stage1.md:22062-22072`). The hot Uranus never had that
discontinuity -- its input sets `Secondary_ionization: Immediate` -- and is
unchanged to five digits.

**Re-solving the hot Uranus on the current tree does not fix it either.** The
root found with section 143's row scales returns `info = 0`, `||R|| = 9.958e-06`,
gate flux spread 2.501e-04, and its own flux is still 3.433e-03 out over
`r >= 1.03`, fourteen times the gate's number. Marched 3000 steps
(`Update_EXHALE_stage1.md:22150-22162`):

| root | p2p at 1.005 | 1.03 | 1.10 | 1.20 | du end |
|---|---|---|---|---|---|
| pre-143 (mass row only) | 1.5476 | 0.5336 | 0.0710 | 0.0125 | 6.387e-04 |
| pre-143 + the restart handling | 1.5473 | 0.5296 | 0.0705 | 0.0126 | 6.389e-04 |
| section 143 root, tree binary | 1.5481 | **0.5333** | 0.0719 | 0.0132 | 3.363e-04 |

Unchanged to 0.06 percent. So on the hot Uranus the cause is neither the restart
nor the row measure.

### Remedy options

1. **A characteristic (non-reflecting) lower boundary.** `TO_BE_DONE.md` (AD)
   names this as the first untried candidate. A subsonic inflow carries one
   outgoing characteristic, so exactly one piece of interior information belongs
   in the ghost; the present closure states two thermodynamic variables (item 4).
   A boundary that does not reflect the outgoing acoustic characteristic would
   remove the trapping the measured period is consistent with. Cost: a rewrite
   of `BC_component_constrho` and a golden refresh of every case. Side effect:
   changes the base state of every planet. Verification: the same reload-and-march
   experiment, peak-to-peak at `r = 1.03` below the flatness the solve achieves;
   plus the full byte-identity matrix showing which goldens moved and why.
2. **Base-cell refinement.** `TO_BE_DONE.md` (AA) records this repairing the
   He/H = 0.3 rung, and it demonstrably halves the cell-1 artifact of item 4.
   It does **not** remove the oscillation: the swing grows by 21 percent at
   He/H = 0.3 and falls by 33 percent at He/H = 30 under a factor-2 refinement
   (`supersonic_molecular_base.md:2724-2738`). Cost: +19 percent per step and
   2.40x wall clock to convergence. Verification: as above.
3. **Not the Shapiro filter.** As configured it is measured 42 times worse than
   doing nothing on the molecular hot Uranus. CETIMB applies it to a base that is
   genuinely ringing; applied to a state that is drifting because its residual is
   not zero it is a low-pass filter on the wrong signal.
4. **Not a smaller time step.** Quartering it changes the excursion by
   3 percent at equal model time.

### Status

**Updated 2026-09-03, late: the update-map question is ANSWERED, and it answers
in favor of the first edition's premise rather than the review's doubt.** The
tool the review asked for is built and landed -- `EXHALE_UPDATE_MAP=<f1,f2,...>`
(`Update_EXHALE_stage1.md` section 145.6, `docs/p55_update_map.md`), with the marching
loop's own body as `Phi_dt` rather than a copy of it. Max over cells of
`|G_dt + R|`:

| `dt/dt_CFL` | hot U mass | momentum | energy | WASP mass | momentum |
|---|---|---|---|---|---|
| 1 | 6.52e-01 | 5.33e-01 | 9.60e-01 | 9.63e-02 | 1.33e-01 |
| 0.1 | 8.93e-02 | 5.27e-02 | 1.40e-01 | 1.09e-02 | 1.69e-02 |
| 0.01 | 9.23e-03 | 5.17e-03 | 2.34e-02 | 1.10e-03 | 1.73e-03 |

**`G_dt` approaches `-R` at first order in `dt` on both planets, so the
production map and the steady residual do share a fixed point**, and the
42 percent discrepancy section 144.1 measured at a CFL step is the splitting
truncation error of an explicit step, not an operator mismatch. The stage
breakdown adds that `Apply_BC`'s conservative-primitive round trip moves the
interior by a relative `4.2e-16` -- so item 10 is not a source of motion here --
and that the semi-implicit energy step is a rate to five digits. The one term
the residual has no counterpart for is the composition pressure projection,
which converges on the hot Uranus and reads `1.124` against a residual of
`1.133` on WASP-121b, where it is not small.

**One exception, and it is provable rather than measured-and-hoped.** With
`Shapiro filter: 0.5 1` the two cannot share a fixed point at all: the filter's
increment is applied per step and not per unit time, so as a rate it grows as
`1/dt` -- `1.187e+01`, `1.189e+02`, `1.188e+03` at the three time steps -- and
`|G_dt + R|` follows it instead of converging. Either the filter's steady
operator goes into the residual or the filter is prohibited in steady-root
validation; not decided.

### Status

Open on the hot Uranus, `TO_BE_DONE.md` (AD). Closed on WASP-121b. **Four
candidates are now excluded by measurement and none is the cause**: the restart
(section 144.2, 144.3), the row scale (section 144.5), the acceptance norm
(item 8; it moves the 3000-step departure by 0.2 percent, 0.5333 against
0.5334), and now an operator mismatch between `assemble_residual` and `RK_rhs`
(section 145.6). **What remains is the base boundary itself** -- review Phase C,
and items (AA), (AB) and (W) -- which is where the work has moved.

**One caveat on this item's own numbers, from section 155.** Item (AD)'s
measurements quote `||R||` of the accepted roots, and that number is evaluated at
the previous iterate's composition: on the five states measured it understates
the written state's own residual by 4 to 56 times, and every one of them is above
`1e-5` at its own composition. This does not change which states were accepted --
they passed on the composition-independent flux gate -- and it does not change
the departure measurements, which are marching excursions and not residuals. It
does mean that "the accepted root's residual was `1e-5`" is not a statement this
item can lean on.

**Disposition.** Review section 5: "Reported measurement; explanation
incomplete. The Euler spatial routines are shared, but the full operators are
not identical. Measure the complete update-map defect. Repair ordering or
splitting discrepancies before attributing the motion to a physical acoustic
mode." **Our judgment: corrected, then settled by the review's own test.** The
first edition's shared-operator inference was withdrawn as unsupported; the
review's B2 diagnostic was then built and it supports the conclusion the
inference had reached, on evidence the inference did not have. The correction
stands as a correction of the reasoning, not of the answer.

---

## 4. The lower boundary is imposed component by component, and a ghost-center radius is used as a face state

### Symptom

Two symptoms in the same cells, and the documents establish they are separate
defects.

(a) The converged hot Uranus carries an alternating velocity in its first cells:
`v` = -248.9, +21.9, -4.91, +3.06, -0.10, +1.06 cm/s over cells 1-6, i.e. a
cell-centered `rho v r^2` of -196, +18.1, -4.24, +2.73, -0.10, +1.00 times the
wind's own flux, while the Riemann face fluxes of the same cells are +1.75 and
+1.92 `F_0`, both positive and both of order the wind flux.

(b) The reservoir state and the resolved atmosphere disagree on the temperature.
Extrapolating the converged interior to the base face agrees with the ghost on
`p` to better than 0.1 percent on every grid and disagrees on `T` by 5-14
percent, growing under refinement, with `rho` taking up the difference.

### Where it arises

* `src/modules/states/Apply_BC.f90:115-175`, `BC_component_constrho`. It pins
  three primitive variables at the ghost:

  ```fortran
  W_in(1,index) = rho_bc                              ! line 123
  ...
  W_in(2,index) = base_flux_const/(rho_bc*r(index)**2)  ! line 128, massflux
  W_in(2,index) = 0.5d0*(W_in(2,1)                      ! line 134, smooth valve
                  + sqrt(W_in(2,1)**2 + valve_eps**2))
  ...
  W_in(3,index) = ntot_bc + dp_bc                       ! line 171, the default
  ```

  Density and pressure together fix the temperature, so with `rho = rho_bc` and
  `p = (ntot_bc + dp_bc) T0` the ghost is isothermal by construction. That is two
  thermodynamic statements.

* `src/modules/states/Apply_BC.f90:225-255`, `Rec_BC`. The reconstructed FACE
  states at the base are overwritten with the same routine:
  `call BC_component_constrho(WR_out,1-Ng)` at line 236 and
  `call BC_component_constrho(WL_out,1-k)` at line 240. But the expressions
  inside it are evaluated at `r(index)`, the ghost CELL CENTER, while the face
  they are used at is `r_edg`. `src/modules/time_step/RK_rhs.f90:60-92` consumes
  them as face states: `call Num_flux(WL(:,j),WR(:,j),Fp,...)` with `rp =
  r_edg(j)`. A cell-centered value used as a face value is a first-order-in-dr
  error, which is exactly the scaling the refinement measurement finds.

* `src/modules/states/Apply_BC.f90:18-91`, the header of
  `base_inflow_mach_number`, states the characteristic accounting the boundary
  is meant to respect: "A subsonic inflow carries one outgoing characteristic,
  so one piece of interior information belongs in the ghost and the other two
  are the reservoir's to state." **The default branch honors that count**: it
  prescribes `rho` and `p` and takes the velocity from cell 1 through the valve.
* The `Base velocity: massflux` branch does not. At `Apply_BC.f90:127-128` it
  prescribes the velocity as well, so all three primitive variables are
  reservoir data while one characteristic is still leaving the domain. That
  branch is the one every `Base velocity: massflux` run uses and it has not been
  checked for characteristic consistency.

### Why it is wrong

**Corrected after review.** The first edition of this catalog said the boundary
was over-specified because it prescribes two thermodynamic variables. That is
wrong, and the review is right to reject it: for outward subsonic flow through
the inner boundary two characteristics enter the domain and one leaves it, so
two reservoir conditions plus one interior compatibility relation is the
expected count, and the default branch supplies exactly that.

What is actually wrong is three separate things.

**(i) The boundary is imposed component by component rather than through the
characteristic invariants.** `BC_component_constrho` writes `rho`, `v` and `p`
independently; nothing in it forms the outgoing acoustic invariant from the
interior, so the interior information reaches the ghost only through whichever
single component the branch happens to copy.

**(ii) A ghost-cell center radius is used in a face state.** `Rec_BC` calls
`BC_component_constrho` to fill reconstructed FACE states, and the velocity
formula inside it evaluates `r(index)` -- the ghost cell center -- while
`RK_rhs` then consumes the result as a state at `r_edg`. This is a
first-order-in-`dr` error at the face, and it is the scaling the refinement
ladder below measures. This is the one part of the first edition's account that
the review confirms in code.

**(iii) The `massflux` branch prescribes all three primitive variables.** For
subsonic inflow that is one condition too many, and it is a different statement
from (i) and (ii). It has not been separately measured.

**What (b) is, and is not.** Prescribing both density and pressure is a
legitimate reservoir model. The measured pressure agreement and temperature
disagreement diagnose an **incompatibility between that reservoir model and the
hydrostatic and thermal state the resolved atmosphere selects** -- not an
incorrect number of boundary conditions.

**What (a) is.** Numerics: the cell-centered velocity of cell 1 is not the
quantity the scheme transports, and any diagnostic that reads `v(1)` -- or the
ghost velocity -- reports a discretization error rather than a flow. The reading
is supported rather than assumed: `|v(1)|` is first order in `dr` to 2 percent
over a 4x refinement (-248.3, -122.2, -60.0 cm/s at `dr` = 1.93e-4, 9.60e-5,
4.78e-5 R_p), its cell count falls 6 -> 4 so it is not a fixed-width physical
boundary layer, and its amplitude falls so it is not an undamped null mode.

### Measured size

**(a) The velocity artifact.** `docs/p44_base_sawtooth.md` section 7 and
`TO_BE_DONE.md` (V). Refinement ladder (`p44:590-592`):

| grid | realized dr [R_p] | v(1) [cm/s] | ratio | cell-centered `rho v r^2` in cell 1 [F_wind] | log10 Mdot | H2 front [R_p] |
|---|---|---|---|---|---|---|
| 1x | 1.930e-4 | -248.3 | -- | -202.8 | 10.28 | 1.0527 |
| 2x | 9.601e-5 | -122.2 | 2.03 | -94.5 | 10.29 | 1.0634 |
| 4x | 4.776e-5 | -60.0 | 2.04 | -43.8 | 10.29 | 1.0636 |

Five boundary prescriptions were tried and none moves cell 1
(`TO_BE_DONE.md` (V); the full table with `||R||` and `info` is
`p44:297-306`). `v(1)` reads -249, -249, -248, -246, -247 cm/s across the
baseline, the Shapiro filter, `Base velocity: massflux`, Low-Mach damping, and
the two combined. Cells 2-12 do respond to dissipation (alternating density
amplitude 2.28e-3 -> 2.84e-4 with the filter while marching, 5.06e-4 -> 1.09e-4
with Low-Mach damping in the converged state); cell 1 does not. Of the arms that
reached `info = 0`, `log10 Mdot` reads 10.29 (A), 10.34 (B), 10.28 (C), 10.29 (E)
and 10.29 (G) against a cycle-to-cycle noise of about 0.01 dex, so only the
Shapiro filter's +0.05 dex is outside the noise, and it is worse on every base
metric. Arms D and F never converged and their 10.35 and 10.47 are marked in the
source document as not quotable.

The layer is not under-resolved: `H(T_eq)/dr = 288` cells, 127 with the run's own
`mu = 2.27 m_H`, against 2.5 for converged HD 189733 b. The mode is not inherited
from the initial condition (all velocities positive and monotone there,
alternating density amplitude 1.5e-13); it appears between steps 24,000 and
30,000 as `du` falls through 0.5 to 0.13, and then holds `v(1) = -248 +- 2 cm/s`
for 150,000 further steps.

The **cell-centered** product of cell 1 is that velocity times `rho_bc`, so it
scales with the base density: -21, -203, -589 times the wind flux at 1, 9 and
28 microbar. Moving the hot-Uranus gate to 1 microbar therefore cut it 9.7x
without touching the artifact, and made the converged molecular layer
reproducible across restart cycles (H2 front to five digits, against a
1.053-1.089 R_p scatter at 9 microbar). The two documents disagree in the
magnitude of that number (-159 to -196 `F_wind` at 9 microbar in `p44`,
-21.1 `F_0` at 1 microbar in `p54`) because they are different states; both are
measured.

**Read the two quantities separately.** `TO_BE_DONE.md` (V) and
`p44_base_sawtooth.md` describe this number as a "base-face mass-flux error" in
places; the review is right that this conflates two quantities, and item 8's
section 11.3 measurement settles which is which. The Riemann FACE mass flux at
the base is `+1.75` and `+1.92 F_0`, positive and of order the wind flux, and
above `r = 1.10` it takes one unique double-precision value over three hundred
cells. The large negative numbers above are the cell-centered product, which is
not the quantity the scheme conserves. In this catalog "cell-centered
`rho v r^2`" and "face mass flux" are never used interchangeably.

**(b) The reservoir/atmosphere thermal incompatibility.** `TO_BE_DONE.md` (W),
`p44_base_sawtooth.md` sections 8.3, 9.2, 9.4. The two closures that let
something float both fail on this planet: `Base ghost temperature: continuous`
converges once to a 1501 K base with the H2 front collapsed onto the boundary
(1.0009 R_p) and does not stay converged on the next restart cycle;
`Hydrostatic base: True` runs the base away to 10,200 K. Raising `T0` to the
1200 K the interior asks for closes the thermal gap (+29.1 -> -8.7 K) and leaves
the cell-1 cell-centered flux at -203.8 against -202.8, which is the measurement
that separates the two defects.

### Remedy options

1. **Solve the lower condition at the face, characteristically** -- the
   review's Phase C and the option this catalog now recommends over the two
   below. Take the outgoing acoustic invariant from the first interior
   reconstructed state, prescribe two independent reservoir quantities with an
   explicit physical meaning (entropy or temperature plus a pressure or density
   datum), solve for the boundary-face primitive state, and construct ghost
   averages that reproduce it to the order of the reconstruction. Handle
   supersonic inflow and flow reversal by an explicit characteristic count
   rather than by passing the interior velocity through a softplus expression.
   This addresses (i), (ii) and (iii) at once and subsumes options 2 and 3.
   Cost: high, and it is a public change to the mathematical boundary-value
   problem, so the reservoir pair and the reversal behavior should be agreed
   before implementation. Verification, from review section 6 Phase C: a
   stationary hydrostatic atmosphere under the same source discretization as
   `RK_rhs`; an outgoing small-amplitude acoustic pulse, measuring the reflected
   amplitude; a smooth manufactured inflow, demonstrating the expected order at
   the boundary face; the hot-Uranus reload-and-march experiment comparing face
   flux, cell-centered flux and the complete update-map defect; then three base
   grids.
2. **Evaluate the boundary expressions at the face radius rather than the ghost
   center where they are used as face states.** The narrowest repair, and it
   removes exactly the first-order term of (ii). Cost: small. Side effect:
   byte identity of every golden breaks. Verification: the refinement ladder --
   if (ii) is the whole of (a), `v(1)` should stop scaling with `dr`. Useful as
   a measurement even if option 1 is adopted, because it isolates how much of
   the artifact (ii) owns.
3. **A staggered or momentum-interpolated velocity at the base face**, the
   standard cure for a collocated odd-even mode, named in `TO_BE_DONE.md` (V) as
   not yet tried. Verification: `v(1)` no longer first order in `dr`; the
   cell-centered and face fluxes at cell 1 agree in sign; `Mdot` and the H2 front
   unchanged within the 0.01 dex and 0.003 R_p cycle noise.
4. **Check the `massflux` branch's characteristic count separately** (iii). It
   is a distinct defect from (i) and (ii) and has never been measured on its own.
   Cost: a measurement, not a change.
5. **For (b), a reservoir model the resolved atmosphere accepts.**
   `TO_BE_DONE.md` (W) proposes a ghost that pins `p` and lets `rho` and `T`
   float together, and records that no such key exists. Note this is a modeling
   choice about what the lower reservoir *is*, not a repair of a wrong condition
   count. Verification: extrapolate the converged interior to the base face and
   require agreement on all three variables, on three grids.

### The prescription, built and measured (section 152, staged)

**Judgment first: the base sawtooth was never a molecular-layer phenomenon; it
was this boundary.** An isolated Euler-plus-gravity probe -- no chemistry, no
radiation, no molecules, no wind -- reproduces the production ladder. The legacy
closure's `R_mom(1)/(rho g)` reads 0.2514, 0.2507, 0.2503, 0.2502 over an 8x
refinement, against an interior truncation error of `4.9e-6`, and drives
`-247.1, -123.2, -61.3, -30.5 cm/s` against the production ladder's
`-248.3, -122.2, -60.0`, to within 1-2 percent.

**The closure.** Two conditions at the face and one from the interior, which is
the characteristic count for subsonic inflow: the lower atmosphere's **pressure
and specific entropy** at the base level, plus the **outgoing acoustic
invariant** of cell 1, `p_b - rho_i c_i v_b = p_i - rho_i c_i v_i` with `c_i` at
cell 1's own `gamma_eff`. The ghosts become volume averages of the hydrostatic
isentrope continued below the face state. `Valve eps`, `Hydrostatic base`,
`Base ghost temperature` and `Base velocity` are retired and refused at startup.

**(i) The sawtooth is gone, not reduced.** On the probe the boundary residual
becomes first order in `dr` (ratios 2.02, 2.02, 2.01) where it was zeroth order
and stuck at 0.2503 -- 1437x to 11762x smaller -- and `v(1)` falls faster than
first order, `-247, -123, -61, -30` becoming `-0.28, -0.04, -0.007, -0.0001`
cm/s. In production the hot-Uranus gate rung's `v(1)` goes **-251.88 -> +14.40
cm/s**, outward, and every one of the ten matrix cases flips the sign of `v(1)`.

**(ii) The reproducibility defect is repaired.** A cold start and a reload land
on the same root: `info = 0`, the same three gate numbers to four figures, the
same `log10 Mdot = 10.34`, agreeing to `3.9e-5` in `rho`, `1.9e-4` in `v` and
`1.3e-5` in `T`, with `v(1) = +14.393` against `+14.395` and `T` at
`r - R_p = 8e-3` equal to 818.01 K in both. That is a direct repair of
`p44_base_sawtooth.md` section 8.4, which had three `info = 0` states of this
configuration putting the H2 front anywhere in 1.053-1.089 R_p and `T` anywhere
in 583-949 K.

**(iii) Two negative results, and they are the honest headline.**

* **It does not repair item 3, and it makes it worse.** Marched 3000 steps, the
  peak-to-peak departure goes 0.4443 -> **0.9284** at `r = 1.005` (2.09x) and
  0.2517 -> 0.3074 at 1.03 (1.22x). The source document offers a mechanism and
  labels it "a hypothesis, not a measurement".
* **A new item, D6: the boundary still reflects sound completely.**
  `|R| = 0.956` against the legacy `0.953`, independent of amplitude
  (`1e-3` to `1e-5` all give 0.956) and nearly of grid. This is not a tuning
  failure: **any pair of static reservoir quantities pins the incoming acoustic
  amplitude**, so `(p, s)`, `(rho, s)` and `(rho, p)` all reflect alike -- the
  pair does not matter, the staticness does. The one form that is transparent
  buys it by corrupting the steady state: weighting toward the invariant takes
  `|R|` to 0.0024 but drives the steady-state test's `R_mom(1)/(rho g)` from
  `1.97e-4` to `1.32e-1`, 460 times the interior truncation error, so it is
  rejected. The earlier prediction that Phase C would damp the base-layer
  ringing of item (AA) is **withdrawn**.

**A third positive result, from section 155.** Phase C also removes most of the
composition inconsistency in the reported residual: the factor by which the
accepted `||R||` understates the written state's own residual falls from **56 to
5.7** on WASP-121b, because the new boundary removes the base cell's dependence
on the previous sweep's composition -- no `n_part_cell1`, no valve. The same
measurement is what retracts the base-cell irreproducibility this catalog
previously attributed to item 4 (item 8): repeating `eval_residual` on an
accepted `Y`, cells 1-3 move by at most a factor two before Phase C and are
reproducible to four significant figures under it.

**What it costs.** Goldens are not byte-identical for any case; this is a model
change. The two WASP cases take 60 percent more steps to reach the same `du`
threshold (11289 -> 18026). The eight molecular matrix rows are fixed-step
snapshots and their `log10 Mdot` shifts (10.43 -> 10.52) are a different point in
a transient, not a different rate -- the one case that is actually solved, the
gate rung, gives 10.34 -> 10.34, and `wasp_full_newton` gives 13.21 -> 13.21.
`lower_profile` is not converged in either build and its 9.56 -> 8.44 should not
be read as a rate at all. `run_fcheck.sh` is clean.

**Not established:** the three-base-grid test (e) was **not run**; every probe
number is from the isolated Euler probe, not the production binary; and whether
the cold start is harder than it would be under section 146 has no baseline.

### Status

Open, `TO_BE_DONE.md` (V) and (W); **the prescription is built, measured and
staged, not applied.** The source document states "Nothing here is applied to
the tree"; the five design decisions D1-D5 were put to the user and answered on
2026-09-03, and D6 was raised by the implementation and is open. Until it lands,
the converged hot Uranus is a solution whose first cell carries a cell-centered
velocity that is not a flux, and every diagnostic reading `v(1)` or the ghost
velocity has to say so. Note in particular that the ghost velocity
printed in the output files of a JFNK-finished run is `eps^2/(4|v_1|)` from the
softplus valve, which is why it reads +0.94 cm/s on a base whose cell-1 velocity
is -249 cm/s, and 0.000 in a marching-only run of the same state.

**Disposition.** Review section 5: "Location mismatch confirmed; interpretation
partly rejected. The ghost-center radius is reused in a face state. Two
reservoir conditions are appropriate for subsonic inflow. Replace component-wise
face overwrites with a characteristic face solve. Treat the cell-centered
velocity and Riemann face flux as distinct diagnostics." **Our judgment:
corrected.** The over-specification reading is withdrawn for the default branch
and replaced by (i), (ii) and (iii) above; the table heading is corrected; the
two flux quantities are separated throughout.

---

## 5. Layer quantities do not converge under grid refinement

### Symptom

Six He/H rungs were re-converged on a base grid refined by a factor 2, and two
of them also by a factor 4. All fourteen states return `info = 0` on both gates
and every one was reproduced by a second restart cycle. On the gate rung -- the
one every published molecular number rests on -- only one quantity settles.

### Where it arises

Not a single code site. This is the observable consequence of item 3: the layer
is not steady, so the quantity being refined is a sample of an oscillation
rather than a converged value.

### Why it is wrong

A grid-refinement sequence measures a discretization error only if the
underlying solution is unique. `docs/supersonic_molecular_base.md` section 14.7
establishes that it is not: restarting each `info = 0` state under pure time
marching swings `rho v r^2` at 1.03 r_base by 16 to 50 percent of the wind value
over 5000 steps on all four states sampled, and 5000 steps is about two of the
2100-step periods P54 measures. So the sequence is measuring the phase of an
oscillation, and the refinement's effect on the *swing* is not systematic.

### Measured size

`docs/supersonic_molecular_base.md:2678-2712`. Convergence criterion: a quantity
counts as converged from the 2x grid on when its 2x-to-4x change is less than a
third of its 1x-to-2x change.

| quantity | He/H = 0.3: 1x / 2x / 4x | ratio | gate rung: 1x / 2x / 4x | ratio |
|---|---|---|---|---|
| log10 Mdot | 10.29 / 10.26 / 10.27 | 0.33 | 10.29 / 10.25 / 10.32 | **1.75** |
| `f(H2) = 0.5` [R_p] | 1.0330 / 1.0308 / 1.0318 | 0.44 | 1.0431 / 1.0375 / 1.0470 | **1.70** |
| `f(H2) = 1e-4` [R_p] | 1.1138 / 1.0922 / 1.0995 | 0.34 | 1.1690 / 1.1390 / 1.1759 | **1.23** |
| n(H3+) peak [cm^-3] | 4.053e5 / 7.535e4 / 7.549e4 | 0.00 | 1.760e5 / 1.803e5 / 1.680e5 | **2.86** |
| T minimum [K] | 312 / 341 / 387 | **1.52** | 424 / 318 / 459 | **1.32** |
| ionization front [R_p] | 2.0788 / 2.0296 / 2.0361 | 0.13 | 2.3762 / 2.3047 / 2.4404 | **1.90** |
| `x2` at cell 1 | 0.98833 / 0.94034 / 0.94027 | 0.00 | 0.95876 / 0.95845 / 0.95837 | 0.29 |
| `\|v(1)\|` [cm/s] | 353.9 / 124.6 / 56.9 | 0.30 | 220.0 / 103.9 / 31.8 | 0.62 |
| closure at 1.03 | 0.39 / 0.97 / 0.75 | 0.38 | 0.52 / 1.00 / 0.39 | **1.25** |

"On the gate rung only `x2` at cell 1 meets it: every other quantity moves
further from 2x to 4x than from 1x to 2x, and fifteen of the seventeen reverse
direction" (`supersonic_molecular_base.md:2696-2705`). The one quantity that
behaves like a discretization error on all six rungs is `|v(1)|`, which is
item 4's artifact.

The first-cell artifact also *relocates* under refinement: at 1x the He/H = 0.3
rung has the cold alternating cell 1, and at 2x that rung is clean while
He/H = 10 has become the alternating one; peak n(H3+) falls by 5.4 on He/H = 0.3
and rises by 1.8 on He/H = 10 (`supersonic_molecular_base.md:2652-2660`).

Cycle-to-cycle scatter of the base energy closure ratio at 1.02/1.03/1.05 r_base
is comparable to or larger than the grid-to-grid difference on several rungs
(`supersonic_molecular_base.md:2747-2756`); on He/H = 0.3 the two 2x cycles read
0.98/0.97/1.00 and 2.55/1.36/1.01.

Cost of the refinement, for whoever decides whether to adopt it: +19 percent per
step, 2.40x wall clock to convergence in the median -- but on the gate rung the
refined run is *cheaper*, 88 minutes to an `info = 0` state against 84 minutes
that reached none at 1x (`supersonic_molecular_base.md:2597-2604`).

### Remedy options

There is no separate remedy: this closes when item 3 closes. What the refinement
does buy on its own, and what it does not, is stated at
`supersonic_molecular_base.md:2793-2800`: it halves the cell-1 artifact on all
six rungs, it converges the gate rung from cold where 1x did not, and it removes
the 37 percent time-averaged base flux deficit on He/H = 0.3; it does not buy a
grid-converged solution, a steady base layer, or a smaller oscillation.

### Status

Open, and the source document is explicit about what was not done: the 4x grid
was run on two rungs only; the oscillation was sampled on four states with
single realizations, no 4x state was sampled, no period was fitted, and the gate
rung itself was not sampled; the regression matrix was not re-run on any refined
grid.

**Disposition.** Review section 5: "Reported measurement. It cannot yet be
assigned uniquely to item 3. Boundary truncation, operator switching, and
sampling of a transient are all possible. Repeat only after items 3 and 4 are
repaired. Compare the same discrete steady definition on each grid, not states
stopped at different oscillation phases." **Our judgment: agree**, and the first
edition's attribution of this item wholly to item 3 is weakened accordingly:
item 4's boundary truncation is first order in `dr` and is a candidate in its
own right, which the `|v(1)|` column of the table above is consistent with.

---

## 6. The base layer's energy budget did not close on the accepted state -- CLOSED

### Symptom (as it was)

In a steady state with conduction and viscosity off,
`heating = radiative cooling + p div v + advection` must hold pointwise. On the
converged 1 microbar He/H = 0.0793 rung (`info = 0`, `||R|| = 2.714e-6`, gate
flux spread 2.499e-3) it held to 0.4 percent from 1.25 to 3.0 r_base and to
4 percent at 1.10, and failed below. The Status section below records how it was
closed; the account in between is the record of the defect.

### Where it arises

* `src/modules/time_step/steady_residual.f90:69-94`, `assemble_residual`, is
  what the accepted state is a zero of; its energy row is
  `R(3,:) = dF(3,:) - S(3,:) - (heat - cool)`.
* `src/modules/time_step/steady_residual.f90:313-372`, `residual_row_scale`,
  and `energy_row_scale` at line 287, is what that row is divided by. Since
  section 143 the divisor is `energy_largest_term(j) = max(|dF_3|, |S_3|,
  |heat|, |cool|, |Sene|)` (stored by `store_row_terms`, `steady_residual.f90:150-170`), i.e. the
  largest term the row itself contains.
* `src/modules/time_step/steady_residual.f90:498-541`,
  `flux_spread_of_state`, is the flux gate, and it measures over `r >= r_flux`
  (`r_flux = 1.2d0` at `parameters.f90:838`), which does not reach the cells in
  question. The comment at lines 509-523 gives the reason and the reason is
  sound; what it does not cover is that a factor-2 gap in the **energy** budget
  over the same cells is not a statement about the launch valve.

### Why it is wrong

Physics. A steady state closes its energy budget everywhere, not above
1.25 r_base. Either the state is not steady there (which item 3 says it is not)
or the terms are being evaluated inconsistently. The measurements have separated
the mass-flux half of this question and not the energy half.

### Measured size

`TO_BE_DONE.md` (AA) and `docs/p23_thermal_budget.md:646-668`. Closure ratio
(heating / losses) on the gate rung: **0.71 at 1.05, 0.52 at 1.03, 0.56 at 1.02**
r_base. The same failure appears in the total-energy form (2.1 at 1.02, 2.0 at
1.03, 1.4 at 1.05, 0.99 at 1.08) and in the integral form over shells (2.19 over
[1.005, 1.05], 1.04 over [1.05, 1.10], 1.005 over [1.10, 1.30], 1.048 over the
whole domain), and it is not removed by using logarithmic derivatives instead of
centered ones.

The cell-center mass flux `rho v r^2` varies by 43 percent over 1.005-1.10 on
the same state and is constant to 0.2 percent above 1.15, while the solver's own
residual reads `||R|| = 2.7e-6`.

**The mass-flux half is closed** (`Update_EXHALE_stage1.md` section 143). The
continuity row's old scale made its residual read the fractional flux error
times the local Mach number -- verified to one percent at ten radii -- so a layer
at Mach 5e-5 was invisible by a factor 2e4; the momentum and energy rows were
mis-scaled the same way by 1e3 to 1e5. With each row divided by its own largest
term, the gate rung converges to a mass flux flat to 1.0e-4 from 1.01 R_p
outward where it was 30 percent out at 1.03, in 14 JFNK iterations against 5,
`log10 Mdot` unchanged at 10.34 and the H2 front at 1.0449 against 1.0425 R_p.
The 30 percent deficit was carrying none of the answer.

**The energy half stands.** Section 143 changes how the energy row is measured,
not what it contains.

**Removed from this item after review: the `Heating_breakdown.txt` diagnostic is
repaired.** The first edition said it under-reports heating by 49 percent in the
first cells and that its patch was not in the repository. Both halves are stale.
`write_heat_breakdown_eq` (`src/modules/radiation/util_ion_eq.f90:1557`) now
passes `f_vibq` and `e_vibq` into `PH_heat_HHe` (lines 1675-1685) and its
`heat_tot` (line 1769) sums photoheating, photoelectron heating, Lyman-alpha
de-excitation, the helium-recombination coupling, both Penning channels,
Lyman-Werner heating, FUV photolysis and molecular reaction heat.
`docs/p23_thermal_budget.md:599-601` marks it fixed and reports agreement with
the solver's own `heat` column to 1.6e-6 on the state used there. Nothing about
the base energy closure changes: the numbers above were computed from the
solver's `heat` and `cool` columns from the start, precisely so that this could
not enter them.

**Identifier note -- collision resolved 2026-09-03.**
`docs/p23_thermal_budget.md` section 9.1 labeled that heating diagnostic
"item P52", but `TO_BE_DONE.md:2385` uses `(P52)` for a different item -- "a run
could not say why it stopped, and 'accepted' had two definitions", closed under
`Update_EXHALE_stage1.md` section 140 -- so the label pointed a reader at the wrong
entry. The heating diagnostic is now called the **heating-dump fix
(`Update_EXHALE_stage1.md` section 140.5)**, which is the subsection that actually
records it being carried in ("A second diagnostic-only defect travels with it").
The scratch file names on disk keep their `p52_` spelling and section 9.1 says
why. `P52` now means only the `TO_BE_DONE.md` item, here and there.

### Remedy options

1. **Establish whether the layer has a steady solution at all** (item 3). If it
   does not, the closure ratio is not a property of a state and there is nothing
   to repair here.
2. **Add an energy-budget gate over the layer.** The item's own point is that
   nothing currently measures the energy closure below the flux window. A
   reported (not gating) closure ratio at 1.02, 1.03, 1.05 and 1.10 r_base
   beside the existing flux spreads costs nothing and would make this visible in
   every run log, exactly as the flux spreads at `r >= 1.03` and `r >= 1.10`
   already do (`EXHALE_main.f90:1099-1109`). Cost: negligible. Verification:
   it reproduces the numbers above on the state they were measured on.
3. **Compare the energy row against the complete production update** rather
   than against the steady residual alone -- the `G_dt` discriminator of item 3.
   A factor-2 energy gap in cells whose momentum row is a near-cancellation of
   two terms worth 1.1e5 times the wind flux may be a splitting or ordering
   defect rather than a physical imbalance, and nothing currently separates the
   two. Cost: it comes free with the item 3 tool.

### Status

**CLOSED 2026-09-03, late** (`docs/p23_thermal_budget.md` section 9.2, marked
SINCE CLOSED, and section 10, which carries the side-by-side). **It was two
errors at once, and neither alone accounts for it.**

1. **The specific heat used to evaluate the advection term.** The budget was
   read off the written profile with a constant `gamma = 5/3`, while the code
   itself now integrates the caloric `gamma_eff` of item (Z). Re-evaluating the
   advection term with `gamma_eff` is the first change.
2. **The state it was read from was not steady in the layer.** The second is to
   read it off a state finished under the face-flux measure of section 145.

Closure `heating/(radiative + adiabatic + advection)`:

| | 1.02 | 1.03 | 1.05 | 1.08 | 1.10 | 1.15 | 1.30 | 1.50 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| section 4 state, `gamma = 5/3` (the original measurement) | 0.56 | 0.52 | 0.71 | 1.01 | 1.04 | 1.02 | 1.00 | 1.00 |
| section 4 state, `gamma_eff` alone | 0.72 | 0.71 | 0.93 | 1.02 | 1.04 | 1.02 | 1.00 | 1.00 |
| P53 cold-converged, `gamma_eff` | 1.09 | 1.42 | 0.75 | 0.93 | 0.97 | 1.03 | 1.00 | 1.00 |
| **face-flux re-finished, `gamma_eff`** | **0.99** | **1.01** | **1.02** | **1.00** | **1.00** | **1.00** | **1.00** | **1.00** |
| face-flux re-finished, `gamma = 5/3` | 0.89 | 0.82 | 0.71 | 0.97 | 1.00 | 1.00 | 1.00 | 1.00 |

**Two caveats, both live.**

* **The section 10 table is marked PROVISIONAL against exactly the change that
  has since landed.** Its header says the state was accepted "as that norm
  stands in the tree at 17:38 KST" and that "if it changes, this state is
  re-finished and the table is re-measured". Section 145 did change it -- to the
  cellwise maximum and the face-flux gate at `2.0e-5`. The state was re-finished
  under the face-flux measure, which is the row that closes; whether it also
  satisfies the cellwise maximum norm at `1e-5` is not stated in that document.
  Until it is re-measured on a state accepted by the landed gates, read the
  closing row as measured-and-provisional.
* **`TO_BE_DONE.md` (AA) has not caught up.** It still reads "the MASS-FLUX half
  is CLOSED ...; the energy half stands." That is the pre-closure status.

The mass half's conclusion is unchanged, and section 145.1 sharpened it: reading
(i) -- "the finite-volume state is steady in its own discrete sense while the
cell-center products are not" -- is now not merely unsupported but inverted. The
conserved face flux of the accepted root is a single double-precision value over
three hundred cells above 1.10 R_p, and the `2.5e-04` the gate used to report was
the difference between two functionals, wrong by six orders.

**Disposition.** Review section 5: "Physical requirement confirmed; cause open.
`assemble_residual` contains the expected energy row, but acceptance can average
a narrow imbalance. The cited heating-output defect is already fixed. Add local
and shell-integrated energy gates. Compare solver terms with the full update map
and then repair the term or boundary causing the mismatch." **Our judgment:
corrected, then closed.** The diagnostic half was already fixed when the review
said so; "acceptance can average a narrow imbalance" was then measured and is
item 8; and the physical cause turned out to be the two errors above rather than
a term or a boundary needing repair. What the review asked for that is still
worth doing regardless is the **shell-integrated and local energy gate**: the
budget now closes on one state, and nothing in the code checks that it closes on
the next one.

---

## 7. Three concepts -- marching stop, solver handoff, scientific acceptance -- share two windows and no names

### Symptom

A run can stop as "converged" on `du < du_th` while its mass flux is still
several percent out over the window the steady solver's own gate uses.

### Where it arises

* `src/EXHALE_main.f90:1024-1032`, the marching stop:

  ```fortran
  mom_max = maxval(abs(mom(j_min:N)))
  mom_min = minval(abs(mom(j_min:N)))
  ...
  du = abs((mom_max-mom_min)/max(mom_min, 1.0d-30))
  ```

  `j_min` is the first cell with `r >= r_esc`
  (`src/modules/init/define_grid.f90:194-198`), and `r_esc` is read from
  `input.inp` (`input_read.f90:181`); the values actually in use across
  `examples/` and `backup/regression/` are 2.00 (47 files) and 1.50 (23 files).

* `src/modules/time_step/steady_residual.f90:498-541`,
  `flux_spread_of_state`, the acceptance gate:

  ```fortran
  do j = j_flux, N
     f    = u(2,j)*r(j)*r(j)
     ...
  spread = (fmx - fmn)/max(abs(fmean), 1.0d-30)
  ```

  `j_flux` is the first cell with `r >= r_flux` (`define_grid.f90:228-233`),
  and `r_flux = 1.2d0` (`parameters.f90:838`).

The two also differ in functional and in sign: `du` is
`(max|F| - min|F|)/min|F|` over absolute values, the gate is
`(max F - min F)/|mean F|` signed; and in threshold, `du_th = 1e-3` against
`flux_spread_th = 5e-3`.

### Why it is wrong

**Corrected after review.** The first edition treated this as one defect --
two definitions of a single statement -- and proposed unifying them. The review
rejects that and it is right to: these are **three different concepts with three
different roles**, and they should not be forced onto one threshold.

* **Marching stop** (`du < du_th`): an economical indication that the time
  evolution has slowed. It is allowed to be cheap and local to the wind window.
* **Nonlinear-solver handoff** (`du < newton_du_switch`): a criterion chosen for
  solver robustness -- when is the state close enough for JFNK to succeed. It
  answers a numerical question and has no physical content at all.
* **Scientific acceptance** (`steady_gates_met`): satisfaction of the physical
  residual and flux gates. This is the only one that licenses calling a state a
  solution.

What is wrong is not that they differ. It is that they are **not named apart,
not calibrated apart, and not reported apart**, so a run stopped only by `du`
prints "converged" and a reader cannot tell which of the three statements was
actually satisfied. The physical error that follows is real: `wasp_full`, a
matrix golden, stops as converged at an `r >= 1.2` flux spread six times the
acceptance gate's own threshold.

The first edition also quoted the project rule "one kind of quantity, one rule"
against this. That rule applies to one kind of quantity; these are three kinds.

### Measured size

`Update_EXHALE_stage1.md` section 143.5 (lines 21902-21946) and
`docs/p55_du_window_census.md`. `wasp_full`, a matrix golden, stops as converged
on `du = 9.999e-4 < 1e-3` at an `r >= 1.2` spread of **3.13e-2**, six times the
flux gate's own threshold, and 1.37 at `r >= 1.01`. P51's transport-on hot Uranus
stopped the same way at 4.4e-2. The two numbers are uncorrelated on the states
these cases stop at: every `armD_*`/`armHeH_*` snapshot ends with `du` between
1.1 and 28 and a gate spread between 0.77 and 35.

A build with the `du` window moved to `[j_flux:N]` was prepared and a 13-case
census run (`p55_du_window_census.md`). Where both builds stop, the answer is the
same: `log10 Mdot` agrees to 0.01 dex on five cases and 0.02 dex on a sixth, and
on the three cases that reach the 60,000-step cap under both builds the state is
identical to four digits. **The window change is a stop criterion change and
moves no physics.** What it costs is the JFNK hand-off, which arms at
`du < newton_du_switch = 1e-2`: on `wasp_full_newton` (which the census
recorded under both its names, the duplicate `solver_newton_cold` having since
been deleted) and `wasp_he23off_newton` the hand-off arms at step ~4242 under
the present window
and at ~55,950 under the wider one, a factor 13.2. And one case gets worse:
`arm_heh1_x2matched` reaches a JFNK finish at a flux spread of 2.7e-4 today, and
under the wider window never arms the hand-off, ending at the cap 218 times less
flat and 0.23 dex higher in the rate. So widening the window without also
re-deriving `newton_du_switch` turns a converged case into an unconverged one.

Caveats the census states about itself: one cap (60,000 steps) with no
re-derivation of `du_th` or `newton_du_switch`; the runs that did not stop are
"had not stopped by then", not "diverged"; and all of it was measured on the
pre-143 mass-flux row scale, not on the measure now in the tree.

### Remedy options

1. **Name the three concepts apart and say which one stopped the run** -- the
   review's Phase F, and the option this catalog now recommends. Give each its
   own name and its own calibrated tolerance, and make the stop message state
   which criterion fired and what the other two read at that moment. A run
   stopped only by `du` must not print "converged" or be labeled a steady
   solution. Cost: small; it is reporting and naming, not a change to any stop.
   Side effect: none on any number. Verification: every matrix case's stop line
   names its criterion, and `wasp_full`'s says "marching stop, `du`", not
   "converged".
2. **Recalibrate the handoff separately, and only for solver robustness.** The
   census shows the handoff is what actually breaks under a window change:
   `arm_heh1_x2matched` reaches a JFNK finish at a flux spread of 2.7e-4 today
   and under a wider window never arms the handoff at all, ending at the cap
   218 times less flat and 0.23 dex higher in the rate. That is a solver
   question, to be tuned against solver success, not against a physical
   criterion.
3. **Require the acceptance gates for any publishable state**, independently of
   how the run stopped. Cost: none in code -- `steady_gates_met` already exists
   and is already the single definition; what is missing is the discipline of
   not quoting a `du`-stopped state.
4. **Withdrawn: unify the window and the functional.** The first edition
   proposed making `du_th` and `flux_spread_th` one number. The review rejects
   it and the census supports the rejection. Recorded because section 143.5
   still carries the proposal.
5. **Already in the tree:** the gate spread and the `r >= 1.03` and `r >= 1.10`
   spreads are printed beside `du` on every diagnostic line
   (`EXHALE_main.f90:1099-1109`), which makes the mismatch visible without
   moving any stop.

### Status

Open; the naming and stop-message work is not done. The unification proposal is
**withdrawn**. `Update_EXHALE_stage1.md:21904-21906` records the tree's present
position: "the tree keeps `du` as it is. It is recorded because the census makes
it unavoidable"; `Update_EXHALE_stage1.md:21942-21946` states why the unification was
left to a user decision: "It is a different change with a different blast
radius: it moves where every marching run stops ... Section 143 changes what a
converged state must satisfy; this would change where an unconverged run may
stop." Note that `TO_BE_DONE.md` (P52), closed under `Update_EXHALE_stage1.md`
section 140, already did part of Phase F's job -- it made a run say why it
stopped, after a residual stop had been printing `converged: momentum constant
(du < du_th)` at `du = 3.156e-3` with `du_th = 1.0e-3`. What section 140 did not
do is give the three concepts three names.

**Disposition.** Review section 5: "Different definitions confirmed, but not
necessarily one defect. They serve termination, handoff, and acceptance roles.
Do not force one threshold onto all roles. Rename and document them by role;
require the steady gates for a publishable state. Recalibrate handoff
separately." **Our judgment: corrected.** The unification proposal is withdrawn
and replaced by naming, separate calibration, and stop-message labeling.

---

## 8. What the residual norm hides, and where the layer's `3.4e-3` comes from -- SETTLED

### Symptom (as it was posed)

The hot-Uranus steady root re-solved with the current tree passes both gates
(`r >= 1.2` spread 2.501e-04, `||R||` 9.958e-06) while its own mass flux was
reported 3.433e-03 out at `r >= 1.03`, fourteen times the gate's number. The
first edition asked whether the volume weighting of the residual norm averages
the base layer away.

### Settled, and not the way it was posed

`docs/p55_base_mode.md` section 11, written after the first edition of this
catalog, answers both halves and neither answer is the expected one.

**(1) The `3.4e-3` is not a conservation failure. It is the difference between
two functionals.** The mass row's scale is
`max(|F_in| r^2, |F_out| r^2)/dV`, so `R_1/s_1` is exactly the fractional change
of the conserved face flux across the cell, and summing it over a window bounds
that flux's total relative variation. Measured on the same accepted states:

| case | window | face-flux spread | `sum R_1/s_1` (the bound) | cell-centered spread (**the gate**) |
|---|---|---|---|---|
| hot Uranus | `r >= 1.03` | 1.21e-08 | 9.33e-08 | **3.433e-03** |
| | `r >= 1.10` | 0 (one unique value) | 4.38e-09 | 2.375e-03 |
| | `r >= 1.20` | 0 | 3.17e-10 | 2.501e-04 |
| WASP-121b | `r >= 1.03` | 0 | 1.29e-09 | 6.751e-05 |
| | `r >= 1.20` | 0 | 3.16e-10 | 5.242e-05 |

Above `r = 1.10` the face mass flux takes **one unique double-precision value
over three hundred cells** (`4.94565219e-05` on the hot Uranus). It varies only
in cells 1-3, where it falls `4.80e-04 -> 3.49e-04 -> 4.95e-05` -- item 4's
artifact. So the wind carries one mass flux to the last bit while the gate
reports `2.5e-04`, because `flux_spread_of_state`
(`steady_residual.f90:498-541`) measures the cell-centered `rho v r^2`, which is
not the quantity the finite-volume scheme conserves. The number it returns on
this state is a reconstruction difference, not a failure of conservation.

This also dissolves an apparent disagreement between two of our own documents:
`docs/p54_base_layer_mass_flux.md`'s `1.043e-04` and section 144.5's
`3.433e-03` are the **face** and **cell-centered** spreads of nearly the same
state. On the conserved functional the newer root is better by four orders
(1.043e-04 -> 1.213e-08 at `r >= 1.03`; 6.077e-06 -> 0 at `r >= 1.10`), and that
improvement is invisible to the gate.

**(2) The dilution is real, and it hides the base cells and the wind, not the
layer.** Dumping the solver's own `F` at the acceptance point and dividing each
cell by its own `residual_row_scale` -- the quantity `rnorm` is a norm of --
gives, for cells above the `1e-5` tolerance:

| state | 1.00-1.30 | of which 1.03-1.30 | whole column |
|---|---|---|---|
| hot Uranus | 1 mass, 0 momentum, 6 energy | **0, 0, 0** | 1, 0, 6 (all in cells 1-3) |
| WASP-121b | 0, 0, 20 | 0, 0, 18 | 0, 0, **70** |

Between `r = 1.03` and `1.30` not one cell of either row exceeds `1e-5` on the
hot Uranus; the layer's own scaled residual there is `1e-11` to `5e-6`. **The
layer is not hidden by the norm.** What is hidden is:

* on the **hot Uranus**, cells 1-3 -- item 4's cells -- where the mass row reads
  `3.10e-05` and the energy row `2.17e-03` against a gate of `1e-5`, i.e. 218
  times the accepted `||R|| = 9.958e-06`. The run's own log even names the cell
  (`max cell j=3 r=1.0006 k=3`) and prints the diluted number beside it;
* on **WASP-121b**, seventy cells, mostly in the **WIND**: 50 in
  `r = 1.30-2.00`, 18 in `1.10-1.30`, 2 at the base, at `1.0e-05` to `6.4e-05`.
  There the volume weighting is diluting a wind, not a layer.

**(3) The review's formulation concern is exactly right, and sharper than our
hypothesis.** Neither branch of `relnorm_over_cells`
(`steady_residual.f90:411-440`) is a cellwise maximum of relative residuals. The
default returns `sum_j |R_kj| V_j / sum_j s_kj V_j`, a scale-weighted average of
the cellwise ratios; the alternative returns `max_j |R_kj| / max_j s_kj`, a
ratio of two maxima taken independently and **not in the same cell**, since
`s ~ 1/dr` is largest in the smallest base cells while `|R|` need not be. The
usual local relative infinity norm, `max_j (|R_kj|/s_kj)`, is what the comment
at `residual_norms` describes and what neither branch computes.

**What the first edition got wrong.** Its arithmetic -- that `scale * V` is
roughly uniform for the mass row, so the weight of each cell is the cell-count
fraction (0.53 for `r < 1.10`) rather than the volume fraction (0.048) -- is
correct as far as it goes and led to the right conclusion for the wrong reason.
The layer is indeed not diluted; but the reason the gate's number looked bad was
the functional, not the weight, and the weight *is* hiding something, three
cells lower down.

### The repair, adopted

Both halves are the user's decision, taken 2026-09-03.

**(a) The acceptance norm becomes the maximum over cells of a cell's own scaled
residual**, and it becomes the only form: `resid_vol` and the `Resid norm` key
are removed with the second form they selected. Measured on both states
(`p55_base_mode.md` section 11.5):

| | hot Uranus, default | hot Uranus, **cellmax** | WASP-121b, default | WASP-121b, **cellmax** |
|---|---|---|---|---|
| JFNK `info` | 0 | **0** | 0 | **0** |
| JFNK iterations | 7 | 9 | 8 | 12 |
| cell-max, mass row | 3.100e-05 | **2.021e-07** | 3.232e-07 | **4.374e-08** |
| cell-max, energy row | 2.174e-03 | **2.002e-06** | 6.397e-05 | **9.764e-06** |
| cells above `1e-5` | **7** | **0** | **70** | **0** |
| `log10 Mdot` | 10.34 | 10.34 | 13.21 | 13.21 |
| 3000-step departure, p2p at `r=1.03` | 0.5333 | 0.5334 | 0.0103 | 0.0104 |

It converges on both, in two to four extra JFNK iterations, leaves no cell of
any row above the tolerance, and changes nothing else: the rate is identical to
the printed digit, the flux spreads agree to five digits, and the 3000-step
departure is the same to 0.2 percent. That last line is a result in its own
right -- **the norm is not what makes a root leave the marching loop**, so the
third candidate for item 3 is excluded.

**(b) The flux gate moves to the faces.** The measurement says the functional
the gate uses is not the conserved one. `docs/p54_base_layer_mass_flux.md`
recommendation 7(a) said so and was never carried in; section 11.3 is a much
sharper argument for it.

### Where the repair stands, verified in the tree

* **(a) is in patch, not in the tree.** `EXHALE_RESID_NORM=cellmax` exists only
  as `p55_cellmax.diff` in the P55 scratch directory (three files, all
  environment-gated, all off by default). The section 145 patch that makes it
  the single form is `apply145_resid.py` and `apply145_main.py` in the same
  scratch tree, built into `t145/`. In the working tree,
  `parameters.f90:798` still reads `logical :: resid_vol = .true.` and
  `relnorm_over_cells` still carries both branches; `Update_EXHALE_stage1.md` has no
  section 145 yet.
* **(b) is decided but not in the patch.** `t145/src/modules/time_step/steady_residual.f90`
  lines 518 and 552 still form `f = u(2,j)*r(j)*r(j)`, the cell-centered
  product, in both `flux_spread_of_state` and `flux_spread_above_radius`. The
  face-flux gate is adopted and not yet written.
* The tolerance a face-flux gate should carry is **not** settled. Every
  threshold in the code is calibrated against the present functional, and
  `p55_base_mode.md` section 11.7 says so explicitly.

### RETRACTED: the base cells are reproducible, and `calc_rho` does not move `rho`

The first and second editions of this catalog carried a further finding from
`p55_base_mode.md` section 11.6 -- that re-assembling the residual on the same
accepted `u` moved cells 1 and 2 of the MASS row by a factor 30,000, attributed
to `ioniz_eq` rewriting `rho` through `calc_rho`, and read as item 4 in a new
form. **Both halves are withdrawn by measurement** (section 155, staged;
`p55_base_cells_155.md`).

**`calc_rho` does not move `rho`.** Over every sweep of two roots on two builds
(346, 86 and 57 sweeps), the density the sweep writes back differs from the
density it was handed by at most `1.7e-15`, with `n_H` and `n_He` at `4.7e-16`.
That is round-off, and it is round-off *for a reason*: the sweep's unknowns are
fractions of the element totals it computes from the state it is handed, and the
code's masses are the nucleus counts exactly (H2 = 2 m_H with 2 H nuclei,
H3+ = 3 with 3, HeH+ = 5 = 4 + 1), so **conserving nuclei conserves mass
identically**. There is one place where that is not closed by construction --
the `max(...,0.0d0)` clamps on `nhi`/`nhei` in the molecular write-back -- and on
these states it does not bite.

**The base cells' residuals are not O(1) either.** Repeating `eval_residual` on
one accepted `Y`, first evaluation against the twentieth, the largest excursion
anywhere in cells 1-3 is **a factor two** (WASP-121b cell 1 before Phase C,
0.545); under Phase C that same cell is reproducible to four significant figures.
So the base cells are not the problem, and to the extent they were one, Phase C's
boundary closed it -- it removed the base cell's dependence on the previous
sweep's composition, with no `n_part_cell1` and no valve.

### What is real, and it is a different defect: the accepted residual is evaluated at the previous iterate's composition

`f_sp` is `intent(inout)` in `ioniz_eq` and the residual mutates it, and it is
not part of the Newton unknown `Y`. So `eval_residual(Y, f_sp)` is a function of
two arguments and the state carries only one of them. It is **deterministic in
both** -- with the composition restored before every repetition, all four states
give bitwise identical residuals, and nothing else (frozen chemical background,
excited-hydrogen populations, stored WENO smoothness factors) contributes
anything at all.

**But the number reported at acceptance is evaluated at the wrong composition.**
In the line search the accepted trial is evaluated as
`eval_residual(Ytry, f_sp_j)` with `f_sp_j` copied from the *previous* iterate,
and `rnorm` is taken from that `F`, while the state written to disk carries the
composition the accepted trial left. Evaluating the same `Y` at its own
composition, and iterating the sweep at fixed `Y` to its fixed point:

| build | state | reported `\|\|R\|\|` | own composition, 1 sweep | at the composition fixed point | factor |
|---|---|---|---|---|---|
| Phase C | hot Uranus | 2.211e-06 | 1.283e-05 | 1.257e-05 | **5.7** |
| Phase C | WASP-121b | 6.287e-06 | 3.865e-05 | 3.593e-05 | **5.7** |
| pre-C | hot Uranus | 4.141e-06 | 1.762e-05 | 1.722e-05 | **4.2** |
| pre-C | WASP-121b | 6.264e-06 | 1.737e-04 | 3.539e-04 | **56** |
| pre-C | `wasp_full_newton` | 8.873e-06 | 1.851e-04 | 3.412e-04 | **38** |

The composition has its own fixed point at fixed `Y`, reached in 3 to 14 sweeps.
**Phase C cuts the WASP-121b factor from 56 to 5.7** -- a second, independent
argument for item 4's boundary.

**What this does and does not change.** It does **not** change which states are
accepted today: only 7 of the 63 regression cases set `Resid tol` at all,
without it the residual gate is off, and every root measured here was accepted on
the **flux gate alone** -- which reads the conserved face mass flux of `u` and is
composition-independent, so it is untouched. What changes is what the reported
number *means*: `||R||` is quoted in run logs, in `docs/`, and in item 3's own
measurements, and on these states it understates the residual of the state that
was written by 4 to 56 times. **Every one of the five states above is above
`1e-5` at its own composition**, so a run that did set `Resid tol: 1.0e-5` would
be accepting a state that does not meet it.

**Adopted (section 155, block I):** `resid_at_own_composition`, which
re-evaluates at hand-back and iterates the sweep at fixed `Y` to the fixed point,
**default off with the recommendation to turn it on** -- it is a decision about
what "converged" means, not a refactor, and turning it on makes every root
accepted under `Resid tol` report a residual 4 to 56 times larger. Cost is 3 to 9
sweeps, one-off. And `rho` becomes an input to the sweep rather than an output:
that changes no number today (the write-back is element-conserving to `1e-15` on
every state tested) but it makes the property structural instead of accidental,
and turns the `max(...,0.0d0)` clamps into a detectable condition rather than a
silent mass adjustment. It is a golden-refreshing change at the `1e-15` level, so
it belongs in a series that is refreshing them anyway.

### Status

**Settled**, with one part of the account retracted and replaced. The face-flux
gate and the cellwise norm landed in section 145; the retraction above removes
the base-cell irreproducibility from this item and puts a real defect in its
place, staged as section 155.

Still open: the face-flux gate's threshold calibration; whether the cellmax norm
passes the regression matrix (it was run on two states, not the matrix); whether
the composition fixed point at fixed `Y` always exists and is unique (reached in
3 to 14 sweeps on all five states, but not proved); whether accepting on the
self-consistent residual changes a mass-loss rate (it cannot change the state
written, but a run that keeps solving instead of stopping ends somewhere else);
and whether the `max(...,0.0d0)` clamps ever bite -- shown not to on five states,
not shown never to.

**Disposition.** Review section 5: "Open, with a verified formulation concern.
Both current norm branches differ from a cellwise maximum of relative residuals.
Add dual integrated/local gates and contribution diagnostics before changing
weights or split radii." **Our judgment: settled by later measurement**, and the
review's formulation concern is confirmed in code and is the thing that was
repaired. One divergence from the review's recommendation is deliberate and is
the user's decision: the review proposes **dual** integrated and local gates
with separately calibrated tolerances, while section 145 makes the cellwise
maximum the **single** form and removes the integrated one. The argument for
removal is that the integrated norm's only demonstrated effect on these two
states was to hide seven cells and seventy cells respectively; the argument for
keeping it is the review's, that an integrated norm measures a column budget
that a cellwise norm does not. That trade is recorded here and was decided in
favor of the single form.

## 9. The fluorescent trapping inside the Lyman-Werner cross section carries a plane-parallel geometry our layer does not have

### Symptom

None visible in a run: the rate is smooth and the table interpolates. The defect
is that a geometry-dependent factor computed for a slab is applied to a
spherical, outward-open layer, so the H2 photodissociation rate at depth carries
an error of the size of the geometry difference.

### Where it arises

* `src/modules/lower_atmosphere/h2_self_shielding_table.f90:93-98`, the module's
  own statement of what is not line by line:

  ```
  ! p_eff/p_single -- the suppression of the branching by re-absorption of
  ! the fluorescent decay photons -- is taken node by node from the same
  ! [CLOUDY runs] ... Reimplementing CLOUDY's escape probabilities would be a
  ! second, unverified version of the same thing.
  ```

* The three tabulated cubes and their accessors:
  `h2_shield_log_sigma_diss` (used by
  `h2_lw_dissociation_cross_section`, line 1027),
  `h2_shield_log_p_single` (line 457, accessor line 1039) and
  `h2_shield_log_p_eff` (line 727, accessor line 1049).
* The one consumer of the trapped cross section is
  `src/modules/lower_atmosphere/lyman_werner.f90:507-520`,
  `lyman_werner_dissociation_rate`:

  ```fortran
  k = F_LW/e_lw_photon_erg                                           &
    * h2_lw_dissociation_cross_section(N_H2, T, n_H)
  if (tau_cont .gt. 0.0d0) k = k*exp(-tau_cont)
  ```

  `src/modules/radiation/util_ion_eq.f90:322` reads `p_single` for the heat
  return.

### Why it is wrong

Physics. An escape probability is a geometric quantity. CLOUDY computes these in
a plane-parallel slab illuminated on one face and closed on the other; our layer
is spherical and open outward. Applying a factor computed in one geometry to
another is a model mismatch whatever its size.

**Corrected after review: the sign is not established.** The first edition
argued that because the thick-inward side is the same and dominates, what is
wrong is the escape the open outward side allows, so the code must trap too much
and under-dissociate. That argument is not supported by the implementation or by
any two-geometry calculation. Angular path lengths, spherical dilution, the
inward optical depth and frequency redistribution all enter the escape
probability, and they do not all push the same way. **Until slab and spherical
transfer are compared with identical molecular data, both the sign and the
magnitude of the error should be reported as unknown.**

### Measured size

**No measurement of the error exists**, in either sign or magnitude. What is
measured is how large the affected factor is, which bounds how large the error
*could* be. `TO_BE_DONE.md` (P47): the trapping enhancement itself is 1.13-1.38
across the tabulated grid, rising to 1.25-3.07 at `N_H2 = 5e21 cm^-2`, so an
error of a few tens of percent in it would be an error of the same size in the
rate at depth. And the factor does depend on the geometry-sensitive quantity: at
fixed `T = 1300 K` the converged face cross section moves 13 percent between
`n_H = 1e12` and `1e14` while the untrapped value is identical to five digits,
i.e. it responds to the column behind the point, which is exactly what changing
the geometry changes.

`p_single`, which the fluorescence heat return uses, is free of this: a trapped
decay and the re-absorption that follows it move one molecule down and another
up, depositing nothing, so the trapping cancels out of the count of decays per
dissociation. The module header and the accessor comment both say so.

Related and separately open: `E_bound` is an unshielded average and is 4-9
percent low inside the shielded layer.

### Remedy options

1. **Label the provenance now.** Keep the present table and give it an explicit
   model name -- the review proposes `slab_surrogate` -- with its geometry
   recorded beside it and a stated uncertainty, so no downstream reader treats
   it as a spherical calculation. Cost: negligible. This should be done before
   anything else on this item.
2. **A matched slab and spherical comparison with identical molecular data.**
   This is what determines the sign and the size, and it must come before any
   replacement table. The replacement must recover the present slab result in
   the slab limit and converge in angle, frequency and radial resolution.
   Cost: high; it is a radiative-transfer calculation, not a table lookup.
   Verification: the slab limit reproduces CLOUDY to the digitization.
3. **Do not reimplement CLOUDY's slab escape probabilities.** The module already
   records why: it would be a second, unverified version of the same thing.
4. **Do not apply a correction of assumed sign.** Withdrawn from the first
   edition.

### Status

**Updated 2026-09-03, late: the two-geometry comparison was done, and the sign
IS one-signed after all -- in the direction the first edition guessed, with a
bound the first edition did not have** (section 150, staged;
`e2lw/docs/e2_lw_geometry.md`).

**The bound.** `f_trap = p_eff(spherical)/p_eff(slab) <= 1` everywhere and
monotonically in height: the spherical layer lets more fluorescent photons out,
so `p_eff` and hence `sigma_diss` are smaller, so **the installed table
over-predicts the Lyman-Werner rate**. By 0.2 percent at the base, **7-10 percent
through the H2-bearing layer**, up to 20 percent at the top where H2 is already
negligible, and **6.5 percent weighted by the dissociation rate**. The slab limit
was validated first: `sigma_pump` reproduces the generator to seven significant
digits and `p_single` to 0.1-0.5 percent.

**Why the table is nonetheless NOT rebuilt, and this is the right call.** The
two independent trapping calculations differ from each other by 20-40 percent,
and not with one sign -- more than the geometry does. Rebuilding the cube in
spherical geometry would move `sigma_diss` by less than the uncertainty already
inside it. The measurement now stands in the headers of
`h2_self_shielding_table.f90` and `lyman_werner.f90` in place of the open-item
statement.

**Two further geometric factors, measured and not adopted.** `f_curv` (curvature
of the pumping path) is 0.98-1.06 through the H2 layer. `f_shell` is 0.20-0.33
there: the code evaluates the substellar ray and applies the full band flux at
the planet, so **its Lyman-Werner rate is 3.8-4.4 times the shell average a 1-D
spherically symmetric H2 budget needs**. A `Lyman-Werner shell average` key
exists and is default False, because it changes what a converged solution means
and is calibrated on one planet.

### A separate defect found in the same work, and fixed: three dayside conventions in one run

`fuv_band_flux` returned all five FUV band fluxes **raw**, with no `appx_mth`
dayside factor, while `set_energy_vectors.f90` halved `F_XUV` on an exact string
match and `lya_rt.f90` halved the stellar Ly-alpha beam on a substring test. The
H2 dissociation rate did not even go through `fuv_band_flux` -- it read
`F_LW_star` directly -- so the band had one owner for the oxygen photolysis and
the photon ledger and another for the rate itself. Under
`2D approximate method: Rate/2 + Mdot/2`, the setting of the whole regression
matrix, **the H2 Lyman-Werner rate stood a factor 2 above the convention the same
run applied to every other band.** Section 150 makes it one definition,
`dayside_dilution()`, read by all three paths.

**It moves results.** Of the ten matrix cases only `mol_lyman_werner` carries an
FUV band flux, and on its 12,000-step snapshot the H2 front at `x(H2) = 1e-3`
moves 1.07299 -> 1.08122 R_p, the H2 column rises 13.4 percent, `x(H2)` rises by
up to 4.28x at `r = 1.081`, and `log10 Mdot` reads 10.47 -> 10.44. Those are
snapshot values on a fixed-step relaxation, not converged rates. `wasp_full` and
`mol_base_handoff` are byte-identical in every output file.

### Status

`TO_BE_DONE.md` (P47) is marked **CLOSED 2026-09-03** in the section 150 draft;
the FUV beam defect is fixed in the same patch. **Both are staged, not applied.**

**Disposition.** Review section 5: "Model mismatch confirmed; sign and size
unverified. Add an explicit geometry provenance flag and uncertainty. Implement
a two-geometry transfer comparison before replacing the table." **Our judgment:
agree, done, and the outcome inverts our own correction.** The first edition
claimed the slab over-traps; the review rightly refused that as unestablished;
the two-geometry comparison was then run and it establishes exactly that, with a
bound. **The claim was right and the reasoning behind it was not** -- the sign
does not follow from "the thick-inward side dominates", it follows from the
measured `f_trap <= 1` monotone in height. The review's remaining
recommendation, provenance and uncertainty rather than a replacement table, is
what was implemented.

**What is still not established:** whether CLOUDY's `p_eff/p_single` or this
one is right in absolute terms; complete frequency redistribution and a
single-flight escape probability are assumed; one planet, one rung; and which
`f_shell` convention EXHALE should adopt is not decided.

Note that `mol_lyman_werner`'s golden is currently **deliberately stale**: the reservoir-weighting correction of
`Update_EXHALE_stage1.md` section 135 moved the H2 front from 1.146163 to 1.144327 R_p
and `make check` reports FAIL on that case by design, the refresh being a
separate decision.

---

## 10. The steady residual fills the base ghost before its ionization sweep

### Symptom

The ghost pressure the reconstruction reads is not the pressure the boundary
condition asked for. The gap vanishes at the fixed point and is largest during
relaxation, which is exactly when the solver is trying to make progress.

### Where it arises

`src/modules/time_step/steady_newton.f90`, `eval_residual` (line 427). The order,
with the tree's line numbers:

```
527   call get_species_densities(...)   ! composition of the PREVIOUS sweep
529   call Apply_BC(u)                  ! fills the ghosts of u FROM that composition
531   call U_to_W(u, W)
533   call get_species_densities(...)
535   call comp_T_from_p(p, n_tot, ne, T)
537   call ioniz_eq(T,rho,f_sp,heat,cool,eta,sweep)   ! composition re-solved HERE
542   call get_species_densities(...)   ! the composition the residual will use
545   call assemble_residual(u, n_tot + ne, heat, cool, R)  ! reads the ghosts of u
```

`Apply_BC` (`src/modules/states/Apply_BC.f90:93-111`) states the base condition
in primitive variables and converts to conservative to store it
(`U_to_W_interior` -> `Apply_BC_W` -> `W_to_U`); the reconstruction downstream
converts back.

### Why it is wrong

Under a constant adiabatic index the primitive-to-conservative round trip returns
the pressure the boundary condition asked for, so the condition is enforced
exactly and the lag cannot show. Under the caloric equation of state now in the
code (`src/modules/states/caloric_eos.f90`) the two conversions use different
heat capacities, because the composition was re-solved between them. The gap is
proportional to how far the ghost composition moved in one sweep.

This is a property of the iteration, not of the solution. The marching loop does
not have it: there `Apply_BC(u)` at `EXHALE_main.f90:940` follows the ionization
sweep and `comp_p_from_T`, so the ghost is written with the composition the next
step reads.

### Measured size

Not measured. `Update_EXHALE_stage1.md` section 141.6 quotes no residual, iteration or
percentage figure for the lag itself, and bounds it only by where it vanishes.
`n_part_cell1`, which the `base_ghost_T_continuous` branch divides the ghost
pressure by (`Apply_BC.f90:169`), has had the same lag all along and the code
already said so; what changed with the caloric equation of state is that the lag
now reaches the ghost pressure itself.

### Remedy options

1. **Move the boundary condition after the ionization sweep.** Direct, and it
   makes the residual's ghost the ghost of its own composition. Cost: it is a
   change to the residual, so it belongs with a golden refresh and not inside
   another change. Verification: the whole matrix, with the moved cases
   explained.
2. **Make `Apply_BC`'s interior round trip exact** -- write the interior back
   untouched instead of reconstructing it -- then a second `Apply_BC(u)` after
   the last composition refresh is free. This is the cleaner repair and it also
   removes a latent non-idempotency. Cost: it breaks byte identity on its own,
   because the round trip today evaluates `0.5*W(1)*W(2)**2`
   (`src/modules/functions/UW_conversions.f90:64,68`) where the input was
   `0.5*u(2)*u(2)/u(1)`, and those are not the same double, so a second call
   perturbs every interior cell at the last bit and the whole atomic regression
   matrix moves. Verification: a golden refresh with the last-bit movement
   documented.

### Status

**Updated 2026-09-03, late: repaired, and it took two changes, not one**
(sections 146 and 147, staged).

**Section 146 fixes the ordering, and it is not what makes `F(Y)` a function of
`Y`.** `Apply_BC` now writes only the ghost indices -- `U_to_W_interior`,
`Apply_BC_W`, then `W_to_U_comp` on the `2*Ng` ghosts alone -- instead of
reconstructing the whole array, and `eval_residual` calls it a second time
*after* its ionization sweep, so the conservative ghost state the fluxes read
carries the composition the fluxes use. The measured reproducibility, three
evaluations of `F(Y0)` with trials in between:

| how the sweep is seeded | entries differing (of 1500) | max, in row-scale units |
|---|---|---|
| one shared workspace, as the interface allowed | 1406 and 1400 | **2.968e-03** |
| the same, with 146's ordering fix | 1412 and 1399 | 4.239e-03 |
| **seeded from one composition every time** | **0 and 0** | **0** |

Against an acceptance tolerance of `1e-5`, `3.0e-03` is 300 times out. **The
ordering fix alone does not help** -- it is section 147's seed interface that
does it.

**Section 147 makes the interface enforce it.** `eval_residual` took the
composition as one `intent(inout)` array that was both the seed and the
destination, so whatever the caller's workspace last held became the next
evaluation's initial guess -- the residual was a function of its argument only
because five call sites remembered to make it one. It now takes
`f_sp_seed (intent(in))` and `f_sp (intent(out))` with `f_sp = f_sp_seed` as its
first statement; Fortran forbids aliasing the two, so the failure mode is no
longer expressible. The five call-site copies are deleted. **0 of 1500 entries
differ.** The control -- the first evaluation of a run against the second --
differs in 6 entries, cells 1 and 2 only, by `8e-12` to `7.8e-11`: one-time
table building, not history.

**Cost.** Section 147 moves nothing: both named Newton cases re-solve to
byte-identical output and all ten matrix cases are byte-identical. Section 146
moves everything at the last digit -- all ten cases differ, worst
`8.2e-04` in `rho` on `mol_metals` -- and **needs a golden refresh**; no quoted
number moves. Isolated on a one-step cold start the `Apply_BC` change alone moves
`p` by `1.26e-03`. Two side effects worth naming: the hot-Uranus JFNK takes 36
percent more iterations and WASP-121b 118 percent more; and `v(1)` goes
`-15.02 -> -251.88 cm/s`, i.e. **section 146 restores the sawtooth that section
145's face-flux measure had been masking** and that section 152 then removes.
The item 3 departure improves by 44 percent at `r = 1.03` and 73 percent at
1.005, and the source document explicitly does not separate how much of that is
the consistency and how much the cell-1 displacement.

**A caveat on our own item 3 text.** Section 145.6 measured `Apply_BC`'s round
trip moving the interior by a relative `4.2e-16` and concluded "item 10 is not a
source of motion here". Both are right and they are about different things: that
was the round trip's contribution to `G_dt` on an accepted state, where the
composition barely moves; section 146's `1.26e-3` is a cold start, where it
does.

### Status

**Repaired, staged, not applied.** `TO_BE_DONE.md` (AB) is drafted as "CLOSED
2026-09-03, `Update_EXHALE_stage1.md` section 146" but the draft is not applied and the
tree still carries the pre-146 ordering. Apply order is 146, then 147, then one
golden refresh. The lag is stated at the head of `caloric_eos.f90`.

**Not settled by either section:** the ghost cells' own composition is still
solved at the pressure the *first* `Apply_BC` wrote, deliberately, because
fixing it would cost an extra ionization sweep; the marching loop's `ioniz_eq`
is untouched; and nothing here makes the sweep's root seed-independent -- how far
the root moves with the seed is not measured.

**Its priority was argued from a finding that has since been retracted, and the
conclusion survives on better evidence.** The second edition cited
`p55_base_mode.md` section 11.6 -- cells 1 and 2 of the mass row moving by a
factor 30,000 on re-assembly. Section 155 withdraws that (item 8): the base
cells reproduce to within a factor two, and to four significant figures under
Phase C. What is real is the same conclusion by a different route: `f_sp` is
`intent(inout)` in `ioniz_eq` and is not part of `Y`, so `eval_residual(Y, f_sp)`
is a function of two arguments and the state carries one. **It is deterministic
in both**, which is what section 147's seed interface then exploits. The correct
acceptance condition for the repair is therefore still the review's -- repeating
`eval_residual` at the same `Y` must give the same residual regardless of prior
trial history -- and section 147 meets it.

**And section 155 adds a second, separate defect on the same interface**: the
residual the gate reports is evaluated at the PREVIOUS iterate's composition,
understating the written state's own residual by 4 to 56 times. That is item 8's
account; it is a different defect from this one and needs its own change
(`resid_at_own_composition`, default off, recommended on).

**Disposition.** Review section 5: "Ordering confirmed; effect unmeasured. Make
residual evaluation a deterministic composition fixed point and fill the
boundary from that same state. Test repeated `F(Y)` evaluations and Jacobian
directional derivatives." **Our judgment: agree**, and both halves are now done
-- the deterministic composition fixed point by sections 146/147, the repeated
`F(Y)` test by `src/tests/residual_determinism/`. (The Jacobian directional
derivative is unchanged by either: `||J r - dFD||_inf / ||J r||_inf = 1.0261e-01`
in both builds.) The review adds one implementation constraint worth
carrying: do not convert every interior cell to primitive and back merely to
fill the ghosts -- preserve the interior bytes and write only the ghost indices,
which removes the non-idempotent round trip without making last-bit
perturbations part of the boundary operation.

---

## 11. Migration and compatibility: restart state transport

**Reclassified after review.** This is not an open numerical defect. The three
repairs are verified in the tree, and what remains is a compatibility policy for
files written before them and a small diagnostic-consistency issue. It is kept
in the catalog because a reader restarting an old file needs to know what the
code will and will not infer.

### Symptom (as it was)

Three separate ways a restart did not reproduce the state it was restarting.

(a) The staged secondary-ionization coupling was not carried. WASP-121b's input
sets no `Secondary_ionization:` key, so the coupling is STAGED; the run that
wrote the state armed it at step 2242 and finished with it on, and every restart
re-staged it. Every marching experiment on WASP-121b in that series was run at a
photoionization rate its root was not converged under.

(b) A restart did not start from its own equilibrium composition. The first
step's `p = (n_tot + n_e) T` rewrite is a rate on any step of a relaxed run and
a **projection of fixed size** on the first step of a restart, because the
composition it corrects came from the file.

(c) The `T` column written from the steady solve belonged to a composition one
sweep behind the species columns beside it.

### Where the repairs are, verified in the tree

* `src/modules/functions/utilities.f90:39-79`, `write_coupling_state_header`,
  emits one `#` line into `Hydro_ioniz.txt`, which is the file a restart is fed
  as `Hydro_ioniz_IC.txt`. Confirmed present in written output: the header of
  `backup/regression/mol_sec_ion/output/Hydro_ioniz.txt` line 4 reads
  `# coupling: sec_ion=T sec_ion_step=0 valve=-1.00000E+00 recon=WENO3
  fluxconst=-1.00000E+00`.
* `src/modules/files_IO/write_output.f90:116` calls it.
* `src/modules/files_IO/load_IC.f90:54, 114-120, 585-623`,
  `ic_coupling_present` and `parse_coupling_header`.
* `src/EXHALE_main.f90:605-625`, the arming rule:

  ```fortran
  sec_ion_active = (use_sec_ion .and. (sec_ion_immediate .or. do_only_pp   &
                    .or. (do_load_IC .and. ic_coupling_present            &
                          .and. ic_sec_ion_active)))
  ```

* `src/EXHALE_main.f90:1741-1842`, `equilibrate_loaded_composition`, called at
  lines 358 and 679; opt-out `EXHALE_RELOAD_EQ=0` at line 346.
* `src/EXHALE_main.f90:1993`, the second `comp_T_from_p` after the post-solve
  sweep, with the reason written above it.
* `src/modules/files_IO/write_setup_report.f90:186-` reports which of the three
  restart cases the run is in.

### Measured size of what was fixed

(a) 3000 steps from the WASP-121b solution, peak-to-peak of `Delta(F r^2)/F_0`
(`Update_EXHALE_stage1.md:22062-22072`): at `r = 1.03`, 4.2961 before, 4.3197 with the
composition equilibration alone, and **0.0174** with the coupling header;
`du` ends at 5.187e-05 against the 5.03e-05 the accepted state carries.

(b) The first step's `dE_chem/E` at four CFLs, before and after: WASP-121b at
`r = 1.10` reads -6.2295e-03, -6.2297e-03, -6.2297e-03, -6.2297e-03 before (four
decades of `dt`, four identical values, i.e. a projection) and -3.2041e-06,
-3.2067e-07, -3.2070e-08, -3.2072e-09 after (one factor of ten per decade, i.e.
a rate). The first-step magnitude falls by 1944 on WASP-121b and by 93 on the hot
Uranus. Every named Load-IC case starts more than four parts in a thousand from
the equilibrium of its own state, and `armD_D2_newton` and `armD_D2_LW_newton`
start 71 and 50 percent away. **The bytes move; the answer does not**:
`armD_D2_LW_newton` gives `log10 Mdot` 8.88 either way, `armD_D2_newton` is
1.2 percent apart at a matched step 40,000 and `jfnk_hd189` 0.34 percent. Writing
more digits would not have helped: the files are already full double precision
and the gap is four to five orders above anything the format could cause.

(c) Only the `T` column moves, by at most 1.770e-05 (at `r = 1.0488`);
`r`, `rho`, `v`, `p`, `heat`, `cool` and all of `Ion_species.txt` are bit
identical and `log10 Mdot` is 10.34 either way.

### Status

**Fixed for files written from now on; restart files written before 2026-09-03
carry no `# coupling:` header and still start staged.** That is deliberate and is
not a formality: on `jfnk_hd189`, whose IC predates the line, arming it blindly
makes the composition mismatch worse (4.12e-03 -> 6.26e-03) and produces a
`du = 1.0e+02` transient by step 5500, so nothing can be inferred about a file
that does not carry the line. `Update_EXHALE_stage1.md` section 144, build `ec1f367c`.
Regression: 10/10 data identical (the header is a comment and the harness
compares with `grep -v '^ *#'`).

**What is left, and it is compatibility and diagnostics, not numerics.**

1. **Legacy files are permanently ambiguous.** A file with no `# coupling:`
   header starts staged, and nothing can be inferred about it. This is the right
   default and it does not expire.
2. **`heat` and `cool` lag the corrected `T`** in a file written from the steady
   solve, by one sweep. They are diagnostic columns. Either recompute every
   diagnostic column from the written state, or label them stale in the file's
   own header -- the review's recommendation, and the second is cheap.
3. **No restart schema version exists.** The `# coupling:` line is a census of
   run state, not a version. A schema version would let a reader -- and
   `load_IC` -- know which repairs a file was written under, instead of
   inferring it from the presence of one comment line.
4. **No executable round-trip test exists**, and the case that was supposed to
   be one cannot run (item 12.4). The review's recommendation is to rebuild the
   fixture from a deliberately short deterministic run rather than from an
   expensive converged product, and to test write/read/write identity and
   one-step continuity. That is now a decided action; see item 12.4.
5. **Not chased:** the reloaded WASP-121b residual improves 14 to 17 times in
   the hydro rows and its energy row still reads 1.0e-05 against the 3.138e-06
   the solve accepted. And the WASP-121b run with the arming was not carried
   past 3000 steps, so "fixed point" there means "did not leave in 3000 steps at
   `du = 5.2e-05`".

**Disposition.** Review section 5: "Main repairs confirmed in code. Legacy files
remain ambiguous and `heat/cool` can lag the corrected output state. Add a
restart schema version and an executable round-trip test; either recompute all
diagnostic columns from the written state or label them stale." **Our judgment:
agree**, and reclassified from open defect to migration and compatibility as
review section 3.5 recommends.

---

## 12. Smaller open items

> **Items 12.1 to 12.3 are superseded by section 151 (staged, not applied),
> which does what the review's Phase E1 asked: four explicit, mutually exclusive
> final states with one cross section each, opacity from their sum, and every
> species, electron and heat source built from a stoichiometric table. The four
> cross sections sum to `sigma_H2` to 4.4e-16 at every energy. Both new channels
> **default ON**, because both reactions are measured and the code was making
> fragments the events do not make; `off`/`False` exist only to reproduce a
> pre-151 result. Measured on `mol_base_handoff`: `log10 Mdot` 10.43 in all five
> configurations, the proton source moves 0.26 percent (against the "at most
> 0.5 percent" bound 12.2 derived), everything is below 1.10 `r_base`, and
> `n(H2+)` moves up to 2.0e-2. Across the matrix `n(H2+)` moves up to 24 percent
> and `n(H2)` 2 to 13 percent; **the two atomic cases are bitwise unchanged and
> the molecular goldens move.** The accounts below are kept as the record of the
> defects. Only 12.1 survives section 151 untouched.**

### 12.1 The H2 photoionization branching above 124 eV is an extrapolation (P31a) -- still open after section 151

Chung, Lee, Masuoka & Samson (1993), J. Chem. Phys. 99, 885, Table II -- the
source of `frac_H2_dissociative_ionization` -- stops at 124 eV, and the code
holds `f_di` at its value there, 0.2048, over the whole X-ray grid. The curve is
flat to 8 percent over its last 50 eV (0.2219 at 74 eV to 0.2048 at 124 eV), so
this is mild. A measurement or calculation above 124 eV would remove it. Note the
related reading error that was caught and fixed: Yan, Sadeghpour & Dalgarno
(1998) section 4 calls Chung et al.'s H+/H2+ column "the ratio to the total
photoionization cross section", and it is not; building on that sentence
overstates the channel by 1.26-1.28 above 70 eV, which the code measured as a
2.1 percent error in `Mdot` and a spurious 0.004 R_p outward shift of the H2
front before it was corrected (`Update_EXHALE_stage1.md` section 130).

### 12.2 Double photoionization of H2 is inside `f_di` instead of beside it (P31b)

`H2 + hv -> H+ + H+ + 2e-` has a vertical threshold of 51.4 eV and by 80 eV
about 20 percent of Chung et al.'s `sigma(H+)` comes from it. Their `sigma(H+)`
is a proton count and a double ionization yields two, so the code's single
dissociative branch stands in for two different reactions. Written out with
`N_d/N_i = 0.0225` above 80 eV: protons -2.2 percent, H2+ -2.2 percent, free
electrons -2.2 percent, and **H atoms +22 percent**, since the code makes one H
per dissociative event and the double channel makes none. Below 51.4 eV all four
errors are identically zero. Size bound: the flux-weighted `f_di` in the shielded
molecular base is 0.0467 against 0.22 above 80 eV, so at most 21 percent of the
dissociative rate there comes from the hard band and the proton-source error is
at most 0.5 percent. It is the H-atom co-product, not the proton count, that is
wrong at a level worth removing. Closing it needs `f_dd(E)` from Yan et al.
section 4 and an H+ row term of its own, since the H-nucleus closure cannot
supply a fragment that is not there.

### 12.3 Over 33-41 eV the H2 cross section is the absorption cross section, so the neutral-dissociation share is handed to H2+ (P31c)

A property of the Yan, Sadeghpour & Dalgarno (1998) fit, which follows Samson &
Haddad, whose photoionization yield drops below unity exactly in that window.
Against Chung et al.'s Table I the neutral share `sigma_n/sigma(abs)` runs
0.8 percent (33 eV), 4.0 (35), 6.6 (36.5), **7.4 (37.5, the peak)**, 6.0 (39),
0.8 (41), and is zero outside 33-41 eV. The H2 destruction rate is right --
absorption destroys the molecule in all three channels -- but up to 7.4 percent
of those events make two H atoms and no ion, and the code gives them an H2+ and
an electron. The resonance that makes this window matter is the same Q1/Q2
Rydberg structure that gives `f_di` its 35 eV peak, so the two corrections sit on
top of each other.

### 12.3b A separate finding of the same work: `sigma_H2` itself was a fit, and the fit was wrong where the wind absorbs

Section 151 replaced the three-piece Yan, Sadeghpour & Dalgarno (1998) fit with
**measurements** over the whole range the wind absorbs in: Backx, Wight & Van der
Wiel (1976), J. Phys. B 9, 315, Table 1 over **15.4-18 eV** (7 rows), and Samson
& Haddad (1994) Table 1 over 18-300 eV (72 rows); only the `E^-7/2` Eq. 19 tail
above 300 eV survives, rescaled.

**What was wrong.** The fit stepped by a factor **1.719 at 18 eV** and 0.619 at
85 eV -- a 72 percent discontinuity inside the band. Worse, **Eq. 17 fell to 0.17
at the 15.4 eV ionization threshold where the measurement is 14.2, a factor 85**:
the analytic fit went to zero at the ionization threshold, which no
photoionization cross section does. The error ran from 85x at threshold to 1.7x
at 18 eV. The new table steps from 0 to 14.20 at 15.4 eV, as an ionization
threshold must. The largest fractional step over a 1e-3 eV scan from 15.5 to
2000 eV is now `1.7e-4`, at 17.501 eV.

**How much any of it is known.** Three independent experiments were compared:

| E [eV] | 18 | 25 | 30 | 33 | 38 | 50 | 60 | 66 |
|---|---|---|---|---|---|---|---|---|
| Backx / Samson & Haddad | 0.984 | 0.998 | 0.996 | 1.080 | 1.210 | 1.104 | 1.005 | 1.207 |
| Lee / Samson & Haddad | 0.921 | 0.932 | 0.935 | 1.036 | 1.045 | 1.080 | 1.234 | 1.431 |

**Below about 31 eV the three agree to a few percent; above 32 eV they diverge by
10-20 percent, and by more at the top** (Lee runs 23 and 43 percent high at 60
and 66 eV). So `sigma_H2` is good to a few percent below 30 eV and to 10-20
percent above it, **whatever any single paper's quoted accuracy says** -- Samson
& Haddad's stated +/-2 to +/-3 percent is the accuracy of one experiment, not of
the quantity. An independent cross-check falls out of the same data: Backx's own
`f(H2+)/f_abs` reproduces Chung's `1 - f_di` to 0.1-0.5 percent over 20-30 eV.

**What it moves.** `sigma_H2` rises 204 percent at 16 eV and 74-79 percent over
17-18 eV, sits within 5 percent over most of 18-85 eV, falls 37.9 percent in the
last few eV below 85, and falls 11.4 percent at 300 eV. **Staged, not applied**;
the molecular goldens move and the atomic ones do not.

**One tension recorded, not resolved:** Samson & Haddad raised their own 300 eV
point from `1.54e-21` to `1.75e-21` cm^2 "to achieve agreement with the sum
rules", and the rescaled tail gives up 12 percent of that correction.

### 12.4 `backup/regression/roundtrip` cannot run: its restart IC is empty (AF) -- fixture rebuilt (section 149, staged)

> **Rebuilt 2026-09-03 (section 149, staged, not applied), the way the review
> asked**: from a deliberately short deterministic run -- 40 steps of the
> hot-Uranus molecular gate -- and not from an expensive converged product.
> Stage A is an ordinary named case; stages B and C reload with
> `EXHALE_DUMP_IC=1` and compare every species column **by label** to `1e-12`,
> plus `rho`/`v`/`p` and the `# coupling:` header. It **stays out of
> `DEFAULT_CASES`**, because its second half is not a bitwise comparison and
> needs its own invocation, and the `run_check.sh` header says so. Exercised on a
> replica: `check roundtrip` reports PASS against its own golden and
> `roundtrip_check.sh` passes; run on a binary without the patch the same case
> fails on all four molecular columns and on the coupling header, so it guards
> exactly what it is for.
>
> **Two defects were found by building it, both in the diagnostic and not in
> `load_IC`:** `EXHALE_DUMP_IC` ran before the first sweep, so it reported a
> molecular restart as atomic -- all four molecular columns exactly zero against
> the `5.34e12 cm^-3` of H2 the file carried -- and it ran ~280 lines before the
> rule that adopts the file's `# coupling:` state, so it wrote
> `sec_ion=F sec_ion_step=-1` over a file that stated `sec_ion=T
> sec_ion_step=40`.

`backup/regression/roundtrip/output/Hydro_ioniz_IC.txt` has zero data rows and
there is no `Ion_species_IC.txt` at all. The case sets `Load IC? True`, so it
aborts in `load_IC` before anything else runs. It is not in `run_check.sh`'s
`DEFAULT_CASES`, so the matrix never notices -- which is worth changing on its
own merits: its name says it is the restart round-trip guard, and the restart
path is exactly what item 11 changed, while all ten default cases are cold
starts and none exercises any of it. Either regenerate the pair of IC files from
a state the code wrote and add the case to the matrix, or retire it.

### 12.5 A stored arm case is refused at startup by the base-composition check (AE)

**QUARANTINED 2026-09-03 (user decision):** moved to
`backup/regression/_quarantined/armD_D1`, inputs untouched so the refusal stays
reproducible, with a `NOTE.txt` recording why and preserving the note the
directory already carried. `run_check.sh`'s header now says `_quarantined/` is
not a case directory. Paths in the rest of this section are the pre-move ones.

`backup/regression/armD_D1` carries `Molecular chemistry: False` together with
`Molecular base: True` and a `base.inp` with `q_H2_base 0.326896922`. The P35
refusal in `src/modules/files_IO/input_read.f90` stops it with `ERROR STOP 1`,
and **the refusal is correct**: the ghost the case would have run with sits at
859 K where its own `input.inp` asks for 1140 K. The case is stale, not the
check. It was the only one of the 65 case directories then present with this
combination, and it is not in the default matrix. Either it is molecular and `Molecular chemistry: True`
belongs in its input -- which changes what it tests and makes its stored outputs
meaningless -- or it is atomic and `Molecular base: True` plus `q_H2_base` come
out. Nothing in the directory records which. The numbers it was cited for were
produced by the pre-P35 binary at a base 25 percent colder than its input asked
for and should be treated as superseded either way.

### 12.5b Three regression case directories held the same cases twice (P53) -- resolved

`crit_warm` and `wasp_hybrid_finish`, `ptc_wasp` and `resid_on`,
`solver_newton_cold` and `wasp_full_newton` each held byte-identical `input.inp`
and `metals.inp` under two names, differing only in the stored outputs and logs
because the two members of each pair were run at different times.
**The first of each pair was deleted on 2026-09-03** (user decision), after the
md5 equality was re-checked immediately beforehand. None of the three was ever
in `DEFAULT_CASES` or in any `golden*/` directory -- `golden/` holds exactly the
ten default cases, and the dated `golden_block*/` snapshots hold subsets of them
-- so no golden entry had to be removed. Each surviving directory's `NOTE.txt`
names its deleted duplicate and repeats the md5, and every document that quoted
a measurement under a deleted name now says which case it is. `TO_BE_DONE.md`
(P53) carries the record; its remaining open half is
`backup/regression/valve_sens/eps4`, which holds two output files and no
`input.inp`, so no harness can run it.

### 12.6 No regression case takes the chemical-equilibrium branch of the base H2 fraction (P25)

`h2_mixing_ratio_base` has two sources -- a lower-atmosphere handoff, and the
Visscher/Koskinen fit at `(p_base, T0)` when no handoff exists -- and the matrix
exercises only the first. The five molecular cases all carry `q_H2_base 0.75`;
`lower_profile` takes its value from the profile and has the molecular base off
anyway; the two `wasp_*` cases have no molecular base. `grep -l "chem.eq. fit"
backup/regression/*/run.log` matches nothing. The branch is not dead --
`examples/18_oxygen_chemistry` and `examples/19_molecular_ir_bands` both run it
-- it is simply outside the matrix, and since the handoff branch now imposes the
base species composition while the fit branch deliberately does not, the two are
different paths through `ionization_equilibrium` and only one is guarded.

### 12.7 The PLM to WENO3 switch changes the discrete operator in one step (S) -- continuation built, default off (section 148, staged)

The review reframed this from a `du` floor to what it is: switching the discrete
operator discontinuously and expecting the state to remain a solution. Section
148 implements the continuation it asked for,
`R_lambda = (1 - lambda) R_PLM + lambda R_WENO3`, behind
`Reconstruction continuation: <dlambda> [<dtu_tol>] [fixed|adaptive]`,
**default off** (`dlambda <= 0`, or the key absent, is the shipped one-step
hand-off; `1.0` reaches it through the continuation code and is byte-identical).

**What the operator difference actually is.** On `mol_sec_ion` at the hand-off
step, the two operators differ by a factor 6.74 of `R_PLM` in the mass row -- and
the whole of that difference is at the base: the `1.00-1.02 R_p` band carries
1.000 of the mass-row and energy-row difference. **In the wind between the base
and the outer edge the two operators agree to a part in 10^4 to 10^2.** The
switch costs a `dtu` jump of 1.763x, and a 500-step ramp reduces it to 0.990x.

| ramp | lambda = 1 at | back-offs | peak `dtu` ratio |
|---|---|---|---|
| shipped one-step switch | step 11861 | -- | **1.763** |
| `dlambda = 0.02`, fixed | 11910 | 0 | 1.248 |
| `dlambda = 0.02`, adaptive | 11919 | 2 | 1.200 |
| `dlambda = 0.002`, adaptive | 12360 | 0 | **0.990** |

`log10 Mdot` is 10.35 in every arm; the ramps cost no wall time (they are
slightly faster). A control that stays on PLM to the cap does **not** converge
(`du = 3.9e-3` at 15000 steps) and lands 0.08 dex higher with a base 400-700 K
hotter, so the switch is buying something real.

**Two corrections to the first edition's account.** The A2 molecular wind is not
what was measured -- the two cases that cross the hand-off are `mol_sec_ion` and
`wasp_full` -- and `wasp_full` crosses at step 1, where the hand-off is
measurably harmless. **Not tested:** the JFNK finish was never exercised through
a ramp, which is the form Phase F actually asks for; and the eight matrix cases
that never leave the PLM stage were not re-run. Byte-identical with the key
absent and at `1.0`; no golden moves.

### 12.8 The steady Picard loop double-counts the oxygen-carrier nuclei into the free-stage columns (T) -- NOT REPRODUCED

**This is the item the review raised above all the others in section 12 -- "fix
this before any coupled A2 result is accepted" -- and section 149 went looking
for it with the assertions the review specified. It did not find it.**

The element assertions were run through the A2 oxygen configuration on the
steady/Picard path item (T) names, including 20 forced passes of
`relax_photochemical_composition`, 42 steady-iterate sweeps and 1230 probe
sweeps, and **no element budget broke**: every `ioniz_eq` sweep, every carrier
write-back, every transport step and the restart equilibration closed inside
`1e-9`, and the finished output closes all five elements at `1.77e-10`. The
original report was carbon `+90.6%` and oxygen `+49.0%` near
`r = 1.003-1.006` -- eight orders larger. A dedicated test reproduces the exact
configuration (a transported CO collapsing by 2.5 decades) and closes at
`1.9e-16`. The code paths were followed and none creates a nucleus.

**This is "not reproduced on the equivalent path", not "proved fixed".** The
likeliest explanation offered is that the carrier rework of `Update_EXHALE_stage1.md`
section 128 (2026-09-02) removed it -- which postdates the item's opening
(2026-08-31). What could not be done is the decisive run: the original state is a
converged A2 wind reached from a ~25,000-step march plus a JFNK finish, it does
not exist anywhere in the tree, and rebuilding it is hours of wall time.

**Status: open.** Settling it needs one state and one run -- a state that reaches
the accepted A2 root, marched with the assertions armed. Until then the item
stays open and the numbers measured on the contaminated state keep their caveat.
The residue that *is* seen -- the steady path closing at `1.77e-10` where the
marching path closes at `2.6e-13`, identically for all five elements -- is the
hydrogen denominator, not the elements, and is reported rather than chased.

### 12.9 `Base ghost temperature: isothermal` does not settle to a single base state (D)

Open, `TO_BE_DONE.md` (D), opened 2026-08-12. Related to item 4(b).

### 12.10 The molecular layer is 100-450 K colder than Koskinen et al. (2022) (X)

Measured, cause not established. The obvious suspect was tested and is **not**
supported: removing the H3+ coolant entirely buys 56 K at 1.001 r_base, 11 K at
1.01, 3 K at 1e-2 and 3 K at 3e-2, against deficits of 359 K and 206 K at the
last two rows, and the temperature minimum does not move; our H3+ formula applied
to Koskinen's own digitized `T`, `n(H2)` and `n(H3+)` reproduces the radiative
cooling they plot to a factor 1.0-2.9 (median 2.0). The 60-100x between the two
codes' radiative cooling is `n(H3+)`, which follows `f(H2)`, which is item 1.
What the trough does sit against is the expansion term: over 1.02-1.05 r_base the
adiabatic `p div v` is 8.6e-8 to 9.1e-8 erg cm^-3 s^-1 against a photoheating of
4.9e-8 to 7.2e-8 and a radiative cooling falling from 1.5e-8 to 5e-10. The
caveat that limits this is item 6: the converged state's own energy budget does
not close below 1.08 r_base.

### 12.11 The paper and poster materials stand on a superseded binary (R)

Open, `TO_BE_DONE.md` (R), opened 2026-08-31, and the physics changes since --
sections 128 through 144 -- have all landed after the figures were made. The
review classifies this as a provenance issue rather than a code defect and asks
for one thing this catalog did not: **mark every stored product with the
executable or source revision and its physics-option metadata**, so that a
figure can be told from its file which binary and which options produced it.
Stored outputs without those fields should not be used as regression references.

### Dispositions for section 12

| item | Review disposition (section 5) | Our judgment |
|---|---|---|
| 12.1-12.3 H2 photochannels | Approximations documented in code; channel stoichiometry incomplete for double ionization and neutral dissociation. Introduce explicit cross sections for mutually exclusive final states and derive every species and electron source from their stoichiometric vectors. | **Agree**, and this is a better design than the three separate patches the first edition proposed. See Phase E1. |
| 12.4 round-trip case | Confirmed. Rebuild the fixture from a deliberately short deterministic run, not from an expensive converged product. Test write/read/write identity and one-step continuity. | **Agree**, and adopted -- the first edition proposed regenerating from a converged run, which is both expensive and the wrong kind of fixture. |
| 12.5 stale arm case | Confirmed configuration conflict; intent cannot be inferred from code. Quarantine it from active regression and label the historical result invalid. | **Agree**, and executed 2026-09-03: the case is moved to `backup/regression/_quarantined/armD_D1` with a `NOTE.txt` recording why. |
| 12.6 base-H2 fit coverage | Confirmed by path analysis; not rerun. Add a small setup/equilibrium test that omits the handoff and asserts the fitted base H2 fraction and the resulting ghost composition. | **Agree.** A setup-level test is cheaper than the run-level case the first edition implied. |
| 12.7 PLM to WENO3 | Open path verified; diagnosis historical. Treat PLM and WENO3 as separate discretizations: use PLM to form an initial guess, then continue in the residual or flux blend before solving the WENO3 system. | **Agree**, and it reframes the item -- the first edition recorded it as a `du` floor, which is a symptom of switching the discrete operator in one step. See Phase F. |
| 12.8 oxygen-carrier double count | Critical reported conservation failure; the exact offending refresh was not isolated. Add element-budget assertions around every ionization sweep, carrier write-back, and outer steady pass. Fix before any coupled A2 result is accepted. | **Agree, and this item is under-ranked in the first edition.** It is a prerequisite for item 2, not a smaller item. See Phase A1. |
| 12.9 isothermal base state | Reported measurement, closely coupled to item 4. Resolve through the characteristic boundary work; do not add another empirical ghost-temperature mode first. | **Agree**, and it supersedes item 4's remedy option 5 as a standalone change. |
| 12.10 cold molecular layer | Symptom, not an independent defect yet. Reassess only after transported H2 and base energy closure are obtained. | **Agree.** |
| 12.11 stale products | Provenance issue, not a code defect. Mark products with executable/source revision and physics-option metadata; do not update scientific figures until the critical physics path is accepted. | **Agree.** |

---

## Program

Review section 6 proposes an implementation program in six phases. It is
reproduced here in outline with the status of each, because most of it now
exists in patch or in scratch builds rather than in the tree, and a reader of
this catalog needs to know which.

**Three status words are used below and they mean different things.**

* **landed** -- in the working tree, verified there.
* **staged for block I** -- written, built, measured, dry-run clean, and **not
  applied**. Every such section's own document says so in its own words.
* **in progress / not started** -- as stated.

**Phase status at a glance.**

| phase | status |
|---|---|
| A1 elemental accounting | **staged for block I** (149), default off; 12.8 not reproduced |
| A2 residual reporting | **landed** (145.3-145.4), as a single cellwise gate with both numbers reported |
| A3 deterministic tests | **partly staged** (147, 149): two of five exist |
| B1 ghost/composition lag | **staged for block I** (146 + 147) |
| B2 update-map diagnostic | **landed** (145.6) -- and it decided item 3 |
| C characteristic boundary | **staged for block I** (152); D1-D5 answered, D6 opened |
| D1 decouple transport from oxygen | **not started** |
| D2 trust region | **staged for block I** (153/154), default off, does not converge |
| D3 replace the references | **not started**, correctly blocked on D2 |
| E1 H2 photochannels | **staged for block I** (151), both channels default on |
| E2 Lyman-Werner geometry | **staged for block I** (150); bound obtained, table deliberately not rebuilt |
| F continuation and stop naming | **partly landed** (145.5); continuation **staged** (148), default off; the three names not given |
| (155) self-consistent acceptance residual | **staged for block I**, default off with a recommendation to turn it on; plus `rho` as an input to the sweep |
| (158) carrier row transported in the conserved variable | **LANDED 2026-09-04.** The premise of item FC(2) did not survive: the two paths already run one operator, and the flat `|G+R|` was the diagnostic measuring `G` in the carrier density, so the hydro stage's term arrived as a constant. The real defect was the ADVECTED VARIABLE -- a particle mixing ratio against the non-conserved `n_tot`, leaving a spurious `n_c v d ln(mbar)/dr` across both fronts. The row now rides the fraction per unit mass; diffusion stays on the particle mixing ratio. Base cell first order after, no trend before. Regression: only `mol_carrier` moves (log Mdot 10.50 -> 10.52). Two side defects fixed with it: the carrier operator was overwriting the LOWER ghost, degenerating the base Dirichlet into zero-gradient from the second relaxation pass on, and `EXHALE_UPDATE_MAP` indexed its step list by the marching counter, reading out of range on a cold start (`dt = 0`, NaN). |

**One ordering fact that matters more than any single status.** The review's
Phase A was meant to come first -- "change diagnostics and tests before it
changes physical solutions". What happened is that A2 and B2 landed first and
**B2 settled item 3**, while A1 landed nowhere and its target, 12.8, could not be
reproduced. The review's ordering held up: the two diagnostics that landed are
the two that decided things.

The review's own ordering argument is worth keeping in front: **items 1 to 6 of
this catalog are not independent.** Items 3, 5 and 6 are different observations
of the base-layer state, and item 12.8 blocks the physically correct resolution
of item 1. Solving them as unrelated patches would make it hard to establish
which equation the final state satisfies.

### Phase A -- establish non-negotiable invariants (diagnostics and tests before solutions)

| step | what it is | status |
|---|---|---|
| **A1** centralize elemental accounting -- **STAGED FOR BLOCK I** | one routine computing H, He, C, O and metal nuclei from `rho` and `f_sp` with exact stoichiometric multiplicities for H2, the molecular ions, HeH+, OH, H2O and CO; debug assertions before and after `carrier_write_back`, `ioniz_eq` in the steady outer loop, `relax_photochemical_composition`, restart equilibration, and acceptance/output. Tolerance from round-off and the cell solver's own tolerance, not a percentage of abundance. On failure report element, cell, radius, before/after totals and each contributing species. | **Staged, section 149.** New module `src/modules/functions/element_census.f90`, reading every multiplicity from `species_table` and writing none of its own. Seven assertion sites, exactly the ones the review listed. **Default off** (`EXHALE_ELEMENT_ASSERT=1` reports, `=2` stops); cost within noise. Tolerance `1e-9`, sitting between the measured accumulated round-off (`3.5e-16` per operator) and the cell solver's own `xtol = 1.5e-8` -- not a percentage of abundance, as asked. Byte-identical on three cases at 300 steps. **Its target was not reproduced: item 12.8.** |
| **A2** dual residual reporting -- **LANDED** | keep the integrated norm as a budget measure, add `max_j \|R_j\|/S_j` for each row and region, report dimensional terms at each worst cell; diagnostic first, then part of acceptance after calibration. **Do not** use `max(\|R\|)/max(S)` as a replacement -- it is not a local relative norm. | **LANDED, diverging in one respect and satisfying the review in another.** 145.3 makes the cellwise maximum the **single** acceptance form and removes `resid_vol`, `residual_norms_vol` and the `Resid norm` key (the key still parses and warns it is retired). The review asked for both norms as gates with separate tolerances; that half is declined and item 8 records the trade. But 145.4 grants the concrete request: `write_residual_breakdown` prints, wherever a residual is reported, each row's cellwise maximum with its cell index and radius, the integrated ratio with numerator and denominator, the contributions of the four radial bands, and for the worst cell the signed dimensional residual and every term its row differences. **Both numbers are visible everywhere; only one gates.** |
| **A3** focused deterministic tests -- **PARTLY STAGED** | elemental closure through one carrier write-back and one ionization sweep; repeated residual evaluation requiring `F(Y)` to be independent of the preceding trial state; an H2 photoevent ledger conserving H nuclei and charge in every channel; restart write/read/write preservation; the chemical-equilibrium base-H2 branch without a handoff. Small module or short-run tests, not regenerated benchmarks. | **Two of five staged, three not started.** `src/tests/residual_determinism/` (147) is the `F(Y)` test and it now passes structurally rather than by discipline. `src/tests/element_census_tests.f90` (149) is 20 tests, 20 passing, covering elemental closure through a carrier write-back, the H2 photoevent ledger (H nuclei and charge to `1e-14` over 15-200 eV) and the chemical-equilibrium base-H2 branch -- so item 12.6's test exists, at setup level. **Not started:** restart write/read/write as a unit test; the round-trip case (12.4) covers it as a regression case instead. |

### Phase B -- make the residual and the production map consistent

| step | what it is | status |
|---|---|---|
| **B1** remove the ghost and composition lag -- **STAGED FOR BLOCK I** | refactor `eval_residual` into a state function: unpack unknowns; solve the composition and temperature fixed point for the physical cells; construct boundary face and ghost states from the converged composition; evaluate columns and radiation consistently with that boundary state; assemble the hydrodynamic and carrier residuals from the same state. Preserve the interior bytes and write only the ghost indices. | **Staged, sections 146, 147 and 155; this is item 10.** Section 155 adds the half the other two do not cover: the residual the gate REPORTS is evaluated at the previous iterate's composition. `resid_at_own_composition` re-evaluates at hand-back and iterates the sweep at fixed `Y` to its fixed point (3 to 9 sweeps, one-off), default off because turning it on makes every `Resid tol` root report 4 to 56 times more residual -- a decision about what "converged" means. The review's acceptance condition -- "repeating `eval_residual` at the same `Y` produces the same residual regardless of the prior trial history" -- is met, and it took the seed interface to meet it: the ordering fix alone left 1412 of 1500 entries differing. The review's implementation constraint was followed exactly: `Apply_BC` now writes only the ghost indices instead of round-tripping every interior cell. 146 needs a golden refresh; 147 is byte-identical. |
| **B2** the complete update-map diagnostic -- **LANDED** | expose `G_dt(q) = [Phi_dt(q) - q]/dt` from one production step without stop logic or file output; compare with the steady residual at the same state for at least three decreasing `dt`, broken down after the Euler, chemistry, energy, carrier, diffusion, boundary and filter stages. | **LANDED (145.6), and it decided item 3.** `EXHALE_UPDATE_MAP` uses the marching loop's own body as `Phi_dt`, not a copy. `G_dt -> -R` at first order in `dt` on both planets, so the two maps share a fixed point and the review's FIRST outcome is the one that obtained. The stage breakdown is there too. The review's third outcome is realized and is now provable rather than suspected: with the Shapiro filter on, the filter's increment is applied per step, so as a rate it grows as `1/dt` and the two cannot share a fixed point. Two implementation traps are recorded in the code, and a third cost the regression outright -- item 12.12. |

### Phase C -- redesign the lower boundary at the face

Implement the lower condition at `r_edg(0)` instead of building a face state by
calling a ghost-cell routine: take the outgoing acoustic invariant from the first
interior reconstructed state, prescribe two independent reservoir quantities with
explicit physical meaning, solve for the boundary-face primitive state, and
construct ghost averages reproducing it to the order of the reconstruction.
Prescribe all incoming characteristics for supersonic inflow, and choose the
characteristic count explicitly for outflow or reversal rather than passing the
interior velocity through a softplus expression.

**Status: STAGED FOR BLOCK I (section 152), not applied.** The design decisions
D1-D5 the first edition said should be agreed first were put to the user and
answered on 2026-09-03; the reservoir pair is `(p, s)`. The probe and production
measurements are in item 4: the base sawtooth is gone rather than reduced, the
restart reproducibility defect is repaired, and two negative results stand --
the item 3 marching departure gets **worse** (2.09x at `r = 1.005`), and the
acoustic-pulse validation **fails and is shown to be unreachable by any static
reservoir pair**, which is the new item D6. Of the review's five validations,
the hydrostatic and manufactured-inflow tests pass, the acoustic pulse fails for
a reason now understood, the reload-and-march comparison was made, and the
**three-base-grid test was not run**. Goldens move for every case; the refresh
belongs at the end of the series. `docs/EXHALE_BC_and_IC.tex` is the
existing boundary and initial-condition reference and does not describe a
characteristic condition. This is item 4's remedy option 1 and item 12.9's
resolution, and it is **a public change to the mathematical boundary-value
problem** -- the reservoir pair and the flow-reversal behavior should be agreed
before implementation. Validation is the five-test list in item 4's option 1.

### Phase D -- finish transported molecular steady states

* **D1 decouple transport from oxygen.** Replace
  `if (.not. carrier_transport_stated) carrier_transport = thereis_oxychem`
  (`input_read.f90:1662`) with one molecular rule: transport is irrelevant for an
  atomic run; a molecular run with transport explicitly false continues only with
  a prominent invalid-closure warning reporting the omitted transport/chemistry
  ratio; a production molecular result requires transport and a carrier residual
  gate. Oxygen chemistry adds the OH, H2O and CO carriers but must not decide
  whether H2 is transported. **Status: NOT STARTED**, and it is now the cheapest
  unstarted item in the program. This is item 1's remedy 3 plus its remedy 2,
  and the review makes them one step.
* **D2 a real trust region in the coupled solver.** Item 2's remedy 1.
  **Status: STAGED FOR BLOCK I (sections 153/154), default off, and it does not
  converge.** Everything the review specified was built -- a dogleg step,
  predicted reduction from the same scaled model, radius updated by the
  actual/predicted ratio, carrier headroom and admissibility enforced inside the
  step, and a finite-difference ray test that refuses a model whose directional
  derivative disagrees. It works as a trust region: reduction ratio 0.998 to
  1.000 on the first iterations, median 0.933, **zero models refused by the ray
  test**, element-budget refusals 18-19 down to **0**. And the D2 success list is
  missed on every physical line. **The obstruction is the base cell**, whose row
  no step can balance because the state there is imposed rather than solved --
  which makes D2 downstream of Phase C, an ordering the review did not
  anticipate. `p51b_dtau_rule.patch` remains the rejected alternative.
* **D3 replace the local-equilibrium molecular references.** Only after D2. The
  acceptance comparison is physical and numerical -- front location, H3+ column,
  temperature, mass-loss rate, energy closure, grid convergence, restart
  stability. **Byte identity with the old local-equilibrium results is not an
  acceptance criterion for this physics correction.**

### Phase E -- independent molecular microphysics

* **E1 H2 photoionization final states. Status: STAGED FOR BLOCK I (section
  151).** Done as specified -- four mutually exclusive channels, one cross
  section each, opacity from their sum, sources from a stoichiometric table,
  closing to `4.4e-16` at every energy. Both channels **default on**. Two things
  the review did not ask for came out of it: `sigma_H2` itself was replaced by
  measurements (item 12.3b), and above 124 eV the branching is still held at its
  124 eV value, which is item 12.1 and is unchanged. The specification was:
  represent the absorption cross section
  as a sum of explicit, mutually exclusive channels -- `H2 + hv -> H2+ + e-`,
  `-> H + H+ + e-`, `-> H+ + H+ + 2e-`, `-> H + H` -- with one cross section
  each, opacity from their sum, and species, electron and heat sources from a
  stoichiometric table. This prevents double ionization and neutral dissociation
  from being hidden inside `f_di` and subsumes items 12.1 to 12.3. Above 124 eV,
  either obtain a defensible model or expose the constant branching as an
  uncertainty option; do not present it as measured. **Status: not started.**
* **E2 Lyman-Werner geometry.** Item 9. **Status: STAGED FOR BLOCK I (section
  150).** The matched slab/spherical comparison was done and the slab limit
  validated to seven significant digits first, as asked. The outcome is a
  one-signed bound (the slab over-traps; 6.5 percent rate-weighted) and a
  decision **not** to rebuild the table, because two independent trapping
  calculations differ from each other by more than the geometry does. The
  provenance and uncertainty the review asked for are what went in. A separate
  defect found in the same work -- three dayside beam conventions in one run --
  is fixed in the same patch and does move `mol_lyman_werner`.

### Phase F -- convergence workflow and reconstruction

Keep marching stop, solver handoff and scientific acceptance separate in name and
in calibration, and never label a `du`-stopped run a steady solution -- item 7.
For the PLM-to-WENO3 cases, do not switch the discrete operator discontinuously
and expect the state to remain a solution; use PLM to form an initial guess and
then continue along `R_lambda = (1 - lambda) R_PLM + lambda R_WENO3` with step
control in `lambda`, evaluating the final gates with the pure WENO3 operator.
This also separates a WENO3 instability from a transient caused only by changing
operators.

**Status: partly landed, partly staged.** Section 145.5 landed the half that
matters most: every `du`, `dtu`, plateau and step-cap stop now prints "THIS IS A
MARCHING STOP, NOT AN ACCEPTED STEADY STATE" together with all three gates
**measured on the state being stopped at**, and three branches that used to say
"converged" now say "stopped". Measured on a hot-Uranus run that used to print
"converged: mass flux constant" at `du = 1.454e-04`: its cellwise residual is
**1.282** in the energy row of cell 7. The unification is formally withdrawn,
with the 13-case census kept as the reason. The PLM-to-WENO3 continuation is
**staged** (section 148, default off) -- see item 12.7, including the caveat that
the continuation ending inside a JFNK solve, which is the form Phase F actually
asks for, is implemented but untested. **Still not done:** the three concepts do
not have three names, and their tolerances are not separately calibrated.

### Definition of completion, and what to say until then

Review section 8 sets the bar for the molecular defect group: one transported
molecular case and one oxygen-chemistry case must show elemental and charge
conservation through every operator; the H2 continuity equation passing
integrated and local gates; mass, momentum and energy passing both gates
throughout the domain with the boundary cells reported separately; the accepted
state a fixed point of the complete production update; a subsonic and
refinement-convergent lower boundary face state; front position, H3+ column,
temperature minimum and mass-loss rate converged over at least three grids or
with a quantified order and error; a restart returning to the same fixed point;
and independence from whether the initial guess came from PLM marching, a
previous WENO3 state, or a compatible restart.

Until those hold, the review's recommendation -- which this catalog endorses --
is that the code should **label default-off carrier transport and
local-equilibrium molecular outputs as experimental or closure-limited.** That is
preferable to treating solver non-convergence as permission to publish a solution
known not to satisfy the transported H2 continuity equation.

### 12.12 Instrumentation that does nothing at run time changed the answer, and the binary had stopped being a function of its source (section 145.8, landed)

Three findings from one regression failure, all worth carrying because they are
about how every other number in this catalog is verified.

**A disabled `if` moved the last bit.** The tree binary failed all ten regression
cases at a relative `3e-14` to `9e-11`. It was not a stale object --
`make distclean && make` reproduced it exactly, and a bisect over clean builds
left one patch. The update-map snapshots had been placed *between the statements
of one RK3 stage*: `a - b*c` is a fused-multiply-add candidate,
`-ffp-contract=fast` is gfortran's default, and inserting a statement changes
whether the FMA is taken -- one ulp per step, `1e-11` over eleven thousand.
**Every inserted statement was inside `if (upmap_n .gt. 0)` and did nothing with
the diagnostic off.** The contraction setting is deliberately not changed: it is
the default the goldens were made under and it is both faster and more accurate.
What changed is where instrumentation may go -- out of line, or between
operator-split stages, never between the statements of one.

**The binary was not a function of its source.** Section 145's first build stamp
compiled the *build time* into a constant, so two `make distclean && make` of
identical source produced different executables, differing in one object in
exactly that string, and three workers comparing md5s read "different program"
from the same program. The stamp now carries the git revision and dirty flag
only, and the run time comes from `date_and_time` in the provenance header.

**And nothing checked that the binary matched the sources.** The stamp rule had
a phony prerequisite, so `make` never reached a fixed point and `make -q`
reported "out of date" immediately after a successful `make` -- the one cheap
check that would catch a stale binary was useless. The Makefile now writes the
stamp at parse time and only when its content changes, and `run_check.sh` runs
`make -q` after its own `make` and refuses the matrix with exit code 2 if
anything is still out of date.

**What an md5 of the binary does and does not identify.** A comment does not
change the code, but a *shifted line number* does: every I/O statement carries
its own source line as an immediate for runtime error reporting, so adding one
comment line above a `write` turns `movl $0x4` into `movl $0x5`. An eight-line
comment moved the binary from `e83b4ab4` to `62cc6c64` with no numerical change,
both bitwise identical to the goldens. **The md5 is a fingerprint of one source
text, not a semantic fingerprint; two binaries with different md5s are not
thereby different programs.**

**Not established:** whether the build system could ever have missed a rebuild
before section 145. Reading it found no mechanism and a recurrence test gives
identical binaries, but no earlier build recorded what it was built from, so the
earlier builds of 2026-09-03 cannot be reconstructed either way.

---

## D6. The lower boundary reflects sound completely, and no static reservoir pair can stop it

### Symptom

An acoustic pulse leaving the domain through the lower boundary is reflected
almost entirely: `|R| = 0.956` under the characteristic closure of section 152,
against `0.953` under the legacy closure it replaces. The reflection is
independent of amplitude (`1e-3`, `1e-4`, `1e-5` all give 0.956) and falls only
slowly with `dr`.

### Where it arises

`characteristic_base_face_state` (section 152, staged). Not a coding defect and
not a tuning failure: the closure states two *static* reservoir quantities --
pressure and specific entropy at the base level -- and takes the third from the
outgoing invariant.

### Why it is wrong

If both reservoir quantities are static then the incoming acoustic amplitude is
pinned, so a wave arriving from inside is reflected. **That is true of `(p, s)`,
of `(rho, s)` and of `(rho, p)` alike: the pair does not matter, the staticness
does.** A reservoir that holds a static quantity is an acoustic node whatever
else it does.

### Measured size

`|R| = 0.956`, amplitude-independent, against a design target of `< 0.1`. The
one form that is transparent is rejected on the steady state: weighting the face
pressure toward the invariant gives `|R| = 0.0024` but drives the steady-state
test's `R_mom(1)/(rho g)` from `1.970e-4` to `1.323e-1` -- at base Mach `1e-3`
and `H/dr = 131` that is 460 times the interior truncation error. **Acoustic
transparency is bought by corrupting the steady state**, so the blend constant
stays at zero.

### Remedy options

The source document names three and takes none.

1. **Leave it.** The reflection is no worse than the closure being replaced.
2. **Add a relaxation condition to the MARCHING map only and forbid it during a
   steady solve**, exactly as the Shapiro filter is now forbidden by the update
   map's measurement (item 3). A relaxation condition is not a function of the
   instantaneous state, so it has no fixed point the steady residual can
   express -- which is precisely why it cannot be in both maps.
3. **Reformulate the steady residual to carry a boundary state of its own.**

### Status

**New, open, raised by the implementation rather than by the review.** It decides
whether item (AA)'s base-layer ringing is in scope at all. One prediction is
withdrawn because of it: Phase C was expected to damp that ringing, and it does
not -- item (AA)'s layer swings and the 2050-2300 step period of
`p54_base_layer_mass_flux.md` section 5 are not addressed by the new boundary.

**Naming caveat.** "D6" is the label the Phase C design document uses for this,
as the sixth of its design decisions. It collides twice: `TO_BE_DONE.md`'s letter
series already has an item **(D)**, and `docs/a2_oxygen_option_design.md` uses
"D6" for an unrelated base-condition arm. If this becomes a `TO_BE_DONE.md` item
it needs a name that is unique in that series -- the same discipline applied to
`P52` earlier in this catalog.

---

## FC. The H2 front creeps outward and never stops

### Symptom

Marched from a cold start with the current physics, the transported hot Uranus
converges on its own criterion at step 243,458 (`du = 9.9996e-4`) -- and the H2
front has not stopped moving at any point along the way. It creeps outward at
**1.0 to 1.5 cells per 20,000 steps** from the first snapshot to the last, and
gains `0.010 R_p` in the final 43,000 steps alone.

**The state the marching calls converged is not the flattest one it passed
through.** The gate window falls monotonically to a minimum of `1.502e-2` at
200,000 steps and then rises by a factor 4.7 to `7.037e-2` by the time `du`
crosses its threshold; the two inner windows turn at the same place and rise by
25 and 61 percent. The flattest state it reached is still three times the flux
gate.

### Where it arises

Not localized to a routine. The measurement is `p23_transport_on_state.md`
sections 8 and 9, on the tree's own binary with `Coupled carrier solve` off and
the stall stop disabled through the ordinary `Stall [tol,N]: 0.0 2000` key --
no modified source.

### What has been excluded

**The base cells are not what carries the front.** Three 60,000-step
continuations from the gate-window minimum, differing only in an environment
probe: a reference arm, one freezing the first five physical cells to their
step-1 state after every step, and one freezing the first two. The freeze
demonstrably works (the frozen cells move by `~7e-16` against `~1e-3` in the
reference). **Every front level is identical across the three arms at both
times**, to the four decimals the trace prints, and the front moves in all three
by the same amount; `du` agrees to three digits and all three converge at 42,458,
42,598 and 42,530 steps with `log10 Mdot = 10.45`. The two inner spread windows
differ between arms by 3 to 7 percent -- the freeze does change the layer, as it
must. **So the base cells move the residual and not the front.**

That separates this from item 4 and from item 2's floor: the residual floor at
the base cell and the front's motion are two different problems, and **the
characteristic boundary of section 152 will not fix this one.**

### What has been decided, in three sentences (section 10)

**(1) The creep is motion in physical time, not a step bias.** Three arms
differing only in `CFL:` -- 0.6, 0.3, 0.15, a factor 4 on `dt`, an input key and
no source change -- traced at step intervals in the same ratio so that matched
rows are matched physical times, put **all three front levels at exactly the
same radius at every matched time**:

| matched time | steps at CFL 0.6 / 0.3 / 0.15 | `f=0.5` | `f=1e-2` | `f=1e-4` |
|---:|---|---|---|---|
| 1 | 2,000 / 4,000 / 8,000 | 1.0566 / 1.0566 / 1.0566 | 1.1389 / 1.1389 / 1.1389 | 1.1518 / 1.1518 / 1.1518 |
| 3 | 6,000 / 12,000 / 24,000 | 1.0566 / 1.0566 / 1.0566 | 1.1407 / 1.1407 / 1.1407 | 1.1518 / 1.1518 / 1.1518 |
| 6 | 12,000 / 24,000 / -- | 1.0574 / 1.0574 | 1.1425 / 1.1425 | 1.1537 / 1.1537 |

A step-count bias would put the small-`dt` arm four times further per unit
physical time; it puts it in the same place, with `du` agreeing to two or three
digits. This also retires an earlier reading: section 134's grid pair, "the same
advance per 20,000 steps on 500 and 1000 cells", cannot be read as a step bias
without the physical times behind it.

**(2) A new defect, found on the way: the carrier marching step does not
converge to the carrier steady residual as `dt -> 0`.** Extending the section
145 update map to the carrier row -- the one row it did not carry, and the row
the front lives in -- and comparing
`G_c = [n(H2)_after - n(H2)_before]/dt` with `-R_carrier` over four decades of
`dt`, **`|G + R|` is `dt`-independent to four significant figures at every
radius**: `3.2125e2, 3.2191e2, 3.2194e2, 3.2195e2` at the front. A consistent
map gives `O(dt)`. As a fraction of the row's own terms it is worst at the
**base cells** (1.9e-2 at `r = 1.0002` and `1.0015`) and 3.4 per cent at the
front, falling to 1.7 per cent by 1.16 `R_p`. Three candidates are **excluded by
bit-identical probes** -- the deferred second-order advection correction
(`EXHALE_CARR_NOCORR`), the frozen limited slope (`EXHALE_CARR_SLOPELIVE`) and
the element-budget clamp (`EXHALE_CARR_NOLIMIT`), all three banners firing. The
remaining candidate is the **base inflow Dirichlet condition, handled differently
by the two paths**; it has no switch, giving it one is a change to the operator
rather than a probe of it, and it is the one consistent with where the
inconsistency is largest. A repair is in progress as section 157, which does not
exist yet.

> **SETTLED 2026-09-04, and the paragraph above is wrong in its diagnosis while
> right in its instinct** (`Update_EXHALE_stage1.md` section 158, numbered 158 because
> 157 became the coupled-solve re-test). The two paths do NOT handle the base
> inflow condition differently -- they run the same assembly sequence in the
> same order and `carrier_base_state` builds it once for both. Measured in the
> variable the carrier residual is written in, the marching step converges on
> that residual at clean first order (`9.5e-5, 7.2e-6, 8.2e-7, 8.3e-8` over
> four decades of `dt`), with the clamp and the write-back contributing
> nothing. **The flat `|G+R|` was the diagnostic's own variable**: it formed
> `G` in the carrier density `n(H2) = f_sp*rho*n0`, so the `dt`-independent
> part was the HYDRO stage's term, `+7.154e-1` against a transport term of
> `-5.470e-2` matching the residual `-5.474e-2` to three digits.
>
> The instinct -- that something upstream of the step control was wrong -- was
> correct, and the defect is a physics one: the carrier was advected as a
> PARTICLE MIXING RATIO against `n_tot`, which ionization and dissociation do
> not conserve, so the composed map carried a spurious `n_c v d ln(mbar)/dr`
> living exactly across the dissociation and ionization fronts. Section 158
> moves the row onto the fraction per unit mass, the scalar the conserved `rho`
> actually carries, keeping the diffusive flux on the particle mixing ratio
> where it belongs. The base cell, the one place with no `dt` trend before, is
> first order after.

**(3) The reading this supports, and it is the uncomfortable one.** If the creep
is physical and is neither the base cells (section 9) nor this inconsistency
(which is largest where the front does not respond, and 3 per cent where it
does), then **the transported steady state's front may sit much further out than
anything measured -- H2 through the whole wind -- and the `du` stop may simply
have cut the relaxation short.** A long marching run judged on physical time
rather than step count is under way. This is an interpretation offered by the
source document, not a measurement.

### What this does not establish

It does not establish that no transported steady state exists; it removes
explanations. The measurement that would settle it is the one this line of work
has always needed: a state at which an additional transport relaxation leaves
the front where it is. Also not measured: the base inflow Dirichlet in
isolation, the 500-against-1000-cell comparison at matched physical time (it was
running from a cold start and had not finished), and any planet but this one.

### Status

**New, open. Whether the creep is a numerical artifact or physics is not
decided** -- but (1) removes the step-bias reading, and (2) is a defect in its
own right whatever the answer to (1) turns out to be.

**Two corrections to how this item was described to us.** First, the front does
**not** creep to the end of the domain: over 40,000 to 243,458 steps the
`f(H2) = 0.5` front moves 1.0485 to 1.0586 `R_p` and `f(H2) = 1e-4` moves 1.1241
to 1.1623. (A separate *Newton finish* under the same physics does put
`f(H2) = 1e-2` at 1.4860 `R_p`, most of the way to 1.5 -- but that state is
`info = 2` and the document says of it "It is not a converged state and none of
it is a result." Reading (3) above is the hypothesis that the marched front
*would* go there, not a measurement that it does.) Second, the `dt` scan of (1)
is a CFL scan at matched physical time, **not** a `dt`/`dr` discriminator: no
`dr` leg was run at matched physical time, and the document lists that
comparison as unfinished.

---

## What is NOT a defect: approximations with stated validity

These are deliberate, are written down at the code site with their range, and
should not be re-opened as defects.

* **The flux gate's window, `r >= r_flux = 1.2 R_p`.**
  `steady_residual.f90:509-523` states the reason: the base boundary
  deliberately admits a small inflow through the smooth valve and the wind is
  still being launched just above it at `|v| ~ 1-10 cm/s`, eight orders below the
  wind speed, so including the launch region would measure the boundary condition
  rather than the wind. Over all cells the spread of a converged hot Uranus is
  9.9e2 and `F` changes sign, so max/min is not even defined there. Item 6's
  point is not that this reasoning is wrong; it is that a factor-2 gap in the
  **energy** budget over the same cells is a different statement and nothing
  measures it.

* **The two-region split of the residual norm, and taking the maximum rather
  than a sum.** `steady_residual.f90:444-469` states it: the wind's rows and the
  layer's rows are made dimensionless by different physics, so the two numbers
  are not addends of a common quantity, and the maximum says that every part of
  the column is steady on the scale its own physics sets. A volume- or
  cell-count-weighted sum over both regions would let the outer region average
  the inner one away, which is how a wind an order of magnitude from steady
  between 1.2 and 2 R_p was once scored better than one that was not.

* **The carrier gate being a third gate rather than a row of `||R||`.**
  `steady_newton.f90:1019-1027`: `Resid tol` is calibrated against a row scale
  that bounds each hydrodynamic row's largest term and converged states sit at
  1e-6 on it, while the carrier row is measured against the terms themselves and
  a converged row sits near 1e-3. One number cannot hold both tolerances.

* **The marching loop passing `carrier_is_unknown = .false.`.**
  `steady_residual.f90:612-619` and `EXHALE_main.f90:1089-1093`: its carriers are
  moved by the operator-split transport step, not solved for, so no carrier
  residual of a solve exists for that state, and gating its stop on the last
  solve's number would report a measurement the state never had.

* **`p_single` carrying no trapping.** `h2_self_shielding_table.f90:24-31` and
  1033-1038: the count of fluorescent decays per dissociation is
  `(1 - p)/p` with the untrapped `p`, because a trapped decay and the
  re-absorption that follows it move one molecule down and another up,
  depositing nothing. Item 9 applies to `sigma_diss` and `p_eff` only.

* **Only H2 carries a molecular heat capacity in the caloric equation of state.**
  `caloric_eos.f90:35-48`: H2+, H3+, HeH+, OH, H2O and CO are counted at
  `(3/2) k`. Measured on the converged He/H = 0.0793 hot Uranus over every cell
  where H2 exceeds 1e-3 of the particle count, the largest ratios to H2 are
  `n(H2+)/n(H2) = 1.7e-7`, `n(HeH+)/n(H2) = 1.2e-7` and
  `n(H3+)/n(H2) = 3.2e-8`; giving all three the full H2 heat capacity would move
  `C_V` by less than 1e-6.

* **A file with no `# coupling:` header starting staged.** `EXHALE_main.f90:598-605`
  and item 11: nothing else can be inferred about such a file, and arming it
  blindly is measured to be worse.

* **The positivity fallbacks in `Apply_BC`.** `Apply_BC.f90:157-167` and 190-218:
  where a linear ghost extrapolation leaves the admissible state the code drops
  to the zero-gradient copy, which is the first-order limit of the same boundary
  condition, and counts every occurrence in
  `n_ghost_cells_positivity_limited` so a run that never needed it is the run an
  unguarded build would have produced.

* **`Escape radius` and `r_flux` clamping to mid-domain on a truncated grid.**
  `define_grid.f90:200-244`: a deep-RLOF case whose Roche lobe is inside the
  default escape radius would otherwise leave an empty convergence range, and
  `maxval`/`minval` over an empty array would give a spurious converged exit at
  step 0. The clamp warns loudly.

---

## Statements in the documentation that this catalog found to be out of date

Recorded here rather than corrected, because another worker is editing those
files.

* **`docs/p23_thermal_budget.md` section 9.3** says
  "`parameters.f90` line 461 fixes `gamma = 5/3` for the whole domain, including
  the layer that is 80 percent H2 by volume at 1000 K, where the rotational
  degrees of freedom make the effective value nearer 7/5 ... This was not
  tested." That is superseded. `src/modules/states/caloric_eos.f90` implements
  the composition- and temperature-dependent `gamma_eff` and is switched on
  unconditionally for any run carrying molecules
  (`caloric_eos.f90:178-205`: `caloric_mixture_active` is set true for any cell
  with a positive particle count and positive H2, and false only when
  `thereis_mol` is false). `TO_BE_DONE.md` (Z) records the item as closed on
  2026-09-03.

* **`docs/p54_base_layer_mass_flux.md` now carries section 10** (this entry
  originally recorded its absence; the section was brought into the repository
  after the catalog was written). Citations of sections 10.4/10.5 resolve.

* **`docs/p44_base_sawtooth.md`** carries three of its own corrections that a
  reader should not miss: section 8.5's "`T0` is applied at the wrong level" is
  **withdrawn** (lines 98-100); section 7.3's "the base thermal structure does
  not converge" is **corrected** by section 8.3 as a sampling artifact of a steep
  profile at a moving point (lines 636-640); and reading (iii), a physical base
  oscillation, is **excluded** on the grid evidence (lines 455-464).

---

## What this catalog did not verify

* Every line number above was read in the tree on 2026-09-03; every measured
  number was taken from the document cited and was **not** re-measured, with one
  exception: the volume and cell-count fractions worked out for the first
  edition of item 8 were computed from
  `backup/regression/mol_sec_ion/output/Hydro_ioniz.txt` (read only, nothing
  re-run). Item 8's settling measurements are `docs/p55_base_mode.md` section 11
  and were not re-measured here either.
* No run was started for this catalog and no output was regenerated. In
  particular the section 145 patch was **not** built or run here; its state is
  reported from reading `apply145_resid.py`, `apply145_main.py` and
  `t145/src/` in the P55 scratch tree.
* Items 12.7, 12.8 and 12.9 are recorded from `TO_BE_DONE.md` headings and the
  review's judgment; the source sites behind them were not read. The review
  itself says of 12.8 that "the exact offending refresh [was] not isolated here",
  so that item has been read by neither of us in the code.
* The review's own scope limits carry over: it reran no production calculation,
  and its conclusions marked "reported measurement" were not re-measured. Where
  this catalog and the review disagree with a third document, no arbitration was
  attempted beyond what is written in each section.
* The claim that the cellmax norm "changes nothing else" rests on two states.
  `p55_base_mode.md` section 11.7 records that it was not run on the regression
  matrix, that no threshold for a face-flux gate has been derived, and that
  nothing establishes whether the base state the cellmax norm forces the solve
  to reach is physical.
