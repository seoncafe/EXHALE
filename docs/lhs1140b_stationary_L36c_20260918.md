# LHS 1140 b, item L36c: the ionization stage row solved in the run

Second increment of item L36 of `docs/PLAN_20260917.md` section 10, which
carries stages B to F of `docs/lhs1140b_stationary_L12a_design_20260913.md`
as `docs/lhs1140b_stationary_L12b_derivation_20260916.md` sequenced them.
It follows `docs/lhs1140b_stationary_L36_20260918.md`, which built the stage
flux as a module with no caller in a run, and closes the state that memo's
section 1 called "a state to leave as soon as possible and not a design".

Every number is MEASURED (run or computed here) unless marked READ.

---

## 0. Verdict, before the detail

1. **There is now ONE spelling of the ionization transport, and it has a
   caller in a run.** The H II row of the transport-chemistry operator is
   `x(H II) = n(H II)/n(H nuclei)`, a fraction per hydrogen NUCLEUS, whose
   face flux, cell divergence and tridiagonal rows come from
   `ionization_stage_transport` written on the element operator's own
   hydrogen nucleus face flux. The bulk-mass spelling is gone: the proton's
   old gradient and drift face coefficients, its bulk-flux advective
   divergence and its nucleus-budget headroom are all switched off for that
   row, and nothing rebuilds a face coefficient of the element operator.

2. **x(He II) and x(He III) are NOT carried, and that is a blocker in
   files this item does not own, not a choice.** The transported fraction
   reaches the ionization balance through `x_hp_fixed`, and the cell state
   has one imposed fraction, not three
   (`ion_cell_state.f90`, `ion_residual_core.f90`
   `impose_transported_ionization_fractions`, READ: "The cell state states
   no imposed helium fraction, so neither row is written"). Carrying a
   helium stage without that would transport it and let the next sweep
   throw it away. Section 2 names the three edits.

3. **The atomic-gas refusal STAYS, for the same reason, and its message now
   says the real one.** `ionization_equilibrium.f90` imposes the
   transported fraction inside a block gated on
   `thereis_mol .and. carrier_transport` (line 1922, READ), so in an atomic
   gas the next sweep would put the proton back on its local root and the
   transport would leave no trace. The refusal is therefore correct and it
   is kept; what changed is that it now states the SWEEP as the reason
   instead of the operator. Two refusals were added, both for
   configurations the operator does not carry: a mixture with no helium
   (the element nucleus flux is the binary H/He operator's and returns
   nothing for one element), and `Coupled carrier solve: True` (its row
   registry carries every carrier unknown as the species mass fraction of
   its `f_sp` column, and a stage unknown is a fraction per element
   nucleus).

4. **Key-off byte identity: MEASURED** on `wasp_full`, `mol_diffusion` and
   `mol_carrier` (section 4).  `hp_front` is a KEY-ON case and belongs in
   the table below, not in this one.

5. **Key-on movement: MEASURED** on the three `hp_*` cases (section 5). The
   largest is 5.3e-3 relative in `hp_zero_seed`, on `n(H I)` at the cell
   that carries the H2 front; `hp_front` moves by 4.1e-5. No golden was
   refreshed. Those goldens are ALREADY stale against the entry text by
   3.9e-3 and 3.3e-3, which is the size of the movement this increment
   adds, so the refresh is one refresh and not two.

6. **The certification reports the stage sum identity** as
   `ionization stage nucleus sum`: at every face the stage fluxes of an
   element must sum, over all its stages, to that element's nucleus flux.
   The entry REPORTS and does not gate, because the tolerance is not
   anchored; section 7 gives the anchoring measurement and the proposed
   anchor.  It reads 2.268e-16 and 2.196e-16 on the two solved states of
   the one fixture that runs the key with the partitioned stationary route
   (section 6).

7. **Stage D cannot be evaluated and stage E is partly owed.** The He 2^3S
   promotion rule of L12a section 3 needs a state whose three ionization
   fractions agree to their own gate; two of the three are not carried, so
   the rule has no precondition to be evaluated on. Section 8.

---

## 1. What the row is, term by term

The unknown of the H II row is

```
    x(H II)(j) = n(H II)(j) / n_H(j) ,
```

`n_H` the hydrogen NUCLEUS density of the cell, taken from the one
stoichiometric map of the species table (`element_nucleus_counts`, now
exported) so that the cell value and the face value the flux divides by are
one count. It is frozen with the rest of the background of a step
(`cbg_nHnuc`): transport moves a nucleus from cell to cell, and which stage
it sits in is what the row solves.

The face flux is equation (1) of `ionization_stage_transport`,

```
    F(f) = x(f) N_H(f)  -  n_H(f) K(f) [x(j+1) - x(j)] / dr(f) ,
```

with `N_H(f)` the hydrogen nucleus face flux of the element operator
(`element_nucleus_face_flux`, through
`hydrogen_and_helium_nucleus_face_flux`), `x(f)` the reconstructed, limited
and upwinded face fraction with the donor side taken from the sign of the
ELEMENT flux, and the eddy term zero at `f = 0` and `f = N`.

**Which halves of `N_H` this run carries.** The stage row rides on the
element flux of THIS run and not on the flux a differently configured run
would have:

| half | present when |
|---|---|
| advective, `F_rho Y` | the row itself carries the material advection, that is the fixed-wind relaxation and the stationary balance; in the marching loop the Runge-Kutta stages carry it and the caller passes a zero face mass flux, exactly as `element_transport_residual` and the molecular carrier rows do |
| the element's diffusive, gradient, eddy and settling flux | `He_diffusion` with helium in the mixture. With that option off no element crosses a face by diffusion in this run, so a stage riding on one would be carried by a flux the elements themselves do not have |
| the stage's OWN eddy term `-n_H K dx/dr` | always: `K_zz` is a bulk mixing coefficient of the gas and is blind to which element it mixes |

The physical face mass flux is handed to the element operator in
g cm^-2 s^-1: the rows' own `F_rho` is the code face mass flux with the
mixture mass still outside it, so the mixture mass enters at the face on
the arithmetic mean every other face coefficient of that operator uses.

The rest of the row follows from the variable:

| term | before | now |
|---|---|---|
| time term | `rho n0 (y - y_old)/dt` | `n_H (x - x_old)/dt`. With the element continuity `dn_H/dt = -div N_H` holding, the nonconservative form is the conservative `d n(H II)/dt + div F` exactly, which is what the element row is written on |
| advective divergence | `div(F_rho Y_HII)` on the BULK face mass flux | inside `F(f)` above, on the element nucleus flux |
| diffusive and eddy flux | `-n_tot (D + K) dX/dr` with `D = 0` and `X = n(H II)/n_tot`, the WHOLE mixing-ratio flux | `-n_H K dx/dr`, the stage term alone; the element's own eddy flux is already inside `N_H` |
| chemistry | `src(ic_Hp) = fv(1)`, a volumetric rate | unchanged; the row's source is a rate and not a fraction |
| Jacobian | donor-cell linearization of the bulk advective face flux plus the diffusive face derivatives | `ionization_stage_face_jacobian` and the same divergence geometry, which reproduces `ionization_stage_row_jacobian` entry for entry |
| admissible set | the shared hydrogen nucleus budget | the SIMPLEX, `stage_simplex_projection`. The shared budget still binds the whole H-bearing carrier set, because `carrier_source` closes atomic hydrogen with the proton inside that sum |
| write-back | `f_sp(H II) = y` | `f_sp(H II) = x n_H / (rho n0)` |

**The boundary rows.** The two end faces carry a live stage flux, because
the material advection crosses them; what vanishes there is the eddy term
alone, and `ionization_stage_face_jacobian` has already set that half to
zero at both. Nothing below the base states an ionization fraction, so the
inner ghosts carry cell 1's own partition, and the outer ghost is cell N
continued; both are copies of the unknown, so the derivative of the end
face folds onto the DIAGONAL of the end row. That folding has its own
acceptance row (section 4) and its own RED.

---

## 2. What the helium stages and the atomic gas still need

Three edits in two files this item does not own, and one in a third:

1. `src/modules/nonlinear_system_solver/ion_cell_state.f90`: two flags and
   two values beside `x_hp_fixed` / `x_hp_fix`, for `x(He II)` and
   `x(He III)` per helium nucleus.
2. `src/modules/nonlinear_system_solver/ion_residual_core.f90`: two lines
   in `impose_transported_ionization_fractions`, which already says where
   they go ("He II and He III per He nucleus are rows 2 and 3 and enter
   here when the helium stages are carried as well", READ). The seven
   `System_*` modules already call that one routine, so nothing else in the
   sweep changes.
3. `src/modules/radiation/ionization_equilibrium.f90` line 1922: the
   ionization part of the imposed-fraction block gated on
   `ionization_transport .and. (bg_ready .or. do_load_IC)` alone instead of
   on `thereis_mol .and. carrier_transport .and. ...`, which is what lets
   an ATOMIC run carry the transported partition.

With those three in place the carrier operator takes the helium rows by
adding two indices to `carrier_set_init` and one line each to
`carrier_is_ionization_stage`, `carrier_nucleus_reference` and
`carrier_stage_face_state` (which already forms `N_He`, `n_Hef` and
`cbg_nHenuc` and throws them away), and the two refusals in `input_read`
are lifted in the same increment. The module carries them already: the
simplex projection, the closing stage and the sum identity are written for
`n_k` stages and not for one.

`src/modules/time_step/steady_newton.f90` `set_transported_species_rows`
(line 2370, READ) is a separate matter and is NOT needed for the
partitioned route: it gates the carrier rows of the COUPLED solve on
`thereis_mol .and. carrier_transport`, and that solve is refused with this
key for the variable reason of section 0.3.

---

## 3. Tests

`src/tests/ionization_stage_flux/run.sh`, private `OBJDIR=build_L36c`,
`OMP_NUM_THREADS=1`: **TOTAL 27 passed, 0 failed** (MEASURED). The 25 rows
of L36 are unchanged and reproduce that memo's values exactly. Two rows are
new, and they are the ones this increment's boundary rule needs:

| row | measured | tolerance |
|---|---|---|
| `boundary_stage_rows_fold_the_ghost_onto_the_diagonal` | 1.8552e-15 | 1e-7 |
| `dropping_the_end_face_entry_changes_the_end_rows` (RED by construction) | 9.0722e-01 | must exceed 1000 eps |

The first takes the two folded diagonals the transport operator assembles
at `j = 1` and `j = N` and compares them with a central difference of the
stage divergence with the ghosts continued the way the operator continues
them. The second is why the fold is needed: a purely diffusive carrier row
drops the end-face entry because its flux vanishes there, and a stage row
that did the same would drop the whole material advection through the
boundary face, which changes the end rows by 0.91 of themselves.

The suites that link the changed objects, all on the measured build:

| suite | result |
|---|---|
| `ionization_stage_flux` | 27 passed, 0 failed |
| `element_operator` | 40 passed, 0 failed |
| `species_face_flux` | 15 passed, 0 failed |
| `steady_species_rows` | every assertion passed |
| `carrier_reference_scales` | 14 passed, 0 failed |
| `carrier_boundary_jacobian` | 12 passed, 0 failed |
| `carrier_outer_boundary` | 16 passed, 0 failed |
| `carrier_constraint_attribution` | 10 passed, 0 failed |
| `carrier_returned_state_acceptance` | 36 passed, 0 failed |
| `ionization_imposed_fractions` | 31 passed, 0 failed |
| `certification` | every assertion passed |
| `attempted_step` | every assertion passed |
| `carrier_retry` | **one row RED**, section 9 |

---

## 4. Key-off byte identity

Control build `EXHALE_L36c_ctl.x`, md5 `d15a4b54fd5fbb8571bb374b7da9270e`,
a build of the ENTRY TEXT of every file of this item (it is the binary the
L36 increment delivered, so the entry text of L36c is the delivered text of
L36). Measured build md5 `1d5e414ad3e21a27bcf21fe03a38d179`; the DELIVERED
build is `91a0dca3daa3e72cf6ea5b61a16c52d6`, which differs from it by two
comment blocks and by `if (ic .eq. ic_Hp)` written as
`if (carrier_is_ionization_stage(ic))`, whose function body is that same
expression.  The two give BYTE-IDENTICAL data rows on the three `hp_*`
cases, which are the cases that enter the changed branch at all; with the
key off `ic_Hp` is outside the solved set and the branch is not reached.
Scratch case copies, `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`, identical
bounds on the two builds.

| case | steps | `Hydro_ioniz.txt` | `Ion_species.txt` |
|---|---|---|---|
| `wasp_full` | 400 | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |
| `mol_diffusion` | 12000 (its pinned count) | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |
| `mol_carrier` | 12000 (its pinned count) | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |

The two builds also agree on every printed step line of the UNCAPPED
`wasp_full`, over the 2920 steps both had taken when that pair was stopped
and re-run at the 400-step bound, and on every step line of the two
molecular cases over all 12000.

Every one of these carries the key OFF. `mol_carrier` and `mol_diffusion`
run the carrier transport operator itself, which is the operator this item
rewrites a row of, so they are the cases that reach the changed code with
`carrier_solved(ic_Hp)` false; `wasp_full` reaches it through the He 2^3S
and metal configuration with no molecules at all.

`EXHALE.x` was never touched and nothing was written under
`backup/regression/`.

---

## 5. The key-on movement

The three `hp_*` cases carry `Ionization transport: True`, and they are the
cases whose row this item re-spells, so they move. The movement is the
largest relative difference over the 504 data rows of each file, taken
column by column.

| case | `Hydro_ioniz.txt` | `Ion_species.txt` | largest mover |
|---|---|---|---|
| `hp_front` | 3.0201e-05 | 4.1478e-05 | `n(H I)` at cell 32, 2.735543e+05 to 2.735429e+05 |
| `hp_zero_seed` | 4.2160e-03 | 5.2914e-03 | `n(H I)` at cell 193, 7.658211e+04 to 7.698949e+04 |
| `hp_trace_seed` | 3.0318e-03 | 3.5450e-03 | `n(H I)` at cell 196, 8.202580e+04 to 8.173502e+04 |

`x(H II)` itself moves by 1.2851e-10 (`hp_front`), 6.0704e-08
(`hp_zero_seed`) and 1.4259e-07 (`hp_trace_seed`); what moves is the
helium partition and the neutral hydrogen at the one or two cells that
carry the H2 and ionization front.

**The goldens are already stale against the entry text by the same size.**
The control build against the current `backup/regression/golden/`:

| case | control against the golden, `Hydro_ioniz.txt` | `Ion_species.txt` |
|---|---|---|
| `hp_front` | 8.8558e-09 | 1.7538e-07 |
| `hp_zero_seed` | 3.2967e-03 | 3.9439e-03 |
| `hp_trace_seed` | 2.7650e-03 | 3.3068e-03 |

so `hp_zero_seed` and `hp_trace_seed` had already moved by 3.9e-3 and
3.3e-3 before this item touched them (the goldens were refreshed on
2026-09-17 for the L7g molecular change, on a build the tree has since moved
past through L30, L30b and L36). No golden was refreshed here.

**Where the movement sits.** In every case the largest movers are
`n(H I)`, `n(He II)` and `n(He III)` at one or two cells around the H2 and
ionization front, and `x(H II)` itself moves by 1.3e-10 (`hp_front`) to 1.4e-7
(`hp_trace_seed`). That is the signature of the change: the three
cases carry `K_zz = 0` and no `He_diffusion`, so with the key on and the
rows marching (no fixed-wind relaxation runs in them) the stage flux has
neither an eddy term nor an element diffusive half, and the H II row is the
same equation in a rescaled variable. What is left is the different
conditioning of the Newton -- the row is now scaled by the hydrogen nucleus
density instead of by `rho n0` -- and the stages the write-back rescales
against a proton that stopped at a slightly different point.

---

## 6. The key-on stationary solve

`backup/regression/carrier_model_a_newton`, the one fixture that carries
`Ionization transport: True` together with `Restart intent: stationary`,
`Solver: Newton 100.0` and `Coupled carrier solve: False`, which is the
partitioned route. Run on a scratch copy with `EXHALE_PTC_DTAU0=1.0`,
`OMP_NUM_THREADS=1`.

**The stage sum identity on the solved state, MEASURED:**

| state | `ionization stage nucleus sum` | face |
|---|---|---|
| the loaded restart, first evaluation | 2.268e-16 | 485 |
| the state the first JFNK handed back | 2.196e-16 | 284 |

against 2.0245e-16 on the suite's synthetic 500-cell column (L36, READ) and
against 1.5e-4 to 3.6e-4, which is what the suite's RED rows measure when
the construction is broken. **The identity holds to the rounding of the
sums on a real solved state**, which is the acceptance this item owed.

**The row itself is not certified, and it was not before either.** The H+
carrier balance is gated at `cert_tol_carrier_at` = 1e-5 and reads, on the
first evaluation of the loaded restart:

| build | gated row measure of `carrier balance H+` | cell |
|---|---|---|
| control (bulk-mass spelling) | 6.392e-03 | 248 |
| measured (stage spelling) | 6.432e-03 | 248 |

so the row is above its gate by the same two and a half decades in both
spellings, and this item neither fixed that nor made it worse. What that
gate is worth for an IONIZATION row is the open question of L12a section
2.4: its five anchors are all measurements of a TRANSPORT row, whose
balance is a cancellation between an advective and a diffusive divergence,
and an ionization row is a cancellation between an advective divergence and
a reaction rate. **No tolerance was changed.** The anchoring measurement
that item owes is a smooth manufactured ionization column, a chosen
`x_eq(r)` with the `P` and `L` that produce it, refined by a factor two at
the production spacing; it is not made here.

---

## 7. The tolerance the new entry would take

`ionization stage nucleus sum` REPORTS and does not gate. The anchoring
measurement for it, unlike the row tolerance above, is not a refinement
study: the identity is exact in exact arithmetic, so what a tolerance has
to clear is the rounding of the sums at the production size.

| where | measured |
|---|---|
| suite, synthetic 500-cell production Mixed grid | 2.0245e-16 (READ, L36) |
| the same after the simplex projection | 2.1653e-16 (READ, L36) |
| the same at the boundary faces | 1.0621e-16 (READ, L36) |
| the run's loaded restart, 500 cells | 2.268e-16 |
| the state the run's first JFNK handed back | 2.196e-16 |
| the suite's RED rows, a broken construction | 1.4761e-04 to 3.6208e-04 (READ, L36) |

**Proposed anchor: 1e-13**, about 500 times the largest rounding measured
and nine decades below the smallest break the RED rows produce. It is
proposed and not set.

---

## 8. Stage D, and what stage E still owes

**Stage D cannot be evaluated.** The rule of L12a section 3 is: if the
post-process moves the He 2^3S column by more than 5 per cent **when the
three ionization fractions already agree to their own gate**, the level
joins the carried set. Two of the three are not carried (section 2) and the
third does not meet its gate on the one fixture that runs the route
(section 6), so the rule has no precondition. The level is NOT carried,
which is the side its own Damkohler measurement supports (`Da(2^3S)`
smallest 15.7 anywhere in the physical domain, READ, L12a section 1.1), and
the systematic it would be judged against was measured by L36 on the
certified atomic fiducial: the He 2^3S radial column differs between the
two closures by a factor 6.15 (READ).

**Stage E is half done.** L12a section 5.2 defines stage E as the
certification entries and the tolerance anchoring of section 2.4. The stage
sum entry is built with its measurement and its proposed anchor (section 7).
The row entry of each carried stage already exists -- `carrier balance H+`
is one --
and their tolerance is the unanchored one of section 6.

**What L12b calls stage F is the `_adv` retirement AS THE SOURCE OF THE
CARRIED FRACTIONS** (its section 8), and that landed with L36 on the
transit side: `EXHALE_TRANSIT_STATE` selects the input pair, the pair is
named in every product header and the measured distance between the two
states travels with it. `post_process_adv` is not retired; it stays the
independent discretization. So stage E is the certification and stage F is
the line, and the two are not the same item.

---

## 9. Noticed outside this item's scope

1. **`src/tests/carrier_retry/` has one RED row on the measured build**, and
   it is a knife-edge one. `one_further_sweep_leaves_the_returned_composition_within_the_sweep_tolerance`
   reads 1.034627e-06 against `ieq_res_tol` = 1e-6; the control build reads
   2.685237e-07. The fixture runs WITH the key on
   (`setup_globals` sets `ionization_transport = .true.`), so this is a
   key-on row and the movement is the intended physics. **The binding
   species is HeH+ at 1.33e-12 of the gas in the base cell** (f_sp column
   37, cell 1, from 1.329162e-12 to 1.329161e-12), which the fixture counts
   as "present" because its threshold is 1e-20. A diagnostic line naming
   the species and the cell was added to the fixture, so the row now says
   what is not closed; **the tolerance was NOT touched**, because no
   tolerance may be chosen to make a snapshot pass. What the row is
   measuring is the sweep's own closure of a trace molecular ion thirteen
   decades below the gas, and the fixture decision -- raise the "present"
   threshold, or accept the row as a statement about the sweep's noise --
   is the advisor's.

2. **`docs/input_schema.md` K15f was stale on two counts and is rewritten.**
   It stated a `Coupled carrier solve` requirement and a fatal `error stop`
   that was removed on 2026-09-12 (L12a section 7 and L12b section 9 both
   record it, READ). The same stale sentence is repeated in
   `Update_EXHALE_stage1.md` section 169.1; that file is the record of a
   landed series and is not edited here.

3. **`EXHALE_main.f90`'s `EXHALE_PTC` refusal carried a stale reason.** Its
   comment said it refuses "for the reason input_read refuses
   `Solver: Newton` with the ionization state carried"; that refusal was
   removed on 2026-09-12 and the partitioned route is supported. The
   comment now says what is refused and what is not. The refusal itself is
   unchanged and is still correct: the direct steady route solves
   `F(Y) = 0` with the composition on its own local root.

4. **`carrier_element_headroom` gained a fifth argument**, the element's
   nucleus density, because an ionization stage is bounded by the simplex
   and not by the reservoir the molecular carriers share. Its one
   production caller and the `steady_species_rows` rows that read it were
   updated, and the row `headroom_Hp_is_H` became
   `headroom_Hp_is_the_H_nucleus_density` with two different numbers in the
   call so that it states which budget the stage reads.
