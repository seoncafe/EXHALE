# D5b-2: one boundary-state operation, with a stated input set

Item D5b-2 of `docs/PLAN_20260918_rev2.md` (section D5, and the D5b-2 row of
"Decisions made"), the second increment of D5b. It answers
`docs/lhs1140b_stationary_D5a_20260918.md` section 9 items 2 to 5 and closes
the two evaluation sites D5b-1 left named
(`docs/lhs1140b_stationary_D5b1_20260918.md` section 8). 2026-09-18.

Every number below is labelled MEASURED (this item ran it) or READ (from a
file, a log or a source line).

---

## 1. Verdict

**The lower boundary is now a function of the physical column, the prescribed
reservoir, the radiation context and the model options, and of nothing else.**
The ghost rows of a restart pair are not read; the ghost's molecular partition
is solved against the ghost's own ionization balance instead of being
prescribed from the composition the state was entered with; the prescribed
reservoir and the ghost's own particle counts are two named quantities and the
solved one is never written back; the cached face state carries the input set
it was derived from and a read that finds the composition moved derives it
again; and the two evaluation sites of a solve install the boundary after the
sweep, as the three evaluation routes of D5b-1 already do.

MEASURED, the four ghost variants of D5a table 3 on the hot-Uranus molecular
reference state, each an admissible ghost supplied through the two lower
species rows of the restart pair, with every physical cell and every declared
input held fixed:

| ghost composition supplied | entry text, cell-1 row | this item |
|---|---|---|
| the state's own | 7.0724E-09 | **1.5582E-08** |
| the archived ghost rows of the L34b pair | 7.0724E-09 | **1.5582E-08** |
| `n(H2)` scaled by 1 + 1e-6 | 1.6505E-09 | **1.5582E-08** |
| the composition of the first physical cell | 2.5418E-04 | **1.5582E-08** |

and on this item's build the four runs' `Hydro_ioniz.txt` and
`Ion_species.txt` are **identical in every data line** (md5 of the non-comment
lines), where the entry text's three differ. Five decades of dependence on a
row nothing in the physical state selects are gone.

**The re-entry ladder is flat from the first re-entry.** MEASURED, the same
state fed its own product four times: 1.5582E-08 at E1, E2, E3 and E4, against
D5a's 15.98 / 9.06 / 7.14 floors and a 2-cycle that never settled. E1 = E2 is
the acceptance item D5a section 9.6 (c) asks for, and the alternation it
recorded was the ghost's iteration history, which the closure removes.

**The physical column still round trips at the last bit.** MEASURED,
conserved density, momentum and ENERGY and every species fraction, cells 1 to
N, through a write and a read: 0, 2.19e-16, 4.03e-16 and 6.07e-16 relative,
which are D5a's numbers unchanged, and the ghost rows of the file are not read
at all.

**A formerly certified state is refused, and that is the item's own
statement.** The corrected boundary moves the cell-1 continuity row of the
three molecular states from 9.06, 3.09 and 0.08 rounding floors to 19.97,
110.0 and 27.13 (section 5). The cause is measured: the ghost's NEUTRAL
columns and its elemental budget are reproduced to 1e-10 and its He/H exactly,
while its trace IONS move by 3e-3 to 5e-2, which is the band the sweep's own
acceptance leaves a ghost cell (a normalized reaction residual of at most
`ieq_res_tol` = 1e-6 admits that much in a trace ion), and the base continuity
row carries a gain of order 1e9 floors per relative unit of the ghost
composition (D5a table 3). **The base continuity row of a molecular state is
therefore not reproducible beyond the ionization solve's own acceptance band
on the ghost cell**, whatever the boundary does; that is D5a section 8's last
line measured against a second prescription, and it is why D5a section 9.6
refuses to let that row be an acceptance criterion. No tolerance was touched
and no golden was refreshed.

**What this item does NOT do.** The inflow/outflow rule is untouched: the
smoothstep weight, the C- relation and the reservoir isentrope are the current
model's, and the model now carries an identity string
(`characteristic_face_ps_reservoir_C_minus_smoothstep_v1`) in the setup report
and in the provenance header of every state file. A physical change of that
rule is D2b's.

---

## 2. The interface

### 2.1 The boundary-state operation

`Apply_BC(u)` (`src/modules/states/Apply_BC.f90`) is the operation, with its
contract stated at the routine:

| | |
|---|---|
| input | the physical conserved state `u(:,1:N)`; the composition installed for it, through the composition-derived state `caloric_cell_mixture` and `n_part_cell1` carry; the PRESCRIBED reservoir of `base_boundary` with its version; the model options |
| not an input | the ghost rows of a restart file (`load_IC` replaces them at the read) |
| output | the ghost conserved state `u(:,1-Ng:0)` and `u(:,N+1:N+Ng)`, and the cached face states `base_face_W`, `base_ghost_W`, `base_face_lower_W`, all of ONE composition |
| the tag | `bc_u_hold`, `bc_W_hold`, `bc_nk_hold`, `bc_xh2_hold`, `bc_mol_hold`, `bc_npart1_hold`, `bc_res_hold`, `bc_res_version_hold`, written by `hold_boundary_state_inputs` |

The ghost COMPOSITION is not written there: it is solved by the composition
sweep and `Apply_BC` turns it into the ghost's conserved and caloric state.
The two together are the boundary state.

New public names of `BC_Apply`:

- `base_boundary_closure_is_current()`, the half of the input set reachable
  without an argument (the composition-derived state of cells `1-Ng` to 1, the
  cell-1 particle count, the reservoir and its version);
- `base_boundary_cache_is_current(u)`, the whole input set, for a consumer
  that holds the state;
- `recompute_base_face_states`, the derivation again from `bc_W_hold` and the
  composition installed now, called by `Rec_BC` when the cache does not belong
  to it;
- `n_base_face_state_recomputed` and `report_base_face_cache_reads(tag)`, the
  count and its report.

### 2.2 The reservoir, its version and the model identity

`base_boundary` (`src/modules/states/base_boundary.f90`):

| name | what it is |
|---|---|
| `base_reservoir_p`, `base_reservoir_T`, `base_reservoir_nhat`, `r_base_level` | the PRESCRIBED reservoir, set once by `set_base_reservoir` at startup and never written again |
| `base_reservoir_prescription_version = 1` | what those four numbers mean: (p, T, particles per unit mass, level radius) in the code's own units, `p = nhat rho T` at the level, the entropy carried by the isentrope through (p, T) at the base composition |
| `base_boundary_model_id = 'characteristic_face_ps_reservoir_C_minus_smoothstep_v1'` | which boundary model the face state is built by |
| `base_ghost_particle_count`, `base_ghost_electron_count` | the SOLVED ghost's counts, written by the sweep through `set_base_ghost_counts`, read by nothing but the report |
| `base_ghost_closure_residual`, `base_ghost_closure_passes` | what the ghost's molecular partition closed to, through `set_base_ghost_closure` |
| `report_base_boundary_model(tag)` | the block that prints all of it |

### 2.3 The ghost composition solve

`ionization_equilibrium.f90`, the base handoff block of the cell sweep. A
lower ghost whose molecular partition the handoff states takes as many passes
of its own cell solve as its closure needs: impose `x_H2 = x2 (1 - x_ion)`,
solve the ghost's ionization balance against it, impose again at the
ionization the solve returned. The radiation field is recomputed on the
self-field passes alone, so the closure passes solve the same cell against the
same field and a run at the default one self-field pass sees exactly the field
it saw.

| | value |
|---|---|
| `base_ghost_closure_tol` | 1.0e-10, absolute in `x_H2` (a fraction of the cell's hydrogen nuclei). As delivered by this item the iteration stopped when one more pass would move the imposed partition by less than this; since D9fix (`docs/lhs1140b_stationary_D9fix_20260919.md`) it stops when the residual below is under it, and the passes after the first two take a secant step on the scalar equation instead of substituting |
| `base_ghost_closure_passes` | 30 |
| the residual reported | `\|x_H2 - x2 (1 - x_ion)\|` at the state the cell returned |
| a failure | an error message naming the cell, the passes, the last move and the residual, and `error stop 'ioniz_eq: base ghost H2 closure failed'`. No fallback |

### 2.4 What a state file carries

`write_output.f90` writes two `#` lines on both halves of the pair:

```
# boundary_model characteristic_face_ps_reservoir_C_minus_smoothstep_v1
# boundary_reservoir version 1 p[p0] T[T0] nhat[n0/rho0] r_level[Rp]  <p> <T> <nhat> <r>
```

the four numbers at ES23.16E3, which round trips a double exactly. `load_IC`
reads the first (`parse_boundary_model_line`) and reports whether the state's
boundary model is this run's, this being informational: the boundary is
rebuilt from the physical column and this run's own reservoir whatever the
file says. `docs/input_schema.md` appendix D carries both lines and the
statement that the lower ghost rows are not read.

---

## 3. Binaries, snapshot and provenance

The working tree is under edit by other items, so the measurement works on a
frozen copy of `Makefile` and `src/`, and the tree's own `EXHALE.x` and
`build/` were never touched.

| | value |
|---|---|
| snapshot taken | 2026-09-18 19:00 KST, from the working tree at git `3c73905` plus its uncommitted files |
| entry-text manifest (md5 of the md5 list of every `.f90`/`.inc`) | `c249ccea3974100847cbde95d042e7fe` |
| control build (the entry text) | `bb4f21b61a59ba794f0d6431b99d3d56` |
| measured build (the entry text plus THIS ITEM'S FILES ONLY) | `bd3fe4981b54ecc1632ed4fdcc1c622a` |
| compiler | conda-forge gfortran 16.2.0, `-O3 -fopenmp`, OpenBLAS of the same prefix |
| every run | `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` |

**The measured build is the control text with this item's six source files
swapped in, and nothing else.** Another item (D3b) landed
`certification.f90` and `ionization_stage_transport.f90` in the live tree
between the snapshot and the first build, so a binary built from the tree
would have carried those as well; every number in this memo comes from the
isolated pair. MEASURED, the two builds were checked against each other on the
four evaluate-route states with both binaries: the tree build and the isolated
build read the same cell-1 row to every printed digit on all four, so D3b's
files do not reach these measurements, but the isolated build is what the
tables report.

---

## 4. The acceptance table

MEASURED, all of it on the isolated build, single threaded.

| the acceptance item (D5a section 9.6, the plan's D5b-2 row) | measured |
|---|---|
| the four ghost variants produce bitwise equal boundaries | **yes**: the four runs' `Hydro_ioniz.txt` and `Ion_species.txt` agree in every data line (one md5 each), and their ghost rows agree too |
| the four ghost variants produce bitwise equal first residuals | **yes**: the cell-1 row reads 1.5582E-08 in all four, and the signed base continuity rows with both face mass fluxes, written by the trace at sixteen digits, are identical |
| a re-assembly of the boundary after any report moves nothing | **yes**: the signed row, both face mass fluxes and the row scale of cell 1 are bitwise identical before and after the face-flux report, and the trace's boundary gap is `0.0000000000000000E+000` at both points. The only value that moves is the printed rounding FLOOR of cell 1, by one unit in the last place (7.8045986180219940E-010 against ...930E-010, one part in 7.8e15) |
| E1 = E2 on the re-entry ladder | **yes**: 1.5582E-08 at E1, E2, E3 and E4 |
| the physical column round trips at 1e-16 | **yes**: 0, 2.19e-16, 4.03e-16, 6.07e-16 on `u(1)`, `u(2)`, `u(3)` and the species fractions, and 4.02e-16 and 6.03e-16 on p and T |
| the ghost rows of the file are not read | **yes**: `load_IC` replaces both lower rows and says so once |
| the ghost composition solve reports its residual below its tolerance | **yes**: 0.00000E+00 in 4 passes on `HeH2.13` and on `.L22/i3_alt`, 1.11022E-16 in 3 on `HeH9.7`, 4.65605E-12 in 2 on the suite fixture, and 0.00000E+00 in 3 passes on the written states of `mol_base_handoff` and of `mol_carrier`, evaluated in their own configurations. `lower_profile` states no molecular network, so it has no ghost partition to close |
| the file schema stays readable by `examples/exhale_io.py` | **yes**: `load_hydro` and `load_ions` return 500 physical cells and the same values from a file of each build, and `physical_cell_rows` reads (2, 502) on both |

---

## 5. The states of D5a and D5b-1, before and after

MEASURED with `EXHALE_MASS_FLOOR_SCAN=1`, `Restart intent: stationary
evaluate`, the cell-1 continuity row and its own rounding floor.

| state | entry text | this item | floor | floors before | floors after |
|---|---|---|---|---|---|
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | 7.0724E-09 | 1.5582E-08 | 7.8046E-10 | 9.06 | **19.97** |
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | 2.3132E-09 | 8.2400E-08 | 7.4894E-10 | 3.09 | **110.0** |
| `.L22/i3_alt` | 6.3728E-11 | 2.1171E-08 | 7.8046E-10 | 0.08 | **27.13** |
| `atomic_scalar_gj1132_kzz1e9/HeH2.13` | 1.2016E-09 | 1.2016E-09 | 7.9044E-10 | 1.52 | **1.52** |

**The physical cells do not move at all.** MEASURED, control against measured
on all four states: `Hydro_ioniz.txt` and `Ion_species.txt` have ZERO
differing values over the physical cells; every difference is in the two lower
ghost rows. The evaluate route takes no step, so the column written is the
column read, and what this item changes is the boundary that column is
measured against.

**The atomic state does not move, and the reason is measurable.** The ghost
continuation of an atomic mixture reads `adiabatic_index_from_state`, which
returns `gamma_ad` with no composition in it (READ, D5a section 4), so the
replaced ghost composition has no path to the base face state.

**What in the ghost moved, and why the row moves by floors.** MEASURED, the
ghost cell 0 of the two molecular states, control against this item, column by
column:

| column | `HeH2.13` | `HeH9.7` |
|---|---|---|
| H I | 1.6e-10 | 4.0e-09 |
| H2 | 1.4e-10 | 3.8e-09 |
| He I | 4.1e-11 | 6.6e-11 |
| H II | 3.1e-03 | 4.6e-03 |
| He II | 3.1e-03 | 4.6e-03 |
| He III | 1.9e-04 | 1.1e-04 |
| He 2^3S | 5.9e-03 | 9.0e-03 |
| H2+ | 5.9e-03 | 4.7e-02 |
| H3+ | 1.3e-03 | 8.4e-03 |
| HeH+ | 5.9e-03 | 9.1e-03 |

and the ghost's ELEMENT budget is reproduced exactly: He/H reads the input's
9.7000000000 and 2.1300000000 in both builds and the hydrogen nuclei per unit
mass agree to every printed digit. So the elemental composition and the
neutral columns of the ghost are the same gas; what moved is the trace
IONIZATION, by 3e-3 to 5e-2.

That band is the sweep's own. A cell's composition is accepted on a
normalized reaction residual of at most `ieq_res_tol` = 1e-6, and in a trace
ion of a nearly neutral base that admits a few parts in a thousand: MEASURED,
two admissible starting points for the same ghost cell (the reservoir's gas
taken neutral, and the reservoir's gas with the first cell's ionization) land
3e-3 and 5e-3 from the archived ghost on the same column. The cell-1 row
carries a gain of order 1e9 floors per relative unit of the ghost composition
(D5a table 3: `n(H2)` scaled by 1e-6 moved it 425 floors), so an ion band of
1e-3 in a column that is 1e-6 of the hydrogen is worth tens of floors.
**The base continuity row of a molecular state is therefore not reproducible
beyond the ionization solve's acceptance band on the ghost cell**, which is
the conditioning D5a section 8 records, now bounded from a second
prescription.

The certification verdicts of the same runs, READ from the logs: the three
molecular states go from CERTIFIED (D5b-1) to NOT CERTIFIED on one refusing
entry each, the hydrodynamic mass row at cell 1; the atomic state is CERTIFIED
and unchanged in every printed measure. A corrected operator refusing a formerly certified state is the case
the plan's D5 section names ("A corrected operator may refuse a formerly
certified state, and that is reported, not absorbed"). Re-solving those states
under the corrected boundary is a catalog item and is not this one.

---

## 6. The ghost's molecular partition, closed

MEASURED on the hot-Uranus molecular reference state, `Restart intent:
stationary evaluate`, from the run's own report:

```
 [base boundary model] the boundary of this evaluation
   model  characteristic_face_ps_reservoir_C_minus_smoothstep_v1
   prescribed reservoir (version 1):
     p [p0]  8.4025828124448965E-01   T [T0]  1.0000000000000000E+00
     particles per unit mass  2.7803267189421238E-01   level radius [Rp]  1.0000000000000000E+00
   solved ghost counts [n0]: heavy  8.7942506135275134E-01   electrons  1.3117459208170167E-08
   ghost H2 partition closed to  0.00000E+00 in at most 4 passes
 [base boundary] evaluate route: face states recomputed at a read because
 the composition had moved since the boundary was derived: 0
```

The closure residual is exactly zero because the iteration reaches a fixed
point bitwise: once the returned `x_ion` repeats to the last bit, the imposed
partition repeats and the solve returns the same root, so
`x_H2 - x2 (1 - x_ion)` is the same expression on both sides. On the
`carrier_model_a_newton` fixture of the suite it is 4.65605E-12 in 2 passes.

**The prescribed count and the solved count are 4.66 per cent apart on this
state** (8.4026e-01 against 8.7943e-01), which is D5a section 1's measurement
of the same two numbers. They now have two names, only the prescribed one is
serialized, and nothing feeds the solved one back: D5a section 6.2 measured
that doing so takes the cell-1 row to 1.0.

**The failure path was exercised.** MEASURED on a throwaway build of the same
text with the tolerance set to zero and the passes to three, so that the
closure cannot succeed:

```
 (ioniz_eq) STOP: the base handoff partition of a lower ghost did not close
 against its own ionization balance
   ghost cell 0   passes 3
   last move of the imposed 2 n(H2)/n_H  1.25218E-08   tolerance  0.00000E+00
   residual of x_H2 - x2 (1 - x_ion) at the returned state  1.47044E-10
ERROR STOP ioniz_eq: base ghost H2 closure failed
```

An unclosed pair is a refused boundary and not a fallback.

---

## 7. What replaces the lower ghost rows of a restart

`load_IC` writes the composition of the gas the reservoir holds at the base
level into both lower ghost rows (`base_reservoir_composition_row`). Two
statements make it, and both are declared inputs of the boundary.

**The elemental abundances are the reservoir's**: He/H and each trace
element's El/H, the numbers the input states and the restart metadata block
carries as its reservoir field. That is the rule this loader already applies
to the base cells of a diffused state, where the column keeps its own
separation and cells `1-Ng` to 1 are carried onto the input's He/H. MEASURED,
it reproduces what the archived states carry: the ghost of the two molecular
reference states reads He/H = 9.7000000000 and 2.1300000000, the input's, and
the hydrogen nuclei per unit mass agree with the control's ghost to every
printed digit.

**The partition within each element is the first physical cell's**: the
ionization stages, the molecular ions and the oxygen carriers as fractions of
their own element's nuclei. The lower-atmosphere model states no ionization,
the ghost's own balance does and the sweep solves it, so what is needed here
is a starting point in the right basin and the cell half a spacing above the
ghost is one. **MEASURED with a neutral seed instead**: on the hot WASP-121b
base of the `grid_and_gates` restart fixture (T = 2358 K at the ghost) the
ghost's ionization solve returned a non-root above the amnesty cap at BOTH
lower cells, both cells kept the neutral seed, and the state's certification
was refused on a chemistry entry that had nothing to do with the boundary.
With the first cell's partition the same run has no non-root and certifies.

The molecular partition is the one exception: the handoff states x2 for the
inflowing gas and the sweep closes `x_H2 = x2 (1 - x_ion)` against the ghost's
own ionization. **The seed takes that form too**, and it must: at the
element-ratio ceiling x2 = 1 and a seed that gave x2 of EVERY nucleus would
ask for hydrogen the cell does not have. MEASURED before that was written that
way, on the helium-rich `HeH9.7` state whose x2 is exactly 1: the seed's own
hydrogen budget was over-subscribed by the first cell's ionized fraction,
5.5e-06 of the nuclei, and the cell-1 continuity row read 9.6e-05 against the
8.2e-08 it reads with the budget closed.

The thermodynamic row of the two ghost cells is the first physical cell's, a
placeholder the first `Apply_BC` overwrites with the ghost cell averages of
the boundary's own hydrostatic isentrope.

The upper ghost rows are still read: they are the free-outflow continuation's
and the composition sweep starts its columns there. The writer still writes
all `N + 2*Ng` rows, so the file layout and every reader of it are unchanged.

**One further path by which the file's ghost rows reached the state was found
and closed.** MEASURED: with the lower ghost rows replaced but the loader's
"did the composition move" test still reading the ghosts' own He/H rescale
factors (`load_IC.f90`, the `restart_composition_move_tol` test), scaling the
H2 column of those two rows by 1 + 1e-6 carried the ghosts' factors past the
tolerance, turned the flag on and sent the loader's density branch the other
way, moving the physical column at 1e-13 and the cell-1 row from 8.86e-09 to
8.26e-09 (measured on the build this was found on, whose ghost seed was an
earlier form; what the measurement shows is the path, not the value). The test
now reads the physical cells alone, which is what the
block's own comment above it already states about the ghosts.

---

## 8. The face cache and its validity

`Apply_BC` records the input set it derived the boundary from. `Rec_BC`, which
is handed reconstructed face states and not cell averages, checks the half it
can reach without an argument and, where the composition has moved since,
derives the three face states again from the interior primitive array the
boundary was derived from and the composition installed now, which is the same
arithmetic `Apply_BC` performs. The count of such reads is reported at the end
of a marching run, of the evaluate route and of a stationary restart; zero is
the call order being right.

**The marching loop needed one derivation of its own.** MEASURED: with the
check in and nothing else changed, a molecular marching run
(`mol_base_handoff`) recomputed at the first reconstruction of its SECOND
step, because the loop's closing composition refresh (the
`get_species_densities` after the adoption boundary) moves the density and the
composition the ghost continuation and the caloric map are evaluated at AFTER
the step's last `Apply_BC`. The boundary is therefore derived at the top of
each attempt, from the state that attempt begins from; it writes the ghosts
and returns the interior bit for bit, so no attempt's starting point moves,
and on an atomic mixture it reproduces the standing boundary exactly. With it
no reconstruction of a step recomputes.

**What still recomputes, and where.** MEASURED: 19 reads on
`mol_base_handoff` and 21 on `mol_carrier` over 12000 steps, none on
`wasp_full` or `lower_profile`, and the first of them is at step 500, at the
PERIODIC CERTIFICATION of the marching state, one per such checkpoint. That
site reads the boundary of the composition the step had before its closing
refresh, and the recompute repairs the face state it reads; the ghost CELL
AVERAGES it also stands on are not repaired, because a read cannot write them.
D5b-1 section 8 recorded `report_marching_stop` as NOT such a site, on the
grounds that the loop ends its step with `Apply_BC` after the composition
update; the closing `get_species_densities` comes later still, so it is one.
Making that site derive the boundary itself is the remaining piece and is not
in this item's list.

---

## 9. The movement on the regression cases

MEASURED, control build against the isolated measured build, single threaded,
on scratch copies; no regression directory was written to and no golden was
refreshed. `mol_base_handoff`, `mol_carrier` and `lower_profile` at their
pinned `maxsteps` of 12000, `wasp_full` at `EXHALE_MAXSTEPS=400`. Data rows
are the non-comment lines of `Hydro_ioniz.txt` and `Ion_species.txt`; the
largest relative change is over the physical cells.

| case | data rows | largest relative change over the physical cells | `Mdot` | `du` (flux spread) | certification | face states recomputed at a read |
|---|---|---|---|---|---|---|
| `wasp_full`, 400 steps (atomic + metals, cold start) | **IDENTICAL**, every value of both files including the ghost rows | 0 | 7.20 both | 3.782838414699551E-01 both | relaxation snapshot, no claim, both | 0 |
| `lower_profile`, 12000 steps (profile handoff, He and metal diffusion) | **IDENTICAL** | 0 | 8.82 both | 3.166556371655431E+00 both | relaxation snapshot, no claim, both | 0 |
| `mol_base_handoff`, 12000 steps | move | 6.19e-09 (`cool`, cell 491) and 1.05e-08 (`H3p`, cell 207) | 10.57 both | 7.980620893678548E-01 -> 7.980620899766524E-01, 7.6e-10 relative | relaxation snapshot, no claim, both | 19 |
| `mol_carrier`, 12000 steps | move | 1.33e-11 (`v`, cell 3) and 1.29e-12 (`HeI`, cell 502) | 10.57 both | 7.963511276718442E-01 -> ...16818E-01, 2.0e-13 relative | relaxation snapshot, no claim, both | 21 |

The two atomic cases are bit for bit what they were: an atomic mixture's ghost
continuation carries no composition, so the cache fingerprint never moves and
the ghost closure has no subject. The two molecular marching cases move at
1e-8 and 1e-11, and they are cold starts, so the loader's ghost replacement
plays no part in them: what moves them is the ghost's own H2 closure and the
face states recomputed at the reads counted in the last column.

**Where those recomputes are.** MEASURED, the first of them is at step 500 of
`mol_base_handoff`, at the periodic certification of the marching state, and
there are 19 of them in 12000 steps, one per such checkpoint. The marching
loop refreshes the composition at the END of its step, after the step's last
`Apply_BC`, so the certification taken there reads a face state of the
composition before that refresh; the recompute repairs the face state it
reads. D5b-1 section 8 recorded `report_marching_stop` as NOT such a site
because the loop ends its step with `Apply_BC` after the composition update;
the closing `get_species_densities` after the adoption boundary makes it one
after all. The ghost CELL AVERAGES of that state are still the pre-refresh
ones, since a read cannot write them, so making that site derive the boundary
itself is the remaining piece and is not in this item's list.

---

## 10. Tests

`src/tests/boundary_state/`, extended with
`boundary_from_the_declared_inputs.sh` (eleven assertions, one or more for
each of the item's six parts) beside D5b-1's
`boundary_of_the_installed_composition.sh`.

| | measured |
|---|---|
| GREEN, the delivered build | 13 of 13 assertions PASS (2 of D5b-1's, 11 of this item's) |
| RED, a build of the entry text | 11 of 11 of this item's assertions FAIL |

The RED rows fail for the reasons the item names: the two runs that differ
only in the file's lower ghost rows write different data and different base
continuity rows; there is no closure line, no reservoir line on the state
file, no count of cache reads and no boundary-model identity; and the trace
carries no gap at either of the two evaluation sites of a solve, because the
sites do not derive a boundary there.

Suites run whole, MEASURED with `EXHALE_OBJDIR`/`EXHALE_EXE` pointing at the
isolated build: `boundary_state` (13 assertions), `carrier_retry` (156),
`carrier_constraint_attribution` (10), `carrier_reference_scales` (16),
`carrier_boundary_jacobian` (12), `attempted_step` (70),
`carrier_returned_state_acceptance` (36), `carrier_helium_inventory` (15),
`coupled_source_step` (30) and `fuv_band_ledger` (24) pass whole, 0 failing.
`grid_and_gates` reads 279 PASS and **4 FAIL**, against 284 PASS and 0 FAIL on
the control build. All four are consequences of this item, in a suite this
item does not own; they are set out in section 11 and are NOT changed here.

---

## 11. Four assertions of `grid_and_gates` that this item moved

**Closed by D5b-3** (the same day): the three round-trip rows now compare the
physical cells and read the lower ghost rows as the boundary output they are,
a new row asserts that those rows ARE the boundary a reader derives from the
physical column of the same file (zero on the delivered build, 4.69e-16 on
a build that reads them; CORRECTED 2026-09-19: the zero was a coincidence of
equal seeds, since the ghost's solve closes to base_ghost_closure_tol and its
last bits follow the seed the loader takes from the file's first cell; after
D2b the two seeds differ and the row reads 4.0e-16, so it now asserts the
round-trip rounding band 1e-14, and whether the file's ghost rows are read is
decided by the four ghost variants of `src/tests/boundary_state/`), and the stagnation fixture's pass budget
is raised from 50 to 100, the ending announcing at pass 89 under the derived
boundary. `grid_and_gates` reads 285 PASS / 0 FAIL. The fourth
closure-stage-mixing site of section 12 is closed there too: the periodic
certification derives the boundary itself and the recompute counts of
`mol_base_handoff` and `mol_carrier` fall from 19 and 21 to 0. What follows is
the record of the four as this item left them.

## 11.1 The four, as D5b-2 left them

MEASURED, the same suite on the control build and on this item's: 284 PASS /
0 FAIL becomes 279 PASS / 4 FAIL. The four are understood, they are
consequences the item intends, and `src/tests/grid_and_gates/` is not a file
this item owns, so they are reported and left alone.

| assertion | measured | why |
|---|---|---|
| `restart_stored_scalar_fields_round_trip` | 1.472811e-02 at row 1 | the round trip writes a state, loads it and rewrites it, and compares the two files ROW BY ROW. Every difference of the first trip is at row 1 or row 2, the two LOWER GHOSTS, which this item makes derived rather than stored: `v` 8.1e-03 at row 2, `p` 1.5e-02 at row 1, `T` 4.2e-03 at row 1. The physical rows are exact, and the SECOND trip is exactly zero in every column, i.e. the round trip is exact once the state's ghost rows are the derived ones |
| `restart_state_and_inventories_within_budget` | 6.422209e-02 | the same rows, through the inventories: `nH` and `nHe` 1.06e-02 and the charge 1.53e-02, all at row 1; `He/H` 1.8e-14 and `x_H2` 3.9e-10 |
| `restart_density_is_the_mass_policy` | 1.055125e-02 | the same rows again: the writer state's own mass closure is 4.6e-16 at row 110, so it is the ghost rebuild and not a serialization loss |
| `outer_iteration_ending_is_the_stagnation_one` | `pass_budget` against `no_progress` | the fixture drives an alternation that stagnates and ends on the no-progress rule at outer pass 44 on the control. With the boundary derived at the joint test (item 6), the joint distance the rule watches keeps moving by a little each pass and the solve runs to the pass budget of 50 instead. The gated species row is 7.40E-02 of 1.0E-05 at cell 290 in both, so what changed is the ending reason and not the state's convergence; the suite's own message says the configuration no longer reaches the ending under test |

What the first three ask for is a comparison over the physical cells, with the
ghost rows read as what they now are: the output of the boundary-state
operation, which a reload re-derives. The fourth asks for a fixture that still
reaches the stagnation ending.

---

## 12. Noticed outside the scope, reported and not changed

- The base continuity row of a molecular state is conditioned at about 4e8 to
  1e9 rounding floors per relative unit of the ghost composition (section 5),
  so its certification entry states as much about the ghost prescription as
  about the state. A gate of ten floors on that row is a gate on the
  prescription; the other entries of the inventory are unaffected.
- The three items D5a reported outside its scope stand: the comment at
  `ionization_equilibrium.f90` around the reservoir-count block defended a
  base pressure boundary condition `base_boundary` replaced (this item removed
  that comment with the write-back it justified); `write_output.f90` writes
  the state files with list-directed output and no declared round-trip format;
  and `.L26/fid_resolve` and `atomic_scalar_gj1132_kzz1e9/HeH2.13` are the
  same state pair.
- `EXHALE_resolved.out` prints `ntot_bc` and `dp_bc` at the end of a run. They
  are now the PRESCRIBED values for the whole run, where before the sweep had
  overwritten them with the solved ghost's counts; the file's own labels do
  not say which, and `write_setup_report.f90` is not this item's file.
- **The ghost's ionization is defined only to the sweep's acceptance band**,
  and the base continuity row amplifies that band to tens of rounding floors
  (section 5). Closing the ghost further means closing the loop the sweep does
  not iterate: the composition sets the particle count, the particle count and
  the boundary's pressure set the temperature the balance is solved at, and
  the balance sets the composition. This item closes the molecular partition
  against the ionization at a fixed temperature, which is what D5b-2 asks for;
  the temperature leg is the next increment and would need the sweep to
  recompute a ghost's T inside its own passes.
- **The suites write their work directories under `build/tests/<suite>/` when
  a driver does not read `EXHALE_TEST_OUT`.** Every suite here was run with
  that variable pointing at a scratch directory and most honored it, but some
  programs of `grid_and_gates` ran in the tree's `build/tests/grid_and_gates/`
  all the same. No object and no `.mod` in `build/` was written; it is the
  default those drivers carry.

---

## 13. How to repeat any of it

The two builds, every run directory and the analysis scripts are under
`/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/59a81fec-9ef0-4992-9a57-6494f229e6bb/scratchpad/D5b2/`
(`tree_ctl/` control, `tree_m/` the isolated measured text, `tree_fail/` the
closure-failure probe, `runs/` the evaluations, `reg/` the regression copies,
`mkrun.sh`, `runone.sh`, `ladder.sh`, `regcase.sh`, `cmp.py`).

- one reading: copy `input.inp`, `base.inp` and the `_IC` pair into a scratch
  directory, set `Restart intent: stationary evaluate`, run with
  `EXHALE_MASS_FLOOR_SCAN=1`.
- the boundary of the evaluation and its closure: the run prints them
  unconditionally; `EXHALE_BOUNDARY_TRACE=1` adds the gap and the signed base
  mass rows in `output/boundary_trace.txt`.
- the ghost variants: `python3 ghostedit.py <src Ion_species> <dst> scaleH2
  <eps>` / `copycell1` / `fromfile <other Ion_species>` (the script is D5a's).
- the suite: `EXHALE_EXE=<binary> EXHALE_TEST_OUT=<scratch>
  src/tests/boundary_state/run.sh`.
