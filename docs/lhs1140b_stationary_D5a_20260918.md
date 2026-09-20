# D5a: where the boundary of an evaluated state comes from, and which of it the residual reads

Item D5a of `docs/PLAN_20260918_rev2.md` (section D5, "the item", the six
measurements), opened by L34b and L37 and shaped by
`docs/PLAN_20260918_review.md` sections 5.2 to 5.6 and
`docs/PLAN_20260918_review1.md` sections 4.4 and 6.3. 2026-09-18.

This item is a DIAGNOSIS. No production source was changed. Everything below
was measured with a private build of a frozen source snapshot, on scratch
copies of the state pairs; no catalog case directory and nothing under
`backup/regression/` was written to.

Every number is labelled MEASURED (this item ran it) or READ (from a file, a
log or a source line).

---

## 0. Binaries, snapshot and provenance

The working tree is under edit by other items, so the measurement works on a
frozen copy of `Makefile` and `src/`.

| | value |
|---|---|
| snapshot taken | 2026-09-18, from the working tree at git `3c73905` plus 90 uncommitted files (`git diff --stat`: 12194 insertions, 1146 deletions) |
| snapshot source manifest (md5 of the md5 list of every `.f90`/`.inc`) | `913131b4f5cd50f71e6ef8daf0ce8903` |
| control build, the snapshot entry text | `44afa4848d143bd6570f7be86ee72e88` |
| measured build, the same text plus the default-off diagnostics of section 0.1 | `1bd647c287c3722aa27a1c51e577e9e1` |
| compiler | conda-forge gfortran 16.2.0, `-O3 -fopenmp`, OpenBLAS of the same prefix |
| every run | `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`, `Restart intent: stationary evaluate` |

**The diagnostics change nothing when they are off.** MEASURED: the measured
build with every key unset reproduces the control build's
`Hydro_ioniz.txt`, `Ion_species.txt`, `Hydro_ioniz_adv.txt`,
`Ion_species_adv.txt`, `Cooling_breakdown.txt`, `Heating_breakdown.txt` and
`Lyman_Werner.txt` byte for byte on the molecular reference state, the only
difference in any file being the run timestamp of the provenance line. With
the keys ON the same seven products are again byte-identical to the keys-off
run.

**The snapshot binary is not the binary the states were written on**
(`3146d11b` for the L34b pair, `c2e9c999` for the earlier ones, the L22
binary for `i3_alt`). MEASURED, it nevertheless reproduces L37's readings
exactly on five of the six states it shares with that memo: the cell-1 mass
row of one no-step evaluation is `1.2475E-08` (`HeH2.13` molecular),
`2.3398E-08` (`HeH9.7`), `1.2016E-09` (`HeH2.13` atomic), `5.0212E-10`
(`HeH0.55` atomic), each to all five printed digits of L37 section 3. The
exception is `.L22/i3_alt`, where this binary reads `1.6332E-08` against
L37's `1.4863E-08`; `1.6332E-08` is what L37 recorded there for its `R0`
route, so the two routes have converged in this tree. Where a same-operator
round trip is claimed below, the state was written by THIS binary (the
re-entry ladder of section 3).

Two entries of the brief are not what they were taken for, and are reported
rather than used: `.L26/fid_resolve/output` and
`atomic_scalar_gj1132_kzz1e9/HeH2.13/output` are the SAME state pair,
`Hydro_ioniz_IC.txt` bitwise identical (MEASURED, `cmp`), so the "three
atomic states" are two; and the atomic reference of the brief,
`atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13`, is not the state L37 named
`at213`, which is `atomic_scalar_gj1132_kzz1e9/HeH2.13`. The latter is used.

### 0.1 The diagnostics, and where they would have to live

All of it is in the snapshot only. If any of it is to be kept, the advisor
decides; each item is behind an environment key that defaults off.

| file | what was added |
|---|---|
| `src/modules/states/boundary_state_trace.f90` (new) | `EXHALE_BOUNDARY_TRACE`: one block per named point with the reservoir, the cached face states, the ghost and base cells' conserved state, caloric state and every species fraction; a signed base mass row writer; a reader for a second state pair |
| `src/EXHALE_main.f90` | five trace points in `stationary_state_of_the_loaded_restart`, two helper routines, and `EXHALE_TRACE_EXPERIMENT` (section 6) |
| `src/modules/states/stationary_operator.f90` | one trace point inside `stationary_face_mass_flux`, after its own `Apply_BC` |
| `Makefile` | one line in `SRC` |

`EXHALE_EVAL_STATE_DUMP` and `EXHALE_MASS_FLOOR_SCAN` already exist in the
tree and were used unchanged.

---

## 1. Verdict

**The evaluate route assembles the residual of a molecular state on a
boundary built from the composition the state was loaded with, while the
sources, the pressure map and the sound speeds of the same residual are
those of the composition one sweep later.** The two are not the same gas: the
ghost composition moves inside the sweep, and under the caloric equation of
state a moved ghost composition moves the base face mass flux. Rebuilding
the boundary after the sweep, which is the second `Apply_BC` the JFNK
residual already performs (READ, `steady_newton.f90` line 4247, with the
reason stated in the comment above it), takes the cell-1 mass row of the
molecular reference state from **15.98 rounding floors to 9.06**, of
`HeH9.7` from **31.24 to 3.09**, and of `.L22/i3_alt` from **20.93 to 0.08**,
all three inside the gate of ten. MEASURED. On the two atomic states the
same rebuild changes the row by 0 and by 0.35 floors.

**The rest is the ghost composition itself, which the file carries and the
loader projects.** With the physical cells and the declared inputs held
fixed and only the two inner ghost composition rows of the file replaced by
another admissible ghost, the cell-1 mass row of the molecular state reads
9.06, 15.98, 425 or 7.4e8 floors (MEASURED, section 4). The reconstructed
boundary is therefore not a function of the physical state and the declared
inputs alone.

**Serialization of the physical column is innocent.** MEASURED: the
conserved density, momentum and ENERGY and every species fraction of cells 1
to N come back through a write and a read at 0, 2.19e-16, 4.03e-16 and
6.07e-16 relative on the molecular state, and at 0, 0, 0 and 6.47e-16 on the
atomic one. The conserved energy is not carried by the file; it is rebuilt
from the pressure column through the caloric map, and that round trip costs
4.03e-16. Five to seven decades below what moves the base mass row.

**The reservoir count is two quantities under one name.** `ntot_bc + dp_bc`
at startup IS the prescribed reservoir: `set_base_reservoir` is called once
with it (READ, `input_read.f90` line 2122) and fixes `base_reservoir_p` and
`base_reservoir_nhat`, which every base face state is continued from. Inside
the sweep `ntot_bc` is overwritten by the SOLVED ghost's particle count
(`n_tot(1-Ng)/n0`) and `dp_bc` by the solved ghost's electron count, and
neither is fed back: MEASURED, on the molecular state `ntot_bc` moves from
8.402582811444896e-01 to 8.794250686178137e-01 (+4.66 per cent) and `dp_bc`
from 1.0e-10 to 1.316e-08 in one sweep, while `base_reservoir_p` and
`base_reservoir_nhat` keep their startup values for the whole run. Restating
the reservoir at the solved count, as an experiment, takes the cell-1 mass
row to 1.0e+00 (section 6): the two are NOT two estimates of one quantity,
and the solved count must not be written back as a reservoir pressure.

**No physical cell drifts because a diagnostic was requested**, and two
processes agree bitwise (section 7).

**What this does not settle.** A one-ulp perturbation is not the smallest
thing the cell-1 row of a molecular base responds to: from the third
re-entry onward the loaded state repeats to 1e-16 in the ghost and the
cached face state repeats bitwise, and the row still alternates between 7.14
and 9.06 floors with a period of two, carried by a base face mass flux that
moves by 1.72e-09 relative (section 3.2). The conditioning of that row, not
its boundary, is what remains.

---

## 2. Table 1. The boundary at each stage of the evaluate route

The four points, on `molecular_scalar_gj1132_kzz1e9/HeH2.13/output` and on
`atomic_scalar_gj1132_kzz1e9/HeH2.13/output`. A is after `Apply_BC` and the
derived rebuild; B is after `ioniz_eq` and the pressure refresh; C is
immediately before `assemble_residual`; D is inside
`stationary_face_mass_flux`, after its own `Apply_BC` on its copy.

MEASURED, relative difference against the previous stage, formed as
`|x - y| / max(|x|, |y|)`; "-" is bitwise equality. (Against the earlier
value as denominator the `ntot_bc` entry reads 4.66e-02 instead of
4.454e-02; both appear below and the convention is stated where they do.)

| quantity | A to B (the sweep) | B to C | C to D (the report rebuilds) |
|---|---|---|---|
| `base_face_W` rho, v, p | - | - | 1.354e-15, **4.350e-11**, 8.228e-16 |
| `base_ghost_W(:,0)` rho, v, p | - | - | 2.939e-16, 4.350e-11, 5.284e-16 |
| `base_face_lower_W` | - | - | 1.436e-16, 4.350e-11, 1.019e-15 |
| ghost 0 conserved `u(1), u(2), u(3)` | - | - | 2.9e-16, 4.350e-11, 9.899e-14 |
| ghost -1 conserved `u(1), u(2), u(3)` | - | - | 1.404e-15, 4.350e-11, **1.809e-11** |
| ghost 0 `nk_per_mass` | 1.545e-13 | - | - |
| ghost 0 `x_h2` | 9.685e-13 | - | - |
| ghost 0 `gamma_eff` | 3.689e-14 | - | - |
| ghost 0 `f(H I)` | **4.836e-08** | - | - |
| ghost 0 `f(H2)` | 8.140e-13 | - | - |
| ghost -1 `f(H I)` | **8.771e-06** | - | - |
| ghost -1 `f(H2)` | 1.477e-10 | - | - |
| cell 1 `n_part_cell1` | 1.193e-16 | - | - |
| `ntot_bc` | **4.454e-02** | - | - |
| `dp_bc` | **9.924e-01** (1.0e-10 to 1.316e-08) | - | - |
| `base_reservoir_p`, `base_reservoir_nhat` | - | - | - |

On the atomic state the ghost caloric columns are identically zero
(`caloric_mixture_active` false), `ntot_bc` does not move at all (the sweep
updates it only `if (thereis_mol)`, READ) and `dp_bc` moves from 1.0e-10 to
1.088e-06.

**Which of them the residual read, and which the report read.** READ from
the source and confirmed by the trace:

- `assemble_residual` at C reads the ghost conserved state and the cached
  `base_face_W` that `Apply_BC` wrote at **A**, that is, of the composition
  the state was LOADED with; `Rec_BC` inserts `base_face_W` straight into
  the left slot of the base face (`Apply_BC.f90` line 343). Its interior
  sources, its pressure map and its sound speeds, including `gamma_eff` at
  the ghost index 0, are those of the composition installed at **B**.
- `stationary_face_mass_flux` at D reads a boundary it rebuilt itself, on a
  copy of the conserved array, after installing the composition of B.
- The caller's conserved array keeps C's ghosts after the report returns
  (the copy is the report's), but the MODULE cache does not: MEASURED,
  `base_face_W` after the report is D's, not C's. A residual assembled after
  the report therefore stands on a face state that does not belong to the
  ghosts of the array it is given. Section 5 measures what that is worth.

The 4.350e-11 in the face velocity is the whole of what one sweep does to
the base face state of this molecular column; it is 4.4e-11 of the row scale
directly, 0.06 rounding floors, and is NOT the size of the defect. The
defect is what the residual does with the ghost STATE, and the ghost's own
conserved energy differs between the two stages by 1.809e-11 at ghost -1.

---

## 3. Table 2. Serialization alone, and the re-entry ladder

### 3.1 The physical column through a write and a read

MEASURED with `EXHALE_EVAL_STATE_DUMP`, comparing the block the evaluation
names `work_before_products` (the state it is about to write) of one run with
the block `loaded` (the state after the loader hands it to the residual) of
the next. Cells 1 to N, worst over the column.

| state | `u(1)` | `u(2)` | `u(3)` | `p` | `T` | species fraction |
|---|---|---|---|---|---|---|
| `HeH2.13` molecular | 0 | 2.193e-16 (cell 403) | 4.031e-16 (cell 493) | 4.023e-16 | 6.028e-16 | 6.073e-16 (cell 57, He II) |
| `HeH2.13` atomic | 0 | 0 | 0 | 0 | 5.849e-16 | 6.467e-16 (cell 17, He II) |

`u(3)` is not a column of the file: the writer emits `rho, v, p, T`
(`write_output.f90` lines 297 to 315, list-directed) and the conserved
energy is rebuilt by `energy_density_from_pressure` at the loaded
composition. Its 4.031e-16 on the molecular state is therefore the
representation of the pressure column PLUS the caloric round trip, and both
together are at the last bit. On the atomic state the same round trip is
exact, which is the contract L9 section 4 states: with no H2 the map carries
no composition.

The products do not move the work state: `work_before_products` against
`work_after_products` is 0 in every column (MEASURED).

**Representation error is thus bounded at 6.1e-16 relative and closure error
is everything else.** The two are told apart by the next table.

### 3.2 The re-entry ladder, and where it stops

Each evaluation handed the state the previous one wrote. Cell 1,
`|R_1|/s_1`, and the same in units of that cell's own rounding floor
7.8046e-10 (the gate is ten floors). MEASURED, `HeH2.13` molecular.

| reading | `|R_1|/s_1` | floors |
|---|---|---|
| E1, the archived state as the route runs today | 1.2475e-08 | 15.98 |
| E1 with the boundary rebuilt after the sweep (section 6) | 7.0724e-09 | 9.06 |
| E2 | 7.0724e-09 | 9.06 |
| E2 with the boundary rebuilt after the sweep | 5.5734e-09 | 7.14 |
| E3 | 5.5734e-09 | 7.14 |
| E4, E6 | 7.0724e-09 | 9.06 |
| E5, E7 | 5.5734e-09 | 7.14 |
| the run's own certification of this state (READ, `run.log` line 2352) | 5.213e-09 | 6.68 |

**One rebuild of the boundary does what one whole re-entry does**, exactly:
E1 with the rebuild is E2 to every printed digit, and E2 with the rebuild is
E3. That is the isolation this item was opened for.

From E3 the ladder does not settle; it alternates with a period of two
between 9.06 and 7.14 floors, 27 per cent apart. MEASURED, what alternates
is the base face mass flux, 5.7370548910887339e-07 against
5.7370549009580230e-07, 1.72e-09 relative; the inner face, between cells 1
and 2 and untouched by any boundary, alternates by 2.2e-10. The loaded
states of the two phases differ by 1.2e-16 to 1.3e-13 in the ghost trace
species, and their cached `base_face_W` are bitwise identical. So the
alternation is a property of the row's conditioning, not of the boundary,
and the rounding floor underestimates it: a state difference of order 1e-13
produces a reading 1.9 floors apart.

L37 reported this ladder settling at 5.5734e-09 by E3 on binary `3146d11b`;
on this binary it is a 2-cycle. Both phases are inside the gate.

---

## 4. Table 3. Boundary independence: the same physical cells, different ghosts

Physical cells, `input.inp`, `base.inp` and every declared input held fixed.
Only the two inner ghost rows of `Ion_species_IC.txt` were replaced. Every
run starts from the state E1 wrote, so the physical cells are bitwise the
same in all four. MEASURED.

| ghost composition supplied | base face mass flux `F(0)` | cell-1 row | floors |
|---|---|---|---|
| the state's own (control) | 5.7370548910886662e-07 | 7.0724e-09 | 9.06 |
| the ARCHIVED ghost rows of the L34b pair | 5.7370550203212283e-07 | 1.2475e-08 | 15.98 |
| `n(H2)` scaled by 1 + 1e-6 | 5.7370526541754155e-07 | 3.3171e-07 | 425.0 |
| the composition of cell 1 | -2.0064666943373528e-05 | 5.7941e-01 | 7.4e+08 |

The archived ghost reproduces the archived reading exactly, which is L37
section 4.3 met again on this binary. The slope in the H2 count is
`d(|R_1|/s_1)/d(eps) = 0.325`, L37's 0.30 to 0.33.

Giving the ghost the composition of the first physical cell, which is what
an "extrapolate the interior" rule would supply and is a perfectly
admissible composition, REVERSES the sign of the base face mass flux. The
boundary is not weakly dependent on the supplied ghost; it is determined by
it.

The same four on the atomic state (`n(H I)` scaled instead of `n(H2)`):

| ghost composition supplied | cell-1 row | floors |
|---|---|---|
| the state's own (control) | 1.1852e-09 | 1.50 |
| the archived ghost rows | 1.1852e-09 | 1.50 |
| `n(H I)` scaled by 1 + 1e-6 | 9.5494e-10 | 1.21 |
| the composition of cell 1 | 1.1852e-09 | 1.50 |

Three of the four are bitwise equal. In an atomic gas
`adiabatic_index_from_state` and `pressure_from_energy_density` return
`gamma_ad` and `p/(gamma_ad - 1)` with no composition in them (READ,
`caloric_eos.f90` lines 690 to 716), so the ghost composition has no path to
the base face flux at all. The 9.5494e-10 is L37's plateau value, met again.

### Why the ghost composition reaches the base face flux

READ, the two paths:

1. `base_boundary.f90` continues the reservoir isentrope with
   `continue_hydrostatic_isentrope(0, ...)` and takes its sound speed from
   `adiabatic_index_at_T(0, ...)` (lines 663 to 670, 802), both at the
   composition of ghost cell 0.
2. The Riemann solve at the base face takes its left sound speed from
   `adiabatic_index_from_state(jL = 0, ...)` (`Num_Fluxes.f90` line 145),
   again at the ghost-0 composition, while `Rec_BC` supplies the left
   primitive state from the cache.

Both read the ghost composition installed at the moment of the call, which
in the evaluate route is the LOADED one for the cache and the SWEPT one for
the flux.

---

## 5. Table 4. Boundary idempotence and installed module state

MEASURED, `HeH2.13` molecular.

| comparison | result |
|---|---|
| `Apply_BC` called twice on the state the first call left (A vs A2) | **bitwise identical** in every traced quantity |
| the conserved array restored and `Apply_BC` called again (A vs A3) | **bitwise identical** |
| a state installed through `stationary_face_mass_flux`, then a second state (`HeH9.7`), then the first again (E1 vs E3) | **bitwise identical** in every traced quantity |
| the residual re-assembled after `stationary_face_mass_flux` left its boundary installed | cell 1 only; see below |

`Apply_BC` is idempotent and leaves the interior alone, as its own comment
claims (READ, `Apply_BC.f90` lines 113 to 132), and a full install through
`stationary_face_mass_flux` leaves no state behind that a later install
cannot overwrite. **What it does leave behind is the cache.** Re-assembling
the residual of the SAME conserved array after the report has run moves one
cell:

| state | `|R_1|/s_1` before the report | after | floors before | floors after |
|---|---|---|---|---|
| `HeH2.13` molecular | 1.2475e-08 | 9.9891e-09 | 15.98 | 12.80 |
| `HeH9.7` molecular | 2.3398e-08 | 2.4733e-08 | 31.24 | 33.02 |
| `.L22/i3_alt` molecular | 1.6332e-08 | 1.6332e-08 | 20.93 | 20.93 |
| `HeH2.13` atomic | 1.2016e-09 | 1.2016e-09 | 1.52 | 1.52 |
| `HeH0.55` atomic | 5.0212e-10 | 7.0852e-10 | 0.75 | 1.06 |

Cells 2 to 12 are bitwise unchanged in every case (MEASURED). The base face
mass flux moves by 2.49e-09, 1.33e-09, 8.3e-14, 0 and 2.06e-10 relative.
So a diagnostic that installs its own boundary changes the next residual of
the base cell by up to 3.2 floors, and the sign of the change is not fixed.

This is the in-run signature of the two certifications that disagree. READ,
`molecular_scalar_gj1132_kzz1e9/HeH2.13/run.log`: the pass that produced
them printed "no composition update was taken: this pass ends the iteration,
so the state handed back is the state certified above", so the two
certifications are of one conserved state AND one composition, with the
carrier and element residual measurements and the flux-spread reports
between them.

---

## 6. Table 5. Same-state operator agreement, and the two experiments

### 6.1 What the in-run record allows

A run prints one number per row per certification: the maximum, its cell,
the verdict cell's value and the tolerance. It does not print the signed
residual, the two adjacent face mass fluxes or the scale of any cell. A
cell-by-cell comparison of signed dimensional quantities between the in-run
certification and the evaluate route is therefore **not possible from the
record**, and making it possible is a requirement on D5b and D8, not a
measurement this item could take. What follows compares the quantities the
record does carry, as absolute dimensional values over each cell's own
recorded scale, with the cell's own rounding floor stated.

| state | in-run (READ) | evaluate route (MEASURED) | evaluate with the boundary rebuilt (MEASURED) | floor |
|---|---|---|---|---|
| `HeH2.13` molecular | 5.213e-09 at cell 1; a second certification of the same state 1.206e-09 at cell 3 | 1.2475e-08 at cell 1 | 7.0724e-09 at cell 1 | 7.8046e-10 |
| `HeH9.7` molecular | 2.319e-09 at cell 1 | 2.3398e-08 at cell 1 | 2.3132e-09 at cell 1 | 7.4894e-10 |
| `.L22/i3_alt` molecular | 3.249e-09 at cell 1 | 1.6332e-08 at cell 1 | 6.3728e-11 at cell 1 | 7.8046e-10 |
| `HeH2.13` atomic | 1.189e-09 at cell 1 | 1.2016e-09 at cell 1 | 1.2016e-09 at cell 1 | 7.9044e-10 |
| `HeH0.55` atomic | 7.414e-10 at cell 1 | 5.0212e-10 at cell 1 | 7.3953e-10 at cell 1 | 6.7147e-10 |

Signed differences over the common scale, in floors of the cell's own floor:
the evaluate route as it runs today stands **+9.30, +28.14, +16.77, +0.02
and -0.36 floors** above the run's own reading on the five states; with the
boundary rebuilt after the sweep the same differences are **+2.38, -0.01,
-4.08, +0.02 and -0.00 floors**. On three of the five the rebuild brings the
two operators to within 0.02 floors of each other.

`.L22/i3_alt` is the state whose provenance is a different binary, and it is
the one whose rebuilt reading falls four floors BELOW the run's.

### 6.2 The two experiments

Both are `EXHALE_TRACE_EXPERIMENT`, default off, and both remove one
candidate cause from one evaluation. Neither is a proposed change.

`bc_after_sweep`: `Apply_BC` once more after the composition sweep and the
pressure refresh, before the residual is assembled. This is literally the
call `steady_newton.f90` line 4247 already makes inside the JFNK residual
for the same stated reason. MEASURED, the effect is table 5's third column;
on the atomic states it is 0 and +0.35 floors.

`resync_reservoir`: `set_base_reservoir` called again with the particle
count the sweep measured, then the boundary rebuilt. MEASURED, the cell-1
row goes to **1.0000e+00** on the molecular state and **8.1001e-01** on the
atomic one, from 1.2475e-08 and 1.2016e-09. The reservoir pressure moves by
4.45 per cent and the base ceases to be the base the state was solved at.

### 6.3 Which quantity `ntot_bc + dp_bc` is, at each stage

READ from the source and confirmed by the trace.

| stage | what the name holds | who reads it |
|---|---|---|
| startup, `input_read.f90` line 2117 to 2123 | the PRESCRIBED reservoir: `ntot_bc` is the base nuclei count per unit `n0` corrected for H2 binding (`comp_ntot_bc`), `dp_bc` is a placeholder 1e-10; the pair is handed to `set_base_reservoir` as `p_level` and `(ntot_bc + dp_bc)/rho_bc` as `nhat_level`, and to `n_part_cell1` as the cell-1 count before the first solve | `base_boundary`, through `base_reservoir_p`, `base_reservoir_T`, `base_reservoir_nhat`, for every face state and every ghost average, for the whole run |
| every sweep, `ionization_equilibrium.f90` lines 3167 to 3195 | the SOLVED ghost: `dp_bc` is the ghost's electron density and `ntot_bc` (molecular runs only) is `n_tot(1-Ng)/n0`, the ghost's own heavy-particle count | nothing in a restart evaluation. `set_base_reservoir` is not called again (MEASURED: `base_reservoir_p` and `base_reservoir_nhat` are bitwise unchanged from A to D). The only remaining readers are `set_IC` (a cold start) and `write_setup_report` / `write_resolved_config` |

MEASURED consequence on the molecular reference state: the reservoir's
particles per unit mass is `base_reservoir_nhat` =
2.780326718942124e-01 while the ghost's own is
`nk_per_mass(0)` = 2.780326798747016e-01, **2.87e-08 apart and never
reconciled**. The ghost pressure that `base_ghost_averages` builds is
`base_reservoir_nhat * rho * T` while the caloric map that turns that ghost
into a conserved energy and a sound speed uses `nk_per_mass`. At the
sensitivity of table 3 (0.325 of the row per relative unit) a 2.87e-08
inconsistency is of order twelve rounding floors of the cell-1 mass row,
which is the same order as the whole refusal; this is an order-of-magnitude
consistency, not an isolated attribution, because the two counts do not
enter the same expression.

The comment at `ionization_equilibrium.f90` lines 3176 to 3190 justifies the
sweep's update by "the ghost pressure is `(ntot_bc+dp_bc) T0` at the pinned
density", which was the base boundary condition before `base_boundary` was
written. It is no longer: the base pressure now comes from
`base_reservoir_p`, frozen at startup. The update the comment defends
therefore no longer reaches the boundary. Reported, not changed (this item
takes no production source change).

---

## 7. Table 6. Repeated evaluations and a fresh process

MEASURED, `HeH2.13` molecular.

| comparison | result |
|---|---|
| the measured build with every diagnostic key unset, against the control build | every product byte-identical but the provenance timestamp |
| the measured build with `EXHALE_BOUNDARY_TRACE`, `EXHALE_MASS_FLOOR_SCAN` and `EXHALE_EVAL_STATE_DUMP` on, against the same build with them off | every product byte-identical |
| two separate processes, same state, same keys | `boundary_trace.txt`, `eval_state_dump.txt`, `Hydro_ioniz.txt` and `Ion_species.txt` all byte-identical |

No physical cell moves because a diagnostic was requested, and the route is
deterministic at one thread across processes.

---

## 8. The causes, with their measured sizes

On the molecular reference state, whose cell-1 mass row refuses the state at
15.98 floors against a gate of ten:

| cause | measured size at cell 1 | where |
|---|---|---|
| **closure-stage mixing**: the residual reads a boundary built from the loaded composition and sources built from the swept one | **6.92 floors** of the 15.98 (15.98 to 9.06), and 28.15 of 31.24 on `HeH9.7`, 20.85 of 20.93 on `i3_alt`; 0 to 0.35 on the atomic states | `EXHALE_main.f90`, the evaluate route; the same order in the steady trial near line 6047 and in the `EXHALE_RESIDUAL` block near line 1253. NOT in the JFNK residual, the marching loop or the checkpoint restore, all of which already apply the boundary after the composition refresh |
| **entry-composition dependence of the ghost**: the file's ghost composition, projected onto the file's ghost density, is not the ghost the state's own sweep would produce, and one sweep does not close it | the remaining 9.06 floors; supplying the archived ghost instead of the state's own is worth 6.92 floors (9.06 to 15.98) on its own, and an admissible ghost as far from it as cell 1's composition is worth 7.4e8 floors | `load_IC.f90` (the projection), `ionization_equilibrium.f90` lines 1999 to 2010 (`x_h2_fix = q_H2,base (1 - x_ion,entry)`, with `x_ion` read from `nhii` and the MODULE array `nmol_eq` of the previous sweep) |
| **the reservoir count**: `base_reservoir_nhat` is 2.87e-08 from the ghost's own `nk_per_mass`, and `ntot_bc` is 4.45 per cent from `base_reservoir_p`, permanently | order twelve floors by the table-3 sensitivity, not isolated; feeding the solved count back as a reservoir pressure gives a row of 1.0 | `input_read.f90` line 2122 (once), `ionization_equilibrium.f90` line 3195 (every sweep), `base_boundary.f90` (the only reader) |
| **serialization** | <= 6.1e-16 relative on every physical conserved variable including the energy, and on every species fraction | none; the physical column is not the problem |
| **stale module state** | `Apply_BC` is idempotent and a full install through `stationary_face_mass_flux` leaves nothing that a later install cannot overwrite; what it leaves is the FACE CACHE, worth up to 3.2 floors to the next residual of the same array | `Apply_BC.f90` line 175 and 343 (the cache and its only reader), `stationary_operator.f90` lines 160 to 174 |
| **conditioning of the row itself** | a 2-cycle of 1.9 floors (27 per cent) from states differing by 1e-13, with the ghost and the cached face state bitwise identical | not a boundary defect; it bounds what any boundary operation can deliver |

The last row is the reason D5b cannot be judged by "the six states certify
on the first reading". Two evaluations of one molecular state agree on the
base mass row only to about two rounding floors, whatever the boundary does.

---

## 9. What D5b's boundary-state operation must take and produce

A proposal for the user's decision. It follows the requirements of
`docs/PLAN_20260918_rev2.md` D5b and adds what the measurements above
settle.

**Input**, and nothing else: the physical conserved state `u(1:3, 1:N)`, the
physical composition `f_sp(1:N, :)`, the PRESCRIBED reservoir (a pressure, a
temperature, a particles-per-unit-mass and the level radius, versioned), the
radiation context, and the model options. The ghost rows of a restart file
are NOT an input: measured, they are the difference between 9.06 floors and
7.4e8, and nothing in the physical state selects one of them.

**Output**: a ghost composition, a ghost conserved state, a caloric state
(`nk_per_mass`, `x_h2`, `molecular_cell`, and the `gamma_eff` and `c_s` that
follow) and the boundary face data `base_face_W`, `base_ghost_W`,
`base_face_lower_W`, all of ONE composition, tagged with the identity of the
state they were built from.

**What the measurements say it must satisfy.**

1. **One composition, one boundary.** The face cache, the ghost conserved
   state and the caloric arrays the flux reads must be outputs of the same
   call. Measured cost of not doing so: 6.92 to 28.15 floors. The immediate
   consequence for the present code is that the evaluate route, the steady
   trial and the `EXHALE_RESIDUAL` block must install the boundary AFTER the
   composition sweep, as the JFNK residual, the marching loop and the
   checkpoint restore already do. That single change is measurable now and
   is the largest term.
2. **The ghost composition is solved, not read.** The prescription
   `x_h2 = q_H2,base (1 - x_ion)` is a condition on the ghost coupled to the
   ghost's own ionization balance; today it is closed by alternation across
   sweeps with the entry composition supplying `x_ion` and the previous
   sweep's `nmol_eq` supplying the molecular ions. The operation must solve
   that pair to a stated tolerance and report its residual and its failure,
   or the ghost carries an iteration history that no amount of re-entry
   removes deterministically.
3. **The reservoir keeps its own name.** `ntot_bc + dp_bc` must be split:
   the prescribed reservoir (fixed, versioned, serialized with its meaning)
   and the solved ghost's counts (diagnostics of the state, never written
   back into `set_base_reservoir`). The measurement is unambiguous: writing
   the solved count back moves the cell-1 row to 1.0. Where the prescribed
   reservoir's particles-per-unit-mass and the ghost's own `nk_per_mass`
   disagree, that 2.87e-08 is a statement about the model and should be
   reported by the operation, not silently split between two expressions.
4. **A cache carries a validity.** `base_face_W` must know which installed
   state it belongs to, and a consumer that finds it stale must recompute it
   or refuse. Measured cost of not doing so: up to 3.2 floors, and the two
   in-run certifications of one state.
5. **The restart pair need not carry the ghost rows at all**, and should
   not carry them as degrees of freedom. It must carry the prescribed
   reservoir with its version, since that is the only boundary input the
   physical column cannot reconstruct.
6. **The acceptance test for D5b is not "the archived states certify".** It
   is: (a) the boundary is a function of the physical state and the declared
   inputs alone, checked by table 3's four ghost variants agreeing bitwise;
   (b) the residual and every report read a boundary of the same installed
   state, checked by table 4's re-assembly moving nothing; (c) one
   evaluation equals the re-entry ladder's first step, checked by E1 = E2;
   (d) the physical column still round-trips at 1e-16. The residual cell-1
   reading itself is reproducible only to about two floors on a molecular
   base, so it cannot be a bitwise acceptance criterion.

**The inflow/outflow composition rule is left exactly as it is** by all of
the above (review1 section 6.3): every measurement here is of the CURRENT
boundary model, and nothing proposed changes which characteristic carries
which datum. A physical change of that rule belongs to D2b.

---

## 10. Noticed outside the scope

- `.L26/fid_resolve/output` and `atomic_scalar_gj1132_kzz1e9/HeH2.13/output`
  are the same state pair (bitwise). Two of the brief's "three atomic
  states" are one.
- The comment at `ionization_equilibrium.f90` lines 3176 to 3190 describes a
  base pressure boundary condition (`p = (ntot_bc + dp_bc) T0` at a pinned
  density) that `base_boundary` replaced; the `ntot_bc` update it justifies
  no longer reaches the boundary. Not changed here.
- `steady_newton.f90` lines 4220 to 4246 already state this item's finding
  in full, for the JFNK residual, and fix it there. The evaluate route, the
  steady trial near `EXHALE_main.f90` line 6047 and the `EXHALE_RESIDUAL`
  block near line 1253 were not brought along.
- `write_output.f90` writes the state files with list-directed output
  (`write(2,*)`), not a declared round-trip format. It happens to round-trip
  the physical column at 1e-16 on gfortran here (measured), but the format
  is not stated anywhere and a different compiler's default width would
  change it silently.

## 11. How to repeat any of it

The snapshot, the two builds, every run directory and the analysis scripts
are under
`/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/59a81fec-9ef0-4992-9a57-6494f229e6bb/scratchpad/D5a/`
(`tree/` measured, `tree_ctl/` control, `runs/` the evaluations,
`an.py`, `cmp.py`, `ser.py`, `ghostedit.py`, `mkrun.sh`, `ladder.sh`).

- one reading: copy `input.inp`, `base.inp` and the `_IC` pair into a scratch
  directory, set `Restart intent: stationary evaluate`, run with
  `EXHALE_MASS_FLOOR_SCAN=1`.
- the four-point boundary trace: add `EXHALE_BOUNDARY_TRACE=1`; the blocks
  are written to `output/boundary_trace.txt`.
- the interleaving: add `EXHALE_TRACE_STATE_B=<dir holding a second pair>`.
- the two experiments: `EXHALE_TRACE_EXPERIMENT=bc_after_sweep` or
  `=resync_reservoir`.
- the ghost variants: `python3 ghostedit.py <src Ion_species> <dst> scaleH2
  <eps>` / `copycell1` / `fromfile <other Ion_species>`.
