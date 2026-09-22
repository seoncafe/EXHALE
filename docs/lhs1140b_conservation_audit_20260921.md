# LHS 1140 b: the discrete conservation audit of two checkpoints

2026-09-21. Phase 5 of `PLAN_20260920_rev9.md` section 15, built to the
contract of that document's sections 13.2 to 13.6. Design 1 of section 13.2
was chosen: a diagnostic export at the hydrodynamic assembly point, off by
default, with an independent consumer that shares no arithmetic with the
assembly.

## 1. The result, first

**The assembly is consistent to the bit.** On both checkpoints, and in all
four momentum branches, every one of the 500 physical cells rebuilt from the
exported faces, geometry, potential and sources gives back the row the
assembly returned with a defect of exactly zero: signed sum zero, sum of
magnitudes zero, largest local defect zero. The same holds for the three
flux differences separately and for the gravitational work the energy flux
difference carries inside it. Check 1 of section 13.3 closes at the
arithmetic floor, and the floor is quoted in section 7.

**The identity of section 13.4 closes at the arithmetic floor.** The
combined statement

```
R_energy + phi_cell R_mass
   = [A_R (F_E,R + phi_R F_m,R) - A_L (F_E,L + phi_L F_m,L)]/V - Q
```

is a different grouping of the same terms, so it is not bitwise. Its largest
local defect is 2.78E-17 on the atomic checkpoint and 1.72E-17 on the
molecular one, in code units of a row whose own terms reach 1E-03 and 1E-05,
and the largest defect measured in last bits of the largest term the identity
differences is 23 and 91.

**Check 2, stationarity, is a separate verdict and it fails on the atomic
checkpoint and passes on two of three molecular rows**, which is what the
stored identity records already say. The audit's contribution is not that
verdict but its localization: on the atomic checkpoint the entire mass
imbalance lives in cells 1 and 2, and above r = 1.1989 the mass flux is
constant to one last bit over 284 cells.

**The defect is therefore a defect of the STATE, not of the assembly.** The
audit found no ledger error, no missing term and no double-counted gravity.

## 2. Classification

Every new number here is **Mode R** (rev9 section 3): the current binary
evaluating an immutable stored generation. No physical step was taken, no
state was published, nothing in `LHS1140b/models/` was written, and each run
was made on a copy of the stored generation in a scratch directory of its
own. Values taken from stored manifests and state headers are READ; values
this work produced are MEASURED; statements about the source are INSPECTED
with `file:line`.

**Both evaluations are Mode R, the current operator on an old state.** For
the atomic checkpoint the boundary model of the stored state
(`..._contact_upwind_v2`) is not the one this binary solves
(`..._ghost_fixed_point_seed_reservoir_row_v3`) and the boundary is rebuilt,
so the residual measured is not the residual the state was written with. For
the molecular checkpoint the stored state carries the `_v3` model the binary
solves. Both facts are READ from the state headers and are the same ones
`docs/lhs1140b_repeatability_20260921.md` section 3.3 records.

## 3. Identities

### 3.1 The binary

| | |
|---|---|
| path | `<scratch>/EXHALE_test.x`, a private build (`make OBJDIR=<scratch>/obj EXE=<scratch>/EXHALE_test.x -j8`) |
| md5 MEASURED | `79057acb816d56a69b21292d4897255f` |
| repository HEAD | `93eed8667483`, working tree dirty (MEASURED, `git rev-parse`, `git status --porcelain`: 50 entries) |
| source digest | `8366b6de4ee8a0b4681b75bc147c592b`, the md5 of the sorted md5 list of `Makefile` and the 245 `.f90`/`.inc` files of `src/` (MEASURED 2026-09-21) |
| compiler strings | `GCC: (conda-forge gcc 16.2.0-5) 16.2.0`; `GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20)` (MEASURED, `strings -a`) |

Every measurement below was made with that binary. One comment block inside
`conservation_budget.f90` was reordered after it was built and before this
document was written; no executable statement changed, and the tree as
delivered rebuilds to a different md5 for that reason alone (an I/O
statement carries its own source line, so a shifted comment moves the
fingerprint).

This is a live tree with concurrent workers, so the source digest records
what stood beside the binary at the time of these runs.

### 3.2 The two checkpoints

| | `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | `molecular_scalar_gj1132_kzz1e9/HeH0.083` |
|---|---|---|
| generation | `g0002_20260919T004843Z_f7485b14` | `g0004_20260920T053734Z_79b42a03` |
| Stage A0 record | `docs/lhs1140b_identity_records_20260921.md` section 3.1 | the same document, section 3.4 |
| `Hydro_ioniz.txt` md5 MEASURED | `25e11562ec3cd1f6d204bf36f8e503eb` | `fd882d8918052069f8ed229dd376f8c4` |
| `Ion_species.txt` md5 MEASURED | `78ac49d7e1ce90a3a4f271785400efda` | `87288b03ad88a92b7b8b8b6c2b1bdd66` |
| `input.inp` md5 MEASURED | `e43c6697786711fb916d3b8eb3c113aa` | `c3916719a26a2f5a7d1d7025a7cd71fa` |
| `base.inp` md5 MEASURED | (the case uses none) | `1e33d197ead9269661dc6bb0d88c0e26` |
| relation to the binary | the current operator on an old state, boundary rebuilt | the current operator on a state of its own boundary model |

Every md5 above equals the one the generation manifest records and the one
the Stage A0 record MEASURED.

### 3.3 How each evaluation was made

```
env -i PATH=/usr/bin:/bin HOME=$HOME OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    EXHALE_RESIDUAL=1 EXHALE_CONSERVATION_BUDGET=1 ./EXHALE_test.x
```

in a directory holding a copy of the case `input.inp` (and `base.inp` where
the case has one), the stored state pair under the names `load_IC` reads,
`output/Hydro_ioniz_IC.txt` and `output/Ion_species_IC.txt`, a copy of the
executable, and an empty `output/`. `EXHALE_RESIDUAL=1` assembles the
stationary residual of the loaded state once and stops
(`src/EXHALE_main.f90:1263-1341`), so exactly one assembly is exported.
Every invocation exited with status 0.

The consumer was then run from the repository root:

```
python3 EXHALE_conservation_budget.py <run dir>/output/conservation_budget_0001.txt
```

The reports are `docs/audit_20260905/conservation_20260921/report_*.txt`, and
`EXPORT_DATA_MD5.txt` beside them carries the md5 of the DATA ROWS of each
export (the header carries a wall-clock line and cannot be compared by md5).

## 4. The contract as implemented

### 4.1 Design 1, at the assembly point, off by default

The export is `src/modules/time_step/conservation_budget.f90`, called from
`assemble_residual` (`src/modules/time_step/steady_residual.f90`) right after
`store_row_terms`, where the three row arrays, the sources, the radiative
terms and the transport sources exist together and the face data the assembly
stored is still the face data of this state.

**The key is `EXHALE_CONSERVATION_BUDGET`**, a positive integer: the number
of further assemblies to export, `1` exporting the first one. It is read once
per process and counted down, so two assemblies of one run can never export
under two different settings. Absent, empty, `0` or unreadable, nothing is
written and no file is opened. It is documented at its read site, as
`EXHALE_RESID_QUAD` (`hydrodynamic_rows.f90:147-177`) and
`EXHALE_MASS_FLOOR_SCAN` (`steady_residual.f90:598-602`) are, and in
`README.md` where the other diagnostic keys and the post-processing tools
are listed.

Each armed assembly writes `output/conservation_budget_<nnnn>.txt` with
`<nnnn>` the index of that assembly in the process, so a route that
assembles more than once does not overwrite its own record.

### 4.2 The gravitational work of the energy row

INSPECTED `src/modules/time_step/RK_rhs.f90:273-275`:

```
dF3p    = dAp*Fp(1)*(Gphi_i(j) - Gphi_c(j))
        - dAm*Fm(1)*(Gphi_i(j-1) - Gphi_c(j))
dF(3,j) = (dAp*Fp(3) - dAm*Fm(3) + dF3p)/dV
```

and INSPECTED `src/modules/states/Source.f90:44-52`, `S(3) = 0` in every
branch: the well-balanced branch returns `S = 0` outright and the ordinary
branch sets `S(3) = 0.0` explicitly. The gravitational work therefore rides
on the MASS face flux inside the energy flux difference, and an export of the
energy faces and `S` alone omits it.

The export carries the face potentials `Gphi_i(j-1)`, `Gphi_i(j)`, the cell
potential `Gphi_c(j)`, and `dF3p/V` as its own column
`grav_work_over_volume`. The file states in its own header, and it is
repeated here: **it is already inside `dF_energy` with a POSITIVE sign**, so
no second gravitational work term is to be subtracted anywhere in a budget
built on this file.

### 4.3 The four momentum branches

INSPECTED `RK_rhs.f90:231-270`, `Source.f90:44-54` and
`Num_Fluxes.f90:527-551`. `A` is the face area, `V` the cell volume, `F` the
stored numerical momentum flux, `p` the face pressure, `q` the face pressure
measured from the cell's own hydrostatic equilibrium, `dr` the cell width.

| branch | `dF_momentum` | explicit `S_momentum` |
|---|---|---|
| `plm_ordinary` | `(A_R F_R - A_L F_L)/V`, and `F` carries the face pressure (`Phys_flux`, `Num_Fluxes.f90:550`) | the weight plus `(A_R - A_L) p_c/V` |
| `weno3_ordinary` | `(A_R F_R - A_L F_L)/V + (p_R - p_L)/dr` | the weight |
| `plm_well_balanced` | `[A_R (F_R + q_up,R) - A_L (F_L + q_dn,L)]/V` | zero |
| `weno3_well_balanced` | `(A_R F_R - A_L F_L)/V + (q_up,R - q_dn,L)/dr` | zero |

Where transport is active every branch additionally subtracts `Smom`.

The export names the branch in force in its header
(`momentum_branch=<name>`), and carries `dr`, both face pressures, both
equilibrium departures, the explicit source and the transport source, so all
four are reconstructible from one file format. The consumer implements the
four and selects on the header, never on the flag pair.

**All four branches were exercised and all four rebuild bitwise**, section 7.

### 4.4 What every exported quantity is

Stated in the file header and repeated here. All columns are in code units:
lengths in `R0`, mass density in `n0 mu`, velocity in `v0`, time in
`t_s = R0/v0`, pressure and energy density in `p0 = n0 mu v0^2`, potential in
`v0^2`. The five normalization constants and `T0` and `b0` are written into
the header at the same precision as the data, so a reader can convert
without the input file.

| column | kind | units |
|---|---|---|
| `j` | cell index over the padded range | - |
| `physical` | 1 for a physical cell, 0 for a ghost | - |
| `r_cell`, `r_face_lo`, `r_face_hi`, `dr` | geometry | `R0` |
| `area_lo`, `area_hi` | geometry, `r_edg^2`, solid-angle factor omitted | `R0^2` |
| `volume` | geometry, `(r_+^3 - r_-^3)/3`, solid-angle factor omitted | `R0^3` |
| `phi_face_lo`, `phi_face_hi`, `phi_cell` | potential | `v0^2` |
| `face_mass_lo`, `face_mass_hi` | flux density | `n0 mu v0` |
| `face_momentum_lo`, `face_momentum_hi` | flux density | `n0 mu v0^2` |
| `face_energy_lo`, `face_energy_hi` | flux density | `p0 v0` |
| `face_p_lo`, `face_p_hi` | face pressure | `p0` |
| `face_q_dn_lo`, `face_q_up_hi` | face pressure measured from the cell equilibrium | `p0` |
| `grav_work_over_volume` | contribution already divided by the cell volume | `p0/t_s` |
| `rho`, `momentum_density`, `energy_density` | the conserved state | `n0 mu`, `n0 mu v0`, `p0` |
| `dF_mass`, `dF_momentum`, `dF_energy` | contribution already divided by the cell volume | `n0 mu/t_s`, `n0 mu v0/t_s`, `p0/t_s` |
| `S_mass`, `S_momentum`, `S_energy` | the same | the same |
| `heat`, `cool` | the same | `p0/t_s` |
| `Smom`, `Sene` | the same, the viscous and conductive sources | `n0 mu v0/t_s`, `p0/t_s` |
| `R_mass`, `R_momentum`, `R_energy` | the assembled rows | as `dF` |
| `momentum_ram`, `momentum_pressure`, `momentum_gravity` | the production attribution of the momentum row, not inputs to the identity | `n0 mu v0/t_s` |
| `equilibrium_pressure_force` | the weight the well-balanced row is read against | `n0 mu v0/t_s` |

**The spherical normalization.** The assembly carries the face area as
`r_edg^2` and the cell volume as `(r_+^3 - r_-^3)/3`, INSPECTED
`RK_rhs.f90:211-217` and `define_grid.f90:356-392`: **the common solid-angle
factor 4 pi is OMITTED from both**, and the export carries the same two
numbers. It cancels in every row, so a reader that restores it on one side of
a balance is wrong by that factor. The carrier row dump is the cautionary
example the plan names, `diffusive_photochemistry.f90:5754-5756` and
`:5860-5865`: its four face entries already carry the area and the cell
volume, and a reader that applies geometry again is wrong by construction.
This export never applies geometry twice: the face columns are flux
densities and the row columns are already divided by the volume, and each
column says which it is.

**Ghost cells.** The lowest ghost cell has no lower face and carries no
equation, and its row is zero by construction (`RK_rhs.f90:160-162`), so it
is not exported. Every other cell of the padded range is, with the
`physical` flag naming the 500 that carry an equation. The consumer takes
every statistic over the physical cells; the ghosts are there so that a face
a physical cell shares with a ghost can be read from both sides.

**The equilibrium pressure departures.** They are written only under the
well-balanced option and are written as NaN otherwise, because the arrays
then hold whatever an earlier call left (INSPECTED `RK_rhs.f90:193-199`, both
assignments sit inside `if (well_balanced)`). A zero there would read as a
measured zero; a NaN cannot. MEASURED: the two columns are `NaN` in the
`weno3_ordinary` and `plm_ordinary` exports and finite in the two
well-balanced ones.

**The transport sources are combined at this point.** INSPECTED
`viscous_conduction.f90:609-613`, `Sene = vel*Smom + qv + Qc`, the viscous
work, the viscous dissipation and the conduction already added together, so
`Sene` is exported as the one number the assembly subtracts and the header
says so. Separating the three needs a capture inside that routine, which
this export does not make and which neither checkpoint would exercise:
MEASURED, both carry `visc=F cond=F` and both `Smom` and `Sene` are exactly
zero in every cell.

### 4.5 Precision

Every real is written `ES25.16E3`: seventeen significant decimal digits,
which round-trips a binary64. The geometry is at the same precision. This is
the precision of the DOUBLE the assembly returns and it is not the precision
the assembly worked in: the kind-generic rows can form the difference in a
wider kind and convert each returned array separately (rev9 section 13.6,
INSPECTED `hydrodynamic_rows.f90:208-213`,
`hydrodynamic_rows_body.inc:1504-1511`), so on that route a budget rebuilt
from these faces can differ from the exported `dF` by the rounding of that
conversion. Three precisions are therefore three separate facts and the
header states the first of them: the internal arithmetic precision (the
`rows_kind`), the returned-array precision (double, always), and the
serialization precision (seventeen digits). MEASURED here: both evaluations
report `rows_kind=production_double`, so the first and the second coincide
and the question does not arise for these numbers.

For contrast, the existing files the plan measures: the carrier writer prints
`ES14.7`, eight significant digits (`diffusive_photochemistry.f90:2829`,
`:2835`, `:2841`), and the loaded-state hydrodynamic profile `ES16.8`, nine
(`EXHALE_main.f90:1305`). The defects reported in section 7 are between 1E-17
and 1E-22 in code units against terms of order 1E-03; neither of those
formats can resolve them.

### 4.6 Which assembly produced the row

The header states, for the evaluation that produced the file: the actual
`rows_kind`, the actual `ieq_sweep_state_kind`, the effective reconstruction
from `assembled_reconstruction_is_plm()` (the single rule,
`parameters.f90:1454-1475`), the momentum branch that follows from it and the
balance option, `well_balanced`, `transport_active`, the flag pair
`use_plm`/`use_weno3` with `recon_lambda_on` and `recon_lambda`, and
`rec_method`. It says beside them that the flag pair does not name the branch
while the continuation is armed and that the branch line does.

MEASURED, both checkpoints:

```
# assembly rows_kind=production_double ieq_sweep_state_kind=marching reconstruction=WENO3 momentum_branch=weno3_well_balanced
# flags well_balanced=T transport_active=F use_plm=F use_weno3=T recon_lambda_on=F
```

Three consequences, each of which the plan asks to be recorded:

1. `ieq_sweep_state_kind=marching` means the assembly selector of section
   13.6 never fires: `assemble_residual` replaces `ROWS_PRODUCTION` only when
   the state kind is NOT the marching one (INSPECTED
   `steady_residual.f90:239-247`). The `EXHALE_RESIDUAL` route leaves it at
   marching, so **a generic comparison was NOT PERFORMED** and
   `EXHALE_RESID_QUAD` would have had no effect on either evaluation. This
   agrees with what `docs/lhs1140b_repeatability_20260921.md` section 9
   measured on the same route.
2. `recon_lambda_on=F`, so the endpoint restriction of section 13.6 is met
   and the state is not inside a reconstruction continuation: section 13.5
   does not apply to these two files, and the departure arrays hold this
   evaluation's own WENO3 departures.
3. The reconstruction is WENO3 although both cases state
   `Reconstruction scheme: PLM` or `PLM+WENO3` in `input.inp`, because every
   stationary route installs WENO3 first
   (`select_stationary_reconstruction`, `stationary_operator.f90:80-101`).
   The header reports what ran, not what the input asked for, which is the
   point of defect 10.5's fourth condition.

## 5. The consumer

`EXHALE_conservation_budget.py`, at the top level of the repository beside
the other analysis tools, Python 3 standard library only.

It reads the export and forms every row again from the face fluxes, the face
areas, the cell volume, the potential and the sources. It never sums the
exported residual and calls that a reconstruction: the exported `R` columns
are only ever the thing compared against. It selects the momentum expression
on the header's branch name, and refuses a branch name it does not know.

What it reports:

1. **Assembly consistency**, per row: the number of cells, the signed sum of
   the local defects, the sum of their magnitudes and the largest of them
   with its cell, because a signed sum alone hides opposing errors. The same
   three statistics for the three flux differences and for the gravitational
   work column.
2. **The arithmetic floor** of each rebuild: `eps` times the largest term the
   row differences, reported as a minimum, a median and a maximum over the
   column, which is the smallest defect a binary64 difference of those terms
   can resolve.
3. **The identity of section 13.4**, with `phi_cell R_mass` KEPT, because
   these checkpoints are not stationary and their mass rows are not small.
   It is reported BESIDE the gas-energy row of item 1 and never in place of
   it: one combined equation can hide compensating errors in the mass row and
   the energy row, and only the two together locate them.
4. **Internal-face agreement of the shared conservative fluxes**: the upper
   face of cell `j` against the lower face of cell `j+1`, for the area, the
   three fluxes and the potential. It deliberately does NOT cover the
   equilibrium pressure departures: INSPECTED `Num_Fluxes.f90:166-168`,
   `q_up = 0.5*(dev_L + dev_R + dp_eq)` and `q_dn = q_up - dp_eq`, so the two
   values at one face are measured from the equilibria of two different cells
   and are constructed to differ by `dp_eq`. Their failure to cancel is the
   intended arithmetic.
5. **Volume-weighted telescoping over subdomains**: the sum of `V R` over a
   window against the flux its two end faces carry, for the mass row and for
   the energy row augmented by the potential energy transport, both of which
   are exact sums of one conservative flux. The windows are the whole
   physical column, the column below and the column above the certification
   gate `r = 1.20` (`cert_regime_wind_r`, `certification.f90:389`), and the
   column up to each row's own failing cell. `--subdomain J` adds more.
6. **Stationarity**, separately: `|R|` divided by the largest term of that
   row, which is the certification's own measure (`residual_row_scale`,
   `steady_residual.f90:815-876`), against the certification's tolerances.
   It is computed from the exported terms alone, so the agreement with the
   measure the binary printed is itself a check on the reading of the scales.

What it does not do, and why. **It does not audit an elemental or carrier
budget.** The three rows this export carries are total mass, radial momentum
and total energy; no species appears in them, so the nucleus-count weighting
the plan asks for does not arise here: the mass row already carries every
nucleus with its own mass. An elemental budget, where H2 alone is not
conserved through a dissociation front, is the business of the carrier row
dump and is a separate export. Its absence here is a bound on this
deliverable, not a claim that it is unnecessary.

## 6. The two tolerances, kept apart

| | check 1, assembly consistency | check 2, stationarity |
|---|---|---|
| the question | does the row rebuilt from the exported terms equal the row the assembly returned? | does that row meet the physical acceptance criterion? |
| the tolerance | arithmetic and serialization: `eps` times the largest term the row differences, section 7.2 | the certification's: momentum 1E-08, energy 1E-06, mass `max(3E-12, the cell's rounding floor)` |
| source of the tolerance | `sys.float_info.epsilon` and the seventeen written digits | `cert_tol_momentum`, `cert_tol_energy`, `cert_tol_mass`, `cert_tol_mass_at`, `certification.f90:249-251`, `:893-923` |
| what a pass means | the export is complete and the reader's ledger is the assembly's | the state is stationary on that row at that tolerance |

A checkpoint that is not stationary passes the first and fails the second,
which is exactly what section 7 measures. And agreement between a writer and
its reader says nothing about whether the flux or the source formula is
physically right: both checks validate the assembly, not the physics of any
source model.

One limitation, stated where it is used. The continuity row's tolerance in
the binary is a function of the cell, `cert_tol_mass_at`, which raises the
fixed `3E-12` to the cell's own rounding floor and needs the caloric sound
speed. The consumer does not rebuild it and reports the fixed floor, so a
continuity refusal it prints is a LOWER BOUND: the identity record for the
atomic checkpoint READS the binary's own value as `1.9E-08` at cell 2, and
the measured row measure there is `9.806E-01`, seven decades above it, so the
refusal stands on either tolerance.

## 7. The measured defects

### 7.1 Check 1, assembly consistency

MEASURED, over the 500 physical cells of each file. "ulp" is the defect
divided by `eps` times the largest term the row differences.

| file | momentum branch | row | signed sum | sum of magnitudes | largest | max ulp |
|---|---|---|---|---|---|---|
| atomic `HeH9.7` | `weno3_well_balanced` | mass | 0 | 0 | 0 | 0 |
| | | momentum | 0 | 0 | 0 | 0 |
| | | energy | 0 | 0 | 0 | 0 |
| molecular `HeH0.083` | `weno3_well_balanced` | mass | 0 | 0 | 0 | 0 |
| | | momentum | 0 | 0 | 0 | 0 |
| | | energy | 0 | 0 | 0 | 0 |

and likewise, exactly zero in every entry, for the three flux differences
`dF_mass`, `dF_momentum`, `dF_energy` and for the gravitational work column
`grav_work_over_volume` against its rebuild from the areas, the mass faces
and the three potentials.

**The identity closes to the arithmetic floor, and here the floor is zero**:
the consumer performs the same operations on the same binary64 values in the
same order, so the difference is not merely small, it is the same number. The
statement that carries information is the one above it: the export carries
every term the row is built from, and the ledger the consumer applies, taken
from the source and written independently in another language, is the ledger
the assembly applies. Nothing is missing from the file and nothing in it is
counted twice.

The three branch-coverage files confirm the same for the other three
momentum branches (section 7.4).

### 7.2 The arithmetic floor

MEASURED, `eps` times the largest term each row differences, over the 500
physical cells:

| checkpoint | row | min | median | max |
|---|---|---|---|---|
| atomic | mass | 2.805E-27 | 9.904E-23 | 4.857E-19 |
| | momentum | 2.898E-29 | 1.693E-27 | 1.193E-16 |
| | energy | 2.046E-26 | 2.736E-21 | 2.293E-19 |
| molecular | mass | 6.097E-26 | 2.152E-21 | 1.284E-19 |
| | momentum | 7.020E-26 | 6.796E-25 | 1.969E-17 |
| | energy | 2.945E-25 | 3.800E-20 | 1.174E-18 |

This is the floor a budget reconstructed from this file can resolve. It is
what the seventeen digits of section 4.5 buy: at eight printed digits the
serialization alone would inject a defect of order 1E-08 of the largest term,
eleven decades above these floors, and any defect below it would be
unresolvable rather than absent.

### 7.3 The identity of section 13.4

MEASURED, over the 500 physical cells:

| checkpoint | signed sum | sum of magnitudes | largest, with cell | largest in ulp, with cell |
|---|---|---|---|---|
| atomic | +4.453807E-17 | 7.443264E-17 | 2.775558E-17 at cell 2 | 23.07 at cell 473 |
| molecular | -3.602485E-18 | 5.915983E-16 | 1.716682E-17 at cell 15 | 90.97 at cell 477 |

This is a genuine regrouping of the same terms: the identity collects the
mass flux under the potential instead of leaving it inside `dF3p`, so it is
not expected to be bitwise and it is not. Tens of last bits of the largest
term is the size a rearranged sum of terms of that ratio gives, and it is
eleven decades below the rows at the failing cells (1E-03 and 1E-05), so the
identity carries no information about the failure and adds none.

The gas-energy row check of section 7.1 stands beside it and not under it.

### 7.4 Branch coverage

Each of the four momentum branches was assembled and audited. The two
ordinary branches were reached by turning the balance option off on the
atomic checkpoint (with `Restart option change: wellbal`, so the change is
declared and the state is loaded rather than refused) and by two evaluations
of the `hp_front` regression case's input in a scratch directory, whose
marching stop and final certification both assemble under PLM.

| file | branch | assembly consistency, all three rows | identity 13.4, largest |
|---|---|---|---|
| atomic `HeH9.7` | `weno3_well_balanced` | 0 | 2.78E-17 |
| atomic, `Well balanced: False` | `weno3_ordinary` | 0 | 2.84E-14 |
| molecular `HeH0.083` | `weno3_well_balanced` | 0 | 1.72E-17 |
| `hp_front` input, marching stop | `plm_ordinary` | 0 | 1.16E-12 |
| `hp_front` input with `Well balanced: True` | `plm_well_balanced` | 0 | 1.25E-12 |

The last two are not LHS 1140 b states and are not offered as physical
results: they exist to exercise the two PLM branches of the export and the
consumer, because no stationary route can produce one
(`select_stationary_reconstruction` installs WENO3 before every stationary
evaluation). In the two ordinary branches the departure columns are NaN, as
the contract requires, and the consumer never reads them there. In
`plm_ordinary` the explicit momentum source is nonzero and carries both the
weight and the geometric pressure term, and the rebuild uses it; in the two
well-balanced branches `S_momentum` is exactly zero in every cell and the
rebuild adds no gravity of its own, which is the trap section 13.4 names.

### 7.5 Check 2, stationarity

MEASURED, the row measure `max_j |R_j| / s_j` computed by the consumer from
the exported terms alone, beside the value the binary printed in the same
run and the value the Stage A0 identity record READS from the stored
certificate.

| checkpoint | row | consumer, from the export | the binary, same run | the stored record | tolerance | verdict |
|---|---|---|---|---|---|---|
| atomic | mass | 9.806145E-01 at cell 2 | 9.806145E-01 | 9.806E-01 above 1.9E-08 at cell 2 | 1.9E-08 (the binary's) | REFUSES |
| | momentum | 8.723527E-06 at cell 1 | 8.723527E-06 | 8.724E-06 above 1.0E-08 at cell 1 | 1.0E-08 | REFUSES |
| | energy | 1.071147E+00 at cell 1 | 1.071147E+00 | 1.071E+00 above 1.0E-06 at cell 1 | 1.0E-06 | REFUSES |
| molecular | mass | 3.004611E-09 at cell 3 | 3.004611E-09 | (no hydrodynamic row refuses) | `max(3E-12, floor)` | see below |
| | momentum | 3.860743E-13 at cell 1 | 3.860743E-13 | (none) | 1.0E-08 | MEETS |
| | energy | 6.728749E-09 at cell 3 | 6.728749E-09 | (none) | 1.0E-06 | MEETS |

The consumer reproduces the binary's three row measures to every printed
digit on both checkpoints, which is the check that its reading of the three
row scales is the assembly's. The atomic checkpoint's three refusing rows are
the three the Stage A0 record reads from the stored certificate, so the
current operator on that old state refuses it where the historical record
says it was refused.

The molecular continuity row reads 3.0E-09 against the fixed floor 3E-12 and
the stored certificate names no refusing hydrodynamic row, so the binary's
cell-dependent tolerance there is above 3.0E-09. The consumer prints
`REFUSES` for it, which is the lower bound of section 6 and not a
contradiction of the certificate.

### 7.6 Subdomains, and where the atomic failure lives

MEASURED, volume-weighted sums, code units, atomic checkpoint:

| window | cells | sum of `V R_mass` | the two end faces | defect |
|---|---|---|---|---|
| the whole physical column | 1..500 | -4.1786256125E-07 | -4.1786256125E-07 | -2.65E-22 |
| below the certification gate `r = 1.20` | 1..216 | -4.1786256125E-07 | -4.1786256125E-07 | 0 |
| above the certification gate | 217..500 | +1.1689881853E-19 | +1.1689881853E-19 | +4.8E-35 |
| up to the failing cell 2 | 1..2 | -4.1786256125E-07 | -4.1786256125E-07 | 0 |

and the same rows for the energy budget augmented by the potential energy
transport close to 1.0E-19 or better on every window.

Two readings, and they are the finding of this audit.

1. **The telescoping closes.** The sum of the volume-weighted mass rows over
   any window equals the mass flux difference of that window's two end faces
   to 1E-22 or exactly. The interior faces cancel, which is the property a
   conservative finite-volume scheme is supposed to have, and the export
   preserves it.
2. **The whole mass imbalance of the state sits in cells 1 and 2.** The
   window 1..2 carries the same -4.1786256125E-07 as the whole column, and
   the window above the gate carries 1.17E-19, twelve decades smaller.
   MEASURED, the area-weighted mass flux `A F_m`:

   | face | atomic | molecular |
   |---|---|---|
   | base, `r = 1.000097` | +4.2300976616E-07 | +1.1186068062E-07 |
   | the certification gate, `r = 1.198946` | +5.1472049116E-09 | +1.1186068072E-07 |
   | the outer boundary, `r = 29.273418` | +5.1472049117E-09 | +1.1186068072E-07 |

   On the atomic checkpoint the mass entering the base face is 82 times the
   mass leaving at the gate and at the outer boundary: 98.8 per cent of it
   does not reach the wind and accumulates in the first two cells. Above the
   gate the flux is constant to one last bit over 284 cells. On the molecular
   checkpoint the same flux is constant from the base to the outer boundary
   to a relative 9E-10, which is why its continuity row reads 3E-09 and not
   1E+00.

   So the atomic checkpoint's refusal is a two-cell base phenomenon and its
   wind is stationary. This localizes the failure to a term and a place: the
   base mass inflow through the lowest face, against the wind the column
   carries. It is the same picture `docs/p44_base_sawtooth.md` describes for
   the base face flux, measured here from the operator's own faces rather
   than reconstructed by integrating the mass row.

## 8. What the key does not move

MEASURED. With `EXHALE_CONSERVATION_BUDGET` unset the export opens no file,
writes nothing and takes one `get_environment_variable` per process. The
private regression, run with the binary of section 3.1 and the key unset:

```
REGRESSION_EXE=<scratch>/EXHALE_test.x backup/regression/run_check.sh check \
    hp_front roundtrip mol_base_handoff
```

MEASURED, the complete output of that run:

```
[build] skipped: REGRESSION_EXE=<scratch>/EXHALE_test.x (79057acb816d)
[hp_front] running (OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=100) ...
            final: count=100  du= 6.6225E-01  dtu= 3.9604E-03
       PASS Hydro_ioniz.txt (data identical)
       PASS Ion_species.txt (data identical)
       PASS Hydro_ioniz_adv.txt (data identical)
       PASS Ion_species_adv.txt (data identical)
[roundtrip] running (OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=40) ...
            final: count=40  du= 2.5794E+00  dtu= 6.4339E-04
       PASS Hydro_ioniz.txt (data identical)
       PASS Ion_species.txt (data identical)
       PASS Hydro_ioniz_adv.txt (data identical)
       PASS Ion_species_adv.txt (data identical)
[mol_base_handoff] running (OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=12000) ...
            final: count=12000  du= 7.9925E-01  dtu= 6.0928E-03
       PASS Hydro_ioniz.txt (data identical)
       PASS Ion_species.txt (data identical)
       PASS Hydro_ioniz_adv.txt (data identical)
       PASS Ion_species_adv.txt (data identical)
==> REGRESSION PASS (all cases byte-identical)
```

The verdict line, verbatim:

```
==> REGRESSION PASS (all cases byte-identical)
```

Twelve files over three cases, byte-identical to their references. The check
does not touch `backup/regression/golden/`, and no golden was refreshed by
this work.

## 9. What this does not establish

- It does not say the discrete equations are the right equations. Both checks
  validate the assembly against its own terms, and two implementations can
  share one formula error.
- It does not cover the elemental and carrier rows (section 5).
- It does not cover the kind-generic assemblies: both evaluations ran
  `rows_kind=production_double`, because the `EXHALE_RESIDUAL` route leaves
  the state kind at marching (section 4.6 item 1). A budget on the quadruple
  route would need the working-kind faces captured before the conversion,
  which the wrapper does not pass (rev9 section 13.6).
- It does not cover a state inside a reconstruction continuation. The export
  records `recon_lambda_on` and `recon_lambda` so that such a file is
  recognizable, and the consumer's branch selection has no entry for a
  blended row, so an intermediate lambda is an unsupported diagnostic
  configuration and not a failure of the physical equations (rev9
  section 13.5, second admissible response).
- The molecular checkpoint's own refusing rows are the carrier H2 balance and
  the elemental He/H partition, READ from its Stage A0 record. Neither is a
  hydrodynamic row and neither is addressed here.
