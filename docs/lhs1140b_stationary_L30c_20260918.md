# L30c: why `oxygen_chemistry` moves by order unity, decided by measurement

Item L30c of `docs/PLAN_20260917.md`, the bounded diagnosis that follows
L30 (one spherical geometry for all transport) and L30b (one spelling of
that geometry). Written 2026-09-18 (KST). Every number is labeled MEASURED
(I ran it) or READ (from the source, a log or a document).

## Verdict

**Case (i): a discrete branch selection amplifying a rounding-level seed on
a snapshot that is not a solution.** The seed is not the carrier diffusion
geometry of L30 and not a growing divergence. It is the restatement of the
spherical shell volume that L30b carried into the hydrodynamic rows, whose
two spellings differ by at most **5.77e-13** relative on this case's grid
(MEASURED, below). The branch is the root ladder of the constrained
molecular ionization solve: on the **first ionization sweep of step 0** the
two builds are handed a state whose printed temperature agrees at all 43
cells both report, and they accept a different composition at **every one of
those 43 cells**, with the simplex violation differing at 39 of them and the
solver status flipping at one (MEASURED). From that point on there is no
growth phase: at a cap of ONE step every number of `Hydro_ioniz.txt` except
the radius already differs, and the maximum movement is 1.98, the same order
it has at 1000 steps.

The attribution is complete and bitwise, not statistical. A build of the
current tree with **only** the shell-volume restatement reverted to the entry
text reproduces the control's 1000-step snapshot **bit for bit in all four
output files**; a build with only the two carrier diffusion denominators
reverted reproduces the **measured** build's snapshot bit for bit
(MEASURED). The two experiments together place the whole of the movement on
the volume respelling and none of it on the transport geometry the item
L30 corrected.

**What a byte comparison of this case can mean.** `oxygen_chemistry` is a
1000-step relaxation snapshot of a cold isothermal start, stopped at
`du = 4.51` while the state still changes by tens of per cent per step
(READ, `backup/regression/oxygen_chemistry/README.md`; MEASURED here,
`du` rises monotonically from 2.62 at step 1 to 4.51 at step 1000 on both
builds). Its chemistry solve reports `no admissible root` at 102 cells on
one build and 983 on the other over the same 1000 steps (MEASURED). A
comparison of this case against a golden therefore certifies **that the same
binary reproduces itself**, and nothing else: it cannot separate a physics
change from a change in the last bits of one geometric factor, because the
root ladder converts either into an order-unity difference within one step.
It is not evidence about the physics of the oxygen network, and its movement
is not apportionable between the changes of a series, which is the same
reading the 2026-09-17 and the preceding golden refreshes recorded
(READ, `docs/Update_EXHALE_stage2.md` section 11). The integral the case
does pin is stable: `log10 Mdot` is 9.61 on both builds at 1000 steps, and
agrees to the two printed decimals at every cap measured (MEASURED).

No production source was changed by this item.

## What was compared

- **Control**: `git archive HEAD` (HEAD `3c73905ca8a7fe2af92a5c2b014c225b2feedff3`)
  unpacked to a scratch tree and built there with bare `make`. All **263**
  source files listed in `LHS1140b/models/BINARY_MANIFEST_59bfdb3fc4d0.txt`
  match that tree by `md5sum -c`, 263 of 263 OK (MEASURED), so the control
  binary is built from the source text of the binary the 2026-09-17 goldens
  were taken with.
- **Measured**: `make OBJDIR=build_L30c EXE=EXHALE_L30c.x` from the current
  tree, md5 `37f16ff6da99cd5b4b315437c1a70920` (MEASURED).
- Both run on scratch copies of `backup/regression/oxygen_chemistry/` with
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` and `EXHALE_MAXSTEPS` set by
  hand, never through `run_check.sh`.

The control reproduces the 2026-09-17 golden at the three cells the matrix
log names, exactly (MEASURED):

| file | row, column | control | measured |
|---|---|---|---|
| `Hydro_ioniz.txt` | 16, 6 (`heat`) | +1.076562e-03 | -1.242380e-03 |
| `Ion_species.txt` | 285, 40 (`H2O`) | 5.421331e+02 | 3.177630e-14 |
| `Hydro_ioniz_adv.txt` | 56, 3 (`v`) | +1.858263e+05 | -2.984002e+05 |

## The step-cap table

Maximum of `|a-b|/max(|a|,|b|)` over every number of the file, with the cell
and column at which it is attained. MEASURED.

| cap | file | max relative | row, column | control | measured | numbers above 1e-12 |
|---|---|---|---|---|---|---|
| 1 | `Hydro_ioniz.txt` | 1.9762e+00 | 49, 3 `v` | 5.253768e+02 | -5.381831e+02 | 3024 / 3528 |
| 1 | `Ion_species.txt` | 1.0000e+00 | 3, 6 `HeIII` | 5.592868e-07 | 0 | 10560 / 20664 |
| 5 | `Hydro_ioniz.txt` | 1.9629e+00 | 143, 3 `v` | -1.665157e+03 | 1.729271e+03 | 3024 / 3528 |
| 5 | `Ion_species.txt` | 1.0000e+00 | 52, 7 `HeITR` | 6.324484e-13 | 0 | 10560 / 20664 |
| 20 | `Hydro_ioniz.txt` | 1.9690e+00 | 234, 3 `v` | -7.896507e+03 | 7.652040e+03 | 3024 / 3528 |
| 20 | `Ion_species.txt` | 1.0000e+00 | 39, 7 `HeITR` | 3.940559e-13 | 0 | 10560 / 20664 |
| 100 | `Hydro_ioniz.txt` | 1.9536e+00 | 89, 6 `heat` | 1.773753e-04 | -1.860022e-04 | 3024 / 3528 |
| 100 | `Ion_species.txt` | 1.0000e+00 | 54, 7 `HeITR` | 1.112523e-12 | 0 | 10560 / 20664 |
| 300 | `Hydro_ioniz.txt` | 1.9076e+00 | 82, 6 `heat` | 1.276519e-04 | -1.158530e-04 | 3024 / 3528 |
| 300 | `Ion_species.txt` | 1.0000e+00 | 37, 3 `HII` | 1.198433e+04 | 0 | 10560 / 20664 |
| 1000 | `Hydro_ioniz.txt` | 1.8665e+00 | 16, 6 `heat` | 1.076562e-03 | -1.242380e-03 | 3024 / 3528 |
| 1000 | `Ion_species.txt` | 1.0000e+00 | 285, 40 `H2O` | 5.421331e+02 | 3.177630e-14 | 10560 / 20664 |

The count in the last column does not move with the cap. 3528 is 504 rows
times 7 columns and 3024 is 504 times 6: every column of `Hydro_ioniz.txt`
except the radius differs at every cell after ONE step, and no further
column joins them by step 1000. The same holds for the ion file, where the
numbers that stay equal are the species that are exact zeros on both builds.

The movement of the profiles themselves, measured against each column's own
largest magnitude rather than against the entry, tells the same story
(MEASURED):

| cap | rho | v | T | heat |
|---|---|---|---|---|
| 1 | 5.09e-02 | 2.67e-02 | 2.14e-01 | 9.99e-01 |
| 5 | 1.33e-01 | 7.04e-02 | 2.25e-01 | 9.95e-01 |
| 20 | 1.23e-01 | 8.26e-02 | 2.06e-01 | 9.94e-01 |
| 100 | 1.51e-01 | 1.66e-01 | 1.83e-01 | 9.88e-01 |
| 300 | 2.34e-01 | 3.94e-01 | 2.38e-01 | 8.03e-01 |
| 1000 | 1.30e-01 | 1.07e+00 | 7.27e-01 | 1.00e+00 |

A five per cent density and a twenty per cent temperature movement after one
step is not an amplified rounding difference; it is a different state.

## Where the first difference is, and that it jumps

The output files are written once, at the end of a run, so the earliest
observable is the marching log. MEASURED, from the two 1000-step runs:

| step | `du` control | `du` measured | relative |
|---|---|---|---|
| 1 | 2.621409866832086 | 2.623100409920833 | 6.445e-04 |
| 2 | 2.624737463182961 | 2.625796913355990 | 4.035e-04 |
| 5 | 2.634613037003779 | 2.633369203780526 | 4.721e-04 |
| 10 | 2.650696444779950 | 2.645744440804977 | 1.868e-03 |
| 100 | 2.890746768103546 | 2.836506986176454 | 1.876e-02 |
| 1000 | 4.507793195151416 | 4.388424972824752 | 2.648e-02 |

There is no growth by a fixed factor at each step to read: the first
marching measure the run prints is already apart by 6.4e-4, and 1000 steps later it is apart by
2.6e-2. A seed of 5.8e-13 does not reach 6.4e-4 in one step by any linear
response.

The jump is inside step 0, in the first ionization sweep. Both logs open
that sweep with the identical line
`the star-ward H2 column reaches 6.38E+21 cm^-2` and then diverge. Of the
cells the sweep reports, 64 appear on the control and 48 on the measured
build, 43 on both; at those 43 the printed temperature is identical at all
43, the accepted iterate's residual differs at all 43, the simplex violation
differs at 39, and the solver status flips at one (MEASURED). The first
rows:

| cell | T [K] | status ctl / msd | viol ctl | viol msd | res ctl | res msd |
|---|---|---|---|---|---|---|
| 150 | 1.448E+03 | 4 / 4 | 9.816E-06 | 9.816E-06 | 1.703E-06 | 2.703E-06 |
| 118 | 1.452E+03 | 4 / 1 | 2.181E-01 | 2.181E-01 | 3.647E-04 | 1.509E-06 |
| 115 | 1.453E+03 | 4 / 4 | 5.145E-02 | 3.641E-05 | 2.009E-03 | 3.688E-03 |
| 103 | 1.453E+03 | 4 / 4 | 2.812E-03 | 2.812E-03 | 4.424E-03 | 1.346E-06 |
| 102 | 1.453E+03 | 4 / 4 | 3.702E-06 | 3.029E-03 | 6.476E-03 | 4.921E-06 |

The sweep's own summary for the same step reads `11 root(s) outside the
simplex` on the control and `5` on the measured build, with
`every molecular equilibrium root left the physical simplex at 10 cell(s)`
against `5 cell(s)` (READ from the two logs). Over the whole 1000 steps the
ladder's outcomes stand in the ratios below (MEASURED, counts of log lines):

| message | control | measured |
|---|---|---|
| `no admissible root` | 102 | 983 |
| `stored state(s) rejected` | 105 | 986 |
| `constrained-continuation root accepted` | 632 | 1701 |
| `NON-ROOT accepted (relaxation amnesty)` | 160 | 174 |
| `NON-ROOT above the amnesty cap` | 6 | 7 |

Two builds whose only numerical difference is 5.8e-13 in one geometric
factor run the chemistry solve through a tenfold different number of root
rejections. That is the branch.

## The seed, measured

The two spellings of the shell volume are

    entry text   dV = (r_+^3 - r_-^3)/3      (also written (A_+ r_+ - A_- r_-)/3)
    L30b         dV = (r_+ - r_-)(r_+^2 + r_+ r_- + r_-^2)/3

They are the same number in exact arithmetic and differ in double precision
by the cancellation of the difference of cubes. Evaluated in IEEE double on
this case's own 500-cell Mixed grid, with `r_edg` formed from the run's own
radius column exactly as `define_grid.f90` line 150 forms it
(`r_edg(j) = 0.5 (r(j) + r(j+1))`, an exact operation), MEASURED:

- maximum relative difference **5.7707e-13**, at the cell whose faces are
  1.0085039078971381 and 1.00869940014035 R_p, i.e. in the 50 uniform base
  cells where `dr/r = 1.94e-4`;
- median over the column 1.1853e-14; nonzero at 99.4 per cent of the cells.

This confirms by measurement the error analysis L30 and L30b gave for it
(READ: "about 4e-13 relative", L30 report, item 1 of its outside-scope list).

## The three builds

| build | what it is | 1000-step snapshot |
|---|---|---|
| control | `git archive HEAD`, bare `make` | reference |
| measured | current tree, `OBJDIR=build_L30c` | max 1.87 against control |
| A | current tree, ONLY the two carrier diffusion denominators of `diffusive_photochemistry.f90` reverted to `r_j^2 max(dr_j, 1 cm)` and the two face areas to `r_edg^2` | **bit for bit the measured build** in all four files |
| C | current tree, ONLY the shell-volume restatement reverted: `species_face_flux.f90`, `RK_rhs.f90`, `Reconstruction.f90`, `steady_residual.f90`, `hydrodynamic_rows_body.inc` taken from HEAD | **bit for bit the control** in all four files |

Build C differs from the control in its log by three lines only: the L31
validity statement of `molecular_reaction_heat.f90`, an inventory count that
reads 30 instead of 28 equations because L36 registered two, and the wall
time. No data line differs. Build A differs from the measured build in its
log by the wall time alone.

`oxygen_chemistry` sets neither `He_diffusion` nor
`Molecular carrier transport` (READ, its `input.inp`), so it enters neither
the element diffusion operator nor the carrier transport operator, which is
why build A can be bitwise the measured build. Of the five files build C
reverts, only `species_face_flux.f90`, `RK_rhs.f90` and `Reconstruction.f90`
are reachable in this configuration: `steady_residual.f90` belongs to the
stationary solve, which never fires here (its hand-off is at `du < 1.0e-2`
against `du = 4.51` at the stop), and `hydrodynamic_rows_body.inc` is the
kind-generic instantiation. The three reachable files carry one and the same
respelling, so the seed is named by the expression and not by the file.

## `mol_sec_ion` is the other class

Measured the same way, on scratch copies, single-threaded, at a cap of 1000
steps and again at the case's own `maxsteps` of 12000. MEASURED.

| step | `du` control | `du` measured | relative |
|---|---|---|---|
| 1 | 2.508183462993202 | 2.508183462993202 | 0 (bitwise) |
| 2 | 2.509994386123802 | 2.509994386100904 | 9.122e-12 |
| 5 | 2.515429787393657 | 2.515429787368377 | 1.005e-11 |
| 10 | 2.524460399650085 | 2.524460399616931 | 1.313e-11 |
| 100 | 2.688664944875988 | 2.688664944762838 | 4.208e-11 |
| 300 | 3.092128710864969 | 3.092128710706079 | 5.139e-11 |
| 1000 | 5.389588808532107 | 5.389588932122480 | 2.293e-08 |

The first step is bitwise identical, the second is apart by 9e-12, and the
separation then rises smoothly by about three decades over 1000 steps. That
is a rounding difference being carried by a trajectory, and it is what a
seed of 5.8e-13 is expected to do. The snapshot at 1000 steps moves by
6.81e-05 at most (`Hydro_ioniz.txt`, cell 59, `v`) and by 8.44e-06 of the
worst column's own scale; `Ion_species.txt` by 4.55e-05 (`H3p`) and
6.30e-06 of the column scale.

Run to the case's own cap of 12000 steps as well, both builds stop at the
same step with the same printed marching line
(`final: count=12000 du= 7.4360E-01 dtu= 6.0495E-03`), and the matrix
movement is reproduced exactly (MEASURED):

| file | max relative | cell, column | control | measured | of the column's own scale |
|---|---|---|---|---|---|
| `Hydro_ioniz.txt` | 3.2102e-04 | 192, 3 `v` | 6.34034658e-01 | 6.34238260e-01 | 6.31e-07 (`heat`) |
| `Ion_species.txt` | 1.2680e-05 | 241, 35 `H2` | 7.195478e+06 | 7.195387e+06 | 3.93e-07 (`HI`) |

The 3.2e-4 sits at a cell whose velocity is 0.63 cm/s, so on the profile it
is 6.3e-07, three decades below the harness tolerance.

The decisive difference from `oxygen_chemistry` is that the root ladder does
not fire at all: over 1000 steps and over the full 12000 alike, both builds
report `no admissible root` at zero cells, no constrained-continuation
acceptance and no relaxation amnesty, and one stored-state rejection each
(MEASURED, counts of log lines, identical on the two builds). With no branch
to take, the case stays in the class the other fifteen are in, and its
movement is the carried rounding of the same volume respelling.

## Consequences, for the advisor to decide

1. The movement of `oxygen_chemistry` between the 2026-09-17 goldens and the
   current binary is **not** a defect of the L30 transport geometry and not
   a sensitivity of the oxygen carrier rows to the shell geometry. Case (ii)
   is refused by a bitwise experiment.
2. The case as configured is a **branch amplifier**, not a regression gate.
   Its `maxsteps` was chosen for wall time, not for a state (READ, its
   `README.md`), and the state it pins is 1000 steps into a cold isothermal
   transient. Whatever else is decided, its golden should carry a `NOTE.txt`
   that says so, so that the next refresh does not spend a session
   apportioning it again.
3. If a gate on the oxygen network is wanted rather than a reproducibility
   snapshot, it has to be a case whose chemistry solve has settled, or a
   row that reads the network directly (a rate, a branching ratio, a
   conservation identity) instead of a 1000-step profile.

## Method notes

- All runs single-threaded, on scratch copies, never through
  `run_check.sh`; nothing was written under `backup/regression/`.
- The tree's `EXHALE.x` was not touched: md5
  `267e0b2e85be71e7d66c221ecc58cf34` before and after (MEASURED).
- `build_L30c/` and `EXHALE_L30c.x` were deleted at the end.
