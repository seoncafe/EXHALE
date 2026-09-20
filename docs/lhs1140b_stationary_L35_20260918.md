# L35: the matched-domain flux and grid comparison, and the local base closure as a bounded study

Item L35 of `docs/PLAN_20260917.md` section 9, written 2026-09-18 (KST).
Two studies, neither of which changes a line of production source. Step 1
holds the physical face boundaries of the 500-cell catalog grid fixed and
solves `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7` with the HLLC and the Roe
flux at 500, 1000 and 2000 cells. Step 2 forms the local characteristic
closure of the base on the production face flux, on the entry point item L27
added, and compares the ghost it implies with the ghost the remote-window
closure of `base_boundary.f90` installs.

Every number below is MEASURED on the host `lart4` on 2026-09-18 unless it is
marked READ. Both recommendations are recommendations; no decision is taken
here and no source, golden or catalog directory was written to.

---

## 1. Verdict

**Recommendation 1, the Roe scope.** The seven-case Roe scope
stands, and the 1000-cell grid is NOT the catalog's answer, for a reason the
matched comparison makes visible and the L25 pair could not. MEASURED: on
matched physical face boundaries the unmodified HLLC flux, which has no
stationary root of this case at 500 cells (mass row 1.801e-01 and energy row
2.406e-01 at CELL 1, base face carrying 0.700 of the wind's own mass flux),
brings every hydrodynamic row inside its tolerance at the FIRST outer pass at
1000 and at 2000 cells and carries the wind's mass flux through every face to
1.3e-07. So refinement does relieve the stall. But each flux then converges
to its OWN base state and the two do not approach each other: between 1000
and 2000 cells each is converged to a few parts in 1e5 in the first cell and
in Mdot, while the difference BETWEEN them is the same at both resolutions to
three digits, 10.3 per cent in T(1), 11.1 per cent in rho(1), 1.0 per cent in
Mdot and 1.8 per cent in the He I 10830 equivalent width. And refining at
fixed domain costs the certification: all four refined runs are refused by a
hydrodynamic energy row of 2.7e-04 to 8.2e-04 against 1.0e-06 at r = 4.3 to
5.0 R_p, at both fluxes, which the 500-cell Roe state does not carry. The one
certified state of the six is Roe at 500 cells. What this study DOES add to
the R37 list is the Roe flux's own grid convergence, which was the missing
measurement: it is converged at the catalog's own resolution.

**Recommendation 2, the local base closure.** A source item is warranted, but
a narrow one, and not the one the exclusion of L26 section R3 would have had
to be reversed for. MEASURED: on both certified states, on a base whose
velocity is cut by a hundred, on a base whose velocity is reversed, on a
perturbation confined to the remote window, and on the 300-step hydrostatic
column, the local characteristic closure built on the production base face
flux and the retained remote-window closure install the SAME ghost cell,
bit for bit in all three primitive components. They part on exactly two kinds
of state, and both are states no converged run holds: a column at rest, where
the window carries no flux at all and the remote closure falls back on the
50/50 entropy blend while the local one reads the reservoir, a 3.79 per cent
difference in the face density of the fiducial and 3.07e-07 on the
hydrostatic column; and a mapped seed on an extended domain, where the
production base face flux is NEGATIVE (-5.90 of the window mean) and the two
closures select different branches, a 7.30 per cent difference in the face
density. So the local closure is not refuted by the states, which is the
finding: the L26 section R3 exclusion, taken on cell-centered readings and on
PLM evaluations of WENO3 states, does not survive the production face flux.
It is also not established by them, because on every state where the two
agree they agree only because the branch weight is saturated at both
discriminants, and the two states where they differ are the two where the
face flux is not the wind's. The bounded source item this supports is the
zero-window limit, not the closure family: at a window of exactly zero flux
the production boundary returns the average of two isentropes that differ by
7.3 per cent on this planet, and a face flux reading would return the
reservoir. See section 6.

---

## 2. The binary, and where it came from

No production source was changed by this item. The tree carries the work of
other items in progress, so the study was made on a SNAPSHOT of the working
tree, taken at `2026-09-18T01:39:21+0900` and built by itself, and not on
the tree binary `EXHALE.x` (md5 `59bfdb3fc4d0104fc2e9c3734596d2f6`), which
is another item's record binary and was not touched.

| | |
|---|---|
| snapshot | `<scratchpad>/L35/tree/` (`Makefile` and `src/`) |
| taken at | 2026-09-18T01:39:21+0900 |
| repository HEAD at the snapshot | `3c73905ca8a7fe2af92a5c2b014c225b2feedff3` |
| working tree at the snapshot | dirty: `32 files changed, 2341 insertions(+), 496 deletions(-)`, plus 13 untracked paths |
| compiler | conda-forge gfortran 16.2.0, `make` with no override, OpenBLAS of the same prefix |
| binary | `<scratchpad>/L35/tree/EXHALE.x`, md5 `16f5345d07ce70cdac09227839d5188b` |

The snapshot compiled at the first attempt, so no re-snapshot was needed. The
files changed at snapshot time, READ from `git diff --stat`, are
`LHS1140b/MODELS.md`, `Makefile`, `docs/input_schema.md`, three memos,
`src/EXHALE_main.f90`, `src/modules/flux/species_face_flux.f90`,
`src/modules/functions/binary_element_diffusion.f90`,
`src/modules/init/define_grid.f90`,
`src/modules/lower_atmosphere/{diffusive_photochemistry,h2_vibrational_relaxation,mol_rates,molecular_reaction_heat}.f90`,
`src/modules/states/Reconstruction.f90`,
`src/modules/time_step/{RK_rhs,certification,hydrodynamic_rows_body.inc,steady_newton,steady_residual}`,
`src/utils/map_state_to_grid.py` and eight test drivers; the untracked
additions include `src/modules/states/stationary_operator.f90`, which is item
L27's entry point and is in `SRC`. **The numbers of this memo are therefore
not the catalog binary's**: they are a build of the tree as items L27, L28,
L30 and L31 had left it at 01:39 KST on 2026-09-18, and the step 1 numbers
are not to be compared digit for digit with the L25 tables, which were made
on `c2e9c9990b9f14f1be8cd77abca68945`.

The probe of step 2 is `src/tests/grid_and_gates/base_closure_candidates_probe.f90`,
a new file. It is NOT registered in `src/tests/grid_and_gates/run.sh`, which
another item owns at the moment; it was compiled by hand with the command the
suite uses for its own probes,

```
gfortran -O0 -g -fbacktrace -fopenmp -J<out> -I<objdir> \
   -o <out>/base_closure_candidates_probe.x \
   src/tests/grid_and_gates/base_closure_candidates_probe.f90 \
   $(ls <objdir>/*.o | grep -vE '(EXHALE_main|_tests|_probe)\.o$') \
   -L/opt/miniconda3/lib -lopenblas -Wl,-rpath,/opt/miniconda3/lib -ldl
```

with `<objdir>` the snapshot's `build/`, and run at `OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1`. Registering it would be one `want` block following
the one `base_continuity` already has, gated on an environment variable
naming the state directory, so that it never runs by default.

---

## 3. Step 1: holding the physical face boundaries fixed

### 3.1 What the grid generator lets a matched refinement be

READ from `src/modules/init/define_grid.f90`. The `Mixed` grid stacks
`N_low` uniform cells of width `dr_base` on `r = 1`, fills the rest of the
domain with `N - N_low` geometrically stretched cells, smooths the widths
with one 1-2-1 pass down to the base, rebuilds the centers from the smoothed
widths, and then rescales the whole column so that the OUTERMOST GHOST center
sits at `r_max`. Two consequences decide what a matched-domain refinement can
be on this grid.

- The second ghost center is set to `r(2-Ng) = 1` before the rescale and the
  rescale fixes `r = 1`, so `r(0) = 1` exactly at every resolution and the
  inner physical face is `r_edg(0) = (1 + r(1))/2`. **Holding `r_edg(0)`
  fixed therefore holds the width of the FIRST CELL fixed**: the inner
  boundary radius and the first cell width are one parameter on this grid,
  and a matched-domain refinement cannot refine the cell the base boundary
  condition speaks to. That is a property of the generator, not a choice
  made here, and it is the first thing this step measures.
- The last physical face `r_edg(N)` is not `r_max`: `r_max` is where the
  outermost GHOST center is put, and `r_edg(N)` follows from the smoothed
  widths, so it moves with `N` at fixed `r_max`. That is the 1.9 per cent
  motion of the outer end the L25 grid study reported (READ, L25 memo
  section 3.5).

Two conditions, two keys. `Base grid [dr,cells]: <dr_base> 50` sets the first
cell width and `Outer radius [R_p]: <r_max>` sets the rescale, so the pair
(`dr_base`, `r_max`) is solved for the pair (`r_edg(0)`, `r_edg(N)`) of the
500-cell catalog grid. The solve was made on an exact replica of
`define_grid`'s `Mixed` branch in Python and CHECKED against the binary: the
replica reproduces all 504 rows of the catalog state's radius column with a
maximum relative difference of **0.0** in double precision.

### 3.2 The keys, and the faces they produce

MEASURED by running the binary once per grid (one step, `Load IC? False`) and
reading `r` from `output/Hydro_ioniz.txt`, whose radius column is written
list-directed and therefore carries the full double-precision value.

| `Grid cells` | `Base grid [dr,cells]` | `Outer radius [R_p]` | `r_edg(0)` | `r_edg(N)` |
|---|---|---|---|---|
| 500 | `0.00019999999494757503 50` | `30.0` | 1.0000966593889191 | 29.2734180458963 |
| 1000 | `0.00019615093524286384 50` | `29.582942490104671` | 1.0000966593889191 | 29.2734180458963 |
| 2000 | `0.00019453512202014908 50` | `29.406628764578247` | 1.0000966593889191 | 29.2734180458963 |

**The first and last physical face radii agree to 0.0, not merely to 1e-10**:
the three grids carry the same bits in both. The 500-cell row is the catalog
grid itself; stating `dr_base` explicitly at its 17 digits reproduces the
default bit for bit (section 8 says why the 17 digits are needed).

What the refinement then is, MEASURED as the cell width at five radii:

| cell width `dr_j` at | 500 | 1000 | 2000 | 500/1000 | 1000/2000 |
|---|---|---|---|---|---|
| `r = 1.00019` (cell 1) | 1.933188e-04 | 1.933188e-04 | 1.933188e-04 | 1.000 | 1.000 |
| `r = 1.01` | 1.984625e-04 | 1.954652e-04 | 1.942337e-04 | 1.015 | 1.006 |
| `r = 1.05` | 8.885798e-04 | 4.897267e-04 | 3.202620e-04 | 1.814 | 1.529 |
| `r = 1.20` | 3.521475e-03 | 1.586813e-03 | 7.923003e-04 | 2.219 | 2.003 |
| `r = 2.0` | 1.750501e-02 | 7.478980e-03 | 3.303821e-03 | 2.341 | 2.264 |
| `r = 10.0` | 1.573908e-01 | 6.631190e-02 | 2.848854e-02 | 2.373 | 2.328 |

So this is a clean factor-2 refinement of the wind above about 1.05 R_p, a
partial one between 1.01 and 1.05, and NO refinement at all of the 50-cell
uniform block that carries the base, whose first cell is identical in all
three by construction. Any statement this step makes about the base is
therefore a statement about the same base cells under three different outer
resolutions, and not about the convergence of the base itself. A study that
refines the base has to move `r_edg(0)`, which is a different comparison.

### 3.3 The recipe

The case is `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7`, the L25 headline case
and one of the seven low-XUV cases the catalog carries `Numerical flux: ROE`
for. Its `input.inp` is the L25 step 3 file `.L25/roe_s010_HeH9.7/input.inp`
(the catalog file with the spectrum path made absolute), with three lines
changed or added per run: `Numerical flux:`, `Grid cells:`, `Base grid
[dr,cells]:` and `Outer radius [R_p]:`. Nothing else differs between the six
runs, so the sources, the boundary model, `Resid tol`, the composition and
`He_Kzz` are the catalog's in all six.

The seed of every run is the CERTIFIED 500-cell Roe state
`.L25/roe_s010_HeH9.7/output` (the catalog case directory carries no state
and no `seed_from_ladder.txt`), mapped onto each grid by

```
python3 src/utils/map_state_to_grid.py .L25/roe_s010_HeH9.7/output \
    .L35/grid/g<N>/output/Hydro_ioniz.txt .L35/g<N>_<flux>/output --ic
```

with the target grid file the one-step probe of section 3.2 wrote. No
extrapolation option was needed: all three matched grids lie inside the
catalog grid's span. The mapper's output states `certified=F
cert_reason=mapped_seed`, which is what it states for every mapped seed.

The solve is the partitioned stationary route with `Restart intent:
stationary` and `EXHALE_PTC_DTAU0=1.0`. The L25 step 3 runs of this case were
taken at `EXHALE_PTC_DTAU0=1.0e8` directly (READ, L25 memo section 3.2); the
plan asks for 1.0 here, and the first attempt went through
`LHS1140b/models/run_case.sh` (`EXHALE_OUTER_PASSES=40`, `NOSEED=1`,
`FORCE=1`, `SEED_ATTEMPTS=1`), whose own continuation to `dtau0 = 1e8` covers
the 1e8 path when the first solve does not return `info = 0`. Section 3.4
says what that attempt measured, why it was stopped and what replaced it.

Hosts and threads. All six ran first on `lart4` (72 cores, load average 3
before the launch), at 8, 12 and 16 threads by cell count. The re-runs put
the 500-cell and the two 1000-cell cases on `lart4` at 8 and 16 threads and
the two 2000-cell cases on `lart3`, which was idle at load 0 throughout
(72 cores), at 24 threads each. **That is a departure from the brief, which
names `lart4`**: the session itself runs on `lart4`, so the six runs and the
two other items' binaries had put it at load 105 of 72 cores while `lart3`
carried nothing, and the 2000-cell pair was moved rather than left to
contend. `OMP_NUM_THREADS` changes the wall clock and not the answer (READ,
`REPRODUCE.md` of every catalog case); the binary was copied to
`.L35/EXHALE_L35.x` so both hosts read one file, md5
`16f5345d07ce70cdac09227839d5188b`, the snapshot build of section 2.

### 3.4 What came out

Five of the six were solved with a BOUNDED outer-pass budget, and the budget
and the reason are part of the measurement.

The first attempt gave all six the catalog's 40 passes. Roe at 500 cells was
ACCEPTED at its first outer pass, certified, post-processed and finished in
31 s; it was never re-run and its record is the runner's own `REPRODUCE.md`.
The other five showed the following within three hours. HLLC at 500 does not
leave its THIRD outer pass: its first two passes take 166 s each, and it then
spends 327 JFNK iterations at `||R||` 2.4e-01 without progress, reporting 500
of 500 mass cells outside their tolerance with the worst at cell 1; it was
still there at 55 minutes, which is the same behaviour L25 measured on this
configuration (252 minutes, 5 passes, stopped; READ, L25 memo section 3.9).
At 1000 and 2000 cells BOTH fluxes drive the hydrodynamic mass and momentum
rows inside their tolerances at the first pass, but leave a hydrodynamic
ENERGY row that rises from 1e-07 at the first passes to a plateau of 2.7e-04
to 8.5e-04 at r ~ 5 R_p by the seventh and then stops falling (the table at
the end of section 3.5). No refined run can therefore certify, and the
remaining thirty passes could only have relaxed the elemental partition,
which halves every pass and refuses nothing the energy row does not already
refuse.

The five were stopped by PID, with `/proc/<pid>/cwd` checked against this
item's own directories, and their logs kept as
`<case>/run_passes40_stopped.log`. Each was then re-solved from the same
mapped seed with 2 outer passes (HLLC 500, whose third pass is the stall) and
8 (the four refined runs, the plateau being reached by the seventh), through
`.L35/solve_and_post.sh`, which is the catalog recipe of `run_case.sh` with
its refusal to make products for a non-certifying solve removed, so that each
writes a state and its advection-corrected profiles and transit spectrum are
made the way the catalog's are. **Only the Roe 500-cell run is a certified
state; the other five are relaxation states and their rows say so.**

#### outcome

| run | outcome | outer passes | `||R||` | rows that refuse the written state | wall clock, threads, host |
|---|---|---:|---:|---|---|
| HLLC 500 | **NOT CERTIFIED** | 2 of 2 | 2.406e-01 | mass 1.801e-01 at cell 1 (tol 1.2e-07); momentum 7.988e-08 at cell 500 (1.0e-08); energy 2.406e-01 at cell 1 (1.0e-06) | 7m, 8, lart4 |
| HLLC 1000 | NOT CERTIFIED | 8 of 8 | 8.215e-04 | energy 8.215e-04 at cell 734, r = 4.98 (1.0e-06); elemental He/H partition 3.075e-05 at cell 635, r = 2.91 (1.0e-05) | 70m, 16, lart4 |
| HLLC 2000 | NOT CERTIFIED | 8 of 8 | 7.767e-04 | energy 7.767e-04 at cell 1375, r = 4.91 (1.0e-06); elemental He/H partition 4.388e-05 at cell 1152, r = 2.91 (1.0e-05) | 77m, 24, lart3 |
| ROE 500 | **CERTIFIED** | 1 of 40 | 4.233e-07 | none | 0m31s, 8, lart4 |
| ROE 1000 | NOT CERTIFIED | 8 of 8 | 3.831e-04 | energy 3.831e-04 at cell 707, r = 4.26 (1.0e-06); elemental He/H partition 3.092e-05 at cell 635 (1.0e-05) | 99m, 16, lart4 |
| ROE 2000 | NOT CERTIFIED | 8 of 8 | 2.689e-04 | energy 2.689e-04 at cell 1345, r = 4.55 (1.0e-06); elemental He/H partition 4.386e-05 at cell 1152 (1.0e-05) | 73m, 24, lart3 |

The wall clocks are not a cost model: the first three shared `lart4` with each
other and with two other items' binaries, and the 2000-cell pair had `lart3`
to itself.

**Why HLLC 500 does not solve, stated by the rows rather than by the norm.**
Its mass and energy rows are at CELL 1, at 1.8e-01 and 2.4e-01 of their own
scales, while its momentum row is four decades smaller and at the outer
boundary. The JFNK reports 500 of 500 mass cells outside their tolerance with
the worst at cell 1. So the stall is the base, and the next table says the
same thing in the conserved flux.

#### the conserved face mass flux budget of the written state

Printed by the certification report of item L27, in units of the wind-window
mean of `rho v r^2` on the same state and under the same operator (WENO3, the
case's own numerical flux, `Well balanced: True`).

| run | window mean [code] | base face | minimum over the faces | maximum over the faces | (max - min)/mean |
|---|---:|---:|---:|---:|---:|
| HLLC 500 | 5.49215e-08 | **0.700465914** | 0.700465914 at face 0 | 1.064386416 at face 4 | **3.639e-01** |
| HLLC 1000 | 5.54293e-08 | 0.999999897 | 0.999999855 at face 5 | 0.999999897 at face 0 | 4.162e-08 |
| HLLC 2000 | 5.54080e-08 | 1.000000359 | 1.000000226 at face 11 | 1.000000359 at face 0 | 1.328e-07 |
| ROE 500 | 5.49209e-08 | 0.999974000 | 0.999973949 at face 1 | 0.999974028 at face 7 | 7.816e-08 |
| ROE 1000 | 5.48775e-08 | 0.999999979 | 0.999999969 at face 27 | 0.999999979 at face 0 | 9.455e-09 |
| ROE 2000 | 5.48576e-08 | 1.000000250 | 1.000000247 at face 95 | 1.000000251 at face 4 | 4.719e-09 |

**Five of the six states carry one mass flux through every face of the
column, to between 5e-09 and 1.3e-07 of the wind's own.** The sixth is HLLC
at 500 cells, whose base face carries 0.700 of the wind's flux and whose
faces spread by 36 per cent of it over the first few cells. That is the stall
seen in the conserved quantity, and it is a base-face statement, not a column
one. It is also the only state of this study on which the withdrawn L26
section R3 reading would have had anything to point at.

#### the integrated column energy balance over the wind, 1.2 <= r <= 25 R_p

The balance is NOT taken over the whole column: the cell-centered product
`rho v r^2` is not the mass flux in the base layer (MEASURED -226 times the
wind's at cell 1 of the ROE 500 state), so a boundary term read there is not
a flux. 1.2 R_p is the radius the certification gates the wind at and 25 R_p
lies inside every one of the three domains, so the same two physical surfaces
are used at all three resolutions. The heating and the cooling are the
`Heating_breakdown.txt` and `Cooling_breakdown.txt` totals integrated over
the spherical shells between them; the three flux terms are
`4 pi r^2 rho v X` differenced between the two surfaces, with `X` the
enthalpy `5p/2rho`, the kinetic energy `v^2/2` and the planetary potential
`-GM/r`. All in erg s^-1.

| run | heating | cooling | advected enthalpy | kinetic | gravitational work | closure residual / heating |
|---|---:|---:|---:|---:|---:|---:|
| HLLC 500 | 9.55954e+18 | 1.28697e+18 | -1.80849e+18 | +8.64113e+14 | +1.00956e+19 | +1.610e-03 |
| HLLC 1000 | 9.64709e+18 | 1.27192e+18 | -1.79534e+18 | +8.75586e+14 | +1.01858e+19 | +1.679e-03 |
| HLLC 2000 | 9.65662e+18 | 1.27107e+18 | -1.79384e+18 | +8.74150e+14 | +1.01811e+19 | +2.658e-04 |
| ROE 500 | 9.55934e+18 | 1.28622e+18 | -1.80836e+18 | +8.64106e+14 | +1.00950e+19 | +1.509e-03 |
| ROE 1000 | 9.54244e+18 | 1.28001e+18 | -1.80563e+18 | +8.64034e+14 | +1.00844e+19 | +1.800e-03 |
| ROE 2000 | 9.55211e+18 | 1.27930e+18 | -1.80416e+18 | +8.62525e+14 | +1.00799e+19 | +3.994e-04 |

The column balance closes to between 2.7e-04 and 1.8e-03 of the heating in
all six, which is the size of the quadrature and of the interpolation onto
the two fixed surfaces and not a physical imbalance. Nothing in this table
separates the two fluxes: heating differs between them by 1.1e-02 at 1000 and
2000 cells and by 2.1e-05 at 500, cooling by 6.4e-03 and 5.9e-04. This is the
wind, and the wind is the same object in all six.

#### the base, the mass-loss rate and the line

`Mdot` at `r_out` is `4 pi rho v r^2` of the outermost physical cell; `Mdot`
at r = 25 is the same product at the cell nearest 25 R_p, a fixed physical
surface, since the outermost cell CENTER is not the same radius at the three
resolutions. The equivalent width is the red-pair EW over the measurement's
vacuum window 10832.60 to 10834.20 A, read the way
`LHS1140b/make_memo_figures.py` reads it. The `log10 Mdot` column is the
binary's own line from `pp.log`, printed to two decimals.

| run | T(1) [K] | rho(1) [mH/cm3] | v(1) [cm/s] | Mdot at r_out [g/s] | Mdot at r = 25 [g/s] | log10 Mdot (binary) | He I 10830 EW [%A] |
|---|---:|---:|---:|---:|---:|---:|---:|
| HLLC 500 | 221.787 | 1.11762e+14 | -4.8548e-01 | 6.41368e+06 | 6.42309e+06 | 6.81 | 0.2355 |
| HLLC 1000 | **247.299** | **1.00622e+14** | -1.0835e+00 | 6.47885e+06 | 6.48327e+06 | 6.81 | 0.2351 |
| HLLC 2000 | **247.308** | **1.00618e+14** | -1.0837e+00 | 6.47893e+06 | 6.48091e+06 | 6.81 | 0.2331 |
| ROE 500 | 221.776 | 1.11767e+14 | -4.8524e-01 | 6.41367e+06 | 6.42308e+06 | 6.81 | 0.2353 |
| ROE 1000 | **221.777** | **1.11767e+14** | -4.8526e-01 | 6.41434e+06 | 6.41872e+06 | 6.81 | 0.2309 |
| ROE 2000 | **221.777** | **1.11767e+14** | -4.8526e-01 | 6.41456e+06 | 6.41653e+06 | 6.81 | 0.2290 |

**HLLC 500's base is not its own answer**: that run never left its seed, so
its first cell is the mapped ROE 500 state's to five digits, and it is in
this table as the state the stall wrote and not as an HLLC solution of the
base.

### 3.5 The resolution trend and the flux difference

**Each flux converges, and they converge to different base states.** The two
tables below are the same numbers read the two ways the decision needs.

#### the resolution trend of each flux

| quantity | HLLC 500 | HLLC 1000 | HLLC 2000 | (1000-500)/500 | (2000-1000)/1000 |
|---|---:|---:|---:|---:|---:|
| T(1) [K] | 221.787 | 247.299 | 247.308 | +1.15e-01 | **+3.9e-05** |
| rho(1) [mH/cm3] | 1.11762e+14 | 1.00622e+14 | 1.00618e+14 | -9.97e-02 | **-3.8e-05** |
| Mdot at r_out [g/s] | 6.41368e+06 | 6.47885e+06 | 6.47893e+06 | +1.02e-02 | **+1.3e-05** |
| Mdot at r = 25 [g/s] | 6.42309e+06 | 6.48327e+06 | 6.48091e+06 | +9.37e-03 | -3.6e-04 |
| He I 10830 EW [%A] | 0.2355 | 0.2351 | 0.2331 | -1.90e-03 | -8.3e-03 |
| heating [erg/s] | 9.55954e+18 | 9.64709e+18 | 9.65662e+18 | +9.16e-03 | +9.9e-04 |
| cooling [erg/s] | 1.28697e+18 | 1.27192e+18 | 1.27107e+18 | -1.17e-02 | -6.7e-04 |

| quantity | ROE 500 | ROE 1000 | ROE 2000 | (1000-500)/500 | (2000-1000)/1000 |
|---|---:|---:|---:|---:|---:|
| T(1) [K] | 221.776 | 221.777 | 221.777 | **+2.5e-06** | **+1.2e-06** |
| rho(1) [mH/cm3] | 1.11767e+14 | 1.11767e+14 | 1.11767e+14 | **-2.4e-06** | **-1.2e-06** |
| Mdot at r_out [g/s] | 6.41367e+06 | 6.41434e+06 | 6.41456e+06 | +1.05e-04 | +3.5e-05 |
| Mdot at r = 25 [g/s] | 6.42308e+06 | 6.41872e+06 | 6.41653e+06 | -6.78e-04 | -3.4e-04 |
| He I 10830 EW [%A] | 0.2353 | 0.2309 | 0.2290 | -1.89e-02 | -8.2e-03 |
| heating [erg/s] | 9.55934e+18 | 9.54244e+18 | 9.55211e+18 | -1.77e-03 | +1.0e-03 |
| cooling [erg/s] | 1.28622e+18 | 1.28001e+18 | 1.27930e+18 | -4.83e-03 | -5.5e-04 |

The Roe column is converged at 500 cells already: doubling and quadrupling
the resolution moves its first cell by 2.5e-06 and 1.2e-06 and its mass-loss
rate by 1e-04 and 3e-05. The HLLC column moves by 11.5 per cent in T(1) and
1.0 per cent in Mdot between 500 and 1000, and then by 3.9e-05 and 1.3e-05
between 1000 and 2000. **But the HLLC 500 row is not an HLLC answer**: that
solve stalled on its seed, which was the Roe 500 state, so the large
"500 to 1000" step of the HLLC table is the distance from the Roe base to the
HLLC base and not a truncation error. Read as a convergence statement, what
the table establishes is the SECOND column of each: between 1000 and 2000
cells both fluxes are converged to a few parts in 1e5 in the base and in the
mass-loss rate.

The equivalent width is the one quantity still moving at 2000 cells, by
-8.3e-03 and -8.2e-03 between 1000 and 2000 at both fluxes, i.e. it has not
converged and is falling at about the same rate in both. Its 1000-to-2000
change is the same size as its flux-family difference, so no statement about
the line is supported at better than about 2 per cent by this study.

#### the flux difference at each resolution, (ROE - HLLC)/HLLC

| quantity | 500 | 1000 | 2000 |
|---|---:|---:|---:|
| T(1) [K] | -4.98e-05 | **-1.032e-01** | **-1.032e-01** |
| rho(1) [mH/cm3] | +4.80e-05 | **+1.108e-01** | **+1.108e-01** |
| Mdot at r_out [g/s] | -2.43e-06 | **-9.96e-03** | **-9.94e-03** |
| Mdot at r = 25 [g/s] | -2.43e-06 | -9.96e-03 | -9.94e-03 |
| He I 10830 EW [%A] | -8.62e-04 | -1.79e-02 | -1.78e-02 |
| heating [erg/s] | -2.08e-05 | -1.09e-02 | -1.08e-02 |
| cooling [erg/s] | -5.87e-04 | +6.36e-03 | +6.47e-03 |

**The difference between the two fluxes does not shrink with resolution; it
is the same number at 1000 and at 2000 cells to three digits.** The 500-cell
column is near zero only because the HLLC solve there never left the Roe
seed. So the Roe base and the HLLC base are two converged discrete answers to
the same continuous problem, differing by 10.3 per cent in the first cell's
temperature, 11.1 per cent in its density and 1.0 per cent in the mass-loss
rate, and refining the grid at fixed physical boundaries does not bring them
together.

Profile against profile at matched radii, (ROE - HLLC)/HLLC, MEASURED at
1000 and 2000 cells (the two columns agree to the digits printed, which is
the same statement again):

| r [R_p] | T | rho | v |
|---|---:|---:|---:|
| 1.0002 | -1.07e-01 | +1.14e-01 | +5.4e-01 |
| 1.0005 | -2.12e-01 | +2.36e-01 | -1.7e-01 |
| 1.0010 | -1.81e-01 | +1.53e-01 | +9.0e-01 |
| 1.0020 | -7.89e-02 | -6.0e-03 | -1.2e-02 |
| 1.0050 | -1.08e-02 | -9.00e-02 | +8.8e-02 |
| 1.0100 | +9.8e-03 | -1.08e-01 | +1.10e-01 |
| 1.0500 | +1.59e-02 | -8.56e-02 | +8.3e-02 |
| 1.2000 | +1.14e-02 | -3.11e-02 | +2.2e-02 |
| 2.0000 | -1.1e-03 | -1.14e-02 | +1.5e-03 |
| 5.0000 | -5.2e-03 | -8.4e-03 | -1.5e-03 |
| 10.000 | -4.6e-03 | -8.1e-03 | -1.9e-03 |
| 25.000 | -4.3e-03 | -8.3e-03 | -1.6e-03 |

The difference is 10 to 24 per cent below 1.01 R_p, still 3 per cent in the
density at 1.2 R_p, and about 0.5 to 0.8 per cent out to the top of the
domain. It is not confined to two or three cells.

#### what neither flux does at 1000 or 2000 cells

Neither certifies. The refusing entry common to all four refined runs is the
hydrodynamic ENERGY row, at 2.7e-04 to 8.2e-04 against 1.0e-06, at a cell at
r = 4.3 to 5.0 R_p in every one of them. Read from the pass records of the
40-pass attempt (`<case>/run_passes40_stopped.log`), that row is at 3e-08 to
4e-07 for the first two or three outer passes and rises to its plateau as the
elemental partition relaxes, then stops moving:

| outer pass | HLLC 1000 | ROE 1000 | HLLC 2000 | ROE 2000 |
|---|---:|---:|---:|---:|
| 1 | 1.08e-07 | 4.69e-08 | 3.86e-07 | 6.32e-08 |
| 2 | 1.29e-07 | 1.98e-04 | 3.49e-07 | 3.25e-08 |
| 3 | 1.60e-07 | 2.97e-04 | 2.98e-07 | 1.40e-04 |
| 4 | 4.10e-04 | 3.44e-04 | 3.86e-04 | 2.09e-04 |
| 5 | 6.25e-04 | 3.66e-04 | 5.89e-04 | 2.42e-04 |
| 6 | 7.35e-04 | 3.76e-04 | 6.94e-04 | 2.58e-04 |
| 7 | 7.92e-04 | 3.81e-04 | - | 2.65e-04 |
| 8 | 8.21e-04 | 3.83e-04 | - | - |
| 9, 10, 11 | 8.37, 8.45, 8.49e-04 | - | - | - |

(the 40-pass attempt; it was stopped at different passes in the four runs.
The bounded 8-pass re-runs end at 8.215e-04, 3.831e-04, 7.767e-04 and
2.689e-04, which is the same plateau.)

The 500-cell Roe run does not do this: it holds 4.2e-07. So at these
resolutions the coupled system has a wind-cell energy imbalance that grows
with the composition relaxation and saturates above the tolerance, at both
fluxes, and it is the reason no refined run of this study is a certified
state. Whether that is a property of the refined grids, of the coupling, or
of the binary this study was built on (section 2: the tree carries item L30's
work in progress on the residual assembly) is NOT isolated here and is the
first thing a follow-up would separate.

### 3.6 The figure

`docs/figures/lhs1140b_L35_flux_grid.png`: T, rho and v of the six runs, the
whole column in log r on the left and the base 1.00 to 1.05 R_p on the right.

Three columns and not two. The two the plan names are the first two: the
whole column in log r, and the base 1.00 to 1.05 R_p. A third was added
because the two fluxes' base difference is 10 per cent of a density that
falls by three decades across that panel, so on the logarithmic axis the six
curves lie on top of one another and the panel shows nothing; the third
column plots each run as a RATIO to the ROE 2000 run over the same interval,
against `r - 1` on a logarithmic axis, which is where the difference is
legible. `LHS1140b/models/.L35/make_figure.py` draws it.

What the figure shows, and what to look for. In the left column the six runs
are one curve: the wind, its temperature maximum of 7.55e+03 K at 1.393 R_p
and its velocity rising to 1.666e+04 cm/s at the top, is the same object at
both fluxes and all three resolutions. In the right-hand ratio column there are
exactly two curves, not six: HLLC 1000 and HLLC 2000 lie on each other, ROE 500,
1000 and 2000 lie on each other, and HLLC 500 lies on the ROE curves because
its solve never left the Roe seed. MEASURED as HLLC 2000 over ROE 2000, the
density ratio runs from 0.799 at r - 1 = 5.8e-04 to 1.122 at 1.3e-02 and the
temperature ratio from 1.290 at 5.8e-04 to 0.983 at 2.3e-02. The velocity ratio panel is noisy below r - 1 = 2e-03 because the
base velocity passes through zero there and the ratio of two small numbers of
either sign is not a useful measure; the velocity itself is in the middle
column, where the base oscillation of the first two or three cells is
visible at both fluxes and does not fall with resolution.

---

## 4. Step 2: the local characteristic closure on the production face flux

### 4.1 The characteristic analysis, and what is actually in question

The radial coordinate increases outward, so at the lower face `r_edg(0)` the
eigenvalues are `v - c`, `v` and `v + c` and a wave ENTERS the domain when its
eigenvalue is positive; gas entering at the base has `v > 0`. On the subsonic
inflow branch `0 < v < c` the two entering waves carry the reservoir pressure
and the reservoir specific entropy and the single outgoing wave `v - c < 0`
carries the linearized left-running acoustic compatibility relation

```
    p_b - rho_i c_i v_b  =  p_i - rho_i c_i v_i .                      (C-)
```

That is the LODI form of Thompson (1987, J. Comput. Phys. 68, 1) and Poinsot
and Lele (1992, J. Comput. Phys. 101, 104), in the statement of Carlson
(2011, NASA/TM-2011-217181, section 2). Those are the sources
`src/modules/states/base_boundary.f90` cites at its head (READ, lines 23 to
50), and the probe re-derives none of them.

**The closure is not in question; the discriminant is.** The branch is told
apart by the SIGN of the face velocity, and at a base face Mach number of
about 3e-7 the cell-centered state does not resolve it. The production
closure therefore reads a mass flux where a cell-centered product IS a flux,
over the remote wind window `r >= r_flux`, and maps it onto the face through
the interior density,

```
    M_wind = F_wind / (rho_i r_b^2) / c_i ,
```

with `F_wind` the window mean of `rho v r^2`
(`base_boundary.f90:characteristic_base_face_state`, READ). The local
alternative is the same relation read at the face itself, from the mass flux
the PRODUCTION Riemann solve puts through that face under the stationary
operator,

```
    M_face = [ r_edg(0)^2 (rho v)_0 ] / (rho_i r_b^2) / c_i ,
```

which is `budget%flux_base` of
`stationary_operator.f90:stationary_face_mass_flux`, item L27's entry point.

The probe calls ONE production routine twice on one installed state:
`characteristic_base_face_state` takes the flux, its relative spread and an
availability flag as arguments, so passing `(F_wind, d_window, have_F)` gives
the row `base_boundary_states` installs, and passing `(budget%flux_base, 0,
.true.)` gives the local row. The spread is zero because one face carries one
flux. Everything else is the same production code in both rows: the reservoir
isentrope, (C-), the isentropic density at `p_b`, the supersonic-outflow
blend, the face Mach cap and the eight-point Gauss-Legendre ghost quadrature.
The difference printed is the difference between the two discriminants and
nothing else. A CHECK in every run confirms it: the ghost `Apply_BC` installs
on the state is the remote row, bit for bit.

### 4.2 The states

| name | what it is |
|---|---|
| `fid_resolve` | the certified fiducial `atomic_scalar_gj1132_kzz1e9/HeH2.13` re-solved under the repaired boundary, `.L26/fid_resolve/output/*_IC.txt` |
| `x010_HeH2.13` | the certified `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13` state, the catalog `output/*_IC.txt` pair |
| `hydrostatic_column` | `backup/regression/hydrostatic_column`, the mechanical column with the radiation off, as its 300-step record leaves it |
| `weak` | `fid_resolve` with the velocity of physical cells 1 to 6 scaled by 1e-2 at fixed mass density and fixed thermal energy |
| `reversed` | `fid_resolve` with the velocity of physical cells 1 to 6 multiplied by -1, same invariants |
| `window` | `fid_resolve` with the mass density of every cell at `r >= r_flux` scaled by 1 + 1e-3, the momentum and the total energy carried with it, so that only the remote window moves |
| `rest` | every cell's momentum set to zero at fixed mass density and fixed thermal energy, on `fid_resolve` and on `hydrostatic_column` |
| `rmax45` | `fid_resolve` mapped onto a 1.5x extended domain, `Outer radius [R_p]: 45.0`, 500 cells, `--extrapolate-beyond 28.0` |

Two of these are not what their names might promise, and the memo says so
rather than letting the table imply it. `hydrostatic_column` as its record
leaves it is a 300-step relaxation snapshot and NOT a column at rest: its
first cell carries `M_i = +0.2257`, an inflow of 8.09e+04 cm/s. That is why
the `rest` perturbation exists, and it is applied to both wind states. And
`rmax45` is a mapped seed on the extended domain, not a solution of it, so it
is a state the operator does not vanish on; that is exactly what makes it the
discriminating case below.

### 4.3 The two closures, side by side

MEASURED. `F_face` is the production base face mass flux under the stationary
operator (WENO3, the case's own numerical flux, well balanced as the case
states it); the density columns are `mH cm^-3`, the velocity `cm s^-1`.

| state | `F_face / F_window` | `M_wind` (remote) | `M_face` (local) | `w_rev` remote | `w_rev` local | face rho remote | face rho local | (local - remote)/remote |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `fid_resolve` | +9.9998437e-01 | +2.991741e-07 | +2.991694e-07 | 0 | 0 | 9.505840969e+13 | 9.505840969e+13 | **0.0** |
| `x010_HeH2.13` | +9.9997483e-01 | +2.677849e-08 | +2.677781e-08 | 0 | 0 | 9.505840969e+13 | 9.505840969e+13 | **0.0** |
| `hydrostatic_column` | +3.5648933e+03 | +6.336941e-05 | +2.259052e-01 | 0 | 0 | 1.216795599e+14 | 1.216795599e+14 | **0.0** |
| `weak` | +2.1658446e+01 | +2.991741e-07 | +6.479646e-06 | 0 | 0 | 9.505840969e+13 | 9.505840969e+13 | **0.0** |
| `reversed` | +4.2734479e+01 | +2.991741e-07 | +1.278505e-05 | 0 | 0 | 9.505840969e+13 | 9.505840969e+13 | **0.0** |
| `window` | +9.9898539e-01 | +2.994733e-07 | +2.991694e-07 | 0 | 0 | 9.505840969e+13 | 9.505840969e+13 | **0.0** |
| `fid_resolve` `rest` | (no window) | 0 | +6.048816e-06 | **0.5** | **0** | 9.158406986e+13 | 9.505840969e+13 | **+3.79361e-02** |
| `hydrostatic_column` `rest` | (no window) | 0 | +1.256476e-04 | **0.5** | **0** | 1.216795225e+14 | 1.216795599e+14 | **+3.07002e-07** |
| `rmax45` | **-5.8969125e+00** | +2.991451e-07 | **-1.764033e-06** | **0** | **1** | 9.506200157e+13 | 8.812646231e+13 | **-7.29581e-02** |

The ghost cells follow the face state. On the six rows that agree, all three
primitive components of both ghosts agree to 0.0. On the three that do not:

| state | ghost 0 rho | ghost 0 v | ghost 0 p | ghost 1 rho | ghost 1 v | ghost 1 p |
|---|---:|---:|---:|---:|---:|---:|
| `fid_resolve` `rest` | +3.864e-02 | -6.816e-04 | +1.142e-03 | +4.001e-02 | -1.993e-03 | +3.334e-03 |
| `hydrostatic_column` `rest` | +3.074e-07 | -3.864e-10 | +6.442e-10 | +3.082e-07 | -1.157e-09 | +1.929e-09 |
| `rmax45` | -7.422e-02 | +1.360e-03 | -2.271e-03 | -7.664e-02 | +3.984e-03 | -6.613e-03 |

(local minus remote, over remote.)

**Reading the table.** Three things are measured here and each says something
different.

1. **On both certified states the base face carries the wind's own mass flux**
   (0.99998437 and 0.99997483 of the window mean, and the whole face budget is
   flat to 1.26e-09 and 5.05e-08 of it). That reproduces what the review's
   probe measured and what item L27's certification budget prints, and it is
   the reason the L26 R3 numbers were withdrawn.
2. **The local discriminant selects the same branch as the remote one on every
   state whose base is moving, whatever the base is doing.** It does so at a
   very different magnitude: 21.7 and 42.7 times the wind's flux on the `weak`
   and `reversed` bases, 3565 times on the hydrostatic column. That magnitude
   is not the advective flux of the face state. On a base at Mach 1e-6 the
   HLLC mass flux is dominated by its dissipation term, `-(S_L S_R / (S_R -
   S_L))(rho_R - rho_L)`, i.e. by the DENSITY JUMP the boundary makes across
   the face, so `M_face` reads that jump and not the sign of `v`. It reads
   its sign correctly on these states, which is what a branch discriminant
   needs, but a memo that called it a local measurement of the face velocity
   would be wrong.
3. **The `window` row separates the two by construction.** A perturbation that
   scales only the wind density above `r_flux` by 1 + 1e-3 moves the remote
   discriminant by exactly +1.01e-03 and leaves the local one at its
   unperturbed value to the last bit. That is the nonlocality of the retained
   closure, measured.

### 4.4 The local row as a map: does it reproduce the flux it read

A local closure is a closed condition only if the ghost it implies puts back
through the face the flux it was closed on. MEASURED by installing the local
row's face state, ghost averages and lower face state in place of the
production boundary (writing `base_face_W`, `base_ghost_W`,
`base_face_lower_W` and the ghost cells, then reconstructing and assembling
the right-hand side under the stationary operator, with `Apply_BC` NOT called,
since it would put the production boundary back), and re-reading
`r_edg(0)^2 (rho v)_0`:

| state | `F_face` read from the production boundary | after the local ghost is installed | at iterations 3 to 6 |
|---|---:|---:|---:|
| `fid_resolve` | 6.3065423448e-07 | 6.3065423448e-07 | unchanged |
| `x010_HeH2.13` | 5.6225295289e-08 | 5.6225295289e-08 | unchanged |
| `weak` | 1.3659204153e-05 | 1.3659204153e-05 | unchanged |
| `reversed` | 2.6951101577e-05 | 2.6951101577e-05 | unchanged |
| `window` | 6.3065423448e-07 | 6.3065423448e-07 | unchanged |
| `hydrostatic_column` | 3.2133484083e-01 | 3.2133484083e-01 | unchanged |
| `fid_resolve` `rest` | 1.2751005420e-05 | 1.3790806400e-05 | unchanged |
| `hydrostatic_column` `rest` | 1.7872521301e-04 | 1.7872519031e-04 | unchanged |
| `rmax45` | -3.7190878694e-06 | -4.9211662680e-06 | unchanged |

**The local closure is a fixed point of its own map on every state tested, and
it is reached in one iteration.** Where the two closures agree the map does
not move the flux at all, because the installed ghost is the production one.
Where they differ the flux moves once, by 8.2 per cent on the fiducial at
rest, by 1.3e-07 on the hydrostatic column at rest and by 32 per cent on the
extended-domain seed, and then stands. In no case does the map change the
branch it selected, so it does not oscillate between inflow and reversal on
any state tested.

### 4.5 The retained remote closure: the window location

MEASURED. `r_flux` is the inner edge of the wind window and is stated by the
second field of `Flux spread tol:`; its production value is 1.20 R_p. The
probe resets `r_flux`, recomputes `j_flux` the way `define_grid` computes it
and re-closes the boundary. On the fiducial:

| `r_flux` [R_p] | `j_flux` | cells in the window | `F_window` | `d_window` | `M_wind` | `w_rev` | face rho [mH/cm3] | face v [cm/s] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.01 | 52 | 449 | 6.306928e-07 | 8.0961e-05 | 2.991877e-07 | 0 | 9.505841e+13 | -7.253964 |
| 1.05 | 139 | 362 | 6.306873e-07 | 8.7135e-05 | 2.991851e-07 | 0 | 9.505841e+13 | -7.253964 |
| 1.10 | 178 | 323 | 6.306760e-07 | 7.3885e-05 | 2.991798e-07 | 0 | 9.505841e+13 | -7.253964 |
| **1.20** | 217 | 284 | 6.306641e-07 | 5.5854e-05 | 2.991741e-07 | 0 | 9.505841e+13 | -7.253964 |
| 1.50 | 270 | 231 | 6.306546e-07 | 5.0352e-05 | 2.991696e-07 | 0 | 9.505841e+13 | -7.253964 |
| 2.00 | 309 | 192 | 6.306488e-07 | 5.0417e-05 | 2.991668e-07 | 0 | 9.505841e+13 | -7.253964 |
| 5.00 | 389 | 112 | 6.306292e-07 | 4.5232e-05 | 2.991576e-07 | 0 | 9.505841e+13 | -7.253964 |
| 10.00 | 435 | 66 | 6.306103e-07 | 3.3499e-05 | 2.991486e-07 | 0 | 9.505841e+13 | -7.253964 |

**On a converged wind the window location does not matter**: over a factor
ten in `r_flux` and a factor seven in the number of cells in the window, the
window mean moves by 1.3e-05 relative, `M_wind` with it, and the face state
is the same to every digit printed, because the weight is saturated
throughout. The same sweep on the 300-step hydrostatic column moves
`F_window` by six decades (5.62e-02 at 1.01 down to 3.85e-08 at 2.00) and
`d_window` from 1.58 to 2.40 and gives the SAME face state at every location,
because the shape weight `C` is zero at every one of them and the first
interior cell answers. On the state at rest the window carries exactly zero
at every location and the answer is again unchanged. So the window location
is not a free parameter the answer depends on, on any state tested; what the
answer depends on is whether the window says anything at all.

### 4.6 The retained remote closure: a domain extension

MEASURED on the fiducial mapped onto `Outer radius [R_p]: 45.0` (1.5 times
the production domain) at 500 cells, with the state extended above 28.0 R_p
by `map_state_to_grid.py --extrapolate-beyond`:

| | 30 R_p domain (the solved state) | 45 R_p domain (the mapped seed) | change |
|---|---:|---:|---:|
| `F_window` | 6.3066408976e-07 | 6.3068391650e-07 | +3.14e-05 |
| `M_wind` | 2.991741e-07 | 2.991451e-07 | -9.69e-05 |
| face rho [mH/cm3] | 9.505840969e+13 | 9.506200157e+13 | +3.78e-05 |
| face v [cm/s] | -7.253964 | -7.451617 | +2.72e-02 |
| `w_rev` | 0 | 0 | 0 |
| `F_face / F_window` | +0.99998 | **-5.8969** | sign |

**The remote closure is insensitive to the domain extension**: its
discriminant moves by 1e-04 and its face density by 4e-05, and the branch
does not change. What DOES change is the production base face flux, which
goes negative at -5.90 of the window mean. The state is a mapped seed and not
a solution of the extended domain, so a face flux that is not the wind's is
the expected reading; the point for this study is that the LOCAL closure
follows it into the reversal branch (`w_rev = 1`, face density 7.30 per cent
below the reservoir) while the remote one does not. On a state the operator
does not vanish on, the two closures are different boundary conditions.

---

## 5. What was not measured

- **A refinement of the base.** Section 3.1: on this grid generator the inner
  face radius and the first cell width are one parameter, so a matched-domain
  refinement leaves the base cells untouched. The base difference between the
  two fluxes is therefore not resolved by this study and cannot be by any
  study that holds `r_edg(0)` fixed on a `Mixed` grid.
- **A third flux, and any planet but this one.** One case, one planet, one
  XUV level.
- **Roe's behavior on the catalog's other 120 certified cases.** Nothing here
  re-solves any of them, and R37 still stands: adopting Roe means re-solving
  and re-judging every one of them on its own residual.
- **A local closure with its own derivation.** The local row of step 2 is the
  production compatibility solve with one argument replaced. A closure derived
  from the face Riemann problem itself, rather than keyed on the flux that
  problem returns, is a different object and was not built.
- **Whether either closure is right at a state at rest.** Section 4.3 measures
  that they differ there by 3.8 per cent on the fiducial. Which of the two is
  the physical answer is a question about the entropy of a stationary contact
  wave and is not settled by a measurement of the two candidates.
- **A certified refined state, at either flux.** All four refined runs were
  stopped at a stated outer-pass budget with their hydrodynamic energy row on
  the plateau of section 3.5. What their state would be after another thirty
  passes of elemental relaxation is not measured; what IS measured is that
  the energy row stops falling by the seventh pass and that the elemental row
  halves every pass, so a longer run relaxes the composition and not the row
  that refuses the state.
- **Whether the wind-cell energy row is a property of the grids, of the
  coupled solve, or of this binary.** Section 2: the snapshot carries item
  L30's work in progress in `RK_rhs.f90`, `hydrodynamic_rows_body.inc`,
  `steady_residual.f90`, `Reconstruction.f90` and `species_face_flux.f90`,
  which is the residual assembly the row is measured by. The 500-cell Roe
  control reproduces the catalog state's first cell to 5e-05 in T and 4e-05
  in rho and certifies at 4.233e-07 against L25's 4.143e-07, so the binary is
  not obviously moving this case; that is a check and not a separation.

---

## 6. The two recommendations, stated as recommendations

### Recommendation 1: keep Roe to the seven low-XUV cases, and do not make the grid the catalog's answer

Recommended, not decided. Three measured statements carry it.

1. **The grid does relieve the stall.** At 500 cells HLLC has no stationary
   root of this case: its mass and energy rows stand at 1.8e-01 and 2.4e-01
   at cell 1 and its base face carries 0.700 of the wind's mass flux. At 1000
   and at 2000 cells, on the SAME physical domain, the same unmodified HLLC
   flux brings the mass row to 2e-08, the momentum row to 3e-14 and the base
   face to 1.0000000 of the wind's flux at the first outer pass. So the
   500-cell stall is a resolution property of that flux, which is what L25
   concluded, now on matched physical boundaries.
2. **The grid does not remove the flux-family difference, and the difference
   is converged.** Between 1000 and 2000 cells each flux is converged to a
   few parts in 1e5 in its first cell and in its mass-loss rate, and the
   difference BETWEEN them is the same at both resolutions to three digits:
   10.3 per cent in T(1), 11.1 per cent in rho(1), 1.0 per cent in Mdot, 1.8
   per cent in the He I 10830 equivalent width, and 0.5 to 3 per cent in the
   density from 1.2 R_p to the top of the domain. Two converged discrete
   answers that differ by that much are a statement about the discretization
   at this base, not about resolution, and choosing between them is not
   something this study can do.
3. **Refining at fixed domain costs the certification.** None of the four
   refined runs certifies: all four carry a hydrodynamic energy row of 2.7e-04
   to 8.2e-04 against 1.0e-06 at r = 4.3 to 5.0 R_p, which appears as the
   elemental partition relaxes and then stops falling, at both fluxes. The
   only certified state of the six is Roe at 500 cells. A catalog re-run on a
   1000-cell grid would therefore, on this binary, replace 120 certified
   states by 120 uncertified ones.

So: the seven-case Roe scope remains the right one, and the answer to "should
the grid be the catalog's answer instead" is no as it stands. What would have
to be measured before either changes: (a) whether the wind-cell energy row of
statement 3 is a property of the refined grids, of the coupled solve, or of
the tree this study was built on (section 2 names the files item L30 has open
in the residual assembly), which is one evaluation of the same state on a
clean build; and (b) if it is not a defect, which of the two converged base
states is the physical one, which needs an independent reference (a
published solution or a second code), since neither the residual, the
conservation budget nor the grid trend separates them.

The measurement that this item CAN add to the R37 list, and that was missing
when R37 was written, is now available: the Roe flux's own grid convergence.
It is converged at the catalog's own 500 cells, to 2.5e-06 in T(1) and
1.1e-04 in Mdot at a doubling. That removes one of the two reasons the
handoff gave for recommending against a wider Roe adoption; the other, that
every certified case would have to be re-solved and re-judged on its own
residual, stands unchanged.

### Recommendation 2: one bounded source item, on the zero-window limit, not on the closure family

Recommended, not decided. The exclusion that L26 section R3 recorded does not
survive: it was taken on cell-centered readings and on PLM evaluations of
WENO3 states, and under the production face flux the local characteristic
closure selects the same branch as the remote-window one on every state whose
base is moving, installs a bit-identical ghost on all six of them, and is a
fixed point of its own map in one iteration. So a local closure is not
refuted. It is also not established: on those six states the two agree only
because the branch weight is saturated at both discriminants, so they are not
being told apart.

Where they do differ is one identified case, and it is the one worth a source
item. **On a column at rest the wind window carries exactly zero flux, the
amplitude weight vanishes, and the production boundary falls back on the
first interior cell, whose face Mach number is also exactly zero, so the
smoothstep returns 1/2 and the base density becomes the average of the
reservoir isentrope and the interior isentrope.** On the fiducial those two
differ by 7.3 per cent, so the boundary states a face density 3.79 per cent
below the reservoir on a state where nothing is flowing, and the ghosts
follow at 3.86 and 4.00 per cent. The production face flux on the same state
is not zero (it is the HLLC dissipation term reading the density jump the
boundary itself makes, +6.0e-06 in face Mach number), so a face-flux
discriminant returns the reservoir there instead. On the hydrostatic column,
where the two isentropes agree to 6e-07, the same 1/2 costs 3.1e-07, which is
why this has not been visible.

The bounded item that follows is therefore not "replace the closure" but
"state what the boundary does when the wind says nothing": at `A = 0` the
present code hands the decision to a discriminant that is identically zero,
and a weight of 1/2 there is a choice, not a limit. Two candidates are the
reservoir (the entropy of a stationary contact wave at the base of a column
at rest is the reservoir's, since nothing has advected out of the interior)
and the present average. Deciding between them is a physics question about a
stationary contact and is not settled by this memo; measuring the cost is,
and it is 3.8 per cent of the base density on this planet and 3.1e-07 on a
well-balanced column. Whatever is chosen, the smoothstep's value at zero
should be stated in the code at the place it is used, which it is not now.

Two further things this study measured and that any such item should carry:

- **The remote closure's nonlocality is real but its answers are not
  sensitive to the window.** Over a factor ten in `r_flux`, on a converged
  wind, the window mean moves by 1.3e-05 and the face state does not move at
  all; on a 1.5x extended domain the discriminant moves by 1e-04 and the face
  density by 3.8e-05. The window location is not a tuning parameter the
  answer depends on.
- **The production base face flux is not the sign of the face velocity.** At
  a base Mach number of 1e-06 the HLLC mass flux at that face is dominated by
  its dissipation term, so `M_face` reads the density jump the boundary makes
  and not the advective direction. It reads the right sign on every state
  tested, which is all a branch discriminant needs, but a closure that called
  it a local measurement of the flow direction would be misdescribing it.

---

## 7. Reproduce

The run directories are `LHS1140b/models/.L35/`:

- `EXHALE_L35.x`: the binary of section 2, md5
  `16f5345d07ce70cdac09227839d5188b`, copied here so that both hosts read one
  file.
- `grid/g{500,1000,2000}/`: the one-step probes that produced the three grid
  files, `Load IC? False`, `EXHALE_MAXSTEPS=1`. `grid_replica.py` is the
  Python replica of `define_grid`'s `Mixed` branch and `solve_grid.py` the
  two-parameter solve of section 3.2; running the latter reproduces the three
  key lines.
- `g{500,1000,2000}_{hllc,roe}/`: the six solves, each with `input.inp`,
  `seed.log`, `run.log`, `pp.log`, `output/`, `transit.log` and the five
  `tpm_*.txt`. `run_passes40_stopped.log` beside `run.log` is the first
  attempt's log, the one with the 40-pass budget, kept for its pass records;
  `.L35/<case>.launch_passes40_stopped.log` is that attempt's runner output
  and `.L35/<case>.run2.log` the re-run's. `g500_roe` has no such pair: it
  certified at its first outer pass and was never re-run, and its
  `REPRODUCE.md` is the runner's own.
- `launch_all.sh`: the first attempt's launcher, with the environment of
  every run in it. `solve_and_post.sh <dir> <passes> <threads>`: the re-run,
  the catalog recipe with the runner's refusal to make products for a
  non-certifying solve removed.
- `collect.py`, `table.py`, `profiles.py`, `make_figure.py`: the tables of
  sections 3.4 and 3.5 and the figure of section 3.6, from the run
  directories.

No catalog directory, no other `.L*` directory and no golden was written to.
The states of step 2 are scratch copies under
`<scratchpad>/L35/states/`, made from `.L26/fid_resolve`, the catalog
`atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13` and
`backup/regression/hydrostatic_column`, with `Load IC? True` and the
`Restart intent:` line removed; the sources were read and never written to.

---

## 8. Noticed outside this item, reported and not fixed

**`dr_base`'s default is a single-precision literal.**
`src/modules/init/parameters.f90` line 143 declares

```fortran
      real*8  :: dr_base     = 2.0e-4   ! uniform base cell size [R_p]
```

and `2.0e-4` without a `d` exponent is a DEFAULT REAL literal, so the
real\*8 variable is initialized to `1.9999999494757503e-04` and not to
`2.0e-4`. `input_read.f90` line 1176 reads the key list-directed into the
same real\*8, so `Base grid [dr,cells]: 2.0e-4 50` in `input.inp` gives
`2.0000000000000001e-04`: **stating the default explicitly produces a
different grid from omitting the key.** MEASURED: the two grids' cell centers
differ by up to 6.58e-09 relative (largest near cell 343), which is the
2.53e-08 difference in `dr_base` damped by a factor 3.8 by the rescale onto
the fixed `r_max`. This is how the
replica of section 3.1 was validated, and it is why the 500-cell row of the
table in section 3.2 states 17 digits: with `2.0e-4` in the key the run would
not have been on the catalog grid. The effect is far below any physical
tolerance and the fix (`2.0d-4`) moves every grid in the tree by 6.6e-09,
which is a golden refresh, so it is reported and not made here. A sweep of
`parameters.f90` for other real\*8 defaults written with an `e` exponent
would be the right scope for it.

**The `hydrostatic_column` regression fixture is not a column at rest.** Its
recorded output is a 300-step relaxation snapshot (`maxsteps` says 300) whose
first cell carries an inflow of 8.09e+04 cm/s at `M_i = +0.2257`, and whose
wind window carries a relative flux spread of 2.40, the largest in the tree.
Its `README.md` states what a pass would be and says the Phase 4 test is not
in that directory, so nothing in the tree claims otherwise; but a reader who
takes the fixture's name for its state will be wrong, and a future well-
balancedness test should build its own state at rest rather than load this
one. Reported, not changed: the fixture pins a regression output and changing
it is a golden refresh.
