# L18: why a certified state did not reproduce its certification from its own files

Item L18 of `docs/PLAN_20260913_lhs_stationary.md`, opened on the side finding
of `docs/lhs1140b_stationary_L17_20260915.md` section 9 and approved by the
user on 2026-09-15.

The measurement it starts from: `LHS1140b/models/.L14/x003_HeH2.13` was written
at outer pass 19 with its hydrodynamic energy row at 8.255e-07 (cell 7), within
its 1.0e-06, and re-entered with `Restart intent: stationary evaluate` it
measured 1.623e-06 (cell 20) and the run exited 2.

## 1. Verdict

**The state the restart evaluated was not the state that was certified, and the
difference was in the mass density.** The pair of files states the density
twice -- as the `rho` column of `Hydro_ioniz.txt`, and as the species densities
of `Ion_species.txt`, which weigh `sum_i m_i n_i`. The two agree only while the
composition closes its own mass, `sum_i f_i A_i = 1`, which is the definition of
`f_sp` and not a tolerance. **It does not close.** The composition of a state
written at the end of a stationary solve departs from its own density by 5.1e-13
on this case, and the loader took the species as the authority and rebuilt the
density from them, so the state it handed the evaluation stood some 2000 units in the
last place from the one the solve had certified. The hydrodynamic rows of this
subsonic base cancel their largest term by 6.2e+05, and at that cancellation
some 2000 ulp of density is the whole row.

Two things follow, and both are measured below.

1. **The reconstruction rule was wrong, and it is fixed.** The conserved
   variable is the mass density: it is what the hydrodynamics advances and what
   the residual is a function of, while the composition is an eliminated
   variable the first sweep re-solves anyway. The loader now takes the density
   from its own column and projects the loaded species onto it, one factor per
   cell, which leaves every element ratio and every ionization split where the
   file put them. With that, `rho` round-trips to the last bit and this state
   re-reads at **mass 3.764e-08 at cell 7, momentum 1.814e-14 at cell 1, energy
   8.255e-07 at cell 7 -- CERTIFIED, exit 0**, against the 3.764e-08 at cell 7,
   1.960e-14 at cell 1 and 8.255e-07 at cell 7 it was written at.

2. **Bitwise reproduction is not available, and the contract has to say so.**
   One unit in the last place of the density moves the cell-wise maximum of the
   energy row of this state by about 11 per cent (measured, section 3), because
   that row is a cancellation of its largest term by 6.2e+05 and the flux
   assembly that forms it is at its rounding floor there. A state whose row
   stands within a factor two of its tolerance can therefore be certified in the
   run and refused on re-evaluation, or the reverse, with no defect anywhere.
   Section 8 states what the contract can promise instead, with the numbers.

## 2. The mass closure of a written state, measured

`sum_i f_i A_i - 1` of the state pair, computed from the files alone: the
species columns are `f_i rho n0` and the density column is `rho n0`, so their
ratio is exactly that sum (`calc_rho`, `utilities.f90`, is the mass policy both
sides use; the metals enter at `melem_A` and He 2^3S is inside the He I column).

| states | n | median | worst |
|---|---|---|---|
| the state pairs of `LHS1140b/models/` (the production catalog) | 136 | 7.8e-16 | 3.8e-14 |
| `.L14` and `.ab`, states written at the end of a stationary solve | 11 | 9.2e-14 | 6.2e-13 |
| written by a stationary EVALUATION (`.L17`, `.L18`) | 20 | 9.8e-16 | 1.4e-15 |

It grows with the number of outer passes and it grows monotonically:

| case | outer passes | `sum_i f_i A_i - 1` |
|---|---|---|
| `.ab/step_1p55` | 2 | 2.5e-14 |
| `.L14/c` | 4 | 3.4e-14 |
| `.ab/new_2p13` | 4 | 3.8e-14 |
| `.ab/arch_1p50` | 7 | 9.1e-14 |
| `.L14/x005_HeH2.13` | 20 | 1.5e-13 |
| `.L14/x005d8_HeH9.7` | 29 | 1.9e-13 |
| `.L14/x002_HeH2.13` | 40 | 2.3e-13 |
| `.L14/x003_HeH2.13` | 19 | 5.1e-13 |
| `.L14/x003_HeH9.7` | 24 | 6.2e-13 |

about 1e-14 per outer pass. Traced with `EXHALE_ELEMENT_ASSERT=2`, which reports
the closure around every `ioniz_eq` sweep, one pass of `.L18/census` takes it
from 3.34e-16 at entry to 1.17e-15 at the state written, in steps of about
2e-16 per sweep and with no single call responsible: the sweep conserves the
element TOTALS to 4.4e-16 and the mass is a weighted sum of them, so the closure
random-walks and ratchets. The element operator's own test
(`element_mass_closure_tol = 1e-10`, `binary_element_diffusion.f90`) compares a
step's closure with the ENTRY closure of that step, four decades above the
increment and blind to a monotone drift by construction. **Nothing in the code
measures the closure of a state.** This is left open as item L19 below; it is
not what the restart contract needed, because a restart that keeps the conserved
density is exact whatever the closure is.

## 3. The fixture, and what the round trip carries

`.L14/x003_HeH2.13` re-entered with `Restart intent: stationary evaluate`, one
thread, control binary `EXHALE_L18_ctl.x` (md5 `922dec0f`, the tree as L17 left
it). The rows reproduce the side finding exactly: mass 3.444e-08 (cell 20),
momentum 4.495e-13 (cell 19), energy **1.623e-06 (cell 20), ABOVE**, exit 2.

**What the file carries and what the run rebuilds, column by column**, comparing
the file read in with the file written back:

| column | is it read? | max relative departure over the 504 rows |
|---|---|---|
| `r` | compared with the grid | 0 (bitwise, all rows) |
| `v` | read | 0 in every physical cell (4.1e-12 in the two lower ghosts, which `Apply_BC` rewrites) |
| `p` | read | 0 in 501 of 504 rows, 1.0e-15 worst |
| `T` | read, then recomputed from `p` and the particle count | 1.8e-14 |
| `heat` | NOT read; rebuilt by `ioniz_eq` | 2.3e-13 |
| `cool` | NOT read; rebuilt by `ioniz_eq` | 3.3e-13 |
| `rho` | NOT read; rebuilt by `calc_rho` from the species | **5.1e-13, and 0 rows bitwise** |

So the radiative terms of this state DO rebuild, to 3e-13; it is the density
that does not, and it does not because the loader was rebuilding it from the
other half of the file.

**The energy row's terms at the cell that binds** (`write_residual_breakdown`):
flux divergence 4.353857e-06, heat 4.406056e-06, cool 5.219202e-08, source 0,
so `R = -7.151968e-12` against a largest term of 4.406056e-06 -- the row is a
cancellation by 6.2e+05. A relative perturbation of 1.6e-06 in either surviving
term is the whole row; and the flux divergence has a cancellation of its own
(the base-layer face states cancel to 3e5-6e6, L17 section 4), so a density
perturbation of 5e-13 reaches the row multiplied by that.

**Window by window**, the volume-integrated numerators of the same evaluation,
control against fixed. The denominators -- the row scales -- are identical to
every printed digit in the two runs, so this is the residual moving and not the
state:

| row, window | numerator, control | numerator, fixed |
|---|---|---|
| mass, r < 1.03 | 1.189e-14 | 9.385e-15 |
| mass, 1.03-1.10 | 4.699e-16 | 4.591e-16 |
| mass, r >= 1.20 | 1.245e-17 | 1.269e-17 |
| momentum, r < 1.03 | 6.884e-14 | 2.930e-15 |
| momentum, r >= 1.20 | 8.165e-17 | 1.274e-17 |
| energy, r < 1.03 | 3.272e-14 | 1.991e-14 |
| energy, 1.03-1.10 | 2.250e-15 | 2.200e-15 |
| energy, r >= 1.20 | 2.430e-16 | 2.441e-16 |

The departure is the base layer. Above 1.03 the two evaluations agree to a few
per cent and above 1.20 to the third digit.

**The sensitivity, measured directly.** Five evaluations in a chain, each of the
state the previous one wrote. Every state in the chain closes its own mass to
5.7e-16, so what separates consecutive states is a unit or two in the last place
of the density:

| step | max relative change of the density evaluated | energy row |
|---|---|---|
| 1 | -- | 1.6232e-06 |
| 1 -> 2 | 1.02e-14 | 1.4586e-06 |
| 2 -> 3 | 3.84e-16 | 1.6215e-06 |
| 3 -> 4 | 3.75e-16 | 1.4586e-06 |
| 4 -> 5 | 5.17e-16 | 1.3175e-06 |

**One ulp of the density moves the cell-wise maximum of the energy row by 10 to
11 per cent.** That is the reproduction limit of this row on this state, and no
change to the file format reaches below it.

(`EXHALE_RESID_QUAD=1`, the quadruple-precision control of
`hydrodynamic_rows.f90`, does not arm on this route: `assemble_residual` selects
it only when `ieq_sweep_state_kind` is not `ieq_state_marching`, and the
stationary evaluation calls `ioniz_eq` directly, which leaves the default tag.
The five runs above repeated with it set are bitwise the five without it. Noted
here because the knob looks available on this route and is not.)

## 4. What changed in the source

`src/modules/files_IO/load_IC.f90`, one block.

- The `rho` column of `Hydro_ioniz_IC.txt` was read into a discard variable and
  thrown away. It is now read.
- After `calc_rho` weighs the loaded species, the two statements of the density
  are compared, cell by cell over the physical column, and the largest departure
  is reported on every restart -- it was never visible before.
- Below `restart_density_agreement_tol = 1.0d-8` the conserved density is taken
  from its own column and the loaded species are multiplied, cell by cell, by
  the single factor that puts them on it. Every element ratio, every ionization
  split and every metal-to-hydrogen ratio is unchanged by a common factor, and
  the composition then closes the density it is a composition of.
- Above that the species stay the authority, as before, and the run says why:
  the blocks above the reconstruction change the composition on purpose -- an
  element the file does not carry is rebuilt at its abundance, a reservoir the
  handoff moved is rescaled, the oxygen carriers of a pre-oxygen-chemistry file
  are seeded, a molecular seed is taken from an atomic state -- and those move
  the mass by 1e-3 and more. The threshold separates two scales that are five
  decades apart on either side; it is not a physical tolerance.
- The comment that claimed "a restart therefore preserves the conserved mass
  exactly" is gone. It was false whenever the closure was broken, which is
  every state written by a stationary solve.

Also in this item, unrelated to it: `binary_element_diffusion.f90` line 861 used
"pipeline" for the run. Replaced.

## 5. The fixture, after

Same fixture, same command, binary `EXHALE_L18.x` (md5 `2b1458a9`):

```
 (load_IC) the species columns weigh the density column of the restart to
   5.146E-13 (worst at cell 36); the conserved density is taken from its own
   column and the composition is projected onto it.
   the loaded composition against the sweep's own root:
   max |d(n_tot+n_e)|/(n_tot+n_e) = 2.100E-14
   hydrodynamic mass row      max= 3.764E-08  cell=7   tol= 2.0E-08  within
   hydrodynamic momentum row  max= 1.814E-14  cell=1   tol= 1.0E-08  within
   hydrodynamic energy row    max= 8.255E-07  cell=7   tol= 1.0E-06  within
   CERTIFIED: every active equation was evaluated and is within its tolerance
```

against the certification it was written with: mass 3.764E-08 at cell 7,
momentum 1.960E-14 at cell 1, energy 8.255E-07 at cell 7. The mass and energy
rows are reproduced at their own cells to every printed digit; the momentum row,
five decades below its tolerance, differs by 7 per cent.

## 6. Every certified state of `LHS1140b/models/`, re-evaluated through its own files

All 83 `REPRODUCE.md` of `LHS1140b/models/` carry `info = 0`. Nine of them are
the closure iterations of the photochemical `kzzprofile` group, whose states are
the `k00`-`k05` subdirectories, so the states are 74 + 46 = **120**. Each was
copied out of `models/` (the catalog is read-only), given
`Restart intent: stationary evaluate`, and run on one thread with each binary.

**120 of 120 re-certify, with both binaries** (exit 0). The production states
close their own mass to 1e-15, so the reconstruction rule moved them by about
one ulp and the fix changes nothing there; their rows stand two to five decades
below their tolerances and survive the rounding-floor spread easily.

How closely the rows reproduce, as the ratio of the re-evaluated cell-wise
maximum to the in-run one over the 120 states:

| row | median | 10th | 90th | worst |
|---|---|---|---|---|
| mass | 1.12 | 0.78 | 1.71 | 3.45 |
| momentum | 9.2 | 1.80 | 134 | 4.5e+03 |
| energy | 1.22 | 0.89 | 1.93 | 3.34 |

The momentum column is the rounding floor seen bare: those rows sit at 1e-14 to
1e-12 against a tolerance of 1e-08, so the ratio of two draws from the floor is
whatever it is and says nothing about the state. The mass and energy rows
reproduce within a factor two.

The states where it matters are the ones whose closure has ratcheted and whose
rows stand near their tolerance -- the XUV ladder:

| state | row | in-run | control | fixed | tolerance |
|---|---|---|---|---|---|
| `.L14/x003_HeH2.13` | mass | 3.764e-08 @7 | 3.444e-08 @20 | **3.764e-08 @7** | 2.0e-08 |
| | momentum | 1.960e-14 @1 | 4.495e-13 @19 | **1.814e-14 @1** | 1.0e-08 |
| | energy | 8.255e-07 @7 | 1.623e-06 @20 | **8.255e-07 @7** | 1.0e-06 |
| | verdict | CERTIFIED | exit 2 | **exit 0** | |
| `.L14/x003_HeH9.7` | mass | 3.708e-08 @3 | 7.368e-08 @6 | **3.708e-08 @3** | 1.7e-09 |
| | momentum | 1.682e-14 @57 | 5.466e-13 @21 | **1.682e-14 @57** | 1.0e-08 |
| | energy | 8.670e-07 @19 | 2.581e-06 @6 | 1.176e-06 @1 | 1.0e-06 |
| | verdict | CERTIFIED | exit 2 | exit 2 | |
| `.L14/x005_HeH2.13` | energy | 5.441e-07 @8 | 6.417e-07 @45 | **5.441e-07 @8** | 1.0e-06 |
| | verdict | CERTIFIED | exit 0 | exit 0 | |
| `.L14/x005d8_HeH9.7` | energy | 4.422e-07 @40 | 1.807e-06 @1 | 6.134e-07 @1 | 1.0e-06 |
| | verdict | CERTIFIED | exit 2 | **exit 0** | |
| `.ab/arch_1p50` | energy | 1.687e-08 @33 | 3.407e-08 @16 | **1.687e-08 @33** | 1.0e-06 |
| `.ab/new_2p13` | energy | 2.880e-08 @22 | 3.308e-08 @22 | **2.880e-08 @22** | 1.0e-06 |
| `.ab/step_1p55` | energy | 4.807e-08 @1 | 4.290e-08 @1 | 4.879e-08 @1 | 1.0e-06 |

Five of the seven now reproduce their in-run row to every printed digit, and two
of the three that the control refused now certify. **`x003_HeH9.7` still does
not**: its energy row re-reads at 1.176e-06 at cell 1 against 8.670e-07 at cell
19, 18 per cent over its tolerance. Its mass and momentum rows are reproduced
exactly, so its density does round-trip; what is left is the rounding floor of
section 3 acting on a row that was certified at 87 per cent of its tolerance.
That is the case the contract statement of section 8 is about.

## 7. What is still not reproduced, and is not this

**A pinned fixture from an older code generation.** The `grid_and_gates`
restart row reads `backup/regression/wasp_full_newton/IC/`, and that state,
re-entered with `stationary evaluate`, reproduces `r`, `rho` (bitwise, after
this item), `v`, `p` and `heat` -- while its `cool` column comes back at 0.35
to 0.76 of what the file carries, everywhere on the grid. The energy row of
that state therefore re-reads at **3.973e-01** against its 1.0e-06, four
decades out, and its `certified=T` header is refused; the three assertions of
`grid_and_gates` that read it fail with the control binary and with the fixed
one alike.

**It is not a defect of the code that writes states today**, and item L20
established that by measurement (`docs/lhs1140b_stationary_L19_L20_20260915.md`).
That pinned pair was written on **2026-09-08 by git `35d9dd5d3ca7`**, an
earlier generation, and `IC/README.md` says it is refreshed only together with
the golden. Every state the CURRENT code writes carries a `cool` column that is
the cooling of the state written beside it: over the golden matrix the column
and the total of `Cooling_breakdown.txt`, which recomputes every channel from
the written (T, rho, f_sp), agree to between 3.4e-15 and 3.0e-10 -- including
the converged `wasp_full_newton` solution in the same case directory, which
agrees to 2.98e-10 where its own pinned seed is out by 1.8x.

So what section 6 of this memo could not do for `.L14/x003_HeH9.7` is one
thing, and this is another: the fixture is stale, refreshing it with the golden
is the action, and item L20 added the assertion
(`cool_column_matches_breakdown_total`) that would have caught it -- the same
statement had been asserted for `heat` since the closure-lag incident and never
for `cool`.

## 8. What the contract can promise

`docs/restart_contract_design_20260909.md` and every `REPRODUCE.md` say a
written state can be re-entered and re-measured. After this item the promise
that holds, with the numbers behind it, is:

- **The conserved state round-trips.** `r`, `v` and `p` are carried in 17
  significant figures and read back bitwise; `rho` is carried and, since this
  item, read back bitwise. What the run evaluates is the state the file names.
- **The composition round-trips to the file's own closure**, which is 1e-15 on a
  production state and up to 6e-13 on a long ladder solve; the projection puts
  that departure on the composition, where no tolerance in the inventory is
  within nine decades of it, instead of on the density, where it was the whole
  hydrodynamic row.
- **The row measures do NOT round-trip bitwise, and cannot.** The first
  equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11 (it
  is one Picard step of a nonlocal coupling, and the state it is handed is not
  exactly its own fixed point), and the flux assembly of a subsonic base is at
  its rounding floor. Measured over the 120 certified states of
  `LHS1140b/models/`: the mass row re-reads within a factor 0.78-1.71 of its
  in-run value (median 1.12), the energy row within 0.89-1.93 (median 1.22), and
  the momentum row, which sits four to six decades below its tolerance, within
  whatever the floor gives. On one state the response to a single ulp of density
  was measured directly and is 11 per cent of the energy row.
- **Therefore a state certified at more than about half of its tolerance may be
  refused on re-evaluation, and the refusal is not a defect of the file.** Of
  the 127 certified states re-evaluated here, one is refused on that ground
  (`.L14/x003_HeH9.7`, certified at 0.87 of tolerance, re-reading at 1.18).

## 9. Reproduce

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
cd $EX && make OBJDIR=build_L18 EXE=EXHALE_L18.x
M=$EX/LHS1140b/models;  L=$M/.L18

# the fixture, before and after (the control binary is the tree as L17 left it)
mkdir -p $L/x003_eval/output
cp $M/.L14/x003_HeH2.13/input.inp              $L/x003_eval/
cp $M/.L14/x003_HeH2.13/output/Hydro_ioniz.txt $L/x003_eval/output/Hydro_ioniz_IC.txt
cp $M/.L14/x003_HeH2.13/output/Ion_species.txt $L/x003_eval/output/Ion_species_IC.txt
sed -i 's/^Restart intent:.*/Restart intent: stationary evaluate/' $L/x003_eval/input.inp
cd $L/x003_eval && OMP_NUM_THREADS=1 $EX/EXHALE_L18.x

# section 2, the closure of a written state (files only, no run)
python3 - <<'PY'
import numpy as np
mHe = 6.6464790722e-24/1.67353284e-24
h = np.loadtxt('output/Hydro_ioniz_IC.txt'); s = np.loadtxt('output/Ion_species_IC.txt')
w = s[:,1] + s[:,2] + mHe*(s[:,3] + s[:,4] + s[:,5])       # metal-free case
print(np.max(np.abs(w - h[:,1])/h[:,1]))
PY

# section 2, where the closure drifts
OMP_NUM_THREADS=1 EXHALE_ELEMENT_ASSERT=2 EXHALE_ELEMENT_ASSERT_TOL=1.0e-3 \
  EXHALE_OUTER_PASSES=1 EXHALE_JFNK_MAXIT=5 EXHALE_PTC_DTAU0=1.0 $EX/EXHALE_L18.x

# section 3, the sensitivity: chain the evaluation, each run on the last one's output
# section 6, every certified state (the catalog is read from, never written)
```

## 10. Tests and regression

All from a copy; `backup/regression/` and the `models/` catalog were read and
never written. Both binaries are the same source apart from this item, so every
line below is control against fixed.

**Assertion suites** (`EXHALE_OBJDIR=build_L18`):

| suite | control | fixed |
|---|---|---|
| `certification` | 84 PASS, 0 FAIL | 84 PASS, 0 FAIL |
| `grid_and_gates` | 195 PASS, 4 FAIL | 195 PASS, **3** FAIL |
| `run_mode` | 29 PASS, 1 FAIL | **31 PASS, 0 FAIL** (see below) |
| `residual_determinism` | not rerun: the residual is a function of its argument and this item does not touch the residual |

- `grid_and_gates` / `output_state_consistency`:
  `outer_iteration_ending_is_the_stagnation_one` FAILS with the control
  (`measured=pass_budget reference=no_progress`) and PASSES with the fixed
  binary.
- `grid_and_gates` / `restart_intent_and_metadata`:
  `stationary_evaluate_mass_density` goes from 4.823e-14 to **0.000e+00** --
  the density round-trips to the last bit. The three that stay are section 7's
  pinned fixture, identical in both binaries and not this item's.
- `run_mode`: `a_phys_header_without_its_clock_is_refused` failed in both, and
  it was the TEST that was stale. It stripped a ` t_phys=` token from the
  `# coupling:` line alone, while the clock is written twice -- there and as
  the `# t_phys[s]` line of the restart metadata block, which the loader reads
  when the first is absent -- so the state reached the binary with its clock
  intact and was accepted. The assertion now strips both, and the metadata line
  from both halves of the pair (a block present in one half and absent in the
  other is refused on that ground instead, which would have made the assertion
  pass while testing something else), and it asserts the refusal names the
  clock. The code was right throughout: it refuses with "its t_phys field is
  missing or is not a finite non-negative number". `run_mode` is now 31 PASS,
  0 FAIL.

**Regression**, from a copy of `backup/regression/`, one thread.

The two cases the item names, `wasp_full_newton` and `atomic_elem_newton`,
are COLD STARTS: both carry `Load IC? False`, both begin at
`(set_IC.f90) Using b0_eff = 3.0 for IC`, and neither enters `load_IC` at all,
so neither can be reached by this change. (Their `IC/` directories are what the
restart assertions of `grid_and_gates` restart FROM, and that pair was run
control against fixed above, which is the reload those two fixtures exercise.)
Both were started anyway and both were stopped: the golden run of
`atomic_elem_newton` takes 70562 marching steps, which is about eighteen hours
at one thread and would have to be paid twice, and `backup/regression/golden/`
holds no reference for that case at all, so the run could only end in
`MISS`; `wasp_full_newton` reached step 3600 of the 4709 its golden run took.
The measurement is the three cases that DO load a state and carry a golden:
`hp_front`, `hp_trace_seed` and `hp_zero_seed`.

The control binary is byte-identical to itself run to run on all three
(`REGRESSION_REL_TOL=0`, 12 of 12 files "data identical"), so a difference
against a control snapshot is the change and not the machine. Against that
snapshot, strictly:

| case | `Hydro_ioniz` | `Ion_species` | `Hydro_ioniz_adv` |
|---|---|---|---|
| `hp_front` | 1.1e-11 | 5.7e-13 | 2.4e-08 |
| `hp_trace_seed` | 2.1e-04 | 3.2e-04 | 2.1e-02 |
| `hp_zero_seed` | 1.9e-04 | 3.0e-04 | 1.7e-02 |

**All three IC pairs close their own mass to 4.862e-16**, one ulp, so what the
fixed binary does to them is move the loaded density by one ulp and nothing
else. That the outcome then moves by 1e-4 is the fixture, and it was measured
directly: with the CONTROL binary, multiplying every species density of
`hp_trace_seed`'s IC by exactly `1 + 2^-52` and marching the same 100 steps
moves the state by **3.76e-03** (the cooling column; 1.57e-05 in the velocity)
-- more than this change does. These two fixtures are the molecular-basin seeds
whose cell solve is known bistable from a zero seed, and 100 marching steps at
a front amplify one ulp by twelve decades. `hp_front`, seeded at the front
itself, moves by 1e-11 for the same one ulp.

By the standing rule the two state files are identical (0.02 to 0.03 per cent,
below 0.1 per cent); the `_adv` profiles, which correct a 100-step relaxation
snapshot, are not. **No golden was refreshed**: the reference of record was
read and never written, and whether the two seed fixtures should carry a golden
at all, given that they amplify one ulp to 1e-3, is a question for the user and
not for this item.


---

## 11. The Codex review of 2026-09-15, applied

**Which half of the file states the density is decided by INTENT, not by the
size of the disagreement.** The first form of this item read the measured
departure against a single threshold, and a magnitude cannot tell a deliberate
change from a damaged file: a mismatched pair would have been accepted as a
deliberate one.

The loader knows its own intent. Each block that changes the composition on
purpose -- the loaded H/He carried onto the input's He/H, an element the file
does not carry rebuilt at its abundance, a reservoir the handoff moved, the
oxygen carriers of a pre-oxygen-chemistry file seeded -- now says so
(`composition_changed_here`, with the block that did it named in the message),
and the flags are set on the FACTORS those blocks apply, not on the fact that a
block was entered: the He/H branch runs on every diffused restart and applies
factors that are 1 to round-off, which is not a change of composition. A
molecular seed is deliberately not in the list: it reads an ATOMIC pair here,
which is self-consistent, and builds the molecular state afterwards.

Where a block acted, the density follows the composition, as before. Where none
did, the conserved density is the authority and the pair is allowed only
rounding:

| departure | what happens |
|---|---|
| <= `restart_density_rounding_tol` = 1e-10 | the density is taken from its column and the composition projected onto it |
| 1e-10 to `restart_density_agreement_tol` = 1e-8 | the same, and the departure is reported as above what a restart of the same equations should carry |
| > 1e-8, and the pair carries a restart metadata block | **REFUSED**: "the two halves of this restart describe two different gases", the pair is damaged or mismatched, and loading it would choose silently between them |
| > 1e-8, and the pair carries no metadata block | loaded with the composition as the authority and the departure reported: the file predates the block, its density column is another generation's number under a mass policy this loader cannot check, and everything the run writes is already marked `provenance_unknown` |

**The last row is not a softening, it is what the measurement forced.** A hard
refusal above 1e-8 was written first and it rejected, on its first run, every
legacy pair in the tree: all forty sampled pairs of
`LHS1140b/archive_20260830` disagree by **1.4e-03 to 7.1e-03**, and so do the
pinned fixtures `backup/regression/atomic_elem_newton/IC` (**1.288e-03**, and
it is not an `eos_metals` mismatch -- leaving the metal mass out makes it
1.1e-02, not smaller) and `wasp_full_newton/IC`. All were written on
2026-09-08/09 by git `35d9dd5d3ca7`. Those are seeds the campaign is entitled
to use and fixtures the test suite owns, so a pair that makes no claim about
its own provenance is reported and loaded, and only a pair that claims to be
this code's own is refused.

On the certified LHS 1140 b state the departure is 7.1e-16 and the conserved
density is taken from its column, as before; the answer is unchanged bit for
bit (`docs/lhs1140b_stationary_L21_20260915.md` section 10.4).
