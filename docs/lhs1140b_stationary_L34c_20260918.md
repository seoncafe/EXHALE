# L34 step (c): the molecular cases the reference solution did not reach, and what the restart contract allows

Item L34 of `docs/PLAN_20260917.md`, part (c): the three molecular cases
L34b left unsolved (`molecular_scalar_gj1132_kzz1e9/HeH0.55` and the two
`molecular_photochem_gj1132_kzzprofile` cases) and the three
`molecular_scalar_gj1132_wellmixed` cases L34b did not start. The reference
solution itself, `molecular_scalar_gj1132_kzz1e9/HeH2.13`, and its companion
`HeH9.7` are L34b's and are quoted here only as the seeds and the controls
they are.

Every number below is MEASURED on this tree unless marked READ. The binary is
`LHS1140b/models/EXHALE_3146d11b.x`, md5
`3146d11b4090306dcea75bb9718edd22`, manifest
`LHS1140b/models/BINARY_MANIFEST_3146d11b4090.txt` (it carries L27, L28, L30,
L30b and L31). No source file was changed, no tolerance was weakened, and no
key that defaults off was turned on beyond the recipe: `Well balanced: True`,
`Secondary_ionization: Immediate`, `Restart intent: stationary`,
`Coupled carrier solve: False`, `EXHALE_PTC_DTAU0=1.0`.

**The uncertainty every molecular number below carries, stated once.** At the
base of a molecular column the chemical heat is set by the helium third body,
whose coefficient is argon's used as a helium bound, and +/-0.3 dex of it
moves the base chemical heat by +42 and -21 per cent (MEASURED by L34b on the
certified reference solution's own cells; the same bracket on L31's archived
cells is +41 / -21, `docs/lhs1140b_stationary_L31_energy_cycles_20260917.md`
section 5.3). Every base temperature, base H2 fraction and base heating rate
quoted here stands inside that bracket, and it is larger than any difference
between the states compared below.

---

## 0. Verdict

**None of the six cases certified, and the reason is now a measured property
of the equations and of the alternation rather than a guess about seeds.**
Four of the six were seeded better than any previous attempt had been, one of
them from a state solved for the purpose; the wind converges in four of the
six and the H2 carrier balance is what refuses in every one of them. Two
seeding routes the plan offered are inadmissible by the restart contract, and
one of those makes the L33 ladder's element-treatment step impossible on a
molecular state at all. The L33 ladder itself was run and its first rung did
not certify, so by its own acceptance the ladder ends there.

| case | seed | passes | wall clock | verdict |
|---|---|---|---|---|
| `molecular_scalar_gj1132_kzz1e9/HeH0.55` | the certified reference solution `kzz1e9/HeH2.13`, carried by `--reservoir He/H 0.55` | 6 | 6 h 01 m | not solved, stopped at the ceiling |
| `molecular_photochem_gj1132_kzzprofile/HeH9` | the certified wind of its atomic pair, `atomic_photochem_gj1132_kzzprofile/HeH9/k04` | 1 at each of two pseudo-time starts | 37 m 41 s | not solved, REFUSED by the element relaxation |
| `molecular_photochem_gj1132_kzzprofile/HeH2.09` | the certified wind of its atomic pair, `.../HeH2.09/k05` | 12 | 6 h 00 m | not solved, stopped at the ceiling |
| `molecular_scalar_gj1132_wellmixed/HeH2.13` | the atomic-to-molecular conversion of `.L34/atomic_wm2.13`, solved for it here | 7 | 6 h 01 m | not solved, stopped at the ceiling |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | the atomic-to-molecular conversion of its certified atomic pair | 31 | 6 h 01 m | not solved, stopped at the ceiling |
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | the same | 40 + 40 | 5 h 19 m | not solved, REFUSED on its own terms, both budgets spent |
| `.L34/rung1`, rung 1 of the L33 ladder | the certified reference solution, `--reservoir He/H 1.411074554655` | 2 + 4 | 6 h 00 m | did not certify, which ENDS the ladder |
| `.L34/atomic_wm2.13`, the atomic well-mixed wind at He/H = 2.13 | `atomic_scalar_gj1132_wellmixed/HeH1`, reservoir carried | 0 | 2 m 44 s | **CERTIFIED**, `info = 0`, `\|\|R\|\| = 2.476e-08` |

The one solution this item adds is the last row: an atomic well-mixed wind at
the reference composition, which the catalog's atomic ladder does not carry
and which the molecular well-mixed case of that composition needs as its only
admissible seed.
## 1. The route, the hosts and the budgets

Every solve used the alternation that certified the reference solution:
`Restart intent: stationary`, `Coupled carrier solve: False`,
`EXHALE_PTC_DTAU0=1.0`, `Well balanced: True`,
`Secondary_ionization: Immediate`, `OMP_NUM_THREADS=8`, through
`LHS1140b/models/run_case.sh`, which also runs the runner's own pseudo-time
continuation at `EXHALE_PTC_DTAU0=1.0e8` when the first solve does not return
`info = 0` (item L4e). Three solves ran on `lart4` and three on `lart3`, both
72-core hosts. The pass budget was 40 for a catalog case and 25 for a ladder
rung, and the wall-clock ceiling was six hours a solve; a run at the ceiling
was stopped by stopping its binary by PID after its `/proc/<pid>/cwd` was
checked to be that case's own directory, which lets the runner write its own
`REPRODUCE.md` before it exits.

The three `molecular_scalar_gj1132_wellmixed` cases carried
`Restart intent: stationary equilibrate` and did not state
`Coupled carrier solve`; both keys were set to the recipe before the run, as
L34b had already done for the cases it re-ran. `equilibrate` puts the loaded
composition on its own fixed point before the measurement, which is a change
of the state and not the route the reference solution was reached by
(`input_read.f90` lines 76 to 81).

## 2. What a molecular case may be seeded from, measured on this binary

Two of the three routes the brief offers were tried first and are refused by
the restart contract, not by a tolerance. Both refusals are properties of the
equation set, both were MEASURED on this binary, and both decide the route the
rest of this item took.

### 2.1 A scalar molecular state cannot seed a profile case

`molecular_photochem_gj1132_kzzprofile/HeH9` was mapped from the certified
`molecular_scalar_gj1132_kzz1e9/HeH9.7` onto the profile grid, with the target
grid header built by `src/utils/profile_match_level.py` (R0 = 1.1593574793765156E+09 cm,
the radius the profile carries at its matching level times R_J) and the
reservoir carried by `--reservoir He/H 9.0102683952787981`, which is the route
`run_case.sh` takes for a case with a `Lower atmosphere profile:`. `load_IC`
refused it, MEASURED:

```
 (load_IC) ERROR: metadata field "reservoir": this run carries a reservoir for C/H
   ( 2.778025321E-04) and the restart files do not.
```

A `scalar` case is metal-free by definition (`MODELS.md` section 2), so its
state file states one ratio, He/H, while a `photochem` case states He/H, C/H,
N/H and O/H at the matching level. `profile_match_level.py --reservoir-options`
can carry only the ratios the donor state states, and no factor can create a
carbon reservoir out of a state that has none. The two are different
compositions and the refusal is correct.

### 2.2 A diffused molecular state cannot seed a well-mixed run

`molecular_scalar_gj1132_wellmixed/HeH2.13` was mapped from the certified
reference `molecular_scalar_gj1132_kzz1e9/HeH2.13` at the same `q_H2_base`,
with `Restart option change: he_diff` so that the option difference itself is
admitted. The option was admitted and the state was still refused, MEASURED:

```
 (load_IC) the restart CHANGES the physics options it was allowed to change:
  he_diff=T -> he_diff=F
 (load_IC) ERROR: the restart file does not carry the He/H the input asks for, and the state carries
   HeH+, whose nucleus of each element cannot be rescaled by one factor.
     He/H in the file   6.042513123E-01   (cell 500, r = 29.031 R_p)
     He/H in the input  2.130000000E+00
     relative departure 7.1631E-01, threshold 1.0000E-06
```

`load_IC.f90` lines 650 to 663 carry the rule: without `He_diffusion` the input
is the authority on the composition and the whole column is rescaled onto it
cell by cell, and that rescaling is refused outright when the state carries
HeH+, because HeH+ holds one nucleus of each element and no single factor per
cell sets both counts. A state produced by the binary element diffusion has a
He/H that varies with radius by construction: on the certified reference it
runs from the reservoir 2.13 at the base to 0.604 at the top cell, a departure
of 0.716.

**So the element-treatment step of the L33 ladder cannot be taken on a
molecular state at all**, whatever the order of the rungs. The ladder's
composition leg stays inside one element treatment; the crossing between
treatments has to be made where there is no HeH+ to conserve, that is on the
atomic state, before the atomic-to-molecular conversion.

## 3. What this does to the L33 ladder

The L33 continuation (memo section 8.2) and the extension L34b section 5 wrote
for the three well-mixed cases both assume that a certified molecular state can
be carried to a neighbouring composition, and, for the well-mixed leg, across
the element treatment. Measured on this binary:

| step | admitted? | what decides it |
|---|---|---|
| a composition step inside `kzz1e9` (element diffusion on) | **yes** | `load_IC.f90` lines 672 to 677: with `He_diffusion` the cell-by-cell element split IS the state, only the base rows are set to the reservoir, and the HeH+ clause is never reached |
| the element-treatment step, `kzz1e9` to `wellmixed` at fixed composition | **no** | the state's He/H varies with radius (2.13 at the base to 0.604 at cell 500 on the certified reference) and a well-mixed run rescales the whole column onto its own He/H, which HeH+ forbids |
| a composition step inside `wellmixed` | **no at the catalog's step size** | the run's He/H is the authority in every cell and the state must already carry it to 1e-6. `map_state_to_grid.py --reservoir He/H` can carry a well-mixed molecular state part of the way: MEASURED on the 2026-09-16 `wellmixed/HeH0.55` state carried to He/H = 1.0, the mapper accepts it (the base rows come out at the new ratio) and the mapped column then holds He/H between 0.999997254 and 1.000000000, a worst departure of 2.746e-06 at row 409 against the loader's threshold of 1.0e-06. The departure is the HeH+ the one factor cannot set, and it falls with the step |

So the composition leg of the ladder runs inside `kzz1e9`, where it is
admitted, and the crossing between element treatments has to be made where
there is no HeH+ to conserve: on the ATOMIC state, before the
atomic-to-molecular conversion. That is the route this item took for the
well-mixed family, and it is the route the catalog itself already used for
`wellmixed/HeH0.55` and `HeH0.083`, whose atomic pairs exist.
`wellmixed/HeH2.13` had no atomic pair, which is the whole of the reason
`MODELS.md` recorded it as never started; one was solved for it here.

---

## 4. The L33 ladder, run and ended at its first rung

Rung 1 is `q_H2_base = 0.261627353`, He/H = 1.411074554655, four steps of
0.0715 short of the catalog case `kzz1e9/HeH0.55` and 0.179 dex in composition
from the certified reference it was seeded from. Its directory is
`LHS1140b/models/.L34/rung1`; `input.inp` and `base.inp` are the reference's
with the three keys the L33 memo names changed and nothing else, so the
`Planet name:` line still carries the reference's name, which is a label and
enters no equation. Rungs 2 to 4 were created the same way and were never run,
because a rung that does not certify ENDS the ladder
(`docs/lhs1140b_stationary_L33_20260917.md` section 8.3).

MEASURED, the first solve at `EXHALE_PTC_DTAU0=1.0`:

| pass | hydro info | carrier balance H2 | cell | r [R_p] | hydrodynamic mass row | wall [s] | how the relaxation ended |
|---|---|---|---|---|---|---|---|
| 1 | 2 | 8.346e-01 | 252 | 1.3705 | 4.402e-02 | 2176.1 | the movement bound |
| 2 | 2 | 2.963e-01 | 225 | 1.2309 | 7.208e-02 | 1234.1 | **no admissible advance**, entry composition restored |

and the continuation at `EXHALE_PTC_DTAU0=1.0e8`, stopped at the ceiling:

| pass | hydro info | carrier balance H2 | cell | r [R_p] | hydrodynamic mass row | wall [s] |
|---|---|---|---|---|---|---|
| 1 | 2 | 5.030e-01 | 225 | 1.2309 | 1.13e-01 | 2405.5 |
| 2 | 0 | 3.670e-01 | 230 | 1.2521 | 4.30e-09 | 1124.4 |
| 3 | 0 | 3.070e-01 | 232 | 1.2611 | 4.62e-09 | 977.6 |
| 4 | 2 | 2.690e-01 | 234 | 1.2704 | 1.82e-01 | 2002.8 |

**Two things the rung measures, and neither is the pass budget.** First, the
outer loop of the first solve ended at pass 2 of a budget of 25, because the
element composition relaxation took 28 steps and found no admissible 29th
(outcome 3 at a step length 8.127e-02 of the composition time scale,
mass-closure departure 1.332e-15) and restored its entry composition, which
makes every later pass a copy of that one. Second, the carrier row does fall,
by 0.65 and then by 0.19 a pass, while the HYDRODYNAMIC rows do not: the mass
row is 4.4e-02 to 1.8e-01 at three of the six passes against a tolerance of
order 1e-12, so most of these passes were not handed a wind. A composition
step of 0.179 dex is already past what the alternation carries on this family,
and the step L34b's extension proposed for the well-mixed leg is inadmissible
for a different reason (section 2.2), so the ladder as designed has no
remaining rung to try.

---

## 5. The six cases, measured

Each case's own `not_solved.md` carries its full pass table; this section is
the comparison the plan asks for, with the archived value beside the new one.
A row reads "--" where the quantity does not exist: a case that did not
certify runs no post-processing pass and no transit synthesis, so it has no
mass-loss rate and no equivalent width, and a run stopped inside a solve
writes no state at all.

| case | log10 Mdot archived / new | red EW [%A] archived / new | base T [K] archived / new | base `2 n(H2)/n_H` archived / new | H2 front, `x2 = 1e-2` [R_p] archived / new |
|---|---|---|---|---|---|
| `kzz1e9/HeH0.55` | 7.640 / -- | 0.1591 / -- | 948.13 / -- | 0.666201 / -- | 13.3555 (cell 453) / -- |
| `photochem/HeH9` | -- / -- | -- / -- | 389.34 / 490.69 | 0.978033 / 0.975686 | 1.04379 (cell 131) / 1.04379 (cell 131) |
| `photochem/HeH2.09` | -- / -- | -- / -- | 237.92 / -- | 0.993609 / -- | 1.06374 (cell 152) / -- |
| `wellmixed/HeH2.13` | -- / -- | -- / -- | 808.33 / -- | 0.355296 / -- | 1.16256 (cell 205) / -- |
| `wellmixed/HeH0.55` | -- / -- | -- / -- | 605.40 / -- | 0.933788 / -- | 1.49874 (cell 269) / -- |
| `wellmixed/HeH0.083` | -- / -- | -- / -- | 497.36 / 606.64 | 0.974332 / 0.972427 | beyond the grid / 19.7739 (cell 477) |

The archived column is the state in each case's `output_pre_L34/`, which for
the four well-mixed and photochemical cases is an uncertified 2026-09-16 state
and not a solution; only `kzz1e9/HeH0.55` carries an archived state its own
binary certified, and its mass-loss rate and equivalent width are the
2026-09-16 post-processing pass's, kept in `output_pre_L34/case_files/`.
**Every base quantity in this table stands inside the +42 / -21 per cent
bracket of the helium third body** stated at the head of this memo, which is
larger than every difference in it.

The two cases whose runs ended on their own terms wrote a state, and those two
rows are the ones worth reading:

- `photochem/HeH9` moved its base temperature from 389.34 K to 490.69 K and
  its base chemical heat from 1.27447e-07 to 9.95519e-08 erg cm^-3 s^-1, with
  the H2 front unmoved at cell 131. It is a one-pass state and not a solution.
- `wellmixed/HeH0.083` moved its base temperature from 497.36 K to 606.64 K
  and its base heat from 3.30202e-07 to 2.00686e-07, and brought the
  `x2 = 1e-2` contour inside the grid for the first time, to 19.77 R_p. Its
  wind is converged (section 5.3 of its `not_solved.md`) and its carrier
  balance is 9.264e-02 against 1.0e-05.

### 5.1 What refuses, case by case

| case | hydrodynamic rows at the end | carrier balance H2, gated | how the relaxation ended |
|---|---|---|---|
| `kzz1e9/HeH0.55` | refusing at all six passes, mass 1.9e-02 to 1.6e-01 | 1.19e-01 at cell 254, after a minimum of 7.09e-02 | the movement bound |
| `photochem/HeH9` | mass 2.848e-01 at cell 2, energy 3.167e-01 at cell 43 | 1.000e+00 at cell 432 | no admissible advance, entry composition restored |
| `photochem/HeH2.09` | refusing at ten of twelve passes | a band, 1.85e-01 to 4.42e-01, ending 3.07e-01 | the movement bound |
| `wellmixed/HeH2.13` | **within** from pass 3 | 4.29e-01 at cell 237, flat to two digits from pass 5 | the movement bound |
| `wellmixed/HeH0.55` | **within** from pass 5 | 6.48e-02 at cell 267, after a minimum of 3.76e-02 at pass 11 | the movement bound |
| `wellmixed/HeH0.083` | **within**, mass 2.005e-08 of 2.8e-08, momentum 5.509e-13, energy 3.110e-08 of 1.0e-06 | 9.264e-02 at cell 1 | the movement bound |

**The pattern.** In the three well-mixed cases the stationary wind is found
and held and the H2 carrier balance alone refuses, with every relaxation of
every pass ending on the composition movement bound and none on its own
residual. In the two photochemical cases and in the diffused `HeH0.55` the
wind is not found either. So the well-mixed group isolates the carrier
transport as the single obstruction, which no earlier attempt had done: it is
not a seed question there, and it is not the hydrodynamics.

---

## 6. The catalog bookkeeping

- `python3 models/status.py --write` from `LHS1140b/` regenerated `MODELS.md`
  sections 7 and 8. Two changes were made to `status.py` itself and are
  described in section 7 below.
- `MODELS.md` section 3 was corrected where it had gone stale: the well-mixed
  molecular row now states the seeding rule the restart contract imposes and
  what the three cases did on it, and the photochemical molecular row states
  why a scalar molecular state cannot seed that group.
- Every `REPRODUCE.md` of this item was written by its own run, which is the
  only mode `models/write_reproduce.py` has; each names the binary md5
  `3146d11b4090306dcea75bb9718edd22`.
- `python3 LHS1140b/make_memo_figures.py` regenerated the 25
  `docs/figures/lhs1140b_*.pdf`.
- `latexmk -pdf docs/lhs1140b_exhale_vs_pwinds.tex` built cleanly:
  **80 pages**, the same page count as before this item.
- Before any case was touched its `output/` was preserved: `output_pre_L34/`
  where L34b had not already made one, and `output_L34b_stopped/` beside the
  three cases whose `output/` held L34b's stopped state. Nothing under
  `output_pre_L34/` was overwritten.

---

## 7. The two changes to `models/status.py`

**(a) `running` is no longer inferred from the absence of a state file.** The
old rule read a case with no `output/Hydro_ioniz.txt` and no solver verdict as
`running`, which is how `molecular_scalar_gj1132_wellmixed/HeH2.13` printed as
running for two days without a process behind it. A case that was refused or
stopped leaves `not_solved.md` beside it, and that file is the verdict: such a
case now prints `not solved`. A case with neither a state nor that file prints
`not started`. A state file with no verdict beside it is still the
flux-criterion stop and still prints `marching-stop`.

**(b) the claim column is renamed and says when its verdict has been
overtaken.** The column `claim on 59bfdb3f` is now `archived claim on
59bfdb3f`, because what it carries is the answer the L29 evaluate route gave
to the claim an ARCHIVED state file made, and printing it beside a current
`info = 0 certified` read as one statement about one state. A row now reads
`superseded` where the state in the case directory has been written since the
L29 table was made. The comparison is stated in the table's own note: the
`run=` stamp in the provenance header of the very file the table read, which
its "claim file" column names, against the modification time of
`models/CLAIMS_59bfdb3fc4d0.md` (2026-09-17 23:27:50 KST), and that file's
modification time where the header carries no stamp. MEASURED after the
change, 14 of the catalog's 95 rows read `superseded`: the seven low-XUV
atomic cases the L34 steps 1 and 2 worker solved, which the L29 table had as
`no_state`, and seven molecular cases L34b and this item re-solved or
re-seeded. Two molecular rows keep their archived verdict, `no_claim`, because
their runs were stopped inside a solve and wrote no new state file. Every other row's archived verdict is still a true statement about
the file on disk, because the L34 atomic re-evaluation rewrote
`Hydro_ioniz.txt` and not the `_IC` copy the table read.

---

## 8. Noticed beside the measurement, not acted on

1. **The pseudo-time continuation can undo a solve's progress and overwrite
   its best state.** MEASURED on `wellmixed/HeH0.083`: the first solve's forty
   passes brought the gated carrier row to 4.13e-02 and the continuation's
   forty took it back up to 7.37e-02, and the state on disk is the
   continuation's. The runner's rule (`run_case.sh`, item L4e) reloads and
   re-solves whenever the first solve does not return `info = 0`, which was
   written for a hydrodynamic solve stuck at a flat merit and is applied here
   to a case whose hydrodynamics were never the problem. A rule that keeps the
   better of the two states, or that does not continue when the hydrodynamic
   rows are already inside their tolerances, is a change to that script and was
   not made.
2. **The composition movement bound never recovers within a run.** It halves
   on a pass that does not progress and `trust_pass` only ever decreases
   (MEASURED in the source, `src/EXHALE_main.f90` lines 7086 to 7087, the only assignment that changes it after its initialization at line 6756; as L33 section 9 also
   records). MEASURED here on `wellmixed/HeH0.083`: 1.0e-02 to 5.0e-03 at pass
   14 of the first solve, held for its remaining 27 passes; the continuation,
   being a fresh run, starts at 1.0e-02 again and the row rises there, so
   restoring the bound is not by itself the remedy.
3. **`kzz1e9/HeH0.55` and `wellmixed/HeH0.55` now hold a state pair whose two
   halves are from different runs.** A solve writes `Hydro_ioniz.txt` only when
   it ends, so a run stopped inside one leaves the previous run's
   `Hydro_ioniz.txt` beside its own `*_IC.txt` seed. The two files are states
   of two different runs and nothing in the directory says so except the
   headers and each case's `not_solved.md`.
4. **The archived `run_dtau0_first.log` of a case is not renamed when a later
   run does not reach the continuation**, so a stopped attempt can leave a
   first-solve log belonging to a run months older beside its own `run.log`.
   MEASURED on `wellmixed/HeH0.55` and `kzz1e9/HeH0.55`, whose
   `run_dtau0_first.log` files are dated 2026-09-16 while their `run.log` is
   this item's.

---

## 9. Where the record is

- Each case's own `not_solved.md` and `REPRODUCE.md`, in the case directory.
- The ladder rung: `LHS1140b/models/.L34/rung1/` with its two logs and its
  `REPRODUCE.md`; `rung2`, `rung3` and `rung4` hold their `input.inp` and
  `base.inp` and were not run.
- The atomic well-mixed wind solved for the well-mixed molecular seed:
  `LHS1140b/models/.L34/atomic_wm2.13/`, `info = 0`, certified.
- The refused element-treatment rung of section 2.2:
  `LHS1140b/models/.L34/rung_wm2.13/`, which holds its `input.inp` with the
  `Restart option change: he_diff` line, its `base.inp` and the `run.log`
  carrying the refusal.
- The run logs of the runs driven from another host:
  `LHS1140b/models/.L34/logs/`.
- No source file was changed except `models/status.py`, whose two changes are
  section 7. No golden was touched and nothing under `backup/regression/` was
  written.
