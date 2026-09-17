# L14: the five low-XUV LHS 1140 b cases that did not certify

Item L14 of `docs/PLAN_20260913_lhs_stationary.md`: the HeH2.13 cases at 0.10,
0.15 and 0.20 of the GJ 1132 spectrum, and the HeH9.7 and HeH2.13 pair at
0.01, which the campaign of 2026-09-14 11:24 left uncertified
(`LHS1140b/models/.stopped/atomic_scalar_gj1132x0.*_lowxuv_202609141124/`).
Every number below is MEASURED on this tree unless it is marked READ.

## 1. Verdict

**The cell-1 hydrodynamic energy row is not what refused these states.** The
line the brief quotes,

```
(JFNK) it 10: the acceptance gate is met and the certified hydrodynamic energy
row is  4.858E-06, above its own  1.0E-06 at cell 1; the iteration goes on
within the remaining budget
```

is a note printed INSIDE a Newton solve that says the loop-top stop was not
taken yet; five iterations later the same solve printed "the acceptance gate
is met and every one of the 3 certified row(s) of the system this solve
carries is within its own tolerance" and returned `info = 0`. Every outer
pass of that run ended with the three hydrodynamic rows inside their own
tolerances (at the last pass: mass 2.30E-08, momentum 3.39E-14, energy
8.49E-07 of 1.0E-06). **What refuses the state is the gated elemental
transport He/H partition row**, 2.12E-02 against 1.0E-05 at cell 332
(r = 2.50 R_p), and the loop ended on its own progress control:

```
(EXHALE_main) outer pass 4: REFUSED -- the joint distance of the state (its
largest entry over its own tolerance) has not fallen in 3 consecutive passes
```

**The obstruction is that progress control, and the response it makes.** The
alternation is "solve the wind at fixed composition, relax the composition at
fixed wind". Its own error is the undamped distance the element relaxation
returns, `max |X_relaxed - X| / X_base`, which is contracted by a factor
`1 - omega` per pass. The control reads a different quantity, the distance of
the state from certification, and where THAT fails to fall it halves `omega`
-- taking the contraction factor from 0.500 to 0.875 -- and after three such
passes ends the loop. On these three cases the two quantities disagree from
the first pass: the composition distance falls every pass while the gated row
rises, so the control slowed the alternation down and then declared it
stagnant.

The certified case at 0.30 of the same spectrum is the control that settles
it: its composition entered at the SAME distance (4.89E-01 against 6.62E-01
at 0.10), its gated row also failed to fall over its first pass, and it
certified at outer pass 13 with `omega` at 0.500 throughout, the distance
halving every pass.

**The 0.01 pair is a different failure and is not this one**: there the
hydrodynamic solve itself does not converge from any seed available in the
tree, and the item treats it by a continuation in the XUV normalization.

## 2. What was measured

### 2.1 Which log is which

`run_case.sh` runs the stationary solve at `EXHALE_PTC_DTAU0=1.0`; if it does
not return `info = 0` it reloads the written state and solves again at
`DTAU0_CONTINUATION` (1e8, item L4e), moving the first log to
`run_dtau0_first.log` and writing the second to `run.log`. So
`run_dtau0_first.log` is the dtau0 = 1 solve and `run.log` the continuation,
and the certified cases went through the same two stages.

### 2.2 The continuation at dtau0 = 1e8, where the state is decided

`atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13`, seeded from the certified 0.25
state of the same composition:

| pass | hydro `info` | gated elemental He/H row | mass | momentum | energy | `omega` | element distance after the update |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 1.28E-02 at cell 348 | 1.72E-08 | 1.59E-14 | 4.67E-07 | 0.500 | 6.62E-01 |
| 2 | 0 | 1.89E-02 at cell 333 | 1.40E-08 | 1.20E-14 | 5.43E-07 | 0.250 | 4.24E-01 |
| 3 | 0 | 2.06E-02 at cell 332 | 3.19E-08 | 1.17E-13 | 8.36E-07 | 0.125 | 3.39E-01 |
| 4 | 0 | 2.12E-02 at cell 332 | 2.30E-08 | 3.39E-14 | 8.49E-07 | 0.125 | REFUSED, no update |

The element relaxation reported "ended on the fixed point of the element
operator" at every one of those passes. The same shape at 0.15 (row 6.90E-03,
1.27E-02, 1.45E-02, 1.51E-02; distance 5.36E-01, 3.76E-01, 3.08E-01) and at
0.20 (row 2.86E-03, 6.10E-03, 7.36E-03, 7.84E-03; distance 3.21E-01,
2.53E-01, 2.16E-01).

### 2.3 The certified cases at the same stage

| case | passes to ACCEPTED | gated row, first pass -> last | element distance, first -> last | ratio per pass | `omega` |
|---|---|---|---|---|---|
| 0.33 / HeH2.13 | 11 | 4.46E-03 -> 8.52E-06 | 2.32E-01 -> 4.39E-04 | 1.98 to 2.02 | 0.500 throughout |
| 0.30 / HeH2.13 | 13 | 3.56E-03 -> 6.22E-06 | 4.89E-01 -> 2.90E-04 | 1.97 to 2.00 | 0.500 throughout |
| 0.25 / HeH2.13 | 6 | 1.63E-04 -> 6.10E-06 | 2.07E-03 -> 1.35E-04 | 1.95 to 2.01 | 0.500 throughout |
| 0.33 / HeH9.7 | 6 | 1.77E-04 -> 6.90E-06 | 1.98E-02 -> 1.22E-03 | 2.00 to 2.01 | 0.500 throughout |
| 0.30 / HeH9.7 | 9 | 1.08E-03 -> 5.15E-06 | 1.07E-01 -> 8.29E-04 | 2.00 to 2.01 | 0.500 throughout |
| 0.25 / HeH9.7 | 9 | 2.03E-03 -> 9.13E-06 | 1.75E-01 -> 1.36E-03 | 2.00 to 2.00 | 0.500 throughout |

The element distance halves at every pass of every certified case, to the
third digit: that is `1 - omega` at `omega = 0.500` with the wind's response
to the composition below the third digit. On the three refused cases the same
ratio at `omega` 0.500 then 0.250 then 0.125 is 1.56, 1.25 (0.10), 1.43, 1.22
(0.15) and 1.27, 1.17 (0.20), against the 2.00 and 1.33 the under-relaxation
alone gives at those two values of `omega`: the wind's response to the
composition leaves 64 to 78 percent of the contraction at `omega` = 0.500 and
88 to 94 percent at 0.250, so the alternation contracts on these states too.
**What removes most of the contraction is the halving of `omega` itself**,
which takes the factor the under-relaxation sets from 2.00 to 1.14 per
pass.

Helium composition and not spectrum shape is what separates the refusals: at
He/H = 9.7 the 0.10, 0.15 and 0.20 cases are all certified on this tree, and
only the He/H = 2.13 column of the same three normalizations refused.

### 2.4 The first solve at dtau0 = 1

On the same three cases the dtau0 = 1 solve never reaches a stationary wind:
the hydrodynamic energy row stands at 6.27E-01, 6.31E-01, 6.32E-01, 6.33E-01
over its four passes (0.10; 4.24E-01 to 4.29E-01 at 0.15, 2.14E-01 to
2.18E-01 at 0.20) against a tolerance of 1.0E-06, with `hydro info = 2` every
pass, while the element distance falls 1.73E-01, 1.13E-01, 8.83E-02. Its
ending is the same progress control, and it is the RIGHT ending there: a
composition closing on the fixed point of a wind that is not a solution says
nothing about the pair, and the passes would only delay the continuation.

### 2.5 Why the low-XUV cases enter so much further out

At the cell that binds the elemental row, r = 2.50 R_p, the bulk speed falls
with the normalization while the temperature falls with it:

| XUV | v [cm/s] at 2.50 R_p | T [K] | He/H at 1.50 R_p | He/H at 2.50 R_p |
|---|---|---|---|---|
| 0.33 | 1.825E+03 | 2592 | 7.12E-01 | 1.71E-01 |
| 0.30 | 1.624E+03 | 2475 | 6.82E-01 | 1.38E-01 |
| 0.25 | 1.259E+03 | 2180 | 5.87E-01 | 6.47E-02 |
| 0.20 | 9.104E+02 | 1967 | 3.14E-01 | 2.60E-02 |
| 0.15 | 6.244E+02 | 1881 | 8.33E-02 | 2.20E-02 |
| 0.10 | 4.031E+02 | 1862 | 4.73E-02 | 2.09E-02 |

with He/H = 2.13 at the base in every case. A weaker wind separates helium
from hydrogen much more strongly and much lower down -- at 0.10 the ratio has
fallen by a factor 45 by 1.5 R_p, against a factor 3.6 at 0.25 -- so a
composition carried over from the 0.25 solution has much further to travel,
and the wind it is relaxed against changes more while it travels. That is why
the gated row rises over the first passes at 0.10 to 0.20 and barely falls at
0.30, and why it falls cleanly at 0.25.

## 3. What changed

### 3.1 `src/EXHALE_main.f90`, the progress control of `steady_wind_with_element_diffusion`

A pass counts as progress when the joint distance from certification fell, OR
when the state's distance from certification is carried by an entry other
than a hydrodynamic row and the composition's own distance to the fixed point
of its operator fell. Only when neither holds are the step lengths shortened
(`omega`, the carrier movement bound) and the consecutive-pass counter
advanced; the ending itself, `outer_no_fall_max = 3` consecutive such passes,
is unchanged.

Both parts of the second condition are needed and both are measured: the
composition distance is the quantity the alternation actually contracts
(sections 2.2 and 2.3), and it only means that where the wind half of the
alternation is not itself what refuses the state (section 2.4). The test is
which entry carries the distance and not whether the hydrodynamic rows are
inside their tolerances; section 4.6 measures why. The distance is read from
the two most recent updates, so the first two passes of a solve carry no
measurement of the contraction and count as progress.

### 3.2 `LHS1140b/models/run_case.sh`

- `EXHALE_OUTER_PASSES` is exported with a default of 40 and named in
  `REPRODUCE.md`. The binary's own default is 20, which the certified 0.30
  case already spent 13 of from a nearer entry, and the 0.10 case of this
  item was ACCEPTED at pass 30 exactly.
- `EXHALE_PTC_DTAU0` is now a DEFAULT of 1.0 rather than a fixed assignment,
  so a caller may state the pseudo-time start of the FIRST solve. Nothing
  about a run that does not state it changes: the variable is set to 1.0
  exactly as before, the continuation at `DTAU0_CONTINUATION` is reached by
  the same condition, and the value is written into `REPRODUCE.md` either
  way. What it makes possible is section 4.6's recipe, which the ladder of
  section 4.8 needs: where the dtau0 = 1 solve cannot leave the ramp floor,
  `EXHALE_PTC_DTAU0=1.0e8` takes the mapped seed straight to the solve that
  works, and the continuation then does not run because the first solve
  returns `info = 0`.
- The walk down `pick_seed.py`'s ranked list is skipped when the refused
  solve left every hydrodynamic row inside its tolerance, read from the
  certification block of the log. Such a state is a stationary wind that the
  composition refuses, and a seed one further XUV step away can only make the
  wind worse: MEASURED on these three cases, the
  attempt from the next seed spent its passes at `lam` 1e-6 and `dtau` 1e-4
  without moving the residual (plan item L4h).

### 3.3 `LHS1140b/sed/`

`lhs1140_sed_gj1132_at_b_xuv0p05.txt` and `..._xuv0p03.txt`, the whole flux
column of `lhs1140_sed_gj1132_at_b.txt` times 0.05 and 0.03, written the same
way the existing grid was: regenerating `..._xuv0p10.txt` by that rule
reproduces the file on the tree in all 29992 rows, no mismatch. They exist
for the continuation of the 0.01 pair and carry no catalog case.

## 4. Verification

### 4.1 The control build

`EXHALE.x` as it stood at md5 `b9419ba7cd8e519ca86896b29c382762`, built
2026-09-14 10:12, is the control: every Fortran source of this working copy
carried an mtime of 10:04 or earlier apart from `src/EXHALE_main.f90`, which
this item changed at 11:52, so that binary is the same tree minus this
change. It is what the refused runs of 11:24 and the `lower_profile` control
of section 4.7 were made with.

The measured build is `EXHALE_L14.x`, md5
`a11038c245050d4f11852af829e13fc7`, built from `build_l14/` with
`make OBJDIR=build_l14 MODDIR=build_l14 EXE=EXHALE_L14.x` so that `EXHALE.x`
was left where the campaign running beside this item expected it. (An
intermediate `dd7ee3188ead1cf54d2f2c7f99246309` carried the first form of the
condition, section 4.6, and the three cases of section 4.5 were certified
with it; they take no `omega` halving under either form, so the delivered
binary reproduces them.) `EXHALE.x` was itself rebuilt from this tree at
2026-09-15 08:23 and now carries the same md5, so the control binary no
longer stands on disk; the control measurements above were all taken while
it did.

`src/EXHALE_main.f90` stands at +165/-29 against `43bc28c`, of which this
item's change is the progress control, the two distances it reads and the
routine header. `LHS1140b/models/run_case.sh` and `LHS1140b/sed/` are not
in the git remote.

### 4.2 The test suites the change touches

Run against `build_l14/` and the DELIVERED `EXHALE_L14.x`
(`a11038c245050d4f11852af829e13fc7`):

| suite | result |
|---|---|
| `certification` | every row PASS (`the_column_is_not_certified_on_its_account`, `the_unresolved_cell_is_where_the_verdict_binds` 4/4); `certification_contexts`: every assertion passed |
| `steady_species_rows` | every row PASS (`no_species_row_writes_no_composition` 0.0, `element_constraint_rows_are_red_without_them` 6/6, `inventory_shared_rows_are_red_without_the_constraint` 3/3); the geometry and reload rows skip for want of the logs and the binary they name, on any build |
| `state_mapper` | PASSED, 0 failures |
| `element_operator` | every row passed |

`certification` is the suite that owns `within_tol`, which the new condition
reads; `element_operator` owns the element relaxation whose returned distance
it reads.

### 4.3 The dtau0 = 1 stage is the control build, pass for pass

On all three cases the first solve reproduces the refused run of 11:24 to
every printed digit -- 0.10: 5.75E-03 at cell 350, 3.07E-03 at 353, 2.23E-03
at 354; 0.15: 3.39E-03 at 346, 1.74E-03 at 348, 1.28E-03 at 349; 0.20:
1.52E-03 at 344, 7.67E-04 at 345, 5.66E-04 at 346, 4.91E-04 at 346 -- and
ends REFUSED at pass 4 with `omega` walked 0.500, 0.250, 0.125 exactly as
before. That is the new condition doing what it is written to do: the
hydrodynamic energy row of that solve stands at 2.1E-01 to 6.3E-01 against
1.0E-06, so the composition's progress is not read and the ending is
unchanged.

### 4.4 `atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13`, the first case through

CERTIFIED at outer pass 19 of the continuation, where the control refused at
pass 4:

```
[atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13] DONE info=0 ||R||=1.228E-07
    certified Mdot=7.14 EW=0.0026 %A depth=0.0089 %
```

| pass | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| gated elemental row | 2.86E-03 | 6.10E-03 | 8.26E-03 | 9.06E-03 | 8.83E-03 | 7.81E-03 | 6.25E-03 | 4.52E-03 | 2.98E-03 | 1.81E-03 |
| element distance | 3.21E-01 | 2.53E-01 | 1.78E-01 | 1.17E-01 | 7.42E-02 | 4.58E-02 | 2.75E-02 | 1.62E-02 | 9.33E-03 | 5.25E-03 |

| pass | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 |
|---|---|---|---|---|---|---|---|---|---|
| gated elemental row | 1.05E-03 | 5.91E-04 | 3.28E-04 | 1.81E-04 | 9.99E-05 | 5.44E-05 | 2.93E-05 | 1.54E-05 | 7.99E-06 |
| element distance | 2.88E-03 | 1.54E-03 | 7.99E-04 | 4.01E-04 | 1.94E-04 | 8.98E-05 | 5.91E-05 | 4.04E-05 | -- |

`omega` was never halved (no such line in the log). The first two passes are
identical to the control's; the row then rises to a peak at pass 4, exactly
where the control ended the loop, TURNS OVER, and falls by 3.1 decades over
the next fifteen passes while the element distance falls monotonically by 3.9
decades. That is the claim of section 1 measured end to end: the rise was the
composition travelling, not a stalled alternation.

### 4.5 The three cases of the item, all CERTIFIED

| case | outer pass ACCEPTED | `||R||` | log Mdot | EW [%A] | red depth [%] | `omega` halvings |
|---|---|---|---|---|---|---|
| `atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13` | 19 | 1.228E-07 | 7.14 | 0.0026 | 0.0089 | 0 |
| `atomic_scalar_gj1132x0.15_kzz1e9/HeH2.13` | 24 | 7.112E-07 | 7.01 | 0.0004 | 0.0012 | 0 |
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13` | 30 | 4.842E-07 | 6.83 | 0.0001 | 0.0003 | 0 |

Each wrote its own `REPRODUCE.md`, its `output/` pair carries `certified=T`,
and each went through the same two stages as the certified cases above: a
dtau0 = 1 solve refused at pass 4, then the continuation. The passes needed
grow as the seed's XUV step grows in decades, and the 0.10 case was ACCEPTED
at pass 30 of a 30-pass budget: `run_case.sh` now asks for 40 for that
reason.

### 4.6 The correction the first ladder rung forced

The 0.01 pair cannot be seeded from anything on the tree: the hydrodynamic
solve does not converge across that step (section 2). The ladder is
`.L14/x005_*` and `.L14/x003_*`, at 0.05 and 0.03 of the spectrum, and the
0.01 cases are then seeded from the nearest rung. Those are working
directories and carry no catalog case.

The first rung MEASURED a defect in the first form of the new condition.
`.L14/x005_HeH9.7` spent its whole 30-pass budget without certifying, at a
gated elemental row pinned at 2.64E-02 and an element distance creeping at a
ratio of 1.10 per pass. The log says why: `omega` was halved at pass 2 and
again at pass 10, and the hydrodynamic energy row of those two passes reads
1.03E-06 and 1.32E-06 against a tolerance of 1.0E-06, while every other pass
of the run is inside it. An "every hydrodynamic row inside its tolerance"
test therefore called a wind solved to three percent of its tolerance an
unsolved wind, and halved the step twice for it.

**The condition is now which entry carries the state's distance from
certification**, `hydro_worst_dist < prog_worst`: the composition's progress
is read when the largest distance belongs to something other than a
hydrodynamic row. It separates the two cases by five decades rather than by
three percent -- at dtau0 = 1 the energy row stands at a distance of 6.3E+05
against the elemental row's 5.8E+02, and in the continuation at 1.03 against
1.9E+03 -- and it is the statement the rule wanted in the first place: the
wind is not what refuses this state. The three cases of section 4.5 are
unaffected, having taken no halving at all under either form.

### 4.7 `mol_diffusion` and `lower_profile`, the two regression cases that could reach this code

| case | `Hydro_ioniz.txt` | `Ion_species.txt` | `Hydro_ioniz_adv.txt` | `Ion_species_adv.txt` |
|---|---|---|---|---|
| `mol_diffusion` | PASS, max rel 1.272e-07 | PASS, 1.814e-06 | PASS, 2.537e-05 | PASS, 3.083e-06 |
| `lower_profile` | FAIL, 1.269e+00 | FAIL, 9.942e-01 | FAIL, 1.269e+00 | FAIL, 1.000e+00 |

**The `lower_profile` failure is not this item's**: the same case run with the
CONTROL binary `EXHALE.x` gives the identical numbers, `du` 3.1666E+00 and
`dtu` 2.3256E-02 at step 12000 and `max rel 1.269e+00 at row 474 col 3
(518500.78714196314 vs -139260.56102646154)` in both. It is a golden that the
tree has moved past, and it is reported here and left alone. Neither case
reaches the routine this item changed in any event: both end on their pinned
12000 marching steps (`du` 7.9E-01 and 3.2E+00, above the 1.0E-02 JFNK
hand-off) and neither log carries a single outer-pass line.

### 4.8 The ladder to 0.01 XUV

Three chains, one per case, each rung solved from the state the previous rung
wrote:

- `LHS1140b/models/.L14/ladder.sh <He/H> <seed> <rung ...>`, the scalar
  chains, driven by `chain_HeH2.13.sh` and `chain_HeH9.7.sh`;
- `LHS1140b/models/.L14/ladder_photochem.sh`, the lower-atmosphere-profile
  case `atomic_photochem_gj1132x*_kzzprofile/HeH9.7`, whose rungs go through
  `run_case.sh` so that the target grid header is built from the profile's own
  matching level (item L13) and each rung keeps its `REPRODUCE.md`.

Every rung is solved at `EXHALE_PTC_DTAU0=1.0e8` for the reason of section
4.6, with the 40-pass budget. The rungs are 0.05, 0.03, 0.02 and 0.015 of the
GJ 1132 spectrum, in `.L14/` only; the 0.01 case of each chain is the catalog
directory, solved by `run_case.sh` so it gets its advection-corrected
profiles and its transit spectrum. The spectra are
`LHS1140b/sed/lhs1140_sed_gj1132_at_b_xuv0p{05,03,02,015,07}.txt`, the whole
flux column scaled, verified against the existing grid (section 3.3).

The photochem case had been given a seed 20 XUV steps away
(`atomic_photochem_gj1132x0.20_kzzprofile`), on which the hydrodynamic solve
returns `info = 2` at every pass; that attempt is kept in
`.stopped/atomic_photochem_gj1132x0.01_kzzprofile_HeH9.7_farseed_202609142135/`.

| rung | outer pass ACCEPTED | `omega` halvings |
|---|---|---|
| `.L14/x005_HeH2.13` | 20 | 0 |
| `.L14/x005d8_HeH9.7` | 29 | 0 |

### 4.9 Every rung of every ladder

| rung | outcome | passes taken | `omega` halvings |
|---|---|---|---|
| `.L14/x005_HeH2.13` | ACCEPTED at pass 20 | 20 | 0 |
| `.L14/x003_HeH2.13` | ACCEPTED at pass 19 | 19 | 0 |
| `.L14/x002_HeH2.13` | REFUSED, see below | 40 (the budget) | 0 |
| `.L14/x005d8_HeH9.7` | ACCEPTED at pass 29 | 29 | 0 |
| `.L14/x003_HeH9.7` | ACCEPTED at pass 24 | 24 | 0 |
| `.L14/x002_HeH9.7` | REFUSED, stopped at pass 13 | 13 | 0 |
| `.L14/p005_HeH9.7` (profile) | ACCEPTED at pass 29 | 29 | 0 |
| `.L14/p003_HeH9.7` (profile) | ACCEPTED at pass 25 | 25 | 0 |
| `.L14/p002_HeH9.7` (profile) | REFUSED, stopped at pass 1 | 1 | 0 |

Not one `omega` halving in any rung of any chain. Where a rung carries its
certification differs by how it was run and both are read that way by the
drivers: a rung solved by the binary directly keeps its solved state in
`output/Hydro_ioniz.txt` and its mapped seed in `output/*_IC.txt`, while a
rung solved through `run_case.sh` keeps the certified state in
`output/*_IC.txt` and lets the advection-corrected pass rewrite
`output/Hydro_ioniz.txt` with a header that makes no stationary claim.

Every chain certifies at 0.05 and at 0.03 and refuses at 0.02. The refusal is
one and the same on all three and is the subject of the next section; on
`x002_HeH2.13`, the one rung that ran its whole budget, the constant
`mass 1.62E-04` of its pass lines is not a truncation floor but the
consequence of `info = 2`, which hands back the state the solve was given, so
the hydrodynamic rows are identical from pass to pass by construction.

### 4.10 The 0.02 wall, on all three chains, and what it hands to L4h

All three ladders certify at 0.05 and at 0.03 and refuse at 0.02, in the same
shape and for the same reason. The last Newton iteration each 0.02 rung
reached:

```
x002_HeH2.13  (JFNK) it  75  ||R||= 1.472E-03  ||Fs||2= 1.19E-08  lam= 4.88E-04
              dtau= 7.93E-05  gm= 1  worst r= 1.000  worst row: mass of cell 1
x002_HeH9.7   (JFNK) it  25  ||R||= 5.273E-03  ||Fs||2= 5.40E-08  lam= 5.00E-01
              dtau= 2.56E+10  gm= 24  worst r= 1.002  worst row: energy of cell 12
p002_HeH9.7   (JFNK) it 120  ||R||= 2.810E-02  ||Fs||2= 4.36E-07  lam= 6.25E-02
              dtau= 3.91E+05  gm= 33  worst r= 1.005  worst row: energy of cell 26
```

`hydro info = 2` at every outer pass of all three, the worst row always
within 0.005 R_p of the base, `lam` and `dtau` walking down from a start of
1e8. **The composition is not what refuses these states**: on
`x002_HeH2.13`, which ran its whole 40-pass budget, the gated elemental He/H
row fell from 2.95E-02 to 4.58E-08, four decades BELOW its 1.0E-05 tolerance,
with `omega` at 0.500 and no halving. The alternation this item repaired
works; the wind half does not reach its own tolerances.

**It is an XUV level and not a step size.** The step 0.03 -> 0.02 is a factor
1.5, gentler than the 0.05 -> 0.03 (factor 1.67) that certified on every
chain, and gentler than the 0.10 -> 0.05 (factor 2) that certified on two. It
is also independent of the helium abundance (2.13 and 9.7 alike) and of the
lower-atmosphere treatment (the scalar `K_zz` and the `kzzprofile` case
alike). What the three have in common is the pseudo-time ramp collapsing on a
run of damped steps at the base, which is plan item L4h and is that item's to
repair.

The three 0.02 rungs were stopped at that point and their states are kept.
`.L14/` is left in place: its certified 0.03 states are the seeds a resumed
ladder starts from.

### 4.11 Where the five cases stand

| case | verdict line of `run_case.sh` |
|---|---|
| `atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13` | `DONE info=0 \|\|R\|\|=1.228E-07 certified Mdot=7.14 EW=0.0026 %A depth=0.0089 %` |
| `atomic_scalar_gj1132x0.15_kzz1e9/HeH2.13` | `DONE info=0 \|\|R\|\|=7.112E-07 certified Mdot=7.01 EW=0.0004 %A depth=0.0012 %` |
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13` | `DONE info=0 \|\|R\|\|=4.842E-07 certified Mdot=6.83 EW=0.0001 %A depth=0.0003 %` |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | NOT CERTIFIED, blocked on L4h; ladder reached 0.03 |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | NOT CERTIFIED, blocked on L4h; ladder reached 0.03 |
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | NOT CERTIFIED, blocked on L4h; ladder reached 0.03 |

The two lower-XUV profile cases the campaign certified beside them, for the
same spectrum grid: `atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7`
`DONE info=0 ||R||=7.440E-07 certified Mdot=7.04 EW=0.5065 %A depth=1.8713 %`
and `atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7`
`DONE info=0 ||R||=1.921E-07 certified Mdot=6.85 EW=0.2823 %A depth=1.0575 %`.

The He 10830 line is already gone at these levels on the scalar cases
(equivalent width 0.0026, 0.0004 and 0.0001 percent-angstrom at 0.20, 0.15
and 0.10 against 0.2823 on the profile case at the same 0.10), so the 0.01
cases are a completeness item of the catalog and not a measurement the
comparison waits on.

## 5. Found beside the item, reported and not acted on

- `src/tests/grid_and_gates/output_state_consistency.sh`, the row
  `outer_iteration_ending_is_the_stagnation_one`, FAILS ON THE CONTROL BUILD
  (`EXHALE.x`, md5 `b9419ba7cd8e51`): its fixture reaches the pass budget,
  `measured=pass_budget reference=no_progress`, where the comment records the
  stagnation ending measured on 2026-09-12. It fails identically on the
  measured build, so this item did not cause it and did not repair it. The
  test says what to do -- "pick a setting that reaches the ending again rather
  than loosening the second assertion" -- and that is an item of its own.
- `backup/regression/` holds fixture directories named `armA_LW`, `armA_noLW`,
  `armD_D2*`, `armHeH_*` and `arm_heh1_x2matched`, which use "arm" as a noun.
  `docs/named_case_audit.md` already carries open proposals on those
  directories awaiting a user decision, and the names are cited in about
  eighteen documents, so they are left where they are and reported here.
  **DONE 2026-09-16** (PLAN_20260916_rev3 section 10): all sixteen were
  renamed; the old-to-new mapping is `docs/named_case_audit.md` section 6.

## 6. The progress control reads a distance, not a displacement (item R2 of the review of 2026-09-15)

### 6.1 Verdict

**The carrier half of `comp_drift` was wrong, and the sign of the error is the
worst one available.** `composition_closing`, the condition this item wrote,
reads

```fortran
composition_closing = (hydro_worst_dist .lt. prog_worst) .and.  &
                      ((n_comp_updates .lt. 2) .or.             &
                       (comp_drift_last .lt. comp_drift_before))
```

and `comp_drift` was, in its carrier half, `carrier_drift` -- the change the
relaxation KEPT. That change is bounded by `trust_pass`, which the same loop
shortens whenever a pass does not help. So the loop shortens the bound, the
displacement falls because the bound fell, and the control reads "the
composition is closing" precisely when the loop has stopped letting it move.

The distance of a composition from the fixed point of its own operator is the
RESIDUAL of that operator at the composition, on a scale no movement bound can
touch, and the carrier module already forms exactly that for the certification
(`carrier_steady_residual`). That is what the carrier half now is, measured on
the state the pass hands back, in the isolated workspace the certification uses
(`save_carrier_module_state` / `restore_carrier_module_state`), so the frozen
background, the photolysis rates and the row terms are the ones the next step
would read.

The ELEMENT half is unchanged: `relax_element_composition` returns an undamped
distance already (`comp_omega` damps what is applied, not what is measured), so
it is a distance and not a displacement.

`EXHALE_CARRIER_DRIFT_IS_DISPLACEMENT=1` restores the old measure.

**The element half is no longer what this section describes.** Item L22 step 1
(`docs/lhs1140b_stationary_L22_20260916.md`) made it the elemental transport
residual of the composition the pass hands back, for the reason the sentence
above does not cover: an undamped distance is still the distance to the
endpoint of a FINITE inner relaxation, which can exit on its step budget, so it
is not the residual of the returned state.
`EXHALE_ELEMENT_DRIFT_IS_MAP_DISTANCE=1` restores the measure this section
was written about. The carrier half is as described here.

### 6.2 What it changes, measured

The fixture is the one the test suite already runs the uncoupled carrier
relaxation on: `backup/regression/carrier_elem_newton` with
`Coupled carrier solve: False` (the `outer_iteration_ending` stage of
`src/tests/grid_and_gates/output_state_consistency.sh`), single-threaded,
binary `EXHALE_L21b.x` (`6463f0fdf181bd1859d59da3226e7819`), the two measures
selected by the flag.

The two logs differ in **three lines and nothing else** -- no number in the
run changes, and the exception summary at exit is the only other difference:

```
>     -> the carrier kept no transport step last pass; the movement bound is left at 1.00E-06
```

three times, present with the DISPLACEMENT measure and absent with the
residual. Those are three passes on which the old measure declared the
composition not to be closing -- and so counted a pass without a fall and
entered the branch that decides what to do with the movement bound -- while the
carrier's own residual was in fact falling. The correction removes three false
stagnation counts.

**The verdict is the same**: both runs end at the pass budget of 8 with
`worst gated species row 7.37E-02 of 1.0E-05 at cell 290 (carrier balance H2)`,
mass 9.48E-12, momentum 7.03E-10, energy 3.49E-09, omega 0.500, trust 1.0E-06,
and the same refusing entry. The change does not rescue this fixture; it stops
the control from reading its own bound.

### 6.2b A fixture where the carrier does move, and the difference is not cosmetic

`backup/regression/carrier_model_a_newton` (`Molecular carrier transport: True`,
`Coupled carrier solve: False`, `Restart intent: stationary`, 20 outer passes,
one thread, same binary, same reloaded state). Here the relaxation keeps
transport steps, so the two measures part.

| | residual measure (default) | displacement measure (old) |
|---|---|---|
| pass 7 | the composition IS closing; the movement bound is left at 1.00E-02 | `neither the joint distance nor the composition distance fell; carrier movement bound = 5.00E-03` |
| pass 7 carrier drift kept | 1.79E-02 in 11 transport steps | 9.18E-03 in 10 transport steps |
| pass 7 worst carrier balance | 7.43E-04 at cell 373 (H2) | 7.98E-03 at cell 131 (H+) |
| pass 20 | worst gated row 5.21E-03 at cell 248 (carrier balance H+), trust 1.0E-02 | worst gated row 5.75E-04 at cell 500 (carrier balance H2), trust 2.5E-03 |
| ending | budget of 20 passes, 2 active equations refuse it | the same |
| written state | differs | differs |

**The ending is the same and the iterate is not.** The old measure halves the
bound at pass 7 because the displacement it reads fell -- and the displacement
fell because the bound had been halved before -- and it ends at a quarter of the
bound with a smaller worst carrier row than the new one. That smaller row is
not an argument for the old measure: a carrier held still by its own bound
perturbs its balance row less, which is the circularity the correction removes.
Neither run certifies, so nothing here decides between them on the answer; what
decides is that a progress control must not read a quantity the same control
shortens.

### 6.3 The L14 fixture itself is not touched by this

The three cases of section 4.5 (`atomic_scalar_gj1132x{0.20,0.15,0.10}_kzz1e9/HeH2.13`)
are ATOMIC: `thereis_mol` is false, and the block this item changed is gated on
`thereis_mol .and. carrier_transport .and. .not. carrier_in_newton`. Their
`comp_drift` is the element half alone, which is unchanged, so their progress
verdict cannot move.

MEASURED, not only argued: `x0.20` re-solved from the seed its own
`REPRODUCE.md` names
(`models_20260914_preL21/atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13/output`,
interpolated by `src/utils/map_state_to_grid.py`) at
`EXHALE_PTC_DTAU0=1.0e8`, `EXHALE_OUTER_PASSES=40`, ONE thread, twice: once
with the residual measure and once with `EXHALE_CARRIER_DRIFT_IS_DISPLACEMENT=1`.
**Every outer-pass line of the two runs is identical in every field but the
wall time** -- pass 1 `worst gated species row 2.17E-04 at cell 217
(elemental transport He/H partition), mass 3.35E-08, momentum 2.79E-14,
energy 2.46E-07, omega 0.500, trust 1.0E-02`, and so on, which is also the
line for line the `REPRODUCE.md` of the certified run records. The runs were
left going past the pass count reported here; through the passes measured
they do not part -- IDENTICAL through outer pass 20 at the time of writing -- and
they cannot, for the reason above.

(At one thread this case does not reach the acceptance its 8-thread campaign
run reached at pass 19. That is item L15, not this item: the JFNK solve is not
reproducible across thread counts, and a control/measured pair is only read
at one thread, which is what was done.)

**NOTE ADDED 2026-09-17 (item L15).** No thread dependence was found on the
catalog binary or on the 2026-09-16 tree, atomic or molecular, at 1, 8 or 16
threads: the traced trajectory of the element reload is bitwise identical
across thread counts, and so is each molecular case. The sentence above is
left as written because it states why the comparison was taken at one thread,
which remains the rule; what it asserts about the solver is not reproducible
today. `docs/lhs1140b_stationary_L15_20260916.md`.

A case with `Coupled carrier solve: True` is likewise untouched: the carrier is
then a Newton unknown and the relaxation block does not run at all.

## 7. Three indices, logged apart (item R3 of the review of 2026-09-15)

One line used to carry one cell index for three different questions. They are
now printed separately, at the end of each pass (the drift line only where the
pass kept a transport step, since a pass that kept none has no drift cell):

| what | where it comes from | example, `carrier_model_a_newton` pass 7 |
|---|---|---|
| the worst carrier drift cell | `carrier_drift_location` | `carrier drift worst at cell 333 r= 1.641 carrier H2` |
| the cell holding the particle-count bound | `bound_last_j`, `bound_last_ic`, `bound_last_dabs` | `particle-count bound last held at cell 327 r= 1.591 carrier H2 measure 1.000E-02 (0 cell(s) limited in the last solve)` |
| the worst-residual cell | `carrier_steady_residual` | `carrier steady residual: volume-weighted 3.02E-04, worst cell 7.43E-04 at 373 r= 2.100 carrier H2` |

They are three different cells in that pass (333, 327 and 373), and the row
that refuses the certification names a fourth, which is the point: one index cannot answer "where is the bound biting",
"where is the composition furthest from its own steady state" and "where does
the certification refuse".

`bound_last_j` and `bound_last_ic` were already public in the carrier module
and already imported by `EXHALE_main.f90`, so nothing had to be asked of the
module that owns the diffusive photochemistry.
