# L13. The elemental-flux closure could not reload its own iterate

2026-09-14. `LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH*/`,
driver `src/utils/element_flux_closure.py`, launched by
`LHS1140b/models/run_closure.sh`.

## 1. The failure, as measured

Every closure rung stopped at its first continuation iteration, `k01`, with
the restart contract refusing the seed (`HeH9/closure.log` and
`HeH9/k01/run.log`, 09:34:37):

```
 (load_IC) ERROR: metadata field "reservoir": the elemental ratio C/H of the
   restart files is 2.777879401E-04
   this run asks for 2.777952293E-04
   relative departure 2.6240E-05, threshold 1.0000E-06
   The state is a solution at its own composition. Restart it there, or start this
   composition cold.
```

## 2. Why

Each iteration of the closure hands the wind a new photochemical column, and
the run takes its elemental reservoirs from that column: `input_read`'s
`apply_lower_atmosphere_profile` calls `lap_element_ratio_at_match`, which
reads the `X_<El>` column at the matching level, for helium and for each of
the ten elements of `species_table`. Every one of those ratios moves from
one iteration to the next. Measured between the reservoir of the `k00`
solution and the `k01` profile of `HeH9`:

| ratio | seed state | `k01` profile | relative move |
|---|---|---|---|
| He/H | 9.010055601335061 | 9.010147506036333 | +1.020e-5 |
| C/H  | 2.7778794005952804e-4 | 2.7779522928812034e-4 | +2.624e-5 |
| N/H  | 8.190311782068561e-5 | 8.190526723261474e-5 | +2.624e-5 |
| O/H  | 8.864446896318066e-7 | 8.863431763672037e-7 | -1.145e-4 |

`solve_escape_wind` carried only helium onto the new reservoir, through
`map_state_to_grid.py --reservoir-HeH`, so the three metal ratios of the
seed's `# reservoir` line stayed at the previous iteration's values.
`load_IC.f90`'s `compare_reservoir_field` compares EVERY element of that
line against the run's at `heh_dev_tol = 1e-6` and stops the run on the
first departure, so C/H at 2.6e-5 was refused.

The loader's own renormalization of a handoff element (`load_IC.f90`, the
"elements the handoff states" loop) does not rescue this: it runs after the
metadata comparison, so a file carrying the metadata block reaches it
already within 1e-6 and its factor is 1 there. That renormalization bites on
a provenance-unknown file, one written before the metadata block, for which
no comparison is possible. The documentation said otherwise and is corrected
(section 5).

## 3. What changed

**`src/utils/map_state_to_grid.py`.** `--reservoir-HeH <value>` is REMOVED
and replaced by `--reservoir <El>/H <value>`, repeatable, one pair per
element:

```
map_state_to_grid.py <src> <target_grid> <out> --ic \
    --reservoir He/H 9.0101475060363327 --reservoir C/H 0.00027779522928812034
```

`El` is helium or one of the ten metals of `species_table`. The element's
ionization-stage columns (He I, He II, He III and the He 2^3S level; El I,
El II, El III) are multiplied in every row by k = value / (the El/H of the
source's `# reservoir` line), which leaves the ionization split and the
shape of the column untouched and moves only the normalization. The mass
density takes the added nuclei at the weights `calc_rho` uses -- 3.9715259
m_H per free helium nucleus, 4.9715259 per HeH+, and
`amu_over_m_H x A_El` per metal nucleus, `A_El` the standard atomic weight
in u of `species_table`'s `melem_A_u` -- the pressure is held, and the
temperature follows the new particle count, T = p/((n_tot + n_e) k_B), a
metal stage El^(k+) counting 1 + k. The `# reservoir` field of both output
files is rewritten with 17 significant digits, and the `# mapped:` line
records each element's factor and the measured ratio of the base rows and of
the column.

The rescaling of an element is refused when the source states no
`# reservoir` line, when that line does not name the element, when the
species file carries no column of it, when the value is not positive and
finite, or when the base rows do not arrive at the new ratio. The last test
is what catches a species holding nuclei of two elements at once: HeH+ holds
one helium and one hydrogen nucleus, and the oxygen-chemistry carriers OH,
H2O and CO hold oxygen and carbon, so no single factor on the ionization
stages moves the element count by the factor asked for. The base rows are
the two inner ghosts and cell 1 for helium, where the element-diffusion
operator holds its Dirichlet He/H, and cell 1 for a metal, the row
`load_IC.f90` forms its own renormalization factor at.

**`src/utils/element_flux_closure.py`.** `reservoir_HeH` became
`reservoir_ratios`, returning the whole `# reservoir` line as a dictionary,
and `profile_match_HeH` became `profile_match_ratios`, which finds every
`X_<El>` column BY NAME on the profile's `# columns:` header (no column
number is hardcoded) and returns the ratios at the matching row. The seeding
step compares the two element by element and carries every element whose
relative move exceeds 1e-6 in ONE mapper call, logging the moves:

```
  seed: reservoir C/H 0.00027778794 -> 0.000277795229, He/H 9.0100556 -> 9.01014751,
        N/H 8.19031178e-05 -> 8.19052672e-05, O/H 8.8644469e-07 -> 8.86343176e-07
```

**Callers and documentation** moved to the new spelling:
`LHS1140b/models/run_case.sh` (replaced atomically while the campaign was
running), `pick_seed.py`, `write_reproduce.py`, `LHS1140b/models/README.md`,
`LHS1140b/MODELS.md`, the `REPRODUCE.md` files that quote the command, and
`src/tests/state_mapper/run.sh`.

## 4. Tests and the reproduction run

`src/tests/state_mapper/run.sh`: **76 assertions, 0 failures** (46/0 before).
The helium assertions are the previous ones under the new spelling; the new
ones are a metal state carrying C and O stage columns and a `# reservoir`
line naming He/H, C/H, O/H and Mg/H:

- one element (C/H): only C I, C II, C III move, by the ratio, to 1e-15; H,
  He and O are bit-identical to the option-free mapping; C/H at the first
  physical cell is the new ratio to 1e-14 and is flat over the column;
  pressure and velocity bit-identical; T = p/((n_tot + n_e) k_B) and rho the
  species mass sum to 1e-13; the `# reservoir` line has C/H rewritten and
  He/H, O/H, Mg/H untouched;
- several elements in one call (C/H, O/H, He/H): each arrives at its new
  ratio to 1e-14, the state stays consistent, the line is rewritten in all
  three places and the `# mapped:` line records three rescalings;
- refusals: `Mg/H`, named by the reservoir line but carrying no column;
  `N/H`, carrying neither; a non-positive value; an element symbol that is
  not one of the eleven. Nothing is written.

The reproduction ran in `LHS1140b/models/.L13/HeH9/k01/`, seeded from
`atomic_photochem_gj1132_kzzprofile/HeH9/k00/output` with the `k01` profile
and `input.inp`, at `OMP_NUM_THREADS=4`. `load_IC` passes, and the loader's
own renormalization factor is exactly 1 for every metal:

```
 (load_IC) He_diffusion: the diffused He/H profile of the restart file is kept (top cell He/H = 2.5474E+00); base cells set to the input He/H = 9.0101E+00.
 (load_IC) C/H at the base cell: restart file  2.777952E-04, handoff  2.777952E-04; column rescaled by  1.000000E+00
 (load_IC) O/H at the base cell: restart file  8.863432E-07, handoff  8.863432E-07; column rescaled by  1.000000E+00
 (load_IC) N/H at the base cell: restart file  8.190527E-05, handoff  8.190527E-05; column rescaled by  1.000000E+00
```

The stationary solve then started (`(EXHALE_main) stationary restart`, JFNK
`||R|| = 8.713E-01`) and was stopped by PID; the run directory is a scratch
one and is not a rung.

## 5. A SECOND refusal on the same path, still open

With the reservoir refusal gone the next metadata field fails, and the first
reproduction attempt measured it:

```
 (load_IC) ERROR: metadata field "grid" differs between the restart files and this run.
   file: grid N 500 R0[cm] 1.1593574831171966E+09 r_min[Rp] 1.0001933187778382E+00 r_max[Rp] 2.9031224061195065E+01 mode Mixed
   run:  grid N 500 R0[cm] 1.1593574812390110E+09 r_min[Rp] 1.0001933187778382E+00 r_max[Rp] 2.9031224061195065E+01 mode Mixed
```

`R0` is the profile's radius at the matching level times R_J
(`apply_lower_atmosphere_profile`, `lap_value_at_match('r')`), and it moves
with the profile exactly as the elemental ratios do: `k00` states
r_match = 1.62166044189167552e-1 R_J and `k01` states
1.62166043926454839e-1 R_J, a relative move of 1.6e-9. The closure hands
`map_state_to_grid.py` the SEED's own `Hydro_ioniz.txt` as the target grid
file, so the seed keeps the seed's `# grid` line, and the loader refuses it.
The cell centers in R_p are unchanged; only the physical scale moves, by
1.6e-9, far below the state's own truncation error.

This is a separate defect from the reservoir one and is NOT fixed here: the
choice of what the closure should hand the mapper as the target grid (a
header copy carrying the run's R0, or a new option, or a relaxed grid
comparison) is a contract decision. The reproduction above was made
conclusive by writing `target_grid_Hydro_ioniz.txt`, the seed's file with
its `# grid` R0 replaced by the run's, and using that as the mapper target;
that file is in the scratch directory. **The closure rungs will not run
until this is decided.**

## 6. Documentation corrected

- `docs/input_schema.md` key 20 (`Load IC?`) and the `<El>_H_base` row: the
  renormalization on a restart is reached only by a provenance-unknown file,
  because the `# reservoir` metadata field is compared first, element by
  element at 1e-6, and a departure stops the run; `map_state_to_grid.py
  --reservoir <El>/H <value>` is what carries a state to another reservoir.
- `src/modules/files_IO/load_IC.f90`, the module header and the comments of
  the "elements the handoff states" loop: the same correction. Comment text
  only, with the line count unchanged.

## 7. Interface change for the user

`--reservoir-HeH <value>` no longer exists. Write
`--reservoir He/H <value>`, and add one `--reservoir <El>/H <value>` pair
per further element. `run_case.sh` and the closure driver add the pairs
themselves; nothing in a case configuration changes.

## 8. A build incident, recorded

The `load_IC.f90` comment edit was verified against the binary by rebuilding.
The first rebuild was made with `PATH=/usr/bin:$PATH make`, which is NOT the
toolchain this tree is built with: `/usr/bin/gfortran` is 13.1.0 and links
the Ubuntu reference LAPACK (`liblapack.so.3`, no rpath, with a
`libgfortran.so.4 may conflict with libgfortran.so.5` warning), while the
`gfortran` on the default PATH is conda-forge 16.2.0 and links
`-L/opt/miniconda3/lib -lopenblas -Wl,-rpath,/opt/miniconda3/lib`, which is
what the campaign binary `b9419ba7cd8e519ca86896b29c382762` was built from.
That build left seven objects compiled by gcc 13 in a gcc 16.2 tree and put
a 3273376-byte executable (`9e38a707fc8296d00e2ccf3dc96e0b2c`) in place of
the 3622336-byte one for about three minutes.

Repaired: the campaign binary was restored from a running process's own
image (`/proc/<pid>/exe`, md5 `b9419ba7`), the seven objects were removed
and rebuilt with the default toolchain, and the relink reproduces
`b9419ba7cd8e519ca86896b29c382762` exactly. `make -q` is clean. The comment
edit is therefore bit-neutral, which a direct test also shows: compiling the
text before and against the text after with the same compiler, flags and
file name gives byte-identical objects (`712b2c4e77e87284e4cefd3614948f83`
both ways).

ONE CASE ran on the wrong executable and its result must not be kept:
`LHS1140b/models/atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13`, started about
10:04. It was still running at the end of this work and could not be
stopped from here.

**Build note for this tree: use plain `make`.** `PATH=/usr/bin:$PATH make`
builds a different executable against a different BLAS.

## 9. Fixed in passing

`LHS1140b/models/make_models.py::profile_match_HeH`, which names a case by
the He/H its stored column carries, read `word[8]` -- the ninth field of the
row -- instead of finding `X_He` on the `# columns:` header. It now finds the
column by name and refuses a file that names none, as the code reading the
profile does. Measured on
`atomic_photochem_gj1132_kzzprofile/HeH9/k01/lower_atmosphere_profile.dat`:
9.010147506036333 before and after, so no case name moves. The comment above
it named a function `_profile_match_HeH` that does not exist; corrected.

## 10. One sweep left to repeat

`run_case.sh` was replaced atomically at 09:45; a case started before that
keeps the previous file open (bash reads a script as it runs) and writes the
old option name into its `REPRODUCE.md`. A retry re-execs the script by path
and takes the new one, so the last such file is written when the last case
started before 09:45 finishes. Repeat once then:

```

> Note added 2026-09-14 15:20: the private builds of this memo were made with `PATH=/usr/bin:$PATH make`, which selects the Ubuntu gfortran 13.1 and its LAPACK, not the conda-forge gfortran 16.2 + OpenBLAS the tree binary `EXHALE.x` is built with (bare `make`). The control and the measured build of the memo share one toolchain, so every comparison here stands; reproduce with bare `make` (`OBJDIR=... EXE=...` as written) to obtain the tree's toolchain.
cd EXHALE_v1.00
for f in $(grep -rl "reservoir-HeH" LHS1140b/models --include=REPRODUCE.md); do
  sed -e 's/`--reservoir-HeH /`--reservoir He\/H /g' \
      -e 's/--ic --reservoir-HeH /--ic --reservoir He\/H /g' "$f" > "$f.new" \
    && \mv -f "$f.new" "$f"
done
```

## 11. The `grid` field, closed

Decision (advisor, 2026-09-14): the restart contract is left as it stands --
the `grid` field is compared as TEXT (`load_IC.f90::refuse_unequal_field`) --
and the closure driver hands the mapper a target header carrying the R0 of
the run that will load the state.

**How the run forms R0, and the same arithmetic here.** With a profile in
use, `input_read::apply_lower_atmosphere_profile` sets `R0 = val` from
`lap_value_at_match('r')`, which is
`lower_atmosphere_profile::value_at_pressure(ic_r, lap_p_match_bar)`: a level
whose own `p` equals `p_match` returns that level's `r` bit for bit, and
otherwise the value is linear in ln p between the bracketing levels. Nothing
else touches `R0` in `src/modules/` until `input_read.f90:1767`,
`R0 = R0*RJ`, with `RJ = 7.1492d9` (`parameters.f90:835`). So R0 in cm is ONE
double multiply on the profile's radius at the matching level, and
`element_flux_closure.profile_values_at_match` + `RJ_CM` reproduce exactly
that, in that order. The metadata writes it through `load_IC::meta_num`, an
`ES23.16` field, which is `'%.16E'` in Python.

**Measured, bitwise.** Reproducing k00's own R0 from k00's profile gives
`1.1593574831171966E+09`, the character string the k00 run itself wrote into
its `# grid` line; reproducing k01's gives `1.1593574812390110E+09`, the
string the k01 run states. No interpolation enters either: the profile's
match row carries `p` equal to the `# p_match_bar` header, so
`value_at_pressure` takes its exact-level branch.

**What the driver does now.** `write_target_grid_header` copies the seed's
`Hydro_ioniz.txt` to `<iter>/target_grid_Hydro_ioniz.txt` with only the R0 of
the `# grid` line replaced, and that copy is the mapper's target argument, so
the mapper (which takes the `# grid` line from its target) writes the run's
line into the seed. N, r_min, r_max and the grid mode stay the seed's: the
grid is built in units of R_p and does not move with R0; were one of them to
move, `load_IC` would refuse the state and print both lines. The mapper is
now called when the reservoir moved OR when R0 moved, so the branch that
used to copy the seed unchanged can no longer leave a stale grid line behind.

**Reproduction** (`LHS1140b/models/.L13/HeH9/k01`, seeded from `k00/output`
through the driver's own helpers, `OMP_NUM_THREADS=4`):

```
  seed: reservoir C/H 0.00027778794 -> 0.000277795229, He/H 9.0100556 -> 9.01014751,
        N/H 8.19031178e-05 -> 8.19052672e-05, O/H 8.8644469e-07 -> 8.86343176e-07
  seed: grid R0[cm] 1.1593574831171966E+09 -> 1.1593574812390110E+09
```

The run loads it with no refusal of any field:

```
 (load_IC) He_diffusion: the diffused He/H profile of the restart file is kept (top cell He/H = 2.5474E+00); base cells set to the input He/H = 9.0101E+00.
 (load_IC) C/H at the base cell: restart file  2.777952E-04, handoff  2.777952E-04; column rescaled by  1.000000E+00
 (load_IC) O/H at the base cell: restart file  8.863432E-07, handoff  8.863432E-07; column rescaled by  1.000000E+00
 (load_IC) N/H at the base cell: restart file  8.190527E-05, handoff  8.190527E-05; column rescaled by  1.000000E+00
```

and the `# grid` line the run WRITES into `k01/output/Hydro_ioniz.txt` is
character-for-character the one the driver put into the seed (`diff` empty):

```
# grid N 500 R0[cm] 1.1593574812390110E+09 r_min[Rp] 1.0001933187778382E+00 r_max[Rp] 2.9031224061195065E+01 mode Mixed
```

The iteration then solved: `Restart intent: stationary -- the stationary
solve returned info = 0; output written.` `state_mapper` is unchanged at
76/0 and the binary is unchanged at
`b9419ba7cd8e519ca86896b29c382762`.

The `REPRODUCE.md` sweep of section 10 was run, skipping the directories of
cases that were running (the cwd of every `EXHALE.x` plus the case argument
of every `run_case.sh`); none of the remaining files carried the old option
name, and none does now.

## 12. The same field again, from the runner: a case-specific target header

Measured 2026-09-14 12:46: `atomic_photochem_gj1132x{0.15,0.10,0.01}_kzzprofile/HeH9.7`
died at `load_IC` on the `grid` field the moment they started:

```
 (load_IC) ERROR: metadata field "grid" differs between the restart files and this run.
   file: grid N 500 R0[cm] 1.1273716464000001E+09 ...
   run:  grid N 500 R0[cm] 1.1592997812185061E+09 ...
```

**The cause, corrected against the first reading.** `run_case.sh` handed
`map_state_to_grid.py` the shared `models/current_grid_Hydro_ioniz.txt` as its
target, always. That file states R0 = 1.1273716464000001E+09 cm, the planet
radius of a case with no lower-atmosphere profile. The mapper copies the
`# grid` line of its TARGET into the state it writes, so the seed carried the
planet radius whatever the seed itself was -- the seed of
`x0.15` was `atomic_photochem_gj1132x0.33_kzzprofile/HeH9.7/output`, itself a
profile case at the right R0. The refusal is therefore not "the seed was a
scalar case": it is every profile case whose seed carries a metadata block.
`x0.33`, `x0.30` and `x0.20` escaped it only because their seeds are ARCHIVED
states (`archive_20260830/exhale/refresh_j96/xuvkzz/C0p33`, `C0p20`), which
carry no metadata block at all, so `load_IC` compares no field; the states
they then wrote state the run's own R0, 1.1592997812185061E+09.

**One rule, one implementation.** The R0, the matching level and the metadata
number format moved out of `element_flux_closure.py` into
`src/utils/profile_match_level.py`, which both the closure driver (importing
it) and `run_case.sh` (running it) use:

```
profile_match_level.py <profile.dat> <grid_file> <out_file>
```

writes `<out_file>`, a copy of `<grid_file>` with the single `# grid` R0 field
replaced by the profile's, and prints that R0 as `meta_num` writes it. In
`run_case.sh` the seed block now reads `Lower atmosphere profile:` out of the
case's `input.inp`; where there is one, it writes
`target_grid_Hydro_ioniz.txt` in the case directory and maps against that,
and where there is none the shared file already states the run's R0 and is
used as it stands. `write_reproduce.py` takes `--grid-target` and `--grid-r0`
and `REPRODUCE.md` names the file and the number.

**Checked before the runs.** All three profiles give
`1.1592997812185061E+09`, character for character the R0 the three runs stated
in their refusal, against the shared file's `1.1273716464000001E+09`. The
refactor is behavior-preserving: rebuilding the L13 reproduction seed through
the moved helpers gives data rows byte-identical to the pair the run of
section 11 loaded, the only difference being the target path the `# mapped:`
provenance line records.

**The reruns** (`OMP_NUM_THREADS=8 ./run_case.sh <case>`, three at once, the
previous `output/`, `run.log`, `REPRODUCE.md`, `seed.log` and
`EXHALE_setup.out` moved to
`.stopped/<case>_gridR0_202609141259/`). All three load with no field refused:

```
 (load_IC) He_diffusion: the diffused He/H profile of the restart file is kept (top cell He/H = 5.4580E-01); base cells set to the input He/H = 9.7110E+00.
 (load_IC) C/H at the base cell: restart file  2.778017E-04, handoff  2.778017E-04; column rescaled by  1.000000E+00
 (load_IC) O/H at the base cell: restart file  9.027749E-07, handoff  9.027749E-07; column rescaled by  1.000000E+00
 (load_IC) N/H at the base cell: restart file  8.190725E-05, handoff  8.190725E-05; column rescaled by  1.000000E+00
```

and the `# grid` line of each seed is

```
# grid N 500 R0[cm] 1.1592997812185061E+09 r_min[Rp] 1.0001933187778382E+00 r_max[Rp] 2.9031224061195065E+01 mode Mixed
```

`pick_seed.py` chose `atomic_photochem_gj1132x0.20_kzzprofile/HeH9.7/output`
for all three this time, that case having certified at 12:59.

**What the solves then did.** The three runs launched here at 12:59 on
`EXHALE.x` (`b9419ba7cd8e519ca86896b29c382762`) reached outer pass 2 or 3 with
the hydrodynamic rows inside their tolerances, the gated elemental row alone
standing and the joint under-relaxation halving at every pass -- the L14 form:

```
x0.15  outer pass 2: hydro info=2, worst gated species row  1.80E-04 of  1.0E-05 at cell 365 (elemental transport He/H partition), mass  2.90E-08, momentum  3.21E-12, energy  2.87E-01, omega 0.250
x0.10  outer pass 3: hydro info=2, worst gated species row  3.28E-04 of  1.0E-05 at cell 368 (elemental transport He/H partition), mass  3.22E-08, momentum  3.21E-12, energy  5.68E-01, omega 0.125
x0.01  outer pass 2: hydro info=2, worst gated species row  2.94E-03 of  1.0E-05 at cell 379 (elemental transport He/H partition), mass  2.81E-08, momentum  3.21E-12, energy  1.01E+00, omega 0.250
```

At 13:35 the coordinating session (this campaign's main session) moved all three aside into
`.stopped/<case>_control_202609141335/` and relaunched them on ANOTHER binary,
md5 `dd7ee3188ead1cf54d2f2c7f99246309` -- neither `EXHALE.x` nor the
`EXHALE_L14.x` that stands now (`a11038c245050d4f11852af829e13fc7`) -- which
is the L14 progress control under measurement. Under it the hydrodynamic solve
returns `info = 0` from the first pass and the gated elemental row falls
monotonically instead of halving omega:

```
x0.15  outer pass 1: ... 7.52E-04 ... pass 9: hydro info=0, worst gated species row  7.49E-06 of  1.0E-05 at cell 368
       outer pass 9: ACCEPTED -- every active equation of this state is within its own tolerance.
x0.10  outer pass 1: ... 2.24E-03 ... pass 12: hydro info=0, worst gated species row  5.74E-06 of  1.0E-05 at cell 379
       outer pass 12: ACCEPTED -- every active equation of this state is within its own tolerance.
```

`x0.15` CERTIFIED at outer pass 9, log10 Mdot 7.04, 57m46s; `x0.10` CERTIFIED
at outer pass 12, log10 Mdot 6.85, 63m23s; both at `OMP_NUM_THREADS=8`,
13:35:23 to 14:33 and 14:38. `x0.01` was moved aside once more at 21:35
(`.stopped/..._farseed_202609142135/`, its target header among the files, so
it used the seeding of this section) after reaching outer pass 9 still refusing
-- gated row 3.02E-02, mass 1.12E-02, omega 0.125 -- and its case directory is
bare at the time of writing; that case is reached by the XUV ladder of item L14 (0.10 -> 0.05 -> 0.03 -> 0.02 -> 0.015 -> 0.01).

So the L14 form seen on `EXHALE.x` was the progress control and not the
seeding, and the certifications belong to L14's binary, not to this work. What
this section establishes is what it set out to: the `grid` field no longer
refuses these cases, they reach the solver, and two of the three then went to
CERTIFIED.

**Not rebuilt.** `make -q` reports the tree out of date at the end of this
work, from `src/EXHALE_main.f90` edited 11:52 -- the L14 progress-control
change, which was built into `EXHALE_L14.x` at the same minute and is being
measured there. Nothing here touched Fortran, and `EXHALE.x` is left at
`b9419ba7cd8e519ca86896b29c382762`, the binary the campaign is running; a
`make` from this session would have put the unfinished L14 change into it.

## 13. The `reservoir` field, from the runner: every element, not He/H alone

Measured 2026-09-16 00:42, in the re-run of the campaign on the corrected base
boundary (item L21). The six XUV-scaled cases above the fiducial photochemical
column died at `load_IC` the moment they started, all six the same way:

```
 (load_IC) ERROR: metadata field "reservoir": the elemental ratio C/H of the
   restart files is 2.778016888E-04
   this run asks for 2.778027137E-04
   relative departure 3.6893E-06, threshold 1.0000E-06
   The state is a solution at its own composition. Restart it there, or start
   this composition cold.
```

**The cause.** `load_IC::compare_reservoir_field` walks the `# reservoir` line
ELEMENT BY ELEMENT and refuses at `heh_dev_tol` = 1e-6 on any of them.
`run_case.sh` carried ONE element onto the seed, He/H, taken from
`He/H number ratio` in `input.inp`. That is the whole composition only where
the case states it in `input.inp`; a case with `Lower atmosphere profile:`
takes EVERY ratio from the profile's `X_<El>` column at the matching level
(`input_read::lap_element_ratio_at_match`), and these six cases read their
column from the flux-closed He/H = 9.7 rung, which the re-run had just
re-solved. He/H was carried and C/H, N/H and O/H were not, so the state
arrived at the run's helium and the previous column's carbon.

It is the same shape of defect as section 12 and the same cause: a field the
restart contract compares was formed in the runner from the case name instead
of from the profile the run actually reads. It had not shown before because
the profile a prescribed case reads had never changed under it.

**The fix.** The rule now lives once, in `src/utils/profile_match_level.py`,
beside the R0 rule that section 12 put there for the same reason:

- `profile_match_ratios(profile)` -- the reservoirs the run will ask for, from
  the `X_<El>` columns at the matching level;
- `reservoir_ratios(state_file)` -- the `# reservoir` line of a state;
- `reservoir_moves(state_file, ratios)` -- what differs by more than the
  contract's own 1e-6, as {El/H: (seed, run)};
- `profile_match_level.py --reservoir-options <state> <profile>` -- the same,
  as the `--reservoir <El>/H <value>` arguments `map_state_to_grid.py` takes.

`run_case.sh` calls that command for a case with a profile and passes every
option it prints; without a profile it keeps reading He/H from `input.inp`,
which is where the composition is stated in that case.
`element_flux_closure.py`, which had grown its own copy of the same rule while
seeding one closure iteration from the previous one, now imports
`profile_match_ratios` and `reservoir_moves` and computes nothing of its own.

**Verified.** The `--reservoir-options` command on the pre-L21 state of
`atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7` against its new column prints
four options, C/H, He/H, N/H and O/H; the mapped seed's `# reservoir` line
then carries all four at the profile's values, and all six cases load and
enter the solver. The nine closure rungs, which ran before this change on the
driver's own copy of the rule, are unaffected: the routine is the one they
were already using, moved.
