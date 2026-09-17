# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`.L4h/p002_HeH9.7`


The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L4h.x` |
| md5 | `9a944cdf573541a47a50ec9599703c07` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_L14_p002_photochem_HeH9.7` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L14/p003_HeH9.7/output (given)
```

The fields are the state directory, the tier it was taken from (`tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

The grid the seed is written on is named by `target_grid_Hydro_ioniz.txt`. This case builds its grid on a lower-atmosphere profile, so its radial scale R0 is the radius that profile carries at its matching level times R_J, 1.1592997812185061E+09 cm, and not the planet radius the shared `current_grid_Hydro_ioniz.txt` states. `load_IC` compares the `grid` metadata field as text, so the seed has to carry this number: the file is that shared one with the single field replaced, written by `src/utils/profile_match_level.py`.

It was interpolated onto the cell centers of the current code by

```
# the target grid header, this case's R0 in place of the shared file's
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/profile_match_level.py \
    lower_atmosphere_profile.dat /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt target_grid_Hydro_ioniz.txt
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L14/p003_HeH9.7/output \
    target_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
# the target grid header, this case's R0 in place of the shared file's
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/profile_match_level.py \
    lower_atmosphere_profile.dat /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt target_grid_Hydro_ioniz.txt
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L14/p003_HeH9.7/output \
    target_grid_Hydro_ioniz.txt output --ic 

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e8 \
    EXHALE_OUTER_PASSES=40 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L4h.x > run.log 2>&1
```

## What came out

The outer passes of the marching plus the JFNK hand-off, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 2 | 5.08E-02 | 2.86E-09 | 2.24E-01 | 2.94E-02 of 1.0E-05 (elemental transport He/H partition) | 371 | 1988.65 |
| 2 | 2 | 5.08E-02 | 2.86E-09 | 2.24E-01 | 2.72E-02 of 1.0E-05 (elemental transport He/H partition) | 288 | 730.85 |
| 3 | 2 | 5.08E-02 | 2.86E-09 | 2.24E-01 | 3.20E-02 of 1.0E-05 (elemental transport He/H partition) | 288 | 759.51 |
| 4 | 2 | 5.08E-02 | 2.86E-09 | 2.24E-01 | 3.25E-02 of 1.0E-05 (elemental transport He/H partition) | 288 | 821.51 |
| 5 | 2 | 5.08E-02 | 2.86E-09 | 2.24E-01 | 3.04E-02 of 1.0E-05 (elemental transport He/H partition) | 287 | 801.44 |
| 6 | 2 | 5.08E-02 | 2.86E-09 | 2.25E-01 | 2.67E-02 of 1.0E-05 (elemental transport He/H partition) | 286 | 789.91 |
| 7 | 2 | 5.08E-02 | 2.86E-09 | 2.25E-01 | 2.23E-02 of 1.0E-05 (elemental transport He/H partition) | 285 | 679.22 |
| 8 | 2 | 5.08E-02 | 2.86E-09 | 2.25E-01 | 1.78E-02 of 1.0E-05 (elemental transport He/H partition) | 284 | 832.22 |

- solver verdict: **info = 2** (last `||R||` = 2.246E-01), by the marching plus the JFNK hand-off
- certification of the state written: **NOT CERTIFIED**
  - hydrodynamic mass row: row measure  4.733E-02 above  2.0E-07 at cell 13
  - hydrodynamic energy row: row measure  2.246E-01 above  1.0E-06 at cell 26
  - elemental transport He/H partition: gated row measure  1.779E-02 above  1.0E-05 at cell 284 (a wind cell)
- wall clock **139m27s** at `OMP_NUM_THREADS=8`, 2026-09-15 14:22:12 to 2026-09-15 16:41:39

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
