# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: the Photochem column handed over as `Lower atmosphere profile:` (`photochem`)
- spectrum: `gj1132x0.25`
- mixing: binary H/He element diffusion with K_zz(p) from the profile (`kzzprofile`)

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## Why this run was made the way it was

The first solve of this re-run was taken at the default pseudo-time start (`EXHALE_PTC_DTAU0=1.0`) and refused: forty outer passes with the gated species row already at 1e-12 and the hydrodynamic ENERGY row of the outermost physical cell above its tolerance (item L4f), which is the ramp collapse of item L4h. It was stopped, and its logs are in `models/.stopped/atomic_photochem_gj1132x0.25_kzzprofile_HeH9.7_dtau0-1_20260916055251/`. This run is what `run_case.sh` prescribes for that situation (item L14): the same mapped seed at `EXHALE_PTC_DTAU0=1.0e8`.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `db87b88d1ce53facf1d61084fa535ca5` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_photochem_gj1132x0.25_kzzprofile_HeH9.7` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models_20260914_preL21/atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7/output  tier0:models_20260914_preL21/atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7  HeH=9.71107  target=9.71107  dlog10=0.0000  candidates=557  (this same case's own most recent certified state); then a continuation from the first solve's written state at dtau0=1.0e8 (its log run_dtau0_first.log)
```

The fields are the state directory, the tier it was taken from (`tier0` this same case's own most recent certified state -- the one this tree carries, else the one a preserved tree such as `models_20260914_preL21/` holds, which is what a re-solve of the whole catalog continues from -- `tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

The grid the seed is written on is named by `target_grid_Hydro_ioniz.txt`. This case builds its grid on a lower-atmosphere profile, so its radial scale R0 is the radius that profile carries at its matching level times R_J, 1.1592997810042915E+09 cm, and not the planet radius the shared `current_grid_Hydro_ioniz.txt` states. `load_IC` compares the `grid` metadata field as text, so the seed has to carry this number: the file is that shared one with the single field replaced, written by `src/utils/profile_match_level.py`.

It was interpolated onto the cell centers of the current code by

```
# the target grid header, this case's R0 in place of the shared file's
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/profile_match_level.py \
    lower_atmosphere_profile.dat /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt target_grid_Hydro_ioniz.txt
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models_20260914_preL21/atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7/output \
    target_grid_Hydro_ioniz.txt output --ic --reservoir C/H 0.00027780271366409571 --reservoir He/H 9.7110732820759225 --reservoir N/H 8.190755142165824e-05 --reservoir O/H 9.0277252491333903e-07 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
# the target grid header, this case's R0 in place of the shared file's
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/profile_match_level.py \
    lower_atmosphere_profile.dat /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt target_grid_Hydro_ioniz.txt
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models_20260914_preL21/atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7/output \
    target_grid_Hydro_ioniz.txt output --ic --reservoir C/H 0.00027780271366409571 --reservoir He/H 9.7110732820759225 --reservoir N/H 8.190755142165824e-05 --reservoir O/H 9.0277252491333903e-07 

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e8 \
    EXHALE_OUTER_PASSES=40 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > run.log 2>&1

# the first solve returned info=1: the written state reloaded and solved
# again with the pseudo-time start raised (L4e)
cp output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
cp output/Ion_species.txt output/Ion_species_IC.txt
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e8 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > run.log 2>&1
```

## What came out

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 7.97E-09 of 1.0E-05 (elemental transport He/H partition) | 221 | 51.06 |
| 2 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 3.99E-09 of 1.0E-05 (elemental transport He/H partition) | 221 | 48.46 |
| 3 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 2.00E-09 of 1.0E-05 (elemental transport He/H partition) | 221 | 47.41 |
| 4 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 1.00E-09 of 1.0E-05 (elemental transport He/H partition) | 221 | 46.39 |
| 5 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 5.03E-10 of 1.0E-05 (elemental transport He/H partition) | 221 | 47.35 |
| 6 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 2.52E-10 of 1.0E-05 (elemental transport He/H partition) | 221 | 47.98 |
| 7 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 1.26E-10 of 1.0E-05 (elemental transport He/H partition) | 221 | 46.56 |
| 8 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 7.81E-11 of 1.0E-05 (elemental transport He/H partition) | 221 | 47.18 |
| 9 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 4.83E-11 of 1.0E-05 (elemental transport He/H partition) | 221 | 46.98 |
| 10 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 2.94E-11 of 1.0E-05 (elemental transport He/H partition) | 221 | 46.15 |
| 11 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 1.75E-11 of 1.0E-05 (elemental transport He/H partition) | 221 | 48.02 |
| 12 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 1.22E-11 of 1.0E-05 (elemental transport He/H partition) | 224 | 49.35 |
| 13 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 9.76E-12 of 1.0E-05 (elemental transport He/H partition) | 224 | 47.83 |
| 14 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 9.76E-12 of 1.0E-05 (elemental transport He/H partition) | 224 | 47.55 |
| 15 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 9.76E-12 of 1.0E-05 (elemental transport He/H partition) | 224 | 46.95 |
| 16 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 9.76E-12 of 1.0E-05 (elemental transport He/H partition) | 224 | 47.91 |
| 17 | 2 | 1.40E+00 | 6.60E-04 | 1.56E+00 | 9.76E-12 of 1.0E-05 (elemental transport He/H partition) | 224 | 48.64 |

- solver verdict: **info = 1** (last `||R||` = 1.561E+00), by the partitioned stationary route
- certification of the state written: **NOT CERTIFIED**
  - hydrodynamic mass row: row measure  4.538E-02 above  3.3E-09 at cell 59
  - hydrodynamic momentum row: row measure  6.604E-04 above  1.0E-08 at cell 500
  - hydrodynamic energy row: row measure  1.561E+00 above  1.0E-06 at cell 4
- wall clock **168m30s** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-16 05:53:44 to 2026-09-16 08:42:14

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

## Re-measuring this state, and how closely it comes back

The two files of `output/` can be handed back to the binary as `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with `Restart intent: stationary evaluate`, which measures the state as it stands, writes it back unchanged and exits 0 only if every active equation is within its tolerance.

What comes back exactly and what does not, MEASURED over the 120 certified states of `models/` (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8). The conserved state round-trips to the last bit: the radius, velocity and pressure columns are read as written, and since item L18 the mass density is read from its own column rather than rebuilt from the species, so the state re-measured is the state certified. The ROW MEASURES do not round-trip bitwise and cannot: the first equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at its rounding floor, where one unit in the last place of the density moves the cell-wise maximum of the energy row by about 11 per cent. Re-evaluated against in-run, the mass row comes back within a factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to 1.93 (median 1.22); the momentum row sits four to six decades below its tolerance and its ratio is the floor itself. A state certified at more than about half of its tolerance may therefore be refused on re-evaluation, and that refusal is a property of the cancellation and not of this file.
