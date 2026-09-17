# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`.L25/roe_p015_HeH9.7`

- numerical flux: `ROE`

The HLLC flux does not reach the stationary root on this grid at this XUV level (L25 step 3), and the same case solved with HLLC on a grid of twice the cells reaches the state the Roe flux reaches on the catalog grid, so the difference is resolution and not a different wind. The Roe flux systematic against HLLC where both solve is Mdot and the He 10830 equivalent width lower by 0.3 to 0.6 per cent, the first cell colder by 6 to 7 per cent and denser by 6 to 8 per cent, and the mass-flux spread unchanged (L25 memo). That systematic is carried wherever the numbers of this case are quoted.

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs two of its keys changed.

## Why this run was made the way it was

L25 step 3: the catalog recipe with `Numerical flux: ROE` as the only change. Run directory under models/.L25/, not the catalog.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `c2e9c9990b9f14f1be8cd77abca68945` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_photochem_gj1132x0.15_kzzprofile_HeH9.7` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/../models_20260914_preL21/atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7/output (given); then a continuation from the first solve's written state at dtau0=1.0e8 (its log run_dtau0_first.log)
```

The fields are the state directory, the tier it was taken from (`tier0` this same case's own most recent certified state -- the one this tree carries, else the one a preserved tree such as `models_20260914_preL21/` holds, which is what a re-solve of the whole catalog continues from -- `tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

The grid the seed is written on is named by `target_grid_Hydro_ioniz.txt`. This case builds its grid on a lower-atmosphere profile, so its radial scale R0 is the radius that profile carries at its matching level times R_J, 1.1592997810038679E+09 cm, and not the planet radius the shared `current_grid_Hydro_ioniz.txt` states. `load_IC` compares the `grid` metadata field as text, so the seed has to carry this number: the file is that shared one with the single field replaced, written by `src/utils/profile_match_level.py`.

It was interpolated onto the cell centers of the current code by

```
# the target grid header, this case's R0 in place of the shared file's
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/profile_match_level.py \
    lower_atmosphere_profile.dat /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt target_grid_Hydro_ioniz.txt
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/../models_20260914_preL21/atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7/output \
    target_grid_Hydro_ioniz.txt output --ic --reservoir C/H 0.00027780271379679837 --reservoir He/H 9.7110732820766632 --reservoir N/H 8.1907551439822828e-05 --reservoir O/H 9.027725250175568e-07 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
# the target grid header, this case's R0 in place of the shared file's
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/profile_match_level.py \
    lower_atmosphere_profile.dat /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt target_grid_Hydro_ioniz.txt
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/../models_20260914_preL21/atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7/output \
    target_grid_Hydro_ioniz.txt output --ic --reservoir C/H 0.00027780271379679837 --reservoir He/H 9.7110732820766632 --reservoir N/H 8.1907551439822828e-05 --reservoir O/H 9.027725250175568e-07 

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
| 1 | 2 | 1.36E+00 | 2.90E-05 | 7.77E-01 | 1.25E-06 of 1.0E-05 (elemental transport He/H partition) | 234 | 226.63 |
| 2 | 2 | 1.36E+00 | 2.90E-05 | 7.80E-01 | 6.48E-07 of 1.0E-05 (elemental transport He/H partition) | 230 | 255.33 |
| 3 | 2 | 1.36E+00 | 2.90E-05 | 7.80E-01 | 3.24E-07 of 1.0E-05 (elemental transport He/H partition) | 230 | 241.08 |
| 4 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 4.53E-07 of 1.0E-05 (elemental transport He/H partition) | 217 | 240.65 |
| 5 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 2.26E-07 of 1.0E-05 (elemental transport He/H partition) | 217 | 201.68 |
| 6 | 2 | 1.35E+00 | 2.90E-05 | 7.80E-01 | 1.05E-06 of 1.0E-05 (elemental transport He/H partition) | 234 | 251.15 |
| 7 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 8.85E-07 of 1.0E-05 (elemental transport He/H partition) | 224 | 256.09 |
| 8 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 6.64E-07 of 1.0E-05 (elemental transport He/H partition) | 224 | 189.74 |
| 9 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 4.98E-07 of 1.0E-05 (elemental transport He/H partition) | 224 | 353.11 |
| 10 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 3.73E-07 of 1.0E-05 (elemental transport He/H partition) | 224 | 286.69 |
| 11 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 2.80E-07 of 1.0E-05 (elemental transport He/H partition) | 224 | 205.68 |
| 12 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 2.10E-07 of 1.0E-05 (elemental transport He/H partition) | 224 | 360.57 |
| 13 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 1.58E-07 of 1.0E-05 (elemental transport He/H partition) | 224 | 213.96 |
| 14 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 1.18E-07 of 1.0E-05 (elemental transport He/H partition) | 224 | 228.38 |
| 15 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 8.88E-08 of 1.0E-05 (elemental transport He/H partition) | 224 | 211.05 |
| 16 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 6.66E-08 of 1.0E-05 (elemental transport He/H partition) | 224 | 228.89 |
| 17 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 5.00E-08 of 1.0E-05 (elemental transport He/H partition) | 224 | 266.42 |
| 18 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 3.75E-08 of 1.0E-05 (elemental transport He/H partition) | 224 | 199.54 |
| 19 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 2.81E-08 of 1.0E-05 (elemental transport He/H partition) | 224 | 221.04 |
| 20 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 2.11E-08 of 1.0E-05 (elemental transport He/H partition) | 224 | 255.56 |
| 21 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 1.58E-08 of 1.0E-05 (elemental transport He/H partition) | 224 | 208.89 |
| 22 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 1.19E-08 of 1.0E-05 (elemental transport He/H partition) | 224 | 237.02 |
| 23 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 8.91E-09 of 1.0E-05 (elemental transport He/H partition) | 224 | 229.87 |
| 24 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 6.69E-09 of 1.0E-05 (elemental transport He/H partition) | 224 | 217.72 |
| 25 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 5.02E-09 of 1.0E-05 (elemental transport He/H partition) | 224 | 201.03 |
| 26 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 3.76E-09 of 1.0E-05 (elemental transport He/H partition) | 224 | 215.80 |
| 27 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 2.82E-09 of 1.0E-05 (elemental transport He/H partition) | 224 | 284.31 |
| 28 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 2.12E-09 of 1.0E-05 (elemental transport He/H partition) | 224 | 215.79 |
| 29 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 1.59E-09 of 1.0E-05 (elemental transport He/H partition) | 224 | 212.61 |
| 30 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 1.19E-09 of 1.0E-05 (elemental transport He/H partition) | 224 | 253.81 |
| 31 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 8.95E-10 of 1.0E-05 (elemental transport He/H partition) | 224 | 196.85 |
| 32 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 6.72E-10 of 1.0E-05 (elemental transport He/H partition) | 224 | 285.86 |
| 33 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 5.04E-10 of 1.0E-05 (elemental transport He/H partition) | 224 | 202.93 |
| 34 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 3.78E-10 of 1.0E-05 (elemental transport He/H partition) | 224 | 204.78 |
| 35 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 2.84E-10 of 1.0E-05 (elemental transport He/H partition) | 224 | 233.84 |
| 36 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 2.13E-10 of 1.0E-05 (elemental transport He/H partition) | 224 | 203.00 |
| 37 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 1.60E-10 of 1.0E-05 (elemental transport He/H partition) | 224 | 209.16 |
| 38 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 1.21E-10 of 1.0E-05 (elemental transport He/H partition) | 224 | 192.89 |
| 39 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 9.11E-11 of 1.0E-05 (elemental transport He/H partition) | 224 | 178.63 |
| 40 | 2 | 1.35E+00 | 2.90E-05 | 7.81E-01 | 7.40E-11 of 1.0E-05 (elemental transport He/H partition) | 224 | 207.60 |

- solver verdict: **info = 1** (last `||R||` = 1.352E+00), by the partitioned stationary route
- certification of the state written: **NOT CERTIFIED**
  - hydrodynamic mass row: row measure  6.807E-05 above  3.1E-12 at cell 303
  - hydrodynamic momentum row: row measure  2.902E-05 above  1.0E-08 at cell 500
  - hydrodynamic energy row: row measure  7.814E-01 above  1.0E-06 at cell 1
- wall clock **499m20s** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-16 23:22:39 to 2026-09-17 07:41:59

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

## Re-measuring this state, and how closely it comes back

The two files of `output/` can be handed back to the binary as `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with `Restart intent: stationary evaluate`, which measures the state as it stands, writes it back unchanged and exits 0 only if every active equation is within its tolerance.

That is the post-processing pass of this record: it is the route the advection-corrected profiles and the mass-loss line above were produced on. Three states are named in it and every product says which one it describes: the LOADED state, the conserved variables the two files carry, which the pass does not write to; the WORK state, which holds those conserved variables fixed and derives the pressure and the temperature from the composition one equilibrium sweep returns, through the caloric equation of state, and which is what `Hydro_ioniz.txt`, `Ion_species.txt`, the breakdowns and the `certified=` pair of their header describe; and the ADVECTION-DERIVED composition the post-process builds from the work state, which is what `Hydro_ioniz_adv.txt` and `Ion_species_adv.txt` carry. The derived pair states the work state's certification pair on a `# derived_from:` line, as provenance of the state it was built from, and makes no certification claim about its own rows. One sweep is not a closed chemical and thermal fixed point, so the work state reports its own closure defect in `pp.log` and "refreshed" never means "closed".

Two answers are kept apart in `pp.log` and in the list above: whether the stationary claim the solved state was written with reproduces when that state is handed back, and what verdict the work state gets. A passing work state is not a confirmation of a claim that did not reproduce.

What comes back exactly and what does not, MEASURED over the 120 certified states of `models/` (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8). The conserved state round-trips to the last bit: the radius, velocity and pressure columns are read as written, and since item L18 the mass density is read from its own column rather than rebuilt from the species, so the state re-measured is the state certified. The ROW MEASURES do not round-trip bitwise and cannot: the first equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at its rounding floor, where one unit in the last place of the density moves the cell-wise maximum of the energy row by about 11 per cent. Re-evaluated against in-run, the mass row comes back within a factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to 1.93 (median 1.22); the momentum row sits four to six decades below its tolerance and its ratio is the floor itself. A state certified at more than about half of its tolerance may therefore be refused on re-evaluation, and that refusal is a property of the cancellation and not of this file.
