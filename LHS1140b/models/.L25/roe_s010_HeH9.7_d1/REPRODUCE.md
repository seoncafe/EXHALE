# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`.L25/roe_s010_HeH9.7_d1`

- numerical flux: `ROE`

The HLLC flux does not reach the stationary root on this grid at this XUV level (L25 step 3), and the same case solved with HLLC on a grid of twice the cells reaches the state the Roe flux reaches on the catalog grid, so the difference is resolution and not a different wind. The Roe flux systematic against HLLC where both solve is Mdot and the He 10830 equivalent width lower by 0.3 to 0.6 per cent, the first cell colder by 6 to 7 per cent and denser by 6 to 8 per cent, and the mass-flux spread unchanged (L25 memo). That systematic is carried wherever the numbers of this case are quoted.

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs two of its keys changed.

## Why this run was made the way it was

L25 step 3: the same case at the default pseudo-time start, seeded from the stopped catalog attempt's written state, which is the configuration step 2 measured the Roe convergence in.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `c2e9c9990b9f14f1be8cd77abca68945` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_scalar_gj1132x0.10_kzz1e9_HeH9.7` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.stopped/atomic_scalar_gj1132x0.10_kzz1e9_HeH9.7_dtau0-1_20260916055251/output (given)
```

The fields are the state directory, the tier it was taken from (`tier0` this same case's own most recent certified state -- the one this tree carries, else the one a preserved tree such as `models_20260914_preL21/` holds, which is what a re-solve of the whole catalog continues from -- `tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

The grid the seed is written on is named by `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt`.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.stopped/atomic_scalar_gj1132x0.10_kzz1e9_HeH9.7_dtau0-1_20260916055251/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.stopped/atomic_scalar_gj1132x0.10_kzz1e9_HeH9.7_dtau0-1_20260916055251/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 \
    EXHALE_OUTER_PASSES=40 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > run.log 2>&1
```

## What came out

The outer passes of the marching plus the JFNK hand-off, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 5.86E-08 | 1.03E-14 | 2.36E-07 | 2.51E-04 of 1.0E-05 (elemental transport He/H partition) | 217 | 320.37 |
| 2 | 0 | 2.04E-09 | 2.54E-14 | 2.10E-08 | 1.98E-04 of 1.0E-05 (elemental transport He/H partition) | 358 | 268.16 |
| 3 | 0 | 1.89E-09 | 3.48E-14 | 2.17E-08 | 1.59E-04 of 1.0E-05 (elemental transport He/H partition) | 361 | 220.25 |
| 4 | 2 | 2.11E-09 | 3.48E-14 | 1.82E-03 | 1.03E-04 of 1.0E-05 (elemental transport He/H partition) | 363 | 1295.47 |
| 5 | 2 | 1.89E-09 | 3.48E-14 | 2.77E-03 | 6.12E-05 of 1.0E-05 (elemental transport He/H partition) | 364 | 1214.52 |
| 6 | 2 | 1.67E-09 | 3.48E-14 | 3.26E-03 | 3.44E-05 of 1.0E-05 (elemental transport He/H partition) | 366 | 2062.27 |
| 7 | 2 | 2.47E-09 | 3.48E-14 | 3.51E-03 | 1.88E-05 of 1.0E-05 (elemental transport He/H partition) | 367 | 1954.35 |
| 8 | 2 | 2.47E-09 | 3.48E-14 | 3.64E-03 | 1.01E-05 of 1.0E-05 (elemental transport He/H partition) | 368 | 1608.93 |
| 9 | 2 | 4.81E-09 | 3.48E-14 | 3.71E-03 | 5.35E-06 of 1.0E-05 (elemental transport He/H partition) | 369 | 1418.62 |
| 10 | 2 | 1.89E-09 | 3.48E-14 | 3.74E-03 | 2.82E-06 of 1.0E-05 (elemental transport He/H partition) | 370 | 1353.88 |
| 11 | 2 | 1.89E-09 | 3.48E-14 | 3.76E-03 | 1.48E-06 of 1.0E-05 (elemental transport He/H partition) | 371 | 1385.03 |
| 12 | 2 | 1.89E-09 | 3.48E-14 | 3.77E-03 | 7.73E-07 of 1.0E-05 (elemental transport He/H partition) | 371 | 1386.00 |
| 13 | 2 | 1.12E-09 | 3.48E-14 | 3.77E-03 | 4.04E-07 of 1.0E-05 (elemental transport He/H partition) | 372 | 1382.76 |
| 14 | 2 | 1.89E-09 | 3.48E-14 | 3.77E-03 | 2.10E-07 of 1.0E-05 (elemental transport He/H partition) | 372 | 1406.54 |
| 15 | 2 | 5.04E-09 | 3.48E-14 | 3.77E-03 | 1.10E-07 of 1.0E-05 (elemental transport He/H partition) | 373 | 1382.73 |
| 16 | 2 | 2.11E-09 | 3.48E-14 | 3.78E-03 | 5.69E-08 of 1.0E-05 (elemental transport He/H partition) | 373 | 1341.33 |
| 17 | 2 | 3.51E-09 | 3.48E-14 | 3.78E-03 | 2.96E-08 of 1.0E-05 (elemental transport He/H partition) | 374 | 1279.08 |
| 18 | 2 | 3.63E-09 | 3.48E-14 | 3.78E-03 | 1.53E-08 of 1.0E-05 (elemental transport He/H partition) | 374 | 1195.56 |
| 19 | 2 | 2.81E-09 | 3.48E-14 | 3.78E-03 | 7.93E-09 of 1.0E-05 (elemental transport He/H partition) | 375 | 1215.70 |
| 20 | 2 | 2.11E-09 | 3.48E-14 | 3.78E-03 | 4.10E-09 of 1.0E-05 (elemental transport He/H partition) | 375 | 1216.34 |
| 21 | 2 | 3.51E-09 | 3.48E-14 | 3.78E-03 | 2.12E-09 of 1.0E-05 (elemental transport He/H partition) | 375 | 1216.33 |
| 22 | 2 | 2.81E-09 | 3.48E-14 | 3.78E-03 | 1.09E-09 of 1.0E-05 (elemental transport He/H partition) | 376 | 1173.49 |
| 23 | 2 | 3.51E-09 | 3.48E-14 | 3.78E-03 | 5.62E-10 of 1.0E-05 (elemental transport He/H partition) | 376 | 1156.22 |
| 24 | 2 | 2.81E-09 | 3.48E-14 | 3.78E-03 | 2.88E-10 of 1.0E-05 (elemental transport He/H partition) | 376 | 1148.42 |
| 25 | 2 | 3.51E-09 | 3.48E-14 | 3.78E-03 | 1.47E-10 of 1.0E-05 (elemental transport He/H partition) | 376 | 1016.32 |
| 26 | 2 | 6.44E-09 | 3.48E-14 | 3.78E-03 | 7.52E-11 of 1.0E-05 (elemental transport He/H partition) | 377 | 739.61 |
| 27 | 2 | 4.33E-09 | 3.48E-14 | 3.78E-03 | 3.83E-11 of 1.0E-05 (elemental transport He/H partition) | 377 | 739.82 |
| 28 | 2 | 6.44E-09 | 3.48E-14 | 3.78E-03 | 1.95E-11 of 1.0E-05 (elemental transport He/H partition) | 377 | 738.47 |
| 29 | 2 | 2.93E-09 | 3.48E-14 | 3.78E-03 | 9.99E-12 of 1.0E-05 (elemental transport He/H partition) | 377 | 739.48 |
| 30 | 2 | 3.63E-09 | 3.48E-14 | 3.78E-03 | 5.14E-12 of 1.0E-05 (elemental transport He/H partition) | 378 | 739.74 |
| 31 | 2 | 5.04E-09 | 3.48E-14 | 3.78E-03 | 2.88E-12 of 1.0E-05 (elemental transport He/H partition) | 378 | 739.35 |
| 32 | 2 | 5.74E-09 | 3.48E-14 | 3.78E-03 | 2.84E-12 of 1.0E-05 (elemental transport He/H partition) | 378 | 740.04 |
| 33 | 2 | 3.63E-09 | 3.48E-14 | 3.78E-03 | 2.80E-12 of 1.0E-05 (elemental transport He/H partition) | 378 | 738.27 |
| 34 | 2 | 3.63E-09 | 3.48E-14 | 3.78E-03 | 2.77E-12 of 1.0E-05 (elemental transport He/H partition) | 378 | 739.55 |
| 35 | 2 | 4.33E-09 | 3.48E-14 | 3.78E-03 | 2.73E-12 of 1.0E-05 (elemental transport He/H partition) | 379 | 737.13 |

- solver verdict: **info = 2** (last `||R||` = 3.776E-03), by the marching plus the JFNK hand-off
- certification of the state written: **NOT CERTIFIED**
  - hydrodynamic energy row: row measure  3.776E-03 above  1.0E-06 at cell 387
- wall clock **641m16s** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-16 23:23:40 to 2026-09-17 10:04:56

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

## Re-measuring this state, and how closely it comes back

The two files of `output/` can be handed back to the binary as `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with `Restart intent: stationary evaluate`, which measures the state as it stands, writes it back unchanged and exits 0 only if every active equation is within its tolerance.

That is the post-processing pass of this record: it is the route the advection-corrected profiles and the mass-loss line above were produced on. Three states are named in it and every product says which one it describes: the LOADED state, the conserved variables the two files carry, which the pass does not write to; the WORK state, which holds those conserved variables fixed and derives the pressure and the temperature from the composition one equilibrium sweep returns, through the caloric equation of state, and which is what `Hydro_ioniz.txt`, `Ion_species.txt`, the breakdowns and the `certified=` pair of their header describe; and the ADVECTION-DERIVED composition the post-process builds from the work state, which is what `Hydro_ioniz_adv.txt` and `Ion_species_adv.txt` carry. The derived pair states the work state's certification pair on a `# derived_from:` line, as provenance of the state it was built from, and makes no certification claim about its own rows. One sweep is not a closed chemical and thermal fixed point, so the work state reports its own closure defect in `pp.log` and "refreshed" never means "closed".

Two answers are kept apart in `pp.log` and in the list above: whether the stationary claim the solved state was written with reproduces when that state is handed back, and what verdict the work state gets. A passing work state is not a confirmation of a claim that did not reproduce.

What comes back exactly and what does not, MEASURED over the 120 certified states of `models/` (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8). The conserved state round-trips to the last bit: the radius, velocity and pressure columns are read as written, and since item L18 the mass density is read from its own column rather than rebuilt from the species, so the state re-measured is the state certified. The ROW MEASURES do not round-trip bitwise and cannot: the first equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at its rounding floor, where one unit in the last place of the density moves the cell-wise maximum of the energy row by about 11 per cent. Re-evaluated against in-run, the mass row comes back within a factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to 1.93 (median 1.22); the momentum row sits four to six decades below its tolerance and its ratio is the floor itself. A state certified at more than about half of its tolerance may therefore be refused on re-evaluation, and that refusal is a property of the cancellation and not of this file.
