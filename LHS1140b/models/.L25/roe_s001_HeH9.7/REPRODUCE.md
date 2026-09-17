# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`.L25/roe_s001_HeH9.7`

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
| run title | `Simulation for LHS1140b_atomic_scalar_gj1132x0.01_kzz1e9_HeH9.7` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.15_kzz1e9/HeH9.7/output (given); then a continuation from the first solve's written state at dtau0=1.0e8 (its log run_dtau0_first.log)
```

The fields are the state directory, the tier it was taken from (`tier0` this same case's own most recent certified state -- the one this tree carries, else the one a preserved tree such as `models_20260914_preL21/` holds, which is what a re-solve of the whole catalog continues from -- `tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

The grid the seed is written on is named by `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt`.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.15_kzz1e9/HeH9.7/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.15_kzz1e9/HeH9.7/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 

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
| 1 | 2 | 1.79E+00 | 4.58E-06 | 1.93E+00 | 2.65E-02 of 1.0E-05 (elemental transport He/H partition) | 400 | 159.92 |
| 2 | 2 | 1.74E+00 | 4.44E-06 | 1.93E+00 | 2.53E-02 of 1.0E-05 (elemental transport He/H partition) | 403 | 281.79 |
| 3 | 2 | 1.74E+00 | 4.44E-06 | 1.93E+00 | 2.50E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 200.52 |
| 4 | 2 | 1.73E+00 | 4.45E-06 | 1.92E+00 | 2.48E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 317.56 |
| 5 | 2 | 1.73E+00 | 4.45E-06 | 1.92E+00 | 2.45E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 175.80 |
| 6 | 2 | 1.73E+00 | 4.45E-06 | 1.92E+00 | 2.39E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 184.06 |
| 7 | 2 | 1.71E+00 | 4.46E-06 | 1.85E+00 | 2.28E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 307.12 |
| 8 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 2.09E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 303.46 |
| 9 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 1.79E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 288.49 |
| 10 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 1.38E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 280.29 |
| 11 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 9.56E-03 of 1.0E-05 (elemental transport He/H partition) | 500 | 267.65 |
| 12 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 5.90E-03 of 1.0E-05 (elemental transport He/H partition) | 500 | 187.44 |
| 13 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 3.34E-03 of 1.0E-05 (elemental transport He/H partition) | 500 | 193.83 |
| 14 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 1.79E-03 of 1.0E-05 (elemental transport He/H partition) | 500 | 200.13 |
| 15 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 9.28E-04 of 1.0E-05 (elemental transport He/H partition) | 500 | 277.60 |
| 16 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 4.73E-04 of 1.0E-05 (elemental transport He/H partition) | 500 | 190.26 |
| 17 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 2.39E-04 of 1.0E-05 (elemental transport He/H partition) | 500 | 292.12 |
| 18 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 1.20E-04 of 1.0E-05 (elemental transport He/H partition) | 500 | 313.09 |
| 19 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 6.01E-05 of 1.0E-05 (elemental transport He/H partition) | 500 | 273.27 |
| 20 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 3.01E-05 of 1.0E-05 (elemental transport He/H partition) | 500 | 261.66 |
| 21 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 1.51E-05 of 1.0E-05 (elemental transport He/H partition) | 500 | 310.56 |
| 22 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 7.53E-06 of 1.0E-05 (elemental transport He/H partition) | 500 | 168.18 |
| 23 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 3.77E-06 of 1.0E-05 (elemental transport He/H partition) | 500 | 273.12 |
| 24 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 1.88E-06 of 1.0E-05 (elemental transport He/H partition) | 500 | 260.05 |
| 25 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 9.41E-07 of 1.0E-05 (elemental transport He/H partition) | 500 | 260.34 |
| 26 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 4.71E-07 of 1.0E-05 (elemental transport He/H partition) | 500 | 195.45 |
| 27 | 2 | 1.71E+00 | 4.45E-06 | 1.91E+00 | 2.35E-07 of 1.0E-05 (elemental transport He/H partition) | 500 | 158.84 |
| 28 | 2 | 1.79E+00 | 4.23E-06 | 1.88E+00 | 1.07E-03 of 1.0E-05 (elemental transport He/H partition) | 277 | 181.54 |
| 29 | 2 | 1.88E+00 | 4.17E-06 | 1.91E+00 | 7.29E-04 of 1.0E-05 (elemental transport He/H partition) | 274 | 255.15 |
| 30 | 2 | 1.88E+00 | 4.17E-06 | 1.91E+00 | 3.83E-04 of 1.0E-05 (elemental transport He/H partition) | 279 | 222.47 |
| 31 | 2 | 1.88E+00 | 4.17E-06 | 1.91E+00 | 1.80E-04 of 1.0E-05 (elemental transport He/H partition) | 279 | 231.01 |
| 32 | 2 | 1.88E+00 | 4.17E-06 | 1.91E+00 | 1.30E-04 of 1.0E-05 (elemental transport He/H partition) | 257 | 230.14 |
| 33 | 2 | 1.90E+00 | 4.21E-06 | 1.85E+00 | 1.32E-04 of 1.0E-05 (elemental transport He/H partition) | 257 | 421.29 |
| 34 | 2 | 1.90E+00 | 4.21E-06 | 1.85E+00 | 1.21E-04 of 1.0E-05 (elemental transport He/H partition) | 262 | 230.19 |
| 35 | 2 | 1.90E+00 | 4.21E-06 | 1.85E+00 | 1.06E-04 of 1.0E-05 (elemental transport He/H partition) | 262 | 227.79 |
| 36 | 2 | 1.90E+00 | 4.21E-06 | 1.85E+00 | 8.78E-05 of 1.0E-05 (elemental transport He/H partition) | 264 | 243.51 |
| 37 | 2 | 1.90E+00 | 4.21E-06 | 1.85E+00 | 6.96E-05 of 1.0E-05 (elemental transport He/H partition) | 264 | 241.26 |
| 38 | 2 | 1.93E+00 | 3.95E-06 | 1.93E+00 | 2.11E-04 of 1.0E-05 (elemental transport He/H partition) | 229 | 150.13 |
| 39 | 2 | 1.93E+00 | 3.95E-06 | 1.93E+00 | 1.22E-04 of 1.0E-05 (elemental transport He/H partition) | 276 | 239.48 |
| 40 | 2 | 1.94E+00 | 3.90E-06 | 1.95E+00 | 1.31E-04 of 1.0E-05 (elemental transport He/H partition) | 280 | 266.11 |

- solver verdict: **info = 1** (last `||R||` = 1.952E+00), by the partitioned stationary route
- certification of the state written: **NOT CERTIFIED**
  - hydrodynamic mass row: row measure  2.460E-05 above  3.0E-12 at cell 495
  - hydrodynamic momentum row: row measure  3.900E-06 above  1.0E-08 at cell 499
  - hydrodynamic energy row: row measure  1.952E+00 above  1.0E-06 at cell 3
  - elemental transport He/H partition: gated row measure  1.309E-04 above  1.0E-05 at cell 280 (a wind cell)
- wall clock **436m23s** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-16 23:22:41 to 2026-09-17 06:39:04

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

## Re-measuring this state, and how closely it comes back

The two files of `output/` can be handed back to the binary as `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with `Restart intent: stationary evaluate`, which measures the state as it stands, writes it back unchanged and exits 0 only if every active equation is within its tolerance.

That is the post-processing pass of this record: it is the route the advection-corrected profiles and the mass-loss line above were produced on. Three states are named in it and every product says which one it describes: the LOADED state, the conserved variables the two files carry, which the pass does not write to; the WORK state, which holds those conserved variables fixed and derives the pressure and the temperature from the composition one equilibrium sweep returns, through the caloric equation of state, and which is what `Hydro_ioniz.txt`, `Ion_species.txt`, the breakdowns and the `certified=` pair of their header describe; and the ADVECTION-DERIVED composition the post-process builds from the work state, which is what `Hydro_ioniz_adv.txt` and `Ion_species_adv.txt` carry. The derived pair states the work state's certification pair on a `# derived_from:` line, as provenance of the state it was built from, and makes no certification claim about its own rows. One sweep is not a closed chemical and thermal fixed point, so the work state reports its own closure defect in `pp.log` and "refreshed" never means "closed".

Two answers are kept apart in `pp.log` and in the list above: whether the stationary claim the solved state was written with reproduces when that state is handed back, and what verdict the work state gets. A passing work state is not a confirmation of a claim that did not reproduce.

What comes back exactly and what does not, MEASURED over the 120 certified states of `models/` (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8). The conserved state round-trips to the last bit: the radius, velocity and pressure columns are read as written, and since item L18 the mass density is read from its own column rather than rebuilt from the species, so the state re-measured is the state certified. The ROW MEASURES do not round-trip bitwise and cannot: the first equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at its rounding floor, where one unit in the last place of the density moves the cell-wise maximum of the energy row by about 11 per cent. Re-evaluated against in-run, the mass row comes back within a factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to 1.93 (median 1.22); the momentum row sits four to six decades below its tolerance and its ratio is the floor itself. A state certified at more than about half of its tolerance may therefore be refused on re-evaluation, and that refusal is a property of the cancellation and not of this file.
