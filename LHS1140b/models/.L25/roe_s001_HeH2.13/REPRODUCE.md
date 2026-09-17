# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`.L25/roe_s001_HeH2.13`


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
| run title | `Simulation for LHS1140b_atomic_scalar_gj1132x0.01_kzz1e9_HeH2.13` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13/output (given)
```

The fields are the state directory, the tier it was taken from (`tier0` this same case's own most recent certified state -- the one this tree carries, else the one a preserved tree such as `models_20260914_preL21/` holds, which is what a re-solve of the whole catalog continues from -- `tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

The grid the seed is written on is named by `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt`.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e8 \
    EXHALE_OUTER_PASSES=40 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > run.log 2>&1

# the advection-corrected profiles: the solved state is handed back, measured
# as it stands with no step taken, and the products are written from it
for f in Hydro_ioniz Ion_species; do cp -f output/$f.txt output/${f}_IC.txt; done
cp -f input.inp input.inp.solved
sed -i -e 's/^Load IC?.*/Load IC? True/' -e '/^Restart intent:/d' input.inp
echo 'Restart intent: stationary evaluate' >> input.inp
OMP_NUM_THREADS=8 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > pp.log 2>&1
mv -f input.inp.solved input.inp

# the transit spectrum, through the WINERED HIRES-Y kernel of the measurement
. /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/winered_hires_y.sh
MPLBACKEND=Agg PYTHONPATH=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00 python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_transit.py > transit.log 2>&1
```

## What came out

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 2.01E-07 | 1.19E-14 | 3.58E-07 | 5.07E-02 of 1.0E-05 (elemental transport He/H partition) | 315 | 279.54 |
| 2 | 0 | 5.76E-07 | 1.28E-14 | 6.52E-07 | 5.00E-02 of 1.0E-05 (elemental transport He/H partition) | 316 | 154.54 |
| 3 | 0 | 1.82E-07 | 2.07E-14 | 5.97E-07 | 4.88E-02 of 1.0E-05 (elemental transport He/H partition) | 318 | 150.78 |
| 4 | 0 | 7.47E-08 | 1.07E-14 | 3.38E-07 | 4.68E-02 of 1.0E-05 (elemental transport He/H partition) | 323 | 155.98 |
| 5 | 0 | 9.86E-08 | 1.34E-14 | 3.93E-07 | 4.41E-02 of 1.0E-05 (elemental transport He/H partition) | 329 | 137.03 |
| 6 | 0 | 3.02E-07 | 1.49E-14 | 5.34E-07 | 4.10E-02 of 1.0E-05 (elemental transport He/H partition) | 337 | 92.59 |
| 7 | 0 | 1.40E-08 | 8.20E-15 | 2.28E-07 | 3.77E-02 of 1.0E-05 (elemental transport He/H partition) | 345 | 74.67 |
| 8 | 0 | 5.67E-08 | 9.60E-15 | 1.95E-07 | 3.68E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 73.86 |
| 9 | 0 | 3.93E-08 | 2.02E-14 | 5.32E-07 | 3.68E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 76.46 |
| 10 | 0 | 2.82E-08 | 1.25E-14 | 3.01E-07 | 3.68E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 91.81 |
| 11 | 0 | 5.58E-08 | 1.07E-14 | 2.75E-07 | 3.68E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 94.44 |
| 12 | 0 | 1.60E-07 | 1.32E-14 | 6.41E-07 | 3.68E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 82.28 |
| 13 | 0 | 4.31E-07 | 1.62E-14 | 7.77E-07 | 3.68E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 55.14 |
| 14 | 0 | 2.87E-07 | 1.21E-14 | 4.12E-07 | 3.68E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 54.12 |
| 15 | 0 | 5.24E-07 | 1.21E-14 | 9.25E-07 | 3.68E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 37.72 |
| 16 | 0 | 2.02E-07 | 1.16E-14 | 9.65E-07 | 3.67E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 38.18 |
| 17 | 0 | 1.59E-08 | 1.01E-14 | 2.61E-07 | 3.67E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 38.84 |
| 18 | 0 | 2.91E-08 | 1.13E-14 | 1.95E-07 | 3.66E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 37.83 |
| 19 | 0 | 5.08E-08 | 1.16E-14 | 2.64E-07 | 3.64E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 37.63 |
| 20 | 0 | 3.67E-08 | 1.47E-14 | 3.33E-07 | 3.60E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 16.79 |
| 21 | 0 | 1.35E-07 | 9.06E-15 | 4.87E-07 | 3.52E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 13.07 |
| 22 | 0 | 1.77E-07 | 9.06E-15 | 8.85E-07 | 3.37E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.59 |
| 23 | 0 | 1.27E-07 | 1.02E-14 | 4.74E-07 | 3.11E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 13.06 |
| 24 | 0 | 1.39E-07 | 1.02E-14 | 5.50E-07 | 2.70E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.52 |
| 25 | 0 | 1.90E-08 | 1.02E-14 | 3.94E-07 | 2.13E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.50 |
| 26 | 0 | 9.87E-08 | 1.02E-14 | 4.73E-07 | 1.50E-02 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.59 |
| 27 | 0 | 2.00E-08 | 1.02E-14 | 5.18E-07 | 9.39E-03 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.54 |
| 28 | 0 | 1.90E-08 | 1.02E-14 | 5.43E-07 | 5.38E-03 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.63 |
| 29 | 0 | 1.90E-08 | 1.02E-14 | 5.57E-07 | 2.90E-03 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.55 |
| 30 | 0 | 3.43E-08 | 1.02E-14 | 5.65E-07 | 1.51E-03 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.44 |
| 31 | 0 | 3.43E-08 | 1.02E-14 | 5.69E-07 | 7.71E-04 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.48 |
| 32 | 0 | 3.43E-08 | 1.02E-14 | 5.72E-07 | 3.90E-04 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.60 |
| 33 | 0 | 3.43E-08 | 1.02E-14 | 5.73E-07 | 1.96E-04 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.61 |
| 34 | 0 | 3.43E-08 | 1.02E-14 | 5.74E-07 | 9.82E-05 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.62 |
| 35 | 0 | 3.43E-08 | 1.02E-14 | 5.74E-07 | 4.92E-05 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.64 |
| 36 | 0 | 3.43E-08 | 1.02E-14 | 5.75E-07 | 2.46E-05 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.68 |
| 37 | 0 | 3.43E-08 | 1.02E-14 | 5.75E-07 | 1.23E-05 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.55 |
| 38 | 0 | 3.43E-08 | 1.02E-14 | 5.75E-07 | 6.15E-06 of 1.0E-05 (elemental transport He/H partition) | 500 | 0.40 |

- solver verdict: **info = 0** (last `||R||` = 5.749E-07), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- wall clock **30m19s** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-16 23:22:39 to 2026-09-16 23:52:58

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

## Re-measuring this state, and how closely it comes back

The two files of `output/` can be handed back to the binary as `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with `Restart intent: stationary evaluate`, which measures the state as it stands, writes it back unchanged and exits 0 only if every active equation is within its tolerance.

That is the post-processing pass of this record: it is the route the advection-corrected profiles and the mass-loss line above were produced on. Three states are named in it and every product says which one it describes: the LOADED state, the conserved variables the two files carry, which the pass does not write to; the WORK state, which holds those conserved variables fixed and derives the pressure and the temperature from the composition one equilibrium sweep returns, through the caloric equation of state, and which is what `Hydro_ioniz.txt`, `Ion_species.txt`, the breakdowns and the `certified=` pair of their header describe; and the ADVECTION-DERIVED composition the post-process builds from the work state, which is what `Hydro_ioniz_adv.txt` and `Ion_species_adv.txt` carry. The derived pair states the work state's certification pair on a `# derived_from:` line, as provenance of the state it was built from, and makes no certification claim about its own rows. One sweep is not a closed chemical and thermal fixed point, so the work state reports its own closure defect in `pp.log` and "refreshed" never means "closed".

Two answers are kept apart in `pp.log` and in the list above: whether the stationary claim the solved state was written with reproduces when that state is handed back, and what verdict the work state gets. A passing work state is not a confirmation of a claim that did not reproduce.

What comes back exactly and what does not, MEASURED over the 120 certified states of `models/` (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8). The conserved state round-trips to the last bit: the radius, velocity and pressure columns are read as written, and since item L18 the mass density is read from its own column rather than rebuilt from the species, so the state re-measured is the state certified. The ROW MEASURES do not round-trip bitwise and cannot: the first equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at its rounding floor, where one unit in the last place of the density moves the cell-wise maximum of the energy row by about 11 per cent. Re-evaluated against in-run, the mass row comes back within a factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to 1.93 (median 1.22); the momentum row sits four to six decades below its tolerance and its ratio is the floor itself. A state certified at more than about half of its tolerance may therefore be refused on re-evaluation, and that refusal is a property of the cancellation and not of this file.
