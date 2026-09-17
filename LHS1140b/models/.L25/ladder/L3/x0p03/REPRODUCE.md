# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`.L25/ladder/L3/x0p03`

- numerical flux: `ROE`

The HLLC flux does not reach the stationary root on this grid at this XUV level (L25 step 3), and the same case solved with HLLC on a grid of twice the cells reaches the state the Roe flux reaches on the catalog grid, so the difference is resolution and not a different wind. The Roe flux systematic against HLLC where both solve is Mdot and the He 10830 equivalent width lower by 0.3 to 0.6 per cent, the first cell colder by 6 to 7 per cent and denser by 6 to 8 per cent, and the mass-flux spread unchanged (L25 memo). That systematic is carried wherever the numbers of this case are quoted.

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs two of its keys changed.

## Why this run was made the way it was

L25 step 4, continuation in XUV: rung x0p03 of ladder L3, the catalog input of the target case with only the spectrum changed, solved from the certified state of the rung above it (/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L25/ladder/L3/x0p05/output). Host lart4, 8 threads.

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
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L25/ladder/L3/x0p05/output (given)
```

The fields are the state directory, the tier it was taken from (`tier0` this same case's own most recent certified state -- the one this tree carries, else the one a preserved tree such as `models_20260914_preL21/` holds, which is what a re-solve of the whole catalog continues from -- `tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

The grid the seed is written on is named by `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt`.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L25/ladder/L3/x0p05/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L25/ladder/L3/x0p05/output \
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
| 1 | 0 | 2.40E-08 | 9.08E-15 | 8.42E-08 | 3.11E-02 of 1.0E-05 (elemental transport He/H partition) | 362 | 107.06 |
| 2 | 0 | 1.01E-07 | 1.30E-14 | 4.19E-07 | 3.70E-02 of 1.0E-05 (elemental transport He/H partition) | 308 | 73.51 |
| 3 | 0 | 1.87E-07 | 2.56E-14 | 5.69E-07 | 4.51E-02 of 1.0E-05 (elemental transport He/H partition) | 308 | 66.15 |
| 4 | 0 | 2.32E-08 | 2.73E-14 | 1.83E-07 | 4.72E-02 of 1.0E-05 (elemental transport He/H partition) | 308 | 51.55 |
| 5 | 0 | 1.56E-08 | 1.19E-14 | 1.65E-07 | 4.56E-02 of 1.0E-05 (elemental transport He/H partition) | 308 | 49.36 |
| 6 | 0 | 1.01E-07 | 1.88E-14 | 2.15E-07 | 4.13E-02 of 1.0E-05 (elemental transport He/H partition) | 307 | 42.73 |
| 7 | 0 | 3.80E-08 | 7.77E-15 | 1.01E-07 | 3.50E-02 of 1.0E-05 (elemental transport He/H partition) | 306 | 57.08 |
| 8 | 0 | 2.44E-08 | 1.48E-14 | 1.13E-07 | 2.76E-02 of 1.0E-05 (elemental transport He/H partition) | 304 | 56.65 |
| 9 | 0 | 1.17E-07 | 1.19E-14 | 3.79E-07 | 2.04E-02 of 1.0E-05 (elemental transport He/H partition) | 303 | 56.36 |
| 10 | 0 | 1.13E-07 | 4.06E-14 | 6.43E-07 | 1.42E-02 of 1.0E-05 (elemental transport He/H partition) | 301 | 40.36 |
| 11 | 0 | 8.23E-08 | 2.64E-14 | 4.83E-07 | 9.48E-03 of 1.0E-05 (elemental transport He/H partition) | 300 | 40.46 |
| 12 | 0 | 4.49E-08 | 1.01E-14 | 1.40E-07 | 6.11E-03 of 1.0E-05 (elemental transport He/H partition) | 298 | 34.24 |
| 13 | 0 | 1.97E-07 | 9.27E-15 | 6.47E-07 | 3.84E-03 of 1.0E-05 (elemental transport He/H partition) | 297 | 33.47 |
| 14 | 0 | 1.31E-07 | 2.68E-14 | 5.77E-07 | 2.37E-03 of 1.0E-05 (elemental transport He/H partition) | 297 | 26.38 |
| 15 | 0 | 3.09E-07 | 1.19E-14 | 5.66E-07 | 1.44E-03 of 1.0E-05 (elemental transport He/H partition) | 296 | 25.13 |
| 16 | 0 | 2.30E-07 | 9.58E-15 | 6.43E-07 | 8.66E-04 of 1.0E-05 (elemental transport He/H partition) | 296 | 22.31 |
| 17 | 0 | 1.45E-07 | 1.73E-14 | 4.40E-07 | 5.13E-04 of 1.0E-05 (elemental transport He/H partition) | 296 | 27.24 |
| 18 | 0 | 1.08E-07 | 1.07E-14 | 2.98E-07 | 3.00E-04 of 1.0E-05 (elemental transport He/H partition) | 296 | 21.34 |
| 19 | 0 | 1.34E-07 | 4.36E-14 | 6.88E-07 | 1.74E-04 of 1.0E-05 (elemental transport He/H partition) | 296 | 17.59 |
| 20 | 0 | 6.58E-08 | 9.44E-15 | 6.04E-08 | 9.91E-05 of 1.0E-05 (elemental transport He/H partition) | 296 | 6.96 |
| 21 | 0 | 1.72E-07 | 4.93E-14 | 5.37E-07 | 5.59E-05 of 1.0E-05 (elemental transport He/H partition) | 296 | 3.72 |
| 22 | 0 | 8.63E-08 | 9.97E-15 | 3.90E-07 | 3.11E-05 of 1.0E-05 (elemental transport He/H partition) | 296 | 5.53 |
| 23 | 0 | 6.72E-08 | 1.73E-14 | 4.99E-07 | 1.70E-05 of 1.0E-05 (elemental transport He/H partition) | 296 | 5.35 |
| 24 | 0 | 6.72E-08 | 1.73E-14 | 7.38E-07 | 9.25E-06 of 1.0E-05 (elemental transport He/H partition) | 296 | 0.24 |

- solver verdict: **info = 0** (last `||R||` = 7.378E-07), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- wall clock **14m33s** at `OMP_NUM_THREADS=8` on `lart4`, 2026-09-17 04:56:20 to 2026-09-17 05:10:53

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

## Re-measuring this state, and how closely it comes back

The two files of `output/` can be handed back to the binary as `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with `Restart intent: stationary evaluate`, which measures the state as it stands, writes it back unchanged and exits 0 only if every active equation is within its tolerance.

That is the post-processing pass of this record: it is the route the advection-corrected profiles and the mass-loss line above were produced on. Three states are named in it and every product says which one it describes: the LOADED state, the conserved variables the two files carry, which the pass does not write to; the WORK state, which holds those conserved variables fixed and derives the pressure and the temperature from the composition one equilibrium sweep returns, through the caloric equation of state, and which is what `Hydro_ioniz.txt`, `Ion_species.txt`, the breakdowns and the `certified=` pair of their header describe; and the ADVECTION-DERIVED composition the post-process builds from the work state, which is what `Hydro_ioniz_adv.txt` and `Ion_species_adv.txt` carry. The derived pair states the work state's certification pair on a `# derived_from:` line, as provenance of the state it was built from, and makes no certification claim about its own rows. One sweep is not a closed chemical and thermal fixed point, so the work state reports its own closure defect in `pp.log` and "refreshed" never means "closed".

Two answers are kept apart in `pp.log` and in the list above: whether the stationary claim the solved state was written with reproduces when that state is handed back, and what verdict the work state gets. A passing work state is not a confirmation of a claim that did not reproduce.

What comes back exactly and what does not, MEASURED over the 120 certified states of `models/` (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8). The conserved state round-trips to the last bit: the radius, velocity and pressure columns are read as written, and since item L18 the mass density is read from its own column rather than rebuilt from the species, so the state re-measured is the state certified. The ROW MEASURES do not round-trip bitwise and cannot: the first equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at its rounding floor, where one unit in the last place of the density moves the cell-wise maximum of the energy row by about 11 per cent. Re-evaluated against in-run, the mass row comes back within a factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to 1.93 (median 1.22); the momentum row sits four to six decades below its tolerance and its ratio is the floor itself. A state certified at more than about half of its tolerance may therefore be refused on re-evaluation, and that refusal is a property of the cancellation and not of this file.
