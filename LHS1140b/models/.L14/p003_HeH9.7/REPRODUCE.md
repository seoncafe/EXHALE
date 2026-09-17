# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`.L14/p003_HeH9.7`


The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L14.x` |
| md5 | `a11038c245050d4f11852af829e13fc7` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_L14_p003_photochem_HeH9.7` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L14/p005_HeH9.7/output (given)
```

The fields are the state directory, the tier it was taken from (`tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

The grid the seed is written on is named by `target_grid_Hydro_ioniz.txt`. This case builds its grid on a lower-atmosphere profile, so its radial scale R0 is the radius that profile carries at its matching level times R_J, 1.1592997812185061E+09 cm, and not the planet radius the shared `current_grid_Hydro_ioniz.txt` states. `load_IC` compares the `grid` metadata field as text, so the seed has to carry this number: the file is that shared one with the single field replaced, written by `src/utils/profile_match_level.py`.

It was interpolated onto the cell centers of the current code by

```
# the target grid header, this case's R0 in place of the shared file's
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/profile_match_level.py \
    lower_atmosphere_profile.dat /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt target_grid_Hydro_ioniz.txt
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L14/p005_HeH9.7/output \
    target_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
# the target grid header, this case's R0 in place of the shared file's
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/profile_match_level.py \
    lower_atmosphere_profile.dat /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt target_grid_Hydro_ioniz.txt
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L14/p005_HeH9.7/output \
    target_grid_Hydro_ioniz.txt output --ic 

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e8 \
    EXHALE_OUTER_PASSES=40 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L14.x > run.log 2>&1

# the advection-corrected profiles: the solved state becomes the IC and one
# time step is taken at a CFL small enough to leave it where it is
for f in Hydro_ioniz Ion_species; do cp -f output/$f.txt output/${f}_IC.txt; done
cp -f input.inp input.inp.solved
sed -i -e 's/^Load IC?.*/Load IC? True/' -e 's/^Do only PP:.*/Do only PP: True/' \
    -e '/^Restart intent:/d' -e '/^Solver:/d' -e '/^du_th /d' input.inp
echo 'CFL: 1.0e-12' >> input.inp
OMP_NUM_THREADS=8 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L14.x > pp.log 2>&1
mv -f input.inp.solved input.inp

# the transit spectrum, through the WINERED HIRES-Y kernel of the measurement
. /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/winered_hires_y.sh
MPLBACKEND=Agg PYTHONPATH=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00 python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_transit.py > transit.log 2>&1
```

## What came out

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 4.44E-08 | 1.82E-14 | 9.62E-07 | 3.10E-02 of 1.0E-05 (elemental transport He/H partition) | 359 | 2213.92 |
| 2 | 0 | 3.63E-08 | 4.18E-14 | 9.59E-07 | 4.07E-02 of 1.0E-05 (elemental transport He/H partition) | 310 | 1900.61 |
| 3 | 2 | 8.95E-04 | 1.85E-10 | 4.27E-03 | 4.90E-02 of 1.0E-05 (elemental transport He/H partition) | 310 | 1308.94 |
| 4 | 0 | 3.75E-08 | 1.34E-14 | 8.58E-07 | 5.12E-02 of 1.0E-05 (elemental transport He/H partition) | 310 | 1141.79 |
| 5 | 2 | 4.36E-04 | 7.99E-11 | 1.53E-03 | 4.97E-02 of 1.0E-05 (elemental transport He/H partition) | 310 | 1238.94 |
| 6 | 0 | 4.46E-08 | 1.45E-14 | 8.45E-07 | 4.54E-02 of 1.0E-05 (elemental transport He/H partition) | 309 | 1279.00 |
| 7 | 0 | 3.07E-08 | 2.02E-14 | 9.70E-07 | 3.89E-02 of 1.0E-05 (elemental transport He/H partition) | 308 | 1257.31 |
| 8 | 2 | 3.55E-06 | 9.64E-13 | 1.90E-05 | 3.09E-02 of 1.0E-05 (elemental transport He/H partition) | 306 | 1651.54 |
| 9 | 0 | 3.09E-08 | 1.85E-14 | 7.58E-07 | 2.29E-02 of 1.0E-05 (elemental transport He/H partition) | 305 | 1253.53 |
| 10 | 0 | 3.00E-08 | 2.63E-14 | 6.66E-07 | 1.59E-02 of 1.0E-05 (elemental transport He/H partition) | 303 | 1282.14 |
| 11 | 2 | 7.94E-07 | 1.91E-13 | 4.56E-06 | 1.05E-02 of 1.0E-05 (elemental transport He/H partition) | 301 | 1324.86 |
| 12 | 0 | 1.26E-07 | 2.73E-14 | 8.42E-07 | 6.74E-03 of 1.0E-05 (elemental transport He/H partition) | 300 | 304.03 |
| 13 | 0 | 2.35E-08 | 2.32E-14 | 7.80E-07 | 4.24E-03 of 1.0E-05 (elemental transport He/H partition) | 298 | 450.45 |
| 14 | 0 | 4.25E-08 | 1.93E-14 | 9.44E-07 | 2.63E-03 of 1.0E-05 (elemental transport He/H partition) | 297 | 513.84 |
| 15 | 0 | 1.27E-07 | 1.23E-14 | 8.99E-07 | 1.61E-03 of 1.0E-05 (elemental transport He/H partition) | 296 | 445.51 |
| 16 | 0 | 1.68E-07 | 2.48E-14 | 9.00E-07 | 9.77E-04 of 1.0E-05 (elemental transport He/H partition) | 295 | 345.80 |
| 17 | 0 | 2.99E-08 | 1.62E-14 | 8.14E-07 | 5.89E-04 of 1.0E-05 (elemental transport He/H partition) | 295 | 489.12 |
| 18 | 0 | 3.31E-08 | 1.80E-14 | 8.36E-07 | 3.51E-04 of 1.0E-05 (elemental transport He/H partition) | 295 | 285.54 |
| 19 | 0 | 8.18E-08 | 2.24E-14 | 9.50E-07 | 2.07E-04 of 1.0E-05 (elemental transport He/H partition) | 295 | 764.54 |
| 20 | 0 | 1.57E-07 | 1.27E-14 | 8.35E-07 | 1.21E-04 of 1.0E-05 (elemental transport He/H partition) | 295 | 204.49 |
| 21 | 0 | 2.02E-08 | 1.25E-14 | 6.88E-07 | 7.02E-05 of 1.0E-05 (elemental transport He/H partition) | 295 | 129.89 |
| 22 | 0 | 2.51E-08 | 3.73E-14 | 7.56E-07 | 4.02E-05 of 1.0E-05 (elemental transport He/H partition) | 295 | 238.59 |
| 23 | 0 | 5.75E-08 | 1.78E-14 | 7.76E-07 | 2.28E-05 of 1.0E-05 (elemental transport He/H partition) | 295 | 89.70 |
| 24 | 0 | 2.88E-08 | 1.52E-14 | 9.81E-07 | 1.27E-05 of 1.0E-05 (elemental transport He/H partition) | 295 | 93.43 |
| 25 | 0 | 2.12E-08 | 1.35E-14 | 7.59E-07 | 7.06E-06 of 1.0E-05 (elemental transport He/H partition) | 295 | 0.15 |

- solver verdict: **info = 0** (last `||R||` = 7.589E-07), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **6.30**, from the post-processing pass
- He I 10830 red-pair equivalent width = **0.0004** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **0.0011** %, FWHM = 0.3067 A (three-Gaussian fit)
- wall clock **337m17s** at `OMP_NUM_THREADS=8`, 2026-09-15 01:55:56 to 2026-09-15 07:33:13

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
