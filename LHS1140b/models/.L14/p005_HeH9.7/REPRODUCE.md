# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`.L14/p005_HeH9.7`


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
| run title | `Simulation for LHS1140b_L14_p005_photochem_HeH9.7` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7/output (given)
```

The fields are the state directory, the tier it was taken from (`tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

The grid the seed is written on is named by `target_grid_Hydro_ioniz.txt`. This case builds its grid on a lower-atmosphere profile, so its radial scale R0 is the radius that profile carries at its matching level times R_J, 1.1592997812185061E+09 cm, and not the planet radius the shared `current_grid_Hydro_ioniz.txt` states. `load_IC` compares the `grid` metadata field as text, so the seed has to carry this number: the file is that shared one with the single field replaced, written by `src/utils/profile_match_level.py`.

It was interpolated onto the cell centers of the current code by

```
# the target grid header, this case's R0 in place of the shared file's
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/profile_match_level.py \
    lower_atmosphere_profile.dat /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt target_grid_Hydro_ioniz.txt
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7/output \
    target_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
# the target grid header, this case's R0 in place of the shared file's
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/profile_match_level.py \
    lower_atmosphere_profile.dat /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt target_grid_Hydro_ioniz.txt
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7/output \
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
| 1 | 0 | 3.97E-08 | 2.43E-14 | 9.82E-07 | 8.85E-03 of 1.0E-05 (elemental transport He/H partition) | 375 | 654.24 |
| 2 | 0 | 2.79E-08 | 3.63E-14 | 6.76E-07 | 1.86E-02 of 1.0E-05 (elemental transport He/H partition) | 360 | 524.65 |
| 3 | 0 | 2.02E-08 | 4.89E-14 | 3.77E-07 | 2.26E-02 of 1.0E-05 (elemental transport He/H partition) | 358 | 247.68 |
| 4 | 0 | 1.93E-08 | 3.49E-14 | 7.15E-07 | 2.42E-02 of 1.0E-05 (elemental transport He/H partition) | 358 | 352.57 |
| 5 | 0 | 2.29E-08 | 1.11E-13 | 6.97E-07 | 2.50E-02 of 1.0E-05 (elemental transport He/H partition) | 358 | 413.78 |
| 6 | 0 | 2.33E-08 | 1.94E-14 | 6.85E-07 | 2.53E-02 of 1.0E-05 (elemental transport He/H partition) | 358 | 237.40 |
| 7 | 0 | 3.40E-08 | 1.99E-14 | 7.11E-07 | 2.55E-02 of 1.0E-05 (elemental transport He/H partition) | 358 | 216.88 |
| 8 | 0 | 2.17E-08 | 3.74E-14 | 5.94E-07 | 2.55E-02 of 1.0E-05 (elemental transport He/H partition) | 358 | 176.00 |
| 9 | 0 | 1.42E-08 | 1.15E-14 | 5.11E-07 | 2.53E-02 of 1.0E-05 (elemental transport He/H partition) | 358 | 172.93 |
| 10 | 0 | 3.63E-08 | 1.47E-14 | 6.64E-07 | 2.49E-02 of 1.0E-05 (elemental transport He/H partition) | 359 | 158.55 |
| 11 | 0 | 1.90E-08 | 2.17E-14 | 4.26E-07 | 2.41E-02 of 1.0E-05 (elemental transport He/H partition) | 361 | 74.30 |
| 12 | 0 | 1.83E-08 | 5.62E-14 | 5.56E-07 | 2.30E-02 of 1.0E-05 (elemental transport He/H partition) | 362 | 150.73 |
| 13 | 0 | 3.63E-08 | 6.34E-14 | 7.93E-07 | 2.13E-02 of 1.0E-05 (elemental transport He/H partition) | 365 | 48.42 |
| 14 | 0 | 1.41E-08 | 1.06E-14 | 3.18E-07 | 1.89E-02 of 1.0E-05 (elemental transport He/H partition) | 368 | 73.75 |
| 15 | 0 | 2.01E-08 | 3.51E-14 | 8.45E-07 | 1.60E-02 of 1.0E-05 (elemental transport He/H partition) | 371 | 63.29 |
| 16 | 2 | 4.76E-08 | 3.80E-14 | 1.03E-06 | 1.27E-02 of 1.0E-05 (elemental transport He/H partition) | 375 | 47.11 |
| 17 | 0 | 4.40E-08 | 1.11E-13 | 9.79E-07 | 9.41E-03 of 1.0E-05 (elemental transport He/H partition) | 378 | 25.15 |
| 18 | 0 | 2.16E-08 | 1.01E-14 | 6.53E-07 | 6.39E-03 of 1.0E-05 (elemental transport He/H partition) | 381 | 41.42 |
| 19 | 0 | 1.65E-08 | 1.94E-14 | 6.96E-07 | 3.97E-03 of 1.0E-05 (elemental transport He/H partition) | 384 | 35.87 |
| 20 | 0 | 2.10E-08 | 9.33E-14 | 8.72E-07 | 2.28E-03 of 1.0E-05 (elemental transport He/H partition) | 386 | 51.09 |
| 21 | 0 | 3.34E-08 | 4.58E-14 | 9.25E-07 | 1.24E-03 of 1.0E-05 (elemental transport He/H partition) | 387 | 11.92 |
| 22 | 0 | 2.69E-08 | 2.53E-14 | 5.78E-07 | 6.49E-04 of 1.0E-05 (elemental transport He/H partition) | 387 | 11.13 |
| 23 | 0 | 1.80E-08 | 1.15E-14 | 5.68E-07 | 3.32E-04 of 1.0E-05 (elemental transport He/H partition) | 388 | 29.21 |
| 24 | 0 | 2.49E-08 | 9.94E-15 | 4.79E-07 | 1.68E-04 of 1.0E-05 (elemental transport He/H partition) | 388 | 3.54 |
| 25 | 0 | 2.26E-08 | 2.65E-14 | 7.83E-07 | 8.47E-05 of 1.0E-05 (elemental transport He/H partition) | 388 | 6.02 |
| 26 | 0 | 1.74E-08 | 1.32E-14 | 4.73E-07 | 4.25E-05 of 1.0E-05 (elemental transport He/H partition) | 388 | 16.02 |
| 27 | 0 | 1.74E-08 | 1.32E-14 | 8.76E-07 | 2.13E-05 of 1.0E-05 (elemental transport He/H partition) | 388 | 0.40 |
| 28 | 0 | 2.06E-08 | 2.63E-14 | 6.37E-07 | 1.06E-05 of 1.0E-05 (elemental transport He/H partition) | 388 | 2.46 |
| 29 | 0 | 2.06E-08 | 2.87E-14 | 6.37E-07 | 5.32E-06 of 1.0E-05 (elemental transport He/H partition) | 388 | 0.14 |

- solver verdict: **info = 0** (last `||R||` = 6.371E-07), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **6.53**, from the post-processing pass
- He I 10830 red-pair equivalent width = **0.0025** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **0.0073** %, FWHM = 0.3183 A (three-Gaussian fit)
- wall clock **64m37s** at `OMP_NUM_THREADS=8`, 2026-09-14 21:37:04 to 2026-09-14 22:41:41

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
