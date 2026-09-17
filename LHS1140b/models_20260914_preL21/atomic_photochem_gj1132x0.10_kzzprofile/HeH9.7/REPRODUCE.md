# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: the Photochem column handed over as `Lower atmosphere profile:` (`photochem`)
- spectrum: `gj1132x0.10`
- mixing: binary H/He element diffusion with K_zz(p) from the profile (`kzzprofile`)

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/../../EXHALE_L14.x` |
| md5 | `dd7ee3188ead1cf54d2f2c7f99246309` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_photochem_gj1132x0.10_kzzprofile_HeH9.7` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.20_kzzprofile/HeH9.7/output  tier3:models/atomic_photochem_gj1132x0.20_kzzprofile  HeH=9.71103  target=9.71103  dlog10=0.0000  candidates=560  (a certified case at another XUV normalization); then a continuation from the first solve's written state at dtau0=1.0e8 (its log run_dtau0_first.log)
```

The fields are the state directory, the tier it was taken from (`tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

The grid the seed is written on is named by `target_grid_Hydro_ioniz.txt`. This case builds its grid on a lower-atmosphere profile, so its radial scale R0 is the radius that profile carries at its matching level times R_J, 1.1592997812185061E+09 cm, and not the planet radius the shared `current_grid_Hydro_ioniz.txt` states. `load_IC` compares the `grid` metadata field as text, so the seed has to carry this number: the file is that shared one with the single field replaced, written by `src/utils/profile_match_level.py`.

It was interpolated onto the cell centers of the current code by

```
# the target grid header, this case's R0 in place of the shared file's
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/profile_match_level.py \
    lower_atmosphere_profile.dat /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt target_grid_Hydro_ioniz.txt
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.20_kzzprofile/HeH9.7/output \
    target_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
# the target grid header, this case's R0 in place of the shared file's
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/profile_match_level.py \
    lower_atmosphere_profile.dat /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt target_grid_Hydro_ioniz.txt
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132x0.20_kzzprofile/HeH9.7/output \
    target_grid_Hydro_ioniz.txt output --ic 

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 \
    EXHALE_OUTER_PASSES=30 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/../../EXHALE_L14.x > run.log 2>&1

# the first solve returned info=1: the written state reloaded and solved
# again with the pseudo-time start raised (L4e)
cp output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
cp output/Ion_species.txt output/Ion_species_IC.txt
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e8 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/../../EXHALE_L14.x > run.log 2>&1

# the advection-corrected profiles: the solved state becomes the IC and one
# time step is taken at a CFL small enough to leave it where it is
for f in Hydro_ioniz Ion_species; do cp -f output/$f.txt output/${f}_IC.txt; done
cp -f input.inp input.inp.solved
sed -i -e 's/^Load IC?.*/Load IC? True/' -e 's/^Do only PP:.*/Do only PP: True/' \
    -e '/^Restart intent:/d' -e '/^Solver:/d' -e '/^du_th /d' input.inp
echo 'CFL: 1.0e-12' >> input.inp
OMP_NUM_THREADS=8 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/../../EXHALE_L14.x > pp.log 2>&1
mv -f input.inp.solved input.inp

# the transit spectrum, through the WINERED HIRES-Y kernel of the measurement
. /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/winered_hires_y.sh
MPLBACKEND=Agg PYTHONPATH=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00 python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_transit.py > transit.log 2>&1
```

## What came out

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 2.35E-08 | 4.51E-14 | 7.59E-07 | 2.24E-03 of 1.0E-05 (elemental transport He/H partition) | 361 | 244.39 |
| 2 | 0 | 2.44E-08 | 1.81E-13 | 7.71E-07 | 2.37E-03 of 1.0E-05 (elemental transport He/H partition) | 364 | 101.48 |
| 3 | 0 | 9.50E-09 | 2.83E-14 | 3.66E-07 | 1.63E-03 of 1.0E-05 (elemental transport He/H partition) | 367 | 58.98 |
| 4 | 0 | 8.27E-09 | 1.99E-14 | 1.94E-07 | 9.80E-04 of 1.0E-05 (elemental transport He/H partition) | 369 | 213.86 |
| 5 | 0 | 1.24E-08 | 8.23E-14 | 4.18E-07 | 5.53E-04 of 1.0E-05 (elemental transport He/H partition) | 371 | 78.09 |
| 6 | 0 | 1.75E-08 | 4.82E-14 | 5.83E-07 | 3.00E-04 of 1.0E-05 (elemental transport He/H partition) | 373 | 26.15 |
| 7 | 0 | 2.28E-08 | 1.98E-13 | 6.17E-07 | 1.59E-04 of 1.0E-05 (elemental transport He/H partition) | 374 | 25.58 |
| 8 | 0 | 1.17E-08 | 1.00E-13 | 4.51E-07 | 8.33E-05 of 1.0E-05 (elemental transport He/H partition) | 375 | 60.32 |
| 9 | 0 | 2.59E-08 | 2.07E-14 | 6.31E-07 | 4.31E-05 of 1.0E-05 (elemental transport He/H partition) | 376 | 25.62 |
| 10 | 0 | 1.32E-08 | 5.08E-14 | 6.03E-07 | 2.21E-05 of 1.0E-05 (elemental transport He/H partition) | 377 | 12.83 |
| 11 | 0 | 9.73E-09 | 1.57E-14 | 2.52E-07 | 1.13E-05 of 1.0E-05 (elemental transport He/H partition) | 378 | 12.61 |
| 12 | 0 | 6.40E-09 | 3.07E-14 | 1.92E-07 | 5.74E-06 of 1.0E-05 (elemental transport He/H partition) | 379 | 12.31 |

- solver verdict: **info = 0** (last `||R||` = 1.921E-07), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **6.85**, from the post-processing pass
- He I 10830 red-pair equivalent width = **0.2823** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **1.0575** %, FWHM = 0.2487 A (three-Gaussian fit)
- wall clock **63m23s** at `OMP_NUM_THREADS=8`, 2026-09-14 13:35:23 to 2026-09-14 14:38:46

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
