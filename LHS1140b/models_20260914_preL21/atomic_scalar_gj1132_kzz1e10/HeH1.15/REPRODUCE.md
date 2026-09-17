# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_scalar_gj1132_kzz1e10/HeH1.15`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: prescribed scalar base, T = 226 K, R_0 = 0.157692 R_J, p = 1 microbar, metal-free (`scalar`)
- spectrum: `gj1132`
- mixing: binary H/He element diffusion with He_Kzz = 1e10 cm^2/s (`kzz1e10`)

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `b261c3287e64149a85d0d2dfd2ee2c74` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_scalar_gj1132_kzz1e10_HeH1.15` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/misc/kzz1e10/output  tier4:archive/refresh_j96  HeH=0.55  target=1.15  dlog10=0.3203  candidates=20  (an archived state)
```

The fields are the state directory, the tier it was taken from (`tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/misc/kzz1e10/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/misc/kzz1e10/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > run.log 2>&1

# the advection-corrected profiles: the solved state becomes the IC and one
# time step is taken at a CFL small enough to leave it where it is
for f in Hydro_ioniz Ion_species; do cp -f output/$f.txt output/${f}_IC.txt; done
cp -f input.inp input.inp.solved
sed -i -e 's/^Load IC?.*/Load IC? True/' -e 's/^Do only PP:.*/Do only PP: True/' \
    -e '/^Restart intent:/d' -e '/^Solver:/d' -e '/^du_th /d' input.inp
echo 'CFL: 1.0e-12' >> input.inp
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
| 1 | 0 | 1.73E-09 | 2.47E-14 | 3.96E-08 | 3.65E-04 of 1.0E-05 (elemental transport He/H partition) | 302 | 411.40 |
| 2 | 0 | 1.47E-09 | 4.17E-14 | 2.59E-08 | 1.22E-04 of 1.0E-05 (elemental transport He/H partition) | 217 | 185.16 |
| 3 | 0 | 5.63E-10 | 2.52E-14 | 2.17E-08 | 1.01E-04 of 1.0E-05 (elemental transport He/H partition) | 292 | 162.60 |
| 4 | 0 | 9.26E-10 | 1.89E-12 | 2.85E-08 | 1.04E-04 of 1.0E-05 (elemental transport He/H partition) | 295 | 151.85 |
| 5 | 0 | 1.26E-09 | 5.36E-13 | 2.44E-08 | 9.21E-05 of 1.0E-05 (elemental transport He/H partition) | 296 | 151.43 |
| 6 | 0 | 1.10E-09 | 6.06E-13 | 2.47E-08 | 8.00E-05 of 1.0E-05 (elemental transport He/H partition) | 296 | 157.05 |
| 7 | 0 | 1.32E-09 | 4.70E-13 | 2.62E-08 | 6.85E-05 of 1.0E-05 (elemental transport He/H partition) | 297 | 129.13 |
| 8 | 0 | 1.00E-09 | 2.33E-13 | 1.98E-08 | 5.79E-05 of 1.0E-05 (elemental transport He/H partition) | 297 | 130.00 |
| 9 | 0 | 7.40E-10 | 8.40E-14 | 3.06E-08 | 4.84E-05 of 1.0E-05 (elemental transport He/H partition) | 298 | 129.44 |
| 10 | 0 | 8.40E-10 | 7.61E-14 | 1.94E-08 | 4.01E-05 of 1.0E-05 (elemental transport He/H partition) | 298 | 130.28 |
| 11 | 0 | 9.01E-10 | 1.23E-13 | 1.84E-08 | 3.30E-05 of 1.0E-05 (elemental transport He/H partition) | 299 | 130.46 |
| 12 | 0 | 1.09E-09 | 6.96E-14 | 3.33E-08 | 2.70E-05 of 1.0E-05 (elemental transport He/H partition) | 299 | 128.60 |
| 13 | 0 | 7.73E-10 | 1.94E-14 | 2.46E-08 | 2.20E-05 of 1.0E-05 (elemental transport He/H partition) | 300 | 128.55 |
| 14 | 0 | 1.04E-09 | 4.07E-14 | 2.58E-08 | 1.78E-05 of 1.0E-05 (elemental transport He/H partition) | 300 | 130.51 |
| 15 | 0 | 8.85E-10 | 2.20E-12 | 2.97E-08 | 1.43E-05 of 1.0E-05 (elemental transport He/H partition) | 300 | 122.92 |
| 16 | 0 | 8.16E-10 | 1.42E-12 | 2.13E-08 | 1.15E-05 of 1.0E-05 (elemental transport He/H partition) | 301 | 120.71 |
| 17 | 0 | 1.02E-09 | 1.26E-12 | 2.68E-08 | 9.21E-06 of 1.0E-05 (elemental transport He/H partition) | 301 | 124.34 |

- solver verdict: **info = 0** (last `||R||` = 2.682E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.87**, from the post-processing pass
- He I 10830 red-pair equivalent width = **1.1135** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **4.1402** %, FWHM = 0.2527 A (three-Gaussian fit)
- wall clock **44m14s** at `OMP_NUM_THREADS=8`, 2026-09-14 02:36:13 to 2026-09-14 03:20:27

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
