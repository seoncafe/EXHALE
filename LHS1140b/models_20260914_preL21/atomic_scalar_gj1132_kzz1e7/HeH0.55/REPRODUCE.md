# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_scalar_gj1132_kzz1e7/HeH0.55`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: prescribed scalar base, T = 226 K, R_0 = 0.157692 R_J, p = 1 microbar, metal-free (`scalar`)
- spectrum: `gj1132`
- mixing: binary H/He element diffusion with He_Kzz = 1e7 cm^2/s (`kzz1e7`)

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
| run title | `Simulation for LHS1140b_atomic_scalar_gj1132_kzz1e7_HeH0.55` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/heh0p55_diff_kzz1e7/output  tier4:archive/original  HeH=0.55  target=0.55  dlog10=0.0000  candidates=34  (an archived state)
```

The fields are the state directory, the tier it was taken from (`tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/heh0p55_diff_kzz1e7/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/heh0p55_diff_kzz1e7/output \
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
| 1 | 0 | 5.17E-10 | 1.27E-14 | 1.75E-08 | 1.58E-03 of 1.0E-05 (elemental transport He/H partition) | 260 | 208.84 |
| 2 | 0 | 1.50E-09 | 2.50E-13 | 4.13E-08 | 4.42E-04 of 1.0E-05 (elemental transport He/H partition) | 217 | 186.35 |
| 3 | 0 | 6.38E-10 | 6.59E-14 | 1.22E-08 | 3.24E-04 of 1.0E-05 (elemental transport He/H partition) | 217 | 168.77 |
| 4 | 0 | 1.06E-09 | 3.60E-12 | 2.38E-08 | 2.53E-04 of 1.0E-05 (elemental transport He/H partition) | 285 | 159.38 |
| 5 | 0 | 1.05E-09 | 4.87E-14 | 1.95E-08 | 2.02E-04 of 1.0E-05 (elemental transport He/H partition) | 287 | 169.00 |
| 6 | 0 | 1.11E-09 | 1.32E-12 | 3.06E-08 | 1.47E-04 of 1.0E-05 (elemental transport He/H partition) | 288 | 154.70 |
| 7 | 0 | 1.37E-09 | 3.92E-13 | 3.70E-08 | 9.86E-05 of 1.0E-05 (elemental transport He/H partition) | 289 | 153.64 |
| 8 | 0 | 8.86E-10 | 9.15E-14 | 1.96E-08 | 6.26E-05 of 1.0E-05 (elemental transport He/H partition) | 291 | 133.85 |
| 9 | 0 | 1.32E-09 | 2.73E-12 | 2.31E-08 | 3.81E-05 of 1.0E-05 (elemental transport He/H partition) | 292 | 124.26 |
| 10 | 0 | 1.04E-09 | 1.02E-12 | 2.31E-08 | 2.23E-05 of 1.0E-05 (elemental transport He/H partition) | 293 | 123.87 |
| 11 | 0 | 7.57E-10 | 2.92E-13 | 2.21E-08 | 1.27E-05 of 1.0E-05 (elemental transport He/H partition) | 293 | 118.33 |
| 12 | 0 | 8.37E-10 | 2.43E-12 | 2.11E-08 | 7.09E-06 of 1.0E-05 (elemental transport He/H partition) | 294 | 93.84 |

- solver verdict: **info = 0** (last `||R||` = 2.112E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.85**, from the post-processing pass
- He I 10830 red-pair equivalent width = **0.1682** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **0.6909** %, FWHM = 0.2290 A (three-Gaussian fit)
- wall clock **30m33s** at `OMP_NUM_THREADS=8`, 2026-09-14 03:30:41 to 2026-09-14 04:01:14

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
