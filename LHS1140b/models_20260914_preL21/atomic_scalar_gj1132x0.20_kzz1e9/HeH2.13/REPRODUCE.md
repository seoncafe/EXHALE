# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: prescribed scalar base, T = 226 K, R_0 = 0.157692 R_J, p = 1 microbar, metal-free (`scalar`)
- spectrum: `gj1132x0.20`
- mixing: binary H/He element diffusion with He_Kzz = 1e9 cm^2/s (`kzz1e9`)

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L14.x` |
| md5 | `dd7ee3188ead1cf54d2f2c7f99246309` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_scalar_gj1132x0.20_kzz1e9_HeH2.13` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.25_kzz1e9/HeH2.13/output  tier3:models/atomic_scalar_gj1132x0.25_kzz1e9  HeH=2.13  target=2.13  dlog10=0.0000  candidates=100  (a certified case at another XUV normalization); then a continuation from the first solve's written state at dtau0=1.0e8 (its log run_dtau0_first.log)
```

The fields are the state directory, the tier it was taken from (`tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.25_kzz1e9/HeH2.13/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132x0.25_kzz1e9/HeH2.13/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 \
    EXHALE_OUTER_PASSES=30 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L14.x > run.log 2>&1

# the first solve returned info=1: the written state reloaded and solved
# again with the pseudo-time start raised (L4e)
cp output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
cp output/Ion_species.txt output/Ion_species_IC.txt
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e8 \
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
| 1 | 0 | 2.57E-08 | 2.97E-13 | 6.19E-07 | 2.86E-03 of 1.0E-05 (elemental transport He/H partition) | 338 | 57.64 |
| 2 | 0 | 7.05E-09 | 3.20E-14 | 1.49E-07 | 6.10E-03 of 1.0E-05 (elemental transport He/H partition) | 321 | 150.92 |
| 3 | 0 | 8.04E-09 | 3.16E-14 | 2.08E-07 | 8.26E-03 of 1.0E-05 (elemental transport He/H partition) | 320 | 93.98 |
| 4 | 0 | 8.05E-09 | 3.75E-14 | 2.46E-07 | 9.06E-03 of 1.0E-05 (elemental transport He/H partition) | 320 | 129.30 |
| 5 | 0 | 1.53E-08 | 7.57E-14 | 3.25E-07 | 8.83E-03 of 1.0E-05 (elemental transport He/H partition) | 320 | 62.01 |
| 6 | 0 | 9.76E-09 | 2.70E-13 | 4.72E-07 | 7.81E-03 of 1.0E-05 (elemental transport He/H partition) | 320 | 37.92 |
| 7 | 0 | 1.49E-08 | 2.04E-13 | 5.47E-07 | 6.25E-03 of 1.0E-05 (elemental transport He/H partition) | 320 | 44.77 |
| 8 | 0 | 5.85E-09 | 3.77E-14 | 1.77E-07 | 4.52E-03 of 1.0E-05 (elemental transport He/H partition) | 320 | 43.73 |
| 9 | 0 | 6.11E-09 | 7.15E-14 | 3.05E-07 | 2.98E-03 of 1.0E-05 (elemental transport He/H partition) | 320 | 75.24 |
| 10 | 0 | 7.35E-09 | 7.41E-14 | 2.36E-07 | 1.81E-03 of 1.0E-05 (elemental transport He/H partition) | 320 | 50.92 |
| 11 | 0 | 2.13E-09 | 3.58E-14 | 1.14E-07 | 1.05E-03 of 1.0E-05 (elemental transport He/H partition) | 318 | 30.79 |
| 12 | 0 | 1.20E-08 | 8.90E-14 | 3.23E-07 | 5.91E-04 of 1.0E-05 (elemental transport He/H partition) | 316 | 23.69 |
| 13 | 0 | 3.10E-09 | 1.29E-14 | 1.12E-07 | 3.28E-04 of 1.0E-05 (elemental transport He/H partition) | 314 | 16.71 |
| 14 | 0 | 1.51E-08 | 1.63E-13 | 5.26E-07 | 1.81E-04 of 1.0E-05 (elemental transport He/H partition) | 311 | 10.70 |
| 15 | 0 | 4.03E-09 | 2.84E-14 | 8.73E-08 | 9.99E-05 of 1.0E-05 (elemental transport He/H partition) | 309 | 10.57 |
| 16 | 0 | 1.48E-08 | 8.91E-14 | 5.59E-07 | 5.44E-05 of 1.0E-05 (elemental transport He/H partition) | 308 | 15.17 |
| 17 | 0 | 3.96E-09 | 5.15E-14 | 1.08E-07 | 2.93E-05 of 1.0E-05 (elemental transport He/H partition) | 308 | 21.47 |
| 18 | 0 | 4.38E-09 | 2.08E-14 | 1.87E-07 | 1.54E-05 of 1.0E-05 (elemental transport He/H partition) | 308 | 11.13 |
| 19 | 0 | 5.54E-09 | 3.02E-14 | 1.23E-07 | 7.99E-06 of 1.0E-05 (elemental transport He/H partition) | 309 | 11.10 |

- solver verdict: **info = 0** (last `||R||` = 1.228E-07), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.14**, from the post-processing pass
- He I 10830 red-pair equivalent width = **0.0026** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **0.0089** %, FWHM = 0.2765 A (three-Gaussian fit)
- wall clock **67m38s** at `OMP_NUM_THREADS=8`, 2026-09-14 11:53:35 to 2026-09-14 13:01:13

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
