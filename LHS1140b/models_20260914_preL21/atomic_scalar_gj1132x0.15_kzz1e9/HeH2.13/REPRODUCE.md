# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_scalar_gj1132x0.15_kzz1e9/HeH2.13`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: prescribed scalar base, T = 226 K, R_0 = 0.157692 R_J, p = 1 microbar, metal-free (`scalar`)
- spectrum: `gj1132x0.15`
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
| run title | `Simulation for LHS1140b_atomic_scalar_gj1132x0.15_kzz1e9_HeH2.13` |

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
| 1 | 0 | 9.10E-09 | 1.74E-14 | 2.94E-07 | 6.90E-03 of 1.0E-05 (elemental transport He/H partition) | 342 | 184.12 |
| 2 | 0 | 1.43E-08 | 1.26E-13 | 5.11E-07 | 1.27E-02 of 1.0E-05 (elemental transport He/H partition) | 325 | 72.83 |
| 3 | 0 | 1.08E-08 | 6.11E-14 | 2.29E-07 | 1.60E-02 of 1.0E-05 (elemental transport He/H partition) | 324 | 107.16 |
| 4 | 0 | 3.12E-08 | 1.77E-13 | 7.59E-07 | 1.75E-02 of 1.0E-05 (elemental transport He/H partition) | 324 | 51.11 |
| 5 | 0 | 1.99E-08 | 3.83E-14 | 5.74E-07 | 1.81E-02 of 1.0E-05 (elemental transport He/H partition) | 324 | 82.36 |
| 6 | 0 | 5.30E-09 | 1.37E-14 | 1.74E-07 | 1.82E-02 of 1.0E-05 (elemental transport He/H partition) | 324 | 127.17 |
| 7 | 0 | 2.81E-08 | 1.98E-13 | 7.47E-07 | 1.80E-02 of 1.0E-05 (elemental transport He/H partition) | 325 | 36.49 |
| 8 | 0 | 8.03E-09 | 4.00E-14 | 1.29E-07 | 1.75E-02 of 1.0E-05 (elemental transport He/H partition) | 326 | 87.90 |
| 9 | 0 | 4.79E-09 | 4.55E-14 | 2.13E-07 | 1.65E-02 of 1.0E-05 (elemental transport He/H partition) | 327 | 35.59 |
| 10 | 0 | 2.06E-08 | 2.08E-13 | 6.79E-07 | 1.48E-02 of 1.0E-05 (elemental transport He/H partition) | 328 | 34.01 |
| 11 | 0 | 1.53E-08 | 7.74E-14 | 5.80E-07 | 1.25E-02 of 1.0E-05 (elemental transport He/H partition) | 331 | 29.45 |
| 12 | 0 | 6.12E-09 | 1.40E-14 | 1.06E-07 | 9.71E-03 of 1.0E-05 (elemental transport He/H partition) | 333 | 58.18 |
| 13 | 0 | 1.46E-08 | 1.52E-14 | 3.23E-07 | 6.89E-03 of 1.0E-05 (elemental transport He/H partition) | 337 | 24.74 |
| 14 | 0 | 1.72E-08 | 3.67E-14 | 4.09E-07 | 4.45E-03 of 1.0E-05 (elemental transport He/H partition) | 340 | 28.25 |
| 15 | 0 | 4.57E-09 | 1.10E-14 | 1.66E-07 | 2.64E-03 of 1.0E-05 (elemental transport He/H partition) | 343 | 48.12 |
| 16 | 0 | 2.26E-08 | 2.03E-13 | 8.46E-07 | 1.47E-03 of 1.0E-05 (elemental transport He/H partition) | 345 | 16.03 |
| 17 | 0 | 4.91E-09 | 1.01E-14 | 2.00E-07 | 7.82E-04 of 1.0E-05 (elemental transport He/H partition) | 346 | 9.95 |
| 18 | 0 | 7.25E-09 | 1.00E-14 | 2.42E-07 | 4.04E-04 of 1.0E-05 (elemental transport He/H partition) | 346 | 10.95 |
| 19 | 0 | 6.90E-09 | 1.50E-14 | 1.75E-07 | 2.06E-04 of 1.0E-05 (elemental transport He/H partition) | 347 | 7.58 |
| 20 | 0 | 6.61E-09 | 1.98E-14 | 1.54E-07 | 1.04E-04 of 1.0E-05 (elemental transport He/H partition) | 347 | 11.09 |
| 21 | 0 | 5.57E-09 | 2.92E-14 | 2.57E-07 | 5.21E-05 of 1.0E-05 (elemental transport He/H partition) | 347 | 3.11 |
| 22 | 0 | 1.23E-08 | 4.08E-14 | 3.80E-07 | 2.61E-05 of 1.0E-05 (elemental transport He/H partition) | 347 | 3.20 |
| 23 | 0 | 3.25E-09 | 5.92E-14 | 1.68E-07 | 1.31E-05 of 1.0E-05 (elemental transport He/H partition) | 347 | 2.99 |
| 24 | 0 | 3.25E-09 | 5.92E-14 | 7.11E-07 | 6.54E-06 of 1.0E-05 (elemental transport He/H partition) | 347 | 0.14 |

- solver verdict: **info = 0** (last `||R||` = 7.112E-07), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.01**, from the post-processing pass
- He I 10830 red-pair equivalent width = **0.0004** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **0.0012** %, FWHM = 0.2860 A (three-Gaussian fit)
- wall clock **73m16s** at `OMP_NUM_THREADS=8`, 2026-09-14 11:53:35 to 2026-09-14 13:06:51

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
