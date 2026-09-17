# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: prescribed scalar base, T = 226 K, R_0 = 0.157692 R_J, p = 1 microbar, metal-free (`scalar`)
- spectrum: `gj1132x0.10`
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
| run title | `Simulation for LHS1140b_atomic_scalar_gj1132x0.10_kzz1e9_HeH2.13` |

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
| 1 | 0 | 1.72E-08 | 1.59E-14 | 4.67E-07 | 1.28E-02 of 1.0E-05 (elemental transport He/H partition) | 348 | 219.03 |
| 2 | 0 | 1.40E-08 | 1.20E-14 | 5.43E-07 | 1.89E-02 of 1.0E-05 (elemental transport He/H partition) | 333 | 139.08 |
| 3 | 0 | 2.37E-08 | 5.43E-14 | 8.94E-07 | 2.21E-02 of 1.0E-05 (elemental transport He/H partition) | 332 | 242.19 |
| 4 | 0 | 2.59E-08 | 2.49E-14 | 9.47E-07 | 2.34E-02 of 1.0E-05 (elemental transport He/H partition) | 332 | 111.06 |
| 5 | 0 | 1.49E-08 | 1.58E-14 | 2.35E-07 | 2.40E-02 of 1.0E-05 (elemental transport He/H partition) | 332 | 114.58 |
| 6 | 0 | 1.47E-08 | 3.65E-14 | 3.32E-07 | 2.43E-02 of 1.0E-05 (elemental transport He/H partition) | 332 | 148.27 |
| 7 | 0 | 1.09E-08 | 9.01E-14 | 5.66E-07 | 2.45E-02 of 1.0E-05 (elemental transport He/H partition) | 332 | 109.24 |
| 8 | 0 | 2.27E-08 | 3.16E-14 | 4.91E-07 | 2.45E-02 of 1.0E-05 (elemental transport He/H partition) | 332 | 44.17 |
| 9 | 0 | 4.46E-09 | 3.76E-14 | 2.02E-07 | 2.44E-02 of 1.0E-05 (elemental transport He/H partition) | 332 | 80.39 |
| 10 | 0 | 7.47E-09 | 2.93E-14 | 3.51E-07 | 2.43E-02 of 1.0E-05 (elemental transport He/H partition) | 333 | 128.15 |
| 11 | 0 | 1.29E-08 | 1.63E-14 | 2.85E-07 | 2.40E-02 of 1.0E-05 (elemental transport He/H partition) | 333 | 83.34 |
| 12 | 0 | 8.95E-09 | 3.05E-14 | 3.99E-07 | 2.35E-02 of 1.0E-05 (elemental transport He/H partition) | 334 | 45.99 |
| 13 | 0 | 1.36E-08 | 5.60E-14 | 3.76E-07 | 2.26E-02 of 1.0E-05 (elemental transport He/H partition) | 335 | 45.52 |
| 14 | 0 | 1.54E-08 | 1.72E-13 | 7.83E-07 | 2.11E-02 of 1.0E-05 (elemental transport He/H partition) | 337 | 23.03 |
| 15 | 0 | 2.03E-08 | 3.10E-14 | 3.14E-07 | 1.89E-02 of 1.0E-05 (elemental transport He/H partition) | 339 | 51.68 |
| 16 | 0 | 1.90E-08 | 1.27E-13 | 8.44E-07 | 1.61E-02 of 1.0E-05 (elemental transport He/H partition) | 343 | 15.57 |
| 17 | 0 | 1.41E-08 | 2.51E-14 | 5.46E-07 | 1.30E-02 of 1.0E-05 (elemental transport He/H partition) | 349 | 7.72 |
| 18 | 0 | 1.72E-08 | 1.43E-13 | 5.56E-07 | 9.87E-03 of 1.0E-05 (elemental transport He/H partition) | 354 | 14.56 |
| 19 | 0 | 7.32E-09 | 9.77E-15 | 2.77E-07 | 6.98E-03 of 1.0E-05 (elemental transport He/H partition) | 359 | 7.61 |
| 20 | 0 | 8.20E-09 | 7.62E-14 | 2.90E-07 | 4.54E-03 of 1.0E-05 (elemental transport He/H partition) | 362 | 2.26 |
| 21 | 0 | 8.72E-09 | 5.18E-14 | 3.61E-07 | 2.71E-03 of 1.0E-05 (elemental transport He/H partition) | 364 | 2.87 |
| 22 | 0 | 1.16E-08 | 9.62E-15 | 2.18E-07 | 1.51E-03 of 1.0E-05 (elemental transport He/H partition) | 366 | 2.70 |
| 23 | 0 | 1.16E-08 | 9.62E-15 | 6.53E-07 | 8.06E-04 of 1.0E-05 (elemental transport He/H partition) | 367 | 0.20 |
| 24 | 0 | 1.07E-08 | 2.18E-14 | 4.12E-07 | 4.17E-04 of 1.0E-05 (elemental transport He/H partition) | 367 | 2.16 |
| 25 | 0 | 1.12E-08 | 2.50E-14 | 3.60E-07 | 2.12E-04 of 1.0E-05 (elemental transport He/H partition) | 367 | 0.20 |
| 26 | 0 | 1.10E-08 | 2.44E-14 | 3.59E-07 | 1.07E-04 of 1.0E-05 (elemental transport He/H partition) | 368 | 0.20 |
| 27 | 0 | 1.10E-08 | 2.44E-14 | 4.06E-07 | 5.37E-05 of 1.0E-05 (elemental transport He/H partition) | 368 | 0.20 |
| 28 | 0 | 1.07E-08 | 2.18E-14 | 4.47E-07 | 2.69E-05 of 1.0E-05 (elemental transport He/H partition) | 368 | 0.20 |
| 29 | 0 | 1.07E-08 | 2.18E-14 | 4.71E-07 | 1.35E-05 of 1.0E-05 (elemental transport He/H partition) | 368 | 0.20 |
| 30 | 0 | 1.07E-08 | 2.18E-14 | 4.84E-07 | 6.74E-06 of 1.0E-05 (elemental transport He/H partition) | 368 | 0.15 |

- solver verdict: **info = 0** (last `||R||` = 4.842E-07), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **6.83**, from the post-processing pass
- He I 10830 red-pair equivalent width = **0.0001** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **0.0003** %, FWHM = 0.2805 A (three-Gaussian fit)
- wall clock **89m00s** at `OMP_NUM_THREADS=8`, 2026-09-14 11:53:35 to 2026-09-14 13:22:35

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
