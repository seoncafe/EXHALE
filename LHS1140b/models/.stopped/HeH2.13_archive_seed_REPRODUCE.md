# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_scalar_gj1132_kzz1e9/HeH2.13`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: prescribed scalar base, T = 226 K, R_0 = 0.157692 R_J, p = 1 microbar, metal-free (`scalar`)
- spectrum: `gj1132`
- mixing: binary H/He element diffusion with He_Kzz = 1e9 cm^2/s (`kzz1e9`)

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `97e10317a710b9ccc63addbedde3586a` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_scalar_gj1132_kzz1e9_HeH2.13` |

## The seed

`models/pick_seed.py` chose, out of the archived states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/basemetals/a_heh/output  refresh_j96  HeH=2.13  target=2.13  dlog10=0.0000  candidates=50
```

The fields are the state directory, the archive generation, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many archived states carried the physics at all.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/basemetals/a_heh/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/basemetals/a_heh/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > run.log 2>&1
```

## What came out

The outer passes of the marching plus the JFNK hand-off, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 2 | 1.23E+00 | 6.42E-04 | 1.78E+00 | 2.98E-04 of 1.0E-05 (elemental transport He/H partition) | 317 | 547.36 |
| 2 | 2 | 1.69E+00 | 6.41E-04 | 1.78E+00 | 3.16E-04 of 1.0E-05 (elemental transport He/H partition) | 306 | 256.27 |
| 3 | 2 | 1.42E+00 | 6.35E-04 | 1.78E+00 | 2.57E-04 of 1.0E-05 (elemental transport He/H partition) | 303 | 385.42 |
| 4 | 2 | 1.70E+00 | 6.21E-04 | 1.78E+00 | 1.80E-04 of 1.0E-05 (elemental transport He/H partition) | 301 | 357.81 |
| 5 | 2 | 1.54E+00 | 6.17E-04 | 1.79E+00 | 1.16E-04 of 1.0E-05 (elemental transport He/H partition) | 300 | 330.37 |

- solver verdict: **info = 2** (last `||R||` = 1.785E+00), by the marching plus the JFNK hand-off
- certification of the state written: **NOT CERTIFIED**
  - hydrodynamic mass row: row measure  3.777E-01 above  1.1E-09 at cell 12
  - hydrodynamic momentum row: row measure  6.167E-04 above  1.0E-08 at cell 500
  - hydrodynamic energy row: row measure  1.785E+00 above  1.0E-06 at cell 3
  - elemental transport He/H partition: gated row measure  1.162E-04 above  1.0E-05 at cell 300 (a wind cell)
- wall clock **37m32s** at `OMP_NUM_THREADS=8`, 2026-09-13 16:10:59 to 2026-09-13 16:48:31

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is an archived state, not a fresh one: it is read from `archive_20260830/`, which is never rewritten, and interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because the archive changed, is still the same solution.
