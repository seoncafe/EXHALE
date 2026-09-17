# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`.L7e/fixed5/HeH2.13`


The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L7e.x` |
| md5 | `e4ebc751bcf16a8fa0d3e4def9d9d799` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_molecular_scalar_gj1132_kzz1e9_HeH2.13` |

## The seed

The binary built this case's state out of a certified ATOMIC solution, which is how a molecular case is seeded (`EXHALE_MOLECULAR_SEED`; `models/pick_seed.py` carries no molecular state to choose from):

```
molecular seed from /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L7d/src_atomic_HeH2.13/output, x2 local; then a continuation from the first solve's written state at dtau0=1.0e8 (its log run_dtau0_first.log)
```

The fields are the atomic state directory and which extension of the base H2 partition the conversion carried: `local` gives every cell its own chemical-equilibrium q_H2(p,T), `handoff` carries the base x2 to the outer boundary, a number is that fraction everywhere.

The conversion, which writes `output/*_IC.txt` and stops:

```
# with "Load IC? True" in input.inp, which this case states
OMP_NUM_THREADS=8 EXHALE_MOLECULAR_SEED=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L7d/src_atomic_HeH2.13/output \
    EXHALE_MOLECULAR_SEED_X2=local \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L7e.x > seed.log 2>&1
```

## The commands

```bash
# the molecular seed, built by the binary from the certified atomic state
# with "Load IC? True" in input.inp, which this case states
OMP_NUM_THREADS=8 EXHALE_MOLECULAR_SEED=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L7d/src_atomic_HeH2.13/output \
    EXHALE_MOLECULAR_SEED_X2=local \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L7e.x > seed.log 2>&1

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 \
    EXHALE_OUTER_PASSES=40 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L7e.x > run.log 2>&1

# the first solve returned info=1: the written state reloaded and solved
# again with the pseudo-time start raised (L4e)
cp output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
cp output/Ion_species.txt output/Ion_species_IC.txt
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e8 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L7e.x > run.log 2>&1
```

## What came out

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 2.86E-07 | 5.23E-11 | 1.27E-06 | 5.21E-08 of 1.0E-05 (carrier balance H2) | 350 | 0.51 |
| 2 | 0 | 1.57E-07 | 2.84E-11 | 6.93E-07 | 4.97E-08 of 1.0E-05 (carrier balance H2) | 351 | 2.30 |
| 3 | 0 | 1.56E-07 | 2.85E-11 | 6.94E-07 | 4.75E-08 of 1.0E-05 (carrier balance H2) | 352 | 1.61 |
| 4 | 0 | 1.55E-07 | 2.85E-11 | 6.97E-07 | 4.55E-08 of 1.0E-05 (carrier balance H2) | 353 | 1.54 |
| 5 | 0 | 1.56E-07 | 2.84E-11 | 6.97E-07 | 4.36E-08 of 1.0E-05 (carrier balance H2) | 354 | 1.54 |
| 6 | 0 | 1.56E-07 | 2.85E-11 | 6.96E-07 | 4.19E-08 of 1.0E-05 (carrier balance H2) | 355 | 1.54 |
| 7 | 0 | 1.56E-07 | 2.85E-11 | 6.95E-07 | 4.04E-08 of 1.0E-05 (carrier balance H2) | 356 | 1.53 |
| 8 | 0 | 1.55E-07 | 2.84E-11 | 6.98E-07 | 3.89E-08 of 1.0E-05 (carrier balance H2) | 356 | 1.53 |
| 9 | 0 | 1.55E-07 | 2.85E-11 | 6.96E-07 | 3.76E-08 of 1.0E-05 (carrier balance H2) | 357 | 1.53 |
| 10 | 0 | 1.57E-07 | 2.84E-11 | 6.94E-07 | 3.64E-08 of 1.0E-05 (carrier balance H2) | 358 | 1.53 |
| 11 | 0 | 1.55E-07 | 2.85E-11 | 6.93E-07 | 3.52E-08 of 1.0E-05 (carrier balance H2) | 359 | 1.53 |
| 12 | 0 | 1.56E-07 | 2.85E-11 | 6.94E-07 | 3.42E-08 of 1.0E-05 (carrier balance H2) | 359 | 2.16 |
| 13 | 0 | 1.55E-07 | 2.85E-11 | 6.95E-07 | 3.32E-08 of 1.0E-05 (carrier balance H2) | 360 | 2.16 |
| 14 | 0 | 1.57E-07 | 2.84E-11 | 6.97E-07 | 3.22E-08 of 1.0E-05 (carrier balance H2) | 361 | 2.16 |
| 15 | 0 | 1.56E-07 | 2.85E-11 | 6.96E-07 | 3.14E-08 of 1.0E-05 (carrier balance H2) | 361 | 2.15 |
| 16 | 0 | 1.56E-07 | 2.85E-11 | 6.95E-07 | 3.05E-08 of 1.0E-05 (carrier balance H2) | 362 | 2.15 |
| 17 | 0 | 1.55E-07 | 2.85E-11 | 6.95E-07 | 2.97E-08 of 1.0E-05 (carrier balance H2) | 362 | 2.16 |
| 18 | 0 | 1.55E-07 | 2.85E-11 | 6.99E-07 | 2.90E-08 of 1.0E-05 (carrier balance H2) | 363 | 2.16 |
| 19 | 0 | 1.56E-07 | 2.84E-11 | 6.95E-07 | 2.83E-08 of 1.0E-05 (carrier balance H2) | 364 | 2.16 |
| 20 | 0 | 1.56E-07 | 2.85E-11 | 6.96E-07 | 2.77E-08 of 1.0E-05 (carrier balance H2) | 364 | 2.14 |
| 21 | 0 | 1.56E-07 | 2.85E-11 | 6.95E-07 | 2.70E-08 of 1.0E-05 (carrier balance H2) | 365 | 2.16 |
| 22 | 0 | 1.56E-07 | 2.85E-11 | 6.95E-07 | 2.64E-08 of 1.0E-05 (carrier balance H2) | 365 | 2.16 |
| 23 | 0 | 1.55E-07 | 2.85E-11 | 6.95E-07 | 2.59E-08 of 1.0E-05 (carrier balance H2) | 366 | 2.16 |
| 24 | 0 | 1.57E-07 | 2.85E-11 | 6.94E-07 | 2.53E-08 of 1.0E-05 (carrier balance H2) | 366 | 2.16 |
| 25 | 0 | 1.57E-07 | 2.84E-11 | 6.91E-07 | 2.48E-08 of 1.0E-05 (carrier balance H2) | 367 | 2.16 |
| 26 | 0 | 1.57E-07 | 2.84E-11 | 6.95E-07 | 2.43E-08 of 1.0E-05 (carrier balance H2) | 367 | 2.16 |
| 27 | 0 | 1.57E-07 | 2.84E-11 | 6.95E-07 | 2.39E-08 of 1.0E-05 (carrier balance H2) | 368 | 2.15 |
| 28 | 0 | 1.56E-07 | 2.85E-11 | 6.95E-07 | 2.34E-08 of 1.0E-05 (carrier balance H2) | 368 | 2.16 |
| 29 | 0 | 1.57E-07 | 2.84E-11 | 6.92E-07 | 2.30E-08 of 1.0E-05 (carrier balance H2) | 369 | 2.16 |
| 30 | 0 | 1.55E-07 | 2.85E-11 | 6.94E-07 | 2.26E-08 of 1.0E-05 (carrier balance H2) | 369 | 2.16 |
| 31 | 0 | 1.56E-07 | 2.85E-11 | 6.98E-07 | 2.22E-08 of 1.0E-05 (carrier balance H2) | 370 | 2.15 |
| 32 | 0 | 1.55E-07 | 2.85E-11 | 6.96E-07 | 2.18E-08 of 1.0E-05 (carrier balance H2) | 370 | 2.16 |
| 33 | 0 | 1.57E-07 | 2.84E-11 | 6.92E-07 | 2.15E-08 of 1.0E-05 (carrier balance H2) | 370 | 2.17 |
| 34 | 0 | 1.56E-07 | 2.84E-11 | 6.95E-07 | 2.11E-08 of 1.0E-05 (carrier balance H2) | 371 | 2.17 |
| 35 | 0 | 1.56E-07 | 2.85E-11 | 6.95E-07 | 2.08E-08 of 1.0E-05 (carrier balance H2) | 371 | 2.16 |
| 36 | 0 | 1.56E-07 | 2.85E-11 | 6.94E-07 | 2.05E-08 of 1.0E-05 (carrier balance H2) | 372 | 2.15 |
| 37 | 0 | 1.57E-07 | 2.85E-11 | 6.94E-07 | 2.02E-08 of 1.0E-05 (carrier balance H2) | 372 | 2.17 |
| 38 | 0 | 1.56E-07 | 2.85E-11 | 6.94E-07 | 1.99E-08 of 1.0E-05 (carrier balance H2) | 372 | 2.18 |
| 39 | 0 | 1.57E-07 | 2.85E-11 | 6.93E-07 | 1.96E-08 of 1.0E-05 (carrier balance H2) | 373 | 2.17 |
| 40 | 0 | 1.56E-07 | 2.85E-11 | 6.94E-07 | 1.93E-08 of 1.0E-05 (carrier balance H2) | 373 | 2.02 |

- solver verdict: **info = 1** (last `||R||` = 1.039E-08), by the partitioned stationary route
- certification of the state written: **NOT CERTIFIED**
  - hydrodynamic mass row: row measure  1.555E-07 above  3.7E-09 at cell 1
- wall clock **55m05s** at `OMP_NUM_THREADS=8`, 2026-09-15 12:28:50 to 2026-09-15 13:23:55

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
