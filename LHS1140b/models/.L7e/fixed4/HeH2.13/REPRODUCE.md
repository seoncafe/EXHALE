# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`.L7e/fixed4/HeH2.13`


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
| 1 | 0 | 1.25E-09 | 4.09E-13 | 1.33E-08 | 4.11E-02 of 1.0E-05 (carrier balance H2) | 500 | 3.29 |
| 2 | 2 | 7.11E-01 | 1.27E-04 | 1.36E+00 | 3.50E-01 of 1.0E-05 (carrier balance H2) | 217 | 373.51 |
| 3 | 0 | 1.27E-07 | 7.75E-12 | 9.77E-07 | 2.52E-01 of 1.0E-05 (carrier balance H2) | 250 | 116.39 |
| 4 | 2 | 9.82E-01 | 1.84E-04 | 1.52E+00 | 3.56E-01 of 1.0E-05 (carrier balance H2) | 218 | 371.49 |
| 5 | 2 | 1.28E+00 | 1.05E-02 | 1.22E+00 | 1.49E-02 of 1.0E-05 (carrier balance H2) | 456 | 319.80 |
| 6 | 2 | 1.81E+00 | 2.69E-03 | 1.46E+00 | 2.16E-01 of 1.0E-05 (carrier balance H2) | 256 | 283.78 |
| 7 | 2 | 1.99E+00 | 7.75E-03 | 1.08E+00 | 6.97E-02 of 1.0E-05 (carrier balance H2) | 222 | 272.26 |
| 8 | 2 | 1.99E+00 | 6.08E-03 | 1.36E+00 | 1.21E-01 of 1.0E-05 (carrier balance H2) | 255 | 182.58 |
| 9 | 2 | 1.97E+00 | 4.62E-03 | 1.14E+00 | 2.07E-01 of 1.0E-05 (carrier balance H2) | 255 | 216.75 |
| 10 | 2 | 1.98E+00 | 3.28E-03 | 1.15E+00 | 2.31E-01 of 1.0E-05 (carrier balance H2) | 255 | 206.38 |
| 11 | 2 | 1.96E+00 | 2.90E-03 | 1.06E+00 | 2.64E-01 of 1.0E-05 (carrier balance H2) | 255 | 225.10 |
| 12 | 2 | 1.97E+00 | 4.27E-03 | 1.06E+00 | 2.54E-01 of 1.0E-05 (carrier balance H2) | 255 | 181.48 |
| 13 | 2 | 1.99E+00 | 3.92E-03 | 1.22E+00 | 2.66E-01 of 1.0E-05 (carrier balance H2) | 255 | 218.74 |
| 14 | 2 | 1.94E+00 | 4.33E-03 | 1.20E+00 | 2.36E-01 of 1.0E-05 (carrier balance H2) | 255 | 186.01 |
| 15 | 2 | 2.00E+00 | 4.29E-03 | 1.18E+00 | 2.37E-01 of 1.0E-05 (carrier balance H2) | 255 | 191.83 |
| 16 | 2 | 1.97E+00 | 4.08E-02 | 1.01E+00 | 1.84E-01 of 1.0E-05 (carrier balance H2) | 225 | 269.63 |
| 17 | 2 | 1.98E+00 | 2.96E-02 | 1.02E+00 | 2.10E-01 of 1.0E-05 (carrier balance H2) | 225 | 179.69 |
| 18 | 2 | 1.98E+00 | 2.01E-02 | 1.07E+00 | 2.13E-01 of 1.0E-05 (carrier balance H2) | 225 | 183.64 |
| 19 | 2 | 1.92E+00 | 1.69E-02 | 1.04E+00 | 2.16E-01 of 1.0E-05 (carrier balance H2) | 225 | 187.02 |
| 20 | 2 | 1.99E+00 | 1.55E-02 | 1.15E+00 | 2.19E-01 of 1.0E-05 (carrier balance H2) | 225 | 186.43 |
| 21 | 2 | 1.94E+00 | 1.46E-02 | 1.18E+00 | 2.22E-01 of 1.0E-05 (carrier balance H2) | 225 | 187.59 |
| 22 | 2 | 1.93E+00 | 1.39E-02 | 1.12E+00 | 2.24E-01 of 1.0E-05 (carrier balance H2) | 226 | 189.95 |
| 23 | 2 | 1.93E+00 | 1.37E-02 | 1.12E+00 | 2.27E-01 of 1.0E-05 (carrier balance H2) | 226 | 189.76 |
| 24 | 2 | 1.93E+00 | 1.33E-02 | 1.12E+00 | 2.30E-01 of 1.0E-05 (carrier balance H2) | 226 | 187.33 |
| 25 | 2 | 1.93E+00 | 1.28E-02 | 1.10E+00 | 2.31E-01 of 1.0E-05 (carrier balance H2) | 226 | 185.75 |
| 26 | 2 | 1.96E+00 | 2.09E-02 | 1.16E+00 | 2.24E-01 of 1.0E-05 (carrier balance H2) | 226 | 259.04 |
| 27 | 2 | 1.99E+00 | 1.73E-02 | 1.07E+00 | 2.25E-01 of 1.0E-05 (carrier balance H2) | 226 | 180.96 |
| 28 | 2 | 1.98E+00 | 2.65E-02 | 1.10E+00 | 1.93E-01 of 1.0E-05 (carrier balance H2) | 225 | 204.66 |
| 29 | 2 | 2.00E+00 | 3.04E-02 | 1.35E+00 | 1.77E-01 of 1.0E-05 (carrier balance H2) | 225 | 196.97 |
| 30 | 2 | 2.00E+00 | 3.47E-02 | 1.14E+00 | 1.38E-01 of 1.0E-05 (carrier balance H2) | 224 | 245.91 |
| 31 | 2 | 2.00E+00 | 3.34E-02 | 1.10E+00 | 1.14E-01 of 1.0E-05 (carrier balance H2) | 253 | 202.42 |
| 32 | 2 | 2.00E+00 | 4.14E-02 | 1.12E+00 | 2.19E-01 of 1.0E-05 (carrier balance H2) | 253 | 197.65 |
| 33 | 2 | 2.00E+00 | 4.38E-02 | 1.17E+00 | 3.23E-01 of 1.0E-05 (carrier balance H2) | 253 | 209.14 |
| 34 | 2 | 2.00E+00 | 7.91E-02 | 1.10E+00 | 4.08E-01 of 1.0E-05 (carrier balance H2) | 253 | 201.11 |

- solver verdict: **info = 1** (last `||R||` = 2.000E+00), by the partitioned stationary route
- certification of the state written: **NOT CERTIFIED**
  - hydrodynamic mass row: row measure  1.997E+00 above  3.2E-12 at cell 139
  - hydrodynamic momentum row: row measure  7.906E-02 above  1.0E-08 at cell 129
  - hydrodynamic energy row: row measure  1.101E+00 above  1.0E-06 at cell 162
  - carrier balance H2: gated row measure  4.081E-01 above  1.0E-05 at cell 253 (a wind cell)
  - elemental transport He/H partition: gated row measure  1.149E-04 above  1.0E-05 at cell 217 (a wind cell)
- wall clock **226m37s** at `OMP_NUM_THREADS=8`, 2026-09-15 11:44:06 to 2026-09-15 15:30:44

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
