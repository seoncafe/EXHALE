# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`.L7d/t_local`


The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/../../EXHALE_L7d.x` |
| md5 | `9e38a707fc8296d00e2ccf3dc96e0b2c` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (Ubuntu 13.1.0-8ubuntu1~18.04) 13.1.0 |
| linear algebra | LAPACK library: not OpenBLAS (no openblas_* entry points); BLAS threads left to the library |
| run title | `Simulation for LHS1140b_molecular_scalar_gj1132_kzz1e9_HeH2.13` |

## The seed

The binary built this case's state out of a certified ATOMIC solution, which is how a molecular case is seeded (`EXHALE_MOLECULAR_SEED`; `models/pick_seed.py` carries no molecular state to choose from):

```
molecular seed from /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/.L7d/src_atomic_HeH2.13/output, x2 local
```

The fields are the atomic state directory and which extension of the base H2 partition the conversion carried: `local` gives every cell its own chemical-equilibrium q_H2(p,T), `handoff` carries the base x2 to the outer boundary, a number is that fraction everywhere.

The conversion, which writes `output/*_IC.txt` and stops:

```
OMP_NUM_THREADS=8 EXHALE_MOLECULAR_SEED=... EXHALE_L7d.x > seed.log
```

## The commands

```bash
(not recorded)
```

## What came out

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 7.24E-07 | 5.36E-11 | 6.44E-06 | 1.00E+00 of 1.0E-05 (carrier balance H2) | 280 | 5.58 |
| 2 | 0 | 7.24E-07 | 5.36E-11 | 6.45E-06 | 1.00E+00 of 1.0E-05 (carrier balance H2) | 280 | 297.24 |
| 3 | 0 | 7.23E-07 | 5.36E-11 | 6.45E-06 | 1.00E+00 of 1.0E-05 (carrier balance H2) | 280 | 276.53 |
| 4 | 2 | 7.24E-07 | 1.09E-08 | 6.46E-06 | 1.00E+00 of 1.0E-05 (carrier balance H2) | 280 | 451.18 |

- solver verdict: **info = 1** (last `||R||` = 1.619E-08), by the partitioned stationary route
- certification of the state written: **NOT CERTIFIED**
  - hydrodynamic mass row: row measure  8.175E-10 above  3.0E-12 at cell 427
  - hydrodynamic momentum row: row measure  1.091E-08 above  1.0E-08 at cell 428
  - hydrodynamic energy row: row measure  6.465E-06 above  1.0E-06 at cell 2
  - carrier balance H2: gated row measure  1.000E+00 above  1.0E-05 at cell 280 (a wind cell)
- wall clock **0m00s** at `OMP_NUM_THREADS=8`, x to y

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
