# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`molecular_scalar_gj1132_wellmixed/HeH0.083`

- chemistry: H2, H2+, H3+ and HeH+ with the H2 carrier transported (`molecular`)
- lower boundary: prescribed scalar base, T = 226 K, R_0 = 0.157692 R_J, p = 1 microbar, metal-free (`scalar`)
- spectrum: `gj1132`
- mixing: no element diffusion; He/H uniform (`wellmixed`)

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L7f.x` |
| md5 | `42ec3a6cbbd4b010300510609fc304f0` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_molecular_scalar_gj1132_wellmixed_HeH0.083` |

## The seed

The binary built this case's state out of a certified ATOMIC solution, which is how a molecular case is seeded (`EXHALE_MOLECULAR_SEED`; `models/pick_seed.py` carries no molecular state to choose from):

```
molecular seed from /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132_wellmixed/HeH0.083/output, x2 local; then a continuation from the first solve's written state at dtau0=1.0e8 (its log run_dtau0_first.log)
```

The fields are the atomic state directory and which extension of the base H2 partition the conversion carried: `local` gives every cell its own chemical-equilibrium q_H2(p,T), `handoff` carries the base x2 to the outer boundary, a number is that fraction everywhere.

The conversion, which writes `output/*_IC.txt` and stops:

```
# with "Load IC? True" in input.inp, which this case states
OMP_NUM_THREADS=8 EXHALE_MOLECULAR_SEED=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132_wellmixed/HeH0.083/output \
    EXHALE_MOLECULAR_SEED_X2=local \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L7f.x > seed.log 2>&1
```

## The commands

```bash
# the molecular seed, built by the binary from the certified atomic state
# with "Load IC? True" in input.inp, which this case states
OMP_NUM_THREADS=8 EXHALE_MOLECULAR_SEED=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132_wellmixed/HeH0.083/output \
    EXHALE_MOLECULAR_SEED_X2=local \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L7f.x > seed.log 2>&1

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 \
    EXHALE_OUTER_PASSES=40 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L7f.x > run.log 2>&1

# the first solve returned info=1: the written state reloaded and solved
# again with the pseudo-time start raised (L4e)
cp output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
cp output/Ion_species.txt output/Ion_species_IC.txt
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e8 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L7f.x > run.log 2>&1
```

## What came out

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 5.97E-06 | 2.33E-10 | 1.34E-05 | 4.79E-02 of 1.0E-05 (carrier balance H2) | 342 | 4.44 |
| 2 | 0 | 1.00E-08 | 5.48E-13 | 1.78E-08 | 4.75E-02 of 1.0E-05 (carrier balance H2) | 343 | 308.17 |
| 3 | 0 | 9.40E-09 | 5.65E-13 | 1.76E-08 | 4.62E-02 of 1.0E-05 (carrier balance H2) | 344 | 322.07 |
| 4 | 0 | 9.69E-09 | 3.26E-13 | 1.38E-08 | 4.42E-02 of 1.0E-05 (carrier balance H2) | 345 | 320.45 |
| 5 | 0 | 7.82E-09 | 5.33E-13 | 1.31E-08 | 4.21E-02 of 1.0E-05 (carrier balance H2) | 346 | 308.84 |
| 6 | 0 | 4.30E-09 | 1.72E-13 | 1.05E-08 | 3.99E-02 of 1.0E-05 (carrier balance H2) | 347 | 320.76 |
| 7 | 0 | 6.42E-09 | 3.53E-13 | 1.34E-08 | 3.77E-02 of 1.0E-05 (carrier balance H2) | 347 | 251.59 |
| 8 | 0 | 2.49E-09 | 8.12E-13 | 1.08E-08 | 3.57E-02 of 1.0E-05 (carrier balance H2) | 348 | 332.90 |
| 9 | 0 | 4.83E-09 | 6.19E-13 | 1.26E-08 | 3.36E-02 of 1.0E-05 (carrier balance H2) | 348 | 358.74 |
| 10 | 0 | 3.81E-09 | 7.45E-13 | 1.37E-08 | 3.16E-02 of 1.0E-05 (carrier balance H2) | 348 | 630.23 |
| 11 | 0 | 6.88E-09 | 4.74E-13 | 1.47E-08 | 2.96E-02 of 1.0E-05 (carrier balance H2) | 348 | 374.27 |
| 12 | 2 | 1.39E-08 | 1.25E-12 | 3.04E-08 | 2.78E-02 of 1.0E-05 (carrier balance H2) | 347 | 353.83 |
| 13 | 0 | 7.54E-09 | 5.74E-13 | 2.27E-08 | 2.60E-02 of 1.0E-05 (carrier balance H2) | 346 | 355.11 |
| 14 | 0 | 1.72E-08 | 8.91E-13 | 1.51E-08 | 2.69E-02 of 1.0E-05 (carrier balance H2) | 500 | 1213.03 |
| 15 | 0 | 7.15E-09 | 9.66E-13 | 2.09E-08 | 2.48E-02 of 1.0E-05 (carrier balance H2) | 500 | 227.27 |
| 16 | 0 | 7.06E-09 | 4.41E-13 | 2.02E-08 | 2.56E-02 of 1.0E-05 (carrier balance H2) | 500 | 222.87 |
| 17 | 0 | 8.98E-09 | 6.72E-13 | 2.48E-08 | 2.44E-02 of 1.0E-05 (carrier balance H2) | 500 | 222.13 |
| 18 | 0 | 4.20E-09 | 5.85E-13 | 1.34E-08 | 2.47E-02 of 1.0E-05 (carrier balance H2) | 500 | 224.66 |
| 19 | 0 | 1.22E-08 | 8.39E-13 | 1.92E-08 | 2.41E-02 of 1.0E-05 (carrier balance H2) | 500 | 258.09 |
| 20 | 0 | 1.15E-08 | 1.03E-12 | 3.60E-08 | 2.42E-02 of 1.0E-05 (carrier balance H2) | 500 | 216.55 |
| 21 | 0 | 1.25E-08 | 1.13E-12 | 3.25E-08 | 2.42E-02 of 1.0E-05 (carrier balance H2) | 500 | 412.41 |
| 22 | 0 | 6.36E-09 | 4.29E-13 | 1.66E-08 | 2.43E-02 of 1.0E-05 (carrier balance H2) | 500 | 440.45 |
| 23 | 0 | 2.49E-09 | 4.50E-13 | 1.15E-08 | 2.45E-02 of 1.0E-05 (carrier balance H2) | 500 | 256.76 |
| 24 | 0 | 1.06E-08 | 1.06E-12 | 3.11E-08 | 2.46E-02 of 1.0E-05 (carrier balance H2) | 500 | 364.11 |

- solver verdict: **info = 1** (last `||R||` = 2.033E-08), by the partitioned stationary route
- certification of the state written: **NOT CERTIFIED**
  - carrier balance H2: gated row measure  2.462E-02 above  1.0E-05 at cell 500 (a wind cell)
- wall clock **318m09s** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-15 22:46:23 to 2026-09-16 04:04:32

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

## Re-measuring this state, and how closely it comes back

The two files of `output/` can be handed back to the binary as `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with `Restart intent: stationary evaluate`, which measures the state as it stands, writes it back unchanged and exits 0 only if every active equation is within its tolerance.

What comes back exactly and what does not, MEASURED over the 120 certified states of `models/` (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8). The conserved state round-trips to the last bit: the radius, velocity and pressure columns are read as written, and since item L18 the mass density is read from its own column rather than rebuilt from the species, so the state re-measured is the state certified. The ROW MEASURES do not round-trip bitwise and cannot: the first equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at its rounding floor, where one unit in the last place of the density moves the cell-wise maximum of the energy row by about 11 per cent. Re-evaluated against in-run, the mass row comes back within a factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to 1.93 (median 1.22); the momentum row sits four to six decades below its tolerance and its ratio is the floor itself. A state certified at more than about half of its tolerance may therefore be refused on re-evaluation, and that refusal is a property of the cancellation and not of this file.
