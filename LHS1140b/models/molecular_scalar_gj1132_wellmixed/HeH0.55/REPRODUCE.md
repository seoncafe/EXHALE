# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`molecular_scalar_gj1132_wellmixed/HeH0.55`

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
| run title | `Simulation for LHS1140b_molecular_scalar_gj1132_wellmixed_HeH0.55` |

## The seed

The binary built this case's state out of a certified ATOMIC solution, which is how a molecular case is seeded (`EXHALE_MOLECULAR_SEED`; `models/pick_seed.py` carries no molecular state to choose from):

```
molecular seed from /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132_wellmixed/HeH0.55/output, x2 local; then a continuation from the first solve's written state at dtau0=1.0e8 (its log run_dtau0_first.log)
```

The fields are the atomic state directory and which extension of the base H2 partition the conversion carried: `local` gives every cell its own chemical-equilibrium q_H2(p,T), `handoff` carries the base x2 to the outer boundary, a number is that fraction everywhere.

The conversion, which writes `output/*_IC.txt` and stops:

```
# with "Load IC? True" in input.inp, which this case states
OMP_NUM_THREADS=8 EXHALE_MOLECULAR_SEED=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132_wellmixed/HeH0.55/output \
    EXHALE_MOLECULAR_SEED_X2=local \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L7f.x > seed.log 2>&1
```

## The commands

```bash
# the molecular seed, built by the binary from the certified atomic state
# with "Load IC? True" in input.inp, which this case states
OMP_NUM_THREADS=8 EXHALE_MOLECULAR_SEED=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132_wellmixed/HeH0.55/output \
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
| 1 | 0 | 2.25E-06 | 8.51E-11 | 4.63E-06 | 1.27E-02 of 1.0E-05 (carrier balance H2) | 233 | 4.84 |
| 2 | 0 | 7.58E-09 | 8.89E-14 | 1.59E-08 | 8.88E-02 of 1.0E-05 (carrier balance H2) | 244 | 879.13 |
| 3 | 0 | 8.65E-09 | 1.70E-13 | 2.01E-08 | 8.32E-02 of 1.0E-05 (carrier balance H2) | 245 | 337.56 |
| 4 | 0 | 1.02E-08 | 2.25E-13 | 1.83E-08 | 8.17E-02 of 1.0E-05 (carrier balance H2) | 247 | 1431.34 |
| 5 | 0 | 1.33E-08 | 2.52E-13 | 1.80E-08 | 7.59E-02 of 1.0E-05 (carrier balance H2) | 248 | 786.96 |
| 6 | 2 | 3.82E-09 | 9.99E-13 | 9.44E-09 | 7.35E-02 of 1.0E-05 (carrier balance H2) | 250 | 1343.01 |
| 7 | 0 | 3.55E-09 | 4.79E-13 | 1.12E-08 | 6.53E-02 of 1.0E-05 (carrier balance H2) | 252 | 330.14 |
| 8 | 0 | 1.40E-08 | 2.53E-13 | 1.73E-08 | 5.68E-02 of 1.0E-05 (carrier balance H2) | 254 | 899.98 |
| 9 | 2 | 6.35E-09 | 1.29E-12 | 1.25E-08 | 4.74E-02 of 1.0E-05 (carrier balance H2) | 256 | 1073.92 |
| 10 | 0 | 3.82E-09 | 9.98E-13 | 1.64E-08 | 3.93E-02 of 1.0E-05 (carrier balance H2) | 258 | 773.16 |
| 11 | 0 | 1.41E-08 | 4.78E-14 | 1.32E-08 | 3.71E-02 of 1.0E-05 (carrier balance H2) | 251 | 408.18 |
| 12 | 0 | 7.17E-09 | 6.89E-13 | 2.34E-08 | 3.67E-02 of 1.0E-05 (carrier balance H2) | 254 | 263.18 |
| 13 | 0 | 6.02E-09 | 8.63E-13 | 1.94E-08 | 3.69E-02 of 1.0E-05 (carrier balance H2) | 257 | 794.70 |
| 14 | 0 | 3.85E-09 | 6.94E-13 | 1.27E-08 | 3.75E-02 of 1.0E-05 (carrier balance H2) | 258 | 710.67 |
| 15 | 0 | 1.61E-08 | 3.04E-13 | 1.37E-08 | 3.08E-02 of 1.0E-05 (carrier balance H2) | 259 | 624.03 |
| 16 | 0 | 1.99E-08 | 3.44E-13 | 1.03E-08 | 2.93E-02 of 1.0E-05 (carrier balance H2) | 261 | 683.09 |
| 17 | 0 | 1.65E-08 | 2.79E-13 | 1.25E-08 | 2.82E-02 of 1.0E-05 (carrier balance H2) | 264 | 445.79 |
| 18 | 0 | 1.96E-08 | 3.51E-14 | 1.27E-08 | 2.76E-02 of 1.0E-05 (carrier balance H2) | 266 | 612.31 |
| 19 | 0 | 9.80E-09 | 8.84E-13 | 3.29E-08 | 2.67E-02 of 1.0E-05 (carrier balance H2) | 267 | 575.48 |
| 20 | 0 | 9.68E-09 | 7.94E-13 | 2.96E-08 | 2.70E-02 of 1.0E-05 (carrier balance H2) | 269 | 676.86 |
| 21 | 0 | 1.91E-08 | 1.77E-13 | 1.16E-08 | 2.43E-02 of 1.0E-05 (carrier balance H2) | 269 | 283.92 |
| 22 | 0 | 2.38E-08 | 3.26E-13 | 1.04E-08 | 2.36E-02 of 1.0E-05 (carrier balance H2) | 271 | 180.61 |
| 23 | 0 | 9.89E-09 | 9.45E-13 | 2.88E-08 | 2.31E-02 of 1.0E-05 (carrier balance H2) | 273 | 415.62 |
| 24 | 0 | 1.36E-08 | 4.99E-13 | 3.91E-08 | 2.29E-02 of 1.0E-05 (carrier balance H2) | 274 | 451.77 |
| 25 | 0 | 2.52E-08 | 1.37E-13 | 9.76E-09 | 2.27E-02 of 1.0E-05 (carrier balance H2) | 276 | 554.22 |
| 26 | 0 | 1.87E-08 | 4.93E-14 | 1.41E-08 | 2.26E-02 of 1.0E-05 (carrier balance H2) | 278 | 501.85 |
| 27 | 0 | 2.18E-08 | 2.94E-13 | 1.58E-08 | 2.23E-02 of 1.0E-05 (carrier balance H2) | 278 | 336.62 |
| 28 | 0 | 2.82E-08 | 1.75E-13 | 1.36E-08 | 2.24E-02 of 1.0E-05 (carrier balance H2) | 281 | 211.30 |
| 29 | 0 | 1.47E-08 | 6.56E-14 | 1.04E-08 | 2.27E-02 of 1.0E-05 (carrier balance H2) | 291 | 170.81 |
| 30 | 0 | 1.28E-08 | 1.39E-13 | 8.37E-09 | 2.27E-02 of 1.0E-05 (carrier balance H2) | 293 | 167.48 |
| 31 | 0 | 3.17E-08 | 4.30E-13 | 1.64E-08 | 2.27E-02 of 1.0E-05 (carrier balance H2) | 294 | 194.13 |
| 32 | 0 | 2.51E-08 | 2.14E-13 | 1.43E-08 | 2.27E-02 of 1.0E-05 (carrier balance H2) | 295 | 218.55 |
| 33 | 0 | 2.16E-08 | 2.40E-13 | 8.44E-09 | 2.29E-02 of 1.0E-05 (carrier balance H2) | 297 | 188.02 |
| 34 | 0 | 1.77E-08 | 1.94E-13 | 9.92E-09 | 2.30E-02 of 1.0E-05 (carrier balance H2) | 298 | 209.94 |
| 35 | 0 | 2.59E-08 | 1.84E-13 | 1.43E-08 | 2.32E-02 of 1.0E-05 (carrier balance H2) | 299 | 165.99 |
| 36 | 0 | 2.22E-08 | 1.19E-13 | 1.25E-08 | 2.35E-02 of 1.0E-05 (carrier balance H2) | 300 | 239.26 |
| 37 | 0 | 1.33E-08 | 2.88E-13 | 2.46E-08 | 2.38E-02 of 1.0E-05 (carrier balance H2) | 302 | 277.01 |
| 38 | 0 | 1.99E-08 | 2.04E-13 | 1.88E-08 | 2.41E-02 of 1.0E-05 (carrier balance H2) | 303 | 253.71 |
| 39 | 0 | 1.74E-08 | 2.45E-13 | 2.86E-08 | 2.44E-02 of 1.0E-05 (carrier balance H2) | 304 | 235.30 |
| 40 | 0 | 2.86E-08 | 2.34E-13 | 1.79E-08 | 2.48E-02 of 1.0E-05 (carrier balance H2) | 306 | 229.08 |

- solver verdict: **info = 1** (last `||R||` = 3.484E-08), by the partitioned stationary route
- certification of the state written: **NOT CERTIFIED**
  - hydrodynamic mass row: row measure  2.858E-08 above  2.2E-08 at cell 1
  - carrier balance H2: gated row measure  2.477E-02 above  1.0E-05 at cell 306 (a wind cell)
- wall clock **506m01s** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-15 22:52:23 to 2026-09-16 07:18:24

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

## Re-measuring this state, and how closely it comes back

The two files of `output/` can be handed back to the binary as `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with `Restart intent: stationary evaluate`, which measures the state as it stands, writes it back unchanged and exits 0 only if every active equation is within its tolerance.

What comes back exactly and what does not, MEASURED over the 120 certified states of `models/` (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8). The conserved state round-trips to the last bit: the radius, velocity and pressure columns are read as written, and since item L18 the mass density is read from its own column rather than rebuilt from the species, so the state re-measured is the state certified. The ROW MEASURES do not round-trip bitwise and cannot: the first equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at its rounding floor, where one unit in the last place of the density moves the cell-wise maximum of the energy row by about 11 per cent. Re-evaluated against in-run, the mass row comes back within a factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to 1.93 (median 1.22); the momentum row sits four to six decades below its tolerance and its ratio is the floor itself. A state certified at more than about half of its tolerance may therefore be refused on re-evaluation, and that refusal is a property of the cancellation and not of this file.
