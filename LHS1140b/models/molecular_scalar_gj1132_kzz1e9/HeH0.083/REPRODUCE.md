# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`molecular_scalar_gj1132_kzz1e9/HeH0.083`

- chemistry: H2, H2+, H3+ and HeH+ with the H2 carrier transported (`molecular`)
- lower boundary: prescribed scalar base, T = 226 K, R_0 = 0.157692 R_J, p = 1 microbar, metal-free (`scalar`)
- spectrum: `gj1132`
- mixing: binary H/He element diffusion with He_Kzz = 1e9 cm^2/s (`kzz1e9`)

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_L7f.x` |
| md5 | `c2e9c9990b9f14f1be8cd77abca68945` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_molecular_scalar_gj1132_kzz1e9_HeH0.083` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083/output (given); then a continuation from the first solve's written state at dtau0=1.0e8 (its log run_dtau0_first.log)
```

The fields are the state directory, the tier it was taken from (`tier0` this same case's own most recent certified state -- the one this tree carries, else the one a preserved tree such as `models_20260914_preL21/` holds, which is what a re-solve of the whole catalog continues from -- `tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

The grid the seed is written on is named by `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt`.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/molecular_scalar_gj1132_kzz1e9/HeH0.083/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 

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
| 1 | 0 | 2.37E-08 | 3.64E-12 | 2.19E-08 | 1.04E-03 of 1.0E-05 (elemental transport He/H partition) | 340 | 5.80 |
| 2 | 0 | 2.99E-09 | 8.43E-13 | 7.65E-09 | 3.21E-03 of 1.0E-05 (carrier balance H2) | 500 | 602.48 |
| 3 | 0 | 5.46E-09 | 9.83E-13 | 9.07E-09 | 1.34E-02 of 1.0E-05 (carrier balance H2) | 500 | 179.02 |
| 4 | 0 | 1.15E-08 | 5.27E-13 | 6.43E-09 | 5.27E-03 of 1.0E-05 (carrier balance H2) | 500 | 357.70 |
| 5 | 0 | 5.63E-09 | 1.55E-12 | 1.13E-08 | 7.33E-03 of 1.0E-05 (carrier balance H2) | 500 | 99.72 |
| 6 | 0 | 1.94E-09 | 9.10E-13 | 5.61E-09 | 4.57E-03 of 1.0E-05 (carrier balance H2) | 500 | 398.75 |
| 7 | 0 | 4.13E-09 | 6.43E-13 | 7.35E-09 | 3.50E-03 of 1.0E-05 (carrier balance H2) | 500 | 168.08 |
| 8 | 0 | 9.89E-09 | 6.88E-13 | 4.61E-09 | 4.23E-03 of 1.0E-05 (carrier balance H2) | 500 | 425.77 |
| 9 | 0 | 2.73E-09 | 7.89E-13 | 6.22E-09 | 2.56E-03 of 1.0E-05 (carrier balance H2) | 500 | 98.87 |
| 10 | 0 | 4.70E-09 | 1.10E-12 | 9.79E-09 | 2.16E-03 of 1.0E-05 (carrier balance H2) | 500 | 248.60 |
| 11 | 0 | 3.11E-09 | 5.52E-13 | 6.52E-09 | 2.07E-03 of 1.0E-05 (carrier balance H2) | 500 | 331.06 |
| 12 | 0 | 1.07E-08 | 5.75E-13 | 6.65E-09 | 2.12E-03 of 1.0E-05 (carrier balance H2) | 500 | 711.47 |
| 13 | 0 | 1.59E-09 | 3.53E-13 | 6.33E-09 | 2.20E-03 of 1.0E-05 (carrier balance H2) | 500 | 385.25 |
| 14 | 0 | 6.75E-09 | 1.12E-12 | 1.01E-08 | 2.28E-03 of 1.0E-05 (carrier balance H2) | 500 | 301.29 |
| 15 | 2 | 1.83E-08 | 5.26E-13 | 5.64E-09 | 2.39E-03 of 1.0E-05 (carrier balance H2) | 500 | 105.63 |
| 16 | 0 | 2.94E-09 | 3.73E-13 | 4.67E-09 | 2.47E-03 of 1.0E-05 (carrier balance H2) | 500 | 119.00 |
| 17 | 0 | 8.11E-09 | 6.08E-13 | 5.08E-09 | 2.58E-03 of 1.0E-05 (carrier balance H2) | 500 | 105.59 |
| 18 | 0 | 2.64E-09 | 7.33E-13 | 5.60E-09 | 2.65E-03 of 1.0E-05 (carrier balance H2) | 500 | 236.47 |
| 19 | 0 | 4.37E-09 | 8.07E-13 | 7.50E-09 | 2.76E-03 of 1.0E-05 (carrier balance H2) | 500 | 105.90 |
| 20 | 0 | 9.90E-09 | 7.37E-13 | 6.27E-09 | 2.86E-03 of 1.0E-05 (carrier balance H2) | 500 | 105.98 |
| 21 | 0 | 5.14E-09 | 4.13E-13 | 5.83E-09 | 2.97E-03 of 1.0E-05 (carrier balance H2) | 500 | 548.35 |
| 22 | 0 | 2.81E-09 | 7.24E-13 | 5.91E-09 | 3.07E-03 of 1.0E-05 (carrier balance H2) | 500 | 104.86 |
| 23 | 0 | 4.33E-09 | 7.56E-13 | 5.50E-09 | 3.17E-03 of 1.0E-05 (carrier balance H2) | 500 | 106.33 |
| 24 | 0 | 4.10E-09 | 7.18E-13 | 6.73E-09 | 3.29E-03 of 1.0E-05 (carrier balance H2) | 500 | 127.44 |
| 25 | 0 | 6.71E-09 | 8.14E-13 | 7.50E-09 | 3.27E-03 of 1.0E-05 (carrier balance H2) | 500 | 138.84 |
| 26 | 0 | 3.43E-09 | 7.11E-13 | 7.61E-09 | 3.27E-03 of 1.0E-05 (carrier balance H2) | 500 | 105.99 |
| 27 | 0 | 4.67E-09 | 8.26E-13 | 1.03E-08 | 3.25E-03 of 1.0E-05 (carrier balance H2) | 500 | 121.19 |
| 28 | 0 | 3.46E-09 | 8.32E-13 | 6.33E-09 | 3.26E-03 of 1.0E-05 (carrier balance H2) | 500 | 234.01 |
| 29 | 0 | 3.55E-09 | 6.59E-13 | 6.32E-09 | 3.21E-03 of 1.0E-05 (carrier balance H2) | 500 | 119.30 |
| 30 | 0 | 1.53E-08 | 6.13E-13 | 6.07E-09 | 3.20E-03 of 1.0E-05 (carrier balance H2) | 500 | 118.89 |
| 31 | 0 | 6.14E-09 | 1.24E-12 | 1.04E-08 | 3.21E-03 of 1.0E-05 (carrier balance H2) | 500 | 509.71 |
| 32 | 0 | 1.88E-08 | 4.68E-13 | 5.46E-09 | 3.15E-03 of 1.0E-05 (carrier balance H2) | 500 | 102.73 |
| 33 | 0 | 1.34E-08 | 2.58E-13 | 7.44E-09 | 3.12E-03 of 1.0E-05 (carrier balance H2) | 500 | 321.36 |
| 34 | 0 | 1.29E-08 | 5.85E-13 | 6.28E-09 | 3.11E-03 of 1.0E-05 (carrier balance H2) | 500 | 102.48 |
| 35 | 0 | 7.03E-09 | 8.20E-13 | 5.80E-09 | 3.07E-03 of 1.0E-05 (carrier balance H2) | 500 | 102.47 |
| 36 | 0 | 6.33E-09 | 7.66E-13 | 7.88E-09 | 3.06E-03 of 1.0E-05 (carrier balance H2) | 500 | 102.70 |
| 37 | 0 | 6.43E-09 | 1.61E-12 | 1.35E-08 | 3.05E-03 of 1.0E-05 (carrier balance H2) | 500 | 92.73 |
| 38 | 0 | 6.49E-09 | 1.05E-12 | 1.20E-08 | 3.01E-03 of 1.0E-05 (carrier balance H2) | 500 | 93.29 |
| 39 | 0 | 8.38E-09 | 5.05E-13 | 5.28E-09 | 2.99E-03 of 1.0E-05 (carrier balance H2) | 500 | 102.22 |
| 40 | 0 | 1.51E-08 | 4.74E-13 | 7.16E-09 | 2.97E-03 of 1.0E-05 (carrier balance H2) | 500 | 98.53 |

- solver verdict: **info = 1** (last `||R||` = 1.369E-08), by the partitioned stationary route
- certification of the state written: **NOT CERTIFIED**
  - carrier balance H2: gated row measure  2.972E-03 above  1.0E-05 at cell 500 (a wind cell)
  - elemental transport He/H partition: gated row measure  1.308E-03 above  1.0E-05 at cell 337 (a wind cell)
- wall clock **233m16s** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-16 07:40:13 to 2026-09-16 11:33:29

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

## Re-measuring this state, and how closely it comes back

The two files of `output/` can be handed back to the binary as `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with `Restart intent: stationary evaluate`, which measures the state as it stands, writes it back unchanged and exits 0 only if every active equation is within its tolerance.

What comes back exactly and what does not, MEASURED over the 120 certified states of `models/` (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8). The conserved state round-trips to the last bit: the radius, velocity and pressure columns are read as written, and since item L18 the mass density is read from its own column rather than rebuilt from the species, so the state re-measured is the state certified. The ROW MEASURES do not round-trip bitwise and cannot: the first equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at its rounding floor, where one unit in the last place of the density moves the cell-wise maximum of the energy row by about 11 per cent. Re-evaluated against in-run, the mass row comes back within a factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to 1.93 (median 1.22); the momentum row sits four to six decades below its tolerance and its ratio is the floor itself. A state certified at more than about half of its tolerance may therefore be refused on re-evaluation, and that refusal is a property of the cancellation and not of this file.
