# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`molecular_photochem_gj1132_kzzprofile/HeH2.09`

- chemistry: H2, H2+, H3+ and HeH+ with the H2 carrier transported (`molecular`)
- lower boundary: the Photochem column handed over as `Lower atmosphere profile:` (`photochem`)
- spectrum: `gj1132`
- mixing: binary H/He element diffusion with K_zz(p) from the profile (`kzzprofile`)
- numerical flux: `HLLC`

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs two of its keys changed.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `c2e9c9990b9f14f1be8cd77abca68945` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_molecular_photochem_gj1132_kzzprofile_HeH2.09` |

## The seed

The binary built this case's state out of a certified ATOMIC solution, which is how a molecular case is seeded (`EXHALE_MOLECULAR_SEED`; `models/pick_seed.py` carries no molecular state to choose from):

```
molecular seed from /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH2.09/k05/output, x2 local
```

The fields are the atomic state directory and which extension of the base H2 partition the conversion carried: `local` gives every cell its own chemical-equilibrium q_H2(p,T), `handoff` carries the base x2 to the outer boundary, a number is that fraction everywhere.

The conversion, which writes `output/*_IC.txt` and stops:

```
# with "Load IC? True" in input.inp, which this case states
OMP_NUM_THREADS=8 EXHALE_MOLECULAR_SEED=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH2.09/k05/output \
    EXHALE_MOLECULAR_SEED_X2=local \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > seed.log 2>&1
```

## The commands

```bash
# the molecular seed, built by the binary from the certified atomic state
# with "Load IC? True" in input.inp, which this case states
OMP_NUM_THREADS=8 EXHALE_MOLECULAR_SEED=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH2.09/k05/output \
    EXHALE_MOLECULAR_SEED_X2=local \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > seed.log 2>&1

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 \
    EXHALE_OUTER_PASSES=40 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > run.log 2>&1
```

## What came out

The outer passes of the marching plus the JFNK hand-off, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 2 | 1.35E-01 | 1.52E-04 | 2.92E-01 | 1.00E+00 of 1.0E-05 (carrier balance H2) | 412 | 10087.24 |
| 2 | 0 | 3.27E-08 | 8.90E-13 | 2.33E-08 | 7.34E-01 of 1.0E-05 (carrier balance H2) | 235 | 14286.47 |
| 3 | 0 | 3.89E-08 | 6.48E-13 | 2.91E-08 | 2.97E-01 of 1.0E-05 (carrier balance H2) | 217 | 1626.95 |
| 4 | 2 | 3.03E-08 | 7.42E-13 | 5.24E-08 | 2.85E-01 of 1.0E-05 (carrier balance H2) | 217 | 711.13 |
| 5 | 2 | 3.27E-08 | 1.94E-12 | 3.48E-07 | 1.59E-01 of 1.0E-05 (carrier balance H2) | 217 | 910.23 |
| 6 | 0 | 3.20E-08 | 8.90E-13 | 1.28E-07 | 1.34E-01 of 1.0E-05 (carrier balance H2) | 235 | 2014.54 |
| 7 | 2 | 2.67E-08 | 8.06E-13 | 3.94E-08 | 2.01E-01 of 1.0E-05 (carrier balance H2) | 233 | 1176.69 |
| 8 | 2 | 2.73E-08 | 4.93E-12 | 9.34E-07 | 2.97E-01 of 1.0E-05 (carrier balance H2) | 225 | 668.86 |
| 9 | 0 | 2.20E-08 | 1.37E-12 | 1.81E-08 | 3.63E-01 of 1.0E-05 (carrier balance H2) | 224 | 919.96 |
| 10 | 2 | 2.40E-08 | 1.66E-12 | 2.91E-07 | 3.87E-01 of 1.0E-05 (carrier balance H2) | 226 | 1464.61 |
| 11 | 2 | 1.91E-08 | 1.43E-12 | 9.55E-08 | 3.52E-01 of 1.0E-05 (carrier balance H2) | 234 | 639.37 |
| 12 | 2 | 2.43E-08 | 9.72E-13 | 1.74E-07 | 2.85E-01 of 1.0E-05 (carrier balance H2) | 242 | 2150.43 |
| 13 | 0 | 2.54E-08 | 1.26E-12 | 7.28E-08 | 2.37E-01 of 1.0E-05 (carrier balance H2) | 219 | 1565.87 |
| 14 | 2 | 2.86E-08 | 8.33E-13 | 1.71E-07 | 1.99E-01 of 1.0E-05 (carrier balance H2) | 222 | 1099.66 |
| 15 | 2 | 2.66E-08 | 1.10E-12 | 1.19E-07 | 1.69E-01 of 1.0E-05 (carrier balance H2) | 237 | 1971.94 |
| 16 | 2 | 2.51E-08 | 1.31E-12 | 3.99E-08 | 2.18E-01 of 1.0E-05 (carrier balance H2) | 245 | 1673.62 |
| 17 | 2 | 1.07E+00 | 2.48E-04 | 4.50E-01 | 2.43E-01 of 1.0E-05 (carrier balance H2) | 245 | 11890.87 |
| 18 | 0 | 1.94E-08 | 1.52E-12 | 1.35E-07 | 3.43E-01 of 1.0E-05 (carrier balance H2) | 242 | 4283.38 |
| 19 | 0 | 1.62E-08 | 1.24E-12 | 1.10E-07 | 3.84E-01 of 1.0E-05 (carrier balance H2) | 242 | 1613.65 |
| 20 | 0 | 2.06E-08 | 1.29E-12 | 9.99E-08 | 4.35E-01 of 1.0E-05 (carrier balance H2) | 238 | 2523.65 |
| 21 | 2 | 1.46E-08 | 1.88E-12 | 9.49E-08 | 5.19E-01 of 1.0E-05 (carrier balance H2) | 234 | 8233.59 |

- solver verdict: **info = 2** (last `||R||` = 4.951E-07), by the marching plus the JFNK hand-off
- certification of the state written: **NOT CERTIFIED**
  - hydrodynamic mass row: row measure  3.945E-09 above  1.6E-09 at cell 25
  - carrier balance H2: gated row measure  5.194E-01 above  1.0E-05 at cell 234 (a wind cell)
  - elemental transport He/H partition: gated row measure  6.614E-05 above  1.0E-05 at cell 228 (a wind cell)
  - the bound was refused at cell 1, r   1.0002: relative change of the cell particle count  1.00E-02, from an entry n_tot + n_e of  1.75E-01
- wall clock **1375m18s** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-16 11:09:38 to 2026-09-17 10:04:56

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

## Re-measuring this state, and how closely it comes back

The two files of `output/` can be handed back to the binary as `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with `Restart intent: stationary evaluate`, which measures the state as it stands, writes it back unchanged and exits 0 only if every active equation is within its tolerance.

That is the post-processing pass of this record: it is the route the advection-corrected profiles and the mass-loss line above were produced on. Three states are named in it and every product says which one it describes: the LOADED state, the conserved variables the two files carry, which the pass does not write to; the WORK state, which holds those conserved variables fixed and derives the pressure and the temperature from the composition one equilibrium sweep returns, through the caloric equation of state, and which is what `Hydro_ioniz.txt`, `Ion_species.txt`, the breakdowns and the `certified=` pair of their header describe; and the ADVECTION-DERIVED composition the post-process builds from the work state, which is what `Hydro_ioniz_adv.txt` and `Ion_species_adv.txt` carry. The derived pair states the work state's certification pair on a `# derived_from:` line, as provenance of the state it was built from, and makes no certification claim about its own rows. One sweep is not a closed chemical and thermal fixed point, so the work state reports its own closure defect in `pp.log` and "refreshed" never means "closed".

Two answers are kept apart in `pp.log` and in the list above: whether the stationary claim the solved state was written with reproduces when that state is handed back, and what verdict the work state gets. A passing work state is not a confirmation of a claim that did not reproduce.

What comes back exactly and what does not, MEASURED over the 120 certified states of `models/` (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8). The conserved state round-trips to the last bit: the radius, velocity and pressure columns are read as written, and since item L18 the mass density is read from its own column rather than rebuilt from the species, so the state re-measured is the state certified. The ROW MEASURES do not round-trip bitwise and cannot: the first equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at its rounding floor, where one unit in the last place of the density moves the cell-wise maximum of the energy row by about 11 per cent. Re-evaluated against in-run, the mass row comes back within a factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to 1.93 (median 1.22); the momentum row sits four to six decades below its tolerance and its ratio is the floor itself. A state certified at more than about half of its tolerance may therefore be refused on re-evaluation, and that refusal is a property of the cancellation and not of this file.
