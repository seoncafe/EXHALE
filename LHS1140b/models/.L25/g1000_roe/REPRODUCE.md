# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`.L25/g1000_roe`

- numerical flux: `ROE`

The HLLC flux does not reach the stationary root on this grid at this XUV level (L25 step 3), and the same case solved with HLLC on a grid of twice the cells reaches the state the Roe flux reaches on the catalog grid, so the difference is resolution and not a different wind. The Roe flux systematic against HLLC where both solve is Mdot and the He 10830 equivalent width lower by 0.3 to 0.6 per cent, the first cell colder by 6 to 7 per cent and denser by 6 to 8 per cent, and the mass-flux spread unchanged (L25 memo). That systematic is carried wherever the numbers of this case are quoted.

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs two of its keys changed.

## Why this run was made the way it was

L25 step 3 item 3, grid convergence: the seed was mapped onto this grid beforehand and NOSEED=1 was set, because run_case.sh maps onto the 500-cell grid of models/current_grid_Hydro_ioniz.txt only.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `c2e9c9990b9f14f1be8cd77abca68945` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_scalar_gj1132x0.10_kzz1e9_HeH9.7` |

## The seed

None: this case starts from the code's own initial condition (`Load IC? False`). The archive holds no state of its physics.

## The commands

```bash
# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e8 \
    EXHALE_OUTER_PASSES=40 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > run.log 2>&1
```

## What came out

The outer passes of the marching plus the JFNK hand-off, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 2 | 1.69E+00 | 1.59E-03 | 1.88E+00 | 3.71E-03 of 1.0E-05 (elemental transport He/H partition) | 621 | 19357.76 |
| 2 | 2 | 1.97E+00 | 3.02E-04 | 1.90E+00 | 5.74E-03 of 1.0E-05 (elemental transport He/H partition) | 634 | 2888.55 |
| 3 | 2 | 1.55E+00 | 2.77E-04 | 1.92E+00 | 7.20E-03 of 1.0E-05 (elemental transport He/H partition) | 627 | 3301.14 |
| 4 | 2 | 1.59E+00 | 2.16E-04 | 1.89E+00 | 7.90E-03 of 1.0E-05 (elemental transport He/H partition) | 627 | 3146.79 |
| 5 | 2 | 1.59E+00 | 2.16E-04 | 1.88E+00 | 8.15E-03 of 1.0E-05 (elemental transport He/H partition) | 622 | 2010.84 |
| 6 | 2 | 1.59E+00 | 2.16E-04 | 1.88E+00 | 8.23E-03 of 1.0E-05 (elemental transport He/H partition) | 622 | 1989.59 |
| 7 | 2 | 1.59E+00 | 2.16E-04 | 1.88E+00 | 8.18E-03 of 1.0E-05 (elemental transport He/H partition) | 622 | 1989.97 |

- solver verdict: **info = 2** (last `||R||` = 1.880E+00), by the marching plus the JFNK hand-off
- certification of the state written: **NOT CERTIFIED**
  - hydrodynamic mass row: row measure  1.481E+00 above  2.2E-09 at cell 211
  - hydrodynamic momentum row: row measure  2.155E-04 above  1.0E-08 at cell 993
  - hydrodynamic energy row: row measure  1.880E+00 above  1.0E-06 at cell 236
  - elemental transport He/H partition: gated row measure  8.183E-03 above  1.0E-05 at cell 622 (a wind cell)
- wall clock **606m48s** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-16 23:58:08 to 2026-09-17 10:04:56

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

## Re-measuring this state, and how closely it comes back

The two files of `output/` can be handed back to the binary as `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with `Restart intent: stationary evaluate`, which measures the state as it stands, writes it back unchanged and exits 0 only if every active equation is within its tolerance.

That is the post-processing pass of this record: it is the route the advection-corrected profiles and the mass-loss line above were produced on. Three states are named in it and every product says which one it describes: the LOADED state, the conserved variables the two files carry, which the pass does not write to; the WORK state, which holds those conserved variables fixed and derives the pressure and the temperature from the composition one equilibrium sweep returns, through the caloric equation of state, and which is what `Hydro_ioniz.txt`, `Ion_species.txt`, the breakdowns and the `certified=` pair of their header describe; and the ADVECTION-DERIVED composition the post-process builds from the work state, which is what `Hydro_ioniz_adv.txt` and `Ion_species_adv.txt` carry. The derived pair states the work state's certification pair on a `# derived_from:` line, as provenance of the state it was built from, and makes no certification claim about its own rows. One sweep is not a closed chemical and thermal fixed point, so the work state reports its own closure defect in `pp.log` and "refreshed" never means "closed".

Two answers are kept apart in `pp.log` and in the list above: whether the stationary claim the solved state was written with reproduces when that state is handed back, and what verdict the work state gets. A passing work state is not a confirmation of a claim that did not reproduce.

What comes back exactly and what does not, MEASURED over the 120 certified states of `models/` (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8). The conserved state round-trips to the last bit: the radius, velocity and pressure columns are read as written, and since item L18 the mass density is read from its own column rather than rebuilt from the species, so the state re-measured is the state certified. The ROW MEASURES do not round-trip bitwise and cannot: the first equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at its rounding floor, where one unit in the last place of the density moves the cell-wise maximum of the energy row by about 11 per cent. Re-evaluated against in-run, the mass row comes back within a factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to 1.93 (median 1.22); the momentum row sits four to six decades below its tolerance and its ratio is the floor itself. A state certified at more than about half of its tolerance may therefore be refused on re-evaluation, and that refusal is a property of the cancellation and not of this file.
