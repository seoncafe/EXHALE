# Partitioned transport–wind experiment archive

Date: September 11, 2026

The [experiment report](../../solver_partition_experiment_20260911.md) explains the measurements, source findings, and limits. No measured state passed the complete selected stationary certification contract. No production source or existing simulation product was changed.

## Files and reproduction

- `transport_wind_experiment.f90`: diagnostic driver linked to production modules, excluding the production main program.
- `build_experiment.sh`: creates a private build. It uses the project Makefile and the OpenBLAS installation under `/opt/miniconda3/lib` on this machine.
- `run_case.sh`: copies one fixture into a new, uniquely named directory, then runs the diagnostic with a 180 s wall-time limit.
- `summarize_runs.py`: prints measurements from existing logs and tables without regenerating simulation products.
- `build.log`: retained initial private production-build output.
- `runs/`: copied inputs, initial conditions, measurements, and available provenance records.

From the project root, build and run with:

```bash
bash docs/audit_20260905/partition_experiment_20260911/build_experiment.sh
bash docs/audit_20260905/partition_experiment_20260911/run_case.sh /absolute/path/from/build/transport_wind_experiment.x carrier_elem_newton partition 3 40 0.5 0.01 8
python docs/audit_20260905/partition_experiment_20260911/summarize_runs.py
```

Runner arguments are executable, fixture, mode, passes, Newton iteration cap, element `omega`, carrier `trust`, and OpenMP threads. Modes are `partition`, `coupled`, and `species`. BLAS uses one thread. The large carrier setting `1e6` in the frozen-background experiment disables the movement stopping trigger along that trajectory; it is not a physical time step.

The retained private build is `/tmp/exhale_partition_20260911.ZxXJEB`. Temporary build retention is not guaranteed; use the script to rebuild if necessary. The script builds the source currently checked out, so reproducing the archived experiment also requires production revision `a0f4292b39e68cb83c1578a63cf4d50bd081242f` and the recorded inputs. The current diagnostic source corresponds to instrumentation version 3. Earlier instrumentation versions did not record every quantity below.

## Completed measurements used in the report

| Directory under `runs/` | Configuration | Outcome |
| --- | --- | --- |
| `carrier_elem_newton_partition.g44qGR` | Version 3, 8 threads, 3 passes, cap 40, trust 0.01 | Completed; not certified |
| `atomic_elem_newton_partition.ouOtBg` | Version 3, 8 threads, 2 passes, cap 40, omega 0.5 | Completed; not certified |
| `carrier_elem_newton_coupled.VaAVVl` | Version 3, 8 threads, cap 40 | Completed; not certified |
| `carrier_elem_newton_partition.TZ326I` | Version 2, 1 thread, 1 pass, cap 40, trust 1e6 | Frozen H2 balance solved accurately; refreshed state not certified |
| `atomic_elem_newton_species.MlxOaF` | Version 3, 1 thread, no hydro update, omega 1 | Raw composition mass closure failed; chemical refresh became nonfinite; exit 1 |
| `atomic_elem_newton_partition.gVmnY3` | Version 3, 8 threads, 1 pass, cap 10, omega 0.5 | Completed; elemental balances degraded; not certified |

## Pilots and incomplete runs

| Directory under `runs/` | Interpretation |
| --- | --- |
| `carrier_elem_newton_partition.fkIQWF` | Completed version-1, one-thread pilot; row table identical to `g44qGR` |
| `atomic_elem_newton_partition.suxeSF` | Version 2, eight threads; two passes completed, third hydro call interrupted by timeout |
| `atomic_elem_newton_partition.KjyB1h` | Version 1, one thread; timeout before first hydro return |
| `atomic_elem_newton_partition.A7QpWW` | Initial-stage-only incomplete pilot; final exit status was not observed; excluded from conclusions |
| `atomic_elem_newton_species.gcGxmR` | Version-1 repetition of the atomic stress case; nonfinite refreshed state |
| `atomic_elem_newton_species.zOWXKb` | Version-2 repetition of the same stress case; adds raw composition snapshot but not raw mass-closure logging |
| `atomic_elem_newton_coupled.vB6B79` | Version 2, one thread; timeout before returned solution |
| `carrier_elem_newton_coupled.Ru63Ch` | Version 2, one thread; timeout before returned solution |
| `atomic_elem_newton_coupled.FK3F9p` | Version 3, eight threads; timeout before returned solution |

The first pilots predate automatic provenance and exit-status logging. Their missing files have not been filled with inferred values. Repetitions of the same stress configuration are not independent physical test cases.

## Interpreting the output

- `run.log` records native frozen-background residuals, solver return flags, warnings, certification, mass closure, and total elapsed time when execution finishes.
- `rows.tsv` records full-state certification entries. `maximum` is the whole-column maximum; `wind`, `reported_region`, and `gated` must not be treated as interchangeable. The `pass` field is the corresponding row's tolerance result, not an outer iteration number.
- `steps.tsv` is whitespace-separated despite its extension. Its columns are stage, operation, info, steps, drift, and seconds. For hydro/coupled operations, `steps` is the requested iteration cap, not the actual count; consult `run.log` for actual iterations. For species operations, `info = 0` is a diagnostic placeholder because the called routine has no success return. It is not a convergence assertion.
- An operating-system exit status of 0 means the diagnostic completed, not that the atmosphere is certified. Status 124 means the external timeout fired. Incomplete solver logs are not returned-state measurements.
- Full-state `*.state` files contain radius, conserved hydro variables, temperature, and species entries, including ghost cells. Raw `frozen_after_*.state` files retain the temperature of the frozen background before the chemical refresh. Species entries are not all probabilities; an entry above one alone is not the mass-closure test.
- Version 3 adds raw composition mass closure using production `get_species_densities` and `calc_rho`, and `*.faces` files containing face index and numerical mass flux for physical faces `0:N`. Face files do not contain radii.
- `provenance.log` and `driver_snapshot.f90`, where present, identify the executable and driver used. Earlier versions remain distinct from the final source rather than having their results regenerated.

The initial state and copied input in the final molecular coupled and partitioned comparisons were identical in direct file comparisons. Some jobs overlapped, so timings are not isolated performance benchmarks. The diagnostic intentionally continues after failed hydro calls to measure composition feedback, whereas the production outer loop exits on a nonzero hydro return. None of these runs establishes physical-time stability, grid convergence, or global elemental and energy budgets.

## What the repository carries

The driver, the two scripts, `summarize_runs.py`, `build.log` and this file are
committed. The `runs/` directory (38 MB of copied inputs, `.state` and `.faces`
dumps and logs) is kept in this working copy only, like `backup/`; the numbers
the report quotes from it are all in the report itself and in
`docs/Update_EXHALE.md` section 8.
