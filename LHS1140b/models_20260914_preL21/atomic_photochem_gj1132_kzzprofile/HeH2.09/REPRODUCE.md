# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_photochem_gj1132_kzzprofile/HeH2.09`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: the Photochem column handed over as `Lower atmosphere profile:` (`photochem`)
- spectrum: `gj1132`
- mixing: binary H/He element diffusion with K_zz(p) from the profile (`kzzprofile`)

The rung has no `input.inp` of its own: `input_template.inp` beside `closure.json` is what each iteration starts from, and the closure driver writes that iteration's own `input.inp` into `kNN/` from it, adding the route's keys (`input_keys`) and the profile the chemistry step just produced.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `b9419ba7cd8e519ca86896b29c382762` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 6; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_photochem_gj1132_kzzprofile_HeH2.09` |

## The seed

None: this case starts from the code's own initial condition (`Load IC? False`). The archive holds no state of its physics.

## The commands

```bash
# the rung: the Photochem column and the wind alternate until the elemental
# flux at the microbar match is its own fixed point.  closure.json states the
# route each iteration's wind is solved by (input_keys) and its environment.
sed -e 's/"omp_num_threads": *[0-9]*/"omp_num_threads": 6/' \
    closure.json > closure_run.json
/opt/miniconda3/bin/python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/element_flux_closure.py . \
    --config closure_run.json \
    --phi0-H 4.7682069882E+06 --phi0-He 2.0433648550E+07 \
    --tol 0.05 --kmax 8 --resume > closure_stdout.log 2>&1

# the transit spectrum of the iterate the rung converged on
. /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/winered_hires_y.sh
cd $(ls -d k[0-9][0-9] | sort | tail -n 1)
MPLBACKEND=Agg PYTHONPATH=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00 python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_transit.py > transit.log 2>&1
```

## What came out

### Iteration `k00`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 1.92E-09 | 3.95E-13 | 1.20E-07 | 2.67E-03 of 1.0E-05 (elemental transport He/H partition) | 266 | 296.70 |
| 2 | 0 | 1.09E-09 | 2.69E-13 | 2.10E-08 | 1.25E-03 of 1.0E-05 (elemental transport He/H partition) | 268 | 191.41 |
| 3 | 0 | 1.01E-09 | 5.97E-13 | 2.26E-08 | 6.00E-04 of 1.0E-05 (elemental transport He/H partition) | 272 | 166.27 |
| 4 | 0 | 1.02E-09 | 1.58E-13 | 2.77E-08 | 2.93E-04 of 1.0E-05 (elemental transport He/H partition) | 276 | 164.96 |
| 5 | 0 | 9.76E-10 | 5.76E-14 | 1.62E-08 | 1.46E-04 of 1.0E-05 (elemental transport He/H partition) | 281 | 137.44 |
| 6 | 0 | 1.08E-09 | 6.12E-14 | 1.45E-08 | 7.44E-05 of 1.0E-05 (elemental transport He/H partition) | 286 | 165.91 |
| 7 | 0 | 1.01E-09 | 1.69E-12 | 2.22E-08 | 3.84E-05 of 1.0E-05 (elemental transport He/H partition) | 290 | 130.23 |
| 8 | 0 | 5.99E-10 | 7.57E-13 | 1.86E-08 | 2.00E-05 of 1.0E-05 (elemental transport He/H partition) | 294 | 130.63 |
| 9 | 0 | 5.16E-10 | 2.95E-13 | 1.58E-08 | 1.05E-05 of 1.0E-05 (elemental transport He/H partition) | 296 | 128.68 |
| 10 | 0 | 7.62E-10 | 1.77E-13 | 1.85E-08 | 5.45E-06 of 1.0E-05 (elemental transport He/H partition) | 298 | 128.51 |

- solver verdict: **info = 0** (last `||R||` = 1.849E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.92**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:25:16 to 2026-09-14 10:50:05

### Iteration `k01`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 7.00E-10 | 7.54E-13 | 1.42E-08 | 5.46E-06 of 1.0E-05 (elemental transport He/H partition) | 298 | 185.78 |

- solver verdict: **info = 0** (last `||R||` = 1.424E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.92**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:25:16 to 2026-09-14 10:50:05

### Iteration `k02`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 6.31E-10 | 3.56E-13 | 1.94E-08 | 5.46E-06 of 1.0E-05 (elemental transport He/H partition) | 298 | 143.89 |

- solver verdict: **info = 0** (last `||R||` = 1.943E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.92**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:25:16 to 2026-09-14 10:50:05

### Iteration `k03`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 8.29E-10 | 1.94E-13 | 1.81E-08 | 5.46E-06 of 1.0E-05 (elemental transport He/H partition) | 298 | 134.73 |

- solver verdict: **info = 0** (last `||R||` = 1.806E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.92**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:25:16 to 2026-09-14 10:50:05

### Iteration `k04`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 1.36E-09 | 1.64E-12 | 2.76E-08 | 5.46E-06 of 1.0E-05 (elemental transport He/H partition) | 298 | 132.02 |

- solver verdict: **info = 0** (last `||R||` = 2.762E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.92**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:25:16 to 2026-09-14 10:50:05

### Iteration `k05`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 1.08E-09 | 7.81E-13 | 2.20E-08 | 5.46E-06 of 1.0E-05 (elemental transport He/H partition) | 298 | 133.02 |

- solver verdict: **info = 0** (last `||R||` = 2.202E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.92**, from the post-processing pass
- He I 10830 red-pair equivalent width = **1.5006** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **5.4077** %, FWHM = 0.2603 A (three-Gaussian fit)
- wall clock **24m49s** at `OMP_NUM_THREADS=6`, 2026-09-14 10:25:16 to 2026-09-14 10:50:05

### The rung

`closure_history.txt`, last row: the converged trial fluxes, the reservoir He/H at the matching level and log10 Mdot.

```
   5   2.5058892177E+07   5.5588453226E+07   2.5713195869E+07   5.6722718267E+07  2.54462E-02  1.99967E-02   0.5000 overlap   1.81912E-03  5.75299E-03  2.0923771E+00    7.920 d90bfac9a413565fdc6cb6107585da3c4bca4e526cd4d913f0bda08c59b18806 0
```

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=6`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
