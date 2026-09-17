# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_photochem_gj1132_kzzprofile/HeH12`

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
| run title | `Simulation for LHS1140b_atomic_photochem_gj1132_kzzprofile_HeH12` |

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
| 1 | 0 | 1.14E-09 | 2.18E-14 | 1.76E-08 | 2.28E-04 of 1.0E-05 (elemental transport He/H partition) | 269 | 49.16 |
| 2 | 0 | 7.38E-10 | 1.37E-14 | 1.95E-08 | 1.12E-04 of 1.0E-05 (elemental transport He/H partition) | 271 | 68.98 |
| 3 | 0 | 1.24E-09 | 4.04E-14 | 3.05E-08 | 5.49E-05 of 1.0E-05 (elemental transport He/H partition) | 274 | 49.85 |
| 4 | 0 | 8.75E-10 | 1.97E-14 | 2.41E-08 | 2.71E-05 of 1.0E-05 (elemental transport He/H partition) | 276 | 49.80 |
| 5 | 0 | 8.03E-10 | 1.64E-14 | 1.79E-08 | 1.34E-05 of 1.0E-05 (elemental transport He/H partition) | 278 | 48.12 |
| 6 | 0 | 1.51E-09 | 1.15E-13 | 5.63E-08 | 6.65E-06 of 1.0E-05 (elemental transport He/H partition) | 281 | 31.97 |

- solver verdict: **info = 0** (last `||R||` = 5.705E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:25:16 to 2026-09-14 10:39:33

### Iteration `k01`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 1.00E-09 | 1.35E-12 | 1.70E-08 | 6.65E-06 of 1.0E-05 (elemental transport He/H partition) | 281 | 167.45 |

- solver verdict: **info = 0** (last `||R||` = 1.699E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:25:16 to 2026-09-14 10:39:33

### Iteration `k02`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 6.74E-10 | 7.01E-13 | 1.59E-08 | 6.65E-06 of 1.0E-05 (elemental transport He/H partition) | 281 | 137.18 |

- solver verdict: **info = 0** (last `||R||` = 1.592E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:25:16 to 2026-09-14 10:39:33

### Iteration `k03`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 9.17E-10 | 2.09E-12 | 1.93E-08 | 6.65E-06 of 1.0E-05 (elemental transport He/H partition) | 281 | 90.21 |

- solver verdict: **info = 0** (last `||R||` = 1.931E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:25:16 to 2026-09-14 10:39:33

### Iteration `k04`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 1.10E-09 | 1.58E-12 | 1.94E-08 | 6.65E-06 of 1.0E-05 (elemental transport He/H partition) | 281 | 27.57 |

- solver verdict: **info = 0** (last `||R||` = 1.942E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- He I 10830 red-pair equivalent width = **2.2729** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **7.1717** %, FWHM = 0.2966 A (three-Gaussian fit)
- wall clock **14m17s** at `OMP_NUM_THREADS=6`, 2026-09-14 10:25:16 to 2026-09-14 10:39:33

### The rung

`closure_history.txt`, last row: the converged trial fluxes, the reservoir He/H at the matching level and log10 Mdot.

```
   4   5.4282626213E+06   7.1933689548E+07   5.4722752785E+06   7.5367012900E+07  8.04284E-03  4.55547E-02   0.5000 overlap   3.63705E-03  4.79417E-03  1.2013647E+01    7.910 4b5b1987e57612b9b060d4feb6acd87906683b5f53af16b82fd08c0d9c946d98 0
```

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=6`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
