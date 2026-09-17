# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_photochem_gj1132_kzzprofile/HeH5`

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
| run title | `Simulation for LHS1140b_atomic_photochem_gj1132_kzzprofile_HeH5` |

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
| 1 | 0 | 5.94E-10 | 1.78E-14 | 1.10E-08 | 2.70E-04 of 1.0E-05 (elemental transport He/H partition) | 275 | 49.14 |
| 2 | 0 | 1.55E-09 | 9.05E-14 | 3.22E-08 | 1.33E-04 of 1.0E-05 (elemental transport He/H partition) | 279 | 39.21 |
| 3 | 0 | 1.04E-09 | 1.53E-14 | 1.87E-08 | 6.59E-05 of 1.0E-05 (elemental transport He/H partition) | 282 | 43.93 |
| 4 | 0 | 4.46E-09 | 3.58E-14 | 6.72E-08 | 3.29E-05 of 1.0E-05 (elemental transport He/H partition) | 285 | 25.62 |
| 5 | 0 | 7.54E-10 | 2.17E-14 | 1.82E-08 | 1.66E-05 of 1.0E-05 (elemental transport He/H partition) | 289 | 41.87 |
| 6 | 0 | 5.59E-10 | 9.64E-15 | 1.78E-08 | 8.46E-06 of 1.0E-05 (elemental transport He/H partition) | 292 | 37.88 |

- solver verdict: **info = 0** (last `||R||` = 1.779E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:39:33 to 2026-09-14 10:57:08

### Iteration `k01`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 8.66E-10 | 4.52E-13 | 1.77E-08 | 8.46E-06 of 1.0E-05 (elemental transport He/H partition) | 292 | 137.67 |

- solver verdict: **info = 0** (last `||R||` = 1.775E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:39:33 to 2026-09-14 10:57:08

### Iteration `k02`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 9.83E-10 | 2.23E-13 | 2.72E-08 | 8.46E-06 of 1.0E-05 (elemental transport He/H partition) | 292 | 142.10 |

- solver verdict: **info = 0** (last `||R||` = 2.717E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:39:33 to 2026-09-14 10:57:08

### Iteration `k03`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 1.69E-09 | 1.82E-12 | 3.51E-08 | 8.46E-06 of 1.0E-05 (elemental transport He/H partition) | 292 | 131.15 |

- solver verdict: **info = 0** (last `||R||` = 3.512E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=6`, 2026-09-14 10:39:33 to 2026-09-14 10:57:08

### Iteration `k04`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 9.92E-10 | 9.34E-13 | 3.09E-08 | 8.46E-06 of 1.0E-05 (elemental transport He/H partition) | 292 | 132.17 |

- solver verdict: **info = 0** (last `||R||` = 3.093E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.91**, from the post-processing pass
- He I 10830 red-pair equivalent width = **2.1099** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **7.0552** %, FWHM = 0.2805 A (three-Gaussian fit)
- wall clock **17m35s** at `OMP_NUM_THREADS=6`, 2026-09-14 10:39:33 to 2026-09-14 10:57:08

### The rung

`closure_history.txt`, last row: the converged trial fluxes, the reservoir He/H at the matching level and log10 Mdot.

```
   4   1.1973128667E+07   6.6050379837E+07   1.2453313397E+07   6.9091611146E+07  3.85588E-02  4.40174E-02   0.5000 overlap   2.53029E-03  5.10623E-03  5.0056766E+00    7.910 45a07df2565c0db36d79c026101ed2f863646877069f6b1cc77a173bf08e5099 0
```

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=6`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
