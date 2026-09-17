# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_photochem_gj1132_kzzprofile/HeH8`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: the Photochem column handed over as `Lower atmosphere profile:` (`photochem`)
- spectrum: `gj1132`
- mixing: binary H/He element diffusion with K_zz(p) from the profile (`kzzprofile`)

The rung has no `input.inp` of its own: `input_template.inp` beside `closure.json` is what each iteration starts from, and the closure driver writes that iteration's own `input.inp` into `kNN/` from it, adding the route's keys (`input_keys`) and the profile the chemistry step just produced.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `c2e9c9990b9f14f1be8cd77abca68945` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_photochem_gj1132_kzzprofile_HeH8` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models_20260915_db87/atomic_photochem_gj1132_kzzprofile/HeH8/k04/output  tier0:models_20260915_db87/atomic_photochem_gj1132_kzzprofile/HeH8/k04  HeH=8  target=8  dlog10=0.0000  candidates=1615  (this same case's own most recent certified state)
```

The fields are the state directory, the tier it was taken from (`tier0` this same case's own most recent certified state -- the one this tree carries, else the one a preserved tree such as `models_20260914_preL21/` holds, which is what a re-solve of the whole catalog continues from -- `tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models_20260915_db87/atomic_photochem_gj1132_kzzprofile/HeH8/k04/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt seed
```

## The commands

```bash
# the seed of iteration 0, on the current cell centers
mkdir -p seed
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models_20260915_db87/atomic_photochem_gj1132_kzzprofile/HeH8/k04/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt seed

# the rung: the Photochem column and the wind alternate until the elemental
# flux at the microbar match is its own fixed point.  closure.json states the
# route each iteration's wind is solved by (input_keys) and its environment.
sed -e 's/"omp_num_threads": *[0-9]*/"omp_num_threads": 8/' \
    closure.json > closure_run.json
/opt/miniconda3/bin/python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/element_flux_closure.py . \
    --config closure_run.json \
    --phi0-H 4.7682069882E+06 --phi0-He 2.0433648550E+07 \
    --tol 0.05 --kmax 8 --seed /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH8/seed > closure_stdout.log 2>&1

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
| 1 | 0 | 4.59E-09 | 2.28E-13 | 2.42E-08 | 6.56E-06 of 1.0E-05 (elemental transport He/H partition) | 217 | 101.77 |

- solver verdict: **info = 0** (last `||R||` = 2.416E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.90**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-16 11:04:55 to 2026-09-16 11:19:35

### Iteration `k01`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 3.75E-09 | 1.50E-13 | 1.70E-08 | 6.56E-06 of 1.0E-05 (elemental transport He/H partition) | 217 | 102.62 |

- solver verdict: **info = 0** (last `||R||` = 1.696E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.90**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-16 11:04:55 to 2026-09-16 11:19:35

### Iteration `k02`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 4.52E-09 | 2.02E-12 | 2.61E-08 | 6.56E-06 of 1.0E-05 (elemental transport He/H partition) | 217 | 96.42 |

- solver verdict: **info = 0** (last `||R||` = 2.611E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.90**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-16 11:04:55 to 2026-09-16 11:19:35

### Iteration `k03`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 1.71E-09 | 6.42E-14 | 1.52E-08 | 6.56E-06 of 1.0E-05 (elemental transport He/H partition) | 217 | 99.71 |

- solver verdict: **info = 0** (last `||R||` = 1.524E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.90**, from the post-processing pass
- wall clock **(see the whole rung)** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-16 11:04:55 to 2026-09-16 11:19:35

### Iteration `k04`

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 2.72E-09 | 5.10E-13 | 1.24E-08 | 6.56E-06 of 1.0E-05 (elemental transport He/H partition) | 217 | 46.80 |

- solver verdict: **info = 0** (last `||R||` = 1.240E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.90**, from the post-processing pass
- He I 10830 red-pair equivalent width = **2.2004** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **7.1146** %, FWHM = 0.2896 A (three-Gaussian fit)
- wall clock **14m40s** at `OMP_NUM_THREADS=8` on `lart3`, 2026-09-16 11:04:55 to 2026-09-16 11:19:35

### The rung

`closure_history.txt`, last row: the converged trial fluxes, the reservoir He/H at the matching level and log10 Mdot.

```
   4   7.9036481653E+06   6.7824097414E+07   8.1126262030E+06   7.0983497005E+07  2.57596E-02  4.45089E-02   0.5000 steady    2.09158E-03  5.03001E-04  8.0091064E+00    7.900 b71a7cfb0161775f369c58aa82f92f2df0f61feba816f5b31ae2b1f148047058 0
```

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.

## Re-measuring this state, and how closely it comes back

The two files of `output/` can be handed back to the binary as `Hydro_ioniz_IC.txt` and `Ion_species_IC.txt` with `Restart intent: stationary evaluate`, which measures the state as it stands, writes it back unchanged and exits 0 only if every active equation is within its tolerance.

What comes back exactly and what does not, MEASURED over the 120 certified states of `models/` (`docs/lhs1140b_stationary_L18_20260915.md`, sections 6 and 8). The conserved state round-trips to the last bit: the radius, velocity and pressure columns are read as written, and since item L18 the mass density is read from its own column rather than rebuilt from the species, so the state re-measured is the state certified. The ROW MEASURES do not round-trip bitwise and cannot: the first equilibrium sweep of the re-entry moves the composition by 1e-14 to 1e-11, and the flux assembly of a subsonic base is at its rounding floor, where one unit in the last place of the density moves the cell-wise maximum of the energy row by about 11 per cent. Re-evaluated against in-run, the mass row comes back within a factor 0.78 to 1.71 (median 1.12) and the energy row within 0.89 to 1.93 (median 1.22); the momentum row sits four to six decades below its tolerance and its ratio is the floor itself. A state certified at more than about half of its tolerance may therefore be refused on re-evaluation, and that refusal is a property of the cancellation and not of this file.
