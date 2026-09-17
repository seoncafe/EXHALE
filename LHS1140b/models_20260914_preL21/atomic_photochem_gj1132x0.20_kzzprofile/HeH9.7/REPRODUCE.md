# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_photochem_gj1132x0.20_kzzprofile/HeH9.7`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: the Photochem column handed over as `Lower atmosphere profile:` (`photochem`)
- spectrum: `gj1132x0.20`
- mixing: binary H/He element diffusion with K_zz(p) from the profile (`kzzprofile`)

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `b9419ba7cd8e519ca86896b29c382762` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_photochem_gj1132x0.20_kzzprofile_HeH9.7` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/xuvkzz/C0p20/output  tier4:archive/refresh_j96  HeH=9.7  target=9.71103  dlog10=0.0005  candidates=556  (an archived state); then a continuation from the first solve's written state at dtau0=1.0e8 (its log run_dtau0_first.log)
```

The fields are the state directory, the tier it was taken from (`tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/xuvkzz/C0p20/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/xuvkzz/C0p20/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > run.log 2>&1

# the first solve returned info=1: the written state reloaded and solved
# again with the pseudo-time start raised (L4e)
cp output/Hydro_ioniz.txt output/Hydro_ioniz_IC.txt
cp output/Ion_species.txt output/Ion_species_IC.txt
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e8 \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > run.log 2>&1

# the advection-corrected profiles: the solved state becomes the IC and one
# time step is taken at a CFL small enough to leave it where it is
for f in Hydro_ioniz Ion_species; do cp -f output/$f.txt output/${f}_IC.txt; done
cp -f input.inp input.inp.solved
sed -i -e 's/^Load IC?.*/Load IC? True/' -e 's/^Do only PP:.*/Do only PP: True/' \
    -e '/^Restart intent:/d' -e '/^Solver:/d' -e '/^du_th /d' input.inp
echo 'CFL: 1.0e-12' >> input.inp
OMP_NUM_THREADS=8 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x > pp.log 2>&1
mv -f input.inp.solved input.inp

# the transit spectrum, through the WINERED HIRES-Y kernel of the measurement
. /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/winered_hires_y.sh
MPLBACKEND=Agg PYTHONPATH=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00 python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE_transit.py > transit.log 2>&1
```

## What came out

The outer passes of the partitioned stationary route, as `run.log` records them:

| pass | hydro info | mass | momentum | energy | worst gated species row | at cell | s |
|---|---|---|---|---|---|---|---|
| 1 | 0 | 7.43E-09 | 2.77E-14 | 3.29E-07 | 2.10E-02 of 1.0E-05 (elemental transport He/H partition) | 263 | 131.34 |
| 2 | 0 | 8.47E-09 | 6.23E-14 | 2.26E-07 | 1.27E-02 of 1.0E-05 (elemental transport He/H partition) | 265 | 149.34 |
| 3 | 0 | 3.52E-09 | 2.71E-14 | 8.76E-08 | 6.99E-03 of 1.0E-05 (elemental transport He/H partition) | 267 | 110.53 |
| 4 | 0 | 2.90E-09 | 3.26E-14 | 1.16E-07 | 3.60E-03 of 1.0E-05 (elemental transport He/H partition) | 268 | 82.92 |
| 5 | 0 | 1.59E-08 | 2.65E-13 | 4.24E-07 | 1.78E-03 of 1.0E-05 (elemental transport He/H partition) | 268 | 89.96 |
| 6 | 0 | 1.12E-08 | 1.18E-13 | 4.22E-07 | 8.67E-04 of 1.0E-05 (elemental transport He/H partition) | 269 | 57.65 |
| 7 | 0 | 5.10E-09 | 4.64E-14 | 1.43E-07 | 4.19E-04 of 1.0E-05 (elemental transport He/H partition) | 269 | 65.47 |
| 8 | 0 | 6.33E-09 | 5.46E-14 | 2.33E-07 | 2.02E-04 of 1.0E-05 (elemental transport He/H partition) | 269 | 48.04 |
| 9 | 0 | 7.48E-09 | 4.77E-14 | 1.53E-07 | 9.75E-05 of 1.0E-05 (elemental transport He/H partition) | 269 | 40.46 |
| 10 | 0 | 3.24E-09 | 3.88E-14 | 1.28E-07 | 4.70E-05 of 1.0E-05 (elemental transport He/H partition) | 270 | 38.97 |
| 11 | 0 | 3.50E-09 | 2.51E-14 | 8.17E-08 | 2.28E-05 of 1.0E-05 (elemental transport He/H partition) | 270 | 41.43 |
| 12 | 0 | 1.45E-08 | 1.32E-13 | 3.40E-07 | 1.10E-05 of 1.0E-05 (elemental transport He/H partition) | 270 | 27.66 |
| 13 | 0 | 2.19E-08 | 1.72E-13 | 6.44E-07 | 5.33E-06 of 1.0E-05 (elemental transport He/H partition) | 270 | 36.36 |

- solver verdict: **info = 0** (last `||R||` = 6.443E-07), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.18**, from the post-processing pass
- He I 10830 red-pair equivalent width = **0.6990** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **2.5319** %, FWHM = 0.2578 A (three-Gaussian fit)
- wall clock **75m11s** at `OMP_NUM_THREADS=8`, 2026-09-14 11:44:35 to 2026-09-14 12:59:46

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
