# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_scalarCNO_gj1132_kzz1e9/HeH2.13`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: the same scalar base plus the C, N and O reservoirs of the photochemical column (base.inp) (`scalarCNO`)
- spectrum: `gj1132`
- mixing: binary H/He element diffusion with He_Kzz = 1e9 cm^2/s (`kzz1e9`)

The `input.inp` in this directory is the file the solution was started from; the runner puts it back after the post-processing pass, which needs three of its keys changed.

## The binary

| | |
|---|---|
| path | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` |
| md5 | `b261c3287e64149a85d0d2dfd2ee2c74` |
| repository HEAD | `43bc28cef58772560bae019559f3592335922212` |
| working tree | dirty (uncommitted changes present) |
| compiler | GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20); GCC: (conda-forge gcc 16.2.0-4) 16.2.0 |
| linear algebra | LAPACK library: OpenBLAS; BLAS threads set to 1 thread(s) (library default was 8; state OPENBLAS_NUM_THREADS to choose) |
| run title | `Simulation for LHS1140b_atomic_scalarCNO_gj1132_kzz1e9_HeH2.13` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/basemetals/d_cno/output  tier4:archive/refresh_j96  HeH=2.13  target=2.13  dlog10=0.0000  candidates=2  (an archived state)
```

The fields are the state directory, the tier it was taken from (`tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/basemetals/d_cno/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/archive_20260830/exhale/refresh_j96/basemetals/d_cno/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic 

# the wind, from the input.inp of this directory
OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 \
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
| 1 | 0 | 1.18E-09 | 1.60E-14 | 2.54E-08 | 2.71E-03 of 1.0E-05 (elemental transport He/H partition) | 265 | 275.61 |
| 2 | 0 | 5.10E-10 | 1.69E-13 | 2.28E-08 | 1.27E-03 of 1.0E-05 (elemental transport He/H partition) | 267 | 188.16 |
| 3 | 0 | 7.90E-10 | 5.02E-13 | 3.06E-08 | 6.07E-04 of 1.0E-05 (elemental transport He/H partition) | 271 | 159.70 |
| 4 | 0 | 9.74E-10 | 2.54E-13 | 1.91E-08 | 2.95E-04 of 1.0E-05 (elemental transport He/H partition) | 276 | 161.53 |
| 5 | 0 | 5.70E-10 | 1.21E-13 | 1.60E-08 | 1.46E-04 of 1.0E-05 (elemental transport He/H partition) | 281 | 161.69 |
| 6 | 0 | 9.76E-10 | 1.02E-13 | 1.96E-08 | 7.46E-05 of 1.0E-05 (elemental transport He/H partition) | 286 | 159.29 |
| 7 | 0 | 1.25E-09 | 2.62E-12 | 2.66E-08 | 3.88E-05 of 1.0E-05 (elemental transport He/H partition) | 291 | 127.12 |
| 8 | 0 | 9.62E-10 | 1.25E-12 | 2.14E-08 | 2.04E-05 of 1.0E-05 (elemental transport He/H partition) | 295 | 126.22 |
| 9 | 0 | 9.64E-10 | 4.80E-13 | 2.81E-08 | 1.08E-05 of 1.0E-05 (elemental transport He/H partition) | 298 | 125.82 |
| 10 | 0 | 1.72E-09 | 2.24E-13 | 3.19E-08 | 5.68E-06 of 1.0E-05 (elemental transport He/H partition) | 300 | 125.40 |

- solver verdict: **info = 0** (last `||R||` = 3.193E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.87**, from the post-processing pass
- He I 10830 red-pair equivalent width = **1.3917** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **5.0182** %, FWHM = 0.2603 A (three-Gaussian fit)
- wall clock **27m34s** at `OMP_NUM_THREADS=8`, 2026-09-14 02:26:26 to 2026-09-14 02:54:00

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
