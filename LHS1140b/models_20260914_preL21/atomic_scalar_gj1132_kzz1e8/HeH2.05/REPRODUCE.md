# How this model was reached

Written by `models/write_reproduce.py`, called by the runner at the end of the run it describes.

## The case

`atomic_scalar_gj1132_kzz1e8/HeH2.05`

- chemistry: H, He, He(2^3S) and electrons; no molecular network (`atomic`)
- lower boundary: prescribed scalar base, T = 226 K, R_0 = 0.157692 R_J, p = 1 microbar, metal-free (`scalar`)
- spectrum: `gj1132`
- mixing: binary H/He element diffusion with He_Kzz = 1e8 cm^2/s (`kzz1e8`)

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
| run title | `Simulation for LHS1140b_atomic_scalar_gj1132_kzz1e8_HeH2.05` |

## The seed

`models/pick_seed.py` chose, out of the states that carry this case's physics:

```
/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH2.13/output  tier2k:models/atomic_scalar_gj1132_kzz1e9  HeH=2.13  target=2.05  dlog10=0.0166  candidates=260  (a certified case of the same physics at another K_zz)
```

The fields are the state directory, the tier it was taken from (`tier1` a certified case of the same group, `tier2` a certified case of another group with the same physics, `tier3` the same at another XUV normalization, `tier4` and `tier5` the archive) with the group or the code generation it belongs to, the seed reservoir He/H, this case's, the distance |log10 He/H_seed - log10 He/H_case| the choice minimizes, and how many states of that tier carried the physics at all.

The seed was solved at He/H = 2.1299999999999999E+00 and this case is at 2.05, so `--reservoir He/H 2.05` carried its helium onto this composition: every helium column of every row was multiplied by 0.962441314554, which sets the base rows -- where the element-diffusion operator holds its Dirichlet He/H -- to 2.05 exactly and leaves the shape of the diffused He/H profile as it was. The pressure is unchanged; the mass density and the temperature were rewritten to follow the new particle count, and the `# reservoir` line of the seed states 2.05, which is what `load_IC` compares against the input.

It was interpolated onto the cell centers of the current code by

```
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH2.13/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic --reservoir He/H 2.05
```

## The commands

```bash
# the seed, interpolated onto the current cell centers
python3 /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/src/utils/map_state_to_grid.py \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH2.13/output \
    /nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/LHS1140b/models/current_grid_Hydro_ioniz.txt output --ic --reservoir He/H 2.05

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
| 1 | 0 | 7.15E-10 | 2.28E-14 | 1.63E-08 | 1.96E-05 of 1.0E-05 (elemental transport He/H partition) | 217 | 885.19 |
| 2 | 0 | 1.14E-09 | 3.88E-13 | 3.59E-08 | 6.40E-05 of 1.0E-05 (elemental transport He/H partition) | 297 | 160.49 |
| 3 | 0 | 7.94E-10 | 5.18E-14 | 1.90E-08 | 6.36E-05 of 1.0E-05 (elemental transport He/H partition) | 298 | 160.05 |
| 4 | 0 | 1.07E-09 | 4.42E-14 | 2.54E-08 | 5.91E-05 of 1.0E-05 (elemental transport He/H partition) | 298 | 161.23 |
| 5 | 0 | 9.00E-10 | 3.09E-14 | 2.42E-08 | 5.28E-05 of 1.0E-05 (elemental transport He/H partition) | 299 | 135.75 |
| 6 | 0 | 1.57E-09 | 2.16E-12 | 2.79E-08 | 4.57E-05 of 1.0E-05 (elemental transport He/H partition) | 299 | 128.65 |
| 7 | 0 | 9.55E-10 | 1.72E-12 | 2.10E-08 | 3.89E-05 of 1.0E-05 (elemental transport He/H partition) | 300 | 128.31 |
| 8 | 0 | 1.21E-09 | 1.13E-12 | 2.50E-08 | 3.25E-05 of 1.0E-05 (elemental transport He/H partition) | 300 | 127.76 |
| 9 | 0 | 1.04E-09 | 7.80E-13 | 2.35E-08 | 2.68E-05 of 1.0E-05 (elemental transport He/H partition) | 301 | 128.71 |
| 10 | 0 | 1.31E-09 | 5.08E-13 | 2.08E-08 | 2.19E-05 of 1.0E-05 (elemental transport He/H partition) | 301 | 128.72 |
| 11 | 0 | 1.15E-09 | 3.55E-13 | 2.82E-08 | 1.78E-05 of 1.0E-05 (elemental transport He/H partition) | 301 | 128.40 |
| 12 | 0 | 8.11E-10 | 2.47E-13 | 1.58E-08 | 1.43E-05 of 1.0E-05 (elemental transport He/H partition) | 302 | 126.01 |
| 13 | 0 | 6.26E-10 | 2.26E-13 | 1.82E-08 | 1.15E-05 of 1.0E-05 (elemental transport He/H partition) | 302 | 127.23 |
| 14 | 0 | 8.91E-10 | 1.25E-13 | 1.84E-08 | 9.15E-06 of 1.0E-05 (elemental transport He/H partition) | 302 | 125.71 |

- solver verdict: **info = 0** (last `||R||` = 1.843E-08), by the partitioned stationary route
- certification of the state written: **CERTIFIED**
- log10 Mdot [g/s] = **7.87**, from the post-processing pass
- He I 10830 red-pair equivalent width = **1.1748** %A over 10832.60 to 10834.20 A (air)
- red-pair depth = **4.3419** %, FWHM = 0.2543 A (three-Gaussian fit)
- wall clock **44m47s** at `OMP_NUM_THREADS=8`, 2026-09-14 03:43:58 to 2026-09-14 04:28:45

## Reproducing it

Run the commands above, in this directory, in the order they are given; each pass reads what the one before it wrote, so the order is the content and not a convenience.

The OpenMP parallelization of the cell ionization sweep is bitwise identical to the serial result (`README.md`, feature list) and the BLAS thread pool is pinned to one thread by the run itself, so `OMP_NUM_THREADS` is expected to change the wall clock and not the answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1` because that is the only setting under which it compares outputs bitwise. This run used `OMP_NUM_THREADS=8`. A different binary -- another compiler, another BLAS -- reproduces the physics and not the last digits.

The seed is another solved state, not a fresh one: a certified case of `models/` where one of the right physics exists, otherwise an archived state out of `archive_20260830/`, which is never rewritten. Either way it is interpolated onto the current grid. The solve does not depend on which seed of the right physics is used -- three seeds four decades apart in the base velocity land on the same fixed point to 1e-8 (`docs/lhs1140b_stationary_L5c_20260913.md`, T9b) -- so a rerun whose `pick_seed.py` chooses differently, because another case has been certified since, is still the same solution.
