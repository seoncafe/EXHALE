# LHS 1140 b catalog refresh, D9 step 3 (2026-09-19)

Item D9 step 3 of `docs/PLAN_20260918_rev2.md`, the catalog part. Every number is MEASURED in this item unless labeled READ.

## Verdict

- Re-evaluation of the 86 indexed cases with a `latest_certified` generation: **83 CERTIFIED, 3 REFUSED** (the three cases at 0.01 of the GJ 1132 XUV: `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13`, `.../HeH9.7`, `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7`). On the 83, the log10 Mdot the catalog reports and the He 10830 red-pair equivalent width are unchanged (largest |d EW| 3.3e-13 %A; the wind-window mass flux ratio gives d log10 Mdot = 0 at its six printed digits wherever both certificates print it).
- Re-solves: **none certified.** The three refused atomic cases reached the 30 min wall ceiling of their class (UNCERTIFIED, INCOMPLETE). The two molecular reference states could not be re-solved: the runner's molecular seed route stops in the binary (`ERROR STOP ioniz_eq: base ghost H2 closure failed`), a code defect under the frozen source, so those cases were stopped as the brief requires. `.L22/i3_alt` was published as a generation (to protect it) and not re-solved.

## The binary of record

`LHS1140b/models/EXHALE_7670f310.x`, md5 `7670f31031fb4db91d27b44cb0da6f70`, a copy (not a rebuild) of the tree's `EXHALE.x` built 2026-09-19 09:35:47 +0900; manifest `LHS1140b/models/BINARY_MANIFEST_7670f31031fb.txt` (header, then the md5 of the Makefile and of every production source, the file list of `BINARY_MANIFEST_74b96cdcf887.txt` plus `src/modules/states/boundary_state_trace.f90` and `src/utils/pin_base_grid.py`). No source file is newer than the binary (MEASURED, `find -newer`); `make -q` clean is READ from the brief. The runner and `write_reproduce.py` take the binary from `EXHALE_BIN` / `--binary` (their default is `<repo>/EXHALE.x`), as the previous catalog binary was named (READ, L34b memo section 6); `publish_state.py` finds the manifest by the md5 prefix. Boundary model of the evaluations: `characteristic_face_ps_reservoir_C_minus_contact_upwind_v2`, reservoir prescription version 1.

## How the re-evaluation was made

`run_case.sh` has no evaluate-only entry. The evaluate block of `run_case.sh` was reproduced by a scratch driver (`<scratchpad>/D9s3/evaluate_generation.sh`): the case's `latest_certified` generation is copied as the `_IC` pair into `runs/<run_id>/eval/output/`, the input is transformed exactly as the runner does (absolute paths, `Load IC? True`, `Restart intent: stationary evaluate`), the binary runs there with `OMP_NUM_THREADS=8` on `lart4`, the written state is published by `publish_state.py` as phase `evaluate` with the parent named, and the transit synthesis runs through the WINERED kernel. Two differences from the runner, both so that nothing is removed: the pass runs in `runs/<run_id>/eval/` instead of `<case>/eval/` (which the runner deletes with `rm -rf`), and the case-level products it replaces were copied to `runs/<run_id>/superseded_case_products/` before the new ones were installed. Each touched case's `REPRODUCE.md` was rewritten by `write_reproduce.py` with a note naming parent, run and generation. For every case, the case directory's `output/Hydro_ioniz.txt` was byte-identical to the parent generation before the pass, so the parent's Mdot and EW columns are those of the parent state. An evaluate pass takes 20 s to a few minutes.

## The 86 re-evaluations

"d log10 Mdot (window flux)" is log10 of the ratio of the wind-window mean of rho v r^2 printed by the new and the parent certificate (n/a where the parent certificate predates that line or the pass refused). The log10 Mdot columns are the catalog's two-decimal line.

| case | verdict | refusing entries | d log10 Mdot (window flux) | log10 Mdot parent -> new | EW [%A] parent -> new | d EW |
|---|---|---|---|---|---|---|
| `atomic_photochem_gj1132_kzzprofile/HeH10/k04` | CERTIFIED | - | 0.00e+00 | 7.9 -> 7.9 | 2.2224 -> 2.2224 | -4.2e-14 |
| `atomic_photochem_gj1132_kzzprofile/HeH12/k04` | CERTIFIED | - | 0.00e+00 | 7.9 -> 7.9 | 2.2289 -> 2.2289 | -6.6e-14 |
| `atomic_photochem_gj1132_kzzprofile/HeH2.09/k05` | CERTIFIED | - | 0.00e+00 | 7.91 -> 7.91 | 1.4765 -> 1.4765 | 1.4e-14 |
| `atomic_photochem_gj1132_kzzprofile/HeH3/k04` | CERTIFIED | - | 0.00e+00 | 7.9 -> 7.9 | 1.7640 -> 1.7640 | 1.4e-14 |
| `atomic_photochem_gj1132_kzzprofile/HeH5/k04` | CERTIFIED | - | 0.00e+00 | 7.9 -> 7.9 | 2.0746 -> 2.0746 | -2.2e-14 |
| `atomic_photochem_gj1132_kzzprofile/HeH7/k04` | CERTIFIED | - | 0.00e+00 | 7.9 -> 7.9 | 2.1768 -> 2.1768 | -1.8e-14 |
| `atomic_photochem_gj1132_kzzprofile/HeH8/k04` | CERTIFIED | - | 0.00e+00 | 7.9 -> 7.9 | 2.2004 -> 2.2004 | -2.3e-14 |
| `atomic_photochem_gj1132_kzzprofile/HeH9.7/k04` | CERTIFIED | - | 0.00e+00 | 7.9 -> 7.9 | 2.2206 -> 2.2206 | -4.3e-14 |
| `atomic_photochem_gj1132_kzzprofile/HeH9/k04` | CERTIFIED | - | 0.00e+00 | 7.9 -> 7.9 | 2.2137 -> 2.2137 | -4.4e-14 |
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | REFUSED | hydrodynamic mass row: row measure  9.910E-01 above  1.1E-08 at cell 2; hydrodynamic momentum row: row measure  1.346E-05 above  1.0E-08 at cell 1; hydrodynamic energy row: row measure  1.049E+00 above  1.0E-06 at cell 1 | n/a | 5.77 -> 5.77 | 0.0000 -> 0.0000 | 0.0e+00 |
| `atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7` | CERTIFIED | - | 0.00e+00 | 6.84 -> 6.84 | 0.2674 -> 0.2674 | -3.0e-15 |
| `atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7` | CERTIFIED | - | n/a | 7.02 -> 7.02 | 0.4822 -> 0.4822 | -2.4e-15 |
| `atomic_photochem_gj1132x0.20_kzzprofile/HeH9.7` | CERTIFIED | - | n/a | 7.16 -> 7.16 | 0.6711 -> 0.6711 | -3.1e-15 |
| `atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7` | CERTIFIED | - | n/a | 7.26 -> 7.26 | 0.8429 -> 0.8429 | -4.4e-15 |
| `atomic_photochem_gj1132x0.30_kzzprofile/HeH9.7` | CERTIFIED | - | 0.00e+00 | 7.35 -> 7.35 | 1.0097 -> 1.0097 | 2.2e-16 |
| `atomic_photochem_gj1132x0.33_kzzprofile/HeH9.7` | CERTIFIED | - | 0.00e+00 | 7.39 -> 7.39 | 1.0979 -> 1.0979 | -9.8e-15 |
| `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3721 -> 1.3721 | -6.7e-16 |
| `atomic_scalar_gj1132_kzz0/HeH0.55` | CERTIFIED | - | 0.00e+00 | 7.84 -> 7.84 | 0.0305 -> 0.0305 | 0.0e+00 |
| `atomic_scalar_gj1132_kzz0/HeH2.6` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 0.9549 -> 0.9549 | 2.2e-15 |
| `atomic_scalar_gj1132_kzz0/HeH3.0` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.0798 -> 1.0798 | 1.1e-14 |
| `atomic_scalar_gj1132_kzz0/HeH3.5` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2175 -> 1.2175 | -5.8e-15 |
| `atomic_scalar_gj1132_kzz0/HeH3.7` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2698 -> 1.2698 | 7.8e-15 |
| `atomic_scalar_gj1132_kzz0/HeH3.9` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3168 -> 1.3168 | 4.7e-15 |
| `atomic_scalar_gj1132_kzz1e10/HeH0.55` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 0.6463 -> 0.6463 | -1.1e-14 |
| `atomic_scalar_gj1132_kzz1e10/HeH1.06` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.0456 -> 1.0456 | 2.2e-16 |
| `atomic_scalar_gj1132_kzz1e10/HeH1.15` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.1032 -> 1.1032 | -1.8e-15 |
| `atomic_scalar_gj1132_kzz1e10/HeH1.29` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.1867 -> 1.1867 | -1.1e-14 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.55` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 0.8065 -> 0.8065 | 5.4e-15 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.795` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.0302 -> 1.0302 | 2.9e-15 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.865` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.0866 -> 1.0866 | -4.2e-15 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.93` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.1364 -> 1.1364 | 6.2e-15 |
| `atomic_scalar_gj1132_kzz1e5/HeH2.6` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 0.9631 -> 0.9631 | -1.7e-14 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.0` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.0880 -> 1.0880 | 1.8e-15 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.35` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1862 -> 1.1862 | -6.7e-16 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.64` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2627 -> 1.2627 | 8.9e-16 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.93` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3314 -> 1.3314 | 6.0e-15 |
| `atomic_scalar_gj1132_kzz1e6/HeH0.55` | CERTIFIED | - | 0.00e+00 | 7.84 -> 7.84 | 0.0591 -> 0.0591 | 0.0e+00 |
| `atomic_scalar_gj1132_kzz1e6/HeH2.4` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 0.9465 -> 0.9465 | -2.0e-15 |
| `atomic_scalar_gj1132_kzz1e6/HeH2.8` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.0804 -> 1.0804 | 4.4e-16 |
| `atomic_scalar_gj1132_kzz1e6/HeH3.19` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1949 -> 1.1949 | 9.5e-15 |
| `atomic_scalar_gj1132_kzz1e6/HeH3.46` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2669 -> 1.2669 | 2.3e-14 |
| `atomic_scalar_gj1132_kzz1e6/HeH3.74` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3388 -> 1.3388 | -4.4e-16 |
| `atomic_scalar_gj1132_kzz1e7/HeH0.55` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 0.1585 -> 0.1585 | 4.3e-15 |
| `atomic_scalar_gj1132_kzz1e7/HeH2.70` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1996 -> 1.1996 | -1.6e-15 |
| `atomic_scalar_gj1132_kzz1e7/HeH2.94` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2724 -> 1.2724 | 6.9e-15 |
| `atomic_scalar_gj1132_kzz1e7/HeH3.18` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3265 -> 1.3265 | -5.8e-15 |
| `atomic_scalar_gj1132_kzz1e8/HeH0.55` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 0.3089 -> 0.3089 | 8.3e-16 |
| `atomic_scalar_gj1132_kzz1e8/HeH2.05` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1513 -> 1.1513 | -2.4e-15 |
| `atomic_scalar_gj1132_kzz1e8/HeH2.23` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2206 -> 1.2206 | -1.2e-14 |
| `atomic_scalar_gj1132_kzz1e8/HeH2.41` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2826 -> 1.2826 | 8.9e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH0.55` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 0.4793 -> 0.4793 | 1.8e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH1.50` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.0996 -> 1.0996 | -8.4e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH1.60` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1511 -> 1.1511 | -5.8e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH1.70` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1970 -> 1.1970 | -8.0e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH2.13` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3749 -> 1.3749 | 0.0e+00 |
| `atomic_scalar_gj1132_kzz1e9/HeH4.0` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.8229 -> 1.8229 | 4.4e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH9.7` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 2.0592 -> 2.0592 | -9.9e-14 |
| `atomic_scalar_gj1132_wellmixed/HeH0.083` | CERTIFIED | - | 0.00e+00 | 7.84 -> 7.84 | 0.3797 -> 0.3797 | 1.7e-16 |
| `atomic_scalar_gj1132_wellmixed/HeH0.40` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.1011 -> 1.1011 | -1.0e-14 |
| `atomic_scalar_gj1132_wellmixed/HeH0.42` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.1324 -> 1.1324 | 4.2e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH0.44` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.1626 -> 1.1626 | 2.2e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH0.55` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.3118 -> 1.3118 | -9.3e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH1` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.6959 -> 1.6959 | 4.4e-16 |
| `atomic_scalar_gj1132_wellmixed/HeH10` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.9347 -> 1.9347 | -3.3e-13 |
| `atomic_scalar_gj1132_wellmixed/HeH100` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.9644 -> 1.9644 | -1.9e-13 |
| `atomic_scalar_gj1132_wellmixed/HeH1000` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.9197 -> 1.9197 | -2.9e-14 |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | REFUSED | hydrodynamic mass row: row measure  9.589E-01 above  3.7E-08 at cell 2; hydrodynamic momentum row: row measure  5.539E-06 above  1.0E-08 at cell 1; hydrodynamic energy row: row measure  1.112E+00 above  1.0E-06 at cell 1 | n/a | 5.79 -> 5.79 | 0.0000 -> 0.0000 | 0.0e+00 |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | REFUSED | hydrodynamic mass row: row measure  9.806E-01 above  1.9E-08 at cell 2; hydrodynamic momentum row: row measure  8.724E-06 above  1.0E-08 at cell 1; hydrodynamic energy row: row measure  1.071E+00 above  1.0E-06 at cell 1 | n/a | 5.78 -> 5.78 | 0.0000 -> 0.0000 | 0.0e+00 |
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13` | CERTIFIED | - | 0.00e+00 | 6.82 -> 6.82 | 0.0001 -> 0.0001 | 0.0e+00 |
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7` | CERTIFIED | - | n/a | 6.81 -> 6.81 | 0.2353 -> 0.2353 | -4.2e-15 |
| `atomic_scalar_gj1132x0.15_kzz1e9/HeH2.13` | CERTIFIED | - | 0.00e+00 | 7.0 -> 7.0 | 0.0003 -> 0.0003 | 0.0e+00 |
| `atomic_scalar_gj1132x0.15_kzz1e9/HeH9.7` | CERTIFIED | - | 0.00e+00 | 7.0 -> 7.0 | 0.4514 -> 0.4514 | 6.7e-16 |
| `atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13` | CERTIFIED | - | 0.00e+00 | 7.13 -> 7.13 | 0.0022 -> 0.0022 | 0.0e+00 |
| `atomic_scalar_gj1132x0.20_kzz1e9/HeH9.7` | CERTIFIED | - | 0.00e+00 | 7.13 -> 7.13 | 0.6294 -> 0.6294 | 0.0e+00 |
| `atomic_scalar_gj1132x0.25_kzz1e9/HeH2.13` | CERTIFIED | - | 0.00e+00 | 7.24 -> 7.24 | 0.0700 -> 0.0700 | -9.0e-16 |
| `atomic_scalar_gj1132x0.25_kzz1e9/HeH9.7` | CERTIFIED | - | 0.00e+00 | 7.23 -> 7.23 | 0.7913 -> 0.7913 | -1.1e-14 |
| `atomic_scalar_gj1132x0.30_kzz1e9/HeH2.13` | CERTIFIED | - | 0.00e+00 | 7.32 -> 7.32 | 0.1894 -> 0.1894 | -3.1e-15 |
| `atomic_scalar_gj1132x0.30_kzz1e9/HeH9.7` | CERTIFIED | - | 0.00e+00 | 7.31 -> 7.31 | 0.9395 -> 0.9395 | 2.4e-15 |
| `atomic_scalar_gj1132x0.33_kzz1e9/HeH2.13` | CERTIFIED | - | 0.00e+00 | 7.36 -> 7.36 | 0.2555 -> 0.2555 | 3.8e-15 |
| `atomic_scalar_gj1132x0.33_kzz1e9/HeH9.7` | CERTIFIED | - | 0.00e+00 | 7.36 -> 7.36 | 1.0231 -> 1.0231 | -1.1e-14 |
| `atomic_scalar_gj699_wellmixed/HeH0.042` | CERTIFIED | - | 0.00e+00 | 8.57 -> 8.57 | 0.9091 -> 0.9091 | 3.9e-15 |
| `atomic_scalar_gj699_wellmixed/HeH0.046` | CERTIFIED | - | 0.00e+00 | 8.57 -> 8.57 | 0.9857 -> 0.9857 | -1.8e-15 |
| `atomic_scalar_gj699_wellmixed/HeH0.050` | CERTIFIED | - | 0.00e+00 | 8.57 -> 8.57 | 1.0608 -> 1.0608 | 3.3e-15 |
| `atomic_scalar_gj699_wellmixed/HeH0.083` | CERTIFIED | - | 0.00e+00 | 8.57 -> 8.57 | 1.6271 -> 1.6271 | -8.2e-15 |
| `atomic_scalar_gj699_wellmixed/HeH1` | CERTIFIED | - | 0.00e+00 | 8.55 -> 8.55 | 4.5250 -> 4.5250 | 8.9e-16 |
| `atomic_scalar_gj699_wellmixed/HeH1000` | CERTIFIED | - | 0.00e+00 | 8.42 -> 8.42 | 2.8725 -> 2.8725 | -1.5e-13 |

## The three refusals

On all three, the physical columns are the parent's; the refusal is the hydrodynamic rows at the base: mass row 0.96 to 0.99 at cell 2, energy row 1.05 to 1.11 at cell 1, momentum 5.5e-6 to 1.3e-5 at cell 1. The face-flux block of `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` reads a base face mass flux of 38.4 times the wind-window mean (5.29e-9 in code units), with the direction read from the window and agreeing with the base face flux (`w_rev = 0`). The parent certificate (old binary, boundary model before D5b-2/D2b) had the mass row at 1.1e-7 of the window. These are ROE, very weak winds (log10 Mdot 5.8); the interpretation that boundary model v2 changes the base face flux of this subsonic Roe base by that amount is likely but was not isolated here.

## Re-solves

| case | seed (pick_seed) | budget | passes reached | ending | last pass (worst gated row, mass, energy) |
|---|---|---|---|---|---|
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | tier0, own generation `g0002_20260919T004834Z_d59ec9ff` | 40 p / 30 m | 26 | wall ceiling, UNCERTIFIED INCOMPLETE | 9.6e-9 of 1e-5; mass 1.0e-5; energy 6.4e-5; hydro info=2 |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | tier0, own generation `g0002_20260919T004843Z_f7485b14` | 40 p / 30 m | 11 | wall ceiling, UNCERTIFIED INCOMPLETE | 1.1e-5 of 1e-5; mass 2.2e-6; energy 2.8e-5; hydro info=2 |
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | tier0, own generation `g0002_20260919T004137Z_b2a43801` | 40 p / 30 m | 8 | wall ceiling, UNCERTIFIED INCOMPLETE | 4.5e-5 of 1e-5; mass 2.0e-5; energy 1.2e-4; hydro info=2 |
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | pick_seed: none compatible (rc 3); runner's molecular route from the certified atomic generation | 40 p / 6 h | 0 | the seed conversion stops: `base ghost H2 closure failed`, ghost cell 0, 30 passes, last move 1.19e-6 against 1e-10, residual 8.8e-7 | not started |
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | same | 40 p / 6 h | 0 | same stop | not started |
| `.L22/i3_alt` | not run | | | its state published as `g0001_20260916T195207Z_38cec4fe` (no certificate attached) | |

The atomic solves spend 60 to 220 s a pass and in every case the hydrodynamic JFNK returns info=2 from pass 2 on with the mass row rising from 1.9e-8 (pass 1) to 1e-5; the 30 min allowance of the atomic class is short of what these three need, and no pass count was reached. The wall ceiling is the campaign's own class table; no override was authorized.

The molecular comparison figures (T, v, n_H2 fraction, x(H II) before and after) were not made: no after-state exists.

## What did not certify, and the evidence

- The three x0.01 cases: refused on re-evaluation (table above), and not re-certified within the class budget. Logs: each case's `run.log`, `not_solved.md` is not written by a wall-ceiling stop; `models/campaign_status.txt`.
- The two molecular references: their certificates stay stale; the D5b-2 ghost H2 closure (`ionization_equilibrium.f90`, `base_ghost_closure_tol = 1e-10`, 30 passes) does not converge on the ghost of a state converted from an atomic seed (`<case>/seed.log`). This is a code defect under the frozen source and is reported, not fixed.

## Defects noticed in the catalog scripts (reported, not fixed)

1. `run_case.sh`'s evaluate block copies the read-only generation files (`chmod a-w` by the publisher) into `eval/output/*_IC.txt` with `cp -f`, which keeps mode 0444, and moves them into `output/`; the next run of that case with a seed then fails in `map_state_to_grid.py` with `PermissionError` on `output/Hydro_ioniz_IC.txt`. MEASURED here: the first launch of the three atomic re-solves failed that way on the files this item had installed the same way; the write bit was restored on 172 files of the 86 case `output/` directories.
2. `pick_seed.py`: `certified()` asks the index for `latest_certified`, but the seed directory comes from `state_directory()` with `latest_complete`, so where the two differ tier 0 offers the latest complete generation labeled "own most recent certified state". Here it offered the three REFUSED evaluation generations (physical cells equal to the certified parents, so the seed values are the same).
3. The index keeps `latest_certified` on the parent generations of the three refused cases although the current binary refuses them; only the molecular cases carry `stale_under`.
4. `write_reproduce.py` says "the run loaded the state already in `output/`" for an evaluate pass given no seed line; the note of each rewritten record says what was loaded.

## Records written

`LHS1140b/MODELS.md` sections 7 and 8 regenerated by `models/status.py --write`; `REPRODUCE.md` of the 86 re-evaluated cases (the rung-level `REPRODUCE.md` of the flux-closure rungs untouched); the three atomic re-solves wrote their own through the runner. Pre-run copies of the case-level files of the two molecular cases and of `.L22/i3_alt` are in `<case>/pre_D9s3/`.

## Second pass, final binary

Added 2026-09-19 (D9 step 3, the catalog part, on the final binary of the series). Every number in this section is MEASURED here unless labeled READ.

### Verdict

- **Re-evaluation: 83 of 83 CERTIFIED.** Every case whose index named a `latest_certified` generation (83 after the demotions of MODELS.md section 9.8) was handed to `models/run_case.sh --evaluate`. None refused. On all 83 the wind-window mean of rho v r^2 printed by the new certificate is character for character the parent's (d log10 Mdot = 0 at its six printed digits), the catalog's two-decimal log10 Mdot line is unchanged, and the largest |d EW| of the He 10830 red pair is 2.1e-14 %A. The evaluated parents were all the evaluate generations the first pass published on `EXHALE_7670f310.x`.
- **The two molecular reference states: seeds built, not certified.** Both seeds now build (H nuclei conserved to 4.0e-16 and 4.2e-16, equation-of-state closure 3.6e-16 and 3.9e-16). `HeH9.7` spent its 40-pass ceiling and ended `composition_refusal` (carrier balance H2 at 4.3e-4 against 1e-5, the hydrodynamic rows inside their tolerances). `HeH2.13` ended its first solve `hydrodynamic_refusal` at the pass ceiling, and the runner's dtau0 = 1e8 continuation was stopped by the 6 h wall ceiling. Neither has a certified state, so neither has an evaluate pass or a transit synthesis.
- **The three x0.01-XUV atomic cases: all stopped at the raised 3 h wall ceiling, UNCERTIFIED and INCOMPLETE.** In each, the hydrodynamic JFNK returned info = 2 at every outer pass, with a hydrodynamic row stuck at order 0.1 to 1 while the species rows fell.

### The binary of record

`LHS1140b/models/EXHALE_2c3b0acc.x`, md5 `2c3b0acc9983aec03bb4f844294fed18`. It is a copy (not a rebuild) of the tree's `EXHALE.x`, built 2026-09-19 12:02:25 +0900. The manifest is `LHS1140b/models/BINARY_MANIFEST_2c3b0acc9983.txt`, made the same way as `BINARY_MANIFEST_7670f31031fb.txt` (header, then the md5 of the Makefile and of every production source, the same 169-file list). No file of `src/` or the Makefile is newer than the binary (MEASURED, `find -newer`). `make -q` clean is READ from the brief. Against the 7670f310 manifest, two production sources differ (MEASURED): `src/modules/radiation/ionization_equilibrium.f90` (the D9fix secant ghost closure and the `bg_ready` restoration) and `src/modules/init/molecular_seed_from_atomic_state.f90` (`transfer_h2`). Every run of this pass named the binary through `EXHALE_BIN`, and every generation it published carries `binary_manifest: BINARY_MANIFEST_2c3b0acc9983.txt`.

### How the pass was run

All on `lart4`, `OMP_NUM_THREADS=8`, at most three runs at a time, from one queue (`<scratchpad>/D9s3b/jobs.txt`, driver `job.sh`). The evaluations went through `run_case.sh --evaluate <case>`. The five re-solves went through `run_campaign.sh 1 8` with `CAMPAIGN_LIST` naming one case and `FORCE=1`, so each carries its budget, its ceiling and its status line in `models/campaign_status*.txt`. The two molecular cases ran on the class table (`molecular_alternation`, 40 passes, 6 h). The three x0.01 cases ran on `EXHALE_OUTER_PASSES=90 CAMPAIGN_WALL=3h` (the status lines read `source=environment`). The override was authorized by the brief of this pass, for this reason: under the v2 contact law the base is a different stationary problem from the one their seeds solved (the D9fix memo, section 6). No `models/budget_overrides.txt` was written, so the override does not outlive this pass. Each evaluation took 20 to 36 s of wall clock.

### Step 1: the 83 re-evaluations

"d log10 Mdot (window flux)" is log10 of the ratio of the wind-window mean of rho v r^2 in the new certificate (`runs/<run>/pp.log`) to the one in the parent's `certification.txt`. The log10 Mdot columns are the catalog's two-decimal line from `pp.log`. The EW is the He 10830 red-pair equivalent width from `tpm_He10830.txt`, as `status.py` reads it. The parent columns were read before the pass, and in every case the case directory's `output/Hydro_ioniz.txt` held the parent's data rows at that point. After the pass, each case's `latest_certified` is its new evaluate child (83 of 83).

| case | evaluated parent | evaluate child | verdict | refusing entries | d log10 Mdot (window flux) | log10 Mdot parent -> new | EW [%A] parent -> new | d EW [%A] |
|---|---|---|---|---|---|---|---|---|
| `atomic_photochem_gj1132_kzzprofile/HeH10/k04` | `g0003_20260919T003952Z_ede9d72d` | `g0004_20260919T032756Z_093629d3` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 2.2224 -> 2.2224 | -2.1e-14 |
| `atomic_photochem_gj1132_kzzprofile/HeH12/k04` | `g0003_20260919T003952Z_84579f58` | `g0004_20260919T032832Z_be7a184f` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 2.2289 -> 2.2289 | -4.9e-15 |
| `atomic_photochem_gj1132_kzzprofile/HeH2.09/k05` | `g0003_20260919T003952Z_20d38273` | `g0004_20260919T032907Z_a39f11a3` | CERTIFIED | - | 0.00e+00 | 7.91 -> 7.91 | 1.4765 -> 1.4765 | -7.8e-15 |
| `atomic_photochem_gj1132_kzzprofile/HeH3/k04` | `g0003_20260919T004026Z_f0c4b5c0` | `g0004_20260919T032941Z_5b58b90c` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 1.7640 -> 1.7640 | 2.2e-15 |
| `atomic_photochem_gj1132_kzzprofile/HeH5/k04` | `g0003_20260919T004028Z_7d1d1800` | `g0004_20260919T033015Z_0f70dd16` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 2.0746 -> 2.0746 | 1.5e-14 |
| `atomic_photochem_gj1132_kzzprofile/HeH7/k04` | `g0003_20260919T004029Z_e86e15cc` | `g0004_20260919T033050Z_9e2b0d4a` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 2.1768 -> 2.1768 | 8.0e-15 |
| `atomic_photochem_gj1132_kzzprofile/HeH8/k04` | `g0003_20260919T004102Z_4ab2c2ef` | `g0004_20260919T033125Z_c0d39208` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 2.2004 -> 2.2004 | 3.1e-15 |
| `atomic_photochem_gj1132_kzzprofile/HeH9.7/k04` | `g0003_20260919T004104Z_07a6c8b9` | `g0004_20260919T033200Z_1878856e` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 2.2206 -> 2.2206 | -7.1e-15 |
| `atomic_photochem_gj1132_kzzprofile/HeH9/k04` | `g0003_20260919T004105Z_cd0ba4b9` | `g0004_20260919T033236Z_fa443a57` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 2.2137 -> 2.2137 | -8.4e-15 |
| `atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7` | `g0003_20260919T004140Z_f44690af` | `g0004_20260919T033311Z_52be5f9b` | CERTIFIED | - | 0.00e+00 | 6.84 -> 6.84 | 0.2674 -> 0.2674 | 8.9e-16 |
| `atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7` | `g0002_20260919T004141Z_6128e1b8` | `g0003_20260919T033344Z_a49f9ac8` | CERTIFIED | - | 0.00e+00 | 7.02 -> 7.02 | 0.4822 -> 0.4822 | -2.2e-15 |
| `atomic_photochem_gj1132x0.20_kzzprofile/HeH9.7` | `g0002_20260919T004211Z_48a64fc4` | `g0003_20260919T033418Z_c3c26194` | CERTIFIED | - | 0.00e+00 | 7.16 -> 7.16 | 0.6711 -> 0.6711 | -6.4e-15 |
| `atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7` | `g0002_20260919T004215Z_6c42f1af` | `g0003_20260919T033451Z_e0e2e60b` | CERTIFIED | - | 0.00e+00 | 7.26 -> 7.26 | 0.8429 -> 0.8429 | 0.0e+00 |
| `atomic_photochem_gj1132x0.30_kzzprofile/HeH9.7` | `g0003_20260919T004216Z_6a7a7eef` | `g0004_20260919T033526Z_cb15d8ca` | CERTIFIED | - | 0.00e+00 | 7.35 -> 7.35 | 1.0097 -> 1.0097 | 8.9e-16 |
| `atomic_photochem_gj1132x0.33_kzzprofile/HeH9.7` | `g0003_20260919T004246Z_7036fabc` | `g0004_20260919T033600Z_261c5da7` | CERTIFIED | - | 0.00e+00 | 7.39 -> 7.39 | 1.0979 -> 1.0979 | 7.5e-15 |
| `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13` | `g0003_20260919T004250Z_fa16692f` | `g0004_20260919T033634Z_85944973` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3721 -> 1.3721 | 2.7e-15 |
| `atomic_scalar_gj1132_kzz0/HeH0.55` | `g0003_20260919T004250Z_fa837b6a` | `g0004_20260919T033707Z_88eb0d69` | CERTIFIED | - | 0.00e+00 | 7.84 -> 7.84 | 0.0305 -> 0.0305 | -8.4e-16 |
| `atomic_scalar_gj1132_kzz0/HeH2.6` | `g0003_20260919T004311Z_622b5f8b` | `g0004_20260919T033731Z_df801154` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 0.9549 -> 0.9549 | 2.9e-15 |
| `atomic_scalar_gj1132_kzz0/HeH3.0` | `g0003_20260919T004319Z_f4d15da3` | `g0004_20260919T033753Z_5507044b` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.0798 -> 1.0798 | -9.8e-15 |
| `atomic_scalar_gj1132_kzz0/HeH3.5` | `g0003_20260919T004322Z_1d768ce0` | `g0004_20260919T033814Z_51875900` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2175 -> 1.2175 | 4.0e-15 |
| `atomic_scalar_gj1132_kzz0/HeH3.7` | `g0003_20260919T004332Z_33447ece` | `g0004_20260919T033836Z_c901733e` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2698 -> 1.2698 | -5.3e-15 |
| `atomic_scalar_gj1132_kzz0/HeH3.9` | `g0003_20260919T004339Z_bc92fd26` | `g0004_20260919T033859Z_79cfb0fc` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3168 -> 1.3168 | 1.3e-14 |
| `atomic_scalar_gj1132_kzz1e10/HeH0.55` | `g0003_20260919T004343Z_d4e5438a` | `g0004_20260919T033921Z_4ef7ddb7` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 0.6463 -> 0.6463 | 8.7e-15 |
| `atomic_scalar_gj1132_kzz1e10/HeH1.06` | `g0003_20260919T004353Z_676a874b` | `g0004_20260919T033943Z_35b873ad` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.0456 -> 1.0456 | -7.1e-15 |
| `atomic_scalar_gj1132_kzz1e10/HeH1.15` | `g0003_20260919T004400Z_49e99e6e` | `g0004_20260919T034004Z_74ef5931` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.1032 -> 1.1032 | 1.8e-15 |
| `atomic_scalar_gj1132_kzz1e10/HeH1.29` | `g0003_20260919T004404Z_790dc4af` | `g0004_20260919T034026Z_58ad4651` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.1867 -> 1.1867 | 2.2e-15 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.55` | `g0003_20260919T004413Z_c9f15592` | `g0004_20260919T034048Z_dbc1f4a8` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 0.8065 -> 0.8065 | -1.8e-15 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.795` | `g0003_20260919T004421Z_a795c719` | `g0004_20260919T034110Z_926dbb2d` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.0302 -> 1.0302 | -1.1e-14 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.865` | `g0003_20260919T004424Z_4a56fcd8` | `g0004_20260919T034132Z_9af4eb31` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.0866 -> 1.0866 | -6.7e-16 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.93` | `g0003_20260919T004434Z_6b4b84ce` | `g0004_20260919T034155Z_8abe4779` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.1364 -> 1.1364 | 2.4e-15 |
| `atomic_scalar_gj1132_kzz1e5/HeH2.6` | `g0003_20260919T004441Z_50269e8d` | `g0004_20260919T034217Z_39ae932b` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 0.9631 -> 0.9631 | 7.9e-15 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.0` | `g0003_20260919T004445Z_6dc72714` | `g0004_20260919T034239Z_fb323ecf` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.0880 -> 1.0880 | -8.9e-15 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.35` | `g0003_20260919T004454Z_e67035c0` | `g0004_20260919T034301Z_99915b3f` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1862 -> 1.1862 | 7.1e-15 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.64` | `g0003_20260919T004501Z_7194de59` | `g0004_20260919T034323Z_3b108aa2` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2627 -> 1.2627 | -1.2e-14 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.93` | `g0003_20260919T004505Z_fc653ea7` | `g0004_20260919T034345Z_7b4d3899` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3314 -> 1.3314 | -8.9e-16 |
| `atomic_scalar_gj1132_kzz1e6/HeH0.55` | `g0003_20260919T004515Z_2e1ab4d2` | `g0004_20260919T034407Z_30c61e49` | CERTIFIED | - | 0.00e+00 | 7.84 -> 7.84 | 0.0591 -> 0.0591 | 0.0e+00 |
| `atomic_scalar_gj1132_kzz1e6/HeH2.4` | `g0003_20260919T004522Z_2137a84b` | `g0004_20260919T034429Z_e43909f6` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 0.9465 -> 0.9465 | 1.3e-15 |
| `atomic_scalar_gj1132_kzz1e6/HeH2.8` | `g0003_20260919T004526Z_d4e9c90f` | `g0004_20260919T034451Z_5798e5d8` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.0804 -> 1.0804 | -7.5e-15 |
| `atomic_scalar_gj1132_kzz1e6/HeH3.19` | `g0003_20260919T004535Z_33de63df` | `g0004_20260919T034513Z_31fca486` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1949 -> 1.1949 | -1.3e-14 |
| `atomic_scalar_gj1132_kzz1e6/HeH3.46` | `g0003_20260919T004543Z_e6325938` | `g0004_20260919T034535Z_f780d699` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2669 -> 1.2669 | -1.7e-14 |
| `atomic_scalar_gj1132_kzz1e6/HeH3.74` | `g0003_20260919T004546Z_455211ef` | `g0004_20260919T034558Z_6c7bd82d` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3388 -> 1.3388 | 1.1e-15 |
| `atomic_scalar_gj1132_kzz1e7/HeH0.55` | `g0003_20260919T004556Z_9af58e20` | `g0004_20260919T034620Z_553ec42c` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 0.1585 -> 0.1585 | -3.9e-15 |
| `atomic_scalar_gj1132_kzz1e7/HeH2.70` | `g0003_20260919T004604Z_d8d9a9cb` | `g0004_20260919T034642Z_e5b7019c` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1996 -> 1.1996 | 8.2e-15 |
| `atomic_scalar_gj1132_kzz1e7/HeH2.94` | `g0003_20260919T004607Z_e662eba8` | `g0004_20260919T034704Z_11578937` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2724 -> 1.2724 | -1.4e-14 |
| `atomic_scalar_gj1132_kzz1e7/HeH3.18` | `g0003_20260919T004616Z_ca2a9830` | `g0004_20260919T034726Z_9028a211` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3265 -> 1.3265 | 8.9e-16 |
| `atomic_scalar_gj1132_kzz1e8/HeH0.55` | `g0003_20260919T004624Z_4318a66e` | `g0004_20260919T034749Z_fd533291` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 0.3089 -> 0.3089 | -2.1e-15 |
| `atomic_scalar_gj1132_kzz1e8/HeH2.05` | `g0003_20260919T004628Z_1a20f250` | `g0004_20260919T034810Z_36bd393a` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1513 -> 1.1513 | 6.4e-15 |
| `atomic_scalar_gj1132_kzz1e8/HeH2.23` | `g0003_20260919T004636Z_0fe62914` | `g0004_20260919T034832Z_f2bb678b` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2206 -> 1.2206 | -6.4e-15 |
| `atomic_scalar_gj1132_kzz1e8/HeH2.41` | `g0003_20260919T004644Z_e208d14e` | `g0004_20260919T034854Z_9a5c1f0c` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2826 -> 1.2826 | -9.8e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH0.55` | `g0003_20260919T004648Z_8cc392ac` | `g0004_20260919T034917Z_098da21a` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 0.4793 -> 0.4793 | -3.4e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH1.50` | `g0003_20260919T004657Z_b6ecf7a4` | `g0004_20260919T034938Z_4274fa5a` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.0996 -> 1.0996 | 1.8e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH1.60` | `g0003_20260919T004705Z_a1e1f0a0` | `g0004_20260919T035000Z_a0851042` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1511 -> 1.1511 | 2.2e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH1.70` | `g0003_20260919T004708Z_5236c813` | `g0004_20260919T035022Z_fc02686a` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1970 -> 1.1970 | 8.9e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH2.13` | `g0003_20260919T003908Z_c1794389` | `g0004_20260919T032647Z_aece8741` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3749 -> 1.3749 | -1.3e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH4.0` | `g0003_20260919T004717Z_8a8c7580` | `g0004_20260919T035044Z_acba83cc` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.8229 -> 1.8229 | -1.2e-14 |
| `atomic_scalar_gj1132_kzz1e9/HeH9.7` | `g0003_20260919T004726Z_cc4f112e` | `g0004_20260919T032647Z_9ca293a8` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 2.0592 -> 2.0592 | 9.3e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH0.083` | `g0003_20260919T004729Z_338777be` | `g0004_20260919T035107Z_5338d0c2` | CERTIFIED | - | 0.00e+00 | 7.84 -> 7.84 | 0.3797 -> 0.3797 | -1.7e-16 |
| `atomic_scalar_gj1132_wellmixed/HeH0.40` | `g0003_20260919T004739Z_ef5bc5d9` | `g0004_20260919T035130Z_bec4f54e` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.1011 -> 1.1011 | 5.3e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH0.42` | `g0003_20260919T004748Z_c0f8bfcc` | `g0004_20260919T035152Z_ab79b2da` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.1324 -> 1.1324 | -1.8e-14 |
| `atomic_scalar_gj1132_wellmixed/HeH0.44` | `g0003_20260919T004750Z_327a1950` | `g0004_20260919T035216Z_ee409f97` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.1626 -> 1.1626 | 8.9e-16 |
| `atomic_scalar_gj1132_wellmixed/HeH0.55` | `g0003_20260919T004800Z_260d42ce` | `g0004_20260919T035238Z_43fb8ba9` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.3118 -> 1.3118 | 4.9e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH1` | `g0003_20260919T004832Z_cdfcb797` | `g0004_20260919T035301Z_ecc54ed6` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.6959 -> 1.6959 | -1.5e-14 |
| `atomic_scalar_gj1132_wellmixed/HeH10` | `g0003_20260919T004821Z_f43e2847` | `g0004_20260919T035323Z_f3e7ada1` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.9347 -> 1.9347 | 2.2e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH100` | `g0003_20260919T004812Z_add0f75c` | `g0004_20260919T035347Z_090396b9` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.9644 -> 1.9644 | 8.7e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH1000` | `g0003_20260919T004809Z_df8c2885` | `g0004_20260919T035411Z_6b8ed8aa` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.9197 -> 1.9197 | -5.3e-15 |
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13` | `g0003_20260919T004853Z_db160327` | `g0004_20260919T035436Z_89f3eb8e` | CERTIFIED | - | 0.00e+00 | 6.82 -> 6.82 | 0.0001 -> 0.0001 | 0.0e+00 |
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7` | `g0002_20260919T004854Z_6783c824` | `g0003_20260919T035456Z_661ac065` | CERTIFIED | - | 0.00e+00 | 6.81 -> 6.81 | 0.2353 -> 0.2353 | 2.0e-15 |
| `atomic_scalar_gj1132x0.15_kzz1e9/HeH2.13` | `g0003_20260919T004903Z_78ae8980` | `g0004_20260919T035518Z_35541992` | CERTIFIED | - | 0.00e+00 | 7.00 -> 7.00 | 0.0003 -> 0.0003 | 0.0e+00 |
| `atomic_scalar_gj1132x0.15_kzz1e9/HeH9.7` | `g0003_20260919T004912Z_d471c75f` | `g0004_20260919T035539Z_9860d32b` | CERTIFIED | - | 0.00e+00 | 7.00 -> 7.00 | 0.4514 -> 0.4514 | 4.5e-15 |
| `atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13` | `g0003_20260919T004915Z_90cacf31` | `g0004_20260919T035601Z_c2dce3ee` | CERTIFIED | - | 0.00e+00 | 7.13 -> 7.13 | 0.0022 -> 0.0022 | 0.0e+00 |
| `atomic_scalar_gj1132x0.20_kzz1e9/HeH9.7` | `g0003_20260919T004923Z_1060ca96` | `g0004_20260919T035622Z_64f988b8` | CERTIFIED | - | 0.00e+00 | 7.13 -> 7.13 | 0.6294 -> 0.6294 | -8.9e-16 |
| `atomic_scalar_gj1132x0.25_kzz1e9/HeH2.13` | `g0003_20260919T004933Z_4b2732b3` | `g0004_20260919T035644Z_508d351b` | CERTIFIED | - | 0.00e+00 | 7.24 -> 7.24 | 0.0700 -> 0.0700 | 0.0e+00 |
| `atomic_scalar_gj1132x0.25_kzz1e9/HeH9.7` | `g0003_20260919T004934Z_4ea21940` | `g0004_20260919T035705Z_c53efa85` | CERTIFIED | - | 0.00e+00 | 7.23 -> 7.23 | 0.7913 -> 0.7913 | 6.7e-15 |
| `atomic_scalar_gj1132x0.30_kzz1e9/HeH2.13` | `g0003_20260919T004944Z_7b58db80` | `g0004_20260919T035728Z_007ecf98` | CERTIFIED | - | 0.00e+00 | 7.32 -> 7.32 | 0.1894 -> 0.1894 | 4.0e-15 |
| `atomic_scalar_gj1132x0.30_kzz1e9/HeH9.7` | `g0003_20260919T004953Z_885760aa` | `g0004_20260919T035749Z_e6f327b3` | CERTIFIED | - | 0.00e+00 | 7.31 -> 7.31 | 0.9395 -> 0.9395 | -7.8e-15 |
| `atomic_scalar_gj1132x0.33_kzz1e9/HeH2.13` | `g0003_20260919T004956Z_85aad69f` | `g0004_20260919T035812Z_76873a9d` | CERTIFIED | - | 0.00e+00 | 7.36 -> 7.36 | 0.2555 -> 0.2555 | 0.0e+00 |
| `atomic_scalar_gj1132x0.33_kzz1e9/HeH9.7` | `g0003_20260919T005006Z_9e11fea8` | `g0004_20260919T035834Z_14c7d6b1` | CERTIFIED | - | 0.00e+00 | 7.36 -> 7.36 | 1.0231 -> 1.0231 | -1.4e-14 |
| `atomic_scalar_gj699_wellmixed/HeH0.042` | `g0003_20260919T005014Z_0d25e57c` | `g0004_20260919T035856Z_f07b2ff2` | CERTIFIED | - | 0.00e+00 | 8.57 -> 8.57 | 0.9091 -> 0.9091 | -2.1e-15 |
| `atomic_scalar_gj699_wellmixed/HeH0.046` | `g0003_20260919T005016Z_1eb348bd` | `g0004_20260919T035920Z_189e1b29` | CERTIFIED | - | 0.00e+00 | 8.57 -> 8.57 | 0.9857 -> 0.9857 | 0.0e+00 |
| `atomic_scalar_gj699_wellmixed/HeH0.050` | `g0003_20260919T005027Z_d476bbee` | `g0004_20260919T035944Z_5b6a0b9f` | CERTIFIED | - | 0.00e+00 | 8.57 -> 8.57 | 1.0608 -> 1.0608 | 1.8e-15 |
| `atomic_scalar_gj699_wellmixed/HeH0.083` | `g0003_20260919T005037Z_302221be` | `g0004_20260919T040008Z_ae236f2e` | CERTIFIED | - | 0.00e+00 | 8.57 -> 8.57 | 1.6271 -> 1.6271 | 1.3e-15 |
| `atomic_scalar_gj699_wellmixed/HeH1` | `g0003_20260919T005050Z_a3f71e49` | `g0004_20260919T040031Z_2091c6ca` | CERTIFIED | - | 0.00e+00 | 8.55 -> 8.55 | 4.5250 -> 4.5250 | -1.3e-14 |
| `atomic_scalar_gj699_wellmixed/HeH1000` | `g0003_20260919T005039Z_d1f5a0f0` | `g0004_20260919T040054Z_b13d83c8` | CERTIFIED | - | 0.00e+00 | 8.42 -> 8.42 | 2.8725 -> 2.8725 | -2.7e-15 |

**The records.** The runner's record step rewrote `REPRODUCE.md` in all 83 case directories. Each now says, under "The seed", that this was an evaluate-only pass which loaded `states/<latest_certified>/` in `runs/<run>/eval/` and published the named `evaluate` child. Before the pass, 89 `REPRODUCE.md` files of `models/` (dot-directories included) carried the old sentence "the run loaded the state already in `output/`" (MEASURED, `grep -rl`; the brief's 74 counted a narrower set). After the pass, **6** carry it, all outside this set:

- `.L35/g1000_hllc/`, `.L35/g1000_roe/`, `.L35/g2000_hllc/`, `.L35/g2000_roe/`: the L35 grid and flux study directories. They are not catalog cases, have no `latest_certified`, and were not evaluated.
- `molecular_scalar_gj1132_kzz1e9/HeH2.13/pre_D9s3/` and `.../HeH9.7/pre_D9s3/`: the pre-run copies the first pass kept. They are archival and are not rewritten.

### Step 2: the two molecular reference states

Route: `run_case.sh`, molecular branch. `pick_seed.py` offers no compatible state for either case (rc 0, "start cold": 3030 candidates refused), so the runner builds the seed from the atomic case of the same name and He/H. It uses that case's `latest_certified`, which is the step-1 evaluate child certified by this binary: `atomic_scalar_gj1132_kzz1e9/HeH2.13/states/g0004_20260919T032647Z_aece8741` and `.../HeH9.7/states/g0004_20260919T032647Z_9ca293a8`, with x2 `local`. Budget: the `molecular_alternation` class, 40 passes and 6 h.

| case | seed conversion | solve 1 (dtau0 = 1) | continuation (dtau0 = 1e8) | ending, ceiling | published | wall |
|---|---|---|---|---|---|---|
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | built: H nuclei 4.0e-16, mass 3.7e-16, EOS closure 3.6e-16 | 40 passes, `hydrodynamic_refusal`: mass row 3.6e-5 above 3.0e-12 at cell 287; gated carrier balance H2 6.8e-2 at cell 227 at pass 40; hydro info = 2 at every pass except 4 and 5 | pass 1 info = 2, pass 2 info = 0 with the carrier row 1.9e-1 at cell 227; stopped by the wall ceiling during pass 3 | wall ceiling, UNCERTIFIED INCOMPLETE | `g0003_20260919T081651Z_3b86c9ed` (solve 1, NOT CERTIFIED; kept whole in `solve_1/`) | 6 h 00 m |
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | built: H nuclei 4.2e-16, mass 4.3e-16, EOS closure 3.9e-16 | 40 passes, hydro info = 0 at every pass from 8 on (2 at passes 2 to 4, 6 and 7); gated carrier balance H2 1.0 (pass 1), 2.0e-2 (13), 2.0e-3 (29), 4.3e-4 (40) at cell 263 | not taken (the continuation addresses `hydrodynamic_refusal` only) | `composition_refusal`, pass ceiling, UNCERTIFIED INCOMPLETE | `g0003_20260919T052326Z_ffa78b36` (NOT CERTIFIED) | 1 h 55 m |

`HeH9.7` is converging, just not within its budget. From pass 20 on, the gated row falls by a factor 0.87 a pass at omega = 0.125 and trust 2.5e-3 (pass 39 to 40: 4.96e-4 to 4.31e-4). At that rate, 1e-5 would take about 27 more passes. That is an extrapolation, not a measurement. `HeH2.13` is not converging: its hydrodynamic rows sit at 1e-2 to 1 for most of solve 1, and its carrier row stays between 7e-2 and 2e-1.

"Then evaluate the result" could not be done. `run_case.sh` runs the evaluate pass only after a solve it accepts, and `--evaluate` needs a `latest_certified`, which neither case has (both indexes: `latest_certified` null, `latest_complete` the new solve generation). The comparison below is therefore of the published solve states, and neither of them is certified.

Figures: `docs/figures/lhs1140b_D9s3b_molecular_HeH2.13_overlay.pdf` and `docs/figures/lhs1140b_D9s3b_molecular_HeH9.7_overlay.pdf`. Each overlays T, v, the fraction of H nuclei in H2 (2 n_H2 / n_H, with n_H = n_HI + n_HII + 2 n_H2 + 2 n_H2+ + 3 n_H3+ + n_HeH+) and x(H II) = n_HII / n_H against r/R_p. "Before" is the generation the catalog published before this pass (`latest_complete` g0002, the stale L34b-era evaluation); "after" is the new solve generation g0003. Read from the two states:

| case | max abs(dT/T) | 2 n_H2/n_H at the first cell | 2 n_H2/n_H at 20 R_p | x(H II) at 20 R_p | T max [K] | v at r_max [cm/s] |
|---|---|---|---|---|---|---|
| HeH2.13, before -> after | 0.64 (at the first cell) | 0.335 -> 0.802 | 3.1e-5 -> 2.8e-3 | 0.529 -> 0.559 | 5841 -> 6045 | 1.20e5 -> 1.12e5 |
| HeH9.7, before -> after | 2.6e-4 | 0.0676 -> 0.0676 | 1.83e-6 -> 1.82e-6 | 0.578 -> 0.578 | 7538 -> 7539 | 1.064e5 -> 1.064e5 |

The HeH9.7 state is the old one to 3e-4 everywhere. The HeH2.13 state differs at the base (more H2, a colder first cell) and carries 90 times more H2 in the outer wind, from 3 R_p outward. That outer H2 lies in the region where the gated carrier row still refuses (6.8e-2 at cell 227), so it is a property of an unconverged state and not a result.

### Step 3: the three x0.01-XUV atomic cases

One re-solve each through `run_campaign.sh` at 90 passes and 3 h (override above), seeded by `pick_seed.py`. Now that tier 0 reads `latest_certified`, which these three no longer have, tier 3 gives a certified case at 0.10 of the XUV. These are the step-1 evaluate children, certified by this binary.

| case | seed (pick_seed) | what ran | where it stopped | the rows at the last pass |
|---|---|---|---|---|
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | tier3, `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13/states/g0004_20260919T035436Z_89f3eb8e` | one solve, 32 outer passes complete | wall ceiling (3 h) in pass 33, UNCERTIFIED INCOMPLETE; nothing new published (`latest_complete` still `g0002_20260919T004834Z_d59ec9ff`) | hydro info = 2 at every pass; mass row 2.1e-1 (pass 1), then 1.52e-1 unchanged from pass 12 on, with the JFNK worst row the mass of cells 231 to 232; energy 1.37e-1; species row (He/H partition, cell 500) 5.1e-2 at pass 1, 3.9e-4 at pass 32, halving a pass at the end |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | tier3, `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7/states/g0003_20260919T035456Z_661ac065` | one solve, 34 outer passes complete | wall ceiling in pass 35, UNCERTIFIED INCOMPLETE; nothing new published (`latest_complete` still `g0002_20260919T004843Z_f7485b14`) | hydro info = 2 at every pass; energy row 0.98 to 0.998 from pass 1, JFNK worst row at cell 1 (energy, mass and momentum); mass 7.3e-8, momentum 8.3e-14; species row 1.3e-2 (pass 1) to 4.6e-12 (pass 34) |
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | attempt 1: tier3, `atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7/states/g0004_20260919T033311Z_52be5f9b`; attempt 2 (the runner's next candidate, `SEED_ATTEMPTS` default 3): `atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7/states/g0003_20260919T033344Z_a49f9ac8` | attempt 1: solve 1, 36 passes, then the dtau0 = 1e8 continuation, 14 passes; attempt 2: 1 pass | wall ceiling in attempt 2, pass 2, UNCERTIFIED INCOMPLETE; published `g0003_20260919T091147Z_279c76bd` (attempt 1, solve 1) and `g0004_20260919T094145Z_259fe9c3` (attempt 1, continuation), both `hydrodynamic_refusal`, NOT CERTIFIED; attempt 1 kept in `attempt_1/` | attempt 1 solve 1: mass row 0.578 at cell 1 from pass 1 to pass 36, species row to 3.8e-12, stopped by the stall rule; continuation: mass 0.692, energy 1.95, species 3.5e-2, stopped by the stall rule at pass 14; attempt 2 pass 1: mass 2.6e-8, energy 1.01, JFNK step length down to 9.5e-7 at dtau 7.4e-5 |

In all three the composition converges and the hydrodynamic rows do not move, at the base for `HeH9.7` and the photochem case and in the interior mass row for `HeH2.13`. This agrees with the D9fix finding that under the v2 contact law the base these weak winds (log10 Mdot 5.8) must reach is not the one their seeds carry. Why the hydrodynamic JFNK makes no progress toward it from either seed was not isolated here. The "more time" hypothesis is refuted for all three: with 6 times the wall clock and 2.25 times the passes of the first re-solves, the stuck rows kept their values.

### Records written

- `LHS1140b/MODELS.md` sections 7 and 8, by `models/status.py --write` (rc 0). Section 7 moves three rows: `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` now reads `info=1 uncertified hydrodynamic_refusal` with its residual 1.95; the two molecular rows read `info=1 uncertified hydrodynamic_refusal` (HeH2.13) and `info=1 uncertified composition_refusal` (HeH9.7). Section 8's wall column for the 83 evaluated cases now reads the evaluation's wall clock (see the note below).
- `REPRODUCE.md` of the 83 evaluated cases and of `molecular_scalar_gj1132_kzz1e9/HeH9.7`, written by the runner.
- New generations: 83 `evaluate` children, the two molecular solve generations, and the two photochem x0.01 solve generations. Nothing was deleted. Every file an evaluation replaced in a case directory is in `runs/<run>/superseded_case_products/`.
- `models/campaign_status.txt` and its dated predecessors `campaign_status_20260919{103836,142326,160120,172327,182756}.txt`: one per campaign of this pass.

### Noticed, not changed

1. A run stopped by the wall ceiling writes no record. `timeout` sends TERM to the runner's process group, and `run_case.sh` has no trap that runs `record()`. So `molecular_scalar_gj1132_kzz1e9/HeH2.13/REPRODUCE.md` and `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7/REPRODUCE.md` still describe the first pass (written 09:54), although both directories now hold newer published generations, a `solve_1/` and an `attempt_1/`. The two x0.01 scalar cases published nothing and keep their earlier records. The campaign status lines and the per-attempt `ENDING` files are the only records of these stops.
2. `status.py`'s section 8 takes the wall clock from `REPRODUCE.md`. After an evaluate-only record, that is the evaluation's 20 to 36 s and not the solve's, so the column no longer says how long the state took to solve. The seed and pass columns of those rows still read "cold start" and the solve's pass count, as before this pass.
3. Within one run, `run_case.sh` tried a second seed (`SEED_ATTEMPTS`, default 3) on the photochem x0.01 case after attempt 1 ended `hydrodynamic_refusal`. That is the runner's policy, but it means the "one re-solve" of the brief spent 2 h 40 m of its 3 h on attempt 1 (16:01 to 18:41) and the last 20 minutes on another seed.

## Runner records and the HeH9.7 continuation (2026-09-19, D9 step 3c)

Every number in this section is MEASURED unless labeled READ. Binary `LHS1140b/models/EXHALE_2c3b0acc.x` through `EXHALE_BIN`, `lart4`, `OMP_NUM_THREADS=8`.

### Verdict

- The three runner defects of "Noticed, not changed" above are fixed and tested (`models/tests/runner_records.sh`, 18 checks: 17 fail on the earlier tools, 18 pass now). The two records those defects left stale are rewritten.
- `molecular_scalar_gj1132_kzz1e9/HeH9.7` SOLVED and CERTIFIED at outer pass 9 of the continuation (generation `g0004_20260919T103629Z_9b8394a3`, info = 0, carrier balance H2 row 8.82e-6 against 1e-5 at cell 262). The runner's own evaluate pass then REFUSED it on the hydrodynamic mass row at cell 1 (8.335e-8 against the cell's rounding-anchored tolerance 7.5e-9; in the solve's own final block the row was 3.873e-9 there). By the refused re-evaluation rule the solve generation is `stale_under` the binary and the case has no `latest_certified`, so `run_case.sh --evaluate` has nothing to evaluate and was not run.

### The runner fixes

1. A stop from outside is recorded. `run_case.sh` runs each binary in the background, waits for it, and traps TERM, INT and HUP. The trap stops the binary by the PID the runner holds, publishes a complete state the interrupted pass wrote itself, and writes `ENDING`, `not_solved.md`, `run.json` and `REPRODUCE.md` with the class `stopped_by_wall_ceiling` (when `run_campaign.sh` told it the ceiling and the ceiling is reached) or `stopped_by_signal`. The campaign's status line uses the same class name (it read `wall_ceiling_reached` before).
2. Section 8 of `MODELS.md` states the solve. The `record` column says whether the published generation is a solve, an evaluation made in the same run, or an evaluate-only record (with that evaluation's wall clock). Seed, outer passes, verdict and wall clock are those of the solve generation the parent chain reaches (`publish_state.solve_ancestor`), read from its `run.json`, from the log that carries its certification block verbatim, or from a solve `REPRODUCE.md` whose run window holds the state's write time. 79 of the evaluated atomic rows end at an imported evaluate generation with no parent and read `not recorded`, with the reason. The solve records of those states were overwritten on 2026-09-18, before the index existed.
3. Seed attempts: one budgeted re-solve is one solve. An attempt after the first is taken only when the one before it ended before its first outer pass (`another_seed_allowed`). The campaign status header states the rule and each line states `seeds=<used>/<allowed>`.

Defects found and fixed on the way:
- `publish_state.py` accepted `--run-id` and dropped it. The manifest now carries `run_id`.
- A run's `REPRODUCE.md` replaced the previous run's record and did not keep it. It is now kept in `runs/<run_id>/superseded_case_products/REPRODUCE.md`.
- A failed or stopped run's record quoted the `pp.log` and `tpm_*` of an earlier run. Those are now quoted only if the run wrote them.
- An evaluate-only record reported the verdict of the case's `run.log`, which belongs to an earlier solve. It now names the solve generation it descends from.
- The record said the runner "puts back" `input.inp`. The runner never writes that file.
- `status.py` read every seed line of the current forms as `cold start`.
- A seed given with `SEED=` appended to the previous run's `seed.log`.
- `tests/seed_compatibility.sh` N2 read the case's `output/` rather than its published generation. It failed after the D9s3b re-solve left a seed pair there.

### The two rewritten records

Both were rewritten with the trap's own function (`write_stop_record`), which was handed the values each run held (READ from `run.json`, the logs, the manifests and the campaign status files). Each record says it was written afterwards, and why. The records they replace are kept in `runs/<run_id>/superseded_case_products/REPRODUCE.md`.

- `molecular_scalar_gj1132_kzz1e9/HeH2.13`: stopped by the 6 h ceiling at 18:27:56 in the dtau0 = 1e8 continuation after its outer pass 2. The first solve's 40 passes and the continuation's 2 are tabulated, and the published generation is `g0003_20260919T081651Z_3b86c9ed` (hydrodynamic_refusal).
- `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7`: attempt 2 was stopped by the 3 h ceiling at 19:01:21 after its outer pass 1. The `run.json` of attempt 1 now names its seed and its two solve generations, `g0003` (7826 s from the start of the run) and `g0004` (9629 s).

### The HeH9.7 continuation

Route: `run_campaign.sh 1 8` with one case, `FORCE=1`, `SEED=<case>/states/g0003_20260919T052326Z_ffa78b36`. Naming a molecular solution replaces the atomic-to-molecular conversion. `EXHALE_OUTER_PASSES=90` was set in the environment for this run alone (status `source=environment`), because at the measured contraction of 0.87 a pass the row needed ln(1e-5/4.309e-4)/ln(0.87) = 27 more passes. The wall ceiling was the class's 6 h. Status line: `exit=0 solved budget=molecular_alternation:90p/6h source=environment spent=427s ceiling=none seeds=1/3`.

| pass | hydro info | carrier balance H2 (gated) | at cell | ratio to previous | mass | momentum | energy | s |
|---|---|---|---|---|---|---|---|---|
| 1 | 0 | 4.31E-04 | 263 | -- | 2.31E-09 | 7.14E-14 | 3.06E-08 | 4.22 |
| 2 | 0 | 1.47E-03 | 263 | 3.411 | 1.70E-09 | 1.62E-13 | 1.99E-08 | 76.26 |
| 3 | 0 | 7.25E-04 | 262 | 0.493 | 1.23E-09 | 4.00E-14 | 2.26E-08 | 72.20 |
| 4 | 0 | 3.50E-04 | 262 | 0.483 | 6.59E-09 | 3.81E-13 | 2.55E-08 | 58.10 |
| 5 | 0 | 1.68E-04 | 262 | 0.480 | 1.00E-09 | 1.95E-13 | 2.03E-08 | 57.82 |
| 6 | 0 | 8.05E-05 | 262 | 0.479 | 1.88E-09 | 1.66E-12 | 2.70E-08 | 50.66 |
| 7 | 0 | 3.86E-05 | 262 | 0.480 | 8.69E-10 | 8.41E-13 | 2.08E-08 | 45.94 |
| 8 | 0 | 1.85E-05 | 262 | 0.479 | 5.71E-09 | 3.24E-13 | 1.80E-08 | 21.65 |
| 9 | 0 | 8.82E-06 | 262 | 0.477 | 3.87E-09 | 9.23E-13 | 1.73E-08 | 17.27 |

The warm restart restarts the relaxation controls (MODELS.md section 9.6), so the contraction is 0.48 a pass, not the 0.87 the first solve had reached by its pass 40. The rise at pass 2 is the first relaxation step taken at the restarted controls.

- Solve (the `run.log` final block): info = 0, `||R||` = 1.998e-8, CERTIFIED in the wind (r >= 1.2 R_p). Published as `g0004`, 405 s after the run started.
- Evaluate pass (`pp.log`): NOT CERTIFIED. The only refusing entry is the hydrodynamic mass row, 8.335e-8 at cell 1 against 7.5e-9 ("the rounding anchor of that cell", distance 11.1). The in-run value of the same row at cell 1 was 3.873e-9, a factor 21.5 lower. That is outside the round-trip range quoted for the atomic catalog (factor 0.78 to 1.71, READ, `docs/lhs1140b_stationary_L18_20260915.md`). Why the base cell's mass row does not come back was not isolated. The re-derived lower ghost rows of a molecular restart (MODELS.md section 9.6) are one candidate, and it is untested. Published as `g0005`, which is `latest_complete`.
- log10 Mdot 7.94 and a red-pair EW of 2.4061 %A, both unchanged from the previous published state at the printed digits.
- Figure: `docs/figures/lhs1140b_D9s3c_molecular_HeH9.7_overlay.pdf` (script `<scratchpad>/D9s3c/overlay.py`, a copy of D9s3b's). It plots T, v, 2 n_H2 / n_H and x(H II) against r. "Before" is `g0003`, "after" is `g0004`, on the same radii. Largest relative differences: T 2.6e-4 (r = 3.53), v 1.8e-4, 2 n_H2 / n_H 9.3e-3 (r = 1.07, where it is 3.5e-5), x(H II) 1.1e-3. The certified state is the 40-pass state to these amounts.

### Records written

`LHS1140b/MODELS.md` sections 7 and 8 by `models/status.py --write` (rc 0), and section 9.2 (`run_id`) and the new 9.9 by hand. `campaign_status.txt`: its predecessor was moved to a dated name by the campaign. Nothing was deleted. The first attempt at the two rewritten records left a `not_solved_20260919192839.md` in each directory, which is kept.

## Third pass, under the boundary composition contract (2026-09-20)

Items 6 and the catalog half of 8 of the order in
[PLAN_20260919_rev1.md](PLAN_20260919_rev1.md) section 4, on the binary that
carries the P1 step 2 and step 2b boundary composition contract and the P6b
CFL restriction. Every number in this section is MEASURED here unless it is
labeled READ.

### Verdict

- **Re-evaluation: 83 of 83 CERTIFIED.** Every case whose index named a
  `latest_certified` generation was handed to `models/run_case.sh --evaluate`.
  None refused, so no case of the 83 entered the re-solve list. The wind-window
  mean of rho v r^2 printed by the new certificate is character for character
  the parent's in all 83 (d log10 Mdot = 0 at its six printed digits), the
  catalog's two-decimal log10 Mdot line is unchanged in all 83, and the largest
  |d EW| of the He 10830 red pair is 2.6e-14 %A.
- **The nine molecular `latest_complete` states: 3 CERTIFIED, 6 REFUSED**, the
  same three that turned under the first form of the contract (the P1 step 2
  memo, section 8), and for the same reason: their refusal was the base
  continuity row that the ghost composition reaches. The six that stay refused
  are refused by rows the ghost does not reach.
- **Re-solves: none certified.** The two molecular reference cases and the
  three x0.01-XUV atomic cases each had one budgeted re-solve at their class
  budget and none reached a certificate.
- No figure was made: the overlay of the earlier pass is a before-and-after of
  a molecular reference case that certifies, and neither did.

### The binary of record

`LHS1140b/models/EXHALE_75d55d9d.x`, md5
`75d55d9d4fd0e748cd01d6e35713e34c`. It is a copy (not a rebuild) of the tree's
`EXHALE.x`, built 2026-09-20 06:12:37 +0900. The manifest is
`LHS1140b/models/BINARY_MANIFEST_75d55d9d4fd0.txt`, made the way
`BINARY_MANIFEST_2c3b0acc9983.txt` was (header, then the md5 of the Makefile
and of every production source, the same 169-file list, MEASURED identical file
for file). No file of `src/` or the Makefile is newer than the binary
(MEASURED, `find -newer`); `make -q` clean is READ from the brief of this pass.
Against the 2c3b0acc manifest, nine production sources differ (MEASURED), and
the manifest header names the memo that names each: `ionization_equilibrium.f90`
and `base_boundary.f90` (P1 step 2 and step 2b), `load_IC.f90` (P1 step 2),
`eval_dt.f90` (P6b), `charge_exchange.f90` (P6d),
`binary_element_diffusion.f90` (P6e), `certification.f90` (P5a),
`boundary_state_trace.f90`, and `EXHALE_main.f90`, for which no memo of the
series names a change. Every run of this pass named the binary through
`EXHALE_BIN`.

### How the pass was run

All on `lart4`, `OMP_NUM_THREADS=8`, at most three runs at a time. The 83
evaluations went through `models/run_case.sh --evaluate <case>`, three at a
time, and took 06:17 to 06:31 KST together. The five re-solves went through
`models/run_campaign.sh 1 8` with `CAMPAIGN_LIST` naming one case and
`FORCE=1`, so each carries its class budget, its ceiling and its status line;
no budget override was set, and every status line reads `source=table`. The
nine molecular states were evaluated READ ONLY on scratch copies, because
`--evaluate` resolves `latest_certified` and no molecular case has one: the
generation was copied as the `_IC` pair, the input's relative paths were made
absolute with `Load IC? True` and `Restart intent: stationary evaluate`, and the
binary ran there. Nothing was published for them and the catalog directories
were only read.

### Step 1: the 83 re-evaluations

"d log10 Mdot (window flux)" is log10 of the ratio of the wind-window mean of
rho v r^2 in the new certificate to the one in the parent's
`certification.txt`. The log10 Mdot columns are the catalog's two-decimal line
from `pp.log`, and the EW is the He 10830 red-pair equivalent width from
`tpm_He10830.txt`, read as `status.py` reads it. The parent columns were read
before the pass. The ghost composition contract prints nothing on any of these
83 and applies to none of them: it is entered only where the network is
molecular and carries helium (READ, `ionization_equilibrium.f90`,
`ghost_contract_on = thereis_mol .and. thereis_He`), and all 83 are atomic. So
for each of the 83 the reached accuracy and the application count are NOT
APPLICABLE, not merely unmeasured. After the pass each case's
`latest_certified` is its new evaluate child, 83 of 83.

| case | evaluated parent | evaluate child | verdict | refusing entries | d log10 Mdot (window flux) | log10 Mdot parent -> new | EW [%A] parent -> new | d EW [%A] |
|---|---|---|---|---|---|---|---|---|
| `atomic_photochem_gj1132_kzzprofile/HeH10/k04` | `g0004_20260919T032756Z_093629d3` | `g0005_20260919T212000Z_d6e7e86d` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 2.2224 -> 2.2224 | 1.3e-15 |
| `atomic_photochem_gj1132_kzzprofile/HeH12/k04` | `g0004_20260919T032832Z_be7a184f` | `g0005_20260919T212000Z_3b2b64e9` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 2.2289 -> 2.2289 | -8.9e-16 |
| `atomic_photochem_gj1132_kzzprofile/HeH2.09/k05` | `g0004_20260919T032907Z_a39f11a3` | `g0005_20260919T212000Z_60ce62bd` | CERTIFIED | - | 0.00e+00 | 7.91 -> 7.91 | 1.4765 -> 1.4765 | -1.7e-14 |
| `atomic_photochem_gj1132_kzzprofile/HeH3/k04` | `g0004_20260919T032941Z_5b58b90c` | `g0005_20260919T212034Z_53fb0237` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 1.7640 -> 1.7640 | -2.2e-16 |
| `atomic_photochem_gj1132_kzzprofile/HeH5/k04` | `g0004_20260919T033015Z_0f70dd16` | `g0005_20260919T212036Z_e190429b` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 2.0746 -> 2.0746 | -1.2e-14 |
| `atomic_photochem_gj1132_kzzprofile/HeH7/k04` | `g0004_20260919T033050Z_9e2b0d4a` | `g0005_20260919T212041Z_3b038d48` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 2.1768 -> 2.1768 | 2.2e-15 |
| `atomic_photochem_gj1132_kzzprofile/HeH8/k04` | `g0004_20260919T033125Z_c0d39208` | `g0005_20260919T212108Z_5353f857` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 2.2004 -> 2.2004 | -1.3e-15 |
| `atomic_photochem_gj1132_kzzprofile/HeH9.7/k04` | `g0004_20260919T033200Z_1878856e` | `g0005_20260919T212110Z_3a5ff9c4` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 2.2206 -> 2.2206 | 2.0e-14 |
| `atomic_photochem_gj1132_kzzprofile/HeH9/k04` | `g0004_20260919T033236Z_fa443a57` | `g0005_20260919T212115Z_8f797adc` | CERTIFIED | - | 0.00e+00 | 7.90 -> 7.90 | 2.2137 -> 2.2137 | 1.3e-15 |
| `atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7` | `g0004_20260919T033311Z_52be5f9b` | `g0005_20260919T212143Z_12918db7` | CERTIFIED | - | 0.00e+00 | 6.84 -> 6.84 | 0.2674 -> 0.2674 | -8.9e-16 |
| `atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7` | `g0003_20260919T033344Z_a49f9ac8` | `g0004_20260919T212145Z_fd1d1b0b` | CERTIFIED | - | 0.00e+00 | 7.02 -> 7.02 | 0.4822 -> 0.4822 | 4.2e-15 |
| `atomic_photochem_gj1132x0.20_kzzprofile/HeH9.7` | `g0003_20260919T033418Z_c3c26194` | `g0004_20260919T212150Z_3e09efeb` | CERTIFIED | - | 0.00e+00 | 7.16 -> 7.16 | 0.6711 -> 0.6711 | 7.8e-15 |
| `atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7` | `g0003_20260919T033451Z_e0e2e60b` | `g0004_20260919T212216Z_36c57760` | CERTIFIED | - | 0.00e+00 | 7.26 -> 7.26 | 0.8429 -> 0.8429 | -1.8e-15 |
| `atomic_photochem_gj1132x0.30_kzzprofile/HeH9.7` | `g0004_20260919T033526Z_cb15d8ca` | `g0005_20260919T212218Z_9fad6d96` | CERTIFIED | - | 0.00e+00 | 7.35 -> 7.35 | 1.0097 -> 1.0097 | 6.4e-15 |
| `atomic_photochem_gj1132x0.33_kzzprofile/HeH9.7` | `g0004_20260919T033600Z_261c5da7` | `g0005_20260919T212223Z_7418ae60` | CERTIFIED | - | 0.00e+00 | 7.39 -> 7.39 | 1.0979 -> 1.0979 | 1.6e-15 |
| `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13` | `g0004_20260919T033634Z_85944973` | `g0005_20260919T212249Z_0e9c67c9` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3721 -> 1.3721 | -1.2e-14 |
| `atomic_scalar_gj1132_kzz0/HeH0.55` | `g0004_20260919T033707Z_88eb0d69` | `g0005_20260919T212252Z_c2112514` | CERTIFIED | - | 0.00e+00 | 7.84 -> 7.84 | 0.0305 -> 0.0305 | 0.0e+00 |
| `atomic_scalar_gj1132_kzz0/HeH2.6` | `g0004_20260919T033731Z_df801154` | `g0005_20260919T212257Z_247bee5f` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 0.9549 -> 0.9549 | -5.1e-15 |
| `atomic_scalar_gj1132_kzz0/HeH3.0` | `g0004_20260919T033753Z_5507044b` | `g0005_20260919T212313Z_2f715b51` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.0798 -> 1.0798 | 5.8e-15 |
| `atomic_scalar_gj1132_kzz0/HeH3.5` | `g0004_20260919T033814Z_51875900` | `g0005_20260919T212319Z_f9d81f2a` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2175 -> 1.2175 | -1.2e-14 |
| `atomic_scalar_gj1132_kzz0/HeH3.7` | `g0004_20260919T033836Z_c901733e` | `g0005_20260919T212323Z_c0c152d5` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2698 -> 1.2698 | 6.7e-15 |
| `atomic_scalar_gj1132_kzz0/HeH3.9` | `g0004_20260919T033859Z_79cfb0fc` | `g0005_20260919T212334Z_23e7d463` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3168 -> 1.3168 | -1.3e-15 |
| `atomic_scalar_gj1132_kzz1e10/HeH0.55` | `g0004_20260919T033921Z_4ef7ddb7` | `g0005_20260919T212341Z_6171f455` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 0.6463 -> 0.6463 | 0.0e+00 |
| `atomic_scalar_gj1132_kzz1e10/HeH1.06` | `g0004_20260919T033943Z_35b873ad` | `g0005_20260919T212345Z_bd8e39c9` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.0456 -> 1.0456 | 6.0e-15 |
| `atomic_scalar_gj1132_kzz1e10/HeH1.15` | `g0004_20260919T034004Z_74ef5931` | `g0005_20260919T212356Z_a9b56894` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.1032 -> 1.1032 | -4.9e-15 |
| `atomic_scalar_gj1132_kzz1e10/HeH1.29` | `g0004_20260919T034026Z_58ad4651` | `g0005_20260919T212402Z_6ab5a3b1` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.1867 -> 1.1867 | -9.8e-15 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.55` | `g0004_20260919T034048Z_dbc1f4a8` | `g0005_20260919T212406Z_19186f05` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 0.8065 -> 0.8065 | -3.7e-15 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.795` | `g0004_20260919T034110Z_926dbb2d` | `g0005_20260919T212418Z_f48aa231` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.0302 -> 1.0302 | 4.4e-15 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.865` | `g0004_20260919T034132Z_9af4eb31` | `g0005_20260919T212424Z_6b6af9e3` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.0866 -> 1.0866 | -3.3e-15 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.93` | `g0004_20260919T034155Z_8abe4779` | `g0005_20260919T212428Z_807d6f7e` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.1364 -> 1.1364 | -6.7e-15 |
| `atomic_scalar_gj1132_kzz1e5/HeH2.6` | `g0004_20260919T034217Z_39ae932b` | `g0005_20260919T212440Z_f4141569` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 0.9631 -> 0.9631 | 0.0e+00 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.0` | `g0004_20260919T034239Z_fb323ecf` | `g0005_20260919T212445Z_c105667b` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.0880 -> 1.0880 | -4.4e-16 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.35` | `g0004_20260919T034301Z_99915b3f` | `g0005_20260919T212449Z_6ba31c56` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1862 -> 1.1862 | 3.1e-15 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.64` | `g0004_20260919T034323Z_3b108aa2` | `g0005_20260919T212501Z_03986912` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2627 -> 1.2627 | 9.3e-15 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.93` | `g0004_20260919T034345Z_7b4d3899` | `g0005_20260919T212507Z_88f69b03` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3314 -> 1.3314 | -2.2e-15 |
| `atomic_scalar_gj1132_kzz1e6/HeH0.55` | `g0004_20260919T034407Z_30c61e49` | `g0005_20260919T212511Z_628c1387` | CERTIFIED | - | 0.00e+00 | 7.84 -> 7.84 | 0.0591 -> 0.0591 | 0.0e+00 |
| `atomic_scalar_gj1132_kzz1e6/HeH2.4` | `g0004_20260919T034429Z_e43909f6` | `g0005_20260919T212523Z_ebcd0771` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 0.9465 -> 0.9465 | -8.9e-16 |
| `atomic_scalar_gj1132_kzz1e6/HeH2.8` | `g0004_20260919T034451Z_5798e5d8` | `g0005_20260919T212529Z_045648ae` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.0804 -> 1.0804 | 8.4e-15 |
| `atomic_scalar_gj1132_kzz1e6/HeH3.19` | `g0004_20260919T034513Z_31fca486` | `g0005_20260919T212533Z_24b9c7c5` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1949 -> 1.1949 | 9.1e-15 |
| `atomic_scalar_gj1132_kzz1e6/HeH3.46` | `g0004_20260919T034535Z_f780d699` | `g0005_20260919T212544Z_2cd6945f` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2669 -> 1.2669 | 2.6e-14 |
| `atomic_scalar_gj1132_kzz1e6/HeH3.74` | `g0004_20260919T034558Z_6c7bd82d` | `g0005_20260919T212550Z_a2d026d7` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3388 -> 1.3388 | -1.1e-15 |
| `atomic_scalar_gj1132_kzz1e7/HeH0.55` | `g0004_20260919T034620Z_553ec42c` | `g0005_20260919T212554Z_6ede98a7` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 0.1585 -> 0.1585 | 4.4e-16 |
| `atomic_scalar_gj1132_kzz1e7/HeH2.70` | `g0004_20260919T034642Z_e5b7019c` | `g0005_20260919T212606Z_ec322f94` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1996 -> 1.1996 | -3.8e-15 |
| `atomic_scalar_gj1132_kzz1e7/HeH2.94` | `g0004_20260919T034704Z_11578937` | `g0005_20260919T212612Z_9be703e3` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2724 -> 1.2724 | 5.1e-15 |
| `atomic_scalar_gj1132_kzz1e7/HeH3.18` | `g0004_20260919T034726Z_9028a211` | `g0005_20260919T212616Z_585264ab` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3265 -> 1.3265 | 7.5e-15 |
| `atomic_scalar_gj1132_kzz1e8/HeH0.55` | `g0004_20260919T034749Z_fd533291` | `g0005_20260919T212628Z_0c8c8c74` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 0.3089 -> 0.3089 | 4.4e-16 |
| `atomic_scalar_gj1132_kzz1e8/HeH2.05` | `g0004_20260919T034810Z_36bd393a` | `g0005_20260919T212634Z_1fbfb43a` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1513 -> 1.1513 | 1.8e-15 |
| `atomic_scalar_gj1132_kzz1e8/HeH2.23` | `g0004_20260919T034832Z_f2bb678b` | `g0005_20260919T212638Z_fb6f4c53` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2206 -> 1.2206 | 7.8e-15 |
| `atomic_scalar_gj1132_kzz1e8/HeH2.41` | `g0004_20260919T034854Z_9a5c1f0c` | `g0005_20260919T212649Z_829714a0` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.2826 -> 1.2826 | 1.6e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH0.55` | `g0004_20260919T034917Z_098da21a` | `g0005_20260919T211716Z_afd550b8` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 0.4793 -> 0.4793 | -3.3e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH1.50` | `g0004_20260919T034938Z_4274fa5a` | `g0005_20260919T212656Z_f46c0b98` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.0996 -> 1.0996 | -8.0e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH1.60` | `g0004_20260919T035000Z_a0851042` | `g0005_20260919T212700Z_8b6147e6` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1511 -> 1.1511 | -5.6e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH1.70` | `g0004_20260919T035022Z_fc02686a` | `g0005_20260919T212711Z_4fc99db6` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.1970 -> 1.1970 | 1.8e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH2.13` | `g0004_20260919T032647Z_aece8741` | `g0005_20260919T212717Z_f9d878de` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.3749 -> 1.3749 | -3.3e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH4.0` | `g0004_20260919T035044Z_acba83cc` | `g0005_20260919T212721Z_3c042f1f` | CERTIFIED | - | 0.00e+00 | 7.87 -> 7.87 | 1.8229 -> 1.8229 | -4.2e-15 |
| `atomic_scalar_gj1132_kzz1e9/HeH9.7` | `g0004_20260919T032647Z_9ca293a8` | `g0005_20260919T212733Z_8b765d17` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 2.0592 -> 2.0592 | -9.8e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH0.083` | `g0004_20260919T035107Z_5338d0c2` | `g0005_20260919T212739Z_6140ce16` | CERTIFIED | - | 0.00e+00 | 7.84 -> 7.84 | 0.3797 -> 0.3797 | 0.0e+00 |
| `atomic_scalar_gj1132_wellmixed/HeH0.40` | `g0004_20260919T035130Z_bec4f54e` | `g0005_20260919T212744Z_7ea7c437` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.1011 -> 1.1011 | 1.6e-14 |
| `atomic_scalar_gj1132_wellmixed/HeH0.42` | `g0004_20260919T035152Z_ab79b2da` | `g0005_20260919T212756Z_346736f9` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.1324 -> 1.1324 | 9.5e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH0.44` | `g0004_20260919T035216Z_ee409f97` | `g0005_20260919T212802Z_af715748` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.1626 -> 1.1626 | 4.4e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH0.55` | `g0004_20260919T035238Z_43fb8ba9` | `g0005_20260919T212806Z_3ffa0560` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.3118 -> 1.3118 | -3.1e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH1` | `g0004_20260919T035301Z_ecc54ed6` | `g0005_20260919T212818Z_6b9c5642` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.6959 -> 1.6959 | 5.3e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH10` | `g0004_20260919T035323Z_f3e7ada1` | `g0005_20260919T212824Z_d27cedc6` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.9347 -> 1.9347 | -2.2e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH100` | `g0004_20260919T035347Z_090396b9` | `g0005_20260919T212828Z_46dfae9a` | CERTIFIED | - | 0.00e+00 | 7.86 -> 7.86 | 1.9644 -> 1.9644 | -4.4e-15 |
| `atomic_scalar_gj1132_wellmixed/HeH1000` | `g0004_20260919T035411Z_6b8ed8aa` | `g0005_20260919T212840Z_6b628859` | CERTIFIED | - | 0.00e+00 | 7.85 -> 7.85 | 1.9197 -> 1.9197 | 1.3e-14 |
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13` | `g0004_20260919T035436Z_89f3eb8e` | `g0005_20260919T212847Z_80bee6f1` | CERTIFIED | - | 0.00e+00 | 6.82 -> 6.82 | 0.0001 -> 0.0001 | 0.0e+00 |
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7` | `g0003_20260919T035456Z_661ac065` | `g0004_20260919T212852Z_de543034` | CERTIFIED | - | 0.00e+00 | 6.81 -> 6.81 | 0.2353 -> 0.2353 | -1.1e-15 |
| `atomic_scalar_gj1132x0.15_kzz1e9/HeH2.13` | `g0004_20260919T035518Z_35541992` | `g0005_20260919T212904Z_2036e2ba` | CERTIFIED | - | 0.00e+00 | 7.00 -> 7.00 | 0.0003 -> 0.0003 | 0.0e+00 |
| `atomic_scalar_gj1132x0.15_kzz1e9/HeH9.7` | `g0004_20260919T035539Z_9860d32b` | `g0005_20260919T212908Z_8d291422` | CERTIFIED | - | 0.00e+00 | 7.00 -> 7.00 | 0.4514 -> 0.4514 | -3.2e-15 |
| `atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13` | `g0004_20260919T035601Z_c2dce3ee` | `g0005_20260919T212913Z_1ece57d5` | CERTIFIED | - | 0.00e+00 | 7.13 -> 7.13 | 0.0022 -> 0.0022 | 0.0e+00 |
| `atomic_scalar_gj1132x0.20_kzz1e9/HeH9.7` | `g0004_20260919T035622Z_64f988b8` | `g0005_20260919T212924Z_78579b62` | CERTIFIED | - | 0.00e+00 | 7.13 -> 7.13 | 0.6294 -> 0.6294 | -2.7e-15 |
| `atomic_scalar_gj1132x0.25_kzz1e9/HeH2.13` | `g0004_20260919T035644Z_508d351b` | `g0005_20260919T212930Z_8c87e392` | CERTIFIED | - | 0.00e+00 | 7.24 -> 7.24 | 0.0700 -> 0.0700 | 0.0e+00 |
| `atomic_scalar_gj1132x0.25_kzz1e9/HeH9.7` | `g0004_20260919T035705Z_c53efa85` | `g0005_20260919T212934Z_5d59402e` | CERTIFIED | - | 0.00e+00 | 7.23 -> 7.23 | 0.7913 -> 0.7913 | -6.0e-15 |
| `atomic_scalar_gj1132x0.30_kzz1e9/HeH2.13` | `g0004_20260919T035728Z_007ecf98` | `g0005_20260919T212946Z_99b5de22` | CERTIFIED | - | 0.00e+00 | 7.32 -> 7.32 | 0.1894 -> 0.1894 | -8.9e-16 |
| `atomic_scalar_gj1132x0.30_kzz1e9/HeH9.7` | `g0004_20260919T035749Z_e6f327b3` | `g0005_20260919T212950Z_7688a5d5` | CERTIFIED | - | 0.00e+00 | 7.31 -> 7.31 | 0.9395 -> 0.9395 | 7.8e-15 |
| `atomic_scalar_gj1132x0.33_kzz1e9/HeH2.13` | `g0004_20260919T035812Z_76873a9d` | `g0005_20260919T212956Z_bab6dc42` | CERTIFIED | - | 0.00e+00 | 7.36 -> 7.36 | 0.2555 -> 0.2555 | -3.3e-15 |
| `atomic_scalar_gj1132x0.33_kzz1e9/HeH9.7` | `g0004_20260919T035834Z_14c7d6b1` | `g0005_20260919T213007Z_ddde9ea9` | CERTIFIED | - | 0.00e+00 | 7.36 -> 7.36 | 1.0231 -> 1.0231 | 5.3e-15 |
| `atomic_scalar_gj699_wellmixed/HeH0.042` | `g0004_20260919T035856Z_f07b2ff2` | `g0005_20260919T213013Z_7bfbfed6` | CERTIFIED | - | 0.00e+00 | 8.57 -> 8.57 | 0.9091 -> 0.9091 | 0.0e+00 |
| `atomic_scalar_gj699_wellmixed/HeH0.046` | `g0004_20260919T035920Z_189e1b29` | `g0005_20260919T213017Z_b3e606eb` | CERTIFIED | - | 0.00e+00 | 8.57 -> 8.57 | 0.9857 -> 0.9857 | 1.0e-14 |
| `atomic_scalar_gj699_wellmixed/HeH0.050` | `g0004_20260919T035944Z_5b6a0b9f` | `g0005_20260919T213030Z_0c93142c` | CERTIFIED | - | 0.00e+00 | 8.57 -> 8.57 | 1.0608 -> 1.0608 | -2.2e-15 |
| `atomic_scalar_gj699_wellmixed/HeH0.083` | `g0004_20260919T040008Z_ae236f2e` | `g0005_20260919T213036Z_621ab39d` | CERTIFIED | - | 0.00e+00 | 8.57 -> 8.57 | 1.6271 -> 1.6271 | -2.2e-15 |
| `atomic_scalar_gj699_wellmixed/HeH1` | `g0004_20260919T040031Z_2091c6ca` | `g0005_20260919T213040Z_8f23234c` | CERTIFIED | - | 0.00e+00 | 8.55 -> 8.55 | 4.5250 -> 4.5250 | -8.9e-15 |
| `atomic_scalar_gj699_wellmixed/HeH1000` | `g0004_20260919T040054Z_b13d83c8` | `g0005_20260919T213053Z_29acb058` | CERTIFIED | - | 0.00e+00 | 8.42 -> 8.42 | 2.8725 -> 2.8725 | 3.1e-15 |

### Step 2: the nine molecular `latest_complete` states

None of the nine has a `latest_certified`, so each was evaluated READ ONLY on a
scratch copy of its `latest_complete` state (MEASURED: of the 137
`state_index.json` files under `LHS1140b/models`, nine name a molecular case,
seven of them the case directories and two the `pre_D9s3` archives, and all
nine carry `latest_certified: null`). Nothing was published and no catalog
directory was written to. The last two columns are the ghost composition
contract's own report line, as the boundary model prints it at the end of the
evaluation: the last application's move in the ghost's species densities, in the
composition's own units, against the 1e-11 of step 2b, and the move in the
ghost's heavy-particle plus electron count against 1e-06.

| case | `latest_complete` evaluated | verdict | refusing entries | ghost composition move, against 1e-11 | its count move, against 1e-06 | applications |
|---|---|---|---|---|---|---|
| `molecular_photochem_gj1132_kzzprofile/HeH9` | `g0003_20260917T222518Z_c55563e3` | NOT CERTIFIED | hydrodynamic mass row: row measure  1.041E+00 above  3.4E-10 at cell 2; hydrodynamic momentum row: row measure  1.450E-03 above  1.0E-08 at cell 1; hydrodynamic energy row: row measure  1.021E+00 above  1.0E-06 at cell 1; carrier balance H2: gated row measure  1.000E+00 above  1.0E-05 at cell 432 (a wind cell); elemental transport He/H partition: gated row measure  1.494E-04 above  1.0E-05 at cell 217 (a wind cell) | 5.84538E-15 | 2.84553E-14 | 3 |
| `molecular_scalar_gj1132_kzz1e9/HeH0.083` | `g0002_20260916T220756Z_f9fca548` | NOT CERTIFIED | carrier balance H2: gated row measure  8.150E-04 above  1.0E-05 at cell 499 (a wind cell); elemental transport He/H partition: gated row measure  2.503E-04 above  1.0E-05 at cell 280 (a wind cell) | 3.62083E-12 | 4.13900E-12 | 2 |
| `molecular_scalar_gj1132_kzz1e9/HeH0.55` | `g0003_20260916T025746Z_47afb99c` | NOT CERTIFIED | hydrodynamic energy row: row measure  4.634E-01 above  1.0E-06 at cell 257; carrier balance H2: gated row measure  8.922E-03 above  1.0E-05 at cell 261 (a wind cell) | 4.15451E-12 | 1.10343E-11 | 2 |
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | `g0003_20260919T081651Z_3b86c9ed` | NOT CERTIFIED | hydrodynamic mass row: row measure  3.623E-05 above  3.0E-12 at cell 287; hydrodynamic momentum row: row measure  4.036E-04 above  1.0E-08 at cell 500; hydrodynamic energy row: row measure  9.565E-02 above  1.0E-06 at cell 6; carrier balance H2: gated row measure  6.786E-02 above  1.0E-05 at cell 227 (a wind cell); elemental transport He/H partition: gated row measure  7.288E-04 above  1.0E-05 at cell 217 (a wind cell) | 1.13196E-13 | 6.13143E-13 | 3 |
| `molecular_scalar_gj1132_kzz1e9/HeH2.13/pre_D9s3` | `g0002_20260917T174914Z_19d7ee4d` | CERTIFIED | - | 8.64453E-15 | 3.92457E-14 | 3 |
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | `g0005_20260919T103629Z_af3b8d33` | CERTIFIED | - | 6.84134E-14 | 4.00014E-13 | 3 |
| `molecular_scalar_gj1132_kzz1e9/HeH9.7/pre_D9s3` | `g0002_20260917T181637Z_552403ab` | CERTIFIED | - | 6.84134E-14 | 4.00258E-13 | 3 |
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | `g0002_20260918T034603Z_258bc5d7` | NOT CERTIFIED | carrier balance H2: gated row measure  7.371E-02 above  1.0E-05 at cell 500 (a wind cell) | 5.24352E-12 | 6.87245E-12 | 2 |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | `g0002_20260915T221824Z_d7347336` | NOT CERTIFIED | hydrodynamic energy row: row measure  1.000E+00 above  1.0E-06 at cell 246; carrier balance H2: gated row measure  2.363E-02 above  1.0E-05 at cell 306 (a wind cell) | 1.38561E-12 | 3.63096E-12 | 2 |

**Three of the nine certify, and they are the same three.** They are the three
the first form of the contract turned (the P1 step 2 memo, section 8), and the
verdicts here reproduce that table case for case. What step 2b changed is the
norm the move is taken in, not which states close: every one of the nine
reaches its fixed point well inside the contract, in two or three applications,
at 1e-15 to 5e-12 against 1e-11, so the contract is not what separates the
three from the six. The six are refused by rows the ghost composition does not
reach: a carrier balance H2 in the wind in all six, an interior hydrodynamic
energy row in three, and in the photochemical HeH9 the base rows at order one.

### Step 3: the two molecular reference cases

Route: `models/run_case.sh` through the campaign, molecular branch,
`molecular_alternation` budget (40 passes, 6 h). `models/pick_seed.py` offers no
compatible molecular state for either, so the runner builds the seed from the
atomic case of the same name and He/H, using that case's `latest_certified`,
which is the step-1 evaluate child this binary certified.

| case | seed | seed conversion | solve 1 (dtau0 = 1) | continuation (dtau0 = 1e8) | ending class | published | wall |
|---|---|---|---|---|---|---|---|
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | `atomic_scalar_gj1132_kzz1e9/HeH2.13`, `g0005_20260919T212847Z_...` converted by the binary | H nuclei 4.00e-16, mass 3.75e-16, equation-of-state closure 3.40e-16 | 40 passes, `hydrodynamic_refusal`: mass row 3.868E-05 above 3.0E-12 at cell 285; gated carrier balance H2 7.75e-2 at cell 224 at pass 40; hydro info = 2 at pass 40 | 2 passes, info = 2 then info = 0 with the carrier row 2.07e-1 at cell 224; stopped by the wall ceiling in pass 3 | `stopped_by_wall_ceiling`, UNCERTIFIED INCOMPLETE | `g0004_20260920T025822Z_7033d8a4` (solve 1, NOT CERTIFIED; kept whole in `solve_2/`) | 6 h 00 m |
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | `atomic_scalar_gj1132_kzz1e9/HeH9.7`, `g0005_20260919T212733Z_8b765d17` converted by the binary | H nuclei 2.21e-16, mass 3.47e-16, equation-of-state closure 3.26e-16 | 40 passes, hydro info = 0 from pass 8 on; gated carrier balance H2 1.32e-2 (pass 17), 3.16e-3 (27), 6.78e-4 (38), 5.13e-4 at the ceiling, at cell 262 | not taken (the continuation addresses `hydrodynamic_refusal` only) | `composition_refusal` at the pass ceiling, UNCERTIFIED INCOMPLETE | `g0006_20260920T000506Z_5b2061aa` (NOT CERTIFIED) | 2 h 32 m |

Neither certified, so neither was evaluated and neither has a before-and-after
overlay figure. `HeH9.7` is the case that converges and does not arrive: its
gated carrier row falls monotonically through the run and stands at 5.13e-4
against 1e-05 when the 40 passes are spent, which is where the second pass left
it (4.31e-4) within a factor 1.2. `HeH2.13` does not converge: its carrier row
bands between 5e-2 and 2e-1 from pass 20 to the ceiling with no trend, and its
hydrodynamic rows stay at 1e-2 to 1.

Note the two `latest_complete` states these cases hold are now the states this
pass wrote, so the step-2 rows above and these are about different bytes: the
step-2 evaluation of `HeH9.7` read `g0005`, the generation the 2026-09-19
continuation published, and CERTIFIED it; the solve of step 3 then published
`g0006`, which its own certificate refuses. A refusal of a solve's own state is
not a refused re-evaluation and takes no certified reference away, because
`HeH9.7` has none.

### Step 4: the three x0.01-XUV atomic cases

One re-solve each at the `atomic_prescribed` class budget, 40 passes and 30 min,
with no override: the six-times budget of the second pass had already refuted
the "more time" reading (READ, the second pass above), so this pass measures
where each stops at its own class budget and no further. `pick_seed.py` gives
each the case at 0.10 of the XUV, tier 3, whose certified state is this pass's
own step-1 evaluate child.

| case | seed (pick_seed) | what ran | where it stopped | the rows at the last pass |
|---|---|---|---|---|
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | tier 3, `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13/states/g0005_20260919T212847Z_80bee6f1` | one solve, 5 outer passes complete | wall ceiling (30 m) in the wind pass after pass 5, `stopped_by_wall_ceiling`, UNCERTIFIED INCOMPLETE; nothing new published (`latest_complete` still `g0002_20260919T004834Z_d59ec9ff`) | hydro info = 2 at every pass; mass 2.46e-1, momentum 3.29e-7, energy 1.22e-1, unmoved from pass 4 to 5; gated species row 4.38e-2 at cell 331 |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | tier 3, `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7/states/g0004_20260919T212852Z_de543034` | one solve, 4 outer passes complete | wall ceiling in the wind pass after pass 4, same class; nothing new published (`latest_complete` still `g0002_20260919T004843Z_f7485b14`) | hydro info = 2 at every pass; energy row 9.98e-1 at the base, mass 7.33e-8, momentum 8.25e-14; gated species row 1.3e-2 down to 9.77e-4 at cell 354 |
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | tier 3, `atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7/states/g0005_20260919T212143Z_12918db7` | one solve, 4 outer passes complete | wall ceiling in the wind pass after pass 4, same class; nothing new published (`latest_complete` still `g0004_20260919T094145Z_259fe9c3`) | hydro info = 2 at every pass; mass 1.16, energy 1.10, unmoved from pass 3 to 4; gated species row 3.29e-2 at cell 387 |

The finding of the second pass stands under the new boundary: in all three the
composition relaxes and the hydrodynamic rows do not move at all between
passes, at the base for the two `HeH9.7` cases and in the interior mass row for
`HeH2.13`. A ceiling is not a solution, and none of the three is called one.

### Records written

- `LHS1140b/MODELS.md` sections 7 and 8, by `models/status.py --write` (rc 0).
  Section 7 moves five rows: the three x0.01 cases and
  `molecular_scalar_gj1132_kzz1e9/HeH2.13` now carry this pass's
  `stopped_by_wall_ceiling` reason, and
  `molecular_scalar_gj1132_kzz1e9/HeH9.7` moves from `evaluated uncertified` to
  `info=1 uncertified composition_refusal` with its residual 2.58e-08. Section 8
  moves the 83 evaluated rows onto their new generations.
- `models/CASE_INVENTORY.md` and `models/CASE_INVENTORY.json`, by
  `models/status.py inventory --write`: MEASURED 2026-09-20 12:33:51 KST, 137
  state indexes and 83 certified, by kind 84 ladder cases (74 certified), 46
  closure rungs (9), 4 `diffusion_check` (0), 2 archives (0), 1 study directory
  (0). The counts are the same as the 2026-09-19 inventory's.
- New generations: 83 `evaluate` children and the two molecular solve
  generations. The three x0.01 cases published nothing: each was stopped in a
  pass that had written no complete state of its own. Nothing was deleted, and
  every file an evaluation replaced in a case directory is in
  `runs/<run>/superseded_case_products/`.
- `REPRODUCE.md` of the 83 evaluated cases and of the two molecular cases.
- `models/campaign_status.txt` and its dated predecessors
  `campaign_status_20260920{070235,073235,080236,090507}.txt`: one per campaign
  of this pass.
- `LHS1140b/models/EXHALE_75d55d9d.x` and
  `LHS1140b/models/BINARY_MANIFEST_75d55d9d4fd0.txt`.

### Noticed, not changed

1. `models/status.py` carries the binary of record as a constant,
   `BINARY_OF_RECORD = 'EXHALE_2c3b0acc.x'`, with no environment override, so
   the inventory this pass wrote names the PREVIOUS binary in its header
   although every run of this pass used `EXHALE_75d55d9d.x`. The counts the
   inventory reports are unaffected, since they are read from the indexes. The
   source of `models/` is frozen for this pass, so the constant was left alone.
2. The evaluate route resolves `latest_certified` and stops where there is
   none, so no molecular case of the catalog can be measured through
   `run_case.sh --evaluate` at all; the nine states of step 2 had to be
   evaluated outside the runner, as the P1 step 2 memo also had to.


## 2026-09-20, fourth pass: the nine molecular states measured through the runner

A certificate is a property of the state, measured: every active equation of
the certification inventory evaluated on the bytes as written and within its
tolerance, under a named binary. The third pass above could not record such a
measurement for any molecular case, because `run_case.sh --evaluate` resolved
`latest_certified` and stopped where the index named none, and a measurement
made outside the runner publishes nothing (the third pass listed this as its
second item under "Noticed, not changed"). That was a defect of the
procedure and not of the states. The entry and the publisher were changed
(`LHS1140b/MODELS.md` section 9.15), and the nine states were measured again,
this time through the runner and the publisher.

Binary of record `EXHALE_75d55d9d.x`, md5 `75d55d9d4fd0e748cd01d6e35713e34c`,
8 threads, on lart4. Every number below is MEASURED in this pass.

**Three of the nine certify, the same three as the third pass, and each now
carries a certified generation.** The catalog goes from 83 certified cases of
137 state indexes to 85; the two are `molecular_scalar_gj1132_kzz1e9/HeH2.13`
and `.../HeH9.7`, the first certified molecular cases of the catalog. The two
`pre_D9s3` archival copies were NOT run and NOT repaired: they cannot be run
where they sit (their `NOTE.md`), and their states are, MEASURED by md5, the
generations `g0002_20260917T174914Z_19d7ee4d` and
`g0002_20260917T181637Z_552403ab` of the two cases one level up, which is
where they were measured from.

| case | generation measured | how it was reached | verdict | refusing entries | ghost composition, against 1e-11 | its count move, against 1e-06 | applications | log10 Mdot | He 10830 EW [%A] |
|---|---|---|---|---|---|---|---|---|---|
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | `g0004_20260919T103629Z_9b8394a3` | `--generation` | CERTIFIED | none | 6.84134E-14 | 4.00014E-13 | 3 | 7.94 | 2.4061 |
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | `g0002_20260917T181637Z_552403ab`, the `pre_D9s3` state | `--generation` | CERTIFIED | none | 6.84134E-14 | 4.00258E-13 | 3 | 7.94 | 2.4061 |
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | `g0002_20260917T174914Z_19d7ee4d`, the `pre_D9s3` state | `--generation` | CERTIFIED | none | 8.64453E-15 | 3.92457E-14 | 3 | 7.90 | 1.5603 |
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | `g0004_20260920T025822Z_7033d8a4` | `latest_complete` | NOT CERTIFIED | mass 3.868E-05 above 3.0E-12 at cell 285; momentum 3.848E-04 above 1.0E-08 at cell 500; energy 7.977E-01 above 1.0E-06 at cell 2; carrier balance H2 7.747E-02 above 1.0E-05 at cell 224 (a wind cell); elemental transport He/H 1.731E-04 above 1.0E-05 at cell 217 (a wind cell) | 1.09344E-13 | 5.92131E-13 | 3 | 7.80 | 1.3078 |
| `molecular_photochem_gj1132_kzzprofile/HeH9` | `g0003_20260917T222518Z_c55563e3` | `latest_complete` | NOT CERTIFIED | mass 1.041E+00 above 3.4E-10 at cell 2; momentum 1.450E-03 above 1.0E-08 at cell 1; energy 1.021E+00 above 1.0E-06 at cell 1; carrier balance H2 1.000E+00 above 1.0E-05 at cell 432 (a wind cell); elemental transport He/H 1.494E-04 above 1.0E-05 at cell 217 (a wind cell) | 5.84538E-15 | 2.84553E-14 | 3 | 7.93 | 2.3332 |
| `molecular_scalar_gj1132_kzz1e9/HeH0.083` | `g0002_20260916T220756Z_f9fca548` | `latest_complete` | NOT CERTIFIED | carrier balance H2 8.150E-04 above 1.0E-05 at cell 499 (a wind cell); elemental transport He/H 2.503E-04 above 1.0E-05 at cell 280 (a wind cell) | 3.62083E-12 | 4.13900E-12 | 2 | 7.39 | 0.0003 |
| `molecular_scalar_gj1132_kzz1e9/HeH0.55` | `g0003_20260916T025746Z_47afb99c` | `latest_complete` | NOT CERTIFIED | energy 4.634E-01 above 1.0E-06 at cell 257; carrier balance H2 8.922E-03 above 1.0E-05 at cell 261 (a wind cell) | 4.15451E-12 | 1.10343E-11 | 2 | 7.64 | 0.1591 |
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | `g0002_20260918T034603Z_258bc5d7` | `latest_complete` | NOT CERTIFIED | carrier balance H2 7.371E-02 above 1.0E-05 at cell 500 (a wind cell) | 5.24352E-12 | 6.87245E-12 | 2 | 6.98 | 0.1164 |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | `g0002_20260915T221824Z_d7347336` | `latest_complete` | NOT CERTIFIED | energy 1.000E+00 above 1.0E-06 at cell 246; carrier balance H2 2.363E-02 above 1.0E-05 at cell 306 (a wind cell) | 1.38561E-12 | 3.63096E-12 | 2 | 7.41 | 0.6807 |

The verdicts, the refusing entries and the ghost contract numbers reproduce
the third pass row for row, with one difference that is not a difference of
state: `HeH2.13` was measured here at its own newer `latest_complete`
(`g0004`, the 2026-09-20 solve) rather than at the `g0003` the third pass
read, and that state is refused by five entries instead of five at other
cells. Every one of the nine reaches the ghost composition fixed point in two
or three applications, at 6e-15 to 5e-12 against 1e-11, so the contract is
again not what separates the three from the six: the six are refused by rows
the ghost does not reach.

`molecular_photochem_gj1132_kzzprofile/HeH2.09` and
`molecular_scalar_gj1132_wellmixed/HeH2.13` publish no index and are not
among the nine; they hold no state.

### A defect this pass found and fixed: a refusal was recorded as `unknown`

The first six evaluations of this pass published their refusal with
`certification.status: unknown`. MEASURED in their logs: the evaluate route
prints its measurement under the heading
`(certification) work state of the loaded restart`, and, where it refuses the
state, adds a further `(certification)` heading at the end of the run saying
how the pair was written (`the state written is a relaxation snapshot ...`),
which measures nothing. `publish_state.certification_from_log` took the LAST
heading, read that note as the certificate, and recorded no verdict. The
consequence is not cosmetic: `demote_refused_parent` fires on
`status == 'NOT CERTIFIED'`, so a refused re-evaluation of a certified state
by this binary would not have taken the certified reference away.

The reader now takes the last heading whose own block states CERTIFIED or NOT
CERTIFIED, preferring the two headings that name the state as written, and a
heading's block ends at the next heading. The six evaluations were then run
again against the same generations, by name, and their refusals are on record;
nothing was deleted, and the first six evaluate generations stay published
beside them. Check `C9` of `models/tests/certified_evaluation.sh` covers it.

### Records written

- `LHS1140b/MODELS.md` section 9.15 (the rule), and sections 7 and 8
  regenerated by `models/status.py --write` (rc 0): the seven molecular rows
  move onto their evaluations, `kzz1e9/HeH2.13` and `kzz1e9/HeH9.7` now read
  `evaluated certified`, and each of the seven carries a log10 Mdot and an
  equivalent width where it carried none.
- `models/CASE_INVENTORY.md` and `models/CASE_INVENTORY.json`, by
  `models/status.py inventory --write`: MEASURED 2026-09-20 14:42:33 KST, 137
  state indexes and 85 certified, by kind 84 ladder cases (76 certified), 46
  closure rungs (9), 4 `diffusion_check` (0), 2 archives (0), 1 study
  directory (0). Its `certified under` column is now `what has been measured`
  and reads `certified by measurement on <date>, binary <md5>` or `not
  measured under binary <md5>`.
- 16 new `evaluate` generations over the seven cases, `REPRODUCE.md` of each,
  and the products each evaluation wrote. Nothing was deleted, and every
  case-level file an evaluation replaced is in
  `runs/<run>/superseded_case_products/`.
