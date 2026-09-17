# Reproducing the LHS 1140 b model tree

Everything under `models/` is produced by the four scripts beside this file.
`MODELS.md` in the directory above states what each model IS -- the groups,
the composition ladders, the recipe and why it is that recipe; this file
states what to run, in what order, to get the tree from nothing.

Each case also carries its own `REPRODUCE.md`, written by the run that
produced it: the binary and its md5, the state that seeded it, the
commands with their environment, the outer passes and the certification, and
the clock. That file is the record of one case; this one is the record of
the tree.

## What has to be installed

| | |
|---|---|
| Fortran compiler | the `gfortran` on PATH; the binary these models were solved with was built by conda-forge gcc 16.2.0 (`readelf -p .comment EXHALE.x` says which one built any given binary) |
| linear algebra | OpenBLAS, through the LAPACK of the compiler's own prefix, rpath'd by the makefile. The runs pin the BLAS pool to one thread themselves |
| Python | the system interpreter, with `numpy` and `matplotlib`, for the runners, the readers and the transit synthesis |
| Photochem | only for the `*_photochem_*` closure rungs: Photochem 0.9.0 in the interpreter `/opt/miniconda3/bin/python3`, which each rung's `closure.json` names so the rung reproduces on the build it was run on |
| the instrument kernel | `LHS1140b/winered_hires_y.sh`, sourced by both runners: the WINERED HIRES-Y resolving power (R = 68000) of the Cherubim et al. (2026) measurement |

Build the binary from the repository root:

```bash
make                      # -> EXHALE.x
```

## The tree, from nothing

```bash
cd LHS1140b

# 1. write every case directory from the group table of MODELS.md section 3.
#    A case that already carries output/Hydro_ioniz.txt is left alone;
#    --force overwrites it.
python3 models/make_models.py            # 88 cases: 73 prescribed, 15 rungs

# 2. the target grid.  Every case is solved on the same cell centers, and
#    models/current_grid_Hydro_ioniz.txt names them: it is a Hydro_ioniz.txt
#    written by the current code on this planet.  Replace it only when the
#    grid construction changes, and re-solve when you do.

# 3. one prescribed-composition case, seed to line
OMP_NUM_THREADS=8 models/run_case.sh atomic_scalar_gj1132_kzz1e9/HeH1.60

# 4. one flux-closure rung (Photochem + EXHALE alternating)
OMP_NUM_THREADS=8 models/run_closure.sh atomic_photochem_gj1132_kzzprofile/HeH9

# 5. read the tree back and write MODELS.md sections 7 and 8
python3 models/status.py --write
```

`run_case.sh` and `run_closure.sh` each take one `<group>/<case>`, do
everything that case needs, and print one line saying what came out. They
skip a case that already carries `tpm_He10830_metrics.txt` unless `FORCE=1`.
A case can always be run on its own: where no case of its group is certified
yet it seeds itself from the archive. The ORDER nevertheless decides which
seed a case gets, because a certified case of the same group is preferred to
any archived state, so a ladder run from its middle outward takes fewer
outer passes than one run from the archive at every rung. The answer is the
same either way (below).

## The scripts

| | |
|---|---|
| `make_models.py` | writes `input.inp` (and `base.inp`, or the profile, or `input_template.inp` + `closure.json`) for every case of the group table. One key is written in one place, so a key that is wrong is wrong once |
| `pick_seed.py` | names the state that seeds one case: this same case's own state of 2026-09-14 in `models_20260914_preL21/` first (`tier0`), then a CERTIFIED case of the same group, then a certified case of another group with the same physics, then the archive (newest generation first); within a tier the same spectrum, lower boundary, diffusion and `K_zz`, nearest in log He/H. The tier is on the line it prints |
| `run_case.sh` | one prescribed-composition case: seed, wind, advection-corrected profiles, transit spectrum, record. A MOLECULAR case is seeded from the certified atomic case of the same name and He/H; `SEED_X2` says what H2 the conversion carries, and its default `local` gives every cell the smaller of the thermochemical fit and the root of that cell's own H2 carrier row, the base layer included; the root is solved for, the row being quadratic in n(H2) (`docs/input_schema.md` appendix D) |
| `run_closure.sh` | one flux-closure rung: the same, with the Photochem column and the wind alternating until the elemental flux at the microbar match is its own fixed point |
| `write_reproduce.py` | writes a case's `REPRODUCE.md`; called by both runners, not by hand |
| `status.py` | reads every case back and writes `MODELS.md` sections 7 (results) and 8 (how each model was reached) |
| `archive_results.sh` | moves the results of the tree into `LHS1140b/models_<tag>/`, keeping the inputs where they are, so the catalog can be solved again on a new binary while the old answers stay readable. `ARCHIVE_EXCLUDE` leaves named cases in place |
| `compare_trees.py` | puts two or more of those preserved trees, and the current one, side by side: one row per case with the mass-loss rate, the He I 10830 equivalent width and depth and the base cell, the movement against the first tree, and each ladder's crossing |
| `compare_archive.py` | overlays one solved case on its 2026-08-30 counterpart: T, v, rho, the hydrogen ionized fraction, `n(He 2^3S)` and the He I 10830 profile |

## Re-solving the catalog

A fix that reaches the equations makes every solved case an answer to
equations the code no longer has, and the catalog is then solved again. The
order is the same each time:

```bash
models/archive_results.sh <tag>          # results aside, inputs stay
python3 models/make_models.py            # every case should report "unchanged"
# name the new tree in pick_seed.py's PRESERVED, so tier 0 finds it
models/run_campaign.sh 6 8 'atomic_scalar*'                 # prescribed cases
models/run_rungs.sh 4 8 atomic_photochem_gj1132_kzzprofile  # the nine rungs
python3 models/make_models.py --only atomic_photochem_gj1132x --force
models/run_campaign.sh 6 8 'atomic_photochem_gj1132x*'      # the frozen-profile grid
python3 models/status.py --write
python3 models/compare_trees.py --markdown                  # what moved
```

Every case is seeded by `pick_seed.py`'s tier 0, its own most recent
certified state, so a re-solve costs a few outer passes and not a campaign
from the archive. The frozen-profile grid is regenerated between the rungs
and its own run because it reads the converged column of the fiducial
He/H = 9.7 rung, which the rungs have just rewritten.

## The re-run of 2026-09-15 and what it left out

The tree was solved on 2026-09-14 and solved again on 2026-09-15, because
the base boundary condition took its entropy branch from the cell-1 velocity
and its handover window sat above this planet's base Mach number (plan item
L21; `MODELS.md` section 5 states the measurement). The 2026-09-14 results
are in `LHS1140b/models_20260914_preL21/` and seed the re-run through
`pick_seed.py`'s `tier0`.

The three cases at 0.01 of the fiducial XUV --
`atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13`, the same group's `HeH9.7`, and
`atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` -- were left out of the
re-run: they never solved (item L17 diagnoses the HLLC contact selection on a
base subsonic to Mach 1e-8, which L21 does not touch), so they carry no
2026-09-14 state to continue from and would enter the archive tiers and fail
the same way. Their directories hold their `input.inp` and nothing else.

## What a rerun is expected to reproduce

The OpenMP parallelization of the cell ionization sweep is bitwise identical
to the serial result (`README.md`, feature list) and the runs pin the BLAS
pool to one thread, so `OMP_NUM_THREADS` changes the wall clock and not the
answer; the regression harness nonetheless fixes `OMP_NUM_THREADS=1`, which
is the only setting it compares outputs bitwise under. A different compiler
or a different BLAS reproduces the physics and not the last digits.

The seed does not decide the answer: three seeds four decades apart in the
base velocity land on the same fixed point to 1e-8
(`docs/lhs1140b_stationary_L5c_20260913.md`, T9b), so a rerun whose
`pick_seed.py` chooses a different state -- because another case of the
group has been certified since -- still reaches the same solution. What the
seed does decide is whether the solve reaches it at all and how long it
takes: the archived state of 2026-08-30 seeds
`atomic_scalar_gj1132_kzz1e9/HeH1.60` into a certified solution and the same
state leaves `HeH2.13` uncertified after five outer passes, while the
certified `HeH1.60` seeds `HeH2.13` in one
(`docs/lhs1140b_stationary_L11_20260913.md`).

A seed solved at another composition is carried onto the case's by
`src/utils/map_state_to_grid.py --reservoir He/H <value>`, which the runner
adds by itself: every helium column is multiplied by the ratio of the two
reservoirs, so the base rows arrive at this case's He/H and the shape of a
diffused He/H profile is kept. The option takes any element of the
`# reservoir` line the same way (`--reservoir C/H <value>`), which is what
the elemental-flux closure uses when its photochemical column moves C/H,
N/H and O/H as well as He/H. `load_IC` would otherwise refuse the state,
its `reservoir` metadata field not being the run's.

## The cost

One prescribed atomic case, end to end at 8 threads, takes about a quarter
of an hour: 15 min 15 s for `atomic_scalar_gj1132_kzz1e9/HeH1.60`, MEASURED
2026-09-13, of which 7 outer passes of 117 to 166 s each are the wind. A
closure rung costs the same per iteration, plus about 70 s for the Photochem
column.
