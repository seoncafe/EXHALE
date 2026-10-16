# EXHALE v1.10

Branch `v1.10` of `git@github.com:seoncafe/EXHALE.git`. The released code is
branch `v1.00`; this branch is the development line that followed it. It was
copied from the v1.00 working tree on 2026-10-12 (the copy is recorded in
`provenance/source_copy.json`), so its history starts here and does not
continue the v1.00 commits.

## Build and run

```bash
make                 # gfortran on PATH (conda-forge) + its own LAPACK; incremental
make FC=ifort        # or make ifort / make ifx
./run_EXHALE.sh      # Tk GUI: writes input.inp, builds, runs
cd <run dir> && /path/to/EXHALE.x   # direct: reads ./input.inp, writes ./output/
```

A run is a directory: `input.inp` and the optional `metals.inp`,
`opacity.inp`, `base.inp` beside it, `output/` written into it. The input
keys, outputs and physics are those of v1.00 (its `README.md`,
`README_HOWTO.md` and `docs/EXHALE_user_manual.pdf` on branch `v1.00`), plus
the additions below. The inputs of the two atomic marching cases of the
first milestone are `examples/regional_convergence/Huang.inp` and
`Salz.inp`; `examples/exhale_io.py` reads every output file.

## What differs from v1.00

- `src/modules/time_step/regional_convergence.f90`: the regional
  stabilization criterion of the time marching (finite-duration stability of
  the mean, range, inventory and mass flow of a radial region), configured by
  the `Regional ...` input keys (`Regional convergence`, `Regional lower
  radius`, `Regional upper radius`, `Regional window seconds`, ... ; unknown
  `Regional` keys stop the run). It is distinct from the stationary
  acceptance gates, which still require the stationary root test.
- `src/modules/states/stationary_boundary_chemistry.f90` and
  `reservoir_composition.f90`: the lower ghost chemistry and thermodynamics
  eliminated at a fixed column, shared by the Newton residual and the
  certification of the stored pair; an unconverged boundary blocks
  certification.
- The regression harness `backup/regression/` of v1.00 is not in this tree, so
  `make check` fails here by design until it is carried over.

## What is not in the repository

As on `v1.00`: `md/` (the development record), `src/tests/`, `backup/`,
build products, executables, run output. In addition `runs/` (every
calculation of this tree), `provenance/` (historical inputs, states, source
copies, executables) and `workflows/` (plans, execution records and their
builds) stay on disk; the links below into `md/` and `workflows/` resolve
only in the working tree.

## Development record (working tree only)

Development copy created from EXHALE v1.00 on October 12, 2026 to implement
the reviewed stationary-solver plan. The original directory is preserved.

The current first milestone is practical time-marching convergence for the
Huang and Salz atomic cases. Simultaneous stationary roots remain the ultimate
objective. See the [active execution record](md/PROGRESS_20261015_marching_first_milestone.md)
for current calculations, measured histories, known discrepancies, and the
preserved conditional p18 result. Neither atomic case has yet been accepted.

Build with `make -j4`. Source files are in `src/`, required model data are in
`inputdata/` and `cooling_data/`. The copied source and data are recorded in
`provenance/source_copy.json`; the initial 2,574 files were compared byte for
byte before changes. Existing simulation products and build objects were not
copied into the active source tree.

The [plan](md/PLAN_20261012_stationary_failures_F1_F3.md) and
[execution record](md/IMPLEMENTATION_20261012_stationary_failures.md)
describe the earlier strict-root work. Historical inputs, states, source and executables
are retained under `provenance/history/`, with checksums in
`provenance/history_archive.json`. New calculations write separate directories
under `runs/`. Workflow plans and durable execution records are in `workflows/`.

The shared stationary lower-boundary chemistry and thermodynamics operation
is implemented. Newton residuals and final/stored-pair certification use the
same closure, and an unconverged boundary blocks certification. The
[approved design](md/DESIGN_20261012_stationary_boundary_closure.md) explains
its physical conditions; the [measured report](md/stationary_boundary_thermochemistry_efficient_analysis.md)
contains seed comparisons, full profiles, and residual replay results.
The [integration report](md/IMPLEMENTATION_20261012_stationary_boundary_chemistry.md)
describes returned-state and chemical-rate restoration corrections.
These corrections do not establish stationary solutions for the historical
Yan, Huang, or Salz states. Their complete certificates still refuse them.

The finite validation sequence is recorded in
`workflows/boundary_consistent_validation/progress.md` and
`workflows/boundary_ptc_review/progress.md`. Diagnostic completion,
a complete refused Newton solve, and stationary certification are separate
outcomes. Earlier failed and timed-out calculations remain under `runs/`.

The preserved production executable has SHA-256
`66be2bb9f5dbda1b5e07f3168ee90c6d7fd6f7be3485774364eb39b4011d1058`.
The current source additionally corrects the whole-column species transport
summary and makes `make check` fail when its regression harness is absent.
These changes were built separately in
`workflows/marching_first_milestone_20261015/`; the build record and focused
runtime check are documented in the active execution record. A successful
diagnostic run is not evidence of general regression or physical acceptance.

Last updated: 2026-10-16 23:25 KST
