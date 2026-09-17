# L15 stage 2: the record of the thread-reproducibility measurement, 2026-09-16

Item L15 of `docs/PLAN_20260913_lhs_stationary.md`, section 2 of
`docs/PLAN_20260916_rev3.md`. The account is
`docs/lhs1140b_stationary_L15_20260916.md` section "Stage 2"; this
directory holds the identities and the raw outputs the account stands on.

## Identities

| what | md5 |
|---|---|
| `EXHALE_L15.x`, the private build carrying the trace, as the measurements were taken | `f23da8d84fafc2f0587b610add56e714` (source `e378cffce5e64f692401106a8a8ce52f`) |
| `EXHALE_L15.x` rebuilt from the delivered source, which differs from the above in three lines of one comment | `1e8afdcdf6b2a84b723562ecdce4b5a6` (source `33cbd5065c7eea7519918be91d5dde5f`); it reproduces the repeat run's trace and its checkpoint exactly (`trace_rep8_rebuilt.txt`, dump md5 `3830aab11ec1018e740db6c5ab4fb11f`) |
| `EXHALE_L15ctl.x`, the same tree with the entry text of `steady_newton.f90` | `7fece24a7561d6413fb30d07ba9e90c8` |
| `EXHALE.x` of the tree, for reference only (built 2026-09-16 10:43, before the day's other items edited `src/`) | `c2e9c9990b9f14f1be8cd77abca68945` |
| `backup/regression/atomic_elem_newton/IC/Hydro_ioniz_IC.txt` | `097ffbb982bfe8f92910910864c557de` |
| `backup/regression/atomic_elem_newton/IC/Ion_species_IC.txt` | `71a72cfcfe7ede661cb595abc25c06ee` |
| `input_reload.inp` (the fixture's `input.inp` with the three lines of the reload recipe) | `4708a683d9ba07e73a24f49548599887` |

Both binaries were built with `make OBJDIR=build_L15 EXE=EXHALE_L15.x` and
`make OBJDIR=build_L15ctl EXE=EXHALE_L15ctl.x` from the tree of 2026-09-16
evening, which carries the uncommitted edits of the other items running that
day in about twenty source files. The two binaries differ in the text of one
file and in nothing else, so the pair is the control the comparison needs;
neither is comparable with `EXHALE.x`, which is a build of an earlier state of
those same files.

## The fixture and the settings

`backup/regression/atomic_elem_newton` reloaded from `IC/` on a scratch copy,
per that case's README: `Load IC? True`, `Coupled carrier solve: False`,
`Restart intent: stationary`, with `EXHALE_OUTER_PASSES=1
EXHALE_DIFF_OMEGA=0.5` and the JFNK cap stated per run. `Restart intent:
stationary` fixes the pseudo-time start at 1.0, so no `EXHALE_PTC_DTAU0` was
set.

| run | threads | binary | JFNK cap | other |
|---|---|---|---|---|
| `t8_a`, `t8_b`, `t8_c` | 8 | `EXHALE_L15.x` | 200 | trace on in `t8_a`, `t8_b` |
| `t8_ctl` | 8 | `EXHALE_L15ctl.x` | 200 | |
| `t1_a` | 1 | `EXHALE_L15.x` | 200 | trace on |
| `t16_a` | 16 | `EXHALE_L15.x` | 200 | trace on |
| `rep8`, `rep1` | 8, 1 | `EXHALE_L15.x` | 30 | trace, `EXHALE_L15_REPEAT=21`, `EXHALE_L15_DUMP_AT=21` |
| `g8_a`, `g8_b` | 8 | `EXHALE_L15.x`, `EXHALE_L15ctl.x` | 80 | `EXHALE_PTC_RAMP_GUARD=0`, the configuration of the L4h memo section 5.0 |

## Files

- `trace_t8_a.txt`, `trace_t8_b.txt`, `trace_t1_a.txt`, `trace_t16_a.txt`:
  the record of the whole solve, 105 outer-iteration records, 105 merit
  records and 713 records of the linear cycle (the right-hand side of each
  cycle, and for each of the 152 products the preconditioner input, the
  preconditioner output, the action of the operator and the Arnoldi column).
  Each record carries FNV-1a 64 hashes of the bit patterns.
- `trace_rep8.txt`, `trace_rep1.txt`: the same with the repeated evaluation
  and the checkpoint at outer iteration 21.
- `checkpoint_iter21_t8.bin`, `checkpoint_iter21_t1.bin`: the state the next
  evaluation reads at outer iteration 21 (the unknowns, the composition, the
  residual, the model base point, the column and row scales, the faces of the
  species box, the active bound set, the pseudo-time), stream access. What is
  NOT in them: the module caches of the residual assembly outside
  `steady_newton` (`ionization_equilibrium`, `util_ion_eq`,
  `constrained_chemical_equilibrium`, `diffusive_photochemistry`), which have
  no accessor this module may read. The comparison that does not need them is
  the repeated evaluation inside one process, which is what `trace_rep*`
  carries.
- `run_*.log`: the run logs, including the two `wasp_he23off` runs of the
  gate-off byte-identity check.
- `comparisons.txt`: the raw output of every comparison, generated from the
  files in this directory.
- `input_reload.inp`, `metals.inp`: the configuration.

## What the comparisons say

The four trace files of 1, 8 and 16 threads have one md5,
`bd12e7c8ad917beda538d789080e10e3`. The run logs of three 8-thread runs, of
the control binary at 8 threads, and of the 16-thread run differ from each
other only in the wall-clock seconds of the pass line; the 1-thread log
differs further in the two lines that report the thread count and in the
exit-time floating-point exception summary. The checkpoint dumps of the
8-thread and the 1-thread run are byte-identical. The repeated residual and
the repeated action of the Jacobian inside one process are bitwise equal at
both thread counts, and equal to each other across them.

The gate-off check on `backup/regression/wasp_he23off` (cold start, one
thread, both binaries, scratch copies): `Ion_species.txt` byte-identical over
the whole file, `Hydro_ioniz.txt` identical in every data line (md5
`c48c86234c18e9a242803ad194a20e10`) and differing only in the `# provenance`
line, which carries the wall-clock time of the run.

## Step 5, the molecular configurations

Added the same day, with a rebuilt private binary `EXHALE_L15.x` md5
`21ba370bf1ae1dd0816c299fc2b7abc0` (the tree's other files had moved between
the two builds, so every comparison of this step is inside that one binary).
Both cases reloaded from their own `output/*_IC.txt` with `Load IC? True` and
`Restart intent: stationary`, keeping each case's own keys
(`Secondary_ionization: Immediate`, `Molecular carrier transport: True`),
`EXHALE_PTC_DTAU0=1.0 EXHALE_OUTER_PASSES=3`. The scratch copies were placed
three directories below a copy of `LHS1140b/sed/` so that the `Spectrum file:
../../../sed/...` of the input resolves.

| case | runs | trace md5 | states, data lines |
|---|---|---|---|
| `molecular_scalar_gj1132_kzz1e9/HeH2.13` | 8, 8, 1 threads | `d5768c1c3c906a65a60fe0e62043e125` | `Hydro_ioniz.txt` `4b7ce7f1d9b86f34f3f6ea189943cac7`, `Ion_species.txt` `0c29e3e53d11b4d4a5f7c8b6297a7de7` |
| `molecular_scalar_gj1132_wellmixed/HeH0.083` (the stalled state) | 8, 8 threads | `e465602dae6f7cfdae3dac4320891a55` | `Hydro_ioniz.txt` `806f3e3a7fba26fb89bf6fe6c9173538`, `Ion_species.txt` `6fa32cde8fb3c735edb72224c6716379` |

Files: `trace_mol_kzz1e9_t8_a.txt`, `trace_mol_kzz1e9_t8_b.txt`,
`trace_mol_kzz1e9_t1_a.txt`, `trace_mol_wellmixed_t8_a.txt`,
`trace_mol_wellmixed_t8_b.txt`, the five `run_mol_*.log`, and
`input_mol_*.inp` with `base_mol_*.inp`.

