# Recent stage-2 code review diagnostics

Date: September 12, 2026

See [the review](../../Update_EXHALE_stage2_review_20260912.md) for the findings, source locations, and limits. Production source was not changed.

## Retained evidence

- `build.log`: fresh private production build.
- `carrier_retry.log`: existing suite, 127 PASS and 0 FAIL.
- `element_operator.log`: existing suite, 24 PASS and 0 FAIL.
- `carrier_state_contract_probe.f90`: diagnostic derived from the current carrier test driver; adds explicit density, EOS, temperature, and mass-ceiling measurements. Original full-line comments were omitted; fixture statements and assertions were retained.
- `carrier_state_contract_final.log`: final synthetic-fixture and cell-verdict measurements. The earlier `carrier_state_contract_probe.log` predates the added entry-mass diagnostic and is retained separately.
- `stationary_retry_history_probe.f90`: independent invocation using the same initialization routines, with mass normalization explicitly corrected before injecting rejected trials.
- `stationary_retry_history_final.log`: final mass-closed history test. The earlier `stationary_retry_history_probe.log` used the original synthetic normalization and is not the result used for the final history finding.
- `molecular_wind_state_probe.f90`: the archived September 11 production-module driver adapted to the current carrier API, with EOS handoff checks.
- `current_probe_run.oY88nt/`: the completed 500-cell molecular probe, copied inputs and initial conditions, state and face-flux records, provenance, and logs.
- `reviewed_source.diff`, `revision.log`, and `source_sha256.log`: reviewed working-tree provenance. The revision alone is insufficient because the source already had uncommitted changes.

`current_probe_run.oY88nt/carrier_state_contract_probe.log` predates the added synthetic entry-mass diagnostic. Its provenance records the executable and source hashes at that execution. Final source versions and later measurements are retained separately rather than rewriting that run's provenance.

## Reproduction

From the EXHALE project root, create a fresh private build:

```bash
audit_build=$(mktemp -d /tmp/exhale_stage2_review_20260912.XXXXXX)
make -j4 OBJDIR="$audit_build/build" EXE="$audit_build/EXHALE.x"
bash docs/audit_20260905/stage2_review_20260912/run_probes.sh "$audit_build"
```

The script compiles all three diagnostic drivers and creates a new unique run directory. It requires the current source API and the compiler/OpenBLAS installation used on this machine. It does not modify the original fixtures. The actual retained build is `/tmp/exhale_stage2_review_20260912.XJbhQc`; rebuild if temporary files are removed.

The existing suites can be repeated independently with private output directories:

```bash
EXHALE_TEST_OBJDIR="$audit_build/carrier_retry" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 bash src/tests/carrier_retry/run.sh
EXHALE_OBJDIR="$audit_build/build" EXHALE_TEST_OUT="$audit_build/element_operator" OPENBLAS_NUM_THREADS=1 bash src/tests/element_operator/run.sh
```

The carrier suite compiles its own dependencies with checking flags. The element suite and custom diagnostics link the private optimized production objects. Driver checking flags do not imply that those production objects were compiled with runtime checks.

## Interpretation

The molecular diagnostic performs one hydro-only production solve and one carrier update, not the complete production outer controller. Its copied input preserves the original coupled key; the diagnostic explicitly changes the Newton registry. The 180 s timeout did not fire in the recorded run, which completed in approximately 9.04 s. Exit 0 means the measurement completed, not that the full atmosphere was certified.

In `steps.tsv`, the hydro count is the requested iteration cap; the actual count belongs to the solver log. The species `info` value is an inherited diagnostic placeholder, not a success return from the current carrier API. The reported five species steps and drift are the routine's actual outputs.

The history test injects rejection through an existing test hook; it is an interface test, not a prediction of rejection frequency. The mass-ceiling counterexample is a call to a cell verdict, not a complete certified atmosphere. The unchanged synthetic carrier fixture has an invalid mass normalization, explicitly measured in the final diagnostic; its results must not be presented as a physically mass-closed atmospheric calculation.

No long convergence runs or full regression suite were performed. No existing simulation products were regenerated.
