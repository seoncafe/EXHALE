# Plan of 2026-09-13: the findings of the two reviews of 2026-09-12

User instruction (2026-09-13): read `docs/code_bug_optimization_review_20260912.md`
and `docs/parallel_code_review_20260912.md`, confirm each finding in the code,
write the correction plan, and carry it out; run the binary in parallel
wherever that shortens the verification.

Every finding was re-read against the source at `e8eb6e6` before this plan
was written. Verdicts:

| finding | verdict | where confirmed |
| --- | --- | --- |
| B1 clipped trapezoid uses the endpoint fluxes | CONFIRMED, wrong integral | `sed_read.f90` `sed_band_integrated_flux`: `0.5*(f_prev+f)*(wb-wa)` after clipping to `[wa,wb]` |
| B2 one covering segment returns zero | CONFIRMED | same routine, `if (nin .lt. 2) F_band = 0` counts segments, not coverage |
| B3 five exhausted cycles return success | CONFIRMED | `diffusive_photochemistry.f90` `equilibrate_chemistry_at_fixed_conserved_state`: `ok = .true.` at entry, the loop falls through without clearing it |
| B4 off-simplex ledger not consumed | CONFIRMED | same routine checks `n_nonfinite` and finite `f_sp` only; `ledger%n_offsimplex` (set at `ionization_equilibrium.f90:2732`) is never read |
| B5 mapper accepts inconsistent or out-of-support grids | CONFIRMED | `src/utils/map_state_to_grid.py`: hydro radii only, `numpy.interp` clamps outside the support, no monotonicity or column checks |
| B6 mapper copies the `# coupling` claims | CONFIRMED | same tool copies every header line but `# rows`; `load_IC` reads `mode`, `t_phys`, `certified` from it |
| P1 dump path and once-only flag read outside the critical region | CONFIRMED | `constrained_chemical_equilibrium.f90` `write_continuation_dump`: `if (already) return` and `get_environment_variable(envname, cce_dump_path)` precede `!$omp critical`; the flag is set before `open` succeeds |
| P2 lazy first-call initializations unsynchronized | CONFIRMED as conditional | `mol_rates.f90` `keq_H_H_to_H2` reads `keq_table_ready` outside the critical; `newton_solver.f90` sets `nt_init`/`nt_force` with no synchronization; both are protected only by the main program's serial first calls |
| P3 OpenBLAS thread pool uncontrolled | CONFIRMED | GNU link uses `libopenblas.so`; nothing in `src/` or the Makefile sets its thread count; `ldd EXHALE.x` resolves it |
| O1, O2, A, B, C, E | proposals that need measurement first | no code change in this plan beyond the counters O2 asks for; see S8 |
| O3 parse the spectrum once | ADOPTED | folded into S1 |
| O4 keep the ghost-count refresh | agreed, no change | `steady_newton.f90` |

Rules: `docs/worker_rules.md`. Items S1 to S9. Scratch `<scratchpad>/S/`. Runs
of the binary use `OMP_NUM_THREADS=8` and are launched concurrently where they
are independent; the regression harness already runs its cases side by side.

## S1: the band integral of a loaded spectrum (B1, B2, O3)

`sed_read.f90`: one subroutine `sed_band_fluxes(n_band, w_lo, w_hi, F, covered,
status)` reads the file once, validates it (finite, nonnegative flux; positive,
strictly increasing wavelength; a malformed row stops the run with the row
named, as `read_sed` does), and for each band integrates the piecewise-linear
spectrum exactly: the flux is interpolated at each clipped edge before the
trapezoid, so a constant and a linear spectrum integrate exactly whatever the
tabulated nodes; a single segment covering the band is a valid integral.
`covered(b)` is true only when the file reaches both edges of band b; a band
not covered gets `F = 0` and `covered = .false.`, and `input_read` then STOPS
with a message naming the band and the key that would state it (a positive
field must never be silently turned into zero; the run's own convention for
the triplet band is the same). `input_read` calls it once for LW, B3, B4.
Tests: a new program in `src/tests/spectrum_type/` (`fuv_band_quadrature.f90`)
writing small spectra to scratch files: constant and linear spectra with both
edges clipped, one edge clipped, edges on nodes, unequal segment widths, one
covering segment, partial coverage, a malformed row; exact answers to 1e-12.
RED on the entry text for the clipped and the one-segment cases.

## S2: the thermochemical closure contract (B3, B4, O2 counters)

`equilibrate_chemistry_at_fixed_conserved_state`: `ok` starts false and becomes
true only when a cycle's temperature increment is below the tolerance AND the
sweep ledger of that cycle reports no nonfinite cell and no off-simplex cell
AND `p`, `T`, `heat`, `cool`, `eta` are finite with `p > 0`, `T > 0` on the
physical cells. Exhaustion returns false with the final increment; the caller
(`relax_photochemical_composition`) already restores the trial on `.not.
chem_ok` and names the ending `carrier_relax_chemistry_refused`. A named reason
(`converged`, `exhausted`, `nonfinite`, `off_simplex`, `state_not_admissible`)
and the last increment are returned and counted; the relaxation report prints
the cycle counts, closures, exhaustions and refused trials (O2's counters).
Tests in `carrier_retry`: a test knob `chemistry_cycles_cap_for_test` forces
exhaustion (refused ending, entry composition and thermal state restored bit
for bit), and the existing convergence rows stand.

## S3: the state mapper (B5, B6)

`src/utils/map_state_to_grid.py`: validate finite, positive, strictly
increasing radii in both source files and the target; the two source files'
radii must agree to 1e-12 relative; the column labels must carry the species
the mapping needs and `n_H > 0` everywhere; every PHYSICAL target center must
lie inside the source's physical support, else refuse (exit 2, message); the
ghost rows are extrapolated linearly in ln r from the nearest two source
physical cells, and the tool says so. Metadata: the output is an
initialization seed: the `# coupling` line is rewritten with `mode=init
t_phys=0 certified=F` (the other tokens kept), and the source line is
retained as `# mapped-from-coupling:`; `restart_schema`/reservoir lines are
kept (they describe the composition, which the mapping preserves). A test
script `src/tests/state_mapper/run.sh` exercises identity mapping, a
2e-4 shifted grid (round trip), mismatched species radii, out-of-support
targets, and the metadata rewrite. The Model A fixture IC is regenerated with
the corrected tool and the 3-pass run repeated.

## S4: the continuation dump (P1)

`constrained_chemical_equilibrium.f90`: the two paths are read once, serially,
by `cce_dump_paths_read` called from `input_read` (after the keys) into two
saved strings; `write_continuation_dump` takes the path by argument, does
every read and write of the once-only flag inside the critical region, sets
the flag only after `open` succeeds, and writes the cell index into the
dump (the selection is the first thread to enter, stated in the file).

## S5: explicit initialization contracts (P2)

`mol_rates.f90`: the lazy initialization is removed; `keq_H_H_to_H2` stops
with a message if the table is not ready (the main program initializes it
serially; test drivers that use it call the initializer). `newton_solver.f90`:
`newton_solver_read_environment` is called serially from `input_read` and the
lazy block is removed. Team policy: `init.f90` calls `omp_set_dynamic(.false.)`
and reports the team size actually obtained; the outdated coverage comment
there is replaced by the present inventory (the parallel review's table).

## S6: the BLAS thread policy (P3)

`init.f90`: at startup, unless `OPENBLAS_NUM_THREADS` is set, the OpenBLAS
thread count is set to 1 through `dlsym` (`openblas_set_num_threads`, resolved
at run time so that an MKL or reference-LAPACK build is unaffected), and the
setup report states the LAPACK library and both thread policies. MEASURED
before adopting the default: wall time of the `wasp_full_newton` reload
(the band factorizations `dgbtrf`) with OpenBLAS at its default 72 threads
and at 1.

## S7: verification, in parallel

Suites: `spectrum_type`, `carrier_retry`, `carrier_returned_state_acceptance`,
`steady_species_rows`, `krylov_and_dogleg`, `certification`, `state_mapper`,
`attempted_step`; the Model A fixture (3 and 12 passes); the oxygen case
with the loaded spectrum (B3/B4 derived fluxes re-measured against the hand
trapezoid, now exact); `wasp_full_newton`; `make check`. Independent runs are
launched together.

## S8: measured-first items, not changed here

O1 (checkpoint copies), O2 beyond the counters, A (fused species
reconstruction), B (block sizes), C (dependent loops), E (thread time vs
elapsed): each needs a profile on a representative molecular case before a
design is chosen; the counters of S2 and the team/BLAS report of S5/S6 are
the instrumentation those profiles need. Recorded in `docs/TO_BE_DONE.md`
with the measurement each requires.

## S9: records

`docs/Update_EXHALE_stage2.md` section 10, this plan's status table,
`docs/ISSUES_20260909.md` section 2, the two review documents left as they
are (they are records), README/manual where a key's behavior changed (the
covered-band stop).

## Status at the close of 2026-09-13

| item | status | record |
| --- | --- | --- |
| S1 | DONE | `sed_band_fluxes`, `band_from_spectrum`, `fuv_band_quadrature` 14/0, `spectrum_type` 170/0; B3 4.686E+02 / B4 4.443E+05 against the exact 4.685913E+02 / 4.443075E+05 |
| S2 | DONE | closure contract with named reasons and counters; `carrier_retry` 142/0 (four new rows, one restated); Model A fixture 3 passes unchanged |
| S3 | DONE | mapper validation and seed metadata; `state_mapper` 17/0; fixture `IC/` regenerated |
| S4 | DONE | paths read once serially; flag inside the critical region; cell recorded; dump format 2 |
| S5 | DONE | explicit initializations at every serial entry point; team policy and report |
| S6 | DONE | `blas_thread_policy.f90`; 16.33 s vs 16.32 s on `wasp_full_newton` |
| S7 | DONE | twelve suites concurrently, all PASS; `make check` PASS, every case byte-identical |
| S8 | RECORDED | `docs/TO_BE_DONE.md`, measured-first items |
| S9 | DONE | update log section 10, this table, ISSUES section 2, README/HOWTO/manual/physics document/input schema |
