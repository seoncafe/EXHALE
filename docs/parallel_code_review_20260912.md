# EXHALE v1.00: Parallel Correctness and Optimization Review

Date: 2026-09-12  
Inspected revision: `e8eb6e6ebbbd5a64d6e3c8ba5e9184f10868a264`

## 1. Conclusions and evidence limits

There is a concrete synchronization defect in the constrained-chemistry diagnostic writer. Two additional initialization routines have unsafe concurrent-first-call behavior, although the normal main-program startup protects important uses of them. These should be distinguished from a demonstrated race in the normal physical solution arrays: this audit did not establish such an array race.

Optimization is worthwhile, particularly through reducing repeated parallel-region entry, controlling the independent BLAS thread pool, and matching work granularity to the number of cells. Increasing the global thread count alone is not a sufficient strategy.

The actual scalar PLM reconstruction and its limiter were extracted unchanged and executed at 1, 2, 4, 8, and 16 threads. All measured output differences from the one-thread result were zero on the manufactured profiles. This validates that narrow path, not the complete radiation/chemistry solver.

This report extends [the general code audit](code_bug_optimization_review_20260912.md). Production source was not modified. No full atmospheric integration, full regression suite, race-detector run, or full-solver scaling experiment was performed. Findings explicitly distinguish source-established races, conditional API hazards, measured kernel behavior, and unmeasured optimization proposals.

## 2. Parallel coverage actually present

The Makefile enables OpenMP with `-fopenmp` for GNU and `-qopenmp` for Intel. Searching the implementation found OpenMP rather than an implemented MPI domain-decomposition path in `src`.

| Area | Implementation | Dependency and protection inspected |
| --- | --- | --- |
| Ionization equilibrium | `radiation/ionization_equilibrium.f90:1432,1662` | Dynamic chunks of eight cells; nonlinear scratch and cell context are thread-local; acceptance totals use reductions |
| Radiation coefficients and cooling | `radiation/util_ion_eq.f90:2442,2646` | Fixed blocks of 32 cells distributed dynamically; output slices are disjoint |
| Carrier residual and Jacobian | `lower_atmosphere/diffusive_photochemistry.f90:3991,4222` | Cell-indexed writes; residual maximum reduction; source coefficients rebuilt in thread-local state |
| Reconstruction | `states/PLM_rec.f90`, `states/Reconstruction.f90` | Independent face-state writes; scalar PLM was tested directly |
| Conserved/primitive conversion | `functions/UW_conversions.f90:68,98` | Molecular EOS loops parallelized over cells |
| Marching hydrodynamic RHS | `time_step/RK_rhs.f90:188` | Face loop and cell loop share one parallel region, with an intervening barrier |
| Stationary hydrodynamic rows | `time_step/hydrodynamic_rows_body.inc` | Parallel reconstruction and flux/divergence loops |

The startup comment in `init/init.f90:72` saying that coverage is limited to radiation and ionization is now outdated. The implementation has significantly more parallel regions, so the historical rationale for a default of at most 16 threads needs new whole-solver measurements. The reported historical step rates in that comment were not remeasured here.

Several protections should be retained:

- `parameters.f90:1373` makes `sys_x`, `sys_sol`, `wa`, and `info` thread-local; the cell sweep lazily allocates each thread's nonlinear scratch.
- `ion_cell_state.f90`, `System_HeH_mol.f90`, and `System_HeH_metals.f90` declare the relevant cell state and coefficients `threadprivate`.
- The H/He sweep copies the initialized metal-layout state with `copyin(cx_metal_base)`.
- `carrier_source` at lines 3770–3772 copies `bg_cell(j)` into the thread-local cell context and resets molecular and oxygen coefficients before using them. Reading only the master's coefficients would be incorrect, but that is not what this path does.
- The face-to-cell barrier in `RK_rhs` is necessary. Removing it with `nowait` would allow cells to consume unfinished face fluxes.
- The carrier maximum residual uses a reduction and then a deterministic location scan. Atomic counters in fluxes, composition diagnostics, and the nonlinear solver avoid lost increments.

## 3. P1: Shared diagnostic path and guard are accessed outside synchronization

**Priority: medium; confirmed by source-level conflicting accesses.**

Locations:

- `nonlinear_system_solver/constrained_chemical_equilibrium.f90:454–456`: shared saved `cce_dump_path`, `cce_dump_written`, and `cce_ok_written`.
- Lines 718 and 727: calls for unsuccessful and successful continuation outcomes.
- Lines 1717–1725: `write_continuation_dump` checks the shared flag and writes the shared path before entering `critical(cce_dump_guard)`.
- `radiation/ionization_equilibrium.f90:2248`: the continuation solve is called inside the parallel cell sweep.

The relevant sequence is:

```fortran
if (already) return
call get_environment_variable(envname, cce_dump_path)
if (len_trim(cce_dump_path) .eq. 0) return
!$omp critical (cce_dump_guard)
if (.not. already) then
   already = .true.
   open(... file=trim(cce_dump_path), status='replace', ...)
```

The critical region protects the file write from another critical-region file write. It does not protect it from another thread's earlier environment lookup writing the same path buffer. Successful and unsuccessful cells can supply different environment-variable names, so a cell can use the other outcome's path. The first unsynchronized read of `already` also conflicts with its synchronized write; putting only the writer inside a critical region is insufficient.

Concurrent calls can write the shared buffer even when diagnostic environment variables are unset: the environment lookup precedes the empty-path return. Thus disabling output reduces the observable consequence but does not make the implementation's shared accesses correctly synchronized. The comment claiming no production path reaches this code is inaccurate: the actual continuation routine invokes it.

OpenMP defines unordered conflicting accesses, including read/write accesses, as data races with unspecified results. This finding follows from the source accesses, not from a nondeterministic-output experiment. See the [OpenMP 5.2 specification, memory model](https://www.openmp.org/wp-content/uploads/OpenMP-API-Specification-5-2.pdf).

**Proposed correction:** read both diagnostic paths once during serial initialization and retain separate immutable paths. Check and change each once-only flag under the same synchronization. Alternatively, use a local path buffer and put every access to the shared flag inside the critical region. Mark a dump successful only after the file operation succeeds; the present flag is set before `open` succeeds.

For reproducible diagnostics, record the cell identifier and explicitly define selection: the first thread to enter is not necessarily the lowest-index failing cell. If deterministic selection is required, select the cell after the sweep and then write its saved data serially.

**Validation to add:** simultaneous successful and unsuccessful cells, distinct output paths, unset paths, failed file opens, and repeated threaded execution. Verify file contents and selection as well as absence of crashes. No diagnostic path was enabled in this audit, so no existing dump was overwritten.

## 4. P2: Concurrent-first-call hazards hidden by startup order

### H2 thermochemistry table

`lower_atmosphere/mol_rates.f90:608–612` reads `keq_table_ready` outside a critical region, then checks it again and initializes the table inside the region. The initial unsynchronized read can overlap initialization's write. A thread that observes true outside the region also has no acquisition of that critical region establishing the intended table-publication order.

However, `EXHALE_main.f90:1011` explicitly calls `h2_thermochemistry_init` serially before the parallel work. Therefore this is a **conditional module-entry defect**, not evidence that normal main-program molecular runs currently race on the table.

Prefer an explicit serial initialization requirement, checked at use, for production and standalone drivers alike. If concurrent lazy initialization must remain supported, synchronize all accesses needed to publish readiness and table contents. An inner critical region alone does not make the outer readiness read safe.

### Nonlinear-solver configuration

`nonlinear_system_solver/newton_solver.f90:49–54` initializes shared saved `nt_init` and `nt_force` without synchronization. The normal ionization path has `if(count > 0)` on its parallel sweeps, so the first serial sweep is intended to initialize this state. The comment relies on that calling order.

A standalone caller that first invokes `solve_ieq` concurrently violates the assumption. Make initialization explicit or correctly synchronized rather than requiring every new caller to reconstruct historical startup order. This audit did not establish a normal-main-program failure from this condition.

### Thread-local scratch lifetime

The code allocates scratch when unallocated, but otherwise relies on thread-local storage persisting across regions. OpenMP guarantees persistence only under specified conditions, including stable team sizes/affinity and disabled dynamic team adjustment. The startup sets a requested thread count but does not establish an explicit dynamic-team policy. This is a portability/validation concern, not a reproduced failure. See [OpenMP 5.2 `threadprivate`](https://www.openmp.org/spec-html/5.2/openmpse24.html).

Document supported settings, report actual team sizes, and test `OMP_DYNAMIC` and team changes before promising support. If the code is later turned into a library that changes `N_eq`, existing thread-local array sizes must also be checked; the current fixed-configuration executable does not itself establish such a resizing bug.

## 5. P3: Independent OpenBLAS thread policy is not controlled

The GNU Makefile selects `/opt/miniconda3/lib/libopenblas.so` on this machine. `ldd EXHALE.x` confirmed that the existing executable resolves `libopenblas.so.0` from that directory. The Intel branch already requests sequential MKL.

A direct query of the same OpenBLAS library, in a fresh Python process with the inherited environment, returned:

```text
OpenBLAS 0.3.34 DYNAMIC_ARCH NO_AFFINITY SkylakeX MAX_THREADS=128
openblas_get_num_threads() = 72
openblas_get_parallel() = 1
```

The process CPU affinity contained 72 logical CPUs; `OPENBLAS_NUM_THREADS` and `OMP_NUM_THREADS` were unset. The query reports a configured thread limit, **not a measurement that each LAPACK call actively uses 72 threads**.

The EXHALE OpenMP default cap of 16 does not control this independent non-OpenMP BLAS thread pool. This can waste resources or worsen latency on unsuitable matrix sizes. Do not infer a guaranteed 16-times-72 nesting: the inspected stationary band solves are outside the cell parallel loops. Also, the `dgesvd` calls found in constrained chemistry are in `cce_probe_from_dump`, not evidence that every production chemical cell runs an SVD.

Use a controlled baseline with `OPENBLAS_NUM_THREADS=1`, then benchmark the stationary band factorization separately before deciding whether a larger BLAS team helps. Record the resolved library and both thread policies in setup output. The [OpenBLAS FAQ](https://www.openmathlib.org/OpenBLAS/docs/faq/) and [runtime-variable documentation](https://www.openmathlib.org/OpenBLAS/docs/runtime_variables/) describe independent library threading and its controls.

## 6. Measured reconstruction behavior

Retained test: [run_probe.py](audit_20260905/parallel_review_20260912/run_probe.py). Raw measurements: [results.json](audit_20260905/parallel_review_20260912/results.json).

```bash
python docs/audit_20260905/parallel_review_20260912/run_probe.py
```

The driver extracts `PLM_rec_scalar` and `minmod_mc_slope` unchanged from current source, supplies a positive nonuniform grid and smooth oscillatory scalar profile, and compiles with GNU Fortran 16.2.0, `-O3 -fopenmp`. It requests fixed teams, core placement, and close binding. Every requested team size was observed. Each timing is the best of three batches of 3000 calls; the machine was not reserved exclusively for this experiment.

| Cells | 1 thread, microseconds | 2 threads | 4 threads | 8 threads | 16 threads |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 500 | 3.876 | 3.253 | 2.687 | 2.669 | 2.840 |
| 8000 | 54.110 | 34.671 | 25.304 | 23.961 | 25.113 |

The maximum absolute difference in either reconstructed array from the one-thread reference was **0 in every case**. This is numerical equality on tested data, not an exhaustive proof of parallel correctness.

The best measured times were approximately 1.45 and 2.26 times faster than the corresponding one-thread OpenMP execution. Sixteen threads were slower than eight in both cases. These are kernel results, not whole-solver speedups or a universal recommendation to use eight threads. The one-thread baseline still executes OpenMP directives; a separately compiled serial implementation was not benchmarked.

The actual caller chain is `species_face_flux.f90:259` → `species_face_fraction` → `Reconstruct_scalar` → `PLM_rec_scalar` when PLM is selected. The scalar reconstruction is invoked repeatedly for different species. This makes reducing region-entry overhead a plausible useful change rather than optimization of an unused routine.

## 7. Optimization priorities

### A. Combine repeated species reconstruction work

Evaluate reconstructing all active species in one parallel region rather than starting a region for each scalar field. Preserve the existing limiter and ghost/face conventions. Cell blocks spanning several species may improve reuse, but Fortran's array layout and the number of active species should decide loop order. Do not merely add an outer parallel species loop while retaining inner parallel regions: nested regions can add overhead, and shared work arrays must be reassessed.

The scalar routine initializes both output arrays in a serial loop before its parallel stencil loop. That serial initialization and memory traffic also limit scaling, especially on larger grids. Test a fused initialization/stencil design with explicitly owned boundary values; do not remove initialization of faces outside the complete stencil.

### B. Tune granularity separately for rate tables and nonlinear chemistry

`util_ion_eq.f90:241` fixes `xuv_rate_block=32`. At the default 500 physical cells and four ghosts, there are exactly 16 work blocks. This limits available block parallelism in these loops regardless of a larger requested team. Dynamic scheduling may be useful for uneven work, but uniform coefficient blocks should be compared with static distribution.

The chemistry sweep uses `schedule(dynamic,8)`, appropriate as a candidate for variable nonlinear iteration counts, but its optimal chunk length is not established by the source. Benchmark several chunk sizes with the physical state and convergence tolerances held fixed.

Do not conflate these scheduling choices with `xuv_field_block_cells`, currently zero in `ionization_equilibrium.f90:86`. Radiation columns are accumulated in a starward sequence and the updated composition is carried between field blocks. Changing field-block size changes the numerical update structure, not merely thread scheduling. Preserve that dependency or validate a deliberately different radiation/chemistry iteration.

### C. Keep dependent operations serial unless their algorithm changes

The carrier block-tridiagonal forward/back substitutions carry spatial dependencies. A plain OpenMP loop over their rows is incorrect. Parallel cyclic reduction would be an algorithmic change requiring conditioning, residual, and performance tests; it is not the first optimization to pursue on a 500-cell problem.

Likewise, do not parallelize the stationary solver's finite-difference residual evaluations by just distributing columns: those evaluations change module state and caches. Thread-local chemical scratch alone does not make the complete residual evaluator reentrant.

### D. Reduce diagnostic contention and configuration lookups

Keep correctness counters, but consider collecting local counts and reducing them once instead of performing an atomic update for every small cell solve. Preserve exact integer totals and rollback semantics. P1's serial diagnostic initialization also avoids repeated environment reads: with diagnostics unset, the current early return never sets `already`, so lookup is repeated on subsequent continuation calls, contrary to the nearby once-only cost comment.

Use `default(none)` where practical to expose new shared-variable decisions at compilation. It does not detect hidden module writes in callees, so maintain a short contract for mutable state below each parallel entry point.

### E. Separate elapsed time from aggregate thread time

`eval_cool` reduces `ec_t` from block timers and explicitly documents several channels as accumulated thread time. Those totals must not be presented as elapsed-time shares of the full run. Measure wall time outside each parallel region and report worker time separately. Include actual team size, CPU affinity, BLAS policy, and accepted/rejected work counts in scaling reports.

## 8. Follow-up validation and deliverables

Recommended order:

1. Correct P1, then stress both diagnostic paths concurrently without touching existing diagnostic products.
2. Establish initialization and thread-setting contracts for P2 and test standalone entry paths.
3. Compare bounded atomic, molecular, and oxygen/carrier fixtures at 1/2/4/8/16 threads. Compare physical arrays, conservation, reaction residuals, acceptance classes, and termination status, excluding nondeterministic log order and elapsed times from numerical comparisons.
4. Measure complete accepted-step and stationary-pass wall times before changing defaults. Test a large grid separately from the default grid.
5. Apply region fusion or scheduling changes only after profiling, then repeat the affected physics checks. Matching an incorrect reference is not a physical acceptance criterion.

New artifacts are this report and `docs/audit_20260905/parallel_review_20260912/`, which retains the diagnostic script, generated Fortran source, compiler products, and raw results. Earlier audit files were preserved. No production fix, data regeneration, or Git staging/history operation was performed.
