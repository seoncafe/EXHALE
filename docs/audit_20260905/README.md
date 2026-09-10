# Physics and numerical audit diagnostics, 2026-09-05

These files preserve the diagnostic code used for the physics audit and the
development-plan reviews at source revision
`35d9dd5d3ca7eebd01ac658cb12f6021a9b9a36c`.

| File | Purpose |
|---|---|
| `audit_probe.f90` | Original projection-energy, LLF, ROE, H3+, H2 thermodynamics, photoevent-energy, and infrared-equilibrium probes. |
| `run_checks.sh` | Compile the original probe and the existing E1 and oxygen diagnostics in a temporary directory. |
| `roe_equal_pressure.f90` | Call the production ROE estimate for a nonvacuum symmetric rarefaction and check the H3+ high-density boundary and emission normalization. |
| `h2_detector_ratio.f90` | Compare the production H2 channel probabilities with the detector-event interpretation of Chung et al. (1993). |
| `photoevent_energy_check.py` | Preserve the arithmetic checks of revision 1's energy equation and the exponent implied by two Miller table entries. |
| `run_review2_checks.sh` | Compile and execute the additional second-review programs and run the arithmetic checks. |
| `roe_flux_vacuum_probe.f90` | Follow the production ROE speed estimates into the final flux and a local Euler update for symmetric nonvacuum and vacuum-producing expansions. |
| `review3_ledger_checks.py` | Check revision 2's energy identity, published H3+ normalization choices, and a two-cell radiation-exchange example. |
| `run_review3_checks.sh` | Compile and execute the additional third-review flux probe and arithmetic checks. |

From the `EXHALE_v1.00` directory:

```bash
bash docs/audit_20260905/run_checks.sh
bash docs/audit_20260905/run_review2_checks.sh
bash docs/audit_20260905/run_review3_checks.sh
```

Requirements: GNU Fortran and Python 3. `FC` may select a compatible Fortran
compiler; the supplied flags use GNU syntax. The original diagnostic build
also enables OpenMP. Builds are isolated in newly created directories under
`/tmp`; the scripts print their locations and retain them for inspection.
No atmosphere is evolved, and no simulation output or reference product is
regenerated.

The programs print diagnostics. A zero exit code indicates that compilation
and execution completed; it does not mean that every physical check passed.
At the reviewed revision, the equal-pressure ROE probe intentionally exposes
NaN sound speeds. The H2 probe exposes an incorrect detector-ratio inversion
even though the channel sum is one. Assertions for a repaired implementation
should be added when its physical specification is settled.

`photoevent_energy_check.py` evaluates the written proposed equation. Its
chemical-energy value is the value measured by the original production
probe at the reviewed revision. It does not execute a coupled energy solver
or automatically follow future changes to chemical constants.

The published papers are in the workspace's `references/` directory, three
levels above this directory. Findings, measured numbers, interpretation, and
the limits of the tests are recorded in
[the second review](../development_plan_20260905_review2.md) and
[the original physics audit](../physics_numerics_audit_20260905.md).

The additional findings are recorded in
[the third review](../development_plan_20260905_review3.md). Its flux probe
applies a production numerical flux to an isolated Euler update. It does
not run the full EXHALE RK stages, positivity repair, gravity, or chemistry,
and does not implement the proposed ROE changes. The normalization and
radiation examples are arithmetic checks, separate from production-code
execution. Historical probes remain in this directory for reference.
