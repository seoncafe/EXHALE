# Focused checks for PLAN_20260920_review

Date: 2026-09-20.

These checks do not solve an atmosphere or change a case. `inspect_records.py`
reads the selected generations and reports their stored certification values,
identities, and freshly calculated component hashes. The atmospheric residuals
in `records.log` are archived values, not fresh residual evaluations. Its two
small mathematical examples are not EXHALE simulations.

From the EXHALE_v1.00 directory:

```bash
python docs/audit_20260905/plan_20260920_review/inspect_records.py \
  | tee docs/audit_20260905/plan_20260920_review/records.log
```

From this directory:

```bash
gfortran -O0 -fopenmp -fcheck=all -J . \
  ../../../src/modules/lower_atmosphere/h3p_cooling.f90 \
  h3p_domain_counters.f90 -o h3p_domain_counters.x
./h3p_domain_counters.x | tee h3p_domain_counters.log
```

**The defect this probe recorded was repaired on 2026-09-21** (defect 10.2 of
`docs/PLAN_20260920_rev9.md`): the membership tests are strict, a nonfinite
argument has its own category, and the counters were renamed to say that they
count evaluations over a run. The probe below therefore no longer compiles
against the module, and it is kept unchanged as the evidence of what was
measured on 2026-09-20. The assertions that hold today are in
`src/tests/physics_probe/h3p_cooling_limits.f90`.

The Fortran probe compiled and invoked the production cooling module as it
stood on 2026-09-20.
It reproduces the current diagnostic behavior: 30 K increments the below-fit
counter, 300 K increments the outside-non-LTE counter, and repeated evaluations
accumulate overlapping counts. Its assertions characterize the defect; they
are not the desired acceptance assertions for a future repair. Rates and solver
convergence are not validated by this test. All compilation products stay here.

Production H3+ source MD5 at execution:
`a0095be5c0a7d30c46453937a45e1ad1`.

An initial compilation attempt used an incorrect relative source path and failed
before compilation. The corrected command above completed with exit status 0.
No production executable or build products were replaced.
