# `carrier_model_a_newton` -- the carrier reload on the matched Model A configuration

## Purpose

The stationary-solver fixture re-pinned on the configuration that matches
Koskinen et al. (2022) Model A (`docs/koskinen2022_model_a_comparison.tex`,
section 5: the earlier fixture `carrier_elem_newton` is a different
configuration, Roche domain, power-law spectrum, `Rate/2 + Mdot/2`, He 2^3S
on, no ionization transport). Here the stationary route and the marching
solution describe ONE configuration, so what the partitioned loop does to
the H2 extent can be read against Model A directly.

## Contents

- `input.inp`: `benchmarks/koskinen2022_model_a/matched_hnu_minus_I/input.inp`
  with the stationary keys appended (`Solver: Newton 100.0`, `Coupled carrier
  solve: False`, `Resid tol: 1.0e-8`, `Restart intent: stationary`) and the
  stage-2 `du_th` set to 10 so that no march runs before the solve.
  `base.inp`: the matched run's.
- `IC/`: the matched run's final state (`matched_hnu_minus_I/output`, mapped onto the current grid by `src/utils/map_state_to_grid.py` as an initialization seed: `mode=init t_phys=0 certified=F`; the original state is `IC_grid_20260905/`;
  2026-09-05), renamed for `Load IC? True`.
- `run.sh`: as for `carrier_elem_newton`.

## Recipe

```
cd backup/regression/carrier_model_a_newton
EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=40 OMP_NUM_THREADS=8 ./run.sh /abs/path/to/EXHALE.x
```
Pinned 2026-09-12. Not in `DEFAULT_CASES`; a diagnostic entry point, not a golden.
