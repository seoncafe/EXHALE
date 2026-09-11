# The named Newton cases under the section 145 gates

Measured 2026-09-03 on tree binary `62cc6c64` (the final section 145 tree), each
case copied out of `backup/regression/` into scratch and run with
`OMP_NUM_THREADS=2`. The acceptance is now the cellwise-max residual norm and the
Riemann face mass-flux spread over `r >= r_flux`, tolerance `2.0e-5`.

| case | `info` | JFNK iterations | log10 Mdot | face spread at accept | gate |
|---|---|---|---|---|---|
| `newton_rsw01`        | 0 | 89  | 13.21 | 8.890e-12 | met |
| `newton_rsw05`        | 0 | 84  | 13.21 | 1.515e-11 | met |
| `wasp_full_newton`    | 0 | 59  | 13.21 | 1.352e-11 | met |
| `wasp_he23off_newton` | 0 | 38  | 13.22 | 1.635e-11 | met |
| `arm_heh1_x2matched`  | 0 | 120 | 10.31 | 4.602e-13 | met |
| `armD_D2_LW_newton`   | 2 | 1   |  8.88 | 2.116e+01 | NOT met |
| `jfnk_hd189`          | - | -   |  -    | -         | never runs JFNK |

**Every case that converged before still converges, at the same iteration count
and the same mass-loss rate.** The new gate does not refuse any of them: the
accepted face spreads are `4.6e-13` to `1.6e-11`, seven to nine orders below the
`2.0e-5` threshold, which is what the threshold was re-derived from.

**A prediction of mine was wrong and is corrected here.** I expected
`arm_heh1_x2matched` to be refused by the face-flux gate. It is not: it is the
tightest of the seven, at `4.602e-13`. What made the cell-centred functional look
bad on that case is the launch region -- its `r >= 1.03` spread is `4.561e-09`,
four orders above its gate window -- and the gate window has always started at
`1.2 R_p`.

**`armD_D2_LW_newton` fails for a reason that predates section 145.** JFNK
reports "no descent along the Newton or the damped Gauss-Newton direction" at its
first iteration, `||grad merit|| = 3.8e+11`, and the marching retry does not
recover it inside its 2100-step cap. Its face spread is `21`, so no threshold
would have accepted it; the old cell-centred gate refused it too. This is the
molecular Lyman-Werner case, item (S)/(P) territory, not an acceptance
question.

**`jfnk_hd189` is misnamed.** Its `input.inp` carries no `Solver:` key at all, so
despite the name it never runs JFNK; it marches with residual-based convergence
(`Resid tol: 1.0e-5`). Stopped by hand at step 76183 with `du = 1.38e-2` and no
sign of a stop -- the same HD 189733 b base behaviour item (AD) and (AA) track.
Either the case gets `Solver: Newton` or it gets a name that says it is a
marching case; not changed here.
