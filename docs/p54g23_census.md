# Named-case census: A (tree) vs G23 (all three row scales)

Generated 2026-09-03 16:49 KST. A dash means that build's run
had not finished when the table was made. 'golden' says whether the
numeric output moved; 'du(r>=r_esc)' is the marching stop metric
recomputed on the written state, beside the gate's sp(1.20).

| case | stop A/G23 | info A/G23 | its A/G23 | evals A/G23 | r_esc | du(r>=r_esc) A | du(r>=r_esc) G23 | sp(1.20) A | sp(1.20) G23 | 1.03 A | 1.03 G23 | 1.10 G23 | 1.01 G23 | logMdot A | logMdot G23 | d(logMdot) | 1.01 gate G23 | golden |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| armA_LW | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 3.498e-02 | 3.498e-02 | 7.226e-02 | 7.226e-02 | 4.346e-01 | 4.346e-01 | 1.377e-01 | 4.527e-01 | 10.52 | 10.52 | +0.00 | FAIL | same |
| armA_noLW | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 9.230e-03 | 9.230e-03 | 2.002e-01 | 2.002e-01 | 5.646e+00 | 5.646e+00 | 2.027e+01 | 5.911e+00 | 10.54 | 10.54 | +0.00 | FAIL | same |
| armD_D1 | startup-error/startup-error | -/- | 0/0 | -/- | - | - | - | - | - | - | - | - | - | - | - | - | - | - |
| armD_D2 | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 1.874e+00 | 1.874e+00 | 7.922e-01 | 7.922e-01 | 8.706e-01 | 8.706e-01 | 8.029e-01 | 9.280e-01 | 10.42 | 10.42 | +0.00 | FAIL | same |
| armD_D2_LW | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 2.392e+00 | 2.392e+00 | 1.059e+00 | 1.059e+00 | 1.531e+00 | 1.531e+00 | 1.193e+00 | 1.407e+00 | 10.20 | 10.20 | +0.00 | FAIL | same |
| armD_D2_LW_newton | maxsteps/maxsteps | 2/2 | 0/0 | 104/104 | 2.00 | 1.546e+01 | 1.546e+01 | 2.225e+01 | 2.225e+01 | 3.283e+01 | 3.283e+01 | 2.644e+01 | 3.718e+01 | 8.88 | 8.88 | +0.00 | FAIL | same |
| armD_D2_newton | du/du | 2/2 | 8/8 | 214/214 | 2.00 | 9.994e-04 | 9.994e-04 | 1.391e-03 | 1.391e-03 | 1.399e-03 | 1.399e-03 | 1.396e-03 | 1.400e-03 | 16.61 | 16.61 | +0.00 | pass | same |
| armD_D2_newton_bigstack | du/du | 2/2 | 8/8 | 214/214 | 2.00 | 9.994e-04 | 9.994e-04 | 1.391e-03 | 1.391e-03 | 1.399e-03 | 1.399e-03 | 1.396e-03 | 1.400e-03 | 16.61 | 16.61 | +0.00 | pass | same |
| armD_D3 | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 1.944e+00 | 1.944e+00 | 9.592e-01 | 9.592e-01 | 9.203e-01 | 9.203e-01 | 9.278e-01 | 9.114e-01 | 10.40 | 10.40 | +0.00 | FAIL | same |
| armD_D4a | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 1.813e+00 | 1.813e+00 | 7.714e-01 | 7.714e-01 | 1.021e+00 | 1.021e+00 | 7.971e-01 | 1.506e+00 | 10.48 | 10.48 | +0.00 | FAIL | same |
| armD_D4b | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 2.767e+01 | 2.767e+01 | 3.471e+01 | 3.471e+01 | 2.590e+01 | 2.590e+01 | 2.971e+01 | 2.438e+01 | 7.74 | 7.74 | +0.00 | FAIL | same |
| armHeH_0p3 | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 1.586e+00 | 1.586e+00 | 7.735e-01 | 7.735e-01 | 8.078e-01 | 8.078e-01 | 7.775e-01 | 8.344e-01 | 10.47 | 10.47 | +0.00 | FAIL | same |
| armHeH_10 | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 2.866e+00 | 2.866e+00 | 8.833e-01 | 8.833e-01 | 8.645e-01 | 8.645e-01 | 8.661e-01 | 8.666e-01 | 10.14 | 10.14 | +0.00 | FAIL | same |
| arm_heh1_x2matched | maxsteps/jfnk | 2/0 | 40/106 | 1127/89 | 2.00 | 4.696e-02 | 2.678e-04 | 3.316e-01 | 2.718e-04 | 7.369e-01 | 8.357e-03 | 1.339e-03 | 1.352e-02 | 10.50 | 10.31 | -0.19 | FAIL | MOVED |
| armHeH_3 | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 2.284e+00 | 2.284e+00 | 8.714e-01 | 8.714e-01 | 1.283e+00 | 1.283e+00 | 9.181e-01 | 1.871e+00 | 10.12 | 10.12 | +0.00 | FAIL | same |
| armHeH_30 | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 3.102e+00 | 3.102e+00 | 9.006e-01 | 9.006e-01 | 9.001e-01 | 9.001e-01 | 9.217e-01 | 8.977e-01 | 10.13 | 10.13 | +0.00 | FAIL | same |
| jfnk_hd189 | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 1.555e-02 | 1.555e-02 | 2.128e-01 | 2.128e-01 | 5.873e+00 | 5.873e+00 | 1.698e+00 | 4.525e+00 | 9.93 | 9.93 | +0.00 | FAIL | same |
| jfnk_hd189_tight | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 1.542e-02 | 1.542e-02 | 1.106e-01 | 1.106e-01 | 4.115e+00 | 4.115e+00 | 8.400e-01 | 3.820e+00 | 10.01 | 10.01 | +0.00 | FAIL | same |
| lower_profile | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 2.823e+01 | 2.823e+01 | 1.858e+00 | 1.858e+00 | 5.618e+00 | 5.618e+00 | 4.870e+00 | 5.893e+00 | 9.56 | 9.56 | +0.00 | FAIL | same |
| mol_base_handoff | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 1.147e+00 | 1.147e+00 | 7.966e-01 | 7.966e-01 | 7.926e-01 | 7.926e-01 | 7.806e-01 | 8.039e-01 | 10.43 | 10.43 | +0.00 | FAIL | same |
| mol_carrier | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 5.904e-01 | 5.904e-01 | 6.615e-01 | 6.615e-01 | 6.982e-01 | 6.982e-01 | 6.620e-01 | 7.118e-01 | 10.42 | 10.42 | +0.00 | FAIL | same |
| mol_diffusion | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 1.147e+00 | 1.147e+00 | 7.967e-01 | 7.967e-01 | 7.926e-01 | 7.926e-01 | 7.807e-01 | 8.039e-01 | 10.43 | 10.43 | +0.00 | FAIL | same |
| mol_ir_bands | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 1.276e+00 | 1.276e+00 | 8.196e-01 | 8.196e-01 | 7.768e-01 | 7.768e-01 | 7.961e-01 | 7.675e-01 | 10.43 | 10.43 | +0.00 | FAIL | same |
| mol_lyman_werner | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 5.504e-01 | 5.504e-01 | 8.423e-01 | 8.423e-01 | 1.106e+00 | 1.106e+00 | 1.073e+00 | 1.033e+00 | 10.47 | 10.47 | +0.00 | FAIL | same |
| mol_metals | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 1.275e+00 | 1.275e+00 | 8.197e-01 | 8.197e-01 | 7.753e-01 | 7.753e-01 | 7.960e-01 | 7.651e-01 | 10.43 | 10.43 | +0.00 | FAIL | same |
| mol_sec_ion | maxsteps/maxsteps | -/- | 0/0 | -/- | 2.00 | 4.096e-01 | 4.096e-01 | 9.750e-01 | 9.750e-01 | 1.191e+01 | 1.191e+01 | 6.050e+00 | 5.265e+00 | 10.38 | 10.38 | +0.00 | FAIL | same |
| newton_rsw01 | jfnk/jfnk | 0/0 | 9/76 | 205/1679 | 1.50 | 5.031e-05 | 5.030e-05 | 1.535e-04 | 5.242e-05 | 2.044e-04 | 6.751e-05 | 5.242e-05 | 8.432e-05 | 13.22 | 13.21 | -0.01 | pass | MOVED |
| newton_rsw05 | jfnk/jfnk | 0/0 | 20/71 | 471/1594 | 1.50 | 5.043e-05 | 5.029e-05 | 1.703e-04 | 5.245e-05 | 3.537e-04 | 6.751e-05 | 5.245e-05 | 8.433e-05 | 13.22 | 13.21 | -0.01 | pass | MOVED |
| ptc_warm | -/maxsteps | -/- | 0/0 | -/- | 1.50 | 4.691e-04 | 2.233e-03 | 4.900e-03 | 2.488e-02 | 8.236e-02 | 3.433e-01 | 1.043e-01 | 3.551e-01 | 13.22 | 13.24 | +0.02 | FAIL | MOVED |
| wasp_full_newton [1] | jfnk/jfnk | 0/0 | 14/44 | 338/1000 | 1.50 | 8.099e-05 | 5.030e-05 | 1.357e-03 | 5.242e-05 | 5.469e-02 | 6.751e-05 | 5.242e-05 | 8.432e-05 | 13.21 | 13.21 | +0.00 | pass | MOVED |
| wasp_full | du/du | -/- | 0/0 | -/- | 1.50 | 9.999e-04 | 9.999e-04 | 3.129e-02 | 3.129e-02 | 9.430e-01 | 9.430e-01 | 2.996e-01 | 1.372e+00 | 13.26 | 13.26 | +0.00 | FAIL | same |
| wasp_full_newton | jfnk/jfnk | 0/0 | 14/44 | 338/1000 | 1.50 | 8.099e-05 | 5.030e-05 | 1.357e-03 | 5.242e-05 | 5.469e-02 | 6.751e-05 | 5.242e-05 | 8.432e-05 | 13.21 | 13.21 | +0.00 | pass | MOVED |
| wasp_he23off | du/du | -/- | 0/0 | -/- | 1.50 | 9.887e-04 | 9.887e-04 | 3.124e-02 | 3.124e-02 | 9.573e-01 | 9.573e-01 | 3.006e-01 | 1.391e+00 | 13.26 | 13.26 | +0.00 | FAIL | same |
| wasp_he23off_newton | jfnk/jfnk | 0/0 | 11/31 | 281/722 | 1.50 | 4.941e-05 | 5.020e-05 | 1.226e-03 | 5.292e-05 | 2.844e-02 | 6.790e-05 | 5.292e-05 | 8.461e-05 | 13.21 | 13.22 | +0.01 | pass | MOVED |

---

[1] Recorded during this census as `solver_newton_cold`. That directory was
byte-identical to `backup/regression/wasp_full_newton` -- same `input.inp` and
`metals.inp`, only the stored outputs differed -- and was deleted on 2026-09-03
(user decision; `TO_BE_DONE.md` (P53)). The row is the same case and the numbers
in it stand; only the directory name it was recorded under is gone. Note that
`wasp_full_newton` also appears in this census under its own name where the run
was made from that directory.
