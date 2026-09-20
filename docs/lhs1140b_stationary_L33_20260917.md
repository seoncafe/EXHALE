# LHS 1140 b stationary series, item L33: slow outer relaxation or a stalled mode, read from the existing logs

Written 2026-09-17 (KST) for item L33 of `docs/PLAN_20260917.md` section 7.
No solve was run for it. Every number below is MEASURED (computed here from
a run log or from a case input file) or READ (quoted from a memo, from the
revised handoff or from a source line), and each is labeled. The reader is
`LHS1140b/models/.L22/outer_pass_history.py`, written for this item; it
parses the outer-pass lines, the certification blocks and the carrier
relaxation lines of a stationary run log, takes the radius of a cell from
the grid of the state the run loaded, and forms the geometric factor of the
worst gated row from pass to pass together with its running mean. It reads
and writes nothing but its own table.

The logs read: `LHS1140b/models/.L22/i4_wm0083/run.log`,
`i4_wm055/run.log`, `i4_kz0083/run.log`, their `run_lart3.log`
predecessors, and `i3_alt/run.log`, the alternation that certified.

---

## 1. Verdict

| run | case | classification |
|---|---|---|
| `i4_wm0083` | `molecular_scalar_gj1132_wellmixed/HeH0.083` | **front relaxation.** The refusing cell is fixed only because it is the last cell of the grid; the two H2 contours move monotonically through the whole run, the outer one inward from 28.55 to 9.13 R_p and the inner one outward from 1.2565 to 1.7719 R_p, and the row decays once the outer contour is inside the domain. |
| `i4_wm055` | `molecular_scalar_gj1132_wellmixed/HeH0.55` | **front relaxation**, unambiguous: the refusing cell itself walks outward from cell 315 (2.1136 R_p) to cell 445 (11.7471 R_p) over sixteen passes and the row decays from pass 5 on. The caveat is that the wind under it has stopped converging: `info = 2` at six of sixteen passes and the hydrodynamic mass row reaches 4.36e-01. |
| `i4_kz0083` | `molecular_scalar_gj1132_kzz1e9/HeH0.083` | **neither of the plan's two**, and this is the correction the logs force. The refusing cell is fixed at 499 from pass 9 to pass 40 and both H2 contours stand still, which reads as a mode; but the row is not bounded, it falls, and its rate is set by the composition movement bound and not by any property of the alternation map. It is a bound-throttled relaxation of a row hosted by the outermost cells. |

The one bounded continuation experiment the plan allows is designed in
section 7 for the two front-relaxation cases. L32's verdict does not apply
to any of the three: none of them entered the coupled block, so none of them
met the linear-solve failure L32 measures.

---

## 2. The control the three cases have to be read against

`i3_alt` solved `molecular_scalar_gj1132_kzz1e9/HeH2.13` by the same
alternation and certified at outer pass 12 (READ,
`docs/lhs1140b_stationary_L22_20260916.md` step 3 section 2; MEASURED again
here from the log). Three of its quantities separate it from all three
refusing runs, and they are the ones this item turns on (MEASURED):

| quantity | `i3_alt`, certified | `i4_wm0083` | `i4_wm055` | `i4_kz0083` |
|---|---|---|---|---|
| how every carrier relaxation ended | the carriers stopped moving at the fixed wind, 11 of 11 passes with an update | the composition movement bound, 37 of 37 | the movement bound, 16 of 16 | the movement bound, 39 of 39 |
| carrier residual of the returned composition | 2.87e-09 to 4.70e-09 at every pass | 2.03e-01 falling to 6.60e-02 at pass 20, then rising to 8.28e-02 at pass 37 | 2.97e-01 falling to about 1.7e-01 and flat from pass 8 | 1.72e-02 to 1.86e-02, flat over all 39 |
| geometric factor of the worst gated row | 0.3526 a pass over passes 1 to 12 | 1.0611 a pass over 1 to 22, then 0.9527 over 26 to 37 | 0.9522 a pass over 1 to 16 | 0.9886 a pass over 20 to 40 |
| both H2 contours | fixed from pass 2 (x2 = 0.5 at 1.0002 R_p, x2 = 1e-2 at 1.1626 R_p) | both moving every pass | both moving every pass | fixed (1.3705 to 1.3770 R_p; the 1e-2 contour never inside the grid) |

The run that certifies is the run whose carrier relaxation reaches the fixed
point of its own operator at the held wind. In all three refusing runs the
relaxation is cut off by the movement bound at every pass, and the
composition handed to the next wind solve is one to five decades away from
that fixed point. That is the single fact the rest of this memo rests on.

---

## 3. The pass tables

`factor` is the ratio of this pass's worst gated row to the previous pass's;
`running mean` is the geometric mean of every factor up to that pass.
`bound cell` and `bound r` are the cell at which the particle-count movement
bound was refused. `trust` is the movement bound in force. `x2 = 0.5` and
`x2 = 1e-2` are the two H2 contours the pass printed; `> r(N)` means the
contour lies beyond the outer boundary at 30 R_p. The worst gated row is the
carrier H2 balance in every pass of every run.


#### `i4_wm0083` (`wellmixed/HeH0.083`)

| pass | worst gated row | cell | r [R_p] | factor | running mean | hydro info | hydro mass row | x2=0.5 [R_p] | x2=1e-2 [R_p] | bound cell | bound r | trust | wall [s] | carrier residual |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2.650e-02 | 500 | 29.0312 | - | - | 0 | 9.12e-09 | 1.2565 | > r(N) | 221 | 1.2153 | 1.0e-02 | 135.9 | 2.030e-01 |
| 2 | 3.850e-02 | 500 | 29.0312 | 1.4528 | 1.4528 | 2 | 3.34e-08 | 1.2657 | > r(N) | 224 | 1.2269 | 1.0e-02 | 930.1 | 1.930e-01 |
| 3 | 4.110e-02 | 500 | 29.0312 | 1.0675 | 1.2454 | 0 | 1.27e-08 | 1.2751 | > r(N) | 227 | 1.2392 | 1.0e-02 | 574.7 | 1.820e-01 |
| 4 | 4.410e-02 | 500 | 29.0312 | 1.0730 | 1.1850 | 0 | 1.50e-08 | 1.2850 | > r(N) | 230 | 1.2521 | 1.0e-02 | 432.6 | 1.720e-01 |
| 5 | 4.650e-02 | 500 | 29.0312 | 1.0544 | 1.1509 | 0 | 2.62e-08 | 1.2951 | > r(N) | 233 | 1.2657 | 1.0e-02 | 1058.6 | 1.620e-01 |
| 6 | 4.930e-02 | 500 | 29.0312 | 1.0602 | 1.1322 | 0 | 1.56e-08 | 1.3056 | > r(N) | 236 | 1.2800 | 1.0e-02 | 273.2 | 1.510e-01 |
| 7 | 5.230e-02 | 500 | 29.0312 | 1.0609 | 1.1200 | 0 | 9.70e-09 | 1.3165 | 28.5489 | 240 | 1.3003 | 1.0e-02 | 642.6 | 1.400e-01 |
| 8 | 5.540e-02 | 500 | 29.0312 | 1.0593 | 1.1111 | 0 | 1.51e-08 | 1.3278 | 27.1453 | 243 | 1.3165 | 1.0e-02 | 659.6 | 1.300e-01 |
| 9 | 5.870e-02 | 500 | 29.0312 | 1.0596 | 1.1045 | 0 | 1.24e-08 | 1.3395 | 25.8132 | 247 | 1.3395 | 1.0e-02 | 665.7 | 1.200e-01 |
| 10 | 6.170e-02 | 500 | 29.0312 | 1.0511 | 1.0985 | 0 | 3.92e-08 | 1.3516 | 24.1420 | 253 | 1.3770 | 1.0e-02 | 193.2 | 1.110e-01 |
| 11 | 6.490e-02 | 500 | 29.0312 | 1.0519 | 1.0937 | 0 | 2.84e-08 | 1.3641 | 22.9629 | 257 | 1.4044 | 1.0e-02 | 180.1 | 1.020e-01 |
| 12 | 6.790e-02 | 500 | 29.0312 | 1.0462 | 1.0893 | 0 | 1.12e-08 | 1.3770 | 21.4837 | 262 | 1.4413 | 1.0e-02 | 170.2 | 9.440e-02 |
| 13 | 7.080e-02 | 500 | 29.0312 | 1.0427 | 1.0853 | 0 | 2.62e-08 | 1.3905 | 20.1041 | 266 | 1.4733 | 1.0e-02 | 158.2 | 8.750e-02 |
| 14 | 7.370e-02 | 500 | 29.0312 | 1.0410 | 1.0819 | 0 | 1.07e-08 | 1.4044 | 18.8174 | 269 | 1.4987 | 1.0e-02 | 141.3 | 8.160e-02 |
| 15 | 7.670e-02 | 500 | 29.0312 | 1.0407 | 1.0789 | 2 | 3.97e-08 | 1.4188 | 17.6173 | 273 | 1.5349 | 1.0e-02 | 129.9 | 7.670e-02 |
| 16 | 7.930e-02 | 500 | 29.0312 | 1.0339 | 1.0758 | 0 | 1.80e-08 | 1.4337 | 16.4981 | 276 | 1.5636 | 1.0e-02 | 141.1 | 7.280e-02 |
| 17 | 8.190e-02 | 500 | 29.0312 | 1.0328 | 1.0731 | 2 | 8.45e-08 | 1.4570 | 15.4543 | 279 | 1.5940 | 1.0e-02 | 128.3 | 6.980e-02 |
| 18 | 8.440e-02 | 500 | 29.0312 | 1.0305 | 1.0705 | 2 | 3.68e-08 | 1.4733 | 14.2478 | 282 | 1.6259 | 1.0e-02 | 96.7 | 6.780e-02 |
| 19 | 8.680e-02 | 500 | 29.0312 | 1.0284 | 1.0681 | 0 | 2.41e-08 | 1.4987 | 13.3555 | 285 | 1.6596 | 1.0e-02 | 360.6 | 6.650e-02 |
| 20 | 8.880e-02 | 500 | 29.0312 | 1.0230 | 1.0657 | 0 | 6.61e-08 | 1.5165 | 12.5233 | 287 | 1.6830 | 1.0e-02 | 86.4 | 6.600e-02 |
| 21 | 9.060e-02 | 500 | 29.0312 | 1.0203 | 1.0634 | 0 | 6.59e-08 | 1.5443 | 11.7471 | 289 | 1.7073 | 1.0e-02 | 121.6 | 6.620e-02 |
| 22 | 9.200e-02 | 500 | 29.0312 | 1.0155 | 1.0611 | 0 | 5.70e-08 | 1.5539 | 11.3789 | 290 | 1.7198 | 5.0e-03 | 87.8 | 6.920e-02 |
| 23 | 8.230e-02 | 500 | 29.0312 | 0.8946 | 1.0529 | 0 | 3.81e-08 | 1.5636 | 11.1995 | 292 | 1.7454 | 5.0e-03 | 64.0 | 6.900e-02 |
| 24 | 8.220e-02 | 500 | 29.0312 | 0.9988 | 1.0504 | 0 | 2.86e-08 | 1.5736 | 10.8501 | 293 | 1.7585 | 5.0e-03 | 63.9 | 6.890e-02 |
| 25 | 8.200e-02 | 500 | 29.0312 | 0.9976 | 1.0482 | 0 | 4.36e-08 | 1.5940 | 10.5125 | 294 | 1.7719 | 5.0e-03 | 73.1 | 6.960e-02 |
| 26 | 8.160e-02 | 500 | 29.0312 | 0.9951 | 1.0460 | 0 | 6.12e-08 | 1.6044 | 10.3481 | 294 | 1.7719 | 5.0e-03 | 64.8 | 7.030e-02 |
| 27 | 8.100e-02 | 500 | 29.0312 | 0.9926 | 1.0439 | 0 | 4.03e-08 | 1.6151 | 10.0278 | 295 | 1.7854 | 5.0e-03 | 79.2 | 7.110e-02 |
| 28 | 8.000e-02 | 500 | 29.0312 | 0.9877 | 1.0418 | 0 | 4.97e-08 | 1.6370 | 9.8718 | 296 | 1.7993 | 5.0e-03 | 78.9 | 7.210e-02 |
| 29 | 7.860e-02 | 500 | 29.0312 | 0.9825 | 1.0396 | 0 | 4.17e-08 | 1.6482 | 9.7185 | 297 | 1.8134 | 5.0e-03 | 64.5 | 7.310e-02 |
| 30 | 7.680e-02 | 500 | 29.0312 | 0.9771 | 1.0374 | 0 | 7.88e-08 | 1.6596 | 9.5678 | 298 | 1.8277 | 5.0e-03 | 63.7 | 7.410e-02 |
| 31 | 7.440e-02 | 500 | 29.0312 | 0.9688 | 1.0350 | 0 | 3.09e-08 | 1.6712 | 9.4197 | 299 | 1.8423 | 5.0e-03 | 68.8 | 7.530e-02 |
| 32 | 7.150e-02 | 500 | 29.0312 | 0.9610 | 1.0325 | 0 | 1.47e-08 | 1.6951 | 9.2742 | 300 | 1.8571 | 5.0e-03 | 78.3 | 7.640e-02 |
| 33 | 6.800e-02 | 500 | 29.0312 | 0.9510 | 1.0299 | 0 | 5.55e-08 | 1.7073 | 9.2742 | 301 | 1.8722 | 5.0e-03 | 68.0 | 7.760e-02 |
| 34 | 6.380e-02 | 500 | 29.0312 | 0.9382 | 1.0270 | 0 | 2.09e-08 | 1.7325 | 9.1312 | 302 | 1.8875 | 5.0e-03 | 84.4 | 7.890e-02 |
| 35 | 5.910e-02 | 500 | 29.0312 | 0.9263 | 1.0239 | 0 | 4.10e-07 | 1.7454 | 9.1312 | 303 | 1.9032 | 5.0e-03 | 54.2 | 8.010e-02 |
| 36 | 5.380e-02 | 500 | 29.0312 | 0.9103 | 1.0204 | 0 | 4.39e-08 | 1.7585 | 9.1312 | 304 | 1.9191 | 5.0e-03 | 4928.5 | 8.140e-02 |
| 37 | 4.790e-02 | 500 | 29.0312 | 0.8903 | 1.0166 | 2 | 3.40e-03 | 1.7719 | 9.1312 | 304 | 1.9191 | 2.5e-03 | 5571.1 | 8.280e-02 |

#### `i4_wm055` (`wellmixed/HeH0.55`)

| pass | worst gated row | cell | r [R_p] | factor | running mean | hydro info | hydro mass row | x2=0.5 [R_p] | x2=1e-2 [R_p] | bound cell | bound r | trust | wall [s] | carrier residual |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 5.920e-02 | 315 | 2.11360 | - | - | 0 | 3.01e-08 | 1.1626 | 1.5539 | 205 | 1.1626 | 1.0e-02 | 743.2 | 2.970e-01 |
| 2 | 3.240e-02 | 326 | 2.34925 | 0.5473 | 0.5473 | 0 | 3.77e-08 | 1.1684 | 1.6044 | 221 | 1.2153 | 1.0e-02 | 1311.3 | 2.620e-01 |
| 3 | 6.130e-02 | 333 | 2.52451 | 1.8920 | 1.0176 | 0 | 2.27e-08 | 1.1744 | 1.6830 | 228 | 1.2434 | 1.0e-02 | 1142.7 | 2.350e-01 |
| 4 | 8.050e-02 | 341 | 2.75283 | 1.3132 | 1.1079 | 0 | 2.82e-08 | 1.1806 | 1.7719 | 241 | 1.3056 | 1.0e-02 | 802.3 | 2.140e-01 |
| 5 | 9.010e-02 | 350 | 3.05077 | 1.1193 | 1.1107 | 0 | 3.72e-08 | 1.1871 | 1.8875 | 287 | 1.6830 | 1.0e-02 | 1362.7 | 2.010e-01 |
| 6 | 9.520e-02 | 359 | 3.39932 | 1.0566 | 1.0997 | 0 | 3.51e-08 | 1.1904 | 2.0206 | 291 | 1.7325 | 1.0e-02 | 220.5 | 1.910e-01 |
| 7 | 8.950e-02 | 369 | 3.85645 | 0.9401 | 1.0713 | 2 | 1.04e-01 | 1.1938 | 2.1735 | 294 | 1.7719 | 1.0e-02 | 3976.1 | 1.930e-01 |
| 8 | 8.160e-02 | 379 | 4.40062 | 0.9117 | 1.0469 | 0 | 3.05e-08 | 1.2007 | 2.3259 | 297 | 1.8134 | 1.0e-02 | 460.4 | 1.800e-01 |
| 9 | 7.130e-02 | 389 | 5.04842 | 0.8738 | 1.0235 | 0 | 4.00e-08 | 1.2042 | 2.5245 | 300 | 1.8571 | 1.0e-02 | 408.7 | 1.780e-01 |
| 10 | 6.380e-02 | 400 | 5.90433 | 0.8948 | 1.0083 | 2 | 3.05e-03 | 1.2079 | 2.7528 | 302 | 1.8875 | 1.0e-02 | 1063.0 | 1.780e-01 |
| 11 | 5.490e-02 | 411 | 6.94115 | 0.8605 | 0.9925 | 0 | 3.96e-08 | 1.2115 | 3.0153 | 304 | 1.9191 | 1.0e-02 | 544.8 | 1.760e-01 |
| 12 | 4.560e-02 | 423 | 8.32366 | 0.8306 | 0.9766 | 0 | 8.52e-08 | 1.2153 | 3.3578 | 306 | 1.9517 | 1.0e-02 | 105.7 | 1.770e-01 |
| 13 | 3.740e-02 | 435 | 10.0278 | 0.8202 | 0.9625 | 2 | 5.72e-03 | 1.2153 | 3.5282 | 307 | 1.9685 | 5.0e-03 | 2114.0 | 1.820e-01 |
| 14 | 3.310e-02 | 441 | 11.0232 | 0.8850 | 0.9563 | 2 | 4.36e-01 | 1.2153 | 3.6179 | 307 | 1.9685 | 2.5e-03 | 1040.2 | 1.790e-01 |
| 15 | 3.060e-02 | 444 | 11.5614 | 0.9245 | 0.9540 | 2 | 1.98e-01 | 1.2153 | 3.7109 | 307 | 1.9685 | 2.5e-03 | 1171.0 | 1.590e-01 |
| 16 | 2.840e-02 | 445 | 11.7471 | 0.9281 | 0.9522 | 2 | 9.89e-02 | 1.2191 | 3.8071 | 308 | 1.9855 | 2.5e-03 | 1667.9 | 1.690e-01 |

#### `i4_kz0083` (`kzz1e9/HeH0.083`)

| pass | worst gated row | cell | r [R_p] | factor | running mean | hydro info | hydro mass row | x2=0.5 [R_p] | x2=1e-2 [R_p] | bound cell | bound r | trust | wall [s] | carrier residual |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2.970e-03 | 500 | 29.0312 | - | - | 2 | 3.69e-03 | 1.3705 | > r(N) | 96 | 1.0232 | 1.0e-02 | 212.7 | 1.720e-02 |
| 2 | 2.520e-02 | 500 | 29.0312 | 8.4848 | 8.4848 | 0 | 5.25e-09 | 1.3705 | > r(N) | 470 | 17.6173 | 1.0e-02 | 610.4 | 1.850e-02 |
| 3 | 7.980e-03 | 500 | 29.0312 | 0.3167 | 1.6392 | 0 | 1.19e-08 | 1.3641 | > r(N) | 461 | 15.2045 | 1.0e-02 | 139.7 | 1.860e-02 |
| 4 | 3.350e-03 | 499 | 28.5489 | 0.4198 | 1.0409 | 0 | 7.48e-09 | 1.3641 | > r(N) | 453 | 13.3555 | 1.0e-02 | 113.3 | 1.850e-02 |
| 5 | 4.100e-03 | 500 | 29.0312 | 1.2239 | 1.0839 | 0 | 8.04e-09 | 1.3641 | > r(N) | 444 | 11.5614 | 1.0e-02 | 128.7 | 1.830e-02 |
| 6 | 5.610e-03 | 500 | 29.0312 | 1.3683 | 1.1356 | 0 | 1.21e-08 | 1.3641 | > r(N) | 445 | 11.7471 | 5.0e-03 | 668.5 | 1.830e-02 |
| 7 | 2.110e-03 | 500 | 29.0312 | 0.3761 | 0.9446 | 0 | 7.88e-09 | 1.3641 | > r(N) | 442 | 11.1995 | 5.0e-03 | 141.3 | 1.820e-02 |
| 8 | 2.260e-03 | 500 | 29.0312 | 1.0711 | 0.9617 | 0 | 8.62e-09 | 1.3641 | > r(N) | 442 | 11.1995 | 2.5e-03 | 151.7 | 1.820e-02 |
| 9 | 1.370e-03 | 499 | 28.5489 | 0.6062 | 0.9078 | 0 | 8.02e-09 | 1.3641 | > r(N) | 440 | 10.8501 | 2.5e-03 | 476.9 | 1.820e-02 |
| 10 | 1.310e-03 | 499 | 28.5489 | 0.9562 | 0.9131 | 0 | 2.36e-08 | 1.3641 | > r(N) | 439 | 10.6798 | 2.5e-03 | 452.6 | 1.820e-02 |
| 11 | 1.250e-03 | 499 | 28.5489 | 0.9542 | 0.9171 | 0 | 1.88e-08 | 1.3641 | > r(N) | 439 | 10.6798 | 2.5e-03 | 188.7 | 1.810e-02 |
| 12 | 1.230e-03 | 499 | 28.5489 | 0.9840 | 0.9230 | 0 | 5.13e-09 | 1.3641 | > r(N) | 438 | 10.5125 | 2.5e-03 | 187.2 | 1.810e-02 |
| 13 | 1.190e-03 | 499 | 28.5489 | 0.9675 | 0.9266 | 0 | 1.56e-08 | 1.3705 | > r(N) | 438 | 10.5125 | 2.5e-03 | 665.6 | 1.810e-02 |
| 14 | 1.150e-03 | 499 | 28.5489 | 0.9664 | 0.9296 | 0 | 4.25e-09 | 1.3705 | > r(N) | 437 | 10.3481 | 2.5e-03 | 222.5 | 1.800e-02 |
| 15 | 1.120e-03 | 499 | 28.5489 | 0.9739 | 0.9327 | 2 | 1.79e-08 | 1.3705 | > r(N) | 436 | 10.1866 | 2.5e-03 | 174.8 | 1.800e-02 |
| 16 | 1.110e-03 | 499 | 28.5489 | 0.9911 | 0.9365 | 0 | 1.11e-08 | 1.3705 | > r(N) | 436 | 10.1866 | 2.5e-03 | 505.6 | 1.800e-02 |
| 17 | 1.080e-03 | 499 | 28.5489 | 0.9730 | 0.9387 | 0 | 9.68e-09 | 1.3705 | > r(N) | 435 | 10.0278 | 2.5e-03 | 200.7 | 1.800e-02 |
| 18 | 1.060e-03 | 499 | 28.5489 | 0.9815 | 0.9412 | 0 | 1.06e-08 | 1.3705 | > r(N) | 435 | 10.0278 | 2.5e-03 | 138.5 | 1.790e-02 |
| 19 | 1.030e-03 | 499 | 28.5489 | 0.9717 | 0.9429 | 0 | 1.26e-08 | 1.3705 | > r(N) | 434 | 9.8718 | 2.5e-03 | 136.0 | 1.790e-02 |
| 20 | 1.030e-03 | 499 | 28.5489 | 1.0000 | 0.9458 | 0 | 1.24e-08 | 1.3705 | > r(N) | 434 | 9.8718 | 2.5e-03 | 136.6 | 1.790e-02 |
| 21 | 9.930e-04 | 499 | 28.5489 | 0.9641 | 0.9467 | 0 | 5.52e-09 | 1.3705 | > r(N) | 433 | 9.7185 | 2.5e-03 | 135.3 | 1.780e-02 |
| 22 | 9.930e-04 | 499 | 28.5489 | 1.0000 | 0.9492 | 0 | 1.07e-08 | 1.3705 | > r(N) | 435 | 10.0278 | 1.3e-03 | 134.9 | 1.780e-02 |
| 23 | 9.790e-04 | 499 | 28.5489 | 0.9859 | 0.9508 | 0 | 1.58e-08 | 1.3705 | > r(N) | 434 | 9.8718 | 1.3e-03 | 170.2 | 1.780e-02 |
| 24 | 9.650e-04 | 499 | 28.5489 | 0.9857 | 0.9523 | 0 | 1.07e-08 | 1.3705 | > r(N) | 433 | 9.7185 | 1.3e-03 | 158.1 | 1.780e-02 |
| 25 | 9.500e-04 | 499 | 28.5489 | 0.9845 | 0.9536 | 0 | 1.42e-08 | 1.3705 | > r(N) | 433 | 9.7185 | 1.3e-03 | 157.7 | 1.780e-02 |
| 26 | 9.330e-04 | 499 | 28.5489 | 0.9821 | 0.9547 | 0 | 4.50e-09 | 1.3705 | > r(N) | 432 | 9.5678 | 1.3e-03 | 155.8 | 1.780e-02 |
| 27 | 9.240e-04 | 499 | 28.5489 | 0.9904 | 0.9561 | 0 | 1.00e-08 | 1.3705 | > r(N) | 432 | 9.5678 | 1.3e-03 | 155.5 | 1.780e-02 |
| 28 | 9.090e-04 | 499 | 28.5489 | 0.9838 | 0.9571 | 0 | 1.02e-08 | 1.3770 | > r(N) | 432 | 9.5678 | 1.3e-03 | 111.6 | 1.770e-02 |
| 29 | 9.070e-04 | 499 | 28.5489 | 0.9978 | 0.9585 | 0 | 8.62e-09 | 1.3770 | > r(N) | 431 | 9.4197 | 1.3e-03 | 112.0 | 1.770e-02 |
| 30 | 8.950e-04 | 499 | 28.5489 | 0.9868 | 0.9595 | 0 | 1.28e-08 | 1.3770 | > r(N) | 431 | 9.4197 | 1.3e-03 | 119.8 | 1.770e-02 |
| 31 | 8.830e-04 | 499 | 28.5489 | 0.9866 | 0.9604 | 0 | 8.64e-09 | 1.3770 | > r(N) | 431 | 9.4197 | 1.3e-03 | 296.2 | 1.770e-02 |
| 32 | 8.730e-04 | 499 | 28.5489 | 0.9887 | 0.9613 | 0 | 6.78e-09 | 1.3770 | > r(N) | 431 | 9.4197 | 1.3e-03 | 126.2 | 1.770e-02 |
| 33 | 8.630e-04 | 499 | 28.5489 | 0.9885 | 0.9621 | 0 | 7.59e-09 | 1.3770 | > r(N) | 430 | 9.2742 | 1.3e-03 | 144.4 | 1.770e-02 |
| 34 | 8.600e-04 | 499 | 28.5489 | 0.9965 | 0.9631 | 2 | 1.47e-08 | 1.3770 | > r(N) | 430 | 9.2742 | 1.3e-03 | 331.2 | 1.760e-02 |
| 35 | 8.490e-04 | 499 | 28.5489 | 0.9872 | 0.9638 | 0 | 8.10e-09 | 1.3770 | > r(N) | 430 | 9.2742 | 1.3e-03 | 257.2 | 1.760e-02 |
| 36 | 8.410e-04 | 499 | 28.5489 | 0.9906 | 0.9646 | 0 | 2.20e-09 | 1.3770 | > r(N) | 430 | 9.2742 | 1.3e-03 | 369.4 | 1.760e-02 |
| 37 | 8.330e-04 | 499 | 28.5489 | 0.9905 | 0.9653 | 0 | 4.75e-09 | 1.3770 | > r(N) | 429 | 9.1312 | 1.3e-03 | 375.0 | 1.760e-02 |
| 38 | 8.280e-04 | 499 | 28.5489 | 0.9940 | 0.9661 | 0 | 7.43e-09 | 1.3770 | > r(N) | 429 | 9.1312 | 1.3e-03 | 120.5 | 1.760e-02 |
| 39 | 8.250e-04 | 499 | 28.5489 | 0.9964 | 0.9669 | 0 | 8.31e-09 | 1.3770 | > r(N) | 429 | 9.1312 | 1.3e-03 | 102.1 | 1.750e-02 |
| 40 | 8.190e-04 | 499 | 28.5489 | 0.9927 | 0.9675 | 0 | 8.80e-09 | 1.3770 | > r(N) | - | - | 1.3e-03 | 233.7 | not measured |

#### `i3_alt` (`kzz1e9/HeH2.13`)

| pass | worst gated row | cell | r [R_p] | factor | running mean | hydro info | hydro mass row | x2=0.5 [R_p] | x2=1e-2 [R_p] | bound cell | bound r | trust | wall [s] | carrier residual |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 7.820e-01 | 243 | 1.31651 | - | - | 2 | 2.47e-06 | 1.0002 | 1.1597 | - | - | 1.0e-02 | 2691.4 | 4.430e-09 |
| 2 | 5.300e-01 | 221 | 1.21527 | 0.6777 | 0.6777 | 0 | 4.44e-09 | 1.0002 | 1.1626 | - | - | 1.0e-02 | 1725.2 | 4.700e-09 |
| 3 | 1.490e-01 | 221 | 1.21527 | 0.2811 | 0.4365 | 0 | 6.82e-09 | 1.0002 | 1.1626 | - | - | 1.0e-02 | 1098.2 | 4.010e-09 |
| 4 | 2.280e-02 | 222 | 1.21908 | 0.1530 | 0.3078 | 0 | 6.49e-09 | 1.0002 | 1.1626 | - | - | 1.0e-02 | 1008.5 | 3.890e-09 |
| 5 | 2.490e-03 | 242 | 1.31102 | 0.1092 | 0.2375 | 0 | 4.17e-09 | 1.0002 | 1.1626 | - | - | 1.0e-02 | 1026.4 | 4.460e-09 |
| 6 | 4.720e-04 | 221 | 1.21527 | 0.1896 | 0.2271 | 0 | 4.04e-09 | 1.0002 | 1.1626 | - | - | 1.0e-02 | 863.2 | 4.650e-09 |
| 7 | 1.920e-04 | 223 | 1.22296 | 0.4068 | 0.2502 | 0 | 7.77e-09 | 1.0002 | 1.1626 | - | - | 1.0e-02 | 859.1 | 3.200e-09 |
| 8 | 5.460e-05 | 234 | 1.27036 | 0.2844 | 0.2548 | 0 | 6.88e-09 | 1.0002 | 1.1626 | - | - | 1.0e-02 | 834.2 | 4.280e-09 |
| 9 | 2.230e-05 | 246 | 1.33357 | 0.4084 | 0.2703 | 0 | 6.79e-09 | 1.0002 | 1.1626 | - | - | 1.0e-02 | 335.4 | 2.960e-09 |
| 10 | 2.160e-05 | 221 | 1.21527 | 0.9686 | 0.3115 | 0 | 3.26e-09 | 1.0002 | 1.1626 | - | - | 1.0e-02 | 233.2 | 3.600e-09 |
| 11 | 1.420e-05 | 221 | 1.21527 | 0.6574 | 0.3357 | 0 | 6.38e-09 | 1.0002 | 1.1626 | - | - | 1.0e-02 | 231.8 | 2.870e-09 |
| 12 | 8.180e-06 | 221 | 1.21527 | 0.5761 | 0.3526 | 0 | 3.25e-09 | 1.0002 | 1.1626 | - | - | 1.0e-02 | 194.4 | not measured |

The two hosts give the same trajectory pass for pass. MEASURED: the three
passes `i4_wm0083` reached on `lart3`, the one of `i4_wm055` and the five of
`i4_kz0083` carry the same worst gated row to every digit printed as the
same passes on `lart4`. The wall clock per pass fell by a factor 8.545
(`i4_wm0083`, 1161.18 s to 135.89 s), 6.425 (`i4_wm055`, 4775.23 s to
743.17 s) and 8.558 (`i4_kz0083`, 1820.51 s to 212.72 s) on the first pass
of each. Total wall clock on `lart4`: 18774 s over 37 passes
(`i4_wm0083`, 507 s a pass), 18135 s over 16 (`i4_wm055`, 1133 s a pass),
9419 s over 40 (`i4_kz0083`, 235 s a pass), and 11101 s over 12 on `lart3`
for `i3_alt` (925 s a pass).

How each run ended (MEASURED from the log tails): `i4_kz0083` reached the
pass cap of 40, returned `info = 1` and wrote its state as a relaxation
snapshot with no stationary claim. `i4_wm0083` and `i4_wm055` were stopped
by the advisor at 10:10 KST inside the hydrodynamic solve of their next
pass (38 and 17), and in both that solve had stopped falling: `i4_wm0083`
at JFNK iteration 183 with `||R||` at 2.061e-02, 2.076e-02, 2.079e-02 over
its last three iterations and all 500 cells outside the mass tolerance;
`i4_wm055` at iteration 255 with `||R||` 1.110e-01 to 1.494e-01 and all 500
cells outside the mass tolerance.

---

## 4. The row's composition: where it is and where it is not

The term decomposition the item asks for, the chemistry against the two
advective face fluxes, is NOT in any of the five logs. None of the I3 or I4
runs set `EXHALE_CARRIER_ROW_TERMS=1` (READ, `LHS1140b/models/.L22/launch.sh`
and `launch_i4_lart4.sh`, whose environment is `OMP_NUM_THREADS`,
`OMP_WAIT_POLICY`, `EXHALE_PTC_DTAU0`, `EXHALE_OUTER_PASSES` and
`EXHALE_JUDGED_ROWS`), and no `carrier_row_terms.txt` exists anywhere under
`LHS1140b/models/.L22/` (MEASURED, `find`). What the logs do carry per pass
is the aggregate: the carrier residual of the returned composition as
`|res|` over the row's own terms, the volume-weighted and worst-cell carrier
steady residual, and the accepted composition displacement. Those are in
the tables above.

The decomposition exists for the three FROZEN states these runs started
from, and is READ from `docs/lhs1140b_stationary_L22_20260916.md` step 2b
section 2, whose dump is the certification of each frozen state (the state
`wellmixed/HeH0.083` carries is the one its campaign left at its pass 24;
rates are volumetric, cm^-3 s^-1, and the row is
`residual = (F_dif,in + F_dif,out) + (F_adv,in + F_adv,out) - net source`):

| case, cell | r [R_p] | T [K] | x_c | F_dif,in | F_dif,out | F_adv,in | F_adv,out | net source | residual | measure |
|---|---|---|---|---|---|---|---|---|---|---|
| `wellmixed/HeH0.083`, 500 | 29.031 | 285.3 | 7.807e-03 | +9.7857e-04 | +0.0 | -1.3366e-02 | +1.3042e-02 | -1.9574e-05 | +6.7473e-04 | 2.4619e-02 |
| `kzz1e9/HeH0.083`, 500 | 29.031 | 216.5 | 3.271e-02 | +3.4783e-03 | +0.0 | -7.7813e-02 | +7.4743e-02 | -5.5024e-05 | +4.6382e-04 | 2.9715e-03 |
| `wellmixed/HeH0.55`, 306 | 1.952 | 4461.5 | 2.664e-06 | -9.5084e-02 | +7.2185e-02 | -1.8544e-02 | +1.6801e-02 | -1.9149e-02 | -5.4924e-03 | 2.4767e-02 |

READ from the same section: in the two cell-500 rows the outer diffusive
flux is exactly zero, because `carrier_face_coefficients` and
`carrier_residual` fill the face arrays over `j = 1, N-1` only, so the
carrier column is a closed box for the diffusive, eddy and settling flux at
both ends; the chemistry contributes 2.9 and 11.9 per cent of the residual
there, and the two advective face fluxes very nearly cancel. Also READ from
it: at cell 500 of `wellmixed/HeH0.083` the whole outer column sits at the
refusing level, cells 339 to 344 reading 2.01e-02 to 2.06e-02 against
2.4619e-02 at cell 500, while at cell 500 of `kzz1e9/HeH0.083` the ratio to
cell 499 is 8.44 and the row IS a single-cell spike.

Those are properties of the entry states, one to forty passes before the
states this item classifies. Whether they still hold at pass 37 and pass 40
was NOT established here and is the first item of section 6.

---

## 5. The progress factors, confirmed and corrected

MEASURED by the reader on the three `run.log` files. The revised handoff's
three figures are confirmed:

| claim (READ, `docs/session_handoff_20260917_rev1.md` section 4.1) | MEASURED here | verdict |
|---|---|---|
| `i4_wm0083` 0.953 a pass over passes 26 to 37 (8.16e-02 to 4.79e-02) | 0.95273, that is 4.73 per cent a pass | confirmed |
| about 175 further passes at that constant factor | 175.0 | confirmed |
| `i4_wm055` 0.912 a pass over passes 13 to 16 (3.74e-02 to 2.84e-02) | 0.91232, that is 8.77 per cent a pass | confirmed |
| about 86 further passes at that constant factor | 86.7 | confirmed |
| `i4_kz0083` about 1 per cent a pass | 1.14 per cent over passes 20 to 40, 1.04 per cent over 23 to 40, 0.88 per cent over 30 to 40; 384, 420 and 497 further passes respectively | confirmed |

Three things the single-window figures leave out, and they change the
reading:

1. **`i4_wm0083` rose before it fell.** Over passes 1 to 22 the worst gated
   row went from 2.65e-02 to 9.20e-02, a factor 1.0611 a pass. The 0.953 of
   the handoff is the second half of a curve with a maximum at pass 22, not
   a rate the run has held. Any extrapolation from it is conditional in that
   second sense as well.
2. **The rate follows the movement bound.** Over the passes in which the
   bound was constant (MEASURED): `i4_kz0083` fell 2.49 per cent a pass over
   passes 10 to 21 at a bound of 2.5e-03 and 1.04 per cent a pass over
   passes 23 to 40 at a bound of 1.25e-03. Halving the bound took the
   decrement to 0.42 of its value, against 0.50 for the bound itself.
   `i4_wm0083` fell 3.22 per cent a pass over passes 23 to 36 at a bound of
   5.0e-03. This is the same proportionality MODELS.md already records for
   the campaign runs, "the bound was cut to its floor 1.0e-03 with the
   movement of a pass falling in proportion" (READ, `LHS1140b/MODELS.md`
   line 106).
3. **The bound is near its floor and the under-relaxation is at its floor.**
   `carrier_trust_floor = 1.0d-3` and the halving is
   `trust_pass = max(0.5d0*trust_pass, carrier_trust_floor)`, fired on a
   pass that counted as no progress (READ, `src/EXHALE_main.f90` lines 6619,
   7074 to 7075); `comp_omega = max(0.5d0*comp_omega, 0.125d0)` on the same
   passes when element diffusion is on (READ, line 7058). MEASURED: the
   bound reached 2.5e-03 in `i4_wm0083` (2 halvings, at passes 22 and 37),
   2.5e-03 in `i4_wm055` (2, at passes 13 and 14) and 1.25e-03 in
   `i4_kz0083` (3, at passes 6, 8 and 22, one halving above the floor), and
   `comp_omega` reached its floor 0.125 in `i4_kz0083` at pass 8.

Two smaller corrections to the revised handoff's table, both MEASURED:

- For `i4_wm0083` it reads "`info = 0` and mass 1e-08 at passes 26 to 36".
  The hydrodynamic mass row at those passes is 4.50e-09 to 9.56e-08 with a
  median near 5e-08, every one of them inside the cell's own tolerance; and
  passes 26 to 36 are `info = 0`, which is right.
- For `i4_wm055` it reads "mass 5.72e-03 to 9.89e-02" at passes 13 to 16.
  The four values are 5.72e-03, 4.36e-01, 1.98e-01 and 9.89e-02, so the
  stated range understates the maximum by a factor 4.4.

---

## 6. The classification, case by case

### 6.1 `i4_wm0083`: front relaxation, refusal hosted by the last cell

The evidence (MEASURED):

- The refusing cell is 500 at every one of the 37 passes, and 500 is the
  outermost physical cell, r = 29.0312 R_p. A fixed index at the edge of the
  grid carries no information about whether anything is moving.
- What is moving is the molecular column. The x2 = 1e-2 contour entered the
  domain at pass 7 at 28.5489 R_p and walked inward monotonically to
  9.1312 R_p by pass 34, where it stopped moving over the last three passes.
  The x2 = 0.5 contour walked outward monotonically from 1.2565 to
  1.7719 R_p over the same 37 passes. The cell at which the movement bound
  was refused walked outward from 221 (1.2153 R_p) to 304 (1.9191 R_p).
- The row's maximum at pass 22 coincides with the x2 = 1e-2 contour passing
  about 11 R_p, and the row falls from there at 3.22 per cent a pass at a
  fixed bound.

So the refusal at cell 500 is the outer end of a tail that is retreating
through the grid, and the classification is front relaxation. The rate is
not the transport time alone: it is throttled by the movement bound, which
is the plan's own second hypothesis met in a form neither of its two names
covers.

What would distinguish it further, if the two readings needed separating:
one pass reloaded from the state the run left, with
`EXHALE_CARRIER_ROW_TERMS=1`, gives the term decomposition at cell 500 at
pass 37 against the entry values of section 4. If the diffusive influx is
still the largest term and the outer diffusive face is still exactly zero,
the row is the closed outer box and what is decaying is the tail's amplitude
at that box; if the advective pair has taken over, the tail has left and the
row is a transient of the advection alone.

### 6.2 `i4_wm055`: front relaxation on a wind that is not a solution

The evidence (MEASURED):

- The refusing cell walks outward monotonically over the sixteen passes:
  315, 326, 333, 341, 350, 359, 369, 379, 389, 400, 411, 423, 435, 441,
  444, 445, that is 2.1136 to 11.7471 R_p. The revised handoff's list of
  the passes 5 to 16 is reproduced exactly.
- The row rises to 9.01e-02 at pass 5 and then decays monotonically to
  2.84e-02 at pass 16, a factor 0.9522 a pass over the whole run.
- The x2 = 1e-2 contour walks outward from 1.5539 to 3.8071 R_p and the
  x2 = 0.5 contour from 1.1626 to 1.2191 R_p.
- The advance of the refusing cell per pass falls with the movement bound:
  in radius, 8.32 to 10.03 R_p between passes 12 and 13, then 0.99, 0.54
  and 0.19 R_p over the next three, while the bound went 1.0e-02, 5.0e-03,
  2.5e-03, 2.5e-03.

A moving cell with a decaying row is front relaxation by the plan's own
criterion, and this case meets it without qualification. The qualification
is elsewhere: **the wind under the front has stopped converging.** The
hydrodynamic solve returned `info = 2` at passes 7, 10, 13, 14, 15 and 16,
with the certified mass row at 1.04e-01, 3.05e-03, 5.72e-03, 4.36e-01,
1.98e-01 and 9.89e-02 against a cell tolerance of order 3e-08, that is six
to seven decades outside. The coupled residual of the whole state rose over
the run from 5.92e+03 to 2.62e+06. The carrier row is falling on a sequence
of states that are moving away from a solution, so the extrapolation of
86.7 passes is not a forecast of certification; it is an extrapolation of
one row.

### 6.3 `i4_kz0083`: neither, and why

Reading it as a mode is what the cell index alone supports. The refusing
cell is 500 at passes 1 to 3 and 5 to 8 and 499 at pass 4 and at every pass
from 9 to 40; the x2 = 0.5 contour stands at 1.3705 or 1.3770 R_p for all
40 passes, and the x2 = 1e-2 contour never enters the domain at all. The
carrier residual of the returned composition is 1.72e-02 at pass 1 and
1.75e-02 at pass 39, flat to two digits over 39 passes: the alternation
never gets the composition near the fixed point of the carrier operator.

Reading it as a mode is refused by the rate. The row is not bounded: it
falls monotonically from pass 9 (1.37e-03) to pass 40 (8.19e-04), it fell
2.49 per cent a pass while the bound stood at 2.5e-03 and 1.04 per cent a
pass after the bound halved to 1.25e-03, and the ratio of those two
decrements is 0.42 against 0.50 for the bound. A mode the alternation cannot
damp would not care what the bound is.

So the honest statement is the third one: the alternation IS damping this
row, at a rate the composition movement bound sets, and the bound has been
cut to within one halving of its floor by three isolated passes that did not
count as progress. Nothing in the state relaxes, because the tail extends
past the outer boundary and the outer diffusive face carries nothing, so the
row at cells 499 and 500 is what it was at the start, walked down slowly.

What would distinguish a bound-throttled relaxation from a mode, and it is
cheap: reload the state the run left and take ten passes with the movement
bound held at 1.0e-02 and not reduced. The proportionality above predicts
about 8 per cent a pass. If the row falls at that rate the throttle is the
whole story; if it falls at 1 per cent whatever the bound, the row is a mode
and L32's verdict is the next item for it. The movement bound is a runtime
quantity of the outer loop and holding it fixed needs a key, so this is a
source change and NOT part of this item.

**Executed 2026-09-18** (`docs/lhs1140b_stationary_D8bound_20260918.md`,
section 5): the prediction above did not hold. From this run's exit state
with the bound held at 1.0e-2 the row went 8.2e-4, 6.3e-3, 1.9e-3 over three
passes and the wind solve returned `info = 2`; the monotone 1 per cent a
pass descent is a property of the small bound, not a slowed version of a
faster descent. Neither reading of this section (throttle or mode) is
supported, and no item was opened on the movement bound.

L32's verdict does not apply to any of the three as they stand: the coupled
block was never entered in any of them (section 7 below), so the linear
solve L32 measures was never the thing that failed here.

---

## 7. The handover was not exercised, and why

`Coupled carrier solve: On stall` was armed in all three I4 runs (READ,
their `input.inp`). The handover fires when the outer loop has ended on
`outer_no_progress`, which asks for `outer_no_fall_max = 3` CONSECUTIVE
passes that did not count as progress, AND `n_bound_endings >= 3`, three
consecutive carrier relaxations that ended on the movement bound (READ,
`src/EXHALE_main.f90` lines 6613, 7086, 7641 to 7646). MEASURED:

| run | passes | carrier relaxations ending on the movement bound | passes that did not count as progress | longest consecutive run of them | handover |
|---|---|---|---|---|---|
| `i4_wm0083` | 37 | 37 of 37 | 22, 37 | 1 | never fired |
| `i4_wm055` | 16 | 16 of 16 | 13, 14 | 2 | never fired |
| `i4_kz0083` | 40 | 39 of 39 that took an update | 6, 8, 22 | 1 | never fired |

This corrects the wording of the revised handoff, which says "neither of the
handover's two counts reached three". The bound-ending count reached three
in all three runs and stood at 37, 16 and 39. The count that never reached
three is the consecutive-no-progress one, and the closest any run came is
`i4_wm055` at two, one pass short. The handover therefore remains an
**unexercised route and nothing more**: no statement, favorable or not,
about the criterion or about what the block would have done follows from
these runs.

---

## 8. The one bounded continuation experiment, designed and not run

The plan allows one experiment on `lart4` for the cases classified as front
relaxation, sections 6.1 and 6.2. This section is the design; nothing was
run.

### 8.1 Which quantity sets the front's position

The plan offers `K_zz` or the base H2 fraction. The logs and the case files
decide it, and the answer is the base H2 fraction.

`K_zz` in these case names is NOT the carrier's eddy coefficient. READ,
`LHS1140b/MODELS.md` lines 68 and 69: `wellmixed` means no element
diffusion, so He/H is uniform, and `kzz1e9` means binary H/He ELEMENT
diffusion with `He_Kzz: 1.0e9` cm^2/s. MEASURED by diffing the four case
directories: `i4_wm0083/input.inp` and `i4_kz0083/input.inp` differ in the
planet name and in exactly two keys, `He_diffusion: True` and
`He_Kzz: 1.0e9`, and their `base.inp` files are identical. The eddy
coefficient the carrier transport reads, `Kzz_base`, is 1.0e9 in all four
cases including the certified one. So `K_zz` cannot be what separates the
certified state from the three refusing ones, and MEASURED at the same base
H2 fraction 0.8576 the element diffusion moves the x2 = 0.5 front only from
1.2565 to 1.3705 R_p, with the x2 = 1e-2 contour beyond the grid in both.

The base H2 fraction moves the front by far more. READ from the four
`base.inp` files and MEASURED from the first pass of each run:

| case | He/H | `q_H2_base` | x2 = 0.5 [R_p] | x2 = 1e-2 [R_p] | outcome |
|---|---|---|---|---|---|
| `kzz1e9/HeH2.13` | 2.13 | 0.190110 | 1.0002 | 1.1626 | CERTIFIED at pass 12 |
| `wellmixed/HeH0.55` | 0.55 | 0.476179 | 1.1626 | 1.5539 | front relaxation, 16 passes |
| `wellmixed/HeH0.083` | 0.083 | 0.857606 | 1.2565 | beyond 30 R_p | front relaxation, 37 passes |
| `kzz1e9/HeH0.083` | 0.083 | 0.857606 | 1.3705 | beyond 30 R_p | bound-throttled, 40 passes |

The base fraction is tied to He/H in this family by the H2 fraction the
photochemical column carries at 1 microbar,
`q_H2_base = (f/2) / ((1 - f) + f/2 + He/H)` with
`f = 0.9999831600354949` (READ, the header of every `base.inp` of the
family), so a continuation in the base H2 fraction IS a continuation in
He/H, and every rung of it is a physical model of the same family rather
than an interpolation with no meaning.

### 8.2 The ladder

Seed: the state `i3_alt` certified,
`molecular_scalar_gj1132_kzz1e9/HeH2.13` at outer pass 12. It is
`LHS1140b/models/.L22/i3_alt/output/Hydro_ioniz_IC.txt` and
`Ion_species_IC.txt`, kept there by `pp_i3_alt.sh` beside the
post-processing rewrite; MEASURED, that file's `# coupling` line reads
`certified=T cert_reason=certified_in_wind recon=WENO3` with the provenance
of the 04:44:16 stationary run, while `Hydro_ioniz.txt` beside it is the
04:52:07 post-processing write and reads `certified=F`. It is the only certified
molecular state there is: `LHS1140b/MODELS.md` line 106 records that
`molecular_scalar_gj1132_wellmixed/HeH2.13` has no certified state to seed
from either. The ladder therefore runs in the `kzz1e9` family.

Rungs, uniform in `q_H2_base` because that is the quantity that sets the
front, with He/H recovered from it by the relation above (MEASURED
arithmetic):

| rung | `q_H2_base` | He/H | note |
|---|---|---|---|
| 0 | 0.190110258 | 2.1300 | the certified fiducial, re-solved to confirm the seed reproduces |
| 1 | 0.261627353 | 1.4111 | |
| 2 | 0.333144448 | 1.0008 | |
| 3 | 0.404661544 | 0.7356 | |
| 4 | 0.476178639 | 0.5500 | the catalog case `kzz1e9/HeH0.55` exactly |

Four steps of 0.0715 in `q_H2_base`. The size is chosen from the measured
response: the whole interval from the fiducial to He/H = 0.55 moves the
x2 = 1e-2 contour from 1.1626 to 1.5539 R_p, so a rung moves it by about
0.1 R_p, which is inside the front's own width on these states: READ from the
L22 memo step 2 section 3, the width of the H2 front on the frozen
`wellmixed/HeH0.55` state, measured as the cells over which x2 falls from
0.9 to 0.1 of its base-cell value, is 97 cells, cells 149 to 245, which is
r = 1.0604 to 1.3278 R_p on this grid (MEASURED). The one continuation the campaign already tried went
from He/H 0.083 to 0.083 across the element treatment in a single step and
did not descend in forty passes (READ, `LHS1140b/MODELS.md` line 106),
which is the reason for small rungs here.

Each rung: a fresh directory under `LHS1140b/models/.L33/rung<k>/`, the
fiducial's `input.inp` and `base.inp` copied with `He/H number ratio`,
`HeH_base` and `q_H2_base` changed and nothing else, `Load IC? True` with
the previous rung's certified state, `Restart intent: stationary`,
`Coupled carrier solve: False` (the route that certified; leaving it armed
would change what is being measured in the middle of a rung, and the
handover stays where section 7 leaves it), `EXHALE_PTC_DTAU0=1.0` (the
fiducial's own value, READ from its `REPRODUCE.md` and used by `i3_alt`),
`OMP_NUM_THREADS=8`, `EXHALE_JUDGED_ROWS=1`.

Host: `lart4`. Pass budget: `EXHALE_OUTER_PASSES=25` per rung, five rungs,
at most 125 passes. At the measured 235 s a pass of `i4_kz0083` on `lart4`
that is about 8 h if every rung uses its whole budget; the fiducial took 12
passes.

### 8.3 Acceptance

A rung is accepted only on the full certification the code prints,
`CERTIFIED: every active equation was evaluated and is within its tolerance`,
together with the admissibility record of the same block: 0 cells without a
chemical root, no active unvalidated physics, no unbudgeted accepted
correction, and the carrier relaxation of the accepted pass ending on its
own residual rather than on the movement bound. The worst gated carrier row
alone is not an acceptance, and no acceleration's own norm is ever one. A
rung that does not certify within its 25 passes ENDS the experiment; its
last state is the record and no later rung is seeded from a state that did
not certify.

### 8.4 What the outcome would mean

If every rung certifies up to He/H = 0.55, the front is reachable by
continuation and the refusal of the direct solves is a matter of the seed,
not of the equations; the same ladder extended to He/H = 0.083 is then the
next experiment and not this one. If a rung fails, the base H2 fraction at
which it fails is the measurement, and the state it leaves is the first
molecular state whose failure is bracketed between two solved neighbors,
which is what every following diagnostic has lacked.

An acceleration of the outer fixed-point map, Anderson mixing over the
composition update, is NOT part of this experiment. If the ladder stalls at
a rung, it becomes the thing to try on that one rung, behind a key that
defaults off, and any accelerated iterate is accepted on the full physical
residual and the admissibility checks alone.

The alternative continuation, in the carrier eddy coefficient `Kzz_base` at
fixed He/H = 0.55, reaches the catalog case with its own He/H untouched and
is the physically direct knob on the front. It is not proposed here for one
reason: there is no certified state at a reduced `Kzz_base` to start from,
so it would have to make its own seed first, and the plan asks for a
continuation from the certified fiducial.

---

## 9. Noticed beside the measurement, not acted on

1. **The advance of the refusing cell and the decay of the row are both set
   by the composition movement bound, not by a transport time.** Sections
   5.2, 6.1, 6.2 and 6.3 carry the numbers. The consequence for the plan is
   that "more passes" and "a removed movement limit" are not the only two
   options the item was written against: a bound that does not ratchet
   downward on isolated non-progress passes is a third, and the ratchet is
   what took `i4_kz0083` from 2.49 to 1.04 per cent a pass. The halving has
   no way back up in the loop as it stands (READ, `src/EXHALE_main.f90`
   lines 7074 to 7075: `trust_pass` only ever decreases within a run). A
   restoration on a pass that does make progress is a one-line change with a
   physical justification, and it is NOT made here because this item owns no
   source file.
2. **The wellmixed runs' hydrodynamic solve broke down before they were
   stopped.** Section 3 has the last iterations. That is a separate failure
   from the carrier row and it is not described anywhere in the revised
   handoff, which reports only the `info` code of the last completed pass.
3. **The certification's outer ghost is data and the relaxation's is a
   copy.** READ from the L22 memo step 2b section 1, which reports it as a
   defect of its own and does not fix it. It bears directly on every cell-500
   refusal in this memo: the row the certification refuses at cell 500 is
   1.12 and 2.7 times the row the relaxation's own trials balance there
   (READ, the same section). Until that is one ghost, a cell-500 refusal is
   partly a difference between two evaluations of the same row.
4. `LHS1140b/MODELS.md` line 106 states the movement-bound proportionality
   and the floor for the campaign runs. It is consistent with everything
   measured here and no correction to it was needed.

---

## 10. Where the record is

- The reader: `LHS1140b/models/.L22/outer_pass_history.py`. Run it as
  `python3 outer_pass_history.py <case directory> [--window FIRST LAST]`;
  a case directory gives both its logs and takes the cell radii from
  `output/Hydro_ioniz_IC.txt`. `--csv <file>` writes the same table as
  comma-separated values.
- The logs it read, unchanged: `LHS1140b/models/.L22/{i4_wm0083,i4_wm055,
  i4_kz0083}/{run.log,run_lart3.log}` and `i3_alt/run.log`.
- No source file, no case directory and no log was edited by this item, and
  no solve was run for it.
