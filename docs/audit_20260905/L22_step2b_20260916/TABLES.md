# L22 step 2b, 2026-09-16: the tables as the analysis produced them

## The frozen-state carrier row under each intervention (MEASURED)

All rates are volumetric [cm^-3 s^-1]. The row measure is the certification
one, the residual over the sum of the row own terms with the two faces of
each flux counted separately. The state is the same in every column
(largest relative difference of n(H2) over the 500 cells 0.0, MEASURED).

### wm083 cell 500
| intervention | worst gated row | worst excluding cell 500 | cell row measure | residual | F_dif,in | F_dif,out | F_adv,in | F_adv,out | net source |
|---|---|---|---|---|---|---|---|---|---|
| ctl | 2.4619e-02 at cell 500 | 2.1713e-02 at cell 499 | 2.4619e-02 | +6.7473e-04 | +9.7857e-04 | +0.0000e+00 | -1.3366e-02 | +1.3042e-02 | -1.9574e-05 |
| copy | 5.0629e-01 at cell 500 | 2.1713e-02 at cell 499 | 5.0629e-01 | -2.8797e-02 | +9.7857e-04 | -2.9471e-02 | -1.3366e-02 | +1.3042e-02 | -1.9574e-05 |
| extrap | 5.3605e-02 at cell 500 | 2.1713e-02 at cell 499 | 5.3605e-02 | +1.5141e-03 | +9.7857e-04 | +8.3938e-04 | -1.3366e-02 | +1.3042e-02 | -1.9574e-05 |
| upwind | 2.7614e-02 at cell 500 | 2.1713e-02 at cell 499 | 2.7614e-02 | +7.5912e-04 | +9.7857e-04 | +0.0000e+00 | -1.3366e-02 | +1.3127e-02 | -1.9574e-05 |
### kz083 cell 500
| intervention | worst gated row | worst excluding cell 500 | cell row measure | residual | F_dif,in | F_dif,out | F_adv,in | F_adv,out | net source |
|---|---|---|---|---|---|---|---|---|---|
| ctl | 2.9715e-03 at cell 500 | 3.5188e-04 at cell 499 | 2.9715e-03 | +4.6382e-04 | +3.4783e-03 | +0.0000e+00 | -7.7813e-02 | +7.4743e-02 | -5.5024e-05 |
| copy | 5.4027e-01 at cell 500 | 3.5188e-04 at cell 499 | 5.4027e-01 | -1.8398e-01 | +3.4783e-03 | -1.8445e-01 | -7.7813e-02 | +7.4743e-02 | -5.5024e-05 |
| extrap | 6.8795e-02 at cell 500 | 3.5188e-04 at cell 499 | 6.8795e-02 | +1.1497e-02 | +3.4783e-03 | +1.1033e-02 | -7.7813e-02 | +7.4743e-02 | -5.5024e-05 |
| upwind | 8.0880e-03 at cell 500 | 3.5188e-04 at cell 499 | 8.0880e-03 | +1.2690e-03 | +3.4783e-03 | +0.0000e+00 | -7.7813e-02 | +7.5548e-02 | -5.5024e-05 |
### wm055 cell 306
| intervention | worst gated row | worst excluding cell 500 | cell row measure | residual | F_dif,in | F_dif,out | F_adv,in | F_adv,out | net source |
|---|---|---|---|---|---|---|---|---|---|
| ctl | 2.4767e-02 at cell 306 | 2.4767e-02 at cell 306 | 2.4767e-02 | -5.4924e-03 | -9.5084e-02 | +7.2185e-02 | -1.8544e-02 | +1.6801e-02 | -1.9149e-02 |
| copy | 1.0894e-01 at cell 500 | 2.4767e-02 at cell 306 | 2.4767e-02 | -5.4924e-03 | -9.5084e-02 | +7.2185e-02 | -1.8544e-02 | +1.6801e-02 | -1.9149e-02 |
| upwind | 2.4767e-02 at cell 306 | 2.4767e-02 at cell 306 | 2.4767e-02 | -5.4924e-03 | -9.5084e-02 | +7.2185e-02 | -1.8544e-02 | +1.6801e-02 | -1.9149e-02 |

## The Jacobian action against a central difference of the full residual (MEASURED)

Frozen `molecular_scalar_gj1132_wellmixed/HeH0.55` state, direction
concentrated on cells 300 to 312, step 1e-06. Files:
`jacobian/jac_wm055_first_order.txt` and
`jacobian/jac_wm055_with_reconstruction.txt`.

| cell | J (first order) | J (with the banded reconstruction) | central difference | relative error, first order | relative error, with it |
|---|---|---|---|---|---|
| 299 | -2.669715e+00 | -2.651803e+00 | -2.651803e+00 | 6.709e-03 | 9.609e-11 |
| 300 | +5.824929e+00 | +5.813977e+00 | +5.813977e+00 | 1.880e-03 | 5.796e-11 |
| 301 | +2.051771e+00 | +2.041758e+00 | +2.041758e+00 | 4.880e-03 | 2.662e-10 |
| 306 | +8.422568e-01 | +8.386133e-01 | +8.425210e-01 | 3.136e-04 | 4.638e-03 |
| 312 | +9.443764e-01 | +9.392641e-01 | +9.411803e-01 | 3.384e-03 | 2.036e-03 |
| 313 | -4.621056e-01 | -4.624001e-01 | -4.606145e-01 | 3.227e-03 | 3.862e-03 |
| **314** | **+0.000000e+00** | **+0.000000e+00** | **+1.679903e-03** | **1.000** | **1.000** |

Cell 314 is the entry a block-tridiagonal matrix has no place for: the
direction stops at cell 312 and the reconstruction stencil of the face at
r_{313+1/2} reaches it, so the action needs d res(314)/d f_c(312). Over the
whole column the largest advective entry outside the band is 9.435e-02 of the
band of its own row (MEASURED, `EXHALE_L22B_JAC_RECON=1`).

Summary over the direction's support (cells 299 to 313, MEASURED): mean
relative error 1.542e-03 with the first-order entries and 3.463e-03 with the
banded ones; worst 6.709e-03 at cell 299 and 5.306e-03 at cell 302.

## The displacement a relaxation pass applies to the primitive state (MEASURED)

`molecular_scalar_gj1132_wellmixed/HeH0.55`, outer pass 1,
`EXHALE_L22B_DISPLACEMENT=1`.

| quantity | control | far-wind set waived |
|---|---|---|
| kept transport steps | 17 | 70 |
| carrier displacement of the pass | 2.904e-02 | 1.599e-01 |
| max abs(dp)/p | 8.862e-03 at cell 230, r 1.2521 | 5.016e-02 at cell 10, r 1.0019 |
| max abs(dT)/T | 4.971e-03 at cell 275, r 1.5539 | 1.870e-02 at cell 1, r 1.0002 |
| max abs(dmbar)/mbar | 9.804e-03 at cell 217, r 1.2007 | 6.448e-02 at cell 3, r 1.0006 |
| max abs(d(n_tot+n_e))/(n_tot+n_e) | 1.000e-02 at cell 221, r 1.2153 | 6.893e-02 at cell 3, r 1.0006 |
