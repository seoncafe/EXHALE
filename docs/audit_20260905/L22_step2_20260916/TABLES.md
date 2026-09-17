### `molecular_scalar_gj1132_wellmixed/HeH0.083`, refusing cell 500 (r = 29.0312 R_p)

Frozen state, the row the certification refuses: residual 6.747E-04, faces scale measure 2.462E-02, x2 7.807E-03, T 285.3 K.

| set / mode | omitted | interval-selection cells | accepted trials, grow min..max | rejections (reason x count) | cell that set or vetoed the bound | kept steps | drift | relaxation ending |
|---|---|---|---|---|---|---|---|---|
| full grid (control), pass 1 | 0 | 500 | 14, 1.485E-03..2.563E+01 | the movement bound x20 | 224 | 14 | 0.01693 | the movement bound, with the last step inside it |
| full grid (control), pass 2 | 0 | 500 | 16, 3.341E-03..2.563E+01 | the movement bound x21 | 226,227 | 16 | 0.01697 | the movement bound, with the last step inside it |
| narrow mask, veto, pass 1 | 131 | 369 | 14, 1.485E-03..2.563E+01 | the movement bound x20 | 224 | 14 | 0.01693 | the movement bound, with the last step inside it |
| narrow mask, veto, pass 2 | 131 | 369 | 16, 3.341E-03..2.563E+01 | the movement bound x21 | 226,227 | 16 | 0.01697 | the movement bound, with the last step inside it |
| narrow mask, waive, pass 1 | 131 | 369 | 16, 3.341E-03..3.844E+01 | the movement bound x21 | 290 | 16 | 0.03529 | the movement bound, with the last step inside it |
| narrow mask, waive, pass 2 | 131 | 369 | 17, 1.253E-03..3.844E+01 | the movement bound x21 | 290 | 17 | 0.03137 | the movement bound, with the last step inside it |
| far wind, veto, pass 1 | 388 | 112 | 14, 1.485E-03..2.563E+01 | the movement bound x20 | 224 | 14 | 0.01693 | the movement bound, with the last step inside it |
| far wind, veto, pass 2 | 388 | 112 | 16, 3.341E-03..2.563E+01 | the movement bound x21 | 226,227 | 16 | 0.01697 | the movement bound, with the last step inside it |
| far wind, waive, pass 1 | 388 | 112 | 26, 1.505E-03..2.919E+02 | the movement bound x27 | 389 | 26 | 0.2132 | the movement bound, with the last step inside it |
| far wind, waive, pass 2 | 388 | 112 | 24, 1.338E-03..6.568E+02 | the movement bound x26 | 389 | 24 | 0.179 | the movement bound, with the last step inside it |

| set / mode | refusing cell row measure, pass 1 entry | after the carrier update (carrier steady residual, volume / worst / cell) | refusing cell row measure after the coupled pass | after that pass's carrier update | hydrodynamic rows of pass 2 (mass / momentum / energy) | H2 front x2=0.5 [R_p], pass 1 -> 2 | safety stops |
|---|---|---|---|---|---|---|---|
| full grid (control) | 2.462E-02 | 7.640E-03 / 1.910E-01 / 1 | 3.345E-02 at cell 500 | 7.440E-03 / 1.830E-01 / 1 | 8.340E-09 / 6.790E-13 / 2.690E-08 | 1.2565 -> 1.2657 | 0 |
| narrow mask, veto | 2.462E-02 | 7.640E-03 / 1.910E-01 / 1 | 3.345E-02 at cell 500 | 7.440E-03 / 1.830E-01 / 1 | 8.340E-09 / 6.790E-13 / 2.690E-08 | 1.2565 -> 1.2657 | 0 |
| narrow mask, waive | 2.462E-02 | 7.350E-03 / 1.750E-01 / 1 | 4.537E-02 at cell 500 | 6.940E-03 / 1.610E-01 / 1 | 2.030E-08 / 2.230E-13 / 1.180E-08 | 1.2704 -> 1.2850 | 0 |
| far wind, veto | 2.462E-02 | 7.640E-03 / 1.910E-01 / 1 | 3.345E-02 at cell 500 | 7.440E-03 / 1.830E-01 / 1 | 8.340E-09 / 6.790E-13 / 2.690E-08 | 1.2565 -> 1.2657 | 0 |
| far wind, waive | 2.462E-02 | 5.240E-03 / 6.070E-02 / 1 | 2.362E-01 at cell 500 | 4.840E-03 / 4.150E-02 / 1 | 2.150E-08 / 1.200E-13 / 1.590E-08 | 1.3905 -> 1.6151 | 0 |

Column-wide change of the row dump against the frozen state (largest relative change over the 500 physical cells; the dump is the certification of the second pass where that pass was reached):

| set / mode | T | x2 | diffusive | advective | production | loss | net source | residual |
|---|---|---|---|---|---|---|---|---|
| full grid (control) | 1.61E-02 | 1.69E-01 | 6.72E+01 | 4.29E+00 | 1.32E-01 | 1.02E-01 | 8.31E-01 | 4.54E+01 |
| narrow mask, veto | 1.61E-02 | 1.69E-01 | 6.72E+01 | 4.29E+00 | 1.32E-01 | 1.02E-01 | 8.31E-01 | 4.54E+01 |
| narrow mask, waive | 3.38E-02 | 3.70E-01 | 1.30E+02 | 9.58E+00 | 2.56E-01 | 1.99E-01 | 1.67E+00 | 9.37E+01 |
| far wind, veto | 1.61E-02 | 1.69E-01 | 6.72E+01 | 4.29E+00 | 1.32E-01 | 1.02E-01 | 8.31E-01 | 4.54E+01 |
| far wind, waive | 2.04E-01 | 4.06E+00 | 4.04E+02 | 7.69E+01 | 8.35E-01 | 7.31E-01 | 7.36E+00 | 4.39E+02 |

### `molecular_scalar_gj1132_wellmixed/HeH0.55`, refusing cell 306 (r = 1.9517 R_p)

Frozen state, the row the certification refuses: residual -5.492E-03, faces scale measure -2.477E-02, x2 2.664E-06, T 4461.5 K.

| set / mode | omitted | interval-selection cells | accepted trials, grow min..max | rejections (reason x count) | cell that set or vetoed the bound | kept steps | drift | relaxation ending |
|---|---|---|---|---|---|---|---|---|
| full grid (control), pass 1 | 0 | 500 | 17, 9.766E-04..5.767E+01 | the movement bound x22 | 219,220,221 | 17 | 0.02904 | the movement bound, with the last step inside it |
| full grid (control), pass 2 | 0 | 500 | 18, 1.879E-03..5.767E+01 | the movement bound x22 | 226,227 | 18 | 0.02898 | the movement bound, with the last step inside it |
| narrow mask, veto, pass 1 | 99 | 401 | 17, 9.766E-04..5.767E+01 | the movement bound x22 | 219,220,221 | 17 | 0.02904 | the movement bound, with the last step inside it |
| narrow mask, veto, pass 2 | 99 | 401 | 18, 1.879E-03..5.767E+01 | the movement bound x22 | 226,227 | 18 | 0.02898 | the movement bound, with the last step inside it |
| narrow mask, waive, pass 1 | 99 | 401 | 18, 9.766E-04..8.650E+01 | the movement bound x22 | 170 | 18 | 0.04328 | the movement bound, with the last step inside it |
| narrow mask, waive, pass 2 | 99 | 401 | 18, 1.879E-03..5.767E+01 | the movement bound x22 | 278 | 18 | 0.03058 | the movement bound, with the last step inside it |
| far wind, veto, pass 1 | 388 | 112 | 17, 9.766E-04..5.767E+01 | the movement bound x22 | 219,220,221 | 17 | 0.02904 | the movement bound, with the last step inside it |
| far wind, veto, pass 2 | 388 | 112 | 18, 1.879E-03..5.767E+01 | the movement bound x22 | 226,227 | 18 | 0.02898 | the movement bound, with the last step inside it |
| far wind, waive, pass 1 | 388 | 112 | 70, 1.000E+00..1.414E+12 | none | (none refused) | 70 | 0.1599 | the carriers stopped moving at the fixed wind |

| set / mode | refusing cell row measure, pass 1 entry | after the carrier update (carrier steady residual, volume / worst / cell) | refusing cell row measure after the coupled pass | after that pass's carrier update | hydrodynamic rows of pass 2 (mass / momentum / energy) | H2 front x2=0.5 [R_p], pass 1 -> 2 | safety stops |
|---|---|---|---|---|---|---|---|
| full grid (control) | 2.477E-02 | 7.010E-03 / 2.620E-01 / 1 | 2.981E-02 at cell 297 | 6.460E-03 / 2.420E-01 / 1 | 1.150E-08 / 4.230E-13 / 2.280E-08 | 1.1626 -> 1.1684 | 0 |
| narrow mask, veto | 2.477E-02 | 7.010E-03 / 2.620E-01 / 1 | 2.981E-02 at cell 297 | 6.460E-03 / 2.420E-01 / 1 | 1.150E-08 / 4.230E-13 / 2.280E-08 | 1.1626 -> 1.1684 | 0 |
| narrow mask, waive | 2.477E-02 | 6.650E-03 / 2.430E-01 / 1 | 3.425E-02 at cell 300 | 6.140E-03 / 2.290E-01 / 1 | 1.850E-08 / 2.440E-13 / 1.200E-08 | 1.1654 -> 1.1714 | 0 |
| far wind, veto | 2.477E-02 | 7.010E-03 / 2.620E-01 / 1 | 2.981E-02 at cell 297 | 6.460E-03 / 2.420E-01 / 1 | 1.150E-08 / 4.230E-13 / 2.280E-08 | 1.1626 -> 1.1684 | 0 |
| far wind, waive | 2.477E-02 | 7.320E-10 / 1.050E-09 / 309 | not reached | not reached | not reached | 1.1542 | 0 |

Column-wide change of the row dump against the frozen state (largest relative change over the 500 physical cells; the dump is the certification of the second pass where that pass was reached):

| set / mode | T | x2 | diffusive | advective | production | loss | net source | residual |
|---|---|---|---|---|---|---|---|---|
| full grid (control) | 4.00E-02 | 3.41E+00 | 2.93E+01 | 3.89E+01 | 1.67E+00 | 2.34E+00 | 4.48E+00 | 4.11E+01 |
| narrow mask, veto | 4.00E-02 | 3.41E+00 | 2.93E+01 | 3.89E+01 | 1.67E+00 | 2.34E+00 | 4.48E+00 | 4.11E+01 |
| narrow mask, waive | 5.73E-02 | 6.43E+00 | 4.05E+01 | 7.19E+01 | 2.64E+00 | 3.91E+00 | 8.01E+00 | 1.15E+02 |
| far wind, veto | 4.00E-02 | 3.41E+00 | 2.93E+01 | 3.89E+01 | 1.67E+00 | 2.34E+00 | 4.48E+00 | 4.11E+01 |
| far wind, waive | 0.00E+00 | 0.00E+00 | 7.81E-08 | 2.58E-03 | 0.00E+00 | 0.00E+00 | 0.00E+00 | 9.56E-08 |

### `molecular_scalar_gj1132_kzz1e9/HeH0.083`, refusing cell 500 (r = 29.0312 R_p)

Frozen state, the row the certification refuses: residual 4.638E-04, faces scale measure 2.972E-03, x2 3.271E-02, T 216.5 K.

| set / mode | omitted | interval-selection cells | accepted trials, grow min..max | rejections (reason x count) | cell that set or vetoed the bound | kept steps | drift | relaxation ending |
|---|---|---|---|---|---|---|---|---|
| full grid (control), pass 1 | 0 | 500 | 34, 9.643E-03..7.482E+03 | the movement bound x31 | 51,53,54,55 | 34 | 0.02234 | the movement bound, with the last step inside it |
| full grid (control), pass 2 | 0 | 500 | 26, 1.505E-03..4.379E+02 | the movement bound x27 | 455,456,457 | 26 | 0.02623 | the movement bound, with the last step inside it |
| narrow mask, veto, pass 1 | 205 | 295 | 34, 9.643E-03..7.482E+03 | the movement bound x31 | 51,53,54,55 | 34 | 0.02234 | the movement bound, with the last step inside it |
| narrow mask, veto, pass 2 | 205 | 295 | 26, 1.505E-03..4.379E+02 | the movement bound x27 | 455,456,457 | 26 | 0.02623 | the movement bound, with the last step inside it |
| narrow mask, waive, pass 1 | 205 | 295 | 34, 1.205E-03..7.482E+03 | the movement bound x31 | 206 | 34 | 0.02277 | the movement bound, with the last step inside it |
| narrow mask, waive, pass 2 | 205 | 295 | 22, 9.514E-03..4.379E+02 | the movement bound x24 | 455,456,457 | 22 | 0.0262 | the movement bound, with the last step inside it |
| far wind, veto, pass 1 | 388 | 112 | 34, 9.643E-03..7.482E+03 | the movement bound x31 | 51,53,54,55 | 34 | 0.02234 | the movement bound, with the last step inside it |
| far wind, veto, pass 2 | 388 | 112 | 26, 1.505E-03..4.379E+02 | the movement bound x27 | 455,456,457 | 26 | 0.02623 | the movement bound, with the last step inside it |
| far wind, waive, pass 1 | 388 | 112 | 34, 2.411E-03..1.122E+04 | the movement bound x31 | 407,415,416,417,418 | 34 | 0.03551 | the movement bound, with the last step inside it |
| far wind, waive, pass 2 | 388 | 112 | 21, 9.766E-04..1.946E+02 | the movement bound x24 | 457,458 | 21 | 0.02575 | the movement bound, with the last step inside it |

| set / mode | refusing cell row measure, pass 1 entry | after the carrier update (carrier steady residual, volume / worst / cell) | refusing cell row measure after the coupled pass | after that pass's carrier update | hydrodynamic rows of pass 2 (mass / momentum / energy) | H2 front x2=0.5 [R_p], pass 1 -> 2 | safety stops |
|---|---|---|---|---|---|---|---|
| full grid (control) | 2.972E-03 | 2.130E-03 / 1.210E-02 / 1 | 1.391E-02 at cell 500 | 2.390E-03 / 1.310E-02 / 2 | 1.660E-01 / 1.210E-03 / 4.270E-02 | 1.3705 -> 1.3641 | 0 |
| narrow mask, veto | 2.972E-03 | 2.130E-03 / 1.210E-02 / 1 | 1.391E-02 at cell 500 | 2.390E-03 / 1.310E-02 / 2 | 1.660E-01 / 1.210E-03 / 4.270E-02 | 1.3705 -> 1.3641 | 0 |
| narrow mask, waive | 2.972E-03 | 2.110E-03 / 1.200E-02 / 1 | 1.542E-02 at cell 500 | 2.110E-03 / 1.450E-02 / 3 | 3.780E-01 / 2.870E-04 / 9.940E-01 | 1.3770 -> 1.3641 | 0 |
| far wind, veto | 2.972E-03 | 2.130E-03 / 1.210E-02 / 1 | 1.391E-02 at cell 500 | 2.390E-03 / 1.310E-02 / 2 | 1.660E-01 / 1.210E-03 / 4.270E-02 | 1.3705 -> 1.3641 | 0 |
| far wind, waive | 2.972E-03 | 1.860E-03 / 1.070E-02 / 1 | 2.174E-02 at cell 500 | 2.300E-03 / 1.370E-02 / 1 | 1.270E-08 / 3.980E-13 / 7.890E-09 | 1.3905 -> 1.3973 | 0 |

Column-wide change of the row dump against the frozen state (largest relative change over the 500 physical cells; the dump is the certification of the second pass where that pass was reached):

| set / mode | T | x2 | diffusive | advective | production | loss | net source | residual |
|---|---|---|---|---|---|---|---|---|
| full grid (control) | 2.81E-01 | 3.56E-02 | 1.23E+01 | 2.44E+04 | 6.27E-01 | 1.42E-01 | 1.04E+02 | 4.58E+02 |
| narrow mask, veto | 2.81E-01 | 3.56E-02 | 1.23E+01 | 2.44E+04 | 6.27E-01 | 1.42E-01 | 1.04E+02 | 4.58E+02 |
| narrow mask, waive | 2.13E-01 | 3.61E-02 | 1.41E+01 | 1.65E+04 | 5.48E-01 | 1.33E-01 | 1.12E+02 | 4.89E+02 |
| far wind, veto | 2.81E-01 | 3.56E-02 | 1.23E+01 | 2.44E+04 | 6.27E-01 | 1.42E-01 | 1.04E+02 | 4.58E+02 |
| far wind, waive | 3.04E-02 | 5.68E-02 | 2.29E+01 | 2.20E+00 | 1.99E-01 | 1.03E-01 | 1.67E+02 | 7.36E+02 |

