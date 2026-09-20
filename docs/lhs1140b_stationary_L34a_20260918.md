# L34a: the archived atomic catalog re-evaluated, and the seven Roe cases solved

Item L34 of `docs/PLAN_20260917.md`, part (a), the ATOMIC part. Part (b), the
molecular reference solution and the molecular catalog, is another worker's and
is not in this memo. MEASURED 2026-09-18 (KST) on the host `lart4` unless a
number is marked READ.

| | |
|---|---|
| binary | `LHS1140b/models/EXHALE_3146d11b.x`, md5 `3146d11b4090306dcea75bb9718edd22`, manifest `LHS1140b/models/BINARY_MANIFEST_3146d11b4090.txt` (items L27, L28, L30, L30b, L31; HEAD `3c73905` plus uncommitted work), checked before every run |
| route, step 1 | `Restart intent: stationary evaluate` (`docs/lhs1140b_stationary_L9_20260916.md` section 3): the state the case carries is measured as it stands, one equilibrium sweep refreshes the composition at the loaded conserved variables, and the certification, the advection-corrected profiles, the mass-loss line and the transit spectrum are all made on that work state. No step and no solve |
| route, step 2 | `models/run_case.sh` as the catalog runs it, at `EXHALE_PTC_DTAU0=1.0e8`, `EXHALE_OUTER_PASSES=40`, `SEED_ATTEMPTS=1`, `FORCE=1`, 8 threads: the seed mapped by `src/utils/map_state_to_grid.py --ic`, the stationary solve, then the same evaluate route and the same products |
| threads | step 1 single-threaded (`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`), eight cases at a time; step 2 at 8 threads, four cases at a time |
| written in place | yes, as item L34 prescribes. The `output/` each case carried before this item is kept beside it in `output_pre_L34/`, together with that case's own previous `pp.log`, `tpm_He10830.txt` and `tpm_He10830_metrics.txt` |

## 1. Verdict

**Every one of the 79 archived atomic states reproduces its own stationary
claim under the new binary, and every one of the seven low-XUV Roe cases
solves and certifies.** No state was refused, so no state had to be re-solved
from its `_IC` copy, and the number of re-solves is zero. The seven Roe cases
were each ACCEPTED at outer pass 1 with `info = 0` and CERTIFIED.

Counts: 79 reproduced, 0 refused, 0 re-solved; 7 of 7 Roe cases certified.

**The states did not move, and the products moved only where the old
post-processing pass was wrong about them.** The evaluation is a measurement,
so the conserved variables it writes back are the ones it read: over the 79
states the mass density is bit-identical in every cell of every state, the
velocity moves in 32 of 79 states by at most 2.219e-16 relative, and the
temperature moves by a median 5.836e-14 and at most 3.807e-12, which is the
composition refresh of the work state and not a change of state. The mass-loss
rate `4 pi rho v r^2` of the outermost physical cell is unchanged to the last
bit in all 79, and the `Log10 of steady-state Mdot` line is unchanged to the
two decimals it carries in all 79. The He I 10830 red-pair equivalent width
moves by a median 5.4e-06 and at most 3.6e-05 relative, which is the size item
L9 measured for the change of post-processing route on the fiducial (5.1e-06,
READ) and comes from the `_adv` profiles, whose mass row the old pass took from
a PLM assembly of a WENO3 state.

**The element rows moved, and the certification held.** The gated elemental
He/H partition row of the 64 states that run element diffusion moved by a
median 0.9 per cent and at most 16.2 per cent between the value the case's own
solve accepted on `c2e9c9990b9f` and the value this evaluation measures on
`3146d11b`, and the cell that attains it is the same cell in 56 of those 64.
That movement is the size item L30 leads one to expect from one spherical
geometry for the transport, but it is NOT attributed to L30 alone here: the two
numbers come from two different binaries and from two different operations (an
outer pass of a solve against a no-step evaluation), and the predecessor
`59bfdb3fc4d0` is not on disk, so the before and after of that row on one
binary was not measured.

**The L27 face-flux budget is flat on every state.** Over the 79, the base face
mass flux stands within 3.5e-05 of the wind-window mean of `rho v r^2`, and the
largest face to face spread of `r_f^2 (rho v)_f` over a whole column is 5.05e-08
of that mean, with a median of 1.96e-09. The operator identity line reads
`WENO3, numerical flux: HLLC, well balanced: T` on all 79 and
`WENO3, numerical flux: ROE, well balanced: T` on the seven Roe cases, so the
budget was formed on the operator the states were solved with.

## 2. What was run, case by case

Nine of the 79 states are the last rung of a flux-closure ladder and their
`input.inp` was the input of the OLD post-processing pass, which
`src/utils/element_flux_closure.py` writes for that pass (`Do only PP: True`,
`CFL: 1.0e-12`, no `Solver` line, `Restart intent` dropped; line 692 of that
file, READ). The evaluate route refuses to start without `Solver`, and
`Do only PP` would skip the residual evaluation that is the whole measurement.
Each of the nine had its `input.inp` repaired to what its own solve had, READ
from that rung's `run.log` (`(input_read) Solver: Newton, JFNK hand-off at
du < 1.00E-02` and `(input_read) Restart intent: stationary`) and from the
rung's `closure.json` `input_keys` (`Well balanced: True`,
`Restart intent: stationary`, `Solver: Newton`, `Secondary_ionization:
Immediate`, `CFL: 1.0e-12`): `Solver: Newton` added back, `Do only PP` set to
`False`, `Restart intent: stationary` added back. `CFL: 1.0e-12` was left where
it is, because the solve carried it too and it is inert on a route that takes
no time step. Each of the nine `REPRODUCE.md` says so.

The nine rungs are the `k04` of `HeH3`, `HeH5`, `HeH7`, `HeH8`, `HeH9`,
`HeH9.7`, `HeH10` and `HeH12` and the `k05` of `HeH2.09`, all of
`atomic_photochem_gj1132_kzzprofile`. No other `input.inp` in the catalog was
changed by this item.

Which file carried the claim: `Hydro_ioniz.txt` reads `recon=PLM` in all 79
(MEASURED), so it is the relaxation snapshot the old pass left and the claim is
in the `_IC` copy, which is the file the evaluate route loads. This is the same
reading `models/reclassify_claims.sh` makes and it agrees with the L29 table in
every row.

## 3. The 79 archived atomic states

The columns: the verdict item L29 recorded for this state on
`59bfdb3fc4d0` (`models/CLAIMS_59bfdb3fc4d0.md`, READ); the verdict of the
evaluate route on `3146d11b`; the worst refusing entry, which is `--` where
nothing refuses; the gated elemental He/H partition row as the case's own solve
accepted it (READ from that case's `run.log`, on `c2e9c9990b9f`) against the
value this evaluation measures, `--` where the case runs no element diffusion;
the base face mass flux in units of the wind-window mean of `rho v r^2` and the
face to face spread in the same units, both from the L27 budget; whether the
`Log10 of steady-state Mdot` line is unchanged from the archived `pp.log`; the
He I 10830 red-pair equivalent width over 10832.60 to 10834.20 A; and its
relative change from the archived curve.

| state | L29 on 59bfdb3f | on 3146d11b | worst refusing entry | gated He/H row, the solve / this evaluation | base face | max-min face | log10 Mdot | He I 10830 EW [%A] | dEW |
|---|---|---|---|---|---|---|---|---|---|
| `atomic_photochem_gj1132_kzzprofile/HeH10/k04` | reproduced | reproduced | -- | 7.24E-06 at cell 217 / 7.242E-06 at cell 217 | 0.999982 | 1.55e-08 | same | 2.222381 | 3.6e-06 |
| `atomic_photochem_gj1132_kzzprofile/HeH12/k04` | reproduced | reproduced | -- | 5.98E-06 at cell 217 / 5.981E-06 at cell 217 | 0.999982 | 9.88e-09 | same | 2.228897 | 2.2e-06 |
| `atomic_photochem_gj1132_kzzprofile/HeH2.09/k05` | reproduced | reproduced | -- | 6.90E-06 at cell 217 / 6.908E-06 at cell 217 | 0.999984 | 1.70e-09 | same | 1.476507 | 5.4e-06 |
| `atomic_photochem_gj1132_kzzprofile/HeH3/k04` | reproduced | reproduced | -- | 6.03E-06 at cell 217 / 6.041E-06 at cell 217 | 0.999984 | 2.99e-09 | same | 1.764001 | 4.5e-06 |
| `atomic_photochem_gj1132_kzzprofile/HeH5/k04` | reproduced | reproduced | -- | 9.60E-06 at cell 217 / 9.601E-06 at cell 217 | 0.999984 | 2.28e-09 | same | 2.074571 | 7.2e-06 |
| `atomic_photochem_gj1132_kzzprofile/HeH7/k04` | reproduced | reproduced | -- | 6.58E-06 at cell 217 / 6.587E-06 at cell 217 | 0.999983 | 3.42e-09 | same | 2.176836 | 5.1e-06 |
| `atomic_photochem_gj1132_kzzprofile/HeH8/k04` | reproduced | reproduced | -- | 6.56E-06 at cell 217 / 6.565E-06 at cell 217 | 0.999983 | 5.38e-09 | same | 2.200367 | 4.5e-06 |
| `atomic_photochem_gj1132_kzzprofile/HeH9.7/k04` | reproduced | reproduced | -- | 7.36E-06 at cell 217 / 7.364E-06 at cell 217 | 0.999982 | 3.36e-09 | same | 2.220553 | 3.6e-06 |
| `atomic_photochem_gj1132_kzzprofile/HeH9/k04` | reproduced | reproduced | -- | 7.68E-06 at cell 217 / 7.678E-06 at cell 217 | 0.999983 | 4.32e-09 | same | 2.213695 | 4.1e-06 |
| `atomic_photochem_gj1132x0.10_kzzprofile/HeH9.7` | reproduced | reproduced | -- | 6.92E-06 at cell 369 / 6.197E-06 at cell 369 | 0.999971 | 4.71e-08 | same | 0.267366 | 0.0e+00 |
| `atomic_photochem_gj1132x0.30_kzzprofile/HeH9.7` | reproduced | reproduced | -- | 7.37E-06 at cell 217 / 7.375E-06 at cell 217 | 0.999974 | 1.74e-08 | same | 1.009657 | 6.9e-06 |
| `atomic_photochem_gj1132x0.33_kzzprofile/HeH9.7` | reproduced | reproduced | -- | 6.59E-06 at cell 217 / 6.592E-06 at cell 217 | 0.999974 | 1.42e-08 | same | 1.097886 | 7.3e-06 |
| `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13` | reproduced | reproduced | -- | 6.46E-06 at cell 217 / 6.469E-06 at cell 217 | 0.999984 | 1.42e-09 | same | 1.372110 | 5.8e-06 |
| `atomic_scalar_gj1132_kzz0/HeH0.55` | reproduced | reproduced | -- | 5.95E-06 at cell 288 / 5.576E-06 at cell 288 | 0.999965 | 2.16e-09 | same | 0.030495 | 3.3e-05 |
| `atomic_scalar_gj1132_kzz0/HeH2.6` | reproduced | reproduced | -- | 9.33E-06 at cell 302 / 9.204E-06 at cell 301 | 0.999983 | 1.18e-09 | same | 0.954900 | 6.3e-06 |
| `atomic_scalar_gj1132_kzz0/HeH3.0` | reproduced | reproduced | -- | 7.70E-06 at cell 302 / 7.591E-06 at cell 302 | 0.999983 | 3.87e-09 | same | 1.079782 | 4.6e-06 |
| `atomic_scalar_gj1132_kzz0/HeH3.5` | reproduced | reproduced | -- | 6.34E-06 at cell 302 / 6.247E-06 at cell 302 | 0.999984 | 1.61e-09 | same | 1.217480 | 5.7e-06 |
| `atomic_scalar_gj1132_kzz0/HeH3.7` | reproduced | reproduced | -- | 9.99E-06 at cell 301 / 9.901E-06 at cell 301 | 0.999984 | 2.41e-09 | same | 1.269788 | 5.5e-06 |
| `atomic_scalar_gj1132_kzz0/HeH3.9` | reproduced | reproduced | -- | 8.50E-06 at cell 302 / 8.416E-06 at cell 301 | 0.999984 | 5.02e-09 | same | 1.316829 | 5.3e-06 |
| `atomic_scalar_gj1132_kzz1e10/HeH0.55` | reproduced | reproduced | -- | 9.06E-06 at cell 290 / 8.904E-06 at cell 290 | 0.999980 | 9.44e-10 | same | 0.646344 | 4.6e-06 |
| `atomic_scalar_gj1132_kzz1e10/HeH1.06` | reproduced | reproduced | -- | 6.72E-06 at cell 217 / 6.728E-06 at cell 217 | 0.999983 | 1.29e-09 | same | 1.045564 | 4.8e-06 |
| `atomic_scalar_gj1132_kzz1e10/HeH1.15` | reproduced | reproduced | -- | 6.91E-06 at cell 217 / 6.916E-06 at cell 217 | 0.999983 | 8.67e-10 | same | 1.103160 | 5.4e-06 |
| `atomic_scalar_gj1132_kzz1e10/HeH1.29` | reproduced | reproduced | -- | 7.08E-06 at cell 217 / 7.093E-06 at cell 217 | 0.999984 | 1.24e-09 | same | 1.186659 | 5.1e-06 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.55` | reproduced | reproduced | -- | 8.59E-06 at cell 217 / 8.580E-06 at cell 217 | 0.999981 | 1.63e-09 | same | 0.806472 | 5.0e-06 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.795` | reproduced | reproduced | -- | 7.18E-06 at cell 217 / 7.174E-06 at cell 217 | 0.999983 | 1.52e-09 | same | 1.030203 | 4.9e-06 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.865` | reproduced | reproduced | -- | 7.06E-06 at cell 217 / 7.056E-06 at cell 217 | 0.999983 | 2.04e-09 | same | 1.086553 | 5.5e-06 |
| `atomic_scalar_gj1132_kzz1e11/HeH0.93` | reproduced | reproduced | -- | 8.82E-06 at cell 217 / 8.810E-06 at cell 217 | 0.999983 | 2.25e-09 | same | 1.136397 | 5.3e-06 |
| `atomic_scalar_gj1132_kzz1e5/HeH2.6` | reproduced | reproduced | -- | 8.75E-06 at cell 302 / 8.621E-06 at cell 301 | 0.999983 | 2.16e-09 | same | 0.963144 | 5.2e-06 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.0` | reproduced | reproduced | -- | 7.29E-06 at cell 302 / 7.178E-06 at cell 302 | 0.999983 | 1.95e-09 | same | 1.087986 | 5.5e-06 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.35` | reproduced | reproduced | -- | 5.99E-06 at cell 302 / 5.889E-06 at cell 302 | 0.999984 | 1.61e-09 | same | 1.186186 | 5.1e-06 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.64` | reproduced | reproduced | -- | 8.76E-06 at cell 301 / 8.665E-06 at cell 301 | 0.999984 | 3.19e-09 | same | 1.262743 | 5.5e-06 |
| `atomic_scalar_gj1132_kzz1e5/HeH3.93` | reproduced | reproduced | -- | 7.86E-06 at cell 301 / 7.772E-06 at cell 301 | 0.999984 | 1.33e-09 | same | 1.331422 | 5.3e-06 |
| `atomic_scalar_gj1132_kzz1e6/HeH0.55` | reproduced | reproduced | -- | 6.58E-06 at cell 290 / 6.222E-06 at cell 290 | 0.999967 | 1.29e-09 | same | 0.059108 | 0.0e+00 |
| `atomic_scalar_gj1132_kzz1e6/HeH2.4` | reproduced | reproduced | -- | 6.60E-06 at cell 301 / 6.475E-06 at cell 301 | 0.999983 | 1.43e-09 | same | 0.946510 | 5.3e-06 |
| `atomic_scalar_gj1132_kzz1e6/HeH2.8` | reproduced | reproduced | -- | 9.26E-06 at cell 301 / 9.149E-06 at cell 301 | 0.999983 | 3.33e-09 | same | 1.080379 | 5.6e-06 |
| `atomic_scalar_gj1132_kzz1e6/HeH3.19` | reproduced | reproduced | -- | 7.05E-06 at cell 301 / 6.950E-06 at cell 301 | 0.999984 | 4.62e-09 | same | 1.194855 | 5.0e-06 |
| `atomic_scalar_gj1132_kzz1e6/HeH3.46` | reproduced | reproduced | -- | 6.34E-06 at cell 301 / 6.246E-06 at cell 301 | 0.999984 | 1.06e-09 | same | 1.266941 | 5.5e-06 |
| `atomic_scalar_gj1132_kzz1e6/HeH3.74` | reproduced | reproduced | -- | 9.09E-06 at cell 300 / 9.006E-06 at cell 300 | 0.999984 | 2.82e-09 | same | 1.338791 | 6.0e-06 |
| `atomic_scalar_gj1132_kzz1e7/HeH0.55` | reproduced | reproduced | -- | 9.66E-06 at cell 289 / 9.357E-06 at cell 288 | 0.999972 | 1.74e-09 | same | 0.158538 | 6.3e-06 |
| `atomic_scalar_gj1132_kzz1e7/HeH2.70` | reproduced | reproduced | -- | 8.45E-06 at cell 292 / 8.366E-06 at cell 292 | 0.999984 | 1.76e-09 | same | 1.199631 | 5.0e-06 |
| `atomic_scalar_gj1132_kzz1e7/HeH2.94` | reproduced | reproduced | -- | 7.04E-06 at cell 217 / 7.049E-06 at cell 217 | 0.999984 | 9.01e-10 | same | 1.272364 | 5.5e-06 |
| `atomic_scalar_gj1132_kzz1e7/HeH3.18` | reproduced | reproduced | -- | 9.41E-06 at cell 302 / 9.498E-06 at cell 302 | 0.999984 | 7.26e-09 | same | 1.326511 | 6.0e-06 |
| `atomic_scalar_gj1132_kzz1e8/HeH0.55` | reproduced | reproduced | -- | 7.14E-06 at cell 292 / 6.893E-06 at cell 292 | 0.999975 | 1.50e-09 | same | 0.308908 | 3.2e-06 |
| `atomic_scalar_gj1132_kzz1e8/HeH2.05` | reproduced | reproduced | -- | 9.29E-06 at cell 297 / 9.190E-06 at cell 297 | 0.999984 | 1.96e-09 | same | 1.151289 | 5.2e-06 |
| `atomic_scalar_gj1132_kzz1e8/HeH2.23` | reproduced | reproduced | -- | 8.86E-06 at cell 297 / 8.769E-06 at cell 297 | 0.999984 | 2.21e-09 | same | 1.220628 | 5.7e-06 |
| `atomic_scalar_gj1132_kzz1e8/HeH2.41` | reproduced | reproduced | -- | 7.10E-06 at cell 217 / 7.108E-06 at cell 217 | 0.999984 | 1.35e-09 | same | 1.282584 | 5.5e-06 |
| `atomic_scalar_gj1132_kzz1e9/HeH0.55` | reproduced | reproduced | -- | 8.18E-06 at cell 289 / 7.989E-06 at cell 289 | 0.999978 | 1.20e-09 | same | 0.479293 | 6.3e-06 |
| `atomic_scalar_gj1132_kzz1e9/HeH1.50` | reproduced | reproduced | -- | 6.44E-06 at cell 217 / 6.456E-06 at cell 217 | 0.999983 | 2.53e-09 | same | 1.099642 | 5.5e-06 |
| `atomic_scalar_gj1132_kzz1e9/HeH1.60` | reproduced | reproduced | -- | 6.89E-06 at cell 299 / 6.787E-06 at cell 299 | 0.999984 | 1.46e-09 | same | 1.151055 | 5.2e-06 |
| `atomic_scalar_gj1132_kzz1e9/HeH1.70` | reproduced | reproduced | -- | 6.92E-06 at cell 217 / 6.934E-06 at cell 217 | 0.999984 | 1.43e-09 | same | 1.196966 | 5.8e-06 |
| `atomic_scalar_gj1132_kzz1e9/HeH2.13` | reproduced | reproduced | -- | 6.67E-06 at cell 217 / 6.679E-06 at cell 217 | 0.999984 | 1.26e-09 | same | 1.374915 | 5.1e-06 |
| `atomic_scalar_gj1132_kzz1e9/HeH4.0` | reproduced | reproduced | -- | 5.00E-07 at cell 217 / 4.955E-07 at cell 217 | 0.999985 | 2.80e-09 | same | 1.822919 | 8.8e-06 |
| `atomic_scalar_gj1132_kzz1e9/HeH9.7` | reproduced | reproduced | -- | 5.27E-06 at cell 217 / 5.269E-06 at cell 217 | 0.999982 | 2.18e-09 | same | 2.059231 | 3.4e-06 |
| `atomic_scalar_gj1132_wellmixed/HeH0.083` | reproduced | reproduced | -- | -- / -- | 0.999970 | 7.18e-10 | same | 0.379690 | 5.3e-06 |
| `atomic_scalar_gj1132_wellmixed/HeH0.40` | reproduced | reproduced | -- | -- / -- | 0.999980 | 1.01e-09 | same | 1.101067 | 5.4e-06 |
| `atomic_scalar_gj1132_wellmixed/HeH0.42` | reproduced | reproduced | -- | -- / -- | 0.999981 | 1.23e-09 | same | 1.132365 | 6.2e-06 |
| `atomic_scalar_gj1132_wellmixed/HeH0.44` | reproduced | reproduced | -- | -- / -- | 0.999981 | 1.17e-09 | same | 1.162608 | 5.2e-06 |
| `atomic_scalar_gj1132_wellmixed/HeH0.55` | reproduced | reproduced | -- | -- / -- | 0.999982 | 8.52e-10 | same | 1.311761 | 1.1e-05 |
| `atomic_scalar_gj1132_wellmixed/HeH1` | reproduced | reproduced | -- | -- / -- | 0.999985 | 9.99e-10 | same | 1.695877 | 8.8e-06 |
| `atomic_scalar_gj1132_wellmixed/HeH10` | reproduced | reproduced | -- | -- / -- | 0.999982 | 1.06e-09 | same | 1.934720 | 0.0e+00 |
| `atomic_scalar_gj1132_wellmixed/HeH100` | reproduced | reproduced | -- | -- / -- | 0.999973 | 1.95e-09 | same | 1.964411 | 3.6e-05 |
| `atomic_scalar_gj1132_wellmixed/HeH1000` | reproduced | reproduced | -- | -- / -- | 0.999969 | 2.87e-09 | same | 1.919701 | 9.9e-06 |
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13` | reproduced | reproduced | -- | 8.38E-06 at cell 258 / 7.496E-06 at cell 257 | 0.999975 | 5.05e-08 | same | 0.000091 | 0.0e+00 |
| `atomic_scalar_gj1132x0.15_kzz1e9/HeH2.13` | reproduced | reproduced | -- | 8.71E-06 at cell 271 / 7.659E-06 at cell 270 | 0.999976 | 7.38e-09 | same | 0.000348 | 0.0e+00 |
| `atomic_scalar_gj1132x0.15_kzz1e9/HeH9.7` | reproduced | reproduced | -- | 7.88E-06 at cell 358 / 7.502E-06 at cell 358 | 0.999972 | 3.89e-08 | same | 0.451436 | 4.4e-06 |
| `atomic_scalar_gj1132x0.20_kzz1e9/HeH2.13` | reproduced | reproduced | -- | 9.81E-06 at cell 294 / 8.458E-06 at cell 293 | 0.999978 | 8.88e-09 | same | 0.002236 | 0.0e+00 |
| `atomic_scalar_gj1132x0.20_kzz1e9/HeH9.7` | reproduced | reproduced | -- | 5.49E-06 at cell 217 / 5.496E-06 at cell 217 | 0.999973 | 2.29e-08 | same | 0.629407 | 6.4e-06 |
| `atomic_scalar_gj1132x0.25_kzz1e9/HeH2.13` | reproduced | reproduced | -- | 6.56E-06 at cell 328 / 5.495E-06 at cell 328 | 0.999980 | 2.37e-08 | same | 0.069983 | 0.0e+00 |
| `atomic_scalar_gj1132x0.25_kzz1e9/HeH9.7` | reproduced | reproduced | -- | 9.15E-06 at cell 217 / 9.147E-06 at cell 217 | 0.999974 | 4.06e-08 | same | 0.791316 | 6.3e-06 |
| `atomic_scalar_gj1132x0.30_kzz1e9/HeH2.13` | reproduced | reproduced | -- | 9.36E-06 at cell 329 / 8.697E-06 at cell 329 | 0.999980 | 1.33e-08 | same | 0.189423 | 5.3e-06 |
| `atomic_scalar_gj1132x0.30_kzz1e9/HeH9.7` | reproduced | reproduced | -- | 7.55E-06 at cell 217 / 7.549E-06 at cell 217 | 0.999975 | 1.26e-08 | same | 0.939533 | 6.4e-06 |
| `atomic_scalar_gj1132x0.33_kzz1e9/HeH2.13` | reproduced | reproduced | -- | 6.02E-06 at cell 329 / 5.480E-06 at cell 328 | 0.999980 | 2.56e-09 | same | 0.255489 | 3.9e-06 |
| `atomic_scalar_gj1132x0.33_kzz1e9/HeH9.7` | reproduced | reproduced | -- | 6.48E-06 at cell 217 / 6.487E-06 at cell 217 | 0.999975 | 2.24e-08 | same | 1.023141 | 5.9e-06 |
| `atomic_scalar_gj699_wellmixed/HeH0.042` | reproduced | reproduced | -- | -- / -- | 1.000013 | 1.65e-10 | same | 0.909107 | 2.2e-05 |
| `atomic_scalar_gj699_wellmixed/HeH0.046` | reproduced | reproduced | -- | -- / -- | 1.000013 | 1.68e-10 | same | 0.985733 | 2.1e-05 |
| `atomic_scalar_gj699_wellmixed/HeH0.050` | reproduced | reproduced | -- | -- / -- | 1.000013 | 1.52e-09 | same | 1.060830 | 2.1e-05 |
| `atomic_scalar_gj699_wellmixed/HeH0.083` | reproduced | reproduced | -- | -- / -- | 1.000016 | 2.19e-10 | same | 1.627076 | 2.0e-05 |
| `atomic_scalar_gj699_wellmixed/HeH1` | reproduced | reproduced | -- | -- / -- | 1.000021 | 5.17e-10 | same | 4.524962 | 5.7e-06 |
| `atomic_scalar_gj699_wellmixed/HeH1000` | reproduced | reproduced | -- | -- / -- | 1.000011 | 1.43e-09 | same | 2.872491 | 4.9e-06 |

## 4. The seven low-XUV Roe cases

These seven carry `Numerical flux: ROE` in their catalog `input.inp` (the
decision of 2026-09-17) and had no `output/` state; their certified states lived
outside the catalog, in `models/.L25/`. Each was solved here from the state item
L25 left for it: three from the certified direct solves of L25 step 3 and four
from the last rung of the ladder its own `seed_from_ladder.txt` names (READ from
those files and from `docs/lhs1140b_stationary_L25_20260916.md` sections 3.3 and
4.5). For the three direct solves the certified `_IC` pair was copied to a
scratch directory and the mapper pointed at that copy, because
`map_state_to_grid.py` prefers `Hydro_ioniz.txt` where one exists and in those
directories that file is the old relaxation snapshot; nothing under
`models/.L25/` was written to. The seeds are already on the catalog cell centers
(MEASURED: the `# grid` line of each seed equals the one its target case is
solved on), and the mapper still ran, so each seed entered the solve stamped
`certified=F cert_reason=mapped_seed`, as it must.

The L25 numbers beside them: `Mdot` is `4 pi rho v r^2` of the outermost
physical cell with the state's own `R0`, which is the definition of L25 section
4.3, and the L25 value is READ from section 3.3 for the three direct solves and
from section 4.6 for the four ladder targets; the equivalent width is READ from
L25 section 3.3 for `x0.01_kzz1e9/HeH2.13` and from section 4.10 for the other
six.

| case | seed | pass | `info` | `||R||` | certification | log10 Mdot | Mdot [g/s], this run / L25 | He I 10830 EW [%A], this run / L25 | base face / window mean | max - min over the faces |
|---|---|---|---|---|---|---|---|---|---|---|
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | `models/.L25/roe_s001_HeH2.13`, the certified direct solve of L25 step 3 | 1 | 0 | 5.767E-07 | CERTIFIED | 5.79 | 6.1644e+05 / 6.1609e+05 (5.6e-04) | 1.3699e-06 / 1.37e-06 (6.1e-05) | 0.9999753 | 4.83e-08 |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | `models/.L25/ladder/L3/x0p01`, the last rung of ladder L3 | 1 | 0 | 8.546E-07 | CERTIFIED | 5.78 | 6.0017e+05 / 6.0016e+05 (1.9e-05) | 3.2459e-05 / 3.2456e-05 (8.3e-05) | 0.9999793 | 3.13e-08 |
| `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7` | `models/.L25/roe_s010_HeH9.7`, the certified direct solve of L25 step 3 | 1 | 0 | 4.362E-07 | CERTIFIED | 6.81 | 6.4138e+06 / 6.4102e+06 (5.6e-04) | 0.23532 / 0.23532 (1.5e-05) | 0.9999740 | 7.87e-08 |
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | `models/.L25/ladder/L4/x0p01`, the last rung of ladder L4 | 1 | 0 | 3.698E-07 | CERTIFIED | 5.77 | 5.8235e+05 / 5.8235e+05 (7.3e-06) | 3.3576e-05 / 3.3576e-05 (1.2e-05) | 0.9999710 | 2.25e-07 |
| `atomic_photochem_gj1132x0.15_kzzprofile/HeH9.7` | `models/.L25/ladder/L1/x0p15`, the last rung of ladder L1 | 1 | 0 | 6.210E-08 | CERTIFIED | 7.02 | 1.0543e+07 / 1.0543e+07 (2.9e-05) | 0.48218 / 0.48218 (8.1e-07) | 0.9999715 | 2.28e-08 |
| `atomic_photochem_gj1132x0.20_kzzprofile/HeH9.7` | `models/.L25/roe_p020_HeH9.7`, the certified direct solve of L25 step 3 | 1 | 0 | 6.628E-08 | CERTIFIED | 7.16 | 1.4347e+07 / 1.4347e+07 (1.4e-06) | 0.67111 / 0.67112 (8.7e-06) | 0.9999725 | 1.98e-08 |
| `atomic_photochem_gj1132x0.25_kzzprofile/HeH9.7` | `models/.L25/ladder/L2/x0p25`, the last rung of ladder L2 | 1 | 0 | 6.245E-08 | CERTIFIED | 7.26 | 1.8194e+07 / 1.8194e+07 (2.5e-05) | 0.84287 / 0.84287 (4.8e-06) | 0.9999736 | 2.36e-08 |

**Agreement with L25.** The largest disagreement in the mass-loss rate is
5.6e-04 relative and the largest in the equivalent width 8.3e-05, over states
re-solved on a binary that carries L27, L28, L30, L30b and L31 and from seeds
that were mapped, so this is the whole effect of those items plus one outer pass
on these seven columns. Both of the 0.01 cases keep an equivalent width of a few
times 1e-05 per cent A, so the statement of L25 section 4.10 that the line is
gone at 0.01 of the fiducial spectrum stands on the catalog states as well.

`atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` is the one case whose history
changes here. L25 step 3 solved it directly under Roe and was refused by the
gated elemental He/H partition row at 2.234e-04 against 1.0e-05 at cell 440
(READ, section 3.3); the L4 ladder then certified the same wind by continuation.
Seeded from that ladder state, the catalog case certifies at pass 1 on this
binary.

## 5. The bookkeeping

- `python3 models/status.py --write` was run from `LHS1140b/` after the
  86 runs of this item and rewrote sections 7 and 8 of `LHS1140b/MODELS.md`.
  The seven Roe rows changed from `none` to `info=0 certified` with their
  mass-loss rate, equivalent width, depth and FWHM filled in; two group
  crossings moved in their fourth decimal
  (`atomic_scalar_gj1132_kzz0` 3.1010 to 3.1011,
  `atomic_scalar_gj1132_kzz1e10` 1.1579 to 1.1580), three equivalent widths and
  two line depths in their fourth decimal, all of it the 5e-06 product movement of section 1.
  **The molecular rows of those sections are a snapshot of item L34b in
  flight** and will have to be regenerated when that item finishes.
- The `claim on 59bfdb3f` column of section 7 still points at
  `models/CLAIMS_59bfdb3fc4d0.md`, which is the L29 table on the previous
  binary and is a historical record. This item did not rewrite that table and
  did not write a `CLAIMS_3146d11b4090.md`; the verdicts of the current binary
  for the atomic states are the table of section 3 above.
- `models/write_reproduce.py` was called by the run that made each case, as the
  runners call it, so every one of the 86 `REPRODUCE.md` names
  `3146d11b4090306dcea75bb9718edd22`. It takes a case directory and was
  therefore run once for each case, not once for the tree.

## 6. What this item did not measure

- **The same rows on the predecessor binary.** `59bfdb3fc4d0` is not on disk
  (the tree's `EXHALE.x` is `3146d11b` and no other copy of the predecessor
  remains), and this item builds nothing, so the element-row movement of
  section 1 is quoted against the value each case's own solve accepted on
  `c2e9c9990b9f` and not against an evaluation of the same state on the
  binary L29 used.
- **The molecular catalog**, which is item L34b.
- **The memo figures and the typeset memo.**
  `LHS1140b/make_memo_figures.py` and `latexmk -pdf` for
  `docs/lhs1140b_exhale_vs_pwinds.tex` were NOT run, because the molecular
  states are not in yet.
- **Any comparison of the Roe cases against HLLC at matched resolution**,
  which is item L35.

## 7. Two defects of the record, found and not fixed here

Both are in `models/write_reproduce.py`, which this item runs but does not own.

1. **A run that took no seed is described as a cold start.** With an empty
   `--seed-line` the file writes "None: this case starts from the code's own
   initial condition (`Load IC? False`). The archive holds no state of its
   physics.", which is false for a case whose `input.inp` states `Load IC?
   True` and which was re-measured in place. The sentence was corrected in the
   79 files this item wrote; the tool would write it again.
2. **A seed named on the command line is attributed to `pick_seed.py`.**
   With `SEED=` set, `run_case.sh` passes `<dir> (given)` and the file writes
   "`models/pick_seed.py` chose, out of the states that carry this case's
   physics:" followed by the paragraph explaining `pick_seed`'s tier columns,
   none of which applies. This affects every case ever run with `SEED=`. The
   sentence and the paragraph were corrected in the seven files this item
   wrote; the tool would write them again.

A third point, not a defect: `models/run_case.sh` still describes its
post-processing pass in the comment at the head of the file as "That pass takes
one time step before it stops, which would move a certified state, so it runs at
CFL 1e-12", which is the route item L9 replaced. The code below it takes the
evaluate route and says so; only the file header is stale.
