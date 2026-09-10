# G23: the own-largest-term measure extended to the momentum and energy rows

**DECIDED AND APPLIED 2026-09-03**, as `Update_EXHALE_stage1.md` section 143. This is
the measurement the decision rests on. In the code the three scales are
`mass_flux_row_scale`, `momentum_row_scale` and `energy_row_scale` in
`steady_residual.f90`, with no other path; "G23" below is that variant and "A"
the measure it replaced.

P55 section 7(1) proposes giving rows 2 and 3 the measure P54b gave row 1, and
records the one thing it could not test: whether the JFNK can reach a state
that passes it. It can, on both planets. What it does not buy is the marching
departure the proposal was aimed at.

## What G23 is

```
   s_1(j) = max( |F_{j-1/2}| r^2 , |F_{j+1/2}| r^2 ) / dV_j      (= G, approved)
   s_2(j) = max( |dF_2(j)| , |S_2(j)| )                          (new)
   s_3(j) = max( |dF_3(j)| , |S_3(j)| , heat(j) , cool(j) )      (new)
```

Terms taken as `RK_rhs` and the ionization sweep produce them, cached at the
same point the face fluxes are, plus the operator-split transport sources where
those are active. **Acceptance measure only**: `cell_state_scales` is untouched,
so the Newton system and the line-search merit are those of variant A, for the
reason P54b measured (rescaling the solver by a quantity that varies by a factor
9 across the first two cells destroyed the solve; `s_2` varies by more).

Fallback: where a row contains nothing at all the legacy `D_k(|v|+c_s)/dr` is
used. No physical state reaches it.

Scratch build `censusG23/`; the cache rebuild counter reported per run
(`(G23) row-scale cache rebuilds:`) is 15 on the hot Uranus and 45 on
WASP-121b, out of order `1e5` scale requests, all from JFNK probe states. The
rebuild cannot recompute `heat`/`cool` from `u`, so those few requests carry a
neighbouring state's radiative terms; it is counted and reported rather than
hidden.

## (i) Does the Newton converge, and what does it cost

Same start, same procedure as the P54b A/B/G comparison. The hot Uranus is the
1 microbar `He/H = 0.0793` gate rung reloaded and Newton-finished; WASP-121b is
`backup/regression/wasp_full_newton` as that case defines it (cold start).

| | hot U **A** | hot U **G** | hot U **G23** | WASP **A** | WASP **G** | WASP **G23** |
|---|---|---|---|---|---|---|
| JFNK `info` | 0 | 0 | **0** | 0 | 0 | **0** |
| `\|\|R\|\|` at the start, in its own units | 2.805e-05 | 2.272e-02 | **4.134e-01** | 5.777e-04 | 1.087e-02 | **1.248e+00** |
| `\|\|R\|\|` accepted | 2.031e-06 | 8.315e-06 | 5.208e-06 | 7.398e-06 | 3.138e-06 | 8.369e-06 |
| iterations | 5 | 9 | **14** | 14 | 17 | **44** |
| residual evaluations | 115 | 206 | **321** | 338 | 406 | **1000** |
| flux spread at the gate | 2.106e-03 | 2.489e-04 | 2.501e-04 | 1.357e-03 | 5.627e-05 | 5.242e-05 |
| `log10 Mdot` | 10.34 | 10.34 | 10.34 | 13.21 | 13.21 | 13.21 |

**The answer to P55's open question is yes.** The state the new measure demands
exists and the existing solver reaches it: the hot Uranus starts at 0.41 on the
own-term measure -- the gate would have refused it by four orders of magnitude
-- and converges to 5.2e-6 in 14 iterations. WASP-121b starts at 1.25 and
converges to 8.4e-6 in 44.

The cost is real but not prohibitive: **1.6x** the residual evaluations of G on
the hot Uranus, **2.5x** on WASP-121b (2.8x and 3.0x the baseline A). The
mass-loss rate does not move on either planet, and the flux gate is not
degraded (2.501e-4 against G's 2.489e-4; 5.242e-5 against 5.627e-5).

## (ii) What the accepted roots look like row by row

Each root measured with the SAME instrument -- the G23 own-term scales -- so
the three columns are comparable. `R_k/s_k`, dimensionless.

**Hot Uranus**

| r | A: R1 / R2 / R3 | G: R1 / R2 / R3 | G23: R1 / R2 / R3 |
|---|---|---|---|
| 1.00019 (cell 1) | -2.15e-01 / 1.42e-01 / **-9.94e-01** | -2.16e-01 / 1.43e-01 / **-9.94e-01** | -2.16e-01 / 1.43e-01 / **-9.94e-01** |
| 1.01005 | 9.70e-04 / -1.54e-04 / -1.89e-01 | -2.47e-07 / -3.32e-08 / -8.20e-04 | **-3.13e-10 / -1.53e-08 / -3.89e-06** |
| 1.03016 | -6.68e-04 / -1.30e-05 / -5.53e-02 | 2.47e-06 / -2.50e-05 / 1.72e-03 | **4.67e-09 / 8.68e-09 / 4.70e-06** |
| 1.05012 | -2.70e-04 / -2.90e-04 / -3.88e-02 | -1.42e-05 / 9.41e-06 / 4.68e-03 | **-1.36e-09 / 5.00e-08 / -4.20e-06** |
| 1.09969 | -5.62e-04 / 2.92e-05 / 7.76e-02 | -2.88e-07 / -3.85e-07 / -1.13e-03 | **2.62e-10 / -1.31e-09 / -1.75e-05** |
| 1.20011 | -5.26e-05 / 4.17e-07 / 2.11e-02 | 5.28e-08 / 3.44e-08 / -2.11e-03 | **-3.15e-12 / -8.93e-12 / -7.63e-06** |

**WASP-121b**

| r | A: R1 / R2 / R3 | G: R1 / R2 / R3 | G23: R1 / R2 / R3 |
|---|---|---|---|
| 1.00020 (cell 1) | 6.86e-03 / 1.38e-03 / **8.03e-01** | -4.07e-03 / 1.87e-02 / **-7.75e-01** | -4.60e-03 / 2.02e-02 / **-8.78e-01** |
| 1.01000 | -5.95e-04 / 1.17e-05 / -4.43e-02 | -2.14e-08 / 3.92e-09 / -9.91e-06 | **8.04e-12 / -2.17e-12 / 3.02e-06** |
| 1.10007 | -1.01e-04 / -2.65e-04 / 2.81e-04 | -2.56e-08 / -3.13e-08 / 1.68e-04 | **-4.94e-12 / 3.66e-10 / 1.26e-05** |
| 1.19951 | -4.03e-05 / 3.87e-04 / 3.80e-03 | -8.49e-09 / -3.29e-07 / 1.71e-03 | **1.96e-12 / 3.88e-09 / 4.28e-04** |

Two readings, and they point in opposite directions.

* **Above the first cell G23 does exactly what it claims.** The hot Uranus's
  energy row at 1.05 goes from 4.68e-3 on the G root to 4.20e-6, a factor 1100;
  its momentum row at 1.03 from 2.50e-5 to 8.68e-9, a factor 2900. WASP-121b's
  energy row at 1.1995 goes from 1.71e-3 to 4.28e-4 and at 1.10 from 1.68e-4 to
  1.26e-5. P55 section 4's headline numbers -- 1.2e-3 and 6.0e-3 on the energy
  row, 2.5e-5 on the momentum row -- are removed.
* **Cell 1 does not move at all, on either planet, under any variant.** Its
  energy row is out by 99 percent on the hot Uranus and 78-88 percent on
  WASP-121b whichever measure is used, and the state passes only because the
  norm is volume weighted and that cell's volume is negligible. The measure
  change is not what stands between the solver and that cell.

## (iii) Does a root that passes it stay put under marching?

Each accepted root reloaded and marched 5000 steps with the JFNK hand-off and
the residual gate removed (`Solver:` and `Resid tol:` deleted,
`du_th [PLM,WENO3]: 1.0e9 1.0e-30`, `EXHALE_P54_TS=25`), **the same marching
binary in every case**, so only the starting state differs.

Peak-to-peak excursion of the face mass flux at `r = 1.03`, in units of the
root's own `F_0`:

| root | hot Uranus | vs A | WASP-121b | vs A |
|---|---|---|---|---|
| **A** (tree as it stands) | 1.127 | 1.00 | 4.277 | 1.00 |
| **G** (approved: mass row) | **0.509** | **0.45** | 4.291 | 1.00 |
| **G23** (rows 1, 2 and 3) | **0.495** | **0.44** | 4.292 | 1.00 |

* **The approved change already buys the reduction.** On the molecular hot
  Uranus, starting from the G root instead of the A root halves the departure
  (1.127 -> 0.509). P55's own baseline for this planet is 0.534 over its window,
  which the G number reproduces.
* **G23 adds 3 percent to that, and nothing at all on WASP-121b.** Driving the
  momentum row at 1.03 down by a factor 2900 and the energy row by 1100 moves
  the departure from 0.509 to 0.495. On WASP-121b all three roots depart
  identically to within 0.3 percent, on an excursion of 4.3.

**So the departure is not controlled by the part of the residual G23 removes.**
P55 section 4 predicts the hot Uranus's early departure rate from `R_2` at
`r = 1.03` to 9 percent; that prediction holds for the *first steps* of a run
from a root with that residual, but the excursion over 5000 steps does not
follow it down when the residual is removed. What is left is root-independent:
the same three WASP roots, with residuals spanning four orders of magnitude,
give the same 4.29. That is the signature of a kick applied at the restart
rather than a drift driven by the state -- P55 section 7(2), the composition
projection the reload injects -- and of cell 1, which no variant here moves.

## What the decision rests on, and what it does not

The case FOR G23, and it is the one that was taken: a row that is out by one
part in a thousand of the physics it balances is out, whether or not a
consequence has been demonstrated downstream, and G23 is the only measure here
under which that is visible at all. Above the first cell it drives every row of
both planets to 5e-12 to 4e-4 of its own terms, at 1.6x to 2.5x the residual
evaluations of the mass row alone and with no change to the rate or the flux
gate. The states the previous measure accepted are refused by four orders of
magnitude, which is the point.

What it is NOT: a fix for the marching departure. That is worth stating plainly
because the proposal (P55 section 7(1)) was aimed there. The mass row's scale
alone halves the hot Uranus's departure (1.127 -> 0.509); adding the momentum
and energy rows moves it to 0.495, and on WASP-121b three roots whose residuals
span four orders of magnitude depart identically (4.277 / 4.291 / 4.292). The
departure is root-independent, which is the signature of the restart kick of
P55 section 7(2) and of cell 1, neither of which any row scale reaches.

**The experiment that follows.** Remove the reload kick (re-equilibrate the
composition at the loaded `T` before the first hydro step of a `Load IC` run),
then re-measure the table above. If the departure becomes root-dependent, the
G23 root should hold where the others do not -- and that is the test this
measure was built to make possible. Run before the kick is removed, it measures
the kick.

## What was not established

* Whether G23 changes any regression case's verdict is measured over the named
  cases in `p54g23_census.md`; the two flagship cases above are the ones with
  the full row-by-row and marching measurement.
* The mechanism by which cell 1 resists every variant. It is item (V), and the
  own-term numbers here are the sharpest statement of it so far: its energy row
  is out by order unity on both planets in every root measured.
* Whether `s_3` should include the operator-split composition energy step. It
  does not (that stage is not in the residual at all), and P55 section 7(3)
  records that as an open question about the residual itself, not about the
  scale.
