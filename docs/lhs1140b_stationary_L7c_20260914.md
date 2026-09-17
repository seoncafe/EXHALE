# LHS 1140 b, item L7c: T-L7-5 on the consolidated tree, and what the H2 carrier row at 1.2 R_p is made of

Item L7c of `docs/PLAN_20260913_lhs_stationary.md`: T-L7-5 repeated on the
consolidated tree (L4d, L4g, L6, L6b, L7, L7b, L10, L12-A) for both the
`local` and the handoff partition, and, since the H2 carrier row still
refuses, the term-by-term account of the cell it refuses at.

Binary: `EXHALE.x`, md5 `b261c3287e64149a85d0d2dfd2ee2c74`, used as delivered
(nothing was built and no source file was changed). Every number below is
MEASURED on that binary or on the files it wrote, except where it is marked
READ (from a source file) or RECONSTRUCTED (arithmetic done here on a
measured state, with the expression it uses named).

Scratch tree: `LHS1140b/models/.L7c/`. The atomic state the seeds are built
from was copied to `.L7c/src_atomic_HeH2.13/` before the runs, so the
campaign writing `LHS1140b/models/atomic_*` cannot move it underneath them.

## 1. Verdict

**Neither partition certifies, and the H2 carrier row that refuses is not a
statement about H2: the H2 it measures was put there by catastrophic
cancellation in `q_h2_equilibrium` (`src/modules/lower_atmosphere/lower_column.f90`
line 116).** The fit is evaluated as

```fortran
tenu = 10.0d0**u
qh2  = (1.9845d0 + tenu - sqrt(tenu*(3.9690d0 + tenu)))/2.3670d0
```

and the exact value of that difference is `A^2/((A + t) + sqrt(t(2A + t)))`
with `A` = 1.9845, i.e. about `A^2/(2 t)` once `t` is large. Formed as
written, the two leading terms cancel: the true difference falls below the
rounding of `t` itself at `t ~ A/sqrt(2 eps)` = 9.4e7, and beyond that the
result is the round-off of `t`, quantized in units of `ulp(t)`. MEASURED on
the LHS 1140 b column: the error passes 1 percent at cell 187 (r = 1.118,
T = 3311 K, `t` = 2.2e7), 186 of the 500 physical cells are wrong by more
than 1 percent, 166 of them receive exactly zero where the fit's value is
positive, and six receive more than twice the true value, the worst by a
factor **1.48e4** at cell 219.

That is the cell the certification refuses at. The `local` seed writes
x2 = 2 n(H2)/n_H = 5.044e-6 at cells 215, 216, 218 and 219 and exactly 0 at
214, 217 and 220, the same number in every one of the four -- the signature
of a difference that has become the rounding of `t` -- while the fit's own
value there is 5.8e-10 falling smoothly to 3.0e-10. At cell 218 the state
therefore carries **1.30e4 times** the H2 the fit allows there, and
**9.86e4 times** the H2 the code's own H + H <-> H2 equilibrium constant
allows (section 5.2); the carrier row measures that.

The fix is one line and is algebraically the same function (section 5): the
conjugate form has no cancellation anywhere, agrees with the present
expression to round-off wherever the present one is not cancelling, keeps the
cold limit 1.9845/2.3670 = 0.838403 exactly, and removes the need for the
`u > 30` branch, which today drops the value discontinuously to zero where
the function is 8.3e-31.

Three further things this item establishes:

- **The row's tolerance is not the obstruction.** By the code's own
  arithmetic bound the carrier row's rounding floor is at most
  `carrier_row_roundoff` = 64 eps = 1.42e-14 of the row's terms (READ,
  `diffusive_photochemistry.f90`), nine decades below the 1e-5 the wind gate
  asks. A cell holding 2.5e-6 of the hydrogen can satisfy 1e-5: the measure
  is relative to the row's own terms and carries its own absolute floor.
  Nothing arithmetic stands between this row and its tolerance.
- **The artifact cannot be removed by running the route longer.** MEASURED
  on the L7b pass table's own state files: over the four outer passes the H2
  mixing ratio at cell 218 moved by 3.9e-9 relative (5.044356915952e-06 ->
  5.044356935412e-06) while the H2 DENSITY fell 6 percent with the
  hydrogen. The carrier is advected with the nuclei and its mixing ratio is
  frozen, so a seeded spike survives every pass and the row reads 1 forever.
- **The consolidated tree does not descend on this route where the L7b tree
  did, and the pseudo-time ramp is why.** MEASURED on one state with one
  binary, the only difference the environment (section 4): with L4g's
  `pseudo_time_doubles_on_an_accepted_step` at its DEFAULT ON, `dtau` runs
  from 2 to its ceiling 1.0e14 in 124 iterations, the line search holds
  lam = 1 in 8 of them, and the residual norm never leaves 1.06 to 1.87; with
  `EXHALE_PTC_RAMP_DOUBLE=0` the same state at the same dtau0 keeps `dtau`
  between 0.09 and 408, holds lam = 1 in 149 of 159 iterations, and takes
  `||Fs||2` from 5.54 to 3.9e-4 and the residual norm to 4.5e-2. Neither
  completes an outer pass inside the hour, so the ramp is the difference
  between descending and not descending on this state; what is left between
  that and L7b's own speed (`||R||` 2.7e-3 by iteration 24) is not explained
  by this control and is a question for the ramp's own item.

**The molecular recipe is therefore NOT confirmed and
`run_case.sh.molecular.patch` is not updated to it.** Section 6 states what
blocks it and the smallest change, and reports one defect in the held patch
found while reading it.

## 2. The two seeds

Source: `LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH2.13/output`, whose
`Hydro_ioniz_IC.txt` carries `certified=T cert_reason=certified_in_wind`
(READ; the run log's closing verdict is `CERTIFIED ... certified IN THE WIND,
r >= 1.20E+00`). Copied to `.L7c/src_atomic_HeH2.13/` and read from there.

Target: `molecular_scalar_gj1132_kzz1e9/HeH2.13`, its own `input.inp` and
`base.inp` (`q_H2_base` = 0.19011025828621544 at He/H = 2.13, ceiling
0.19011406844106465), with `Load IC?` set to True for the conversion.
Invariant `p` for both, as the item asks.

MEASURED, the conversion (`seed.log` of each directory):

| | `local` | handoff |
|---|---|---|
| x2 requested | per cell | 9.9998316003549481e-01 |
| largest x2 actually transferred | 9.9999437572717320e-01 | 9.9998316003549481e-01 |
| cells clipped at the element-ratio ceiling | 234 | -- |
| cells capped by their own neutral hydrogen | 235 of 504 | 480 of 504 |
| max relative change of H nuclei / He nuclei / mass | 2.15e-16 / 0 / 2.22e-16 | 1.47e-16 / 0 / 2.53e-16 |
| largest relative move of T | 1.8169e-01 | 2.7580e-01 |
| equation-of-state closure, p round trip | 3.3647e-16 | 8.3748e-16 |
| `# coupling` of the written pair | `certified=F cert_reason=molecular_seed mode=init` | the same |

Both reproduce the pilot's conversion measurements (L7 section 7.2, READ:
234 cells clipped and 18.2 percent for `local`, 480 capped and 27.6 percent
for the handoff), so the conversion on the consolidated tree is the
conversion the pilot measured.

## 3. T-L7-5 on the consolidated tree

Route, the recipe of the plan's section 5.3, identical for both (the file is
byte-identical to `.L7/route213_local/input.inp`): `Load IC? True`,
`Reconstruction scheme: PLM` with no `du_th`, `Well balanced: True`,
`Solver: Newton`, `Restart intent: stationary equilibrate`,
`Secondary_ionization: Immediate`, `Molecular carrier transport: True`,
`EXHALE_PTC_DTAU0=1.0`, 8 threads, a 60 minute cap on each. `CFL: 1.0e-12`
belongs to the post-processing pass and is not in this run.

### 3.1 The state as loaded

MEASURED, the certification block the run prints before the first solve:

| row | `local` | handoff | tolerance |
|---|---|---|---|
| hydrodynamic mass | 1.001 at cell 1 | 1.001 at cell 1 | 3.0e-12 (anchored per cell) |
| hydrodynamic momentum | 1.800e-01 at cell 1 | 1.800e-01 at cell 1 | 1.0e-08 |
| hydrodynamic energy | 1.027 at cell 2 | 1.989 at cell 203 | 1.0e-06 |
| carrier balance H2, gated at r >= 1.2 | 9.953e-01 at cell 218 | 9.216e-01 at cell 224 | 1.0e-05 |
| elemental transport He/H, gated | 2.403e-04 at cell 409 | 2.251e-02 at cell 372 | 1.0e-05 |
| level balance He 2^3S | 1.114e-19 | 1.571e-11 | 1.0e-06 |
| eliminated-species closure `System_HeH_mol` | 4.080e-17 | 8.677e-07 | 1.0e-06 |
| cells without a chemical root | 0 | 0 | -- |
| verdict | NOT CERTIFIED, 5 entries | NOT CERTIFIED, 5 entries | |

The `local` seed's carrier row reproduces L7b's 9.953e-01 at cell 218 to four
digits, and the 163 cells reported as absent with the first at cell 194 are
L7b's count and first cell exactly. The seed is the same seed.

### 3.2 The solve, 60 minutes each, 8 threads

MEASURED. **Neither partition completes a single outer pass, and neither
reaches a certified state**; both runs were ended by their 60 minute cap
(exit 124).

| | `local` | handoff |
|---|---|---|
| JFNK iterations taken in 60 min | 124 | 54 |
| residual norm at the start / lowest / last | 1.867 / 1.062 / 1.673 | 2.166 / 2.021 / 2.055 |
| `dtau` first iteration -> last | 2.00 -> 1.00e+14 | 2.00 -> 1.00e+14 |
| line-search `lam` at the last iteration | 7.81e-03 | 1.95e-03 |
| outer passes completed (`done info=`) | 0 | 0 |
| cells outside their row's tolerance at the last iteration | mass 500, momentum 425, energy 500 | -- |
| the worst row throughout | energy/momentum of cells 3 to 4, r = 1.001 | mass of cell 379 and 411, r = 4.4 to 6.9 |

For the handoff partition this is the third independent measurement of the
same outcome (L7 section 7.2: no descent, `||R||` 1.93 to 2.73, zero outer
passes in 14 minutes; L7b section 7.1: 62 iterations, `||R||` 1.769 to 2.894,
no outer pass), and it is the partition's own doing: carrying x2 = 0.99998 to
30 R_p removes the H I photoionization heating the wind is held up by.

For `local` it is **new and it is worse than L7b**. On the L7b tree the same
route from the same seed reached `info = 0` on its first pass -- MEASURED in
that run's own log, `||R||` 1.867 at iteration 1 falling to 2.68e-03 by
iteration 24, with `lam` = 1 at every iteration and 55 iterations over the
whole pass series. Here the FIRST iteration is identical to the last digit
printed (`||R||` 1.867E+00, `||Fs||2` 5.54E+00, `lam` 1.00E+00), so the two
runs start from one state and one residual; they part as soon as `dtau`
leaves 1. Section 4 is the control that asks whether the ramp is the cause.


## 4. The pseudo-time ramp, measured against itself

Same seed, same route, same binary, same `EXHALE_PTC_DTAU0=1.0`, 8 threads,
60 minutes; the only difference is `EXHALE_PTC_RAMP_DOUBLE`, which L4g
documents as restoring the arithmetic of every result before 2026-09-14.
MEASURED:

| | ramp doubling ON (the default) | `EXHALE_PTC_RAMP_DOUBLE=0` |
|---|---|---|
| JFNK iterations in 60 min | 124 | 159 |
| residual norm, start -> last | 1.867 -> 1.673 (lowest 1.062) | 1.867 -> 4.487e-02 (lowest 4.261e-02) |
| linear merit, start -> last | 5.540 -> 3.120e-02 | 5.540 -> 3.910e-04 |
| `dtau`, start / last / range | 2.00 / 1.00e+14 / 2.00 to 1.00e+14 | 7.42 / 20.8 / 9.2e-02 to 4.08e+02 |
| iterations at full step, lam = 1 | 8 of 124 | 149 of 159 |
| outer passes completed | 0 | 0 |
| judged distance at the last iteration | 3.214e+08 | 6.184e+06 |

The ramp doubles on every ACCEPTED step, and on this state the accepted steps
are line-search steps of lam = 8e-3: the pseudo-time is raised by a step that
barely moved, so within about forty iterations the shift I/dtau is gone and
the solve is the unshifted Newton on a state it is not near. That is the
opposite regime from the one L4e and L4g were measured on, where the state
was already close and dtau = 1 was the floor that would not lift. Both
observations are consistent with the same reading -- the merit is a poor
ramp signal on these states -- but the clip that repairs one breaks the
other, so the growth rule needs a signal that tells "this step was accepted
because it was good" from "this step was accepted because it was tiny".
Reported here as a measurement; it is the ramp's item, not this one's.

## 5. What the refusing row is made of

### 5.1 The H2 at cell 218 is round-off, not chemistry

READ, `src/modules/lower_atmosphere/lower_column.f90` lines 116-128. With
`u = -23672/T - log10(p_bar) + 6.2645` and `t = 10^u`,

```fortran
qh2 = (1.9845d0 + tenu - sqrt(tenu*(3.9690d0 + tenu)))/2.3670d0
```

Write `A` = 1.9845 and `B` = 2.3670. Multiplying by the conjugate,

    A + t - sqrt(t(2A + t)) = A^2 / ( (A + t) + sqrt(t(2A + t)) )   exactly,

so the exact value falls like `A^2/(2Bt)` for large `t` while the two terms
that are differenced both grow like `t`. The difference drops below
`ulp(t) = t eps` at `t = A/sqrt(2 eps)` = 9.4e7, and beyond that the computed
value is the rounding of `t`, quantized in steps of `ulp(t)/B`. RECONSTRUCTED
in the two forms, in double precision:

| `t` | as coded | conjugate form | relative error |
|---|---|---|---|
| 1e5 | 8.318894e-06 | 8.318889e-06 | 6.0e-07 |
| 1e6 | 8.319244e-07 | 8.319038e-07 | 2.5e-05 |
| 1e7 | 8.341377e-08 | 8.319053e-08 | 2.7e-03 |
| 1e8 | 1.259076e-08 | 8.319054e-09 | 5.1e-01 |
| 10^8.5 | 0.000000e+00 | 2.630716e-09 | 1.0 |
| 1e9 | 5.036303e-08 | 8.319054e-10 | 60 |
| 1e10 | 8.058085e-07 | 8.319054e-11 | 9.7e+03 |

MEASURED over the 500 physical cells of the atomic state the seeds are built
from (`.L7c/src_atomic_HeH2.13`, T 418 to 5838 K, p 1.6e-15 to 9.5e-7 bar):

| | |
|---|---|
| first cell wrong by more than 1 percent | cell 187, r = 1.1184, T = 3311 K, `t` = 2.2e7 |
| cells wrong by more than 1 percent | 186 of 500 |
| cells given exactly 0 where the fit is positive | 166 |
| cells given more than twice the fit's value | 6 |
| worst | 1.480e+04 too large, at cell 219 |

and this is what the `local` seed writes into the state (MEASURED, the
written `Ion_species_IC.txt`, against the conjugate form of the same fit at
the same p and T):

| cell | r [R_p] | T [K] | x2 written | x2 the fit actually gives | ratio |
|---|---|---|---|---|---|
| 190 | 1.1249 | 3475 | 9.85227e-08 | 1.00566e-07 | 0.98 |
| 191 | 1.1271 | 3531 | 7.88181e-08 | 7.60854e-08 | 1.04 |
| 192 | 1.1293 | 3586 | 7.88181e-08 | 5.79638e-08 | 1.36 |
| 193 | 1.1316 | 3642 | 7.88181e-08 | 4.44620e-08 | 1.77 |
| 194 | 1.1340 | 3697 | 0 | 3.43370e-08 | 0 |
| 202 | 1.1542 | 4140 | 3.15273e-07 | 5.47449e-09 | 57.6 |
| 214 | 1.1904 | 4757 | 0 | 6.75319e-10 | 0 |
| 215 | 1.1938 | 4804 | 5.04436e-06 | 5.84413e-10 | 8.6e+03 |
| 216 | 1.1972 | 4850 | 5.04436e-06 | 5.07780e-10 | 9.9e+03 |
| 217 | 1.2007 | 4896 | 0 | 4.42923e-10 | 0 |
| **218** | **1.2042** | **4940** | **5.04436e-06** | **3.87821e-10** | **1.30e+04** |
| 219 | 1.2079 | 4984 | 5.04436e-06 | 3.40831e-10 | 1.48e+04 |
| 220 | 1.2115 | 5027 | 0 | 3.00614e-10 | 0 |

Four cells carrying the SAME x2 to twelve digits at four different
temperatures and pressures, with exact zeros between them, is not a
chemical-equilibrium profile; 5.04436e-06 is 7.88181e-08 times 64 and
3.15273e-07 is the same number times 4, which is what a difference quantized
in `ulp(t)` looks like as `t` doubles.

The same defect is what the run's own chemistry warnings report: with the
`local` seed loaded, the first sweeps print `every molecular equilibrium root
left the physical simplex at 32 cell(s), clamped onto the element budget`
(MEASURED, `run.log`). `q_h2_equilibrium` is the H2 partition of the
molecular starting point in `dissociation_ionization_balance_at_fixed_ne`
(`ionization_equilibrium.f90` line 4223) and of `molecular_limit_stages`
(`constrained_chemical_equilibrium.f90` line 1510), READ, so those starting
points carry the same rounding.

**Reach.** Both call sites are molecular-layout routines, and
`h2_mixing_ratio_base` (`composition.f90` line 292) and the cold molecular
start (`set_IC.f90` line 216) evaluate the fit at the base, T = 226 K and
p = 1e-6 bar, where `u` = -92.5 and there is no cancellation. **No atomic run
is affected**, so the LHS 1140 b atomic campaign is untouched. On the hot
Uranus regression state (`backup/regression/golden/mol_base_handoff`, T 886
to 2615 K) only 18 of 500 cells are wrong by more than 1 percent and the
worst is 4.9 percent, which is why the molecular matrix never showed it.

### 5.2 The term budget of the row at cell 218

The code carries no diagnostic that splits a carrier row into its terms: the
only numbers printed for a row are its residual, its scale
(`row_terms`/`row_terms_phys`) and the cell they are attained at.
`EXHALE_CARRIER_DEBUG=1` adds one line per step (Newton iterations, the
relative residual, limiter touches, base-partition movement) and
`EXHALE_DIFFUSION_CHECK=1` reports the ELEMENT diffusion, not the carrier
rows. The budget below is therefore RECONSTRUCTED from the written state with
the expressions the code uses, named at each line; it is an order-of-size
account and not a reading of the assembled row.

READ, `carrier_residual` (`diffusive_photochemistry.f90` lines 3953-4185).
The row of carrier `ic` in cell `j` is

    res = nrho (f - f_old)/dt                     the time term (killed in the steady evaluation)
        + K_j ( s_R J_j - s_L J_{j-1} )           the diffusive and drift face fluxes
        + adv_j                                   the divergence of the mass row's own face mass flux
        - src_ic                                  the chemistry, carrier_source -> mol_heh_rows

and the measure is `|res|` over the sum of the magnitudes of those terms
(without the time term) plus `1e-20 n_element sigrate`.

MEASURED at cell 218 of the `local` seed: r = 1.20424 R_p, T = 4940 K,
n(H I) = 1.3733e9, n(H2) = 3.5043e3 cm^-3, x2 = 5.044e-6, v = 223.7 cm/s,
dr = 4.040e6 cm, dr/v = 1.806e4 s, r/v = 6.069e6 s. Its neighbour at
cell 217, which is the upwind cell, holds n(H2) = 0 exactly.

RECONSTRUCTED, with first-order upwinding for the advective term and the
code's own expressions for the chemistry: R15 `k3b = 2.8e-31 T^-0.6`
[cm^6 s^-1] and R12 `k_diss = k3b/K_eq` (`mol_rates.f90`, READ), with `K_eq`
evaluated from the code's own definition -- the Roueff et al. (2019) bound
ladder `h2_lev_T`/`h2_lev_g` of `molecular_infrared_data.f90` (302 levels,
nuclear spin inside the weights), `D0_H2_cm` = 36118.11 and the constants of
`parameters.f90`, in the expression `h2_thermochemistry_init` builds the
table with. That gives `Z_H2(4940 K)` = 190.01 and
`K_eq(4940 K) = 1.8839e-20` cm^3.

| term of the H2 row at cell 218 | value [cm^-3 s^-1] | |
|---|---|---|
| material advection, div(n(H2) v r^2)/r^2 | +1.95e-01 | the whole of it: the upwind cell 217 holds n(H2) = 0 exactly |
| collisional dissociation R12, H2 + M -> 2H + M | -4.40e-01 | time constant 7.97e+03 s |
| three-body formation R15, H + H + M -> H2 + M | +4.46e-06 | **five decades** below the dissociation |
| Lyman-Werner photodissociation, and the diffusive and drift faces | not evaluated here | both add to the imbalance |

so the row is advection against dissociation with formation five decades out
of the account. The cell's equilibrium H2 is

    n(H2)_eq = K_eq n(H I)^2 = 1.8839e-20 x (1.3733e9)^2 = 3.55e-02 cm^-3,

against the 3.5043e+03 cm^-3 the state carries: **9.86e4 times** its own
equilibrium by the code's own equilibrium constant. (The Visscher/Koskinen
fit read in its conjugate form puts the same cell at 1.30e4 times
equilibrium; the two chemical-equilibrium statements differ by 7.6 at this
temperature, which is a separate matter and is not what this row refuses on.)

**The cell cannot hold that H2 and cannot have made it.** Its dissociation
time, 7.97e3 s, is less than half the cell crossing time dr/v = 1.81e4 s, so
no stationary state of this column carries it; building 3.5e3 cm^-3 at the
formation rate above would take 7.9e8 s, against a column flow time
r/v = 6.07e6 s; and the gas upwind of it carries none. The measure the code
reports, 9.953e-01, is what such a row reads.

### 5.3 The artifact is frozen into the solution

MEASURED on the L7b route's own files (`.L7/seed213_local/output` as written
against `.L7b/green_cert_L7_5/output` after four outer passes):

| cell | x2 of the seed | x2 after four passes | relative move |
|---|---|---|---|
| 191 | 7.881813933165e-08 | 7.881813938281e-08 | 6.5e-10 |
| 215 | 5.044356915952e-06 | 5.044356932433e-06 | 3.3e-09 |
| 218 | 5.044356915952e-06 | 5.044356935412e-06 | 3.9e-09 |
| 219 | 5.044356915952e-06 | 5.044356936483e-06 | 4.1e-09 |
| 400 | 4.059890826423e-03 | 4.059903480331e-03 | 3.1e-06 |
| 450 | 6.101959246127e-01 | 6.102033582639e-01 | 1.2e-05 |

The H2 DENSITY at cell 218 falls from 3.504e3 to 3.281e3 over those passes,
which is the hydrogen falling with the expansion; the MIXING RATIO does not
move. The carrier is carried with the nuclei and its partition is not
relaxed, so a seeded spike stays and no number of outer passes can clear it.

### 5.4 The tolerance is reachable at that cell

The question the item asks -- whether 1e-5 relative is arithmetically
possible in a cell whose H2 is 2.5e-6 of the hydrogen -- is answered by the
code's own bound. READ, `carrier_row_roundoff = 64 eps` = 1.42e-14
(`diffusive_photochemistry.f90`): that is the fraction of the row's FULL
terms below which an exact solve of the row leaves nothing, and the
certification's measure divides by the row's own terms, with an absolute
floor `1e-20 n_element sigrate` inside the scale. So the arithmetic floor of
the measure at any cell is at most 1.4e-14, **nine decades below the 1e-5 of
the wind gate**, and it does not depend on how little of the element the
carrier holds: the smallness of x2 scales the residual and the scale
together. The L3 estimate for the hydrodynamic rows -- eps times the largest
term formed, over the row's scale -- gives the same statement here.

The refusal at cell 218 is therefore not arithmetic and not a tolerance that
is too tight. It is a state whose H2 is 1.3e4 above its own equilibrium,
placed there by a cancelling difference and held there by the partition's
not being relaxed.
## 6. The change this asks for, and its reproduction

`src/` was not touched (the item is a diagnosis). What the measurement asks
for, in `src/modules/lower_atmosphere/lower_column.f90`, is the conjugate
form of the same function:

```fortran
      double precision function q_h2_equilibrium(p_bar, T) result(qh2)
      ! q_H2 = n_H2/(n_H2 + n_H + n_He) of the Visscher/Koskinen
      ! chemical-equilibrium fit.  The fit's own expression,
      !     (A + t - sqrt(t (2A + t)))/B,   A = 1.9845, B = 2.3670,
      ! differences two terms that both grow like t while their difference
      ! falls like A^2/(2t), so above t = A/sqrt(2 eps) = 9.4e7 it returns
      ! the rounding of t.  Multiplying by the conjugate gives the SAME
      ! number with no cancellation at any t, the cold limit A/B = 0.838403
      ! exactly, and the correct A^2/(2Bt) tail.
      real*8, intent(in) :: p_bar, T
      real*8 :: u, tenu
      u    = -23672.0d0/T - log10(max(p_bar, 1.0d-30)) + 6.2645d0
      tenu = 10.0d0**min(u, 300.0d0)
      qh2  = 1.9845d0**2                                                 &
           / (2.3670d0*((1.9845d0 + tenu)                                &
              + sqrt(tenu*(3.9690d0 + tenu))))
      if (qh2 .gt. 1.0d0) qh2 = 1.0d0
      end function q_h2_equilibrium
```

The `u > 30` branch goes with it: at u = 30 the function is 8.3e-31 and the
branch returns 0, a step at the branch's own threshold of the same kind the
comment in that file records having removed from the cold side in 2026-08-31.
The `qh2 < 0` clamp goes too -- the expression is positive by construction.
`min(u, 300)` keeps `tenu` finite; at that value q is already 1e-300.

MEASURED equivalence of the two forms in double precision, over the range
where the coded one is not cancelling:

| `t` | as coded | conjugate | relative difference |
|---|---|---|---|
| 1e-30 | 0.8384030418 | 0.8384030418 | 0 |
| 1e-10 | 0.8383946252 | 0.8383946252 | 1.3e-16 |
| 1e-2 | 0.7583547490 | 0.7583547490 | 0 |
| 1 | 0.3191274233 | 0.3191274233 | 3.5e-16 |
| 1e2 | 0.0081579474 | 0.0081579474 | 2.1e-13 |
| 1e4 | 0.0000831740 | 0.0000831740 | 8.2e-09 |

and the cold limit of the conjugate form at t = 0 is 0.8384030418250951,
which is 1.9845/2.3670 to the last bit.

**What it would move.** No atomic run: both cancelling call sites are
molecular-layout routines (section 5.1). Of the molecular regression cases,
only the cells above about 3000 K at nanobar pressures change, and on the
`mol_base_handoff` golden state that is 18 cells of 500 with a worst change
of 4.9 percent in the STARTING POINT of the chemical solve, not in a
converged state; whether any golden moves at all has to be measured by
running the seven `mol_*` cases with the change, which this item did not do
(no build was allowed). What it would move on LHS 1140 b is the whole
`local` partition above cell 187 and the molecular starting point of every
sweep above 3300 K.

## 7. The recipe, and the held patch

**The molecular recipe is not confirmed**, so
`LHS1140b/models/run_case.sh.molecular.patch` was NOT rewritten. The recipe
of the plan's section 5.3 is the right route -- it loads, it equilibrates,
every cell has a chemical root, the element and level rows are inside -- but
with the seed the tree writes today the H2 column above 1.13 R_p is round-off
and the carrier row cannot be satisfied by any solve. The order this puts
the remaining work in:

1. the conjugate form of `q_h2_equilibrium` (section 6), with the seven
   `mol_*` cases measured against their goldens;
2. the pseudo-time ramp's growth rule (section 4). Until it is settled, a
   molecular run of this route has to carry `EXHALE_PTC_RAMP_DOUBLE=0`, which
   is the only one of the two that descends on this state;
3. T-L7-5 again on the `local` partition with both of those, which is the
   first run in which the H2 column above 1.13 R_p is the fit's own value and
   the hydrodynamic rows are being driven down while it is measured;
4. only then the recipe, and the `input.inp` of the molecular cases below.

**One defect in the held patch, found while reading it and NOT fixed** (it is
not this item's file to change, and the patch is held unapplied): after the
seed pass the patch restores the case's own `input.inp` and runs the wind
pass with it, but `molecular_scalar_*/HeH*/input.inp` still carries the
cold-march configuration -- `Load IC? False`, `Reconstruction scheme:
PLM+WENO3`, `du_th [PLM,WENO3]`, and no `Restart intent:` or
`Secondary_ionization:` line. `run_case.sh` reads `Load IC?` from that file
(line 137), so the wind pass would start a cold march from `set_IC` and throw
the seed away. For the patch to run the route it documents, each molecular
case's `input.inp` has to carry what the atomic cases carry -- `Load IC?
True`, `Reconstruction scheme: PLM`, no `du_th` -- with `Restart intent:
stationary equilibrate` (the atomic cases say plain `stationary`) and
`Secondary_ionization: Immediate`.

## 8. How to reproduce

`EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00`, from
`$EX/LHS1140b/models`. The binary is the delivered `$EX/EXHALE.x`
(md5 `b261c3287e64149a85d0d2dfd2ee2c74`); nothing was built.

```bash
# the certified atomic state, copied out of the campaign's reach
mkdir -p .L7c/src_atomic_HeH2.13
\cp -f atomic_scalar_gj1132_kzz1e9/HeH2.13/output/*.txt .L7c/src_atomic_HeH2.13/
SRC=$PWD/.L7c/src_atomic_HeH2.13

# the two seeds, invariant p (the default); "local" asks for each cell's own
# chemical-equilibrium fit, the absent variable asks for the handoff's x2
for p in local handoff; do
  d=.L7c/HeH2.13_$p; mkdir -p $d/output
  \cp -f molecular_scalar_gj1132_kzz1e9/HeH2.13/base.inp $d/
  sed 's/^Load IC?.*/Load IC? True/' \
      molecular_scalar_gj1132_kzz1e9/HeH2.13/input.inp > $d/input.inp.seed
  \cp -f $d/input.inp.seed $d/input.inp
done
( cd .L7c/HeH2.13_local   && OMP_NUM_THREADS=8 EXHALE_MOLECULAR_SEED=$SRC \
      EXHALE_MOLECULAR_SEED_X2=local $EX/EXHALE.x > seed.log 2>&1 )
( cd .L7c/HeH2.13_handoff && OMP_NUM_THREADS=8 EXHALE_MOLECULAR_SEED=$SRC \
      $EX/EXHALE.x > seed.log 2>&1 )

# the route of the plan's section 5.3, in the same directory
for p in local handoff; do d=.L7c/HeH2.13_$p
  sed -e 's/^Reconstruction scheme:.*/Reconstruction scheme: PLM/' \
      -e 's/^Load IC?.*/Load IC? True/' -e '/^du_th /d' \
      $d/input.inp.seed > $d/input.inp
  printf 'Restart intent: stationary equilibrate\nSecondary_ionization: Immediate\n' >> $d/input.inp
  ( cd $d && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 timeout 3600 $EX/EXHALE.x > run.log 2>&1 ) &
done; wait

# the pseudo-time ramp control: the same route with the arithmetic of
# every result before 2026-09-14 (L4g's option off)
d=.L7c/HeH2.13_local_rampoff; mkdir -p $d/output
\cp -f .L7c/HeH2.13_local/{input.inp,base.inp} $d/
\cp -f .L7c/HeH2.13_local/output/Hydro_ioniz_IC.txt \
       .L7c/HeH2.13_local/output/Ion_species_IC.txt $d/output/
( cd $d && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 EXHALE_PTC_RAMP_DOUBLE=0 \
      timeout 3600 $EX/EXHALE.x > run.log 2>&1 )
```

The cancellation itself needs no run: it is arithmetic on `(p, T)` of the
atomic state, in the two forms of section 5.1, and every number of the tables
there was produced that way. The frozen mixing ratio of section 5.3 is a
comparison of `.L7/seed213_local/output/Ion_species_IC.txt` with
`.L7b/green_cert_L7_5/output/Ion_species.txt`, both already in the tree.

## 9. Noticed outside this item's scope, reported and not fixed

- `q_h2_equilibrium`'s cancelling difference (section 5.1 and 6). It is the
  finding of this item, and it is not confined to the seed: it is the H2
  partition of the molecular chemical starting point in
  `ionization_equilibrium` and `constrained_chemical_equilibrium` as well.
- The `u > 30` branch of the same function returns 0 where the function is
  8.3e-31, a step at the branch's own threshold. The conjugate form removes
  the need for the branch.
- The pseudo-time ramp doubles on a line-search step of lam = 8e-3 as
  readily as on a full step, and on a far state it therefore removes its own
  shift within about forty iterations (section 4). L4g's default is what
  stops this route descending.
- `run_case.sh.molecular.patch` would run its wind pass as a cold march, and
  the molecular cases' `input.inp` still carries the cold-march
  configuration (section 7).
- The code carries no diagnostic that splits a carrier row into its terms in
  a named cell; `EXHALE_CARRIER_DEBUG=1` gives one summary line per step.
  Every carrier-row account in this document had to be reconstructed from the
  written state. A dump of a row's terms cell by cell, on the model of the
  hydrodynamic rows', would have made section 5.2 a measurement instead of a
  reconstruction.

Processes: the three runs of this item were `EXHALE.x` in
`.L7c/HeH2.13_local` (PIDs 34891, 34893), `.L7c/HeH2.13_handoff` (34892,
34894) and `.L7c/HeH2.13_local_rampoff` (62751, 62752); all three ended at
their own 60 minute cap and none is left running.
