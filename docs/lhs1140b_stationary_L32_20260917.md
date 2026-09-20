# L32: one frozen coupled Newton system, measured

Item L32 of `docs/PLAN_20260917.md` (section 6). A measurement item: no
production numerical change was made, no default moved, and no
preconditioner was designed. Every number below is MEASURED on the frozen
checkpoint of section 2 unless it carries READ, in which case its source is
named beside it.

The binary of the measurements is a private build of the working tree,
`EXHALE_L32.x` md5 `858b34cb7d0d4d600bd425301ffbe30a`, source
`src/modules/time_step/steady_newton.f90` md5
`21db26649b6253bdf2275c786a642842`, which is the tree's text plus the one
print of section 10. The control binary of that print is
`EXHALE_L32_control.x` md5 `c350cc5a391c07949e46801e6fa21d94`, the same
tree with `steady_newton.f90` at the text this item found it in. Both were
deleted with their object directory at the end of the item, as the private
build of an item is; the md5 in `identities.txt` is what the numbers below
were made with. Every
measurement was made single-threaded, `OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1`, on scratch copies. The record is
`docs/audit_20260905/L32_20260917/`.

---

## 1. Verdict

**Of the three explanations the plan names, the measurements support (i):
the directional action is inaccurate. It is not a linear map of its
direction.** That a wider subspace does not repair it is MEASURED
(section 8b). That no preconditioner repairs it is an argument from what a
preconditioner is, not an experiment: right preconditioning replaces
`A x = b` by `A M^-1 y = b`, and where `A(.)` is not linear neither system
is one a Krylov method has a theory for. A preconditioner experiment would
say how far the cycle gets before the nonlinearity binds; it cannot say
that the cycle converges.

Four independent measurements say so. The first three are made on the
operator the solve itself builds, with its own scaling, its own
preconditioner and its own directions; the fourth is on a fifty-cell copy
of the same configuration.

1. **The residual the block reports is not the residual of its step.** In
   the configuration every production block solve runs, the Krylov cycle at
   the checkpoint reports a relative residual of `2.589e-01`. One fresh
   product of the operator on the step it returned gives a TRUE relative
   residual of `5.113e-01`, and the image the Arnoldi relation hands the
   trust region differs from the operator's image of the same step by
   `4.210e-01`, 42 per cent. The basis is orthogonal to `2.164e-13`, so
   neither number is a loss of orthogonality. The `2.597e-01` of the L22
   memo and of the handoff is the reduced least-squares residual of a
   projection, not a statement about the step.

2. **The true residual of the Krylov iterate turns UPWARD with the
   subspace size.** On one operator and one right-hand side, the cycle at
   subspace sizes 40, 80, 160 and 320 stops at 70 products in every case
   above 40 and returns `2.816e-01`; a second orthogonalization pass moves
   the fifth digit. For a linear operator the GMRES residual is
   monotonically non-increasing in the subspace dimension and cannot turn
   upward. Widening the subspace is measured to buy nothing.

3. **The action is not additive.** `A(v1 + v2) - A(v1) - A(v2)` is
   `1.867e+02` of `||A(v1) + A(v2)||` on the first two directions of the
   solve's own Arnoldi basis, and `1.871e-03` on two directions off it. The
   defect falls as the arc to the power `2.29` and `1.96`, which is
   curvature of the residual over the probe arc and neither rounding
   (power -1) nor a decision of the feasible set (power 0). The action is
   homogeneous to `1.40e-14`, so the failure is across directions and not
   in the step length.

4. **The exact solve of the assembled operator does not solve the
   matrix-free one.** On a 50-cell copy of the same configuration, the
   Jacobian assembled column by column from the same action and solved by
   `dgesv` satisfies its own assembled system to `9.41e-12`, and the
   matrix-free action of that step differs from the assembled matrix's
   product by `1.20e+04` relative, the same at unit step length.

**Explanation (iii), a long-range coupling the banded model omits, is
refuted at this checkpoint.** On the Ritz vectors of the three smallest
Ritz values of the preconditioned operator and on the first two Arnoldi
directions, `||(A - A_band) v||/||A v||` is `2.8e-02` to `3.9e-02`. The
diagnostic then splits six columns, the carrier H2 columns of cells 1, 2,
223, 274, 276 and 286, which are the cells those directions live on: the
part of each of the six lying OUTSIDE the band is exactly zero, and with
the radiation frozen it is zero at cells 1 and 2 and `6e-08` to `1.5e-07`
of `1e-02` at the outer four. The whole of the `2.8e-02` to `3.9e-02` is
inside the band, where the band's entries differ from the true column by
`0.09` to `3.5` per cent. The band is not missing a coupling on these
directions; it is a slightly wrong approximation inside its own support.
Six columns are not every column, so this is a statement about the
directions the spectrum and the cycle named and not about the whole
operator.

**Explanation (ii), a stronger preconditioner, is not what the measurements
support, and the measurement 6 comparison is the reason.** Two cycles of 40 and
one cycle of 80 cost 1.5 times the residual evaluations of one cycle of 40
and return a true residual no better than `2.816e-01`; the nonlinear
progress of the outer iteration is unchanged in all three. A preconditioner
changes which linear system the cycle solves; it does not make the sampled
map linear, and the cycle here stops because the map is not linear, not
because the subspace ran out.

**Where the nonlinearity lives.** Nine tenths of the squared additivity
defect sits in 2 cells of 500, cells 1 and 2, and 87 per cent of it in the
momentum rows and 11 per cent in the energy rows. Nine tenths of the band
difference sits in the same two cells. Beside that, the run's own
closure-map line reports the base carrying `3.563e+04` of the closure error,
so the composition elimination of the first two cells is the most strongly
amplified part of the column; whether that amplification is what makes their
residual curve over the probe arc is a plausible reading of the two numbers
together and was not measured here. The binding row of the solve is
elsewhere (carrier H2 of cell 243, `r = 1.317`), so the row that refuses and the
rows that break the linear model are not the same rows.

**What a base-only block would need, and whether it has it.** It would need
the problematic modes to be confined to the base. They are not, and the two
measurements disagree about which statement is being made. The modes that
break the LINEAR MODEL are confined to the base: cells 1 and 2 carry nine
tenths of the additivity defect and nine tenths of the band difference. The
modes the SOLVE is held by are not: the smallest Ritz values of the
preconditioned operator are carried 98.4 to 99.3 per cent by the carrier
rows and their nine tenths is spread over 65, 39 and 4 cells, with cells 1
and 2 among the largest but cells 223 to 287 beside them, and the row that
refuses the state is the carrier row of cell 243. A base-only block would
therefore remove the source of the model's nonlinearity from the coupled
system without touching the row that refuses. That is worth knowing, and it
is not a justification for the base-only block as a cure.

**What this item says the next item is.** Not a preconditioner. The item is
the directional action: the residual of the first two cells is nonlinear
over the arc the action probes at, `sqrt(eps)(1 + ||Y||) = 3.15e-07`, and
that is where a fix would have to act (a shorter arc for those rows, a
row-dependent arc, an analytic derivative of the base closure, or a
formulation of the base rows that is not amplified by `3.6e4`). Whether to
open that is the user's decision, as the plan states.

---

## 2. The checkpoint: what it is and how it was obtained

**The item's own words place the checkpoint at "trust-region iteration 14,
`||R||` 4.241e-01, thirteen linear solves exhausting 40 of 40 at relative
residual 2.597e-01". Those are two different iterates, and the log says so
(READ, `LHS1140b/models/.L22/i3_block/run.log`):**

| trust-region iteration | `\|\|R\|\|` | relative residual the Krylov leg reached | log line |
|---|---|---|---|
| 1 | 4.2666e-01 (entry) | **2.597e-01** | 154, 155 |
| 13 | 4.242e-01 | 1.509e-01 | 263 |
| 14 | 4.241e-01 | 1.337e-01 | 271 |
| 17 (the last the run reached) | 4.237e-01 | 1.339e-01 | 298 |

`2.597e-01` is the FIRST linear solve's, at the entry state, and it is the
largest of the seventeen. Iteration 14 reached `1.337e-01`, and it was an
ACCEPTED step (ratio 1.001, `accepted=T`), not a rejected one; the rejected
steps of that run are iterations 13 and 16.

**The checkpoint of this item is therefore the entry state**, the coupled
block's iterate at outer iteration 1 on the certified
`molecular_scalar_gj1132_kzz1e9/HeH2.13` state. It is the state at which
the `2.597e-01` the plan quotes was measured, it is the worst of the
seventeen, and it is on disk, so no regeneration was needed: it is the
`_IC` pair the run loaded,
`LHS1140b/models/.L22/i3_block/output/{Hydro_ioniz_IC.txt,Ion_species_IC.txt}`.
A second run reaching iteration 14 is reported in section 9.

Frozen copy and identities (`docs/audit_20260905/L32_20260917/`):

```
1b7fe1d8e44da593f0b198949bd4dc31  checkpoint/input.inp
ac3939005cda413a0cd305f989fbbb2c  checkpoint/base.inp
551bff8b9a4c7d87095c0ecde216daf3  checkpoint/output/Hydro_ioniz_IC.txt
050bc77c2e59f8a322ef63b099618d2d  checkpoint/output/Ion_species_IC.txt
```

The exact command every diagnostic run of this item used, on a scratch copy
of that directory:

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 EXHALE_PTC_DTAU0=1.0 \
  EXHALE_OUTER_PASSES=1 EXHALE_JFNK_MAXIT=1 <diagnostic keys> ./EXHALE_L32.x
```

**The system the checkpoint poses** (MEASURED, the run's own report): N =
500 cells, 5 unknowns in each (mass, momentum, energy, the H2 carrier and
the He element row, since this configuration carries `He_diffusion: True`),
2500 rows; band half-widths `kl_jac = ku_jac = 14`, 29 colors, 87 residual
evaluations for one preconditioner build plus the cycle's products;
`||Y|| = 2.0167e+01`; the arc the action probes at `3.1541e-07`;
`||F/Drow||` of the checkpoint `1.3526e+00`; `||R|| = 4.2666e-01`; no
species unknown outside its box and none on a bound; `dgbcon rcond` of the
shifted band `4.180e-06`; zero unresolved colors, zero zero columns.

**What is frozen, and why the measurements are made inside one process.**
The operator of one block iteration is the state, the composition seed, the
column and row scales, the species box, the active bounds, the pseudo-time
shift and the frozen reconstruction weights. The radiation and chemistry
caches that the residual assembly keeps outside `steady_newton` have no
accessor a second process could be restored through (the `l15_dump_state`
comment says so in the source), so every sample that has to share one
operator is taken between one initialization and one exit, with the
composition seed named at every call.

**The pseudo-time shift of the linear solve is zero.** READ,
`steady_newton.f90` lines 15327 to 15329: `idtau_leg = 0` unless
`EXHALE_TR_SHIFTED_LEG=1`. The trust region's Krylov leg therefore solves
the UNSHIFTED scaled system `A_z v = D^-1 J (D v) / Drow` with
`b = -F/Drow`, and that is the operator every measurement below uses.

---

## 3. Measurement 1: is one operator being sampled?

MEASURED by `src/tests/coupled_block_jacobian/coupled_block_linear_system.f90`
at the checkpoint, single process.

| what was repeated | spread | reference |
|---|---|---|
| `F(Y)` twice, nothing in between | **0.0** (bitwise) | `1.6e-12` in row-scale units, READ, the control of `src/tests/residual_determinism` (`docs/TO_BE_DONE.md` entry (AG)) |
| `F(Y)`, an unrelated state, `F(Y)` again | **0.0** (bitwise) | the same |
| `A(v)` twice, nothing in between | **0.0** (bitwise)  | not taken |
| `A(v)`, an unrelated state, `A(v)` again | **0.0** (bitwise)  | not taken |

The unrelated state is the same column with every cell's internal energy
raised by a thousandth, evaluated through the same caches between the
samples.

**The operator is one map and it is bitwise reproducible.** Nothing is
carried between evaluations, and no part of what follows is a
reproducibility defect. The solve's own repeat check agrees: the merit
reproducibility at the iterate is `0.000E+00` and `3.553e-15` at later
iterates (READ, the `[TR] merit reproducibility` lines of
`i3_block/run.log`).

---

## 4. Measurement 2: the perturbation ladder, by row family

MEASURED. For each direction and each arc the forward difference the solve
uses is compared with a CENTRAL difference over the same arc, and the error
is reported as a maximum and as an RMS separately for the mass, momentum,
energy, elemental and carrier rows, with the cell of the maximum. A row
enters the statistics only where its central difference reaches `1e-3` of
the largest central difference of its own family, which is the threshold
the I1 report used; the count of rows that entered is in the last column of
the program's output. The step is cut to nine tenths of what the species
box and the element budget leave on both sides, so both samples are
feasible and the two sides are the same length; the room along each
direction is quoted so that a rung the box cut is visible as such. It was
not cut at any rung here: the box leaves `1.9e+06` to `3.1e+07` arcs in
every direction measured.

Arc `x1` is the arc the route itself probes at, `3.1541e-07`.

**Direction: the scaled right-hand side** (the direction the Krylov cycle
starts from), maximum relative error:

| arc | mass | momentum | energy | elemental | carrier |
|---|---|---|---|---|---|
| x0.01 | 6.12e-05 | 2.99e-04 | 3.61e-05 | 6.01e-05 | 6.18e-05 |
| x0.1 | 4.90e-06 | **7.19e-01** | 4.51e-06 | 4.84e-06 | 4.98e-06 |
| **x1** | 3.60e-07 | **8.16e-02** | 3.54e-07 | 3.56e-07 | 3.55e-07 |
| x10 | 7.44e-02 | 1.00e+01 | 7.33e-02 | 7.38e-02 | 7.42e-02 |

RMS at the route's own arc: mass `1.11e-07`, momentum `1.04e-02`, energy
`1.09e-07`, elemental `1.19e-07`, carrier `1.10e-07`. Worst cells at that
arc: 11 for every family but the momentum rows, 13 for those.

**Direction: the carrier rows**, maximum relative error:

| arc | mass | momentum | energy | elemental | carrier |
|---|---|---|---|---|---|
| x0.01 | 3.17e-03 | 1.52e-03 | 4.46e-03 | **3.50e-01** | 3.32e-03 |
| x0.1 | 4.13e-04 | 2.74e-04 | 9.50e-04 | 2.43e-03 | 4.29e-04 |
| **x1** | 4.13e-05 | **7.15e-01** | 2.98e-05 | **4.71e+00** | 4.22e-05 |
| x10 | 2.55e-01 | 8.84e-02 | 2.58e-01 | 2.55e-01 | 2.57e-01 |

**Direction: the energy rows**, maximum relative error:

| arc | mass | momentum | energy | elemental | carrier |
|---|---|---|---|---|---|
| x0.01 | 2.99e-04 | 1.47e-04 | 4.24e-04 | **1.05e-02** | 2.76e-04 |
| x0.1 | 4.48e-05 | **7.24e-01** | 6.59e-05 | **3.55e+01** | 4.23e-05 |
| **x1** | 1.61e-01 | 1.27e-01 | 3.87e-01 | 1.00e+00 | 1.61e-01 |
| x10 | 3.22e-01 | 2.03e-01 | 3.23e-01 | **4.75e+00** | 3.22e-01 |

**What the ladder says.** For four of the five families along the
right-hand side direction the error behaves as a difference quotient
should: it falls as the arc grows from `x0.01` to `x1` (rounding, error
proportional to the reciprocal of the arc) and rises again at `x10`
(curvature), with its minimum, `3.6e-07`, at the arc the route chose. The
route's arc is well chosen for those rows.

**The momentum rows are the exception and they do not behave like a
difference quotient at all**: `8.2e-02` at the route's own arc, `7.2e-01`
one decade below it, `1.0e+01` one decade above, non-monotone, on 131 of
the rows that the direction moves. The elemental rows behave the same way
along the carrier and the energy directions, reaching `4.71e+00` and
`3.55e+01`. The mean of `3.692e-04` that I1 measured on another state does
not bound these, which is what the review said it would not.

**How far this table may be pushed.** The filter admits a row at `1e-3` of
its family's largest response, so a relative error quoted here can carry up
to three decades of cancellation in its own denominator, and a large entry
is not by itself a large error in the operator. That is why the verdict of
section 1 rests on the norm statements of sections 5, 6 and 8 and not on
this table; what this table adds is WHICH rows and WHICH cells are worst,
and that the momentum rows do not behave like a difference quotient at any
of the four arcs.

The worst cells are not a single place: 11 to 14 in the base layer, 42,
123 to 130, and 216 to 223. Cells 123 to 130 and 216 to 223 are inside the
H2 front, which sits at `r = 1.1597` (cell 204) on this case; cells 11 to
14 are the base layer.

---

## 5. Measurement 3: additivity, homogeneity, and the Arnoldi image

MEASURED. The first three rows are the solve's own hook,
`EXHALE_JV_ADDITIVITY=1`, on the checkpoint's own directions; the last two
are the probe program's on deterministic directions.

| | first two directions of the Arnoldi basis | two directions off the basis |
|---|---|---|
| `\|\|A(v1+v2) - A(v1) - A(v2)\|\|` relative to `\|\|A(v1) + A(v2)\|\|` | **1.867e+02** | **1.871e-03** |
| `\|\|A(v1) + A(v2)\|\|` | 1.473e-03 | 3.760e+03 |
| `\|\|defect\|\|` | 2.750e-01 | 7.037e+00 |
| the defect at a tenth of the arc / at the arc / at ten times it | 1.469e-05 / 1.867e+02 / 5.619e-01 | 3.475e-07 / 1.871e-03 / 2.861e-03 |
| the power of the arc the defect scales as | **2.291** | **1.958** |
| share of the squared defect: mass / momentum / energy / element He / carrier H2 | 0.011 / **0.868** / 0.114 / 0.007 / 0.000 | **0.706** / 0.290 / 0.001 / 0.000 / 0.003 |
| nine tenths of the squared defect sits in | **2 cells of 500, cells 1 and 2** | **1 cell of 500, cell 1** |
| components blocked by the feasible set | 0 | 0 |
| the probe took its nominal step | yes | yes |

**What the exponent says, and what it does not.** The hook's legend is `+1`
for curvature, `-1` for rounding or a residual that is not smooth on the
arc, `0` for a decision of the feasible set. The measured `2.291` and
`1.958` are on the curvature side and far from the other two, and the
blocked counts (0) rule out the feasible set directly, so the defect FALLS
with the arc and is not a floor. It is steeper than the `+1` the first
curvature term of a forward quotient predicts, which the second-order term
alone does not explain and which this item does not explain either. The
second-difference rows of the same hook are ambiguous on those rows and are
reported as such: at half the arc, the arc and twice it the second
difference is `3.320e-02`, `3.574e-01`, `4.224e-01`, ratios `10.765` and
`1.182` against `4` and `4` for a residual smooth on the arc and `1` and
`1` for a floor, so they are neither. The row cancels its own largest term
by a factor `3.555`, which is the cancellation those second differences are
read through.

**The consequence that matters for a remedy.** A defect that falls with the
arc is one a shorter arc would reduce; a floor would not be. The ladder of
section 4 says where the room is: for the mass, energy, elemental and
carrier rows along the right-hand side direction the error is already at
its minimum at the route's own arc, so there is none there, while for the
momentum rows the error is non-monotone in the arc and its minimum is not
at the route's arc.

The defect on the basis directions is `1.867e+02` RELATIVE because the two
images nearly cancel (`||A(v1) + A(v2)|| = 1.473e-03` against a defect of
`2.750e-01`). That is the case that matters: GMRES builds its subspace out
of exactly these directions, and the linear part of their images cancels
while the curvature does not.

**Homogeneity and the Arnoldi image are both exact.**

| | measured |
|---|---|
| `A(3v)` against `3 A(v)`, relative (probe program) | **1.400e-14** (500 cells), **1.170e-16** (50 cells) |
| worst column gap of the Arnoldi image, `A_z M^-1 V_k` against `V Hbar_k` (`EXHALE_GM_IMAGE_CHECK=1`) | **3.523e-16** over 40 columns |
| loss of orthogonality of the basis (`EXHALE_GM_ORTHO=1`) | **2.164e-13**, worst of the solve `3.269e-12` |

So the recursion is faithful column by column and the basis is orthogonal.
The image of the ASSEMBLED step is nevertheless 42 per cent wrong
(section 6), which is only possible because the map is not linear: each
column of the recursion is a correct product of the map on ITS direction,
and the sum over the columns is not the map's product on the sum.

---

## 6. Measurement 4: the recurrence residual against the recomputed one

MEASURED at the checkpoint with `EXHALE_ELEM_DIAG=1`, which makes the solve
recompute `||b - A x||/||b||` by one fresh product of the operator on the
step it returns. Same scaling, same (zero) pseudo-time shift, same
projection.

| configuration | what the cycle reports | TRUE `\|\|b - A x\|\|/\|\|b\|\|` | relative gap of the Arnoldi image | `\|\|sN\|\|` |
|---|---|---|---|---|
| the default cycle (what production runs) | **2.589e-01** | **5.113e-01** | **4.210e-01** | 9.192e-01 |
| the same cycle returning its step by the true residual (`EXHALE_GM_TRUE_RESIDUAL=1`) | 3.189e-01 | 3.189e-01 | 0.000e+00 | 8.900e-01 |

**The number the block reports is smaller than the residual its step
actually has, by a factor 1.98, and the model image the trust region reads
its predicted decrease from is 42 per cent wrong.** The default cycle stops
on the reduced least-squares residual of the projected problem, which keeps
falling because the projection is consistent; the operator's own residual
does not follow it. That is the same shape N27 recorded on the atomic
element reload, where the Arnoldi image was 15 per cent off (READ,
`docs/ISSUES_20260909.md`); here it is 42 per cent.

The source already states this property in the comment of
`gm_step_by_its_true_residual` (READ, `steady_newton.f90` lines 1760 to
1780, item N35). This item measures it on the coupled block at the
checkpoint the L22 campaign refused, and adds the size scan of section 8.

---

## 7. Measurement 5: orthogonality, the residual history, the chemistry tolerance

| what | measured |
|---|---|
| loss of orthogonality of the Arnoldi basis, worst of the solve | **3.269e-12** |
| the level at which a second orthogonalization pass is worth its inner products (READ, `gm_reorthogonalization_loss_level`) | 1e-8 |
| the cycle repeated with the basis orthogonalized twice, subspace 160 and 320 | `2.81571e-01` against `2.81563e-01`, the fifth digit |

Orthogonality is four decades below the level at which it could matter and
reorthogonalization moves nothing. The reduced problem is not losing rank:
where the step is chosen by its true residual the reduced and the true
residual are equal to every digit printed (section 8).

**The preconditioned residual history of one cycle** (`EXHALE_GM_HISTORY=1`
with the size scan, `||b|| = 5.728e+00`), which is the clearest single
statement of the item:

| product | reduced residual | TRUE residual, one product of the operator |
|---|---|---|
| 1 | 9.85110e-01  | not taken |
| 20 | 3.72110e-01 | 3.72989e-01 |
| 30  | not taken | 3.18872e-01 |
| 40 | 2.58916e-01 | **5.11329e-01** |
| 50  | not taken | 2.81563e-01 |
| 60 | 1.75644e-01 | 3.45380e-01 |
| 70  | not taken | risen at two consecutive checks; the cycle stops |

**The reduced residual falls monotonically and the true residual does
not.** At product 40, where the default cycle stops, the reduced residual
has fallen to `2.589e-01` and the operator's own residual of that same
iterate stands at `5.113e-01`, above where it was at product 30. The
subspace is doing what a Krylov subspace does; the map it is built on is
not the map the residual is taken with.

**The chemistry stopping tolerance does not move the attainable linear
residual floor.** `EXHALE_RESID_EQ_TOL` was run at `1e-6`, the default
`1e-8`, `1e-10` and `1e-12` on the checkpoint, everything else equal:

| `EXHALE_RESID_EQ_TOL` | products | true relative residual | residual samples admitted |
|---|---|---|---|
| 1e-6 | 40 of 40 | 3.189e-01 | 80 |
| 1e-8 (default) | 40 of 40 | 3.189e-01 | 80 |
| 1e-10 | 40 of 40 | 3.189e-01 | 80 |
| 1e-12 | 40 of 40 | 3.189e-01 | 80 |

Identical to four digits. The elimination is not what limits the linear
solve here, and the closure-map line of the run says why it cannot be: the
sweep exits at increment `5.003e-14` after one pass where the reference
runs four, and is asked for `2.806e-13`.

---

## 8. Measurement 6: subspace size and restart cycles on one operator

MEASURED, all on the checkpoint, one outer iteration, single-threaded.
**The cost measure is the counted residual evaluations and not the wall
clock**: the host carried other work of this tree throughout, so wall
clocks are not comparable between the runs and none is quoted.

**(a) The three configurations the item names**, each a separate run:

| configuration | products of the operator | true relative residual of the step | residual evaluations of the whole iteration | the nonlinear step |
|---|---|---|---|---|
| one cycle of 40 (`EXHALE_GM_TRUE_RESIDUAL=1`) | 40 of 40 | **3.189e-01** | 80 | accepted, ratio 0.975, `\|\|s\|\| = 1.536e-04`, `\|\|R\|\|` 4.2666e-01 to 4.267e-01 |
| two cycles of 40 (`EXHALE_GM_CYCLES=2`) | 80 of 40 | **3.473e-01** | 123 | the same `\|\|R\|\|` |
| one cycle of 80 (`EXHALE_GM_M=80`) | 70 of 80 | **2.816e-01** | 113 | the same `\|\|R\|\|` |

The outcome line of the last two is "the true residual of the step turned
upward and the best one seen was returned". Memory is the subspace times
the row count: 40, 40 and 80 vectors of 2500 doubles, so 0.8, 0.8 and 1.6
MB; it is not a constraint at this size and it is not what decides.

**(b) The ladder of subspace sizes on ONE operator and ONE right-hand side**
(`EXHALE_KRYLOV_SIZE_SCAN=1`, `||b|| = 5.728e+00`), which is the cleaner
comparison because nothing else can differ between the rungs:

| subspace size | orthogonalized twice | products spent | reduced residual | TRUE residual | `\|\|x\|\|` |
|---|---|---|---|---|---|
| 40 | no | 40 | 3.18872e-01 | 3.18872e-01 | 8.900e-01 |
| 80 | no | **70** | 2.81563e-01 | 2.81563e-01 | 8.667e-01 |
| 160 | no | **70** | 2.81563e-01 | 2.81563e-01 | 8.667e-01 |
| 320 | no | **70** | 2.81563e-01 | 2.81563e-01 | 8.667e-01 |
| 160 | yes | 70 | 2.81571e-01 | 2.81571e-01 | 8.667e-01 |
| 320 | yes | 70 | 2.81571e-01 | 2.81571e-01 | 8.667e-01 |

**Eight times the subspace buys nothing.** The cycle stops at 70 products
whatever it is allowed, because the true residual of its iterate has risen
at two consecutive checks, and it saturates at `2.816e-01` against the
`1.00e-01` asked for; 80, 160 and 320 return the same step to six digits.
The statement that cannot hold for a linear operator is the one inside the
cycle: the true residual of the GMRES iterate is `3.730e-01` at product 20,
`3.189e-01` at 30, `5.113e-01` at 40, `2.816e-01` at 50 and `3.454e-01` at
60 (section 7). A GMRES iterate minimizes the residual over a Krylov space
that only grows, so on a linear operator that sequence is non-increasing.

**The call site is the restart wrapper.** READ, `steady_newton.f90` line
15503: the trust-region step calls `pgmres_with_restarts`, so
`EXHALE_GM_CYCLES` is honored there, which the `products 80 of 40` line of
the two-cycle run confirms MEASURED. The two plain `pgmres` call sites,
lines 17226 and 17401, are the cycle that re-initializes the trust radius
when it is not a positive finite length and the Newton leg of the
non-trust-region branch; neither inherits the wrapper, and the block takes
neither at this checkpoint (`use_tr` is on and the radius is finite).

---

## 9. Measurement 7: the assembled Jacobian and a direct solve, on a small grid

MEASURED on a 50-cell copy of the same configuration (`Grid cells: 50`,
`Base grid [dr,cells]: 1.0e-3 10`), the checkpoint state mapped onto that
grid by `src/utils/map_state_to_grid.py`. 250 rows, 5 unknowns in each. A
diagnostic reference, not a proposal for production: it costs one residual
evaluation per unknown.

Every column is `A_z e_i = D^-1 J (D e_i) / Drow`, sampled with the same
`jacobian_action_of_direction` the Krylov cycle calls, with the same zero
pseudo-time shift, so the assembled matrix and the matrix-free operator are
the same object wherever the action is linear.

| | measured |
|---|---|
| `\|\|b - (assembled A) x\|\|/\|\|b\|\|` after `dgesv` | **9.408e-12** |
| `\|\|A(x) - (assembled A) x\|\|/\|\|(assembled A) x\|\|` | **1.196e+04** |
| `\|\|b - A(x)\|\|/\|\|b\|\|`, the direct step put back through the action | **1.196e+04** |
| the same at unit step length | **1.196e+04** |
| `\|\|x\|\|` of the direct step | 1.300e+02 |
| additivity of the action on this grid | 1.719e-01 |
| homogeneity of the action on this grid | 1.170e-16 |
| repeated residual and action, with an unrelated state in between | 0.0, bitwise |

The direct solve solves the system it was given to eleven digits and the
map it was assembled from does not agree with it at all. Because the defect
is the same at unit length, it is not a property of the step's length: the
action is homogeneous, and what it is not is additive.

**One caveat on the size of `1.196e+04`, stated because it bears on how the
number may be used.** The column probe displaces one unknown by the arc in
that unknown's own units, while a combined direction distributes the same
arc over every unknown; the two therefore sample the residual at different
relative displacements, and the number measures the departure of the action
from ANY single linear map rather than an error against a known Jacobian.
The statements that do not depend on that are the ones of sections 5, 6 and
8, which are taken on the solve's own directions.

---

## 10. Step 0: the iterate record of a loop-top stop

Found by L28 and repaired here. In `solve_steady_jfnk` the
`(JFNK) it N ||R||=` line is printed at the END of the iteration body while
the acceptance test sits at the loop top, so a solve accepted on its first
test printed a completion line and no iterate record, and a reader scoped
to one solve then had no iterate to compare the hand-back measurement with.
One record is now printed at the loop-top stop, in the same format as the
iteration line and carrying the `solve=<n>` token.

**RED and GREEN**, on the reader `src/tests/steady_selfconsistent_residual/run.sh`
and two logs that differ only in that record:

| log | verdict |
|---|---|
| the control format: a solve accepted at its loop top, no iterate record | `FAIL selfconsistent_residual_of_dir measured=incomplete_evidence`, reason `no_iterate_record_in_that_solve` |
| the same log with the record this item prints | `PASS handback_matches_the_accepted_iterate_of_dir measured=1 reference=1 tol=2`, read from the identified scope of solve 2 |

The other four rows of that suite pass on both. The two logs are
synthetic: a first-iteration loop-top stop appears in none of the archived
logs this item read (`.L22/i3_alt/run.log` has eleven loop-top stops and
every one of them follows iteration lines of its own solve), and synthetic
logs are what L28's plan text prescribes for this reader.

**The print cannot move a number, and the control says it does not.**
`backup/regression/carrier_model_a_newton` on a scratch copy,
`EXHALE_OUTER_PASSES=3 EXHALE_JFNK_MAXIT=5`, single-threaded, run with
`EXHALE_L32_control.x` (the entry text) and with `EXHALE_L32.x`: the data
lines of `output/Hydro_ioniz.txt` and `output/Ion_species.txt` are
byte-identical (md5 `bbfa8936ec2e3cc6701b36c3c1a807d4` and
`e00f35f9466f5760fef170e2957b2e00` in both), and the two logs differ only
in the wall clock of the three outer-pass lines.

**A second control of the same print.** The eight-row suite of this
directory, `run.sh`, was run with both binaries on
`LHS1140b/models/molecular_scalar_gj1132_wellmixed/HeH0.55`: every row
reads the same digits on both (`4.472282e-09`, `1.721303e-11`,
`3.6918105E-04`, `3.7921027E-03`, `0_of_13`, `0_of_2`), and the one row
that FAILS, `the_species_rows_scale_is_the_certifications` at `2_of_13`,
fails identically on both. That row is a pre-existing failure of this
tree, not of this item; it compares `diffusive_photochemistry.f90`, which
item L30 is editing.

Suites run on the measured build, being every suite whose driver links
`steady_newton`: `steady_selfconsistent_residual` (all rows pass on its own
run log), `steady_species_rows` with `EXHALE_OBJDIR=build_L32` (195 rows
pass, none fails), `krylov_and_dogleg` (334 rows pass), `element_operator`
(40 rows pass).

---

## 11. What this item did not settle

- **Iteration 14 was not reached by the second reading, and the verdict
  does not rest on it.** A rerun of the block route from the checkpoint was
  started single-threaded on this build with the diagnostics armed at
  iteration 14; after 66 minutes of processor time it stood at iteration 6,
  the host carrying the tree's regression matrix throughout, and it was
  stopped rather than held open. What it did reach is a reproduction
  statement and is in the record as `it14_partial.run.log`: the archived
  run was made at eight threads on an earlier binary, and this
  single-threaded run on the present one walks the same trajectory.

  | trust-region iteration | this run, one thread, `EXHALE_L32.x` | the archived run, eight threads, `EXHALE_L22e.x` |
  |---|---|---|
  | 1 | 2.589e-01, `\|\|R\|\|` 4.267e-01 | 2.597e-01, 4.267e-01 |
  | 2 | 1.679e-01, 4.266e-01 | 1.477e-01, 4.266e-01 |
  | 3 | 1.189e-01, 4.265e-01 | 1.198e-01, 4.265e-01 |
  | 4 | 1.242e-01, 4.265e-01 | 1.252e-01, 4.264e-01 |
  | 5 | 1.327e-01, 4.264e-01 | 1.330e-01, 4.261e-01 |

  Every leg of both spends 40 products of 40 and every iterate of both is
  bound by the carrier H2 row of cell 243. The verdict rests on iteration
  1, which is the state the `2.597e-01` of the plan belongs to and the
  worst leg of the seventeen the archived run reached. All seventeen of
  those iterations carry that same binding row, and sixteen of them end
  with the subspace exhausted at `1.2e-01` to `2.6e-01` against `1.0e-01`
  (the eighth reached its tolerance in 37 products), so nothing in
  section 1 depends on which of them is chosen; a measurement that did
  depend on it would have to say so.
- **No remedy was tried.** Naming what a fix would have to act on
  (section 1) is not a design, and this item proposes none.
- **The Schur-complement question is untouched.** The measurements say a
  stronger preconditioner is not what the evidence supports; they do not
  measure one.
- **Nothing here is a statement about other states.** The checkpoint is one
  iterate of one configuration. The three refusing cases of I4 were not
  measured.

---

## 12. The record

`docs/audit_20260905/L32_20260917/`:

- `identities.txt`, the md5 of the checkpoint files, of both binaries and
  of the two source files this item touched;
- `checkpoint/`, the frozen state, `input.inp` and `base.inp`;
- `d_default.run.log`, the default cycle with the true-residual and image
  checks (section 6);
- `d_linear.run.log`, the additivity hook, the image check, the
  orthogonality measurement (sections 3 and 5);
- `d_scan.run.log`, the subspace-size ladder (section 8b);
- `d_band.run.log`, the Ritz values and the band difference (section 1);
- `d_cyc1`, `d_cyc2`, `d_m80`, the three configurations of section 8a;
- `d_eqtol_m6`, `d_eqtol_m10`, `d_eqtol_m12`, the chemistry tolerance
  ladder (section 7);
- `t1.run.log`, the bare run of the checkpoint with no diagnostic armed;
- `it14_partial.run.log`, the second reading that was stopped at iteration
  6 (section 11);
- `cbj_suite_measured.log` and `cbj_suite_control.log`, the eight-row suite
  of this directory on both binaries (section 10);
- `d_hist.run.log`, the residual history of one cycle (section 7);
- `probe1.log` and `small_probe.log`, the probe program on 500 and on 50
  cells.

The program is `src/tests/coupled_block_jacobian/coupled_block_linear_system.f90`
with `run_linear_system.sh` beside it; the eight-row suite `run.sh` of that
directory is unchanged.
