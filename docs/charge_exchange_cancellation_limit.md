# Cell 387: the trajectory of a diagnosis, and the actual cause

*Kept as a diagnostic record. The file name is unchanged so that earlier
references resolve, but the title it carried -- "H--He charge-exchange
cancellation in the constrained molecular solve" -- names a mechanism that was
measured and found not to operate here. The sections below are preserved in
the order they were written, because the sequence is the point: three
successive explanations were proposed for one cell, two of them numerical,
and the one that closed it was neither.*

## The outcome first

Cell 387 at step 0 of the `He/H = 1` molecular arm no longer exists as a
failure. It closed when a **missing reaction channel** was added to the
network -- the non-radiative He+ + H electron capture of Kingdon & Ferland
(1996), whose source is the same Zygelman et al. (1989) calculation the code's
radiative rate already came from. The He+ balance had been incomplete, and no
amount of numerical work on an incomplete balance was going to produce a root
of it. With the channel adopted the He/H = 1 arm runs to 12,000 steps with
**zero** non-root acceptances.

The three explanations, in the order they were advanced:

1. **An unavoidable floating-point limit** from subtracting two nearly equal
   H--He charge-exchange rates -- the original verdict, and the title of this
   file. **False, and measured to be false**: at the cell's own temperature
   the two coded rates differ by 36 orders of magnitude (section 1a), so there
   is no cancellation to lose digits to, and deleting the pair outright moves
   neither the residual nor the condition number.
2. **A finite-difference step mismatched to the logarithmic unknowns** -- the
   rebuttal that replaced it, reasoned from the code rather than from the
   numbers, and reproduced in the section below as it was written. **Real, and
   fixed**: measured against a central-difference reference, MINPACK's own
   rule got the HeH+ Jacobian column wrong by 86% and the H3+ column by 70%.
   But it was **not** the cause of this cell: with the step corrected, cell 387
   still failed, and it still failed when the step was varied by a decade.
3. **Missing physics in the He+ balance** -- what actually closed it.

The methodological lesson is recorded deliberately. Explanations 1 and 2 were
both about arithmetic, and both were argued from plausible mechanisms before
the failing state itself had been dumped and examined. The state was
diagnostic as soon as it was captured: the H+ row is carried by hydrogen
recombination, ten orders above either charge-exchange term. A residual that
will not go to zero is first evidence that the equations are not the right
equations, and only then evidence about how they are being solved.

What did not change is the acceptance rule. The tolerance was never loosened,
the class-4 reporting and the persistent-streak stop stayed exactly as
section 113 of `Update_EXHALE.md` set them, and the cell stayed classified as
unresolved for as long as it was one. That is what kept the question open
until the physics answer was found, and it is why an incomplete network was
caught at all.

## 1. Observed failure

The source memo records the following step-0 acceptance:

```text
(ioniz_eq) NON-ROOT accepted (relaxation amnesty): cell 387 step 0
           r 1.943414  T[K] 1.140E+03  info 1  viol 7.856E-02  res 2.072E-02
```

It also records that the constrained continuation reduced the largest scaled
reaction residual to `1.7e-5`, with the H+ balance as the limiting row, while
the element residual was about `1.4e-9`. The requested root tolerance is
`ieq_res_tol = 1e-6`.

These values are retained as quoted measurements. The run log and a saved
cell-state dump containing the two charge-exchange terms, all other H+ terms,
the final `u`, and the Jacobian were not found beside this document. The
numeric attribution to charge exchange therefore could not be reproduced from
repository artifacts during this review. The code paths and formulas below
were verified directly from the implementation.

## 1a. Measured outcome (2026-08-31): the charge-exchange attribution is false

The diagnostic this memo asked for was built and run on the actual cell-387
step-0 state, captured from the He/H = 1 arm and re-evaluated without the
hydrodynamics. **The premise that the H+ row is a near-cancellation of the two
H--He rates is wrong, and section 2.1's expectation below - that the residual
measures "a small net rate against a scale containing the two gross rates" -
does not describe this cell.**

At T = 1139.868 K the two coded rate coefficients are

    kcx(He0 + H+) = 1.6991e-60      kcx(He+ + H0) = 1.7452e-15

so the two gross terms in the H+ balance are **-7.75e-47 and +5.74e-11**. They
differ by 36 orders of magnitude: there is no cancellation between them, and
both sit about 10 orders below the terms that actually carry the row. The
dominant terms, before scaling, are

| term | value |
|---|---|
| recombination of H+ | **-2.1743e+01** |
| R17 He+ + H2 | +3.0153e-03 |
| R10+R13 H+ + H2 | -6.6962e-05 |
| photoionization of H0 | +6.4515e-06 |
| R9 H2+ + H0 | +1.1825e-07 |
| Penning He(2^3S) + H0 | +8.9564e-08 |

Deleting the H--He pair outright changes the condition number of the scaled
Jacobian not at all (1.17130e11 against 1.17131e11) and the full residual from
2.2255e-04 to 2.2261e-04.

The other two candidates in section 3 were also measured. Residual precision
is **not** limiting: summing every signed term of the H+ and He+ rows in
real*16 and comparing with binary64 gives relative summation errors of
1.6342e-16 and 1.3416e-16, confirming section 3.1. The difference-step
mismatch of section 3.3 **is** real - against a central-difference reference
MINPACK's rule got the HeH+ column wrong by 86% and the H3+ column by 70% -
and was corrected. Conditioning is present (smallest singular value 1.60e-11,
five of ten singular values below 1e-8) but does **not** by itself discriminate,
because a cell that converges in seven field solves was measured at condition
number 2.03e16, three orders worse.

What finally resolved cell 387 was neither numerical: the He+ balance was
**incomplete**. The non-radiative He+ + H channel was missing from the network
(see `molecular_chemistry_audit_he_rich.md` section 4.3.1). With it adopted,
the He/H = 1 arm runs 200 steps with **zero** non-root acceptances and cell 387
no longer appears at all, against two non-roots and a res = 2.072e-02 event at
that cell before.

One number from the earlier review round could **not** be reproduced and is
therefore not carried here: a claim that on the `mol_diffusion` golden this
pair supplies 11% of the H+ recombination sink at 1213 K with He/H = 0.0793.
Recomputing at that golden's base cell gives He/H = 221.9 and a ratio far above
11%, so either a different cell or a different convention was meant. Treat it
as unverified.

## 2. What the implementation establishes

### 2.1 The H--He pair is assembled with opposite stoichiometric signs

`src/modules/radiation/charge_exchange.f90` evaluates

```text
R1 = k1 n(He I) n(H+)      He I + H+ -> He+ + H I
R2 = k2 n(He+) n(H I)      He+ + H I -> He I + H+

H+ row  += R2 - R1
He+ row += R1 - R2
```

The constrained residual calls the same routine and then multiplies every
reaction row by a fixed reciprocal turnover scale. The turnover definitions
for the molecular H+ and He+ rows include both H--He rate coefficients. Thus
the reported normalized residual is intended to measure a small net rate
against a scale containing the two gross rates.

*Measured correction (2026-08-31): that is the intent of the normalization,
but it is not what happens at cell 387, where the two gross rates are 36
orders of magnitude apart and both are negligible against the recombination
term. See section 1a. The two directions are also not each other's reverse -
B2 is radiative charge transfer and B1 is collisional - so the "pair" framing
of this section is itself misleading; see `molecular_chemistry_audit_he_rich.md`
section 4.3.*

### 2.2 The constrained variables differ from the existing analytic-Jacobian variables

The constrained solve uses

```text
n_i = n_i,ref exp(u_i)
```

for selected species densities, appends explicit element-conservation rows,
and computes electron density from charge neutrality. In contrast,
`he_h_cx_jac` is written for the atomic fraction layout

```text
x1 = n(H+)/n_H
x2 = n(He+)/n_He
x3 = n(He++)/n_He.
```

It supplies only the H--He contribution for that three-column layout. It
cannot be inserted unchanged into the constrained solver. A correct Jacobian
for the constrained system must also apply the exponential chain rule,
include the dependence of all reaction rows on electron density, include
element-conservation derivatives, and follow the active species partition.

### 2.3 `hybrd` does not accept a partial analytic Jacobian

The constrained path calls `hybrd`, which rebuilds the complete Jacobian via
`fdjac1`. Supporting an analytic or mixed Jacobian requires one of the
following concrete changes:

- call the MINPACK routine that accepts a complete user Jacobian, if that
  routine is present and linked;
- add a complete Jacobian callback and use the existing analytic Newton
  driver after adapting it to the constrained system; or
- modify the local solver interface so selected columns or terms can be
  replaced in a controlled way.

Calling `he_h_cx_jac` before `hybrd` would have no effect because `fdjac1`
would overwrite the matrix.

## 3. Problems with the original numerical diagnosis

### 3.1 Cancellation in the residual is not yet at the binary64 floor

Let the turnover scale `S` be comparable to or larger than `|R1| + |R2|`, as
the implemented molecular scale intends. The floating-point error of forming
`R2 - R1`, divided by `S`, is then normally of order machine epsilon, roughly
`2e-16`, not of order the cancellation ratio. A net rate that is `1e-8` of
the gross rates loses about eight decimal digits but still retains roughly
eight decimal digits in binary64. That is more than required by a normalized
tolerance of `1e-6`.

This does not prove that the complete residual is accurate: the H+ row is a
sum of many positive and negative terms, and their ordering matters. It does
show that `1e-8` cancellation alone cannot explain a plateau at `1.7e-5`.
The actual gross terms and their summation error must be measured.

### 3.2 A `sqrt(eps)` Jacobian error does not imply a `sqrt(eps)` residual floor

For a smooth function evaluated to binary64 accuracy, a forward difference
with a step near `sqrt(eps)` usually has a derivative error of order
`sqrt(eps)`. A derivative accurate to about eight decimal digits is normally
adequate to reduce a residual to `1e-6`. A residual plateau requires additional
evidence such as an ill-conditioned Jacobian, a badly scaled column, an
incorrect active-row set, a nonsmooth residual, or solver stagnation.

No singular values, condition estimate, derivative comparison, or step-size
sweep was recorded in the original memo. The statement that no solver tuning
can cross the plateau was therefore stronger than the evidence.

### 3.3 The actual difference step is poorly matched to logarithmic departures

`fdjac1` uses `h = eps*abs(x(j))`, with the special value `h = eps` only when
`x(j)` is exactly zero. Here `x` is `u`, a logarithmic *departure* from the
rung reference, not the logarithm of an absolute density. Values close to
zero are expected near every successful continuation state.

Consequently:

| `u(j)` | current `h` with `eps ~= 1.5e-8` |
|---:|---:|
| `0` | `1.5e-8` |
| `1e-4` | `1.5e-12` |
| `1e-8` | `1.5e-16` |
| `1e-12` | `1.5e-20` |

The last two perturbations may not produce a resolvable density or residual
change. This behavior is a better first hypothesis than an intrinsic
charge-exchange subtraction limit. Merely changing `epsfcn` does not fix the
underlying scale definition, although a step sweep is useful diagnostically.

A suitable rule for these dimensionless logarithmic departures is based on
an absolute floor, for example

```text
h_j = eta max(1, |u_j|)
```

with `eta` selected by a convergence study. A central difference should be
used in the diagnostic because it separates truncation error from residual
noise more clearly. Production use should account for its doubled residual
evaluation cost.

### 3.4 Positivity is not absolute under the current transform

The code comment says that the exponential cannot produce a nonpositive
density, but only the upper side of `u` is capped. A sufficiently negative
`u` can underflow `exp(u)` to exact zero. This may not occur in cell 387, but
the statement and invariant are not exact. Either bound both sides of `u` or
explicitly accept nonnegative, rather than strictly positive, densities and
handle zero derivatives accordingly.

### 3.5 The physical rate pair has a separate consistency issue

The repository's molecular chemistry audit records that the independently
fitted forward and reverse H--He rates do not satisfy detailed balance. That
is not a numerical explanation for failure to solve the equations as written,
but it matters physically: a highly accurate numerical root can still be the
root of a thermodynamically inconsistent rate pair.

Before interpreting the resulting composition, choose and document one of
these physical policies:

- retain both published fits and state that the network is a kinetic fit
  outside thermodynamic consistency;
- select the better-supported direction in the relevant temperature range and
  derive its reverse from detailed balance; or
- replace both rates with a mutually consistent evaluated pair.

This decision requires checking the published validity ranges and is separate
from fixing the nonlinear solve.

## 4. Assessment of the proposed remedies

### 4.1 Complete analytic or automatic Jacobian: feasible and definitive

This is the strongest long-term solution if derivative diagnostics confirm a
Jacobian problem. It removes finite-difference ambiguity and permits direct
condition estimates. It is more work than the original memo implied because a
complete constrained Jacobian is required, not only the H--He terms.

For a reaction `R = k n_a n_b` and logarithmic unknowns, the useful identities
are

```text
dR/du_a = R
dR/du_b = R
```

followed by the fixed row scaling. Conservation derivatives are the number of
element nuclei in species `i`, multiplied by `n_i` and divided by the element
total. Electron-dependent rates require the derivative of the charge-neutral
electron density. These expressions are simple enough to assemble from the
reaction stoichiometry rather than maintaining handwritten derivatives in
several layouts.

Feasibility: **high**, but the affected path is the entire constrained network
and needs focused derivative tests for every active network configuration.

### 4.2 Scale-aware finite differences: feasible first correction

Replace the difference step for `u` with an absolute-scale rule and make it a
local feature of the constrained solver, so legacy fraction-layout paths are
unchanged. Compare forward and central differences over several values of
`eta`, such as `1e-4`, `1e-5`, `1e-6`, and `1e-7`. The solution and dominant
Jacobian entries should be stable over an interval, not only at one selected
step.

Feasibility: **very high**. This is the smallest implementation experiment and
directly addresses a verified mismatch between `fdjac1` and the chosen
variables.

### 4.3 Row recombination: mathematically possible, but underspecified

Adding the H+ and He+ balance rows cancels the H--He pair exactly because its
stoichiometric contributions are opposite. A nonsingular left transformation
of the complete independent equation set preserves the root. However, replacing
one row with the sum is useful only if the remaining independent row is chosen
and scaled carefully. The exchange information cannot simply be removed; it
must remain in another independent balance.

The phrase “carry the exchange balance in a variable scaled to the difference”
does not define a square system. Before implementation, construct the
stoichiometric matrix for the active network, select an independent row basis,
and check its rank and condition number at the failing state. A QR- or
SVD-selected row basis is safer than a special case tied only to rows 1 and 2.

Feasibility: **moderate**. It may improve conditioning, but only after the
independent equation basis is demonstrated.

### 4.4 Higher precision: useful as a diagnostic, not the first production fix

Evaluating the residual and numerical Jacobian in higher precision is a good
experiment: if the plateau moves by the expected number of digits, residual
or derivative roundoff is implicated. If it does not move, precision was not
the limiting factor.

Compensated summation can reduce error in a sum of many terms. It does not
recover digits lost when only two nearly equal already-rounded values are
subtracted. Pairwise or sign-separated accumulation may help the complete H+
row, but it must be tested against a higher-precision reference.

Feasibility: **high for diagnosis**, **low priority for production** until an
actual binary64 error floor is measured.

### 4.5 Continuation and partition changes: not ruled out

The original memo ruled out more rungs and different holdout thresholds
because all attempted rungs stalled on one row. That observation makes a
simple large continuation jump unlikely, but it does not exclude:

- a singular or inaccurate first weak-field seed;
- a trace species whose small density has a large derivative contribution;
- a change in the selected independent rows when the active partition changes;
- loss of sensitivity caused by `n_i,ref` and `dn_i/du_i = n_i` for a very
  small reference density.

Record the active species, held species, row-to-species map, and column norms
at every failed rung before excluding the continuation and partition logic.

## 5. Focused diagnostic and implementation plan

### Phase A: make the failure reproducible

Add an opt-in diagnostic for the failing cell and rung. Save:

- temperature, element totals, radiation-field scale, and all rate
  coefficients;
- `reference_density`, final `u`, active and held species, and selected rows;
- every signed term in the H+ and He+ balances before scaling;
- row scales, the full residual, MINPACK status, and accepted/rejected steps;
- the numerical Jacobian before QR factorization.

The diagnostic must be sufficient for a small standalone driver to evaluate
the same residual without running the hydrodynamics.

### Phase B: distinguish residual error, derivative error, and conditioning

At the saved state:

1. evaluate the residual in binary64 and higher precision;
2. compute Jacobians with the current rule, scale-aware forward differences,
   central differences, and independent analytic or automatic derivatives;
3. compare every column using a norm and identify zero or disagreeing columns;
4. compute singular values of the independently verified scaled Jacobian;
5. repeat after removing only the H--He pair, as a diagnosis rather than a
   physical solution.

Interpretation:

| Result | Conclusion | Next action |
|---|---|---|
| Binary64 residual disagrees with higher precision near `1e-6` | residual evaluation is limiting | improve accumulation or precision locally |
| Current finite differences disagree, scale-aware differences agree | step selection is limiting | adopt the constrained-specific step rule |
| All derivatives agree but the smallest singular value is tiny | formulation is ill-conditioned | choose a better independent row basis or variables |
| Removing H--He exchange fixes conditioning | exchange coupling exposes the bad basis | reformulate the H/He balance basis |
| Removing H--He exchange does not help | original attribution is false | inspect the dominant row terms and other columns |

### Phase C: implement the narrowest demonstrated fix

Prefer the first remedy supported by Phase B:

1. constrained-specific scale-aware finite differences;
2. a complete analytic or automatic Jacobian;
3. a rank-verified row transformation;
4. higher-precision row evaluation only if a measured residual floor requires
   it.

Do not loosen `ieq_res_tol` to accept the present state. The tolerance is a
physical balance requirement, whereas the suspected problem is derivative or
equation conditioning.

## 6. Validation criteria

The fix is accepted only if all of the following hold:

1. Cell 387 reaches a physical state with normalized reaction residual at or
   below `1e-6` and element residual at or below the existing conservation
   criterion.
2. The result is stable across at least a decade of admissible numerical
   derivative steps or agrees with the independent Jacobian result.
3. The H+ and He+ balances, including their gross H--He terms, reproduce in a
   standalone residual evaluation.
4. The constrained solve reaches the same root from nearby seeds and does not
   rely on projection.
5. The low-He/H molecular regression and the He/H = 1 focused arm show no new
   accepted non-roots. Atomic-only paths need no rerun if they are untouched.
6. If a row basis changes, numerical rank is verified for every active network
   configuration: molecular H/He, triplet, oxygen carriers, and metals.

Until these checks pass, class 4 remains a diagnostic containment mechanism,
not a resolution of the chemistry state.

## 7. Relevant implementation paths

- `src/modules/radiation/charge_exchange.f90`: `he_h_cx_fvec` and the
  atomic-fraction `he_h_cx_jac`.
- `src/modules/nonlinear_system_solver/constrained_chemical_equilibrium.f90`:
  logarithmic density variables, continuation, active-row selection, scaled
  residual, and the call to `hybrd`.
- `src/modules/nonlinear_system_solver/fdjac1.f90`: the current
  `sqrt(eps)*abs(x)` difference rule.
- `src/modules/nonlinear_system_solver/System_HeH_mol.f90`: molecular balance
  rows and turnover scales, including both H--He rates.
- `src/modules/radiation/ionization_equilibrium.f90`:
  `normalized_reaction_residual`, acceptance classes, and the persistent
  non-root stop.
- `docs/molecular_chemistry_audit_he_rich.md`: physical detailed-balance
  caveat for the independently fitted H--He rates.
- `docs/supersonic_molecular_base.md`: the broader hydrodynamic failure and
  required chemistry acceptance criteria.
