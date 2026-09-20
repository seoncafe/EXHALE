# Item D3a: the transported He II row measured, source channel by channel

The measurement side of item D3 of `docs/PLAN_20260918_rev2.md` (steps 1 to
4 of its item text). D3b, the decision on what the certification should ask
of a transported ionization stage, is the user's and is NOT taken here.
Nothing in the acceptance interface was changed: no gate, no tolerance and
no denominator moved, and `certification.f90` was not touched.

Every number is MEASURED (run or computed here) unless marked READ.

Binaries, all built with the conda-forge gfortran on PATH into private
object directories, single threaded wherever results are compared.

| binary | md5 | what it is | used for |
|---|---|---|---|
| `EXHALE_D3a_ctl.x` | `b39379a638109fb20952f9e89d47ad5c` | the tree as this item entered it | the solves that produced the two hot-Uranus states |
| `EXHALE_D3a.x` (first) | `3e7e3a7db0ea1bf58458556c6918ae36` | that tree with the channel diagnostic | the evaluations that produced the records |
| `EXHALE_D3a_ctl2.x` | `237fb66c722be016a14c2e9fa6a5ced6` | the tree at the close of this item with THIS item's two production files at their entry text | the control of the byte-identity pair |
| `EXHALE_D3a.x` (final) | `7bcd4f59ded4c6e069f6ce897ad652e3` | the same tree with them as delivered | the measured build of that pair, and the suites |

The byte-identity pair is the last two, built back to back from one tree
with a source manifest taken before and after and verified unchanged
between them, because two other items edited their own files in this tree
while this one ran. Their difference is exactly this item's two files. The
final binary reproduces the certification of the state the records were
taken on to every digit quoted below (`carrier balance He+` gated 8.075e-1,
`H+` 2.831e-4, `He++` 2.193e-4 at cell 248).

---

## 0. Verdict, before the detail

1. **The gross channels of the three ion-stage rows are now an optional
   output of both production kernels**, `mol_heh_rows` in molecular gas and
   `heh_tr_rows` (and `heh_rows`) in atomic gas, written from the same
   factors as the rows themselves. The signed channel sum of a stage
   reproduces that stage's assembled source, after the caller's He <-> H
   charge-exchange pair, to 3.0e-16 of the larger of its two gross sides,
   in both gas branches, on the fixtures and on synthetic states. Default
   off; the key-off cases are byte-identical (section 7).

2. **The He II row of the hot Uranus is stiff by three decades and its
   arithmetic is clean.** At the gated cell 248 of
   `carrier_model_a_newton`, P = 7.6156e2 and L = 7.6034e2 cm^-3 s^-1
   against a net of 1.2153, so |P - L|/max(P, L) = 1.60e-3. The
   cancellation floor of that net, from the number and size of the summed
   channels, is 3.72e-12: **the net source stands 3.3e11 above the floor**,
   and a quadruple-precision re-summation of the same channels reproduces
   the kernel's double-precision net to the last bit. The refusal is not a
   floating-point cancellation.

3. **The refusing residual is an abundance error of 0.16 per cent.** The
   row's local relaxation is d(P - L)/dx(He II) = -1.4472e9 cm^-3 s^-1 per
   unit fraction, linear over six decades of perturbation, so the residual
   -1.2007 corresponds to a correction of 1.58e-3 of x(He II) at cell 248
   and 1.34e-3 at cell 107. That is the same 0.2 per cent composition
   difference L36d measured between the two builds (READ), now read off the
   row itself.

4. **The certification denominator admits a dimensional residual that
   varies by ten decades over one column, and the turnover denominator
   admits a fixed relative abundance error.** At a tolerance of 1e-5 the
   present scale `row_terms_phys` admits an abundance error of 6.0e-13 of
   x(He II) at one cell of the locally equilibrated hot Uranus and 2.1e-2
   at cell 500 of the certified atomic fiducial. The scale with the gross
   production and loss in place of the net source admits **2.0e-5 of the
   stage fraction at every stiff cell measured, on both fixtures**, and
   coincides with the present scale where the stage is not stiff
   (section 4). This is a measurement, not a recommendation.

5. **The transported state is farther from the balance point than the
   locally equilibrated one.** The abundance correction that would close
   the transport is 1.58e-3 of x(He II) on the transported state at cell
   248 and 2.35e-5 on the locally equilibrated state evaluated through the
   same operator, a factor 67. The carrier relaxation at 25 passes under
   its `trust 2.5e-3` bound has not reached the root the row defines
   (section 3.4); this is the one candidate solver defect the measurement
   exposes and it is stated, not acted on.

6. **The recorded normalized order 2.886 is a DIMENSIONAL order of 1.92.**
   The same manufactured column, measured without the normalization, gives
   3.7493e-2, 1.0076e-2 and 2.6089e-3 cm^-3 s^-1 at 250, 500 and 1000
   cells, an observed dimensional order of 1.923, and the stiff column
   measures 1.928 with an acceptance-normalized 2.886 on the same three
   grids. The face-magnitude scale grows as 1/h, exactly as the second
   review said it would, so the normalized order stands one above the
   truncation's own. L36d section 4's reading of 2.886 as "the WENO3
   reconstruction's own order" is corrected (section 9).

7. **The stage-sum entry's floating-point bound is 3.11e-15 for hydrogen
   and 5.33e-15 for helium** times an eddy-cancellation factor g >= 1
   measured at 1.00 to 1.50. Largest entry measured on any production state
   here: 2.203e-16 (H) and 4.195e-16 (He), 0.071 and 0.079 of their own
   bounds (section 6). The entry is left REPORTING and not gating, which is
   D3b's to change.

8. **The sweep and the transport operator read one radiation field.**
   Bitwise equality of the temperature, the four photoionization rates, the
   recombination, collisional and He 2^3S coefficients and the two
   charge-exchange coefficients over every sweep/operator pair in the
   record, on the evaluation runs and on a one-pass solve; the densities
   differ by 1.5e-15 at the evaluation and by the transported trial
   elsewhere, which is the difference this item measures.

---

## 1. What was added, and where

### 1.1 The channels

`ion_residual_core.f90` now declares the channel layout shared by both
kernels: 26 channels, each a single term of the H II, He II or He III
balance, with the stage each belongs to (`stage_chan_row`), the sign it
enters that stage's source with (`stage_chan_sign`) and a name.
`stage_source_from_channels` and `stage_gross_of_channels` are the two
readers.

| channel | stage | sign | reaction |
|---|---|---|---|
| 1 `HII_photo_HI` | H II | + | photoionization of H I, secondaries included |
| 2 `HII_collisional_HI` | H II | + | electron-impact ionization of H I |
| 3 `HII_penning_He23S` | H II | + | He(2^3S) + H Penning branch |
| 4 `HII_H2_photo_diss_ion` | H II | + | H2 + hv -> H + H+ + e- |
| 5 `HII_H2_photo_double` | H II | + | H2 + hv -> H+ + H+ + 2e-, twice |
| 6 `HII_H2p_plus_HI` | H II | + | H2+ + H (R9) |
| 7 `HII_HeII_plus_H2` | H II | + | He+ + H2 (R17) |
| 8 `HII_recombination` | H II | - | H+ + e- radiative recombination |
| 9 `HII_plus_H2` | H II | - | H+ + H2 (R10 and R13) |
| 10 `HII_cx_gain_HeII_HI` | H II | + | He+ + H0 -> He0 + H+ |
| 11 `HII_cx_loss_HeI_HII` | H II | - | He0 + H+ -> He+ + H0 |
| 12 `HeII_photo_HeI_singlet` | He II | + | photoionization of He I (1^1S) |
| 13 `HeII_photo_HeI_triplet` | He II | + | photoionization of He(2^3S) |
| 14 `HeII_coll_HeI_singlet` | He II | + | electron impact on He I (1^1S) |
| 15 `HeII_coll_HeI_triplet` | He II | + | electron impact on He(2^3S) |
| 16 `HeII_rec_from_HeIII` | He II | + | He++ + e- recombination |
| 17 `HeII_rec_to_HeI` | He II | - | He+ + e- recombination, both series |
| 18 `HeII_photo_to_HeIII` | He II | - | photoionization of He+ |
| 19 `HeII_coll_to_HeIII` | He II | - | electron impact on He+ |
| 20 `HeII_plus_H2` | He II | - | He+ + H2 (R17 and R23) |
| 21 `HeII_plus_CO` | He II | - | He+ + CO (D1) |
| 22 `HeII_cx_gain_HeI_HII` | He II | + | He0 + H+ -> He+ + H0 |
| 23 `HeII_cx_loss_HeII_HI` | He II | - | He+ + H0 -> He0 + H+ |
| 24 `HeIII_photo_HeII` | He III | + | photoionization of He+ |
| 25 `HeIII_coll_HeII` | He III | + | electron impact on He+ |
| 26 `HeIII_recombination` | He III | - | He++ + e- recombination |

Channels 10, 11, 22 and 23 are NOT filled by the kernels. The He <-> H
charge-exchange pair is applied by the caller of the rows
(`charge_exchange::he_h_cx_fvec`) with a reservoir and an orientation that
are the caller's: the ground singlet at `he_row_sign = +1` in molecular gas,
the whole neutral helium at `-1` in atomic gas. They are filled by whoever
applies the pair, and the measurement here fills them by calling that one
routine twice, each time with one of its two rate coefficients set to zero.
No rate is written a second time anywhere.

### 1.2 The record

`EXHALE_STAGE_CHANNELS` names up to eight grid cells as a comma-separated
list. Unset, which is every run by default, nothing is reached and no file
is opened. Named, every evaluation of the rows at one of those cells
appends a line to `./output/stage_channels.txt` carrying the gas branch and
the call site, the cell, a call counter, **every argument the kernel was
given** (27 entries common to both branches, 30 more for the molecular
network) and the 26 channels. Recording the arguments is what makes the
record complete: the rows can be re-evaluated afterwards at the very state
the run gave them, which is how every sensitivity below is taken.

The call site is resolved for the molecular kernel and not for the atomic
one: the transport operator is the only caller that asks `mol_heh_rows` for
its production and loss split, so `present(p_Hp)` tells its evaluations
(`MC`) from the local sweep's (`MS`). `heh_tr_rows` has no such argument
and its records are written `A-`. Resolving the atomic call site needs an
argument from `carrier_source`, which is `diffusive_photochemistry.f90` and
belongs to item D4b this increment; the lines are 5020 to 5035
(`mol_heh_rows`, with the optional split already passed) and 5073 to 5089
(`heh_tr_rows` and the pair that follows it). REPORTED, not done.

The append sits in an OpenMP critical region, so the record is safe at any
thread count; the ORDER of its lines is deterministic only at one thread,
which is how every record below was taken.

---

## 2. The ledger identity (step 1)

`src/tests/stage_row_balance/` links the production objects, calls the two
kernels, applies the pair through `he_h_cx_fvec` and checks that the signed
channel sum of each stage is the assembled source of that stage. In
molecular gas the sources are the rows taken whole,
`src(H II) = fvec(1)`, `src(He II) = fvec(2)`, `src(He III) = fvec(3)`; in
atomic gas row (2) is the SUMMED He I balance written He-I-gain positive, so
`src(He II) = -fvec(2) - fvec(3)`.

| row | state | departure of the channel sum from the source |
|---|---|---|
| `the_channel_sums_are_the_assembled_source_molecular` | synthetic molecular cell, 1800 K, metastable tracked | 1.334e-16 |
| `the_channel_sums_are_the_assembled_source_atomic` | synthetic atomic cell, 6000 K, metastable tracked | 0 |
| `the_channel_sums_are_the_assembled_source_atomic_no_metastable` | the same with the level untracked | 0 |
| `the_standard_form_channels_are_its_own_He_II_row` | `heh_rows`, the two-level form | 1.973e-16 |
| `the_recorded_channels_are_the_rows_of_the_recorded_state` | 24 records, hot Uranus, both states | 0, bitwise |
| `the_recorded_channel_sums_are_the_recorded_source` | the same | 0 |
| `the_recorded_channel_sums_are_the_recorded_source` | 36 records, atomic fiducial | 1.533e-16 |
| `the_recorded_channel_sums_are_the_recorded_source` | 10303 records, hot Uranus one-pass solve | 2.961e-16 |

The departure is measured against the larger of the stage's gross
production and gross loss, which is the quantity the identity is a
statement about. The bitwise row is the stronger one: re-calling the kernel
with the recorded arguments returns the recorded channels to the last bit,
so the record carries the whole state the rows were evaluated at.

**The helium `sprod`/`sloss` of `carrier_source` are the positive and the
negative parts of the NET source**, and item D6 documented them as such where
a reader meets them (READ, D6 report and `diffusive_photochemistry.f90` near
the `sprod` block); the proton ledger that D6 repaired reproduces the
assembled source after the charge-exchange pair. Nothing was relabeled here.

---

## 3. The transported He II row, term by term (step 2)

Three states, each measured through the same operator with
`Ionization transport: True`:

- **hot Uranus, transported**: `backup/regression/carrier_model_a_newton` on
  a scratch copy, `Restart intent: stationary`, `EXHALE_PTC_DTAU0=1.0`,
  `EXHALE_OUTER_PASSES=25`, eight threads, then the resulting state
  re-measured through `Restart intent: stationary evaluate` at one thread.
  Final certification: `carrier balance He+` 1.000 at cell 107, gated
  8.075e-1 at cell 248 (its worst gated row), `carrier balance H+` gated
  2.831e-4, `carrier balance He++` gated 2.193e-4, all against 1e-5.
- **hot Uranus, local root**: the same fixture with the key OFF, solved to
  `info = 0` in 17 passes, then evaluated with the key ON. Its composition
  satisfies P = L at every cell, and the same rows read gated 5.329e-2
  (He+ at cell 257), 3.630e-2 (H+) and 7.554e-2 (He++).
- **atomic fiducial, certified**: the L36e 90-pass state of
  `LHS1140b/models/.L26/fid_resolve` with the key on, re-evaluated here.
  `carrier balance H+` 8.620e-7 at cell 150, `He+` 7.147e-7 at 177, `He++`
  2.360e-6 at 183, `elemental transport He/H partition` 8.679e-6 at 500, all
  within; CERTIFIED, reproducing L36e's table (READ).

### 3.1 The rows, hot Uranus transported state

Volumetric rates in cm^-3 s^-1. `D` is the stage flux divergence (the
advective term is zero for an ionization stage, which carries its advection
inside its own face flux), `R = D - (P - L)` the unscaled residual, and
`scale` is `row_terms_phys` as the certification forms it.

| cell | stage | x | P | L | P - L | \|P-L\|/max(P,L) | D | R | scale | measure |
|---|---|---|---|---|---|---|---|---|---|---|
| 248 | H II | 3.980e-4 | 1.07515e3 | 9.92753e2 | 8.2401e1 | 7.664e-2 | 8.1640e1 | -7.608e-1 | 2.6873e3 | 2.831e-4 |
| 248 | He II | 5.254e-7 | 7.61558e2 | 7.60343e2 | 1.2153e0 | 1.596e-3 | 1.4595e-2 | -1.2007e0 | 1.4869e0 | 8.075e-1 |
| 248 | He III | 1.449e-9 | 2.8708e-4 | 2.3018e-4 | 5.6901e-5 | 1.982e-1 | 5.6725e-5 | -1.762e-7 | 8.0320e-4 | 2.193e-4 |
| 107 | H II | 6.077e-7 | 2.63922e2 | 2.64028e2 | -1.0507e-1 | 3.979e-4 | 1.8240e-1 | 2.8747e-1 | 8.2776e1 | 3.473e-3 |
| 107 | He II | 7.492e-11 | 2.27282e2 | 2.26977e2 | 3.0495e-1 | 1.342e-3 | 5.0672e-6 | -3.0494e-1 | 3.0495e-1 | 1.0000e0 |
| 107 | He III | 2.941e-15 | 1.2324e-8 | 1.2458e-8 | -1.3424e-10 | 1.078e-2 | -1.3181e-10 | 2.429e-12 | 5.5667e-9 | 4.364e-4 |

The components of `row_terms_phys` at those rows: the two face
contributions by their magnitudes, the advective magnitudes (zero for a
stage), the absolute NET source and the 1e-20 free-element floor. At cell
248 the He II row's scale 1.4869 is 1.2153 of net source, 2.716e-1 of face
magnitudes and 7.825e-13 of floor; at cell 107 it is 3.0495e-1 of net
source, 7.4e-7 of face magnitudes and 4.610e-10 of floor. **The scale of a
fast stage is its own net source**, which is what the module's header says
of the diagnostic split and what the second review said of the denominator.

### 3.2 The same rows on the locally equilibrated state

| cell | stage | P | L | P - L | D | R | scale | measure |
|---|---|---|---|---|---|---|---|---|
| 248 | H II | 1.13925e3 | 1.13925e3 | 3.337e-10 | 1.1254e2 | 1.1254e2 | 3.5030e3 | 3.213e-2 |
| 248 | He II | 7.86255e2 | 7.86255e2 | 6.439e-11 | 1.8447e-2 | 1.8447e-2 | 3.4615e-1 | 5.329e-2 |
| 248 | He III | 3.2225e-4 | 3.2225e-4 | 1.827e-17 | 8.4648e-5 | 8.4648e-5 | 1.1205e-3 | 7.554e-2 |
| 107 | He II | 2.35725e2 | 2.35725e2 | 3.794e-13 | 7.8489e-6 | 7.8489e-6 | 1.4254e-5 | 5.506e-1 |

Here P = L to the arithmetic floor, so the row is pure transport and the
scale is the face magnitudes alone. This is the conditioning statement of
the second review in production numbers: **near local equilibrium the
denominator collapses onto the same cancellation the numerator sits on**,
and a row whose transport is real reads 5.3e-2 to 5.5e-1 on a state whose
chemistry is exactly balanced.

### 3.3 The certified atomic fiducial

| cell | r [R_p] | stage | P | L | \|P-L\|/max(P,L) | D | R | scale | measure |
|---|---|---|---|---|---|---|---|---|---|
| 150 | 1.043 | He II | 5.79669e2 | 5.78814e2 | 1.475e-3 | 8.5486e-1 | 7.083e-5 | 3.2010e2 | 2.213e-7 |
| 177 | 1.086 | He II | 4.30499e2 | 4.22404e2 | 1.880e-2 | 8.0950e0 | 1.157e-5 | 1.6190e1 | 7.147e-7 |
| 183 | 1.110 | He II | 3.96069e2 | 3.86301e2 | 2.466e-2 | 9.7672e0 | -2.550e-6 | 1.1193e2 | 2.278e-8 |
| 500 | 29.031 | He II | 7.0702e-4 | 1.7094e-5 | 9.758e-1 | 6.8993e-4 | 6.268e-9 | 3.7822e-1 | 1.657e-8 |

Cell 150 is as stiff as the hot Uranus cell 248 (1.5e-3 against 1.6e-3) and
its row certifies with four decades to spare, because the state IS the root
there: the residual is 7.1e-5 against a transport of 0.855. **Stiffness
alone does not refuse a row**; what refuses it is a state that is not the
root of a stiff row.

### 3.4 Sensitivity, the abundance correction, and smoothness

The perturbation respects the budgets: a helium nucleus moved into He II
comes out of the ground singlet and takes its electron into `n_e`; the same
for H II out of H I and He III out of He II.

| state, cell | d(P-L)/dx(He II) | d(P-L)/dln(n_e) | R/(dS/dx), as a fraction of x |
|---|---|---|---|
| hot Uranus transported, 248 | -1.4472e9 | -2.338e-2 | 1.579e-3 |
| hot Uranus transported, 107 | -3.0295e12 | -8.840e-5 | 1.343e-3 |
| hot Uranus local root, 248 | -1.3527e9 | -2.776e-2 | -2.346e-5 |
| hot Uranus local root, 107 | -2.7401e12 | -1.073e-4 | -3.330e-8 |
| fiducial certified, 150 | -6.2787e6 | -1.747e2 | -4.357e-8 |
| fiducial certified, 177 | -8.2339e5 | -1.897e2 | -2.482e-8 |
| fiducial certified, 183 | -5.2746e5 | -1.847e2 | -1.010e-7 |
| fiducial certified, 500 | -9.6571e-4 | -8.384e-7 | -1.596e-7 |

**Smoothness.** The derivative was taken at relative steps 1e-2, 1e-3, 1e-4,
1e-5, 1e-6 and 1e-8. For the He II row of the hot Uranus it reads
-1.447190e9 at ALL SIX steps, to every digit printed; for the H II row of
the same cell it moves from -2.993882e6 to -2.988956e6 over the same range,
converging monotonically. The row is linear in its own unknown over six
decades of budget-respecting perturbation and shows no roughness at any
step: **the residual is a smooth function of the composition** and its
largeness is not a nondifferentiable numerator.

**And the transported state is farther from the balance point than the
locally equilibrated one.** At cell 248 the correction that would close the
transport is 1.58e-3 of x(He II) on the transported state and 2.35e-5 on the
local root, a factor 67; at cell 107 it is 1.34e-3 against 3.33e-8, a factor
4.0e4. The carrier relaxation ran 25 outer passes with the movement bound
`trust 2.5e-3` and did not reach the root; on a stiff row the distance
that has to be closed is 1e-3 of the fraction, and a bounded alternation
that moves at most 2.5e-3 of the state per pass has no difficulty of size
with it. This is a candidate solver statement, and it is REPORTED, not
acted on: the pass budget, the bound and the route are D9's and D3b's.

### 3.5 The cancellation floor

The floor is formed from the channels the kernel returned: each channel is a
product of two or three doubles, so its own relative error is at most 2 eps,
and the sum of m of them carries at most (m - 1) eps of their magnitudes;
the floor reported is (m + 2) eps times the sum of the channel magnitudes.
A quadruple-precision re-summation of the same channels isolates the
assembly's own rounding.

| state, cell, stage | \|net\| | floor | \|net\|/floor | quad-sum departure |
|---|---|---|---|---|
| hot Uranus transported, 248, H II | 8.2401e1 | 5.510e-12 | 1.50e13 | 0 |
| hot Uranus transported, 248, He II | 1.2153e0 | 3.717e-12 | 3.27e11 | 0 |
| hot Uranus transported, 248, He III | 5.6901e-5 | 5.743e-19 | 9.91e13 | 0 |
| hot Uranus transported, 107, He II | 3.0495e-1 | 1.110e-12 | 2.75e11 | 0 |
| hot Uranus local root, 248, He II | 6.439e-11 | 3.841e-12 | 1.67e1 | 0 |
| hot Uranus local root, 107, He II | 3.794e-13 | 1.10e-12 | 0.35 | 0 |

**On the transported state the net source stands eleven to fourteen decades
above its own arithmetic floor.** On the locally equilibrated state it sits
at the floor, which is what a converged local root is. A
higher-precision evaluation of the ROWS is not available without writing
every rate a second time (the kernels are double precision throughout), and
that is forbidden; what is available and was taken is the higher-precision
SUM of the kernel's own channels, which is zero to the last bit, so the
assembly contributes nothing and the whole of the rounding sits inside the
individual products, which the floor above bounds.

### 3.6 One radiation field

`the_sweep_and_the_transported_source_read_one_radiation_field`: over the
sweep/operator pairs in the record, the largest relative difference of the
temperature, `P_HI`, `P_HeI`, `P_HeII`, `P_HeITR`, the four recombination
coefficients, the four collisional ones, the five He 2^3S rates, the two
charge-exchange coefficients and the five molecular photo rates is **0,
bitwise**, over 2 pairs on the evaluation of the transported state, 2 on the
locally equilibrated one and 7 on a one-pass solve. The densities differ by
1.550e-15 and 4.031e-16 at the two evaluations (the operator evaluates the
row at the composition the sweep left) and by 5.625e-16 on the first pass of
the solve. The atomic fiducial's record holds no such pair: a stationary
evaluation of an atomic gas runs the operator alone, and the record cannot
name the sweep's call site there in any case (section 1.2).
This is the source reading confirmed by execution: `carrier_source` sets
`ieq_cell = bg_cell(j)`, one stored cell state for both.

---

## 4. What each candidate norm admits, dimensionally

The comparison the second review asks for: `tolerance x scale`, the
dimensional residual a norm admits, and what that residual is as an error in
the stage fraction, which is the observable. Tolerance 1e-5 throughout, the
value `cert_tol_carrier_wind` carries (READ).

| state, cell (He II row) | scale, present | scale, gross | admitted \|R\| present | admitted \|R\| gross | admitted dx/x present | admitted dx/x gross |
|---|---|---|---|---|---|---|
| hot Uranus transported, 248 | 1.4869e0 | 1.5222e3 | 1.487e-5 | 1.522e-2 | 1.956e-8 | 2.002e-5 |
| hot Uranus transported, 107 | 3.0495e-1 | 4.5426e2 | 3.050e-6 | 4.543e-3 | 1.344e-8 | 2.001e-5 |
| hot Uranus local root, 248 | 3.4615e-1 | 1.5729e3 | 3.461e-6 | 1.573e-2 | 4.402e-9 | 2.000e-5 |
| hot Uranus local root, 107 | 1.4254e-5 | 4.7145e2 | 1.425e-10 | 4.714e-3 | 6.047e-13 | 2.000e-5 |
| fiducial certified, 150 | 3.2010e2 | 1.4777e3 | 3.201e-3 | 1.478e-2 | 5.302e-6 | 2.447e-5 |
| fiducial certified, 177 | 1.6190e1 | 8.6100e2 | 1.619e-4 | 8.610e-3 | 3.473e-7 | 1.847e-5 |
| fiducial certified, 183 | 1.1193e2 | 8.8453e2 | 1.119e-3 | 8.845e-3 | 2.585e-6 | 2.043e-5 |
| fiducial certified, 500 | 3.7822e-1 | 3.7825e-1 | 3.782e-6 | 3.783e-6 | 2.096e-2 | 2.096e-2 |

"Scale, gross" is `row_terms_phys` with the gross production plus the gross
loss in place of the absolute net source, the other terms unchanged.

Two readings, both MEASURED:

- **The gross scale admits 2.0e-5 of the stage fraction wherever the stage
  is stiff**, on both fixtures and over four decades of x, and the number is
  twice the tolerance. That is arithmetic, not coincidence: for a stage whose
  loss dominates its own fraction, d(P-L)/dx is about -(P + L)/x, so
  `tol (P + L)/|d(P-L)/dx|/x` is about `2 tol`. A turnover gate is a gate on
  the relative abundance error, uniformly.
- **The present scale admits between 6.0e-13 and 2.1e-2 of the stage
  fraction over the same set of rows**, ten decades, because near local
  equilibrium it collapses onto the residual's own cancellation and far from
  it (cell 500, where |P-L|/max(P,L) = 0.976) it is the transport alone.
  Where the stage is NOT stiff the two scales coincide to four digits
  (cell 500: 3.7822e-1 against 3.7825e-1), as they must.

Neither reading selects a norm. What they establish is that the two norms
are statements about different quantities: the present one bounds the
residual against the terms the row actually contains, the gross one bounds
the abundance error the residual implies.

---

## 5. The manufactured cases (step 3)

### 5.1 The smooth column, restated dimensionally

`test_manufactured_ionization_column` is unchanged in what it asserts; two
diagnostics were added beside it.

| cells | dr [R_p] | normalized row | dimensional \|R\| [cm^-3 s^-1] | \|R\|/max\|S\| |
|---|---|---|---|---|
| 250 | 4.000e-3 | 7.2833e-6 | 3.7493e-2 | 4.568e-5 |
| 500 | 2.000e-3 | 1.0135e-6 | 1.0076e-2 | 1.228e-5 |
| 1000 | 1.000e-3 | 1.3327e-7 | 2.6089e-3 | 3.179e-6 |

`max|S| = 8.2079e2` cm^-3 s^-1 is the largest exact source of the continuous
problem, a scale that does not move with h. Observed orders 250 to 1000:
**normalized 2.886, dimensional 1.923, on the grid-independent scale 1.923**.

**The recorded 2.886 is therefore restated conditionally and corrected.**
The acceptance scale contains the two face MAGNITUDES divided by the cell
volume; for a smooth nonzero flux that is about 2|F|/dr and grows as 1/h, so
the normalized measure carries one power of h more than the truncation. The
dimensional order is 1.92, which is the conservative divergence's second
order, not a third-order reconstruction. What L36d's extrapolation does is
extrapolate the NORMALIZED measure with the NORMALIZED order, which is
internally consistent because the normalized measure is the quantity the
certification gates; the anchored value 1e-5 is unaffected and was not
touched. What changes is the reading of where the extra power comes from.

### 5.2 The stiff column

`test_manufactured_stiff_stage_column`, new. The same manufactured solution
is given a reaction that is stiff and exact by construction:

```
    P(r) = Lambda + D(r)/2 ,     L(r) = Lambda - D(r)/2 ,
```

with `D(r)` the exact divergence of the continuous stage flux, so
`P - L = D` exactly and `x(r)` is still the exact stationary solution, while
P and L are each about Lambda. Written as an ionization-recombination pair
`P = alpha (1-x) n_el`, `L = beta x n_el`, the local relaxation is
`dS/dx = -(P/(1-x) + L/x)`, which is what turns a residual into an abundance
error. Lambda is set to 1, 1e2, 1e4 and 1e6 times the largest |D| of the
column.

**Nothing in the discretization changes with Lambda.** The residual is the
same truncation at every Lambda because the source is exact. What changes is
every normalized measure built on it:

| grid | dimensional \|R\| | Lambda/max\|D\| | Damkohler | acceptance | turnover | abundance error |
|---|---|---|---|---|---|---|
| 250 uniform | 4.0502e-2 | 1e0 | 3.907e3 | 7.2833e-6 | 3.5313e-6 | 9.5033e-7 |
| 250 uniform | 4.0502e-2 | 1e2 | 3.902e5 | 7.2833e-6 | 2.3279e-7 | 9.6710e-9 |
| 250 uniform | 4.0502e-2 | 1e4 | 3.902e7 | 7.2833e-6 | 2.4658e-9 | 9.6727e-11 |
| 250 uniform | 4.0502e-2 | 1e6 | 3.902e9 | 7.2833e-6 | 2.4673e-11 | 9.6727e-13 |
| 500 uniform | 1.0827e-2 | 1e6 | 8.202e9 | 1.0135e-6 | 6.5956e-12 | 2.5855e-13 |
| 1000 uniform | 2.7965e-3 | 1e6 | 1.683e10 | 1.3327e-7 | 1.7035e-12 | 6.6779e-14 |
| Mixed, 500 | 3.6695e-1 | 1e6 | 1.427e9 | 5.6394e-5 | 2.2354e-10 | 3.6349e-11 |

"acceptance" is `|R|/(faces + |P-L|)`, "turnover" `|R|/(faces + P + L)`,
"abundance error" `|R|/|dS/dx|` as a fraction of the stage fraction. The
Damkohler number is the chemical relaxation rate over the transport rate,
read where the manufactured profile varies.

Observed orders 250 to 1000: **dimensional 1.928, on the grid-independent
scale 1.928, acceptance-normalized 2.886, abundance error 1.928**. The
acceptance-normalized order stands exactly one above the dimensional one,
which is the row `the_acceptance_normalized_order_is_one_above_the_dimensional_one`.

**The Mixed grid** is the production grid the catalog runs on (50 uniform
base cells of 2e-4 R_p under a geometric stretch, out to 30 R_p); the
manufactured profile is defined on [1, 2] and flat above it, which leaves it
twice differentiable and its source zero in the far field. Its interior
dimensional residual is 3.6695e-1 cm^-3 s^-1, 1.4e2 times the 1000-cell
uniform value, because its cells near r = 2 are far wider than 1e-3 R_p.

**The boundary operator.** The two end cells carry the operator's own
boundary rule (no eddy flux through the end faces) and are reported
separately: dimensional |R| 7.9321e-1, 4.2293e-1 and 2.1818e-1 at 250, 500
and 1000 cells, an observed order of 0.93, against 1.92 in the interior, and
an acceptance measure of 7.9645e-5, 2.1189e-5 and 5.4601e-6. **The boundary
rows are first order and stand a decade above the interior**; this is the
same statement the existing test makes with its "with the end cells" column
(8.6047e-5, 2.1983e-5, 5.5589e-6) and it is now attached to the operator
that produces it.

The stiff column's interior dimensional residual (4.0502e-2 at 250 cells)
differs by 8 per cent from the smooth column's (3.7493e-2) because the two
continue the profile differently into the two ghost cells: the smooth column
extends the quintic, the stiff one holds it flat above r = 2 so that the
Mixed grid is the same problem. Their observed orders agree to 0.005.

### 5.3 What the solution error is, and is not, here

The discrete solution error was NOT obtained by solving the discrete
stationary system: the row is assembled through a limited WENO3
reconstruction and solving it would be a second solver written in the test.
What is reported instead is `|R|/|dS/dx|`, the linearization of that error,
which in the stiff limit is the error to first order because the local
balance and not the transport coupling sets the correction. Its observed
order is 1.928, the dimensional one, and it does not depend on Lambda
(9.67e-13, 2.59e-13 and 6.68e-14 at the three grids at Lambda/max|D| = 1e6,
falling as h^1.93). That independence is the row
`the_abundance_error_is_the_grid_independent_measure`.

---

## 6. The stage-sum entry's floating-point bound (step 4)

The identity `sum_k F_k(f) + F_close(f) = N_el(f)` is algebraic: it has no
truncation error, so what a tolerance on it clears is the rounding of the
sums. Counting the rounding events on the path from the face fractions to
the measured difference (`ionization_stage_face_fractions` and
`ionization_stage_face_flux`, READ), with nk carried stages:

```
   nk   forming xclose = 1 - sum_k x_k
    1   xclose * N
   nk   subtracting each eddy term from the closing flux
   2nk  each carried stage's product and its sum with the eddy term
   nk   summing the nk + 1 fluxes in the measure
    1   the final difference
```

that is `5 nk + 2` events, each bounded by eps times the magnitude of the
quantity it rounds. Every such quantity is bounded by

```
   S = |N| + |F_close| + sum_k |F_k| + sum_k |E_k| ,
```

the eddy terms appearing because the same `E_k` is added to a carried stage
and subtracted from the closing one (its own rounding therefore cancels in
the sum exactly, and only the terms above leak). So

```
   |sum F - N|  <=  (5 nk + 2) eps S .
```

The entry divides by `sc = max(|N|, |F_close| + sum_k|F_k|)`, at least half
of `S' = |N| + |F_close| + sum_k|F_k|`, so with `g = S/S' >= 1`

```
   d  <=  2 (5 nk + 2) eps g .
```

With `eps = epsilon(1.0d0) = 2.220446e-16` (conservative by two, since each
rounding is bounded by the unit roundoff eps/2) and g = 1:

| element | nk | bound at g = 1 |
|---|---|---|
| hydrogen | 1 | **3.109e-15** |
| helium | 2 | **5.329e-15** |

**Validated by execution.** `the_stage_sum_identity_is_within_its_rounding_bound_*`
runs the production face-flux routine on columns whose eddy coefficient
spans `K_0 = 0` to `1e15` cm^2 s^-1 and checks every face against its own
face's bound: the largest ratio of an entry to its own bound is 0.0703
(hydrogen) and 0.0580 (helium), with g measured at 1.000 to 1.50 (it
saturates at 1.5 because a dominant eddy term enters `S'` as well).

**Against the production states**, the largest entry any state measured here
attains:

| state | `ionization stage nucleus sum H` | of its bound | `... He` | of its bound |
|---|---|---|---|---|
| hot Uranus, transported | 2.203e-16 | 0.071 | 4.195e-16 | 0.079 |
| hot Uranus, local root | 2.144e-16 | 0.069 | 3.049e-16 | 0.057 |
| atomic fiducial, certified | 2.161e-16 | 0.070 | 4.106e-16 | 0.077 |

(the bound taken at g = 1; a production face with a dominant eddy term
raises it by at most 1.5). Against the broken constructions L36 measured,
1.48e-4 to 3.62e-4 (READ), the bound separates the two by ten decades
(1.48e-4 / 5.33e-15 = 2.8e10).

**The entry still REPORTS and does not gate.** Turning it into a gate is an
acceptance-interface change and belongs to D3b with `certification.f90` in
its ownership. The bound above is what such a gate would take; the proposal
of L36c section 7, 1e-13 (READ), sits 1.27 decades above the helium bound
and 9.2 decades below the smallest broken-construction reading, so it would
also separate them, with less of the margin on the rounding side.

---

## 7. Tests, and the key-off impact

New suite `src/tests/stage_row_balance/run.sh`: 6 rows without a record, 8
on a record with no sweep/operator pair in it (the atomic fiducial) and 9
where the record holds one, 0 failed on every record taken here. The
extended `src/tests/ionization_stage_flux/run.sh`: 45 rows, 0 failed (42
before this item; no existing row changed and none was retuned).

RED before, GREEN after, for the rows this item adds: the channel arrays do
not exist in the entry text, so the suite does not compile against it; the
RED that is meaningful is the one the measurement replaces, namely that the
gross production and loss of a helium stage were not available at all
(`carrier_source`'s comment says so in the entry text, READ), and the GREEN
is the ledger identity above at 3.0e-16.

**Key-off byte identity**, control binary against measured, one thread:

`EXHALE_D3a_ctl2.x` against `EXHALE_D3a.x` (final), one thread:

| case | key | steps | `Hydro_ioniz.txt` | `Ion_species.txt` |
|---|---|---|---|---|
| `wasp_full` | off | 400 | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |
| `carrier_model_a_newton` | on | 25 outer passes | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |

`wasp_full` is atomic with helium, the He 2^3S level and metals, so it
reaches `heh_tr_rows` through the local sweep with no carrier at all;
`carrier_model_a_newton` reaches `mol_heh_rows` and the transported stage
rows through both call sites, the local sweep and the transport operator, at
25 outer passes of the partitioned stationary route. Together they cover
both kernels and both call sites with the diagnostic compiled in and
switched off. `mol_carrier` at its pinned 12000 steps was started and
stopped before it finished; it adds the molecular carrier operator with the
ionization key off, which the key-on case above already runs. No golden was
refreshed and none needed to be.

---

## 8. What the measurements support for D3b, without choosing

The decision is the user's. What D3a establishes, and what each option would
have to carry:

**(a) The gate as it stands.** Supported by: the certified atomic fiducial,
where the present scale is four decades above the residual at the binding
cell of every row and the state IS the root. Against it: the scale it
divides by is the row's own net source wherever the stage is fast, so the
dimensional residual it admits varies by ten decades over one column
(section 4), and on a locally equilibrated state whose chemistry balances to
the arithmetic floor the He II row reads 5.3e-2 to 5.5e-1. Keeping it means
accepting that a stiff stage is judged on a denominator that vanishes with
the numerator.

**(b) Complementary criteria.** The measurements give three that are ready:
the chemical backward error (the turnover norm, which admits a uniform
2.0e-5 of the stage fraction at every stiff cell measured), the abundance
error the residual implies (`|R|/|dS/dx|`, measured here at every cell and
grid independent), and the stage-sum identity with the executed bound of
section 6. A turnover norm alone does NOT certify a transport balance: the
review's arithmetic counterexample stands, and the measurement adds that on
the hot Uranus's transported state the turnover measure of the He II row is
7.888e-4, which still refuses a 1e-5 gate by a factor 79, while on the
locally equilibrated state it is 1.173e-5, which refuses it by 1.17. Both
states are wrong for the transported problem, and only the abundance
correction of section 3.4 says by how much.

**(c) The stage-sum gate.** Its bound is derived and validated
(section 6): 3.11e-15 (H) and 5.33e-15 (He) times g, with the largest
production reading at 0.08 of it and the broken constructions ten decades
above. Enabling it changes `certification.f90` and is an
acceptance-interface change.

**(d) A solver defect the measurement exposed.** One candidate, stated in
section 3.4 and not acted on: on the hot Uranus the transported state sits
1.58e-3 of x(He II) from the root the row defines while the locally
equilibrated state sits 2.35e-5 from it, so 25 outer passes of a relaxation
bounded at `trust 2.5e-3` per pass moved the state AWAY from the balance
point rather than onto it. Nothing here shows the row is wrong: the channels
close to 3e-16, the source is eleven decades above its arithmetic floor, the
row is linear in its own unknown over six decades of perturbation, and the
sweep and the operator read one radiation field bitwise. What is not
established is whether the relaxation can reach that root on this fixture at
any pass budget, which is a measurement D9 could take.

One further observation that belongs to whichever option is chosen: the
present scale and the turnover scale COINCIDE where the stage is not stiff
(cell 500 of the fiducial, 3.7822e-1 against 3.7825e-1), so a criterion
built on the larger of the two changes nothing outside the stiff regime.

---

## 9. Noticed outside this item's scope

1. **L36d section 4's reading of the observed order 2.886 was wrong** and is
   corrected there: it is the order of the NORMALIZED measure, and the
   dimensional order of the same residual on the same three grids is 1.923.
   The anchored value 1e-5 is unaffected, because the extrapolation
   extrapolates the normalized measure with the normalized order, which is
   the quantity the certification gates. The one sentence was edited; no
   number in that memo changed.

2. **The boundary rows of the stage operator are first order**, observed
   0.93 against 1.92 in the interior on the manufactured column (section
   5.2). The existing suite already reported the boundary cells' larger
   measure; what is new is the order. D5 owns the boundary operator.

3. **The atomic kernel's call site is not resolvable inside the kernel.**
   Telling the local sweep's evaluation of `heh_tr_rows` from the transport
   operator's needs one argument from `carrier_source`
   (`diffusive_photochemistry.f90` lines 5073 to 5089), which is D4b's file
   this increment. The molecular branch resolves it for free through the
   production/loss split the operator already asks for.

4. **`row_terms_phys` carries no advective term for an ionization stage.**
   Measured zero in the advective column at every stage row of every state
   here, which is correct (a stage carries its advection inside its own face
   flux, `carrier_residual`, READ) but means the stage rows and the
   molecular carrier rows are scaled on different sets of terms. Reported
   only.

---

## 10. What this item does not do

- It does not change any tolerance, gate, denominator or acceptance
  interface, and it does not touch `certification.f90`.
- It does not enable the stage-sum entry as a gate.
- It does not add the metal charge-exchange terms to the stage sources
  (D7b), so a metal-bearing configuration is outside what the ledger
  identity above covers; both fixtures measured here are metal-free.
- It does not refresh a golden.
- It does not solve the relaxation question of section 8 (d).
