# Item D8 step 8: the composition movement bound held at 1e-2, measured

The one bounded experiment the user's D3b decision opened (READ,
`docs/PLAN_20260918_rev2.md`, the D3b row of "Decisions made on 2026-09-18",
and D8 step 8). It asks one question of two cases: is the monotone halving
of the composition movement bound what holds them.

Nothing in the acceptance interface moved: no gate, no tolerance, no
denominator, no default. The only source change is a diagnostic environment
key that defaults off, `EXHALE_CARRIER_TRUST_HOLD`.

Every number is MEASURED (run here) unless marked READ.

| binary | md5 | what it is | used for |
|---|---|---|---|
| `EXHALE_bnd_ctl.x` | `39de2ffc837ba2c03320aee18db3b8e0` | the tree at 17:12 with the ENTRY TEXT of `EXHALE_main.f90` compiled in place of the delivered one and relinked | the control of the byte-identity pair and the RED run of the test row |
| `EXHALE_bnd.x`, first | `e9200686c285942ac2517568eaa5ac28` | the same tree with the key as delivered | sections 3 and 4, the byte-identity pair, the GREEN run of the test row |
| `EXHALE_bnd.x`, second | `1abf66c7024f0c2a066814e9749f2f21` | the tree at 19:17 with the key as delivered | section 5, both pairs |

All built with the conda-forge gfortran on PATH into private object
directories, from the working tree at `git=3c73905ca8a7`, which carries the
uncommitted work of the other D-items as well. The two `EXHALE_bnd.x` builds
differ because that tree moved between them: they carry the same key and
different states of the other items' files. Every comparison reported here is
between two runs of ONE binary, so no comparison crosses that difference; the
byte-identity pair is the 17:12 build against its own entry text.

---

## 1. Verdict

1. **On `backup/regression/carrier_model_a_newton` the movement bound is not
   the binding constraint and holding it changes nothing.** Ten further
   passes from the 25-pass transported state with the bound held at 1.0e-2
   and ten with it halving as today produce **bitwise identical states**.
   Every relaxation of both runs ends on "the carriers stopped moving at the
   fixed wind" with an accepted composition displacement of 9.43e-4 falling
   to 4.71e-4, one decade below either bound, and the particle-count bound
   is not attained at any cell in any pass. The halving the control takes
   between passes 7 and 8 therefore has nothing to shorten.
2. **The He II row's distance from the balance point FALLS, and crosses it.**
   The abundance correction that would close the transport at cell 248,
   `R/(d(P-L)/dx)` as a fraction of `x(He II)`, is 1.579137e-3 on the
   25-pass state (reproducing D3a to every digit, READ and re-MEASURED here)
   and **-7.222846e-4** after ten further passes: the magnitude falls by
   0.4574 over ten passes, 7.52 per cent a pass, and the sign changes, so
   the state passed through the root of the row. `P - L` at that cell goes
   from +1.215281 to -5.356961e-1. The candidate solver defect D3a reported,
   twenty-five bounded passes moving the transported state away from the
   balance point, is **not a stall of the alternation and not an effect of
   the movement bound**: the same alternation closes the distance at 7.5 per
   cent a pass, with the bound never binding.
3. **On `molecular_scalar_gj1132_wellmixed/HeH0.083` the bound IS binding,
   and holding it does not make the row fall: it makes it rise.** Every
   relaxation of both runs ended on the movement bound with the bound
   attained. With it held at 1.0e-2 the gated H2 row rises monotonically
   over all ten passes, 7.37e-2 to 9.35e-2, +2.68 per cent a pass. This
   case is not the one L33 section 6.3 designed the experiment for, and it
   does not test the rate proportionality (its bound was already 1.0e-2 and
   stayed there for eight of the ten passes in both runs). The control's
   bound halves once, after pass 8, and that is the
   only pass in which the two differ: the halved step moved the composition
   by 8.26e-3 against the held step's 1.65e-2 and took the row from 9.27e-2
   to 8.30e-2, while the held step took it to 9.35e-2. Shortening the step
   made the row fall and the full step made it rise, the opposite sign to
   the throttle reading.
4. **On `molecular_scalar_gj1132_kzz1e9/HeH0.083`, the case L33 section 6.3
   actually designed the experiment for, the ordering L33 read as a throttle
   runs the other way.** Two corrections first, both MEASURED: the catalog
   case's published pair is byte for byte the SEED of L33's forty passes
   (md5 `1229b961323d7ac710146650da72454a`
   = `.L22/i4_kz0083/output/Hydro_ioniz_IC.txt`), not the state they left,
   so the experiment was run from `.L22/i4_kz0083/output/*.txt`, gated row
   8.188e-04 at cell 499 with the bound at 1.3e-3 (READ); and the ten passes
   were NOT reached, three completing in both runs and a fourth in the
   control before both solves collapsed into a line search (`lam` 9.5e-7 and
   1.9e-6) at a cost that had grown from 7 s to 3660 s a pass. Over what ran,
   with the bound held at 1.0e-2 the row went 8.15e-04, 6.29e-03 (a factor
   7.72 UP in one pass) and 1.93e-03: **at 1.0e-2 there is no decay whose
   rate could be compared with the 1.04 per cent a pass measured at 1.3e-3.**
   L33's own log says the same where it ran at that bound (READ, its passes 1
   to 5: 2.97e-03, 2.52e-02, 7.98e-03, 3.35e-03, 4.10e-03), and its monotone
   decay begins only at passes 10 and 23, at 2.5e-3 and 1.3e-3. The monotone
   descent is a property of the SMALL bound, not a slowed version of a faster
   one. Holding the bound also cost the wind: pass 3 returned `info = 2` in
   both runs with the mass row at 1.65e-01, the control recovered at the
   halved bound and the held run did not reach a pass 4.
5. **No item is opened on the movement bound's monotone halving, and no
   default is recommended for change.** In the case D3a raised the bound
   never binds at all; in the two molecular cases it binds at every pass and
   holding it at 1.0e-2 makes the row worse, not better, and in the L33 case
   it also breaks the wind solve. Nothing measured here supports the reading
   that the halving is throttling a descent, so there is nothing to open.

---

## 2. The key, and what it does not change

`EXHALE_CARRIER_TRUST_HOLD=<value>` is read once at the entry of the
stationary outer iteration, beside the existing `EXHALE_CARRIER_TRUST` and
`EXHALE_OUTER_PASSES`. Set to a positive value it puts the composition
movement bound at that value and takes the halving branch of the progress
control out of the pass; unset or empty or not positive, nothing in the new
code is reached. It is documented at both sites: at the declaration of the
bound's floor, which it defeats, and at the halving it replaces. It is a
diagnostic, not an input key: `input_read` does not know it.

The progress control is otherwise untouched. `comp_omega` is still halved on
the same passes, the no-fall counter still ends the loop on
`outer_no_fall_max`, and the branch for a relaxation that kept no transport
step still prints what it printed.

**Byte identity with the key absent** (MEASURED): `carrier_model_a_newton`
entered from its pinned IC pair, `EXHALE_OUTER_PASSES=3
EXHALE_JFNK_MAXIT=5`, one thread, the control build against the measured
build.

| file | md5, control | md5, measured |
|---|---|---|
| `output/Hydro_ioniz.txt` | `7b0e0300966472bb6444cc3782d95514` | `7b0e0300966472bb6444cc3782d95514` |
| `output/Ion_species.txt` | `254306f0cdce3210329b093647ee7dd1` | `254306f0cdce3210329b093647ee7dd1` |

The two run logs are identical line for line once the wall-clock fields are
removed.

---

## 3. `carrier_model_a_newton`: ten passes from the 25-pass transported state

**The state.** The 25-pass state D3a produced and measured
(`Restart intent: stationary`, `EXHALE_PTC_DTAU0=1.0`,
`EXHALE_OUTER_PASSES=25`, `Ionization transport: True`), copied to scratch
as the restart pair of two fresh case copies. The fixture's own `input.inp`
is used unchanged, as D3a used it, so that the entry state's measurement can
be compared with D3a's. Both runs: `EXHALE_PTC_DTAU0=1.0`,
`EXHALE_OUTER_PASSES=10`, one thread, `EXHALE_bnd.x`; one with
`EXHALE_CARRIER_TRUST_HOLD=1e-2` and one with the key unset.

**The entry state re-measured** (MEASURED, through
`Restart intent: stationary evaluate` with `EXHALE_STAGE_CHANNELS=248` and
the `stage_row_balance` driver of D3a): `carrier balance He+` gated
8.075116e-1 at cell 248, `H+` 2.831155e-4, `He++` 2.193445e-4, and at cell
248 P = 7.615580e2, L = 7.603427e2, P - L = 1.215281, D = 1.459467e-2,
R = -1.200686, scale 1.486896, d(P-L)/dx = -1.447190e9. Every one of these
reproduces D3a's table to the digits it prints (READ), so the tree's other
uncommitted work has not moved this fixture.

### 3.1 The two pass tables

`trust` is the bound in force. `displacement` is the accepted composition
displacement of the pass. The worst gated row is `carrier balance He+` at
cell 248, r = 1.2955 R_p, in every pass of both runs.

| pass | gated He+ row | cell | trust, held | trust, halved | displacement | hydro info | worst mass row |
|---|---|---|---|---|---|---|---|
| entry | 8.075e-01 | 248 | - | - | - | - | - |
| 1 | 8.08e-01 | 248 | 1.0e-02 | 1.0e-02 | 9.43e-04 | 0 | 3.15e-11 |
| 2 | 8.11e-01 | 248 | 1.0e-02 | 1.0e-02 | 8.63e-04 | 0 | 3.56e-11 |
| 3 | 7.79e-01 | 248 | 1.0e-02 | 1.0e-02 | 7.91e-04 | 0 | 2.95e-11 |
| 4 | 7.82e-01 | 248 | 1.0e-02 | 1.0e-02 | 7.25e-04 | 0 | 3.97e-11 |
| 5 | 7.47e-01 | 248 | 1.0e-02 | 1.0e-02 | 6.65e-04 | 0 | 2.76e-11 |
| 6 | 7.52e-01 | 248 | 1.0e-02 | 1.0e-02 | 6.10e-04 | 0 | 3.05e-11 |
| 7 | 7.13e-01 | 248 | 1.0e-02 | 1.0e-02 | 5.60e-04 | 0 | 3.54e-11 |
| 8 | 7.18e-01 | 248 | 1.0e-02 | 5.0e-03 | 5.13e-04 | 0 | 4.57e-11 |
| 9 | 6.77e-01 | 248 | 1.0e-02 | 5.0e-03 | 4.71e-04 | 0 | 3.35e-11 |
| 10 | 6.82e-01 | 248 | 1.0e-02 | 5.0e-03 | 0.00e+00 | 0 | 3.20e-11 |

The two runs print the same gated row, the same displacement, the same
hydrodynamic rows and the same relaxation ending at every pass; the trust
column is the only entry that differs, from pass 8, and the four output
files of the two runs are bitwise equal. Pass 10 takes no update, which is
the loop's own invariant: the pass that ends the iteration hands back the
state its certification was evaluated on.

**Why the bound cannot matter here.** All twenty relaxations ended on
"the carriers stopped moving at the fixed wind" in 32 or 33 transport steps.
The accepted displacement falls geometrically, 9.43e-4 to 4.71e-4, 8.31 per
cent a pass, and no pass printed a "particle-count bound last held" line at
all: the bound was not attained at any cell. The relaxation reaches the
fixed point of its own operator at the held wind, which is the property L33
found separates the run that certifies from the three that do not (READ,
`docs/lhs1140b_stationary_L33_20260917.md` section 2).

### 3.2 The He II row's distance from the balance point

MEASURED through the D3a route at both endpoints, at cell 248,
T = 1.777e3 K. The state after ten passes is the same for both runs, so one
measurement answers for both.

| quantity | entry, 25 passes | after 10 more |
|---|---|---|
| x(He II) | 5.253926e-07 | 5.259604e-07 |
| P | 7.615580e+02 | 7.613384e+02 |
| L | 7.603427e+02 | 7.618741e+02 |
| P - L | 1.215281e+00 | -5.356961e-01 |
| \|P-L\|/max(P,L) | 1.595782e-03 | 7.031294e-04 |
| D, the stage flux divergence | 1.459467e-02 | 1.459412e-02 |
| R = D - (P - L) | -1.200686e+00 | 5.502902e-01 |
| scale, `row_terms_phys` | 1.486896e+00 | 8.073249e-01 |
| measure, \|R\|/scale | 8.075116e-01 | 6.816217e-01 |
| d(P-L)/dx(He II) | -1.447190e+09 | -1.448540e+09 |
| **R/(dS/dx) as a fraction of x** | **1.579137e-03** | **-7.222846e-04** |
| cancellation floor of the net | 3.71723e-12 | 3.72043e-12 |
| \|net\|/floor | 3.27e+11 | 1.44e+11 |

The derivative is the same at every relative step from 1e-2 to 1e-8, as D3a
found, so the row is still linear in its own unknown over six decades and
the abundance correction is a clean reading.

**The row's distance from its own root falls 7.52 per cent a pass and
changes sign.** The gated measure falls more slowly, 1.87 per cent a pass,
because the denominator falls with the numerator: the scale of this stiff
stage is its own net source (READ, D3a section 3.1), and `|P - L|` fell by
more than half over the same ten passes. **The measure and the distance are
not the same quantity here**, and the measure is the slower of the two.

The other two transported rows fall with it (MEASURED, the same two
evaluations): `carrier balance H+` gated 2.831e-4 to 1.293e-4, `He++`
2.193e-4 to 1.008e-4, and the volume-weighted `He+` entry 5.504e-3 to
2.532e-3.

### 3.3 What this says about D3a's candidate defect

D3a reported, and did not act on, that twenty-five bounded passes left the
transported state 1.58e-3 of x(He II) from the row's root while the locally
equilibrated state sits 2.35e-5 from it (READ, D3a section 3.4). Ten further
passes take that 1.58e-3 to -7.22e-4: the alternation is closing the
distance, and it crosses the root between those two states. It is not a
stalled relaxation and the movement bound is not throttling it, since the
bound is never attained. What the measurement leaves open is the RATE: at
7.5 per cent a pass the distance would need about thirty more passes to
reach the locally equilibrated state's 2.35e-5, and the sign change says the
approach is oscillatory rather than monotone, so an extrapolation of ten
passes is a rate and not a forecast.

---

## 4. `molecular_scalar_gj1132_wellmixed/HeH0.083`: ten passes from the state L34c left

**This case is NOT the one L33 section 6.3 designed the experiment for.**
That design is on `molecular_scalar_gj1132_kzz1e9/HeH0.083` (L33's
`i4_kz0083`), whose bound had been cut to 1.3e-3, one halving above its
floor, and whose row was falling at 1.04 per cent a pass there; the
proportionality it predicts is tested in section 5. The case measured here
entered its ten passes with the bound already at 1.0e-2 and held it for
eight of them in both runs, so it does not test that proportionality. It is
kept because it is a measurement of a case that was named and because the
bound binds in every one of its passes.

**The state.** The published pair of the case directory,
`output/Hydro_ioniz.txt` and `output/Ion_species.txt` of 2026-09-18 12:46,
the state the L34c continuation wrote (READ, `not_solved.md` and
`REPRODUCE.md` of that case). The case directory holds no `ENDING` file;
what the run recorded of its own ending is the log's last two lines, "the
stationary solve returned info = 1; output written" and "the state written
is a relaxation snapshot: the run made no stationary claim, so it is written
certified=F with the reason no stationary claim" (READ). Its final no-step
certification (READ, `not_solved.md`) has every hydrodynamic row within
tolerance and `carrier balance H2` at 9.264e-2 of 1.0e-5, worst at cell 1,
gated 7.371e-2 at cell 500. The `_IC` pair beside it is the 40-pass state
that seeded the continuation, and it was NOT used.

Both runs were copied to scratch; nothing in the case directory was written.
The route is the one the brief names and the case's own `input.inp` carries:
`Restart intent: stationary`, `Coupled carrier solve: False`,
`Well balanced: True` (the resolved default, READ,
`EXHALE_resolved.out`), `Secondary_ionization: Immediate`,
`EXHALE_PTC_DTAU0=1.0`, `EXHALE_OUTER_PASSES=10`, eight threads,
`EXHALE_bnd.x`. The only edit to the copied input is the spectrum path,
made absolute. Wall clock 48 min (bound held) and 49 min (bound halving), the two in
parallel on `lart4`, 17:15 to 18:04 KST.

### 4.1 The two pass tables

The worst gated row is `carrier balance H2` at cell 500, r = 29.0312 R_p,
the outermost physical cell, in every pass of both runs. The bound was
attained in every relaxation of every pass of both runs: each ended on the
movement bound and each printed the bound as held at a cell that walks
outward, 271 (1.516 R_p) to 292 (1.745 R_p).

| pass | gated H2 row, bound held | gated H2 row, bound halving | trust, held | trust, halving | displacement, held | displacement, halving | hydro info | worst mass row | x2=0.5 [R_p] | x2=1e-2 [R_p] |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 7.37e-02 | 7.37e-02 | 1.0e-02 | 1.0e-02 | 1.73e-02 | 1.73e-02 | 0 | 1.15e-07 | 1.4044 | 18.8174 |
| 2 | 7.67e-02 | 7.67e-02 | 1.0e-02 | 1.0e-02 | 1.73e-02 | 1.73e-02 | 0 | 5.95e-08 | 1.4261 | 17.6173 |
| 3 | 7.95e-02 | 7.95e-02 | 1.0e-02 | 1.0e-02 | 1.72e-02 | 1.72e-02 | 0 | 5.79e-08 | 1.4413 | 16.4981 |
| 4 | 8.21e-02 | 8.21e-02 | 1.0e-02 | 1.0e-02 | 1.71e-02 | 1.71e-02 | 0 | 3.94e-08 | 1.4570 | 15.4543 |
| 5 | 8.46e-02 | 8.46e-02 | 1.0e-02 | 1.0e-02 | 1.70e-02 | 1.70e-02 | 0 | 1.52e-08 | 1.4733 | 14.2478 |
| 6 | 8.71e-02 | 8.71e-02 | 1.0e-02 | 1.0e-02 | 1.69e-02 | 1.69e-02 | 0 | 4.59e-08 | 1.4987 | 13.3555 |
| 7 | 8.93e-02 | 8.93e-02 | 1.0e-02 | 1.0e-02 | 1.68e-02 | 1.68e-02 | 0 | 1.19e-08 | 1.5165 | 12.5233 |
| 8 | 9.12e-02 | 9.12e-02 | 1.0e-02 | 1.0e-02 | 1.67e-02 | 1.67e-02 | 0 | 1.68e-08 | 1.5443 | 11.7471 |
| 9 | 9.27e-02 | 9.27e-02 | 1.0e-02 | 5.0e-03 | 1.65e-02 | 8.26e-03 | 0 | 5.30e-08 | 1.5636 / 1.5539 | 11.1995 / 11.3789 |
| 10 | **9.35e-02** | **8.30e-02** | 1.0e-02 | 5.0e-03 | 0.00e+00 | 0.00e+00 | 0 | 8.07e-08 / 7.60e-08 | 1.5636 / 1.5539 | 11.1995 / 11.3789 |

Where two values are printed, the first is the held run and the second the
halving one. Passes 1 to 8 are the same run: the control's first halving
fires after pass 8, and until then the key changes nothing. Pass 10 takes
no update in either run, so its row is the certification of the state pass 9
produced.

### 4.2 What the held bound did

**The row does not decay at a held bound of 1.0e-2.** MEASURED: with the
bound held the row RISES monotonically over all ten passes, 7.37e-2 to
9.35e-2, a factor 1.0268 a pass, +2.68 per cent. L33's rate proportionality
was measured on other cases and its 8 per cent a pass prediction is for
`kzz1e9/HeH0.083` (section 5), not for this one; what this case shows is
that its own bound, binding at every pass, does not set a decay rate here
because there is no decay.

**The one pass in which the two runs differ says the opposite of the
throttle reading.** Across pass 9's update the held bound moved the
composition by 1.65e-2 and the row went 9.27e-2 to 9.35e-2, +0.86 per cent;
the halved bound moved it by 8.26e-3, exactly half, and the row went 9.27e-2
to 8.30e-2, -10.46 per cent. Shortening the step made the row FALL and the
full step made it rise. One pass is one pass and not a rate, but it is the
opposite sign to the hypothesis that the halving is throttling a descent.

**Nothing in the run has reached a fixed point.** Every relaxation of both
runs ended on the movement bound with the bound attained, so the composition
handed to each wind solve is not the fixed point of the carrier operator at
that wind, which is L33's single distinguishing fact (READ, its section 2).
The H2 column is still moving: the x2 = 1e-2 contour walks inward
monotonically from 18.8174 to 11.7471 R_p over the eight shared passes and
the x2 = 0.5 contour outward from 1.4044 to 1.5443 R_p, continuing the
motion L34c recorded over its eighty passes. The wind under it stays
converged, `hydro info = 0` at every pass with the mass row between 1.19e-08
and 1.15e-07 against a cell tolerance near 3e-08.

So on this case the bound IS binding, and holding it at 1.0e-2 does not
recover a descent: the row rises at the held bound, and the only pass in
which the halving acted is the pass in which the row fell. The bound's
monotone halving is not what holds this case either. What the case is doing
is what L34c measured over twice as many passes: a molecular column still in
motion whose outermost cells carry a rising carrier balance, with the
composition never reaching the fixed point of its own operator.

---

## 5. `molecular_scalar_gj1132_kzz1e9/HeH0.083`: the case L33 section 6.3 designed the experiment for

### 5.1 Which state L33 meant, and which state the catalog publishes

L33's `i4_kz0083` table is `LHS1140b/models/.L22/i4_kz0083/run.log`, forty
passes (READ, its section "Read, unchanged"; `run_lart3.log` beside it is the
first five passes of the same run). MEASURED: the catalog case's published
`output/Hydro_ioniz.txt` has md5 `1229b961323d7ac710146650da72454a`, which is
byte for byte `.L22/i4_kz0083/output/Hydro_ioniz_IC.txt`, the **seed** of
that run. **The catalog's published pair is the state L33's forty passes
STARTED from, not the state they left.** The state they left is
`.L22/i4_kz0083/output/{Hydro_ioniz,Ion_species}.txt` (md5
`1ee5430a3c8560c40f017bc5eb618a8b`), whose final certification reads
`carrier balance H2` gated 8.188e-04 at cell 499 and `elemental transport
He/H partition` 2.515e-04 at cell 280, with the movement bound at 1.3e-3 in
its last passes (READ, that log). That is the state section 6.3 names, and it
is the one used here.

Its ending, READ from the same log: "the stationary solve returned info = 1;
output written", written `certified=F` with the reason "no stationary claim".

### 5.2 What was run, and how far it got

Configuration: `.L22/i4_kz0083/input.inp` unchanged except
`Coupled carrier solve: On stall` set to `False`, so that both runs stay on
the alternation route for every pass (the halving of the bound lives in the
branch guarded by `.not. block_now`, so a handover would end the comparison).
L33 records that the handover never fired in that run (READ, its section 7).
The case carries `He_diffusion: True` and `He_Kzz: 1.0e9`, so `comp_omega` is
also under the progress control; MEASURED, it fell to 0.250 at pass 3 in BOTH
runs, once each, so the key changes the bound and nothing else. Route:
`Restart intent: stationary`, `Well balanced: True` (resolved default),
`Secondary_ionization: Immediate`, `EXHALE_PTC_DTAU0=1.0`,
`EXHALE_OUTER_PASSES=10`, eight threads, `EXHALE_bnd.x`.

**A restart re-enters with `trust_pass = carrier_trust` = 1.0e-2 whatever the
state's history**, so the control here is the default restart: it also begins
at 1.0e-2 and halves only on a pass that counts as no progress. The two runs
are therefore identical until the control's first halving.

**The ten passes were not reached.** MEASURED: passes 1, 2 and 3 completed in
both runs and pass 4 in the run whose bound halves; both runs were then in a
line search that had collapsed (`lam` 9.5e-7 and 1.9e-6, `dtau` 5.6e-5) at
JFNK iterations 649 and 177 of a 3000 cap, at a cost that had grown from
7 s to 435 s to 1845 s to 3660 s a pass. At that rate the remaining passes
would have taken hours each, and the wind solve they were spending it on was
no longer converging. Both were stopped at 21:21 KST. What is reported below
is what ran; nothing is extrapolated.

### 5.3 The pass table

The worst gated row is `carrier balance H2` throughout, at cell 499 on entry
and cell 500 afterwards, r = 28.5489 and 29.0312 R_p. The x2 = 0.5 contour
stands at 1.3705 to 1.3905 R_p and the x2 = 1e-2 contour never enters the
domain, as L33 recorded for this case (READ).

| pass | H2 row, held | H2 row, halving | cell | r [R_p] | trust, held | trust, halving | displacement, held | displacement, halving | hydro info | worst mass row | x2=0.5 [R_p] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| entry | 8.188e-04 | 8.188e-04 | 499 | 28.5489 | 1.3e-03 | 1.3e-03 | - | - | - | - | 1.3770 |
| 1 | 8.15e-04 | 8.15e-04 | 499 | 28.5489 | 1.0e-02 | 1.0e-02 | 2.74e-02 | 2.74e-02 | 0 | 1.68e-08 | 1.3837 |
| 2 | 6.29e-03 | 6.29e-03 | 500 | 29.0312 | 1.0e-02 | 1.0e-02 | 2.80e-02 | 2.80e-02 | 0 | 4.34e-09 | 1.3837 |
| 3 | 1.93e-03 | 1.93e-03 | 500 | 29.0312 | 1.0e-02 | 5.0e-03 | 2.78e-02 | 1.41e-02 | 2 | 1.65e-01 | 1.3705 / 1.3770 |
| 4 | not reached | 1.28e-02 | 500 | 29.0312 | 1.0e-02 | 5.0e-03 | not reached | 1.16e-02 | 0 | 2.14e-08 | 1.3905 |

The entry row is READ from `.L22/i4_kz0083/run.log`; every other number is
MEASURED here. Where two values are printed, the first is the held run.
Passes 1 to 3 are the same run in both: the control's first halving fires in
pass 3's progress control, so the bound differs from pass 3's UPDATE onward
and the printed row of pass 3 is still common. The halved update is exactly
half the held one, 1.41e-2 against 2.78e-2.

**The bound binds at every pass of both runs.** All relaxations ended on the
movement bound with it attained (measure 1.000E-02, then 5.000E-03 in the
control), at cells 422, 426, 419 and 227, that is 8.197, 8.717, 7.830 and
1.239 R_p.

### 5.4 The reading L33 asked for

L33's test: about 8 per cent a pass of decay at a held 1.0e-2 means the
throttle is the whole story; about 1 per cent whatever the bound means the
row is a mode (READ, its section 6.3).

**Neither. At a bound of 1.0e-2 there is no decay whose rate could be
measured.** MEASURED here from the state L33 left: 8.15e-04, then 6.29e-03,
a factor 7.72 UP in one pass, then 1.93e-03, a factor 0.307 down, with the
bound held at 1.0e-2 throughout. The row is an order of magnitude above where
the small bound had left it after a single pass at the restored bound.

**L33's own log says the same thing where it ran at 1.0e-2**, which is the
strongest available check and is READ, not inferred: its passes 1 to 5 stood
at 2.97e-03, 2.52e-02, 7.98e-03, 3.35e-03 and 4.10e-03, factors 8.48, 0.317,
0.420 and 1.224, while the bound was 1.0e-2. The monotone decay L33 measured,
2.49 per cent a pass and then 1.04, begins only at passes 10 and 23, at
bounds of 2.5e-3 and 1.3e-3.

So the ordering L33 read as a throttle is real but its direction is the other
one: **the monotone descent is a property of the SMALL bound, not a slowed
version of a faster descent that a larger bound would restore.** Restoring the
bound to 1.0e-2 does not resume a descent at 8 per cent a pass; it puts the
row back into a band an order of magnitude wide.

**And it costs the wind.** MEASURED: pass 3 returned `info = 2` in both runs
with the hydrodynamic mass row at 1.65e-01 against a cell tolerance near
3e-08; the control recovered at pass 4 (`info = 0`, mass 2.14e-08) at the
halved bound, and the held run did not reach a pass 4 at all, its JFNK
standing at `||R|| = 2.225e-01` with `dtau` collapsed to 5.64e-05 when it was
stopped, against `||R|| = 1.010e-01` and the same `dtau` in the control one
pass further on. The composition step the held bound allows, 2.78e-2 against
1.41e-2, is what the wind solve then has to absorb.

The claim this measurement does NOT support is the one L33 offered as the
alternative: nothing here shows the row to be a mode the alternation cannot
damp. What it shows is that the bound at 1.0e-2 is too large for this state,
not too small.

---

## 6. The test row

`src/tests/grid_and_gates/carrier_movement_bound_hold.sh`, wired into that
suite's `run.sh` and README as `carrier_bound_hold`. Five assertions on
scratch copies of `backup/regression/carrier_model_a_newton` entered from
its pinned IC pair: one outer pass with the key unset against one with it
empty, and fourteen outer passes with the key unset against fourteen with it
held at 1.0e-2.

MEASURED, the measured build, one thread:

```
PASS carrier_bound_absent_is_inert measured=bitwise_equal reference=bitwise_equal tol=0
PASS carrier_bound_absent_is_silent measured=0 reference=0 tol=0
  bound by pass, control: 1.0E-02 x13 5.0E-03
  bound by pass, held:    1.0E-02 x14
PASS carrier_bound_default_halves measured=1 reference=ge_1 tol=0
PASS carrier_bound_key_holds measured=held_every_pass reference=held_every_pass tol=0
PASS carrier_bound_passes_before_agree measured=13_passes_equal reference=equal tol=0
```

MEASURED, the control build of the entry text, the same command: the same
five rows with `carrier_bound_key_holds` FAILING at
`passes_off_the_value=1,halvings=1`, since without the key the held run
halves with the control at pass 14. The other four are green there, which is
what "the key's absence changes nothing" means in the suite.

**Changed 2026-09-19.** `carrier_bound_default_halves` is no longer an
assertion. Whether the control halves within fourteen passes turns on the
mass row moving at its own rounding (MEASURED after item D9fix: 3.13e-11 to
7.16e-11 on one build, 4.06e-11 to 4.03e-11 on the next, so one halved at
pass 14 and the other never did), so it is printed as a DIAGNOSTIC and the
row `carrier_bound_hold_is_announced` (the held run states the held bound)
takes its place. Where the control does not halve, `carrier_bound_key_holds`
cannot tell the key from its absence, and the announcement row is what
shows the key reached the loop.

The entry costs about four minutes at one thread with the two fourteen-pass
runs in parallel, which is stated in its header and in the README row;
`EXHALE_BOUND_HOLD_PASSES` moves the pass count for a fixture whose first
halving falls elsewhere.

---

## 7. Noticed outside scope, reported and not fixed

1. **The case directory has no `ENDING` record.** The brief asked for the
   published generation's `ENDING`;
   `LHS1140b/models/molecular_scalar_gj1132_wellmixed/HeH0.083/` holds
   `not_solved.md`, `REPRODUCE.md` and the two logs, and the ending exists
   only as the last two lines of `run.log`. That is D8 item 1's own subject
   (the ending reason of each generation recorded) and is not fixed here.

2. **The certification measure and the distance from the root move at
   different rates on a stiff stage.** On `carrier_model_a_newton` the He II
   measure at cell 248 fell 1.87 per cent a pass while the abundance
   distance behind it fell 7.52 per cent a pass and changed sign, because
   the scale of a stiff stage is its own net source and that source shrank
   with the residual. A reader watching the gated measure alone would judge
   this state to be nearly still. This is a reading of D3a's conditioning
   statement in a trajectory rather than at one state, and it belongs to
   D3b, not here.

3. **None of the three is a state the acceptance can certify yet, for
   reasons that have nothing to do with the bound.** `carrier_model_a_newton`
   is closing on the root at 7.5 per cent a pass from 1.58e-3 of the stage
   fraction; `wellmixed/HeH0.083` has a molecular column still in motion
   whose outer carrier balance rises; `kzz1e9/HeH0.083` loses its wind solve
   within three passes of the bound being restored. All three are recorded
   here as measurements; no item is opened on the movement bound.

4. **A catalog case publishes its own seed as its state.**
   `molecular_scalar_gj1132_kzz1e9/HeH0.083`'s published
   `output/Hydro_ioniz.txt` is byte for byte the `_IC` file of the run that
   produced it, so the forty passes of `.L22/i4_kz0083` are not represented
   in the catalog at all and a reader taking "the state the run left" from
   the case directory gets the state it started from. That is D8 items 1, 3
   and 4 in one example, and it is reported, not fixed.

5. **`src/EXHALE_main.f90` was edited by another item while this one held
   it.** MEASURED against the entry snapshot: eight hunks that are not this
   item's stand in the file at lines 89, 1695, 2236, 4212, 5660, 5788, 5792
   and 6948, beginning with a `report_base_boundary_model` import from
   `base_boundary`, written at 19:35. This item's three hunks are intact and
   the other eight are untouched by it. Reported, not acted on; it is why two
   measured builds exist.

6. **A JFNK line-search collapse ends these solves without ending the pass.**
   Both `kzz1e9` runs spent hours at `lam` near 1e-6 and `dtau` near 5.6e-5
   with the outer iteration cap at 3000, so a pass that cannot converge costs
   the full cap before the outer loop can judge it. The cost grew 7 s, 435 s,
   1845 s, 3660 s over four passes. Nothing here diagnoses it; it is recorded
   because it is what stopped this experiment.
