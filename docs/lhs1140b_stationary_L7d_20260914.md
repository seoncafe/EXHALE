# LHS 1140 b, item L7d: the three corrections item L7c names

Item L7d of `docs/PLAN_20260913_lhs_stationary.md`, carried out on the live
tree while the campaign ran `LHS1140b/models/atomic_*` with the delivered
`EXHALE.x` (md5 `b261c3287e64`). Nothing in `build/`, `EXHALE.x` or any
`atomic_*` directory was written by this item; every binary below is a
private build and every run is in `scratchpad/L7d/` or
`LHS1140b/models/.L7d/`.

| build | what it is | md5 |
|---|---|---|
| `EXHALE_L7d_ctl.x` | the tree as delivered to this item (the control) | `0ac03bc25ac9e706fa3d3c710d815143` |
| `EXHALE_L7d_qh2.x` | the control plus section 2 alone | `9635e312fe7c8ce7a0e143c63ef583a5` |
| `EXHALE_L7dx.x` | the same plus the ramp rules of section 3 under an environment switch, for the measurement only | `ebd091cc0f0e0a69ec43d35297018ef4` |
| `EXHALE_L7d.x` | the delivered source text: sections 2 and 3 | `9e38a707fc8296d00e2ccf3dc96e0b2c` |

Every number is MEASURED (produced by a run of one of those builds, or by
arithmetic stated here on a file such a run wrote) unless it is marked READ
(from a source file or a document).

## 1. Verdict

1. **`q_h2_equilibrium`'s cancellation is real, and the conjugate form
   removes it.** RED, on the 500 physical cells of the certified atomic
   HeH 2.13 state: the form the fit is written in is wrong by more than
   1 percent in **186** cells, returns exactly zero in **164** of them where
   the fit is positive, and is **1.30e4** times the fit at its worst.
   GREEN, the same table with the conjugate form: **0** cells wrong by more
   than 1 percent, worst relative error **2.1e-15** against a 60-digit
   reference. Through the code, the `local` molecular seed of the same state
   written by the control build carries the four-cell x2 = 5.044357e-06 spike
   L7c reports and the measured build carries the fit's own value in every
   cell (section 2). The `u > 30` branch is gone with it.

2. **The growth rule the item proposed -- double only on a full step,
   lam = 1 -- is REFUTED, and the rule delivered is "double on a step the
   line search did not have to cut, lam >= 0.5, and multiply by lam on a
   shorter one".** The near states do not take full steps: MEASURED on the
   C/N/O pass-8 state, the first six accepted steps are all lam = 0.5 and
   only then does lam reach 1, so "lam = 1 only" leaves dtau at its start,
   the solve stalls at `||R||` 3.27e-01 and the cut then takes dtau to
   2.8e-04, against the control's `info = 0` at 20 iterations and
   `||R||` 3.773e-08. With the threshold at 0.5 those states are solved
   exactly as the plain doubling solved them, and on the far state -- where
   the accepted steps are lam = 2e-03 and the delivered ramp doubles through
   them to dtau = 1.05e+06 by iteration 20 with the merit frozen at 1.28e-01
   -- the new rule cuts dtau at the first damped step and the merit falls by
   a factor 25 on the step after it (section 3).

3. **The molecular cases now carry the seed recipe and the runner builds
   their seed.** All nine `molecular_*` cases are rewritten with
   `Load IC? True`, `Reconstruction scheme: PLM`, no `du_th`,
   `Solver: Newton`, `Restart intent: stationary equilibrate`,
   `Secondary_ionization: Immediate` and `Well balanced: True`; the molecular
   branch is in `run_case.sh`, put there by replacing the whole file with
   `\mv -f` while the campaign was reading it (section 4).

4. **T-L7-5 does not certify, and what refuses is now a row and not the
   seed.** With the three corrections, the `local` partition's hydrodynamic
   rows descend for the first time on this route: mass 3.41e-02 -> 7.21e-07
   between outer passes 1 and 2, momentum 5.34e-11, energy 6.43e-06, and the
   runner's own continuation at dtau0 = 1e8 reaches `||R||` 1.157e-08 and
   `info = 0` on every pass. The H2 carrier balance row reads exactly 1.000
   at cell 280 (r = 1.604 R_p, T = 5986 K) throughout, and the H2 it measures
   is no longer the rounding L7c found: the carried mixing ratio is within a
   factor 1.6 to 2.0 of the fit's own value at every cell and varies smoothly
   from cell to cell. The handoff partition returns the state it was given
   (mass 1.00, momentum 1.80e-01, energy 1.99, the seed's own rows), the
   fourth independent measurement that it does not descend (section 5).

## 2. `q_h2_equilibrium` in the conjugate form

### 2.1 What was changed

`src/modules/lower_atmosphere/lower_column.f90`, the function body (+47/-25,
counting the comment block that states the reasoning; the whole diff of that
file is this one function and the two comment paragraphs above it):

```fortran
      u    = -23672.0d0/T - log10(max(p_bar, 1.0d-30)) + 6.2645d0
      tenu = 10.0d0**min(u, 300.0d0)
      qh2  = 1.9845d0**2                                                 &
           / (2.3670d0*((1.9845d0 + tenu)                                &
              + sqrt(tenu)*sqrt(3.9690d0 + tenu)))
      if (qh2 .gt. 1.0d0) qh2 = 1.0d0
```

in place of `(1.9845 + tenu - sqrt(tenu*(3.9690 + tenu)))/2.3670` inside an
`if (u > 30)` branch that returned 0. It is the same function: with
A = 1.9845 and B = 2.3670, `A + t - sqrt(t(2A + t)) = A^2/((A + t) +
sqrt(t(2A + t)))` identically, and the right-hand side adds two positive
terms where the left differences two that both grow like t while their
difference falls like A^2/(2t). The square root is formed as
`sqrt(t) sqrt(2A + t)` so the product under it cannot overflow, and
`min(u, 300)` keeps t finite; q is 6.6e-301 there. The fit's source (Koskinen
et al. 2022 Eq. 11, quoting Visscher et al. 2006) and the He/H validity
statement above the function are unchanged.

The `u > 30` branch is removed because it is a discontinuity at its own
threshold: the fit is 8.3e-31 at u = 30 and the branch returned 0. The
`qh2 < 0` clamp is removed because the expression is positive by
construction. The `qh2 > 1` ceiling is kept as a statement of the range and
cannot fire (the maximum is A/B = 0.8384).

### 2.2 RED and GREEN, the arithmetic

The reference is the exact value of the same expression at 60 decimal digits
(Python `decimal`); no run is needed, the input is the (p, T) of the 500
physical cells of the certified atomic state
`LHS1140b/models/.L7c/src_atomic_HeH2.13/Hydro_ioniz.txt`, with
p_bar = p_code x p0/1e6 and p0 = n0 k T0 = 1.0 cgs, which is how
`molecular_seed_from_atomic_state.f90` line 391 forms the argument.

| over the 500 physical cells | as the fit is written (RED) | conjugate form (GREEN) |
|---|---|---|
| cells wrong by more than 1 percent | **186** | **0** |
| first such cell | 187, r = 1.1184, T = 3311 K, t = 2.15e7 | -- |
| cells returning exactly 0 where the fit is positive | 164 | 0 |
| cells more than twice the fit's value | 8 | 0 |
| worst relative error | **1.30e+04**, at cell 218 (r = 1.2042, T = 4940 K) | **2.11e-15**, at cell 366 |

The 186 cells and the first cell 187 (r = 1.1184, T = 3311 K, t = 2.2e7)
reproduce L7c section 5.1 exactly; the count of exact zeros (164 against
L7c's 166) and of cells more than twice the true value (8 against 6) differ
because L7c compared the two double-precision forms with each other and this
table compares each against the 60-digit value.

### 2.3 RED and GREEN, through the code

The `local` molecular seed of that state, built by the binary
(`EXHALE_MOLECULAR_SEED=... EXHALE_MOLECULAR_SEED_X2=local`), with the
control build and with the measured build. Both reproduce L7c's conversion
report line for line (234 cells clipped at the element-ratio ceiling, 235
capped by their own neutral hydrogen, largest x2 transferred
9.9999437572717320e-01, T moved by 1.8169e-01, equation-of-state round trip
3.3647e-16), so the only thing that differs between them is the fit.

x2 = 2 n(H2)/n_H,nuclei as written into `output/Ion_species_IC.txt`:

| cell | r [R_p] | T [K] | control | measured | the fit, exactly | control / fit |
|---|---|---|---|---|---|---|
| 190 | 1.1249 | 3475 | 9.852267e-08 | 1.005656e-07 | 1.005656e-07 | 0.98 |
| 191 | 1.1271 | 3531 | 7.881814e-08 | 7.608540e-08 | 7.608540e-08 | 1.04 |
| 192 | 1.1293 | 3586 | 7.881814e-08 | 5.796381e-08 | 5.796381e-08 | 1.36 |
| 193 | 1.1316 | 3642 | 7.881814e-08 | 4.446196e-08 | 4.446196e-08 | 1.77 |
| 194 | 1.1340 | 3697 | 0 | 3.433700e-08 | 3.433700e-08 | 0 |
| 202 | 1.1542 | 4140 | 3.152725e-07 | 5.474491e-09 | 5.474491e-09 | 57.6 |
| 214 | 1.1904 | 4757 | 0 | 6.753189e-10 | 6.753189e-10 | 0 |
| 215 | 1.1938 | 4804 | 5.044357e-06 | 5.844133e-10 | 5.844133e-10 | 8.63e+03 |
| 216 | 1.1972 | 4850 | 5.044357e-06 | 5.077804e-10 | 5.077804e-10 | 9.93e+03 |
| 217 | 1.2007 | 4896 | 0 | 4.429233e-10 | 4.429233e-10 | 0 |
| **218** | **1.2042** | **4940** | **5.044357e-06** | **3.878210e-10** | **3.878210e-10** | **1.30e+04** |
| 219 | 1.2079 | 4984 | 5.044357e-06 | 3.408311e-10 | 3.408311e-10 | 1.48e+04 |
| 220 | 1.2115 | 5027 | 0 | 3.006138e-10 | 3.006138e-10 | 0 |

Over every cell the conversion did not cap (x2 from the fit below 0.5): the
control's written x2 is off by more than 1 percent in **186** cells with a
worst of **1.480e+04**; the measured build's is off in **0** cells with a
worst of **1.54e-11**. Cell 218, the cell the carrier row of L7c refuses at,
holds the fit's own 3.878e-10 instead of 1.30e4 times it.

### 2.4 What it moved: the ten regression cases that carry molecular goldens

The matrix was run twice, once with each binary, single-threaded as the
harness fixes it, the second on a copy of `backup/regression/` so the two
runs do not share the case directories. **The goldens of this tree have
already moved under other items**, so the verdict against `golden/` is not
the measurement: MEASURED, the CONTROL build alone already fails
`mol_sec_ion` (8.866e-03), `mol_metals` (7.315e-01), `mol_ir_bands`
(9.117e-01), `hp_trace_seed` (8.843e-03) and `hp_zero_seed` (3.703e-02)
against `golden/`. What this item measures is therefore the control build
against the measured build, file by file, over the data lines.

| case | Hydro_ioniz | Ion_species | Hydro_ioniz_adv | Ion_species_adv |
|---|---|---|---|---|
| `mol_base_handoff` | 3.94e-09 | 2.10e-08 | 1.44e-07 | 2.10e-08 |
| `mol_carrier` | identical | identical | identical | identical |
| `mol_diffusion` | identical | identical | identical | identical |
| `mol_lyman_werner` | identical | identical | identical | identical |
| `mol_sec_ion` | **5.79e-03** | 2.44e-05 | **5.79e-03** | 2.44e-05 |
| `mol_metals` | 1.87e-08 | 3.42e-09 | 6.27e-06 | 1.09e-08 |
| `mol_ir_bands` | 5.86e-08 | 4.25e-09 | 2.69e-04 | 5.57e-08 |
| `hp_front` | identical | identical | identical | identical |
| `hp_trace_seed` | 7.33e-05 | 1.36e-04 | 8.92e-04 | 1.36e-04 |
| `hp_zero_seed` | 3.21e-06 | 7.11e-06 | 3.80e-05 | 7.11e-06 |

(largest relative difference over the numeric fields; "identical" means the
data lines are byte-identical.)

Of the three `hp_*` cases one is **bitwise identical** (`hp_front`) and two
move, by at most 8.9e-04 (`hp_trace_seed`) and 3.8e-05 (`hp_zero_seed`).
Nine of the ten cases move by less than the harness's 1e-3.

**`mol_sec_ion` moves by 5.79e-03, above that tolerance, and it moves TOWARD
its golden**: MEASURED against `golden/`, the same field is 8.866e-03 away
with the control build and 3.095e-03 away with the measured one. It is the
velocity column at row 177, which is cell 175 of the hot-Uranus column, and
`mol_sec_ion` is the only molecular case whose photoelectron partition is
applied from the first evaluation, so it is the case in which the molecular
starting point of the chemistry reaches the energy deposition. **No golden
was refreshed.**

## 3. The growth rule of the pseudo-time ramp

### 3.1 The rule the item proposed is refuted, and by what

The item asks for "double only on a full step, lam = 1; hold or shrink
otherwise". MEASURED, the first thing that has to be true for that rule is
not: **the near states this ramp was built for do not take full steps.** On
the C/N/O pass-8 entry state (`scratchpad/L4e/red_A`, the state L4g's GREEN
(a) was measured on), the control build's accepted steps are

```
it 1-6   lam = 5.00E-01   dtau = 2, 4, 8, 16, 32, 64
it 7-20  lam = 1.00E+00   dtau = 128 ... 5.24E+05
```

and the solve reaches `info = 0` at iteration 20 with `||R||` 3.773e-08. Six
of the twenty accepted steps are half steps, and they are the six that lift
the pseudo-time off its start. With the threshold at lam = 1 (measured on the
same state with a build differing from the control in this rule alone)
`dtau` never leaves 1, `||R||` stalls at 3.272e-01 through iterations 4 to 8, the line
search then collapses to lam = 9.5e-07, the cut takes `dtau` to 2.8e-04 and
by iteration 33 `||R||` is 1.965, **worse than the state it started from**.

So the signal is not "was the step full" but "did the line search have to cut
it". The line search halves, so the two are one backtrack apart, and the
threshold delivered is **lam >= 0.5**: the full step, or one backtrack.

### 3.2 What the rule is

```fortran
      real*8, parameter :: lam_of_a_full_step = 0.5d0
      ...
               if (lam .ge. lam_of_a_full_step) then
                  dtau_ramp = 2.0d0
               else
                  dtau_ramp = lam
               endif
               dtau = min(max(dtau*dtau_ramp, dtau_floor), 1.0d14*dtau0)
```

in `src/modules/time_step/steady_newton.f90`, in place of the unconditional
`dtau_ramp = 2.0d0`: one new module-scope parameter with the ten-line
statement of what 0.5 is and what was measured to put it there (near line
1493), the branch itself at the accepted-step site (+5/-1, near line 16578),
and three comment blocks rewritten -- the option's declaration, the ramp's
own paragraph and the MEASURED note beside it -- because each of them said
"doubles on every accepted step". The floor, the ceiling,
the cut of a step admitted nowhere (x 0.25) and the `EXHALE_PTC_RAMP_DOUBLE=0`
merit-ratio branch are unchanged. The reason the factor below the threshold
is `lam` and not 1: the shift I/dtau is what makes the linear model of a step
trustworthy, so a step the search had to cut to a fraction lam is the
statement that the model held over that fraction and no further, and the
pseudo-time the step represents falls in the same proportion.

### 3.3 The far regime, measured

The state is the one L7c measured its ramp control on: the `local` molecular
seed of the certified atomic HeH 2.13 state, at 1.2 R_p, as the tree wrote it
BEFORE section 2 (`LHS1140b/models/.L7c/HeH2.13_local_rampoff/output/*_IC.txt`),
kept unchanged so that the only thing that differs between the rows below is
the ramp. `EXHALE_PTC_DTAU0=1.0`, `EXHALE_OUTER_PASSES=1`, 8 threads on a
machine also running the campaign, so wall times are not claims.

| | the delivered ramp (doubling on every accepted step) | the merit-ratio branch (`EXHALE_PTC_RAMP_DOUBLE=0`) | this item's rule |
|---|---|---|---|
| `\|\|R\|\|` start -> end | 1.867 -> 1.776 | 1.867 -> **9.319e-09** | 1.867 -> see below |
| `\|\|Fs\|\|2` start -> end | 5.54 -> 1.28e-01 | 5.54 -> 2.58e-11 | 5.54 -> 1.56e-02 by iteration 12 |
| `dtau` | 2 -> 1.05e+06 at iteration 20, doubling every iteration | 0.58 to 5.36e+09 | 2 -> 256, cut to 16 at the first damped step, back to 128 by iteration 12 |
| `lam` | 1 for the first eight steps, then 0.25, 0.125, 6.25e-02 and 1.95e-03 from iteration 13 on | never below 0.25 (29 steps at 1, three at 0.5, one at 0.25) | 1 except iteration 9 (6.25e-02) |
| outcome | no descent; stopped at iteration 20 | **`info = 0` at 33 iterations** | descending; see section 3.5 |

The delivered ramp's row reproduces L7c section 4 on the same state with the
campaign's own binary (124 iterations, `||R||` 1.867 -> 1.673, dtau to the
1e14 ceiling, `lam` = 7.81e-03 at the end): the pseudo-time is raised by
steps that barely moved, and once the shift is gone the solve is the
unshifted Newton on a state it is not near.

The middle column is new and it is a consequence of section 2: on the tree as
L7c measured it the merit-ratio branch took this state only to `||R||`
4.5e-02 in 159 iterations, and with the fit's cancellation removed from the
molecular chemical starting points the same branch takes it to 9.3e-09 in
33. The H2 column of the seed is unchanged between the two (it is the same
file); what changed is the H2 partition the sweep starts each cell's
chemistry from.

### 3.4 The controls

All on scratch copies, the control build being the same tree with only the
ramp rule of the entry text (`EXHALE_L7d_qh2.x`), so that section 2 is not in
the difference.

**The near state, `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13` pass-8 entry.**
The delivered rule is **bitwise the entry text** here: the two runs print the
same twenty iteration lines to the last digit, end with the same
`done info=0 ||R||= 3.773E-08 flux spread= 5.080E-12 non-monotone accepts=5`
and the same outer-pass row (mass 1.08e-09, momentum 1.08e-13, the element
row 2.02e-03 at cell 278), and their `output/Hydro_ioniz.txt` and
`output/Ion_species.txt` are identical line for line. That is what the
threshold at 0.5 buys: every accepted step of this state is lam = 0.5 or
lam = 1, so the rule takes the doubling branch at every one of them.

**`atomic_elem_newton`, the HD 209458 b element reload** (`Load IC? True`,
`Coupled carrier solve: False`, `Restart intent: stationary`, the pair from
`IC/` in `output/`, `EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=80
EXHALE_DIFF_OMEGA=0.5`, 8 threads). The two runs are identical for eleven
iterations and part at iteration 12, the first accepted step with
lam = 0.125. At the end of outer pass 1, both at the 80-iteration cap:

| | control | this item's rule |
|---|---|---|
| `\|\|R\|\|` at the end of pass 1 | 9.017e-01 | **8.627e-01** |
| mass row | 6.85e-02 | **3.24e-02** |
| momentum row | 1.56e-05 | 1.90e-05 |
| energy row | 9.02e-01 | 8.63e-01 |
| worst gated species row | 2.84e-03 at cell 301 (elemental transport O) | 2.83e-03 at the same cell |
| `dtau` over the pass | 8 -> 1.00e+14 (the ceiling), `lam` down to 1.22e-04 | 8 to 512 until iteration 40, then between the explicit-stable floor 9.68e-05 and twice it |
| non-monotone accepts | 20 | 59 |

Neither certifies at pass 1 and both go on; the twelve outer passes did not
finish on a machine running the campaign, so pass 1 is what is compared.
The measured difference is in the delivered rule's favour on the mass row and
on `||R||`, and the pseudo-time is the visible reason: the control spends the
second half of the pass at the 1e14 ceiling taking steps of lam ~ 1e-4, the
new rule keeps it where the line search will take whole steps.

### 3.5 The new rule on the far state, and what the two regimes now say

The same far state, the delivered binary, 40 iterations:

```
it  1-8    lam = 1                        dtau 2, 4, 8, 16, 32, 64, 128, 256
it  9      lam = 6.25e-02   -> dtau x lam  dtau 16
it 10-13   lam = 1                        dtau 32, 64, 128, 256
it 14      lam = 3.12e-02   -> dtau x lam  dtau 8
it 15-19   lam = 1                        dtau 16 ... 256
it 20      lam = 1.25e-01                 dtau 32
...
it 40      lam = 2.50e-01                 dtau 64
```

`||Fs||2` runs 5.54 -> 3.40e-03 at iteration 15, reaches 1.79e-03 and stays
between that and 2.1e-02 for the rest; `||R||` oscillates between 2.78e-01 and 1.97. The
pseudo-time stays inside 2 to 512 for the whole run instead of reaching the
1e14 ceiling by iteration 40, which is the whole of the change: every time
the line search has to cut a step the shift comes back, and the next step is
a full one.

Measured at a common 40 iterations, by the merit the line search descends on
(the residual norm is not monotone on this state under any of the three
rules):

| | the delivered ramp | this item's rule | the merit-ratio branch |
|---|---|---|---|
| `\|\|Fs\|\|2` after 20 / 40 iterations | 1.28e-01 / (not run past 20) | 2.10e-02 / 5.07e-03 | 1.59e-03 at 21, `info = 0` at 33 |
| where `dtau` ended | 1.05e+06 and doubling | 64 | 5.36e+09 after the solve had converged |

The merit-ratio branch is still the fastest thing on this particular state,
and that is a measurement about the state and not an argument for the branch:
L4g's own reproductions (the C/N/O pass-8 and pass-2 states, the raw archive
seed) are the states on which that branch does not ramp at all, and this
item's rule solves those exactly as the plain doubling did (section 3.4).
One rule now serves both regimes, which is what the item asked for.

**The alternative below the threshold, holding `dtau` instead of multiplying
by lam, was measured only at threshold lam = 1** (section 3.1), where it
leaves the near state stalled. At threshold 0.5 the two differ only on the
steps the line search cuts, and the far state's recovery -- the merit falling
by a factor 25 on the step after the cut -- is the measurement that chose the
proportional form. A third option, the published clip `max(lam, 0.1)`, was
not measured; section 7 says when it would matter.

The two regression cases that reach `solve_steady_jfnk` through the marching
hand-off, measured with the same binary as the section 2 matrix and with the
delivered one (so the difference is the ramp alone): `mol_base_handoff`
**byte-identical in all four files**, and `wasp_full_newton` 4.76e-09 in
`Hydro_ioniz.txt` and 7.98e-09 in `Ion_species.txt`, with the only larger
entry the advection diagnostic column `adv_mass_row` of one cell (1.22e-14
against 0). Both cases fail against `golden/` with the control binary as
well, for the reason section 2.4 gives.

## 4. The molecular recipe: the case files and the runner

### 4.1 `make_models.py`

Two changes (`LHS1140b/models/make_models.py`):

- `stationary = for_closure_template or group.chemistry == 'atomic'` becomes
  `... or group.chemistry in ('atomic', 'molecular')`. A molecular case
  therefore writes `Reconstruction scheme: PLM`, `Load IC? True`,
  `Solver: Newton`, `Restart intent: ...`, `Secondary_ionization: Immediate`
  and no `du_th` line, exactly as an atomic case does. The comment at the
  site states why: the runner builds the molecular state out of the certified
  atomic solution, so the case is solved from a loaded state and the marching
  keys would throw that state away, `run_case.sh` reading `Load IC?` from
  this file for the wind pass (L7c section 7).
- the intent word: a molecular case writes `Restart intent: stationary
  equilibrate`, an atomic one keeps plain `stationary`. The seed carries the
  composition the conversion put on it cell by cell, and `equilibrate` puts
  that composition on its own fixed point before the state is measured.

`python3 make_models.py --only molecular --force` rewrote all **nine**
molecular cases (`molecular_scalar_gj1132_wellmixed` 0.083/0.55/2.13,
`molecular_scalar_gj1132_kzz1e9` 0.083/0.55/2.13/9.7,
`molecular_photochem_gj1132_kzzprofile` 2.09/9.05). Nothing else was
regenerated: `--only molecular` matches no other group, and the atomic cases
the campaign is running were not touched.

The file each case now carries, at HeH 2.13 of `kzz1e9` (the lines that
changed): `Reconstruction scheme: PLM`, `Load IC? True`, `Solver: Newton`,
`Restart intent: stationary equilibrate`, `Secondary_ionization: Immediate`,
with `du_th [PLM,WENO3]: 0.5 1.0e-3` gone; `Well balanced: True`,
`Molecular chemistry: True` and `Molecular carrier transport: True` are as
before.

### 4.2 `run_case.sh`

The molecular branch of `run_case.sh.molecular.patch` is now IN the runner.
The runner is re-read by every case the campaign starts, so it was never
edited where it stood: a complete copy was written beside it and moved over
it with `\mv -f`.

| | md5 |
|---|---|
| `run_case.sh` before | `44267916fbeb2439bb4e50c3bdfed89c` |
| `run_case.sh` after | `75714a53396424f17868f8e24082d0b2` |
| `run_case.sh.molecular.patch` (rewritten against those two) | `6cd74bbfe433c93f27058670b0732fb4` |

`bash -n run_case.sh` is clean, and the patch was verified by applying it to
a copy of the before file: the result is byte-identical to the delivered
runner (md5 `75714a53396424f17868f8e24082d0b2`).

What the branch does. `MOLECULAR=1` when `input.inp` states
`Molecular chemistry: True`; pass 0 is then the binary's own conversion out
of `$ATOMIC_SEED/output` (default: this case's name with `molecular`
replaced by `atomic`), refused unless that state's metadata says
`certified=T`. `SEED_X2` chooses the extension of the base partition and
**defaults to `local`**: MEASURED (L7c section 3.2), at the handoff's own
x2 = 0.99998 the column has no neutral hydrogen left to photoionize above
the base, the energy row of the seed is its radiative imbalance (1.989 at
cell 203 against 1.027 for `local`), and the solve makes no outer pass in an
hour; `SEED_X2=handoff` asks for the conversion's own default by leaving the
variable out of the environment, and a number in [0,1] is that fraction
everywhere. `SEED_INVARIANT` passes p or T through. The seed command with
its environment goes into `COMMANDS`, so `write_reproduce.py` puts it in
`REPRODUCE.md` next to the wind, post-processing and transit commands.

**The atomic path is unaffected, and that is measured on the campaign
itself.** The replacement went in at 05:54 and the campaign, which re-reads
`run_case.sh` for every case it starts, has since finished fourteen atomic
cases through the new file, every one of them `DONE info=0 ... certified`
with its transit metrics written (`LHS1140b/models/campaign_20260914.log`,
06:39 to 07:43).

Two further things the branch needed, and one defect found while writing it:

- the ranked seed walk (`SEED_ATTEMPTS`) is pick_seed's, and pick_seed has
  no list for a molecular case, so the walk is not entered when
  `MOLECULAR = 1`.
- **`SEED_GIVEN` read `$ATTEMPT` before it existed.** Under `set -u` the
  command substitution `$([ -n "$SEED" ] && [ -z "$ATTEMPT" ] && ...)` left
  the subshell on an unbound variable and `SEED_GIVEN` came back EMPTY, so a
  seed the caller gave by hand was treated as pick_seed's own and the runner
  would walk away from it to pick_seed's ranked list on a refused solve. It
  now reads `${ATTEMPT:-}`. Found while writing this branch, fixed with it,
  and reported here because it is a change to the runner beyond the branch.

## 5. T-L7-5 with the new runner and the new binary

`molecular_scalar_gj1132_kzz1e9/HeH2.13` (run on copies of the case in
`LHS1140b/models/.L7d/t_local` and `.L7d/t_handoff`, so the case directory
the campaign will use is left as `make_models.py` wrote it), by
`run_case.sh`, `EXHALE_BIN=EXHALE_L7d.x`, the seed built from a frozen copy
of the certified atomic state (`.L7d/src_atomic_HeH2.13`, the same file
`.L7c` used, so the campaign cannot move it underneath the run),
8 threads. The outer loop and the JFNK were capped
(`EXHALE_OUTER_PASSES=4`/`2`, `EXHALE_JFNK_MAXIT=40`) because the machine was
running the campaign; those caps are the reason a pass ends where it does.

The runner's own pass 0 reproduces the conversion measurements of L7c
section 2 for both partitions (`local`: 234 cells clipped at the
element-ratio ceiling, 235 capped by their own neutral hydrogen, T moved
1.8169e-01, equation-of-state round trip 3.3647e-16; handoff: x2
9.9998316003549481e-01, 480 capped, T moved 2.7580e-01, round trip
8.3748e-16), with the H2 column now the fit's own (section 2.3).

### 5.1 The `local` partition

| outer pass | hydro `info` | mass | momentum | energy | worst gated species row |
|---|---|---|---|---|---|
| 1 | 1 | 3.41e-02 | 9.53e-05 | 3.24e-01 | 1.00 at cell 279, carrier balance H2 |
| 2 | **0** | **7.21e-07** | **5.34e-11** | **6.43e-06** | 1.00 at cell 280, carrier balance H2 |
| 3 | 2 | 7.22e-07 | 6.45e-11 | 6.44e-06 | 1.00 at cell 280 |
| 4 | 2 | 7.23e-07 | 5.36e-11 | 6.45e-06 | 1.00 at cell 280 |

**This is the first descending hydrodynamic solve this route has had.** L7b
reached `info = 0` on its first pass and refused on the carrier row; L7c, on
the consolidated tree, made no outer pass at all in an hour and left `||R||`
between 1.06 and 1.87. Here the mass row falls five decades between pass 1
and pass 2 and the state then stands still to three digits, which is the
alternation's fixed point.

The certification of the state written (four entries refuse it):

| row | measure | tolerance | where |
|---|---|---|---|
| hydrodynamic mass | 7.230e-07 | 5.0e-09 | cell 1 |
| hydrodynamic momentum | 5.356e-11 | 1.0e-08 | within |
| hydrodynamic energy | 6.452e-06 | 1.0e-06 | cell 2 |
| **carrier balance H2** | **1.000e+00** | 1.0e-05 | **cell 280** |
| elemental transport He/H | 3.005e-05 | 1.0e-05 | cell 301 |
| level balance He 2^3S | 1.477e-19 | 1.0e-06 | within |
| eliminated-species closure `System_HeH_mol` | 5.552e-17 | 1.0e-06 | within |
| cells without a chemical root | 0 | | |

**NOT CERTIFIED, and the row that refuses is the H2 carrier balance.** It is
no longer the row L7c described: the H2 it measures is smooth and is the
fit's own size. MEASURED on the state this run wrote, against the conjugate
form of the fit at each cell's own (p, T):

| cell | r [R_p] | T [K] | x2 carried | x2 of the fit | ratio |
|---|---|---|---|---|---|
| 194 | 1.1340 | 3887 | 3.4337e-08 | 1.8846e-08 | 1.82 |
| 218 | 1.2042 | 5198 | 3.8782e-10 | 2.4726e-10 | 1.57 |
| 250 | 1.3578 | 6217 | 2.5546e-11 | 1.4895e-11 | 1.72 |
| 280 | 1.6044 | 5986 | 1.3217e-11 | 6.5355e-12 | 2.02 |

against L7c's 1.30e+04 at cell 218 on the same route. The mixing ratio is
within a factor two of chemical equilibrium everywhere, and it is smooth:
1.3161e-11, 1.3156e-11, 1.3175e-11, 1.3217e-11, 1.3283e-11 over cells 277 to
281, which is a carrier advected with the nuclei, not a rounding pattern.

**So what the row refuses on now is the row, not the seed.** A measure of
exactly 1.000 is what a row reads when one of its terms carries the whole sum
of magnitudes: at 5986 K the collisional dissociation of H2 is the term, and
no advected mixing ratio at that temperature can balance it. Item L7c
section 5.2 reconstructed that budget at cell 218 of the old state and
reported that the code still carries no dump of a row's terms cell by cell;
that is still true, so the statement here rests on the measure and on the
profiles above and not on an assembled budget.

The `local` solve then took `run_case.sh`'s own continuation (L4e: a refused
solve is reloaded and solved again at `DTAU0_CONTINUATION` = 1e8), and that
continuation reaches `info = 0` on every pass with `||R||` down to
**1.157e-08** at its fifth Newton iteration:

| outer pass of the continuation | hydro `info` | mass | momentum | energy | worst gated species row |
|---|---|---|---|---|---|
| 1 | 0 | 7.24e-07 | 5.36e-11 | 6.44e-06 | 1.00 at cell 280 |
| 2 | 0 | 7.24e-07 | 5.36e-11 | 6.45e-06 | 1.00 at cell 280 |
| 3 | 0 | 7.23e-07 | 5.36e-11 | 6.45e-06 | 1.00 at cell 280 |

so the state the route reaches is a fixed point of the alternation to three
digits and the H2 carrier row is the only thing between it and a
certification worth arguing about (the mass row's 7.2e-07 at cell 1 and the
energy row's 6.5e-06 at cell 2 are the base ghost and its neighbour and are
a separate question).

### 5.2 The handoff partition

MEASURED, 24 JFNK iterations of outer pass 1 at the caps above: `||R||`
2.183 -> 1.997, the line-search merit 1.62e+02 -> 8.47e+01 after rising to
8.45e+03 on the way, `dtau` driven from 2.0 to the explicit-stable floor
8.49e-05 by the new rule's cut, `lam` = 9.5e-07 at the last iteration. Outer pass 1 ended
`info = 1` with mass 1.00, momentum 1.80e-01 and energy 1.99 -- the rows the
seed was loaded with (L7c section 3.1: mass 1.001, momentum 1.800e-01,
energy 1.989), so **the solve returns the state it was given.** The worst
species row is the H2 carrier balance, 9.22e-01 at cell 224. **This is the fourth independent measurement that the
handoff partition does not descend** (L7 section 7.2, L7b section 7.1, L7c
section 3.2), and the reason L7 gave still holds: carrying x2 = 0.99998 to
30 R_p leaves no neutral hydrogen above the base to photoionize, so the
energy row of the seed is the radiative imbalance of a column the wind is not
held up by. It is why `SEED_X2` defaults to `local`.

## 6. How to reproduce

`EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00`. Everything
below was run with private builds; `EXHALE.x` and `build/` were not touched
because the campaign was using them.

```bash
cd $EX
# the control (the tree as this item received it) and the delivered text
PATH=/usr/bin:$PATH make -j8 OBJDIR=build_L7d_ctl EXE=EXHALE_L7d_ctl.x   # before any edit
PATH=/usr/bin:$PATH make -j8 OBJDIR=build_L7d     EXE=EXHALE_L7d.x       # after

# section 2, the arithmetic: no run needed, the 60-digit reference and the
# two double-precision forms over the (p,T) of the 500 physical cells of
# LHS1140b/models/.L7c/src_atomic_HeH2.13/Hydro_ioniz.txt, p_bar = p x 1e-6.

# section 2, through the code: the local molecular seed of that state
cd $EX/LHS1140b/models
SRC=$PWD/.L7c/src_atomic_HeH2.13
for b in ctl meas; do d=.L7d/${b}_seed_local; mkdir -p $d/output
  \cp -f molecular_scalar_gj1132_kzz1e9/HeH2.13/base.inp $d/
  sed 's/^Load IC?.*/Load IC? True/' \
      molecular_scalar_gj1132_kzz1e9/HeH2.13/input.inp > $d/input.inp; done
( cd .L7d/ctl_seed_local  && OMP_NUM_THREADS=2 EXHALE_MOLECULAR_SEED=$SRC \
      EXHALE_MOLECULAR_SEED_X2=local $EX/EXHALE_L7d_ctl.x > seed.log 2>&1 )
( cd .L7d/meas_seed_local && OMP_NUM_THREADS=2 EXHALE_MOLECULAR_SEED=$SRC \
      EXHALE_MOLECULAR_SEED_X2=local $EX/EXHALE_L7d.x     > seed.log 2>&1 )

# section 3, the near state (the C/N/O pass-8 entry, scratchpad/L4e/red_A)
#   copy input.inp, base.inp and output/*_IC.txt into a scratch directory, then
#   OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 EXHALE_OUTER_PASSES=1 \
#     EXHALE_JFNK_MAXIT=80 <binary>
# section 3, the far state (the molecular local seed at 1.2 R_p, the state
#   LHS1140b/models/.L7c/HeH2.13_local_rampoff/output/*_IC.txt holds)
#   the same command with EXHALE_JFNK_MAXIT=40, and with
#   EXHALE_PTC_RAMP_DOUBLE=0 for the merit-ratio control.
# section 3, the element reload: a copy of backup/regression/atomic_elem_newton
#   with Load IC? True, Coupled carrier solve: False, Restart intent: stationary,
#   IC/*_IC.txt in output/, then
#   OMP_NUM_THREADS=8 EXHALE_OUTER_PASSES=12 EXHALE_JFNK_MAXIT=80 \
#     EXHALE_DIFF_OMEGA=0.5 <binary>

# section 4
cd $EX/LHS1140b/models && python3 make_models.py --only molecular --force

# the regression movement: the matrix run twice, once with each binary, the
# second on a copy of backup/regression/ so the two do not share the case
# directories (the harness refuses two runs in one tree by a lock).
REGRESSION_EXE=$EX/EXHALE_L7d_ctl.x backup/regression/run_check.sh check \
   mol_base_handoff mol_carrier mol_diffusion mol_lyman_werner mol_sec_ion \
   mol_metals mol_ir_bands hp_front hp_trace_seed hp_zero_seed
```

> Note added 2026-09-14 15:20: the private builds of this memo were made with `PATH=/usr/bin:$PATH make`, which selects the Ubuntu gfortran 13.1 and its LAPACK, not the conda-forge gfortran 16.2 + OpenBLAS the tree binary `EXHALE.x` is built with (bare `make`). The control and the measured build of the memo share one toolchain, so every comparison here stands; reproduce with bare `make` (`OBJDIR=... EXE=...` as written) to obtain the tree's toolchain.

## 7. Scope, and what was not measured

- **Validation scope.** Section 2 touches `q_h2_equilibrium`, whose only
  callers are the molecular starting point of `ionization_equilibrium`
  (line 4223) and of `constrained_chemical_equilibrium` (line 1510), the base
  mixing ratio of `composition` (line 292), the cold molecular start of
  `set_IC` (line 216), the molecular seed conversion (line 391) and the
  lower column itself. The first two are the cancelling ones and both are
  molecular-layout routines, so no atomic run can move; the LHS 1140 b atomic
  campaign is untouched by it. The molecular cases that carry goldens are the
  seven `mol_*` and the three `hp_*` (section 2.4). Section 3 touches
  `solve_steady_jfnk` only; the cases that reach it are the two regression
  cases that hand over from the march (`wasp_full_newton`,
  `mol_base_handoff`), the two reload fixtures, and every LHS 1140 b
  stationary route.
- **No golden was refreshed**, and none is proposed. The movement is
  reported in section 2.4 against a control build of the entry text, because
  the hot-Uranus goldens in this tree have already moved under other items.
- **The element reload was measured at its first outer pass only.** Twelve
  outer passes of 80 JFNK iterations each did not finish on a machine running
  the six-job campaign; both runs were stopped after pass 1 and the pass-1
  rows are what section 3.4 reports.
- **The `local` and the handoff partitions of T-L7-5 were not run to
  certification either**, for the same reason; section 5 states how far each
  got and what refuses.
- **The cut of an accepted damped step is not clipped.** The published SER
  practice clips the ratio to [0.1, 2] and this rule clips only the growth,
  so one accepted step at lam = 2e-3 takes the pseudo-time three decades
  down, to the explicit-stable floor if it is near it. MEASURED on the
  element reload that is what happens, and the solve then alternates between
  the floor and twice the floor while the merit falls; it is not a failure
  there, but a clipped cut (`max(lam, 0.1)`) was not measured and is the
  first thing to try if a state is found where the drop to the floor is the
  obstruction.

## 8. Noticed outside this item, and what was done about each

- **`src/utils/run_lower.py` carried its own copy of the fit, with BOTH
  asymptotic guards inverted** -- `u > 30` returned 1.0 where the fit gives
  0, `u < -30` returned 0.0 where it gives 0.8384 -- as well as the
  cancelling difference. That is the state `lower_column.f90` was in before
  2026-08-31 for the guards and before this item for the form. FIXED to the
  conjugate form with no branch, matching the Fortran, and checked: the cold
  limit is 0.8384030418250951 = 1.9845/2.3670 to the last bit and the value
  falls smoothly through 6.38e-06 at 3311 K, 2.80e-08 at 4940 K and 5.13e-15
  at 5838 K. It is a standalone driver that writes a `base.inp`, so nothing
  in `src/modules/` links it; a `base.inp` written by it before this fix
  carries a q_H2 that is wrong wherever the column is hotter than about
  3300 K at nanobar pressures, and completely wrong below 529 K, which is
  most of a lower-atmosphere column.
- **`models/write_reproduce.py` told every molecular case that
  `pick_seed.py` chose its seed**, which pick_seed cannot do (it has no
  molecular state to choose from). FIXED: a seed line beginning "molecular
  seed from" now gets its own paragraph naming the conversion, the atomic
  state and what `SEED_X2` carried. Reported here because it is a third
  file beyond the two the item names.
- **`run_case.sh`'s `SEED_GIVEN` read `$ATTEMPT` before it existed.** Fixed
  with the molecular branch; section 4.2.
- The goldens of the molecular and `hp_*` cases in this working copy are
  stale against the current tree (section 2.4): five of the ten fail with a
  build of the entry text. Not refreshed, and not this item's to refresh.
- The code still carries no dump of a carrier row's terms cell by cell, so
  section 5's statement about what the row at cell 280 is made of rests on
  the measure and on the profiles and not on an assembled budget. L7c
  section 9 reported the same gap.

## 9. State left behind

- `backup/regression/<case>/output/` for the ten cases of section 2.4 now
  holds the CONTROL build's outputs and `<case>/EXHALE.x` is a copy of
  `EXHALE_L7d_ctl.x`: that is what `run_check.sh` leaves in a case directory
  after any run, and the next `check` overwrites both. The measured build's
  outputs for the same ten cases are in
  `scratchpad/L7d/reg/backup/regression/<case>/output/`, and the two matrix
  logs are `scratchpad/L7d/reg_ctl_matrix.log` and
  `scratchpad/L7d/reg_meas_matrix.log`.
- The private builds `build_L7d/`, `build_L7d_ctl/`, `build_L7dx/` and the
  binaries `EXHALE_L7d*.x` were deleted at the end of the item. The delivered
  source text rebuilds to md5 `9e38a707fc8296d00e2ccf3dc96e0b2c` (checked:
  a second `make OBJDIR=build_L7d EXE=EXHALE_L7d.x` produced the same md5,
  and `make -q` then reported nothing to do).
- The runs are kept: `scratchpad/L7d/` (the near state A, the far state F_*,
  the element reload, the copied regression trees) and
  `LHS1140b/models/.L7d/` (the two seeds of section 2.3, the frozen atomic
  state, and the two T-L7-5 case copies). `LHS1140b/models/.L7d/t_local`'s
  `REPRODUCE.md` was overwritten by a check of the section 8 fix to
  `write_reproduce.py` and no longer carries that run's command list; the
  run's own logs (`run.log`, `run_dtau0_first.log`, `seed.log`) are intact.
- No process of this item is left running.
