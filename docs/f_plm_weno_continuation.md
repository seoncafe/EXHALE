# Phase F: the PLM -> WENO3 hand-off as a continuation, not a switch

Scratch investigation of `docs/open_defects_20260903_review.md` item 12.7 and
its Phase F recommendation, and of `TO_BE_DONE.md` item (S).

Everything below was built and run **outside the working tree**, in
`<scratch>/fweno`, from a copy of the tree's `src/` and `Makefile` taken at
tree revision `6d07d48` **with the working tree's uncommitted changes in
place** (the tree is not clean; see "Provenance"). No tree file was modified.

## 1. What the two-stage recipe actually changes

`Reconstruction scheme: PLM+WENO3` runs stage 1 with `rec_method = 'PLM'` and
flips to `'WENO3'` in a single step when `du < du_th_plm`
(`src/EXHALE_main.f90`, the `in_plm_stage` branch). The flip sets
`use_plm = .false.`, `use_weno3 = .true.`, and those two flags -- not the
`rec_method` string -- are what the discretization reads. They are read in
exactly five places:

| file | what changes |
|---|---|
| `states/Reconstruction.f90` | the face states: PLM slopes vs ESWENO3 nonlinear weights |
| `flux/Num_Fluxes.f90` (`Phys_flux`) | `use_plm`: the pressure is INSIDE the momentum flux |
| `states/Source.f90` (`source`) | `use_plm`: the geometric pressure source `(dA_p - dA_m)/dV * p_C` is added |
| `time_step/RK_rhs.f90` | `use_weno3`: the momentum row instead gets the face-pressure difference `(p_R - p_L)/dr` |
| `states/Apply_BC.f90` (`Apply_BC_W`, `Rec_BC`) | `use_weno3`: the outer free-outflow ghost is linearly extrapolated instead of copied |

So the switch is not "a different limiter". It simultaneously changes the
reconstruction, the *form* of the momentum equation (conservative
pressure-in-flux plus geometric source, versus a non-conservative face-pressure
gradient), and the outer boundary operator. Three changes, one step, no
measurement of what any of them costs.

`steady_residual.f90::assemble_residual` -- the residual the JFNK solver and
every acceptance gate are built on -- calls the same `Reconstruct` + `RK_rhs`
pair, so it inherits whichever scheme the flags currently select. The steady
gates after a two-stage run are therefore WENO3 gates, which is correct; the
point of item 12.7 is that the marched state handed to them was produced by a
PLM operator until one step earlier.

## 2. The instrument: a fixed-state operator difference

`EXHALE_OPDIFF=1` (this patch) takes the state at the moment the hand-off
fires and assembles the steady residual `R = du/dt` twice on the SAME interior
cells -- once with the PLM flags and its zero-gradient ghost, once with the
WENO3 flags and its extrapolated ghost -- and writes both, row by row and cell
by cell, to `output/opdiff_profile.txt` together with `u`, `dt_loc` and the
scale each residual row is measured against. `EXHALE_OPDIFF=2` stops the run
right after writing, so the measurement can be repeated cheaply.

`R_WENO3 - R_PLM` on a fixed state is the instantaneous kick the switch
applies to the right-hand side. It is a property of the state and the two
operators alone; it contains no time integration, so it separates "the
equations changed" from "the solution then moved".

## 3. The continuation

`reconstruction_continuation_rhs` (`steady_residual.f90`) assembles

    R_lambda(u) = (1 - lambda) R_PLM(u) + lambda R_WENO3(u)

by evaluating the whole right-hand side twice, once under each set of flags,
and combining `dF`, `S`, the face states and the stored interface fluxes
linearly. `Apply_BC_W` blends the two outer ghosts under the same lambda. The
three RK stages of the marching loop and `assemble_residual` both call it, so
the marched equation and the residual the gates read are the same equation at
every lambda.

Three properties were designed in, and two of them are measured in section 5:

- **The endpoints are exact and cost one evaluation.** `lambda <= 0` and
  `lambda >= 1` take a single branch with the corresponding flags set, so
  nothing but the pure scheme runs there. Only the open interval costs two
  right-hand sides.
- **A marching stage at intermediate lambda is the convex combination of the
  two schemes' forward-Euler stages.** The admissible set `rho > 0`,
  `rho e > 0` is convex, so a stage is admissible whenever both single-scheme
  stages are; the continuation cannot manufacture a positivity failure that
  neither scheme has. (`positivity_limited_fluxes` was given the matching
  lambda weight on its face-pressure term. Its first-order replacement flux is
  still built with the flags currently set; it fires on none of the runs
  below, and the counter says so.)
- **The final gates are pure WENO3.** When lambda reaches 1 the continuation
  disarms itself and sets the flags to WENO3, so the JFNK finish, the steady
  gates and every stop test after the ramp see the production operator with no
  blending left in it.

The environment variables named in this section were the INTERFACE OF THE
EXPERIMENT, and sections 4 to 7 are the measurements made through them. They
do not exist in the adopted form: section 148 (section 8) replaces them with
one input key, and keeps only `EXHALE_OPDIFF`, which measures a state and
changes no result.

Step control (`EXHALE_RECON_LAMBDA_STEP`, `..._DTUTOL`, `..._MODE`): lambda
advances by the current step each marching step while
`dtu <= tol * dtu(ramp start)`; a step that exceeds the tolerance halves the
lambda step and does not advance lambda; five consecutive accepted steps
double it back, capped at the initial value. `EXHALE_RECON_LAMBDA_MODE=fixed`
disables the control and ramps unconditionally.
`EXHALE_RECON_LAMBDA_FIX=<lambda>` pins lambda at the hand-off with no ramp,
which is how the endpoints are tested.

Everything is off unless an environment variable is set; with none set the
code takes the shipped branch.

## 4. Measured: the operator difference on a fixed state

Two cases, both from `backup/regression/`:

- **`mol_sec_ion`** (hot Uranus, molecular chemistry, `Secondary_ionization:
  Immediate`, `du_th [PLM,WENO3]: 0.5 1.0e-3`, 12000-step cap). The hand-off
  fires **mid-run at step 11861**, on a state that has been marching for
  11860 steps with `dtu` down to 5.9e-4. This is the case that can say
  anything about the switch.
- **`wasp_full`** (WASP-121b, He 2^3S + metals, same two-stage line). Its
  hand-off fires at **step 1**: a cold start already has du = 0.370, below
  the 0.5 stage-1 threshold. The state measured there is the initial
  condition after one step, in the middle of a `dtu = 2.6` transient.

### `mol_sec_ion`, step 11861 (dt = 1.386e-4)

Residual rows are divided by `residual_row_scale`, the largest term that row
itself contains.

| row | max abs R_PLM / scale | max abs R_WENO3 / scale | max abs dR / scale | r at that cell | dR relative to R_PLM |
|---|---|---|---|---|---|
| mass | 1.519 | 10.245 | 10.232 | 1.00019 | 6.74 |
| momentum | 1.678 | 2.379 | 1.036 | 4.68095 | 0.618 |
| energy | 1.924 | 1.925 | 1.202 | 1.00077 | 0.625 |

**The difference is not spread over the wind; it sits at the two boundaries.**
As a fraction of each row's own maximum:

| radial band | mass | momentum | energy |
|---|---|---|---|
| 1.00 - 1.02 | **1.000** | 0.058 | **1.000** |
| 1.02 - 1.05 | 1.1e-4 | 0.013 | 2.6e-3 |
| 1.05 - 1.10 | 1.1e-4 | 0.013 | 5.0e-2 |
| 1.10 - 1.30 | 6.5e-5 | 0.013 | 3.7e-2 |
| 1.30 - 1.60 | 5.2e-7 | 0.018 | 2.4e-4 |
| 1.60 - 2.00 | 3.5e-7 | 0.031 | 5.5e-5 |
| 2.00 - outer | 3.5e-3 | **1.000** | 0.554 |

The mass and energy rows differ almost entirely in the first two or three
cells above the base; the momentum row differs almost entirely in the
outermost cells (its worst cells are r = 4.59-4.73, against a domain ending
near 4.9). Those are exactly the two places where PLM and WENO3 differ by more
than the interior reconstruction -- the base ghost and pressure split at one
end, the extrapolating free-outflow ghost at the other. **In the wind between
them the two operators agree to a part in 10^4 to 10^2.**

### The measured transient is the operator difference, not an instability

One forward-Euler step of each operator, over the cells `dtu` is measured on
(`r >= r_esc = 2.0`):

| | predicted from the fixed state | measured in the run |
|---|---|---|
| dtu of a PLM step | 5.979e-4 | **5.948e-4** (step 11861, the last PLM step) |
| dtu of a WENO3 step | 1.095e-3 | **1.015e-3** (step 11862, the first WENO3 step) |
| ratio | 1.831 | **1.706** (peak over the next 400 steps: 1.763) |

The PLM prediction is right to 0.5 percent, and the WENO3 prediction
overshoots by 8 percent -- which is what a forward-Euler estimate of an
SSP-RK3 step with an ionization sweep inside it should do. **The jump in dtu
at the hand-off is quantitatively the operator difference on a fixed state.**
Nothing in the following 400 steps needs a WENO3 instability to explain it:
`du` never rises at all in this case, it keeps descending monotonically
(0.4998 at the switch, 0.4194 a hundred steps later, 0.3572 at step 14000),
and no positivity counter fires.

### `wasp_full`, step 1 (dt = 1.013e-4)

| row | max abs R_PLM / scale | max abs R_WENO3 / scale | max abs dR / scale | r at that cell | dR relative to R_PLM |
|---|---|---|---|---|---|
| mass | 0.294 | 0.329 | 0.0528 | 1.29681 | 0.180 |
| momentum | 1.996 | 1.996 | 0.222 | 1.29897 | 0.111 |
| energy | 1.224 | 1.591 | 0.367 | 1.00059 | 0.300 |

Over the `dtu` window (`r >= 1.5`) the predicted step sizes are equal to four
digits and the added kick is 7.1e-5 against a step that is itself moving the
state by `dtu = 2.55`. **On `wasp_full` the hand-off is measurably harmless**,
because it happens inside a cold-start transient four orders of magnitude
larger than the operator difference. This case cannot be used to argue either
way about item 12.7; it is included because it is the atomic two-stage case
the review's validation matrix asks for.

## 5. Measured: the continuation

### 5.1 The endpoints reproduce the shipped code

| test | result |
|---|---|
| build with no environment variable set, `wasp_full` | stop step 11289, du 9.9989e-4, log10 Mdot 13.26; `Hydro_ioniz.txt` and `Ion_species.txt` **data lines bitwise identical to `backup/regression/golden/wasp_full`** |
| same, `mol_sec_ion` (12000-step cap) | stop step 12000, du 4.0963e-1, log10 Mdot 10.38; both files **data lines bitwise identical to the golden** |
| `OMP_NUM_THREADS=1` vs `4`, `wasp_full` | outputs bitwise identical |
| `EXHALE_RECON_LAMBDA_FIX=1`, `wasp_full` (blend path pinned at lambda = 1 for the whole run) | stop step 11289, du 9.9989e-4, Mdot 13.26; **every data line bitwise identical to the unpinned run** |
| `EXHALE_RECON_LAMBDA_STEP=1`, `mol_sec_ion` | both output files bitwise identical to the baseline run, and the whole du/dtu trace matches digit for digit |

The only byte that moves under `EXHALE_RECON_LAMBDA_FIX=1` is the provenance
string in the output header, which reads `recon=PLM` because `rec_method` is
left on `'PLM'` while the operator is WENO3 through the blend. **If any of
this is adopted, that header has to report the effective lambda**; a run whose
header says PLM and whose equations are WENO3 is exactly the provenance defect
the review's section 7 asks to close.

(The goldens themselves are one header line behind this working tree: the tree
has an uncommitted `write_output.f90` that adds a `# coupling:` line the
goldens do not carry. That is pre-existing and unrelated to this work; the
comparisons above are of the data lines.)

### 5.2 `mol_sec_ion`: what the ramp does to the transient

All runs identical up to step 11861 (checked at steps 1000 and 4000, and by
the bitwise `lambda`-step-1 result above). `dtu` immediately before the
hand-off is 5.948e-4.

| run | ramp | lambda = 1 at | back-offs | peak dtu after the hand-off | ratio to pre |
|---|---|---|---|---|---|
| baseline (shipped switch) | - | step 11861 | - | 1.060e-3 | **1.763** |
| `LAMBDA_STEP=1` | - | step 11861 | - | 1.060e-3 | 1.763 (identical trace) |
| `LAMBDA_STEP=0.02 MODE=fixed` | 49 steps | 11910 | 0 | 7.501e-4 | **1.248** |
| `LAMBDA_STEP=0.02` adaptive, tol 1.2 | 58 steps | 11919 | 2 | 7.213e-4 | **1.200** |
| `LAMBDA_STEP=0.002` adaptive, tol 1.2 | 499 steps | 12360 | 0 | 5.948e-4 | **0.990** |

A 500-step ramp removes the excursion entirely: `dtu` never rises above its
pre-hand-off value at any point. A 50-step ramp cuts the excursion from 76
percent to 20-25 percent. The adaptive control fires (2 back-offs at
`tol = 1.2`, none at 500 steps) and the tolerance is what it is set to.

`du` descends monotonically in every case, so there is no "steps to recover"
to report -- there is nothing to recover from in this case. The ramp does put
the run slightly behind: at 400 steps past the hand-off, du is 0.38105
(baseline), 0.38129 (50-step ramp), 0.38561 (500-step ramp); by step 12660 the
three are 0.35512, 0.35532, 0.35656.

**Cost.** Wall time for the same 14000 steps, four threads, on a loaded
machine: baseline 10 m 37 s, 50-step ramp 9 m 34 s, 500-step ramp 9 m 36 s.
The doubled right-hand side is paid on 49 to 499 of 14000 steps and is not
measurable against the scatter of a shared machine.

**The state at the fixed 14000-step cap** (this case never converges; it is a
relaxation snapshot):

| run | du at 14000 | log10 Mdot | max rel move of the profile vs baseline |
|---|---|---|---|
| baseline | 0.35715 | 10.35 | - |
| `LAMBDA_STEP=1` | 0.35715 | 10.35 | 0 (bitwise) |
| 50-step ramp, fixed | 0.35768 | 10.35 | rho 2.7e-3, T 3.5e-3, p 6.0e-3 |
| 50-step ramp, adaptive | 0.35782 | 10.35 | rho 3.4e-3, T 4.4e-3, p 7.5e-3 |
| 500-step ramp | 0.36217 | 10.35 | rho 2.4e-2, T 2.1e-2, p 3.2e-2 |

Every one of those maxima sits at r = 1.25-1.32, which in this snapshot is
where the velocity passes through zero and reverses (v = +735 cm/s at
r = 1.291, -2240 cm/s at r = 1.305). The wind above and the base below move
far less. The moves are the ramp arriving at the cap a few hundred steps
"younger", not a different solution.

### 5.3 `wasp_full`: the atomic case

| run | stop step | du | log10 Mdot |
|---|---|---|---|
| baseline | 11289 | 9.9989e-4 | 13.26 |
| 50-step ramp | 11288 | 9.9821e-4 | 13.26 |
| 500-step ramp | 11584 | 9.9358e-4 | 13.25 |
| **held on PLM** (`LAMBDA_FIX=0`, 15000-step cap) | 15000, **not converged** | 3.9114e-3 | 13.34 |

The ramps cost nothing and change nothing: the 50-step ramp lands on the same
step count to within one and moves the converged profile by at most 5.6e-5 in
rho and 3.8e-5 in T; the 500-step ramp stops 295 steps later at a slightly
tighter du and moves it by up to 2.4e-2 in rho, which is the known
path-dependence of a du stop, not a change of solution.

**The PLM-held run is the interesting one.** It does not reach du < 1e-3 in
15000 steps, and its state is not a perturbation of the WENO3 one: at
r = 1.06-1.30 its density is up to 2.4 times higher, its base is 400-700 K
hotter, its base velocity is +2.6e3 cm/s where the WENO3 run has -7.3e2, and
its mass-loss rate is 13.34 against 13.26. **PLM and WENO3 do not share a
discrete solution near the base of this planet.** That is the real content of
"do not switch the operator and expect the existing state to remain a
solution": the continuation is a path between two genuinely different discrete
problems, not a cosmetic smoothing. (Caveat: the PLM run is not converged, so
part of that gap is its own non-convergence; the direction and the base
structure are not.)

### 5.4 How far the goldens would move if a ramp became the default

`mol_sec_ion` re-run at its own 12000-step golden cap. The baseline column is
the run with no environment variable set, whose two output files are bitwise
identical to `backup/regression/golden/mol_sec_ion` (data lines).

| run | du at 12000 | log10 Mdot | max rel move vs baseline | where |
|---|---|---|---|---|
| baseline | 0.40963 | 10.38 | 0 (golden) | - |
| 50-step ramp (`STEP=0.02`, tol 1.2) | 0.41166 | 10.38 | rho 6.1e-3, T 2.1e-3, cool 1.2e-2 | r = 4.814 |
| 500-step ramp (`STEP=0.002`, tol 1.2) | 0.42444 | 10.38 | rho 6.9e-2, T 1.8e-2, cool 1.3e-1 | r = 4.814 |

Both ramps leave `log10 Mdot` at 10.38 and move `du` by 0.5 and 3.6 percent.
The largest profile move is at r = 4.81, the outermost cells -- the same place
the momentum-row operator difference of section 4 peaks. `wasp_full`'s golden
would move by at most 5.6e-5 (50-step ramp) or 2.4e-2 (500-step ramp, which
also stops 295 steps later); its Mdot stays 13.26 and 13.25. The six cases
that never leave the PLM stage are untouched.

## 6. Verdict

**On the two cases measured, the shipped discontinuous hand-off is (a), a
harmless transient -- and the reason it is harmless is measurable, not
assumed.** The `dtu` excursion it causes is quantitatively the fixed-state
operator difference (1.83 predicted, 1.71 measured on the first step), it
decays within about fifty steps (dtu is back below its pre-hand-off value by
step +50), `du` never stops descending, and no positivity guard fires.

**That verdict does not extend to the case item (S) reported**, and this work
did not reproduce it. Item (S) is HD 189733 b with molecular and oxygen
chemistry, transport and IR bands, restarted and marched 23896 steps to
du = 0.5; there the same hand-off is recorded as raising `dtu` 5.5e-5 ->
4.1e-4 in one step and `du` to ~1e5 within 2000 steps. That configuration is
not in the regression matrix and a single trajectory to it costs tens of
thousands of steps, so it was out of reach here. What this work does supply is
the instrument to settle it: `EXHALE_OPDIFF=2` measures `R_WENO3 - R_PLM` on
that state and stops, and the answer is decidable in one run. If the ratio
there is the ~7.5 that the reported `dtu` jump implies, the switch is the
cause; if it is ~1.8 like here, the blow-up is a WENO3 instability that the
switch merely uncovers, and the continuation will not fix it.

**What the measurement does say about (b).** The operator difference is not
uniform: over the wind the two schemes agree to a part in 10^2 to 10^4, but in
the first two or three cells above the base the WENO3 mass-row residual is
6.7 times the largest PLM residual anywhere on the grid, and in the outermost
cells the momentum rows differ by 62 percent. A `wasp_full` run held on PLM
does not converge in 15000 steps and settles on a base that is 400-700 K
hotter and up to 2.4 times denser than the WENO3 one, with the opposite sign
of base velocity and a mass-loss rate 0.08 dex higher. **The two
discretizations do not share a solution near the boundaries**, so the
two-stage recipe is a continuation between two different discrete problems.
Treating it as one -- which is what the review asks for -- is right on the
physics even where the transient it produces is benign.

**If it were adopted.** With `lambda = 1` the code is the shipped code to the
bit, so a default-off continuation moves no golden at all. A default-ON ramp
would move only the two matrix cases whose runs cross the hand-off inside
their step budget, by the amounts tabulated in 5.4 (`mol_sec_ion`: du +0.5 to
+3.6 percent, Mdot unchanged at 10.38, profile at most 0.6 to 6.9 percent, all
of it in the outermost cells) and 5.3 (`wasp_full`: Mdot 13.26 -> 13.26 or
13.25). The six cases that never leave the PLM stage stay byte-identical. So
a default-on ramp is a deliberate golden refresh of exactly two cases, of that
size.

**Recommendation, and what was decided.** Keep the continuation off by
default and expose it as an input key rather than an environment variable,
with the ramp length the user's choice; make the output header report the
effective lambda; and use `EXHALE_OPDIFF` as the standing measurement whenever
a two-stage run is suspected of losing its state at the hand-off. Adopting a
default ramp is not justified by these two cases -- it buys a 20-25 percent
smaller `dtu` excursion that decays on its own within fifty steps -- and would
be justified by the (S) configuration if the operator difference there turns
out to be what drives it.

That is the form section 148 implements; section 8 below is the adopted
version, its input key, and its measurements.

## 7. What was NOT checked

- **Item (S)'s own configuration was not run.** See section 6. Everything said
  about it here is inference from its recorded numbers plus this measurement
  on a different case.
- **The JFNK finish was never exercised through a ramp.** `mol_sec_ion` stops
  at its 12000/14000-step cap with du ~ 0.36, far above the 1e-2 Newton
  hand-off, and `wasp_full` has no `Solver:` line. The residual blend is
  reached only by the periodic residual monitor. **A continuation that ends
  inside a JFNK solve -- which is the form the review's Phase F actually asks
  for -- is implemented but untested.**
- **The operator difference was not decomposed into its three causes.** It is
  measured as one number per row per cell. Which part of it is the
  reconstruction, which the pressure split, and which the ghost extrapolation
  is not separated; the radial localization (base and outer boundary) points
  at the ghosts and the pressure split, but that is an inference, not a
  measurement.
- **`positivity_limited_fluxes` at intermediate lambda is only partly on the
  homotopy.** Its face-pressure term carries the lambda weight; its
  first-order replacement flux is built with the flags currently set. The
  counter is zero on every run in this document, so it never mattered here.
- **No other regression case was run.** The eight remaining matrix cases were
  not touched; with the continuation off the code path is the shipped one, and
  the two cases that were run reproduce their goldens bitwise, but that is not
  a full `make check`.
- **The eight non-hand-off matrix cases were not re-run under section 148.**
  With the key absent the code path is the shipped one and the two cases that
  do cross the hand-off reproduce bitwise, but a full `make check` at
  application time is still the thing that says so.
- **The machine was shared.** Other sessions of the same user were running
  EXHALE and Cloudy throughout; wall times in section 5.2 are therefore lower
  bounds on precision, not clean timings. Step counts and all output numbers
  are unaffected.

## 8. Section 148: the adopted form

Section 148 is the shipped form of everything above: the continuation kept,
**off by default**, driven by an input key rather than an environment
variable, with the file header made to report the operator the state was
actually produced under.

### 8.1 The input key

```
Reconstruction continuation: <dlambda> [<dtu_tol>] [fixed|adaptive]
```

- `dlambda` -- step in lambda per marching step. `<= 0`, or an absent key, is
  the shipped one-step hand-off. `1.0` puts lambda at 1 at the hand-off and is
  the one-step hand-off, reached through the continuation code.
- `dtu_tol` -- the step control, default `1.2`: lambda advances while
  `dtu <= dtu_tol * dtu` at the start of the ramp; a step that exceeds it
  halves the lambda step instead of advancing lambda.
- `fixed` / `adaptive` -- disable or enable that control; `adaptive` is the
  default.

The two optional fields are told apart by content, not by position, because
the mode is a word and the tolerance is a number; positional-only parsing
would have made `... 0.02 fixed` illegal for no reason a user could see.

**Why this name and this shape.** `input_read.f90`'s dominant form for a
tuning line is numeric-first with optional trailing fields on the same line --
`Shapiro filter: <eps> [<every>]`, `Low-Mach damping: <eps4> [<M_th>]`,
`du_th [PLM,WENO3]: <du1> [<du2>]` -- and the two-word key puts the first
value at `get_word(line, 3)`, which is what the surrounding branches expect.
The key sits next to `du_th` in the parse chain because both are settings of
the same two-stage recipe. `lbl_match` anchors on the full two-word key, so
`Reconstruction continuation` and `Reconstruction scheme` cannot be confused,
and the duplicate-key check attributes a line to the longest matching key.

**It is refused rather than ignored when it cannot apply.** The continuation
walks a hand-off; a run without one has nothing to walk. Resolved at the end
of `input_read`, once `Reconstruction scheme:` and both `du_th` values are
known: if the run is not two-stage, or its PLM stage is empty
(`du_th_plm <= du_th`), the key is cleared and a WARNING is printed. Measured
on a single-stage WENO3 input with the key present:

```
 (input_read) PLM -> WENO3 continuation: dlambda = 2.00E-02, dtu tolerance = 1.50E+00, fixed
 (input_read) WARNING: "Reconstruction continuation" needs a two-stage run with a
 non-empty PLM stage ("Reconstruction scheme: PLM+WENO3" and two du_th values): ignored.
```

`EXHALE_RECON_LAMBDA_STEP`, `..._DTUTOL`, `..._MODE` and `..._FIX` are gone.
`EXHALE_OPDIFF` stays an environment variable, with the other `EXHALE_*`
diagnostic probes (`EXHALE_RESIDUAL`, `EXHALE_JAC_TEST`, `EXHALE_PTC`,
`EXHALE_PARSE_DUMP`): it measures a state and changes no result.

### 8.2 The provenance defect, fixed at its single definition

`rec_method` is not the name of the operator a run is solving with. While the
continuation is armed the run marches `R_lambda`, which is neither scheme, and
`rec_method` is left on `'PLM'` -- so a header built from it says PLM of a
state produced by neither. The tree wrote `recon=trim(rec_method)` in **two**
places, `write_provenance_header` and `write_coupling_state_header`
(`utilities.f90`), which is the duplicate this had to avoid multiplying.

`reconstruction_operator_label()` is now the one definition, in
`global_parameters` next to the state it reads, and both writers call it:

Measured, by stopping `mol_sec_ion` inside a slow ramp
(`Reconstruction continuation: 0.001`, `EXHALE_MAXSTEPS=11900`, hand-off at
step 11861, so lambda = 0.04 at the write):

```
# coupling: sec_ion=T sec_ion_step=0 valve=-1.00000E+00 recon=PLM+WENO3(lambda=0.0400) fluxconst=-1.00000E+00
# provenance: recon=PLM+WENO3(lambda=0.0400) base_bc=isothermal carrier=local restart_schema=3 resid_def=145 N=500
```

Both lines, from the one function. A snapshot written before the hand-off says
`recon=PLM` and one written after the ramp finishes says `recon=WENO3`, both
verified in the same run's files.

### 8.3 Rebase status

The tree moved while this was being written: **section 145 landed in the
working tree at 19:00**, and `utilities.f90`, `input_read.f90`,
`parameters.f90`, `steady_residual.f90`, `steady_newton.f90` and
`EXHALE_main.f90` all moved with it. Sections 146 (AB) and 147 (AF) have
**not** landed -- `Apply_BC.f90` and `RK_rhs.f90` are still untouched -- and
the scratch trees that hold them (`p_ab`, `p_af`) were branched from the tree
BEFORE section 145, so no "section 147 applied" source exists to rebase onto.

Section 148 is therefore developed against the **current tree**, and its
relation to the pending sections was measured rather than assumed, by
test-applying the patch to the section-147 source:

| file | result on the section-147 source |
|---|---|
| `parameters.f90` | applies (offset -5) |
| `input_read.f90` | applies (offsets +1, -19) |
| `write_setup_report.f90` | applies |
| `Apply_BC.f90` | applies (offsets +25) |
| `RK_rhs.f90` | applies |
| `steady_residual.f90` | 3 of 5 apply (offsets -5) |
| `EXHALE_main.f90` | 7 of 8 apply (offsets up to +88) |

Every hunk that fails is a hunk whose context section 145 created -- the
`# provenance: recon=` line, the `public ::` list and third residual call site
of `steady_residual.f90`, and the `use steady_residual_mod` list of
`EXHALE_main.f90` -- and section 145 is already in the tree. **No hunk of
section 148 overlaps a hunk of section 146 or 147**: the Apply_BC change of
146/147 rewrites `Apply_BC` (the conservative wrapper), section 148 edits
`Apply_BC_W` (the primitive ghost fill), and they are disjoint text. So the
patch is expected to apply after 146, 147 and block H with offsets only; it
will be re-applied and re-checked against the tree at signal time rather than
assumed.

### 8.4 Goldens

Measured against an unpatched build of the same tree source, four threads, the
regression comparison's own rule (`grep -v '^ *#'`, so the build-time stamp of
the section-145 provenance line is excluded as it always was):

| case | run | stop | data lines vs the unpatched control |
|---|---|---|---|
| `wasp_full` | control (unpatched) | 11289, du 9.9989e-4 | - (and bitwise equal to the stored golden) |
| `wasp_full` | key absent | 11289, du 9.9989e-4 | **identical** |
| `wasp_full` | `Reconstruction continuation: 1.0` | 11289, du 9.9989e-4 | **identical** |
| `mol_sec_ion` | control (unpatched) | 12000, du 4.0963e-1 | - (and bitwise equal to the stored golden) |
| `mol_sec_ion` | key absent | 12000, du 4.0963e-1 | **identical** |
| `mol_sec_ion` | `Reconstruction continuation: 1.0` | 12000, du 4.0963e-1 | **identical** |
| `mol_sec_ion` | `Reconstruction continuation: 0.02` | 12000, du 4.1166e-1 | differs, as it must |
| `mol_sec_ion` | key absent, final binary (repeat) | 12000, du 4.0963e-1 | **identical** |

`Reconstruction continuation: 1.0` is the sharper of the two identity tests:
it goes through the whole continuation -- the blended right-hand side at
lambda = 1, the blended ghost fill at lambda = 1, the lambda-weighted
face-pressure term of the positivity flux correction, the controller and its
disarm -- and lands on the shipped result to the bit, on both an atomic and a
molecular case. Its du/dtu trace matches the control digit for digit at every
step.

**The parse corpus does not move.** `parse_dump.txt` gains the three
continuation values only when the key is armed, so an input without the key
dumps exactly what it dumped before -- verified directly (`mol_sec_ion`
parse dump byte-identical between the unpatched and patched builds; with
`Reconstruction continuation: 0.02 1.3` it gains exactly
`recon_lambda_step0`, `recon_lambda_dtu_tol`, `recon_lambda_adaptive`). None
of the 1122 files in `backup/regression/parse_golden/` needs a refresh.

**No golden moves.** Section 148 is off by default and byte-identical when
off, so `make check` and the parse corpus are expected to pass unchanged; the
two cases that cross the hand-off were the ones at risk and both reproduce.
The eight matrix cases that never leave the PLM stage were not re-run -- with
the key absent the code path is the shipped one, which is what the two cases
above measure -- and the full `make check` should still be run once at
application time.

**The continuation itself is unchanged by the rewrite.** With
`Reconstruction continuation: 0.02` the `mol_sec_ion` ramp is 58 steps with 2
back-offs, the dtu excursion is 1.200 times its pre-hand-off value against
1.763 for the one-step switch, and the 12000-step snapshot moves by at most
6.1e-3 in rho at r = 4.81 with `log10 Mdot` unchanged -- every number
identical to the environment-gated measurement of sections 5.2 and 5.4.


## 9. Provenance of the measurements in sections 4 to 7

- Tree revision `6d07d48afd414dad229e53d4fa441c050623b243`, **working tree not
  clean** -- the scratch copy was taken from `src/` and `Makefile` as they
  stood, including uncommitted changes.
- Built with `PATH=/usr/bin:$PATH make` (gfortran 13, `-O3 -fopenmp`).
- Runs at `OMP_NUM_THREADS=4` except the two golden comparisons, which used
  `OMP_NUM_THREADS=1` as `run_check.sh` does. `wasp_full` was verified
  bitwise identical between 1 and 4 threads.
- Scratch root `<scratch>/fweno`; run directories under `runs/`, patch in
  `docs/f_plm_weno_continuation.patch`, raw operator-difference summaries in
  `docs/opdiff_mol_sec_ion.txt` and `docs/opdiff_wasp_full.txt`.
- Three binaries were built as the patch grew; all three take the shipped
  code path when no environment variable is set, and each was checked
  against a golden on that path:
  - `f689667466b0cf94d5beb7edd20e897a` -- continuation + `EXHALE_OPDIFF`.
    Ran `wasp_gold`/`wasp_base` (bitwise vs the `wasp_full` golden),
    `gold12000` (bitwise vs the `mol_sec_ion` golden), and the 14000-step
    ramp series.
  - `248a65ea844c7b0ff777046bfacf2e2c` -- adds `EXHALE_RECON_LAMBDA_FIX`.
    Ran the three pinned `wasp_full` runs; `wasp_pinref` (no environment
    variable) stopped at the golden step 11289 with the golden du/dtu.
  - `861727394787caba4ad600edccca2f76` -- adds `u`, `dt_loc`, `dt` to the
    operator-difference dump. Ran both `EXHALE_OPDIFF=2` probes and the
    12000-step golden-move series; `g12_base` (no environment variable) is
    bitwise identical to the `mol_sec_ion` golden. This is the state the
    patch file carries.
- Source md5 of the patched files (final state):

```
e00f95a15debc8dd256178b8bd30475e  src/modules/init/parameters.f90
535a7339a1915c6dabdd43c13194b56f  src/modules/states/Apply_BC.f90
97012cdbef97d5fda9971a6a13c531ff  src/modules/time_step/RK_rhs.f90
87f64806309d08ea00facb47eecb1bd1  src/modules/time_step/steady_residual.f90
b9b43378a0e15bf24a6b5e01d79b8a23  src/EXHALE_main.f90
```
