# P50. Why the carrier/wind alternation does not converge, and what a stable one costs

Scratch campaign, 2026-09-03. Everything below was run from
`.../scratchpad/picard`; nothing in the repository tree was read except as
documentation, and nothing in it was changed. The source is a copy of
`.../scratchpad/frontmode/src` (the P49 working copy: continuous positivity
repair plus the refusal of an uncertified composition), built with
`PATH=/usr/bin:$PATH make`.

## 0. Judgment

**The alternation does not converge because, on this state, there is nothing
for it to converge to.** Three things are measured below and each is
independent of the others.

1. **The drift is real and it is not a small-number artefact.** At the first
   hand-off the fixed-wind carrier relaxation drives `x2 = 2n(H2)/n_H` to the
   hydrogen element ceiling, exactly 1.0, in 123 of 500 cells between 1.047
   and 1.227 R_p. In the late hand-off the largest change sits at
   r = 1.795 R_p, where the H2 mixing ratio goes from 0.252 to 0.507 -- a
   quarter of the gas by number. Replacing the drift measure with a relative
   one carrying an absolute floor (option (a)) is therefore **not indicated**:
   the measure in the code is already normalized by the largest H2 mixing
   ratio anywhere on the grid, which is an absolute measure, and the cells it
   flags are cells where H2 is a major constituent.

2. **The code's own carrier convergence measure understates the imbalance by
   about 4000.** `carrier_row_scale` is `n_H(|v|+c_s)/dr`. Measured against
   the sum of the row's own term magnitudes -- the quantity
   `carrier_residual` already computes as `dsc`, and the quantity section 133
   made the principle for every hydrodynamic row -- that scale is 10^2.7 to
   10^6.4 too large. The state the JFNK accepted at `info = 0` reports a
   carrier steady residual of 1.6e-4 and carries a residual of **0.60 of the
   row's own largest terms** at r = 1.241, above 0.1 in 52 cells between
   1.195 and 1.436 R_p. It is not a carrier steady state and the reported
   number cannot say so.

3. **The alternation can be made stable, and it still does not converge.**
   Handing the wind a bounded transport advance instead of the fixed-wind
   steady state blended at omega keeps **every** JFNK at `info = 0` --
   20 passes out of 20, `||R||` between 2.4e-4 and 3.1e-4 -- while the H2
   front walks outward one cell in each pass, 1.203 -> 1.250 R_p, at
   essentially the local gas speed. Larger bounds move it faster and no
   further per unit of composition spent: 1.210 -> 1.464 R_p in 19 passes at
   a bound of 0.05, and 1.227 -> 1.447 in 8 passes at 0.2. The alternation is
   time integration wearing a fixed-point loop's clothes, and the front it
   integrates is the propagating front sections 128.6 and 134.4 already
   report.

The recommendation in section 5 follows from those three, and it does not
include a scheme that makes the loop converge on this state.

## 1. Reproduction, made cheap

`E1d_fix2` reaches its accepted hand-off at marching step 64,120 after about
half an hour. It also has an EARLIER one that the P49 write-up does not
separate out: the very first JFNK of the run, at marching step 2, already
returns `info = 0` at `||R|| = 5.028e-4` (`E1d_fix2/run.log` lines 38-57), and
the pass that follows it throws `||R||` to 8.166e-2. So the failure is
reproducible in **three seconds**, not thirty minutes.

The state is captured by one addition to the scratch source: in
`steady_wind_with_element_diffusion`, after `if (jfnk_info .ne. 0) exit` and
before any carrier pass, `EXHALE_PICARD_DUMP=1` writes the accepted wind with
`write_output` and stops. Copying `output/{Hydro_ioniz,Ion_species}.txt` to
`output/{Hydro_ioniz,Ion_species}_IC.txt` of a new directory restarts from
exactly that wind (`Load IC? True` is already set). That state is
`.../scratchpad/picard/cap/output/`; every experiment below starts from it and
reaches its first post-pass JFNK in under a minute.

Restarted, the wind re-converges to `||R|| = 1.454e-4`, `info = 0`, and the
carrier pass that follows reproduces the failure: drift 0.553, worst at
cell 271, r = 1.227, `||R||` handed on 6.745e-2, `info = 2`.

## 2. Anatomy of one pass

`EXHALE_PICARD_PROBE=1` writes `output/carrier_pass_NN.txt` at the end of each
pass: for every cell, the entry and exit carrier mixing ratio, `x2`, the
relaxation step, and the carrier residual on both scales. It also traces the
inner relaxation.

### 2.1 What the pass changes

Pass 1 from the accepted wind (`picard/A1`), undamped, at the frozen wind of
the `info = 0` solve:

| r [R_p] | v [cm/s] | T [K] | x2 entry | x2 exit | fc(H2) entry | fc(H2) exit |
|---|---:|---:|---:|---:|---:|---:|
| 1.050 | 1.6e+1 | 711 | 0.9425 | **1.0000** | 7.75e-1 | 8.22e-1 |
| 1.100 | 1.8e+2 | 663 | 0.7602 | **1.0000** | 5.44e-1 | 7.15e-1 |
| 1.150 | 1.8e+3 | 873 | 0.4447 | **1.0000** | 2.60e-1 | 5.84e-1 |
| 1.181 | 4.7e+3 | 1078 | 0.2241 | **1.0000** | 1.16e-1 | 5.17e-1 |
| 1.200 | 7.4e+3 | 1202 | 0.0963 | **1.0000** | 4.67e-2 | 4.85e-1 |
| 1.227 | 1.2e+4 | 1413 | 6.3e-4 | **1.0000** | 2.90e-4 | 4.63e-1 |
| 1.250 | 1.7e+4 | 1565 | 8.1e-6 | 0.2488 | 3.77e-6 | 1.15e-1 |
| 1.400 | 6.0e+4 | 2282 | 5.9e-5 | 4.4e-5 | 2.75e-5 | 2.05e-5 |
| 1.835 | 2.2e+5 | 2740 | 2.7e-5 | 1.1e-5 | 1.24e-5 | 5.07e-6 |

`x2 = 1.0000` is not a rounding: **123 cells between 1.047 and 1.227 R_p sit
exactly on the hydrogen element ceiling** that `limit_to_element_budget` and
`hydrogen_available_to_carriers` impose (entry: none). A state pinned by the
limiter over a quarter of the grid is not a solution of the carrier equation;
it is the equation's clamp.

The far wind moves by 7e-6 in fc, four to five decades below the front's move.
The drift of 0.553 is entirely the front region.

### 2.2 The inner relaxation converges; the pass is not bounded

The relaxation exits on its own criterion at step 63 (`dstep` 2.5e-10 against
`relax_tol` 1e-10), with `grow` at 2.4e10, i.e. `dt` about 1e10 cell crossing
times. **Option (b) of the brief is therefore already in force**: the pass
already solves the steady carrier equation at the fixed wind, and the 0.553 is
what that steady solution is. There is nothing to gain by asking for it a
second way.

The trace shows where the move happens -- steps 1 to 9, while `grow` is still
of order 1 to 30:

| k | grow | drift of the step | cumulative drift | worst r |
|---:|---:|---:|---:|---:|
| 1 | 1.0 | 1.13e-2 | 1.13e-2 | 1.152 |
| 3 | 2.3 | 2.52e-2 | 5.34e-2 | 1.160 |
| 5 | 5.1 | 5.46e-2 | 1.44e-1 | 1.177 |
| 7 | 11.4 | 1.07e-1 | 3.20e-1 | 1.205 |
| 9 | 25.6 | 1.25e-1 | 5.44e-1 | 1.218 |
| 20 | 2.2e3 | 2.72e-3 | 5.56e-1 | 1.232 |
| 60 | 2.5e10 | 2.50e-10 | 5.56e-1 | 1.232 |

**The first step already moves the composition by 1.13e-2.** That step is at
`dt = min(t_diff, t_adv)`, the cell's own crossing time, the smallest step the
operator takes. So the Picard drift gate of 1e-3 is BELOW the smallest move
this operator can make, and no amount of damping can make the gate fire: the
gate tests the undamped move, by design (the routine header says so, and the
design is right). This is measured, not inferred -- `EXHALE_PICARD_TRUST=1e-3`
and `=1e-2` produce bit-identical runs (`picard/E1`, `picard/E2`).

### 2.3 The drift measure is already absolute, and the late pass proves it

`relax_photochemical_composition` computes

```
drift = max_{j,ic} |fc_now(j,ic) - fc_entry(j,ic)|
        / max( max_j fc_entry(j,H2), max_j fc_now(j,H2) )
```

The denominator is the largest H2 mixing ratio anywhere on the grid, one
number, not the cell's own value. It is an ABSOLUTE measure with a global
normalization, so the failure mode option (a) protects against -- a cell where
H2 is 1e-5 reporting a large relative change -- cannot occur.

The late pass confirms it directly. `picard/A0` reproduces the whole
`E1d_fix2` sequence with the probe on, and its second hand-off is the one
P49 section 8 reports (400 relaxation steps, worst drift in the wind rather
than at the front). Measured, at that hand-off the front has already run out
through the whole domain during the 44,000 marching steps that followed the
first failure:

| r [R_p] | x2 entry | x2 exit | fc(H2) entry | fc(H2) exit | share of the drift |
|---|---:|---:|---:|---:|---:|
| 1.150 | 0.9942 | 0.9400 | 8.11e-1 | 7.67e-1 | 0.05 |
| 1.302 | 0.8362 | 0.9994 | 6.32e-1 | 7.56e-1 | 0.14 |
| 1.500 | 0.6399 | 0.9921 | 4.21e-1 | 6.53e-1 | 0.27 |
| **1.795** | 0.4349 | 0.8736 | **2.52e-1** | **5.07e-1** | **0.296** |
| 2.204 | 0.3042 | 0.7337 | 1.64e-1 | 3.96e-1 | 0.27 |

At the worst cell H2 holds 25% of the gas before the pass and 51% after. The
"relative drift in a place where H2 is 1e-5" reading of the P49 note is
**refuted**: there is no such place left in that state. The wind is molecular
to the top of the domain, which is the runaway itself, not a measure artefact.

### 2.4 The carrier convergence measure cannot see any of this

`carrier_steady_residual` divides by `carrier_row_scale = n_H(|v|+c_s)/dr`.
Section 128.5 argues correctly that `n(H2)` cannot be its own scale. But the
scale chosen is not the row's own largest term either, and section 133 made
that -- "a bound on that row's own largest term" -- the principle for every
hydrodynamic row. The transport Newton already computes the right quantity,
`dsc`, the sum of the row's term magnitudes. On the accepted wind:

| r [R_p] | reported (`n_H(|v|+c_s)/dr`) | on the row's own terms | dsc / rowscale |
|---|---:|---:|---:|
| 1.050 | 1.49e-7 | 3.63e-3 | 4.1e-5 |
| 1.150 | 5.26e-5 | 3.48e-2 | 1.5e-3 |
| 1.181 | 1.19e-4 | 6.79e-2 | 1.8e-3 |
| 1.200 | 1.59e-4 | **1.33e-1** | 1.2e-3 |
| 1.241 | -- | **6.04e-1** (max) | -- |
| 1.500 | 3.97e-7 | 8.53e-2 | 4.7e-6 |

Reported maximum 1.59e-4; on the row's own terms 6.04e-1, and above 0.1 in
52 cells. Two further points belong with it.

* The own-term measure has **no floor in this evaluation**, because the
  absolute floor `1e-20 n_element/dt` that `carrier_residual` carries is
  divided by the `dt_big = 1e30` that `carrier_steady_residual` substitutes to
  kill the time term. That is why the measure reads 0.55 at r = 1.25 where
  n(H2) is 2.3e4 cm^-3. Adding a floor built as a rate rather than from `dt`,
  `eps n_H (|v|+c_s)/dr`, cleans the far wind and **does not change the
  verdict**: at eps = 1e-8 the front still reads 3.5e-2 at 1.15 and 1.33e-1 at
  1.20, and 51 cells still exceed 0.1.
* So the transported front region of the state the solver accepts is 3 to 13
  per cent out of carrier balance on the row's own terms, and 26 to 60 per
  cent where H2 is trace. This appears to be the cleanest available statement
  of why the pass then moves so far.

## 3. How the carriers reach `||R||`

The chain is short and it closes exactly.

**Carriers -> particle count.** Replacing 2 H by 1 H2 halves the hydrogen's
share of `n_tot`. Per H nucleus the particle count is `1 - x2/2`.

**Particle count -> temperature.** The pass holds the wind, so `p` is the
`p` the JFNK solved; `T = p/(n_tot k)` then rises by the inverse ratio.
Predicted from the table of 2.1 against what `ioniz_eq` reports on the next
sweep (the omega = 0.5 state):

| cell | r [R_p] | x2 entry -> exit | particle ratio | T predicted [K] | T reported [K] |
|---:|---:|---|---:|---:|---:|
| 253 | 1.1810 | 0.2241 -> 1.0000 | 0.563 | 1915 | **1915** |
| 254 | 1.1833 | 0.2084 -> 1.0000 | 0.558 | 1958 | **1961** |

The entry values are 1078 and 1093 K. **The front heats by a factor 1.78 from
the particle count alone**, before any chemistry or radiation responds.

**Temperature -> the energy row.** The next JFNK starts at
`||R|| = 6.745e-2` with below-escape rows mass 2.595e-2, momentum/gravity
2.772e-3, energy 6.745e-2, and its worst cell is **j = 253, r = 1.1810, row
k = 3** -- the cell of the table, the energy row. The energy row exceeds the
mass row by 2.6 and the momentum row by 24. So the answer to "which row, and
where" is: the energy row, inside the region the carrier pass just
molecularized, and it is there because the pass changed how many particles
carry the pressure the wind was solved with.

## 4. The designs, measured

All from the same accepted wind (`picard/cap`, `||R|| = 1.454e-4`,
`info = 0`), same source, same input, `OMP_NUM_THREADS` 4 or 8 (the numbers
below do not depend on it). "passes at `info = 0`" counts consecutive outer
passes whose steady solve was accepted; the outer loop stops at 20.

| id | design | drift of pass 1 | relaxation steps | `\|\|R\|\|` handed to the next solve | that solve | passes at `info = 0` | H2 front, first -> last pass |
|---|---|---:|---:|---:|---|---:|---|
| `A1` | **baseline**: relax to the fixed-wind steady state, blend at omega = 0.5 | 0.553 | 63 (converged) | 6.745e-2 | `info = 2`, abort | 0 | -- |
| `B1` | (own) background (n_tot, T, ionization) re-solved at every relaxation step | 0.860 | 400 (not converged) | 4.762e-1 | `info = 2` | 0 | front to 1.43 |
| `C2` | (c) omega set from the drift, omega = min(0.5, 1e-3/drift) = 1.8e-3 | 0.553 | 63 | 2.031e-1 | `info = 2` | 0 | -- |
| `E6` | the same full relaxation, accepted WHOLE (no blend) | 0.553 | 63 | 1.925e-2 | `info = 2` | 0 | -- |
| `E5` | (e) trust region on the pass, 0.4 | 0.441 | 8 | 6.701e-4 | `info = 0` | 1 | 1.227 -> 1.247 |
| `E4` | (e) trust region 0.2 | 0.220 | 6 | 2.787e-4 | `info = 0` | 8+ | 1.227 -> 1.447 |
| `E3` | (e) trust region 0.05 | 0.053 | 3 | 2.372e-4 | `info = 0` | 19 | 1.210 -> 1.464 |
| `E1` | (e) trust region 0.01 | 0.0113 | 1 | 2.394e-4 | `info = 0` | **20 of 20** | 1.203 -> 1.250 |
| `E2` | (e) trust region 0.001 | 0.0113 | 1 | 2.394e-4 | `info = 0` | **20 of 20** | identical to `E1` |

Option (a) is absent from the table because section 2.3 removes its premise:
the measure it would replace is already absolute, and the cells it flags carry
25 to 100 per cent of the hydrogen.

Option (d) -- the carrier row inside the JFNK unknowns -- was **not built**;
section 6 records what the measurements say about its cost, which is the only
part of it this campaign can speak to.

Reading the table:

* **The under-relaxed blend is the single most destructive ingredient, and it
  is worse than either state it interpolates.** Carrier steady residual on the
  reported scale: entry 1.59e-4, undamped relaxed exit 2.17e-3, omega = 0.5
  blend **8.63e-3**. Accepting the same relaxed state whole (`E6`) leaves
  `||R|| = 1.925e-2` instead of 6.745e-2, a factor 3.5 less damage, and its
  solve reaches 4.119e-3 instead of aborting at 4.285e-2. A blend of two
  states is a solution of neither, and the ionization sweep and the residual
  both notice.
* **A smaller omega does not help; it hurts** (`C2`, 2.031e-1 against
  6.745e-2). Reducing omega by 276 while leaving the pass at the clamped
  ceiling state makes the blend worse, not better, so option (c) is not
  indicated on its own.
* **Re-solving the frozen background makes the pass move further, not less**
  (`B1`, drift 0.860 against 0.553; the front reaches 1.43 R_p in one pass).
  Section 128.6's statement that the frozen background is the conservative
  side of the feedback is confirmed, on the state P49 produced.
  `B1` is not without value: at a LATER hand-off, where the state is better
  relaxed, two consecutive passes give drifts 0.136 then 0.0186, a factor 7
  fall, and the steady solve between them returns `info = 0` -- the only
  contraction of the drift seen anywhere in this campaign. The solve after the
  second pass then fails (`info = 2`, 4.017e-3), so the contraction is not
  carried through, but it is the one sign that a self-consistent background
  helps once the state is close enough.
* **The trust region is what makes the alternation stable.** Bounding how far
  one pass may move the composition and handing over that TRANSPORT state --
  not a blend, not the clamped steady state -- keeps the steady solve at
  `info = 0` up to a bound of about 0.2 per pass. At 0.4 it survives one pass
  and fails on the second; at 0.55, the full relaxation, it fails at once.
* **And it does not converge.** In `E1` the drift creeps 1.13e-2 -> 1.34e-2
  over 20 passes and `||R||` 2.39e-4 -> 3.14e-4, while the front advances one
  cell in each pass, monotonically. Each pass is one transport step at the
  cell's own advective crossing time -- measured, `dt_phys` equals `dr/|v|`
  exactly at every front cell, 3706 s at 1.15 R_p and 799 s at 1.227 -- and
  0.047 R_p (1.65e8 cm) in 20 such steps is **1.03e4 cm/s against a gas speed
  of 1.24e4 cm/s** there. The front is being advected at 83 per cent of the
  gas speed, with no sink holding it. The
  bigger bounds do the same thing faster and no differently -- 1.464 R_p in 19
  passes at 0.05, 1.447 in 8 at 0.2 -- which is what a time integration looks
  like and not what a fixed-point iteration looks like.

## 5. Recommendation

1. **Do not attempt the alternation from a state whose carrier rows are not
   balanced, and measure that on the row's own terms.** The gate that would
   have refused this hand-off exists in `carrier_residual` already: `dsc`.
   The change is to `carrier_steady_residual` -- divide by
   `max(dsc, eps n_H (|v|+c_s)/dr)` instead of `carrier_row_scale`, with the
   floor written as a rate so that it survives `dt -> infinity`, and keep
   `carrier_row_scale` for nothing. On the accepted wind this reads 0.60
   rather than 1.6e-4, and would say plainly that the hand-off is not a
   carrier steady state. This appears to be the single highest-value change
   available and it is small.
2. **Replace the "relax to the fixed-wind steady state, then under-relax the
   blend" pass with a bounded transport advance, handed over whole.** Measured
   above: it converts an immediate abort into 20 consecutive accepted solves.
   A bound of 1e-2 is safe with margin; 0.2 is where the margin runs out. The
   under-relaxation should go with it -- it is measurably the worst part of
   the present pass, and the trust region is the control it was standing in
   for. Note that section 128.6 rejected a bound written in FLOW TIMES for a
   good reason and that reason does not reach a bound written in COMPOSITION:
   the flow-time bound could not be made small enough to matter (0.1 flow
   times still gave 0.533), while a composition bound of 0.01 is met by the
   first step of the operator.
3. **Drop the drift gate of 1e-3, or state what it means.** It is below the
   smallest move the transport operator can make in one step (1.13e-2 at this
   state), so it can only ever be met when the front has stopped moving --
   which is the condition sections 128.6 and 134.4 already name as the one to
   wait for. Keeping a gate that cannot fire reads as a control to the next
   reader. If a gate is wanted, the carrier residual of recommendation 1 is
   the one that means something.
4. **`Molecular carrier transport` should stay default off.** Nothing here
   changes the standing rule. What is added is that the alternation, in the
   form it has, cannot be the thing that flips it: the front's motion is what
   blocks it, and a Picard loop is a poor integrator of a propagating front.

## 6. On option (d), the carrier row inside the JFNK

Not built, so this is an estimate and it is labelled as one. Section 128
rejected it on cost. Two measurements bear on that.

* What one residual evaluation costs would barely move. `eval_residual`
  already runs the ionization sweep, which is the expensive part; the carrier
  residual is one assembly over 500 cells with a 4x4 chemistry block and no
  root find.
* The cost is in the banded Jacobian's colouring (17 colours now) and in
  `block_thomas`'s block size. That is the part that has not been measured
  here, and it is what a re-evaluation would have to measure.

What the campaign does say is that the reason to want (d) is stronger than the
cost argument makes it look: the two halves are coupled through the particle
count with a gain of 1.78 in T at the front (section 3), and a splitting that
evaluates each half at the other's old state is being asked to work outside
the radius where that is reasonable. But it also says something (d) cannot
fix: **if the front is genuinely propagating, a coupled Newton has no steady
state to find either.** Section 134.4's marching result -- the front still
drifting at 0.0033 R_p per 20,000 steps at 150,000 steps -- and this
campaign's alternation, which advects the front at the gas speed with no sink
holding it, point the same way. Settling whether a steady transported state
exists on this planet is prior to choosing a solver for it.

## 7. What was NOT measured

* Whether the front stops at all. No long marching run was made here; section
  134.4's is the reference and it does not settle it.
* Option (d). Not built; section 6 is an estimate.
* Any planet but the He/H = 0.0793 hot Uranus of `E1d_fix2`, and any state but
  its two hand-offs.
* The regression matrix. Nothing in this campaign is a candidate for the tree
  yet, and the scratch source additionally carries the P43 copy's
  `input_read.f90`, which still requires `Log10 lower boundary` and cannot run
  six of the nine cases (P49 section 8 records this).
* `run_fcheck.sh`.
* Whether the recommendation 1 measure changes the verdict on any state that
  IS converged. It has only been evaluated on states that are not.

## 8. Runs and files

Source (scratch only): `.../scratchpad/picard/src`, binaries `EXHALE.x`
(capture), `EXHALE_v5.x` (probe), `EXHALE_v7.x` (all design options).

Additions to the scratch source, all gated by an environment variable and
inert without it:

* `EXHALE_main.f90` -- `EXHALE_PICARD_DUMP=1` (write the accepted wind and
  stop); `EXHALE_PICARD_REFRESH=<M>`, `EXHALE_PICARD_TRUST=<d>` and
  `relax_carriers_with_background`; `EXHALE_PICARD_OMEGA_AUTO=<target>`.
* `diffusive_photochemistry.f90` -- `EXHALE_PICARD_PROBE=1` (the pass dump and
  the relaxation trace), the saved `pct_res_rel`, `pct_row_scale` and
  `pct_dsc`, and the two helpers `carrier_relax_timestep` and
  `carrier_mixing_ratios` that let the loop be driven from the main program.

| directory | what it is |
|---|---|
| `picard/cap` | the capture run; `output/{Hydro_ioniz,Ion_species}.txt` IS the accepted wind every experiment starts from |
| `picard/A0` | baseline, full sequence with the probe on; `output/carrier_pass_02.txt` is the late pass of section 2.3 |
| `picard/A1` | baseline from the accepted wind; `output/carrier_pass_01.txt` is the table of section 2.1 |
| `picard/B1` | background re-solved every relaxation step |
| `picard/C1`, `picard/C2` | omega from the drift, starting at 0.5 and at 1.81e-3 |
| `picard/E1` ... `picard/E6` | trust region 1e-2, 1e-3, 5e-2, 0.2, 0.4, 1.0 (the last never fires, i.e. the full relaxation accepted whole) |

---

# P50 part two. The carrier row inside the Newton (2026-09-03)

The judgment of sections 0 to 5 was accepted and option (d) adopted: the wind
and the carriers are solved together. What follows is the implementation, what
it cost, and what it found. **Nothing below is a converged solution and
nothing below is quoted as a rate.**

## 9. Judgment

**The splitting was not merely slow -- it was pointing the wrong way.** The
alternation moves the H2 front OUTWARD, about one cell in each pass, and runs
away. The coupled Newton, from the same state, moves it **INWARD**:
`f(H2) = 0.5` from 1.1421 to 1.0618 R_p, `f(H2) = 1e-2` from 1.2192 to 1.1488.
The reason is the feedback section 3 measured: H2 costs particles, fewer
particles at the wind's own pressure means a hotter cell, and a hotter cell
destroys H2. A splitting that holds the wind fixed cannot see that, and the
sign of the front's motion is what it gets wrong.

**The coupled solve does not converge either, but it fails somewhere else and
much further in.** It reaches `||R|| = 1.532e-3` against a tolerance of 1e-3 --
a factor 1.5 short, from 1.221e-2 -- and drives the carrier row's worst cell
from 0.604 to 0.044 of its own terms. It ends on the iteration limit
(`info = 1`), never on an abort. What stops it is the line search circling, not
a missing direction and no longer the element budget.

**Three defects were found in the coupled formulation, each by a measurement,
and each is recorded because each looked plausible before it was measured.**
They are in section 11.

## 10. What was built

`Coupled carrier solve: True|False`, default false, and meaningless without
`Molecular carrier transport: True`.

* **Unknowns.** Four per cell, `(rho, rho v, E, n(H2))`, the last in the code's
  own density units `n(H2)/n0` -- not in cm^-3. The Krylov step size is
  `sqrt(eps)(1 + ||Y||)/||v||`, a norm of the whole vector, so a component
  thirteen decades from the others would set that step by itself.
* **The carrier row** is the steady continuity equation the transport operator
  already assembles, with the time term off: `carrier_steady_residual` returns
  it cell by cell and `eval_residual` packs it as row 4, converted to the
  code's time unit so that the merit compares one clock.
* **The chemistry coupling is the one that was already there.** With the
  carrier transport on, `ioniz_eq` holds H2 fixed at the partition it is handed
  (`x_h2_fix`) and solves everything else around it. The coupled solve changes
  only where that number comes from: the Newton unknown instead of the last
  transport step. So the particle count, the temperature `p/(n_tot+n_e)` that
  follows from it, and every heating and cooling term are functions of the
  carrier unknown, and the Jacobian sees it. **No new chemistry was written.**
* **The limiter is frozen inside the model**, on `weno_mode` -- the switch
  section 121 already uses for the WENO3 weights, and for the same reason: the
  limited slope reaches `j-2` and `j+2`, which the banded stencil does not
  carry. Mode 2 (Jacobian columns, Krylov products) reuses the slopes of the
  outer iterate; modes 0 and 1 rebuild them, so the marching path and the line
  search see the slopes of the state they are evaluating.
* **The element budget is a CONSTRAINT, not a row.** A trial whose `n(H2)`
  leaves the hydrogen its cell has is refused where `admissible` is asked for,
  the way the positivity of rho and p is refused in the line search. It is not
  clamped: a clamp puts a kink in the residual, which is the defect section 126
  spent a campaign on, and it is exactly what made the fixed-wind relaxation
  return a state pinned on the ceiling in a quarter of the grid.
* **Band geometry** is set at run time by `set_carrier_unknown`:
  `|row-col| <= 4*2 + 3 = 11`, so `kl = ku = 11` and 23 colors, against 8 and
  17 for three unknowns. The carrier row itself reaches only +/-1 -- its
  second-order correction is a frozen source -- so the WENO3 reach still sets
  the bandwidth and the fourth variable only widens the stride.
* **Acceptance.** The carrier row is one more row of `||R||`, combined by the
  maximum like the others, on its own volume-weighted two-region measure. The
  flux gate of section 133 **cannot be evaluated on this source**: it is in the
  repository tree and not in the P43 working copy this is built on, which
  carries the section 127/131 residual measure alone. That limitation is
  P49's, unchanged.
* **`solve_steady_ptc` refuses the coupled unknown** and says so. It builds its
  Jacobian from `frozen_residual`, which does no equilibrium sweep and so has
  no carrier row; giving it one would mean assembling that row from a
  background nothing in that routine refreshes.

## 11. Three defects of the coupled formulation, each measured

**(i) One scale vector cannot scale both the carrier row and its unknown.**
The solver scales rows and columns with the same `D` (`D^-1 J D`), which is
right when the residual has the units of the unknown over a time -- true of
every row here. It is still not enough, because the carrier row is STIFF: its
natural rate is the chemical one, about 1e3 per code time at the front against
the hydrodynamic rows' 1e-6.

| carrier scale used for both | merit at the hand-off | outcome |
|---|---:|---|
| `n(H2) + row_terms/((\|v\|+c_s)/dr)` (state-like) | **89.8** | abort at iteration 2 |
| `row_terms * (R0/v0)` (the row's own terms) | 2.51 | abort at iteration 2 |
| **row terms for the ROW, state-like for the COLUMN** | 2.51 | runs, `||R||` 1.221e-2 -> 1.532e-3 |

With the state-like scale on the row, the line search sees nothing but the
carrier and the damped Gauss-Newton escape lands on a state whose chemistry
cannot be certified. With the row scale on the column, the Newton proposes
carrier steps **84 times the cell's own `n(H2)`** and the line search cuts to
`lam = 4.9e-4`. The two-sided scaling is the fix: `Drow` beside `D`, identical
for the three hydrodynamic rows -- which is what keeps every run without the
carrier unknown bit-for-bit unchanged -- and different only for the carrier.
The pseudo-transient term is then `diag(D/Drow)` rather than the identity, and
is written that way.

**(ii) The constraint was moving with the trial it constrained.** The element
headroom depends on the very unknown it bounds (`n_H_free` counts H2's own two
nuclei). Refreshed at every state that is not a Jacobian probe, the last
LINE-SEARCH trial leaves its own headroom behind, and the Gauss-Newton escape
that follows finds the ITERATE outside it: every backtracked step is refused
down to `lam = 5e-4` and the solve aborts holding a descent direction. It is
now frozen at the outer iterate (`weno_mode = 1`), and the test is a
COMPARISON against how many cells the iterate itself leaves outside -- the same
rule, and the same reason, as section 138's chemistry gate: a rule that rejects
the neighbourhood of the point it stands on cannot be used to leave that point.
A negative density has no such excuse and is refused outright.

**(iii) The Gauss-Newton escape had no line search, and with the carrier row it
needs one.** The escape used to evaluate the full step and nothing else, which
suffices while the only way it can fail is to ascend -- the direction is built
to descend, so a shorter one descends too. With the carrier unknown a second
failure mode appears and it is not about descent: the full step LOWERS the
merit (2.28 -> 1.50, measured) and puts 73 cells' molecular equilibrium outside
the physical simplex, so the trial is refused and the solve aborts with a
descent direction in hand. It now backtracks up to 12 times. **The three-row
route keeps the single evaluation it has always had**, so no atomic solve
moves.

## 12. What the coupled solve does, measured

From `picard/cap` -- the wind the alternation's first solve accepted at
`info = 0`, `||R|| = 1.454e-4` -- with `Coupled carrier solve: True`:

| | at entry | at exit (best iterate) |
|---|---:|---:|
| `\|\|R\|\|` (all four rows) | 1.221e-2 | **1.532e-3** |
| carrier row, worst cell / own terms | 6.04e-1 | **4.41e-2** |
| carrier row, volume-weighted | 7.72e-3 | 1.25e-2 |
| trials refused for the element budget | -- | 38 in 20 iterations |
| `info` | -- | 1 (iteration limit) |

`||R||` at entry is 1.221e-2 rather than 1.454e-4 because the carrier row is
now in it: the wind that solve accepted was steady in three rows and not in the
fourth, which is the whole finding of section 2.4 restated by the solver
itself.

**The front, and the three layers.** The state the solve returns against the
state it was handed. **Neither is converged** -- the returned state is 1.5
times the residual tolerance and its carrier row is 12 times it -- so this is a
direction and a size, not a result:

| | hand-off | coupled solve |
|---|---:|---:|
| `f(H2) = 0.5` at | 1.1421 | **1.0618** |
| `f(H2) = 1e-2` at | 1.2192 | **1.1488** |
| `f(H2) = 1e-4` at | 1.2300 | 1.2327 |
| `log10 Mdot` | 10.413 | 10.412 |

| r [R_p] | T [K] hand-off -> coupled | f(H2) | n(H3+) | n_e |
|---|---|---|---|---|
| 1.05 | 711 -> 867 | 0.943 -> 0.886 | 9.7e4 -> 3.7e4 | 1.18e6 -> 2.60e6 |
| 1.15 | 871 -> 766 | 0.445 -> **9.8e-3** | 8.5e2 -> 6.6e0 | 7.17e7 -> 8.29e7 |
| 1.25 | 1567 -> 1562 | 8.1e-6 -> 5.9e-5 | 7.9e-4 -> 5.7e-3 | 6.74e7 -> 6.74e7 |
| 1.50 | 2538 -> 2538 | 5.4e-5 -> 2.5e-5 | 1.3e-3 -> 6.2e-4 | 6.47e7 -> 6.48e7 |
| 2.00 | 2691 -> 2693 | 1.8e-5 -> 7.6e-6 | 4.7e-5 -> 1.9e-5 | 4.43e7 -> 4.43e7 |
| 3.01 | 2206 -> 2211 | 4.1e-6 -> 2.4e-6 | 1.8e-6 -> 1.1e-6 | 1.87e7 -> 1.87e7 |

The wind above 1.25 R_p barely moves and the rate does not move at all
(10.413 -> 10.412). **Everything happens between 1.05 and 1.20 R_p**, which is
the molecular layer and the front, and it happens in the direction the
alternation could not go.

## 13. Cost

Both from the same hand-off, `EXHALE_JFNK_MAXIT=20`, `OMP_NUM_THREADS=8`, same
binary:

| | 3 unknowns (the alternation's wind solve) | 4 unknowns (coupled) |
|---|---:|---:|
| colors | 17 | 23 |
| outer iterations used | 6 (converged) | 20 (limit) |
| full residual evaluations | 64 | 778 |
| wall | 1.00 s | 14.51 s |
| result | `info = 0`, `\|\|R\|\| = 1.454e-4` | `info = 1`, `\|\|R\|\| = 1.532e-3` |

Per outer iteration that is 11 residual evaluations against 39. Six of the
extra come from the wider colouring (23 against 17); the rest are the damped
Gauss-Newton escape, which now backtracks and which fires in almost every
iteration. **Section 128's cost objection is not what stands in the way**: on
this state the alternation spends about 1 s in each of 20 passes and does not
converge, and the coupled solve spends 14.5 s and reaches a residual 1.5 times
the tolerance with a qualitatively different front. The two are within a factor
of order one of each other.

## 14. Where it stops, and what that says

The solve reaches its best iterate at about iteration 10 and then circles for
the rest of its budget: `||R||` 1.53e-3 at its best, 4.4e-3 where it plateaus,
61 non-monotone acceptances in 500 iterations. It never aborts.

That is the Grippo non-monotone line search accepting sideways steps
(`steady_newton.f90:1682`, the test against the worst merit of the last five
iterates) -- the pathology section 126 identified, which here is not fatal
because the solve keeps and returns its best iterate. The measure that says so:
392 of 501 reported worst cells are the carrier row and 108 the energy row, and
the carrier row's own worst cell falls by a factor 14 while the volume-weighted
number does not fall at all. The residual stops being concentrated at the front
and becomes spread, and the search has nothing sharp left to follow.

**So the honest reading is that the coupled route removes the obstruction the
alternation had and exposes the one section 126 left.** What it does NOT settle:

* whether a converged coupled steady state exists on this planet. The solve got
  within 1.5 of the tolerance and stopped improving; that is not evidence
  either way.
* whether the inward front is where a converged solve would put it. The
  direction is measured and the size is not.

## 15. What was checked, and what was not

* **Both parts are byte-identical on the repository's own regression matrix,
  and part two had to be made so.** `make check` in the tree, ten cases --
  `wasp_full`, `wasp_he23off`, `mol_base_handoff`, `mol_metals`,
  `mol_lyman_werner`, `mol_diffusion`, `mol_ir_bands`, `mol_sec_ion`,
  `mol_carrier`, `lower_profile` -- run twice: **10/10 byte-identical with
  part one alone, and 10/10 again with the coupled route on top of it.**

  Part two did NOT start out that way. `wasp_full_newton` and
  `wasp_he23off_newton` moved by 1.2e-9 and 4.0e-8 relative, with `log10 Mdot`
  unchanged to six decimals and, on the first, an identical JFNK iteration
  log. The cause was the two-sided row scaling: `sum((F/Drow)**2)` in place of
  `sum((F/D)**2)` is the same arithmetic when `Drow` holds `D`'s values, but
  not the same instruction sequence, and a reduction summed in a different
  order differs in its last bits. **The run-time band geometry is NOT the
  cause and that was tested**: a build with `nvar_jac`, `kl_jac`, `ku_jac` and
  `ncolor_jac` restored to compile-time parameters is bit-identical to the
  run-time build (`picard/reg/wasp_full_newton_p2c`). Making `Drow` a local of
  the solve did not restore identity either. What did was branching every site
  that divides by the row scale on `nvar_jac` and keeping the ORIGINAL
  expression, character for character, in the three-unknown branch -- plus
  allocating `Drow` last and on its own so that every array the three-unknown
  solve has always had keeps its address.

* **Six of the nine regression cases still cannot run on this source**, before
  and after: the P43 copy's `input_read.f90` requires `Log10 lower boundary`
  where the tree has made it optional, and `mol_base_handoff`, `mol_metals`,
  `mol_lyman_werner`, `mol_diffusion`, `mol_ir_bands` and `mol_sec_ion` carry
  that key in `base.inp`. This is P49's limitation, unchanged; the tree's
  `input_read.f90` was not copied in, so the molecular matrix is untested here.
* `run_fcheck.sh` was not run.
* The coupled route on any planet but the He/H = 0.0793 hot Uranus, on any
  state but this hand-off, and with helium diffusion also on (the outer loop
  falls through to the helium relaxation in that case and nothing exercises
  it).
* Whether the band half-width 11 is exact rather than sufficient. It is derived
  (`4*2 + 3`), and no column-by-column check against a dense Jacobian was made.

## 16. Files and runs (part two)

Source: `.../scratchpad/picard/src` (both parts). The P49 copy it is a change
to is `.../scratchpad/picard/src_p49base`, and **part one alone is kept as a
separate source and a separate diff**, because it goes to the tree whether or
not part two does:

| file | what it is |
|---|---|
| `picard/src_p49base` | the P49 working copy, unchanged |
| `picard/src_part1` | part one applied to it -- 473 diff lines, two files |
| `picard/src` | part two applied to that -- 1398 further diff lines, five files |
| `picard/p50_part1.diff` | `src_p49base` -> `src_part1` |
| `picard/p50_part2.diff` | `src_part1` -> `src` |

Binaries: `EXHALE_p1b.x` from `src_part1`, `EXHALE_p2.x` from `src`. The
reconstructed part-one source was checked against the part-one run the
measurements were made on and reproduces it line for line (`picard/P1b`
against `picard/P1`).

**Part one** -- `diffusive_photochemistry.f90` (the rate-converted absolute
floor, `row_terms`, `carrier_steady_residual` on the row's own terms with the
volume-weighted measure, the composition trust region in place of the blend),
`EXHALE_main.f90` (the acceptance test on the carrier residual rather than on
the drift). It goes to the tree whether or not part two does.

**Part two** -- `steady_newton.f90` (the runtime band geometry, `pack_carrier`
/ `unpack_carrier`, the carrier row in `eval_residual`, `Drow`, the headroom
constraint, the backtracked escape, the cost line), `diffusive_photochemistry.f90`
(`carrier_column_scale`, `carrier_row_scale_H2`, `carrier_headroom`, the frozen
limiter), `parameters.f90` and `input_read.f90` (`Coupled carrier solve`),
`EXHALE_main.f90` (the call, and the outer loop running once).

| directory | what it is |
|---|---|
| `picard/P1` | part one alone from the accepted wind: 20 outer passes, all `info = 0` |
| `picard/D1`, `D2` | the coupled solve with the two single-scale choices of section 11(i) |
| `picard/D4`, `D5`, `D6` | with the headroom fixed, then with two-sided scaling, then at a 4000-iteration budget |
| `picard/D7` | the state quoted in section 12, dumped at the solve's return |
| `picard/C20`, `U20` | the cost pair of section 13 |
| `picard/SC271` | the single-cell carrier row scan (`EXHALE_CARRIER_SCAN`) |
| `picard/reg/*_{base,p1,p2}` | the byte-identity matrix of section 15 |
| `picard/reg/wasp_full_newton_p2c` | part two built with the band geometry back at compile time, which decides that the geometry is not what moved the atomic cases |
| `picard/src_paramtest` | the source of that build |
