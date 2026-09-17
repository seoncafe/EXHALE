# LHS 1140 b, item L5b: why the element relaxation refuses every pass

Item L5b of `docs/PLAN_20260913_lhs_stationary.md`, opened by L5a's finding
that all six seeds end at outer pass 1 with `element_step_out_of_bounds`
after 24 to 26 steps, where the HD 209458 b control reaches the fixed point
of the same operator in 30 steps at every pass.

This item is a diagnosis. **No source file was changed.** Every number is
MEASURED on the binary of L5a (`<scratch>/L5a_bin/EXHALE.x`, md5
`97e10317a710b9ccc63addbedde3586a`, newer than every file under `src/`, so
it is a build of the current tree) unless marked READ.

## 1. Verdict

**The bound that fails is the upper one on the helium mass fraction, and it
fails because the FIXED POINT of the element operator on the wind it is
given sits exactly ON that bound: the He/H partition the operator relaxes to
is pure helium over a band of cells. Nothing about the step control can
reach a fixed point that is not in the interior of the admissible set.**

The reason the fixed point is there is measured, not inferred. At the fixed
point the element row makes the helium elemental flux
`F_He = 4 pi r^2 (rho v X + J)` radially constant, so
`X = (F_He - 4 pi r^2 J)/F_rho` with `F_rho = 4 pi r^2 rho v` the state's own
total mass flux. On the LHS 1140 b state the hydrodynamic solve hands over,
`F_rho` falls 27 percent between 1.36 and 26.9 R_p while `F_He` falls 9
percent, so that ratio climbs from 0.85 to **1.000000 at r = 13.9 R_p**
(table T3), where the hydrogen elemental flux passes through zero
(`F_H = 6.1e+03 g/s` at cell 450, `-2.5e-02` at cell 455, against
`8.1e+07 g/s` of total). On the HD 209458 b control the same ratio is
0.199145 flat to 2e-5 over 1.30 to 3.90 R_p (table T4) and no bound is
approached.

So the answer to the item's three candidates is:

- **(ii), the state the hydrodynamic solve hands over, is the cause.** Its
  discrete mass flux is not radially constant; the element row, which is
  conservative in the elemental flux at a density the hydrodynamics owns,
  converts that into a composition and there is no composition below X = 1
  that it can converge to.
- **(i), the bound definition, is correct and is not the cause, but the
  helium-rich base makes the margin 5.6 times smaller.** `X_base` is
  MEASURED at 0.865918 here against 0.245935 on HD 209458 b, so the headroom
  to the upper bound is 0.134 against 0.754. The same relative excursion of
  the partition reaches the bound here and not there.
- **(iii), the step control, is the proximate refusal and cannot be the
  cure.** The excursion above 1 falls LINEARLY with the step length
  (measured factor-2 sequence, table T2), and the retry control halves the
  step at most `he_relax_retry_max = 20` times, i.e. it can divide the
  excursion by 1e6, while clearing an excursion of 5e-3 to the bound
  tolerance `element_fraction_bound_tol = 1e-12` needs 33 halvings. Even
  with an unlimited budget the limit of the halving sequence is the entry
  composition, not an advance.

**Two knobs the item named do nothing here, by code and by measurement.**
`omega` (`EXHALE_DIFF_OMEGA`) damps the blend applied AFTER the relaxation
returns (`binary_element_diffusion.f90` line 1264, READ), so a pass whose
relaxation never returns a relaxed composition is untouched by it: at
omega 0.125 the pass is refused at the same step count and the same step
length as at 0.500, with every hydrodynamic row identical to the printed
digits (table T5). `EXHALE_CARRIER_TRUST` bounds the molecular carrier
relaxation, which is entered only under `thereis_mol .and.
carrier_transport` (`EXHALE_main.f90` line 6529, READ); this configuration
is atomic (`carrier_transport F`, `EXHALE_resolved.out` line 19, READ), so
the knob is never read.

**The element relaxation is not the binding obstruction of the route.**
With `EXHALE_OUTER_PASSES=1` the outer loop takes no composition update at
all (`EXHALE_main.f90` line 6466: an update is taken only while
`it_diff < pass_cap`), so the pass is a pure hydrodynamic solve at the
loaded composition. The state it produces is still refused, by the
hydrodynamic rows: mass 1.075 (worst refusing entry 3.236e-3 against
3.018e-12 at cell 224, r = 1.2269), momentum 2.643e-2, energy 1.842
(table T5). Removing the element relaxation therefore does not produce a
certifiable state; it only moves the refusal to the equation L5c is about.

**What DOES get the route past pass 1 is `Restart intent: stationary
equilibrate`,** and not for the reason the item supposed. That second word
is the ionization-equilibrium sweep of the loaded composition
(`equilibrate_loaded_composition`, `EXHALE_main.f90` line 5307, READ), NOT
the element fixed point; there is no entry point in the code that relaxes
the elements on the loaded state without a hydrodynamic solve. It
nevertheless removes 5.06e-2 of particle-count inconsistency from the seed
(MEASURED, the run's own line: "restart composition equilibrated in 9
sweeps; max |d(n_tot+n_e)|/(n_tot+n_e) removed = 5.063E-02, left =
5.259E-11"), and on the state that follows the element relaxation reaches
its fixed point at EVERY pass of the run: 30, 31, 30 and 32 steps with the
drift falling 3.98e-1, 1.94e-1, 9.83e-2, 7.41e-2, while the hydrodynamic
rows fall from mass 1.78, energy 1.47 to mass 1.27, energy 0.864 and then
hold (table T6). The run ends at pass 5 on the joint progress control, not
on an element refusal: "the joint distance of the state has not fallen in 3
consecutive passes", standing at 1.896e-4 against 1.0e-5. This is the first
LHS 1140 b run of the partitioned route that gets past pass 1, and the
first whose refusal names the equation that is actually stuck.

**Minimal change proposed: none in the element operator's bound test.** The
three proposals in section 6 are a diagnostic message that names the cell
and the bound, a retry budget sized against the tolerance it has to clear,
and one factual correction to a comment. None of them makes this planet
converge; the obstruction that remains is L5c's.

## 2. What `element_step_out_of_bounds` tests, and where it fails

READ, `src/modules/functions/binary_element_diffusion.f90`:

| what | line | text |
|---|---|---|
| the measure | 741-742 | `Xover = maxval(Xhe(1:N)) - 1`, `Xunder = -minval(Xhe(1:N))`, taken on the raw solve output BEFORE the range clip of lines 757-758 |
| the tolerance | 452 | `element_fraction_bound_tol = 1.0d-12` |
| the test | 882 | `else if (max(Xover, Xunder) .gt. element_fraction_bound_tol) then outcome = element_step_out_of_bounds` |
| the restore | 887 | `if (outcome .ne. element_step_accepted) f_sp = f_entry` |
| the name printed | 931-932 | outcome 3 reads `helium fraction left [0,1]` |

`Xhe` is the helium MASS fraction `m_He n_He/(sum_i m_i n_i)` of cells 1..N
(line 682). There is one element in this test: the trace metals carry their
own diagnostic (`trace_ratio_under_zero`) and no bound of their own, and
this run holds no metals (`metal/H dep 0.000E+00` on every diagnostic line,
MEASURED).

**Which bound, which cells, which radius** (MEASURED, `EXHALE_DIFFUSION_CHECK=1`
on seed S1c through the same route as L5a, `<scratch>/L5b_diag_S1c`):

- the UPPER bound, every time. `X_min` sits at 0.8386 for the whole
  refusing part of the relaxation and never approaches 0.
- the band `X_face > 0.9999` runs from cell 449 to cell 461, and the
  maximum `X_face = 0.99999993` is at face 455, **r = 13.91 R_p**
  (`output/element_flux_profile.txt` of the refused candidate, whose
  `X_face` column is `0.5(X_j + X_{j+1})` after the clip).
- the excursion at the first refusal is `Xover = 5.234e-03`, at the last
  `3.999e-07` (table T2), against a tolerance of 1e-12.

The radius depends on the state and not on the planet: on a state whose
mass flux decays much faster (the `EXHALE_JFNK_MAXIT=1` run, whose
`4 pi r^2 rho v` falls 2.90 percent per cell through the base layer) the
same runaway happens at **r = 1.052 R_p**, cells 140-141, and the refusal
comes after 15 steps instead of 24.

## 3. The trajectory, and why 24 to 26 steps

The relaxation grows its step geometrically: `grow` starts at 1, is
multiplied by 1.5 after each accepted step and halved after each refused one
(lines 1194-1225, READ), and the step handed to the operator is
`dt_code(j)*grow` with `dt_code(j) = min(dr^2/(D_12 + K_zz), dr/max(v, w_s))`
per cell. Because multiplication commutes, the step length the refusal
message prints is exactly `1.5^(accepted) x 0.5^(refused)`: the 24-step case
prints 1.605e-2 = `1.5^24 x 0.5^20`, the 25-step case 2.408e-2 and the
26-step case 3.612e-2, which is how L5a's table T8 reads. So every LHS
refusal spent the full budget of 20 retries; the other exit,
`grow <= he_relax_grow_floor = 1e-6`, was never reached.

### T1. The relaxation of pass 1 on seed S1c, full route (`<scratch>/L5b_diag_S1c`)

`grow` reconstructed from the acceptance sequence; `state` is the verdict
the bound test returns (REFUSE iff `Xover > 1e-12`).

| step | grow | state | X min | X max | X over 1 | Newton passes |
|---|---|---|---|---|---|---|
| 1 | 1.00e+00 | accept | 0.565815 | 0.865985 | -1.340e-01 | 14 |
| 5 | 5.06e+00 | accept | 0.566135 | 0.865954 | -1.340e-01 | 18 |
| 10 | 3.84e+01 | accept | 0.590503 | 0.865918 | -1.341e-01 | 10 |
| 13 | 1.30e+02 | accept | 0.690661 | 0.865918 | -1.341e-01 | 11 |
| 14 | 1.95e+02 | accept | 0.740656 | 0.865918 | -1.341e-01 | 7 |
| 15 | 2.92e+02 | accept | 0.809260 | **0.876757** | -1.232e-01 | 9 |
| 16 | 4.38e+02 | accept | 0.840150 | **0.934407** | -6.559e-02 | 7 |
| 17 | 6.57e+02 | accept | 0.838678 | **0.995492** | -4.508e-03 | 12 |
| 18 | 9.85e+02 | accept | 0.838611 | 0.998336 | -1.664e-03 | 6 |
| 19 | 1.48e+03 | **REFUSE** | 0.838511 | 1.000000 | **+5.234e-03** | 30 |
| 20 | 7.39e+02 | REFUSE | 0.838561 | 1.000000 | +9.283e-04 | 30 |
| 21 | 3.70e+02 | accept | 0.838585 | 0.999398 | -6.018e-04 | 6 |
| 24 | 1.39e+02 | accept | 0.838576 | 0.999796 | -2.043e-04 | 11 |
| 30 | 1.95e+01 | REFUSE | 0.838571 | 1.000000 | +8.624e-07 | 11 |
| 36 | 2.74e+00 | REFUSE | 0.838571 | 1.000000 | +1.451e-06 | 5 |
| 39 | 1.03e+00 | REFUSE | 0.838571 | 1.000000 | +3.999e-07 | 4 |

24 steps were accepted and 20 refused. Five of the 20 refusals left no
diagnostic line: a step whose Newton row did not solve returns before the
diagnostic is written (line 750, READ), so those are
`element_step_solve_failed`, not out-of-bounds; the last one is outcome 3,
which is what the message prints.

The shape is the whole story. Through step 14 the maximum is pinned at
`X_base = 0.865918` and only the MINIMUM moves: the operator is
redistributing helium in the wind. Between steps 15 and 19 the maximum
leaves the base value and runs 0.8768, 0.9344, 0.9955, 0.9983, 1.0000 (the
partition in the outer wind is being driven to pure helium) and from there
the maximum never comes back below 0.9994 whatever the step length.

### T2. The excursion is linear in the step length

From the `EXHALE_JFNK_MAXIT=1` run (`<scratch>/L5b_maxit1_S1c`), whose
runaway cell is in the base layer and whose refusal sequence is a clean
halving chain:

| step | grow | X over 1 | over / grow |
|---|---|---|---|
| 32 | 6.682e-03 | 1.510e-05 | 2.26e-03 |
| 33 | 3.341e-03 | 7.198e-06 | 2.15e-03 |
| 34 | 1.670e-03 | 3.250e-06 | 1.95e-03 |
| 35 | 8.352e-04 | 1.276e-06 | 1.53e-03 |

The predicted slope is the non-conservative source of the element row
carried over one step. The row's advective term is the divergence of the
state's own face mass fluxes, and its coefficients sum to `div(F_rho)`
rather than to zero (the donor-cell linearization, lines 2268-2270, READ):
at `X = 1` the drift flux `X(1-X)` has shut off and what is left is
`rho dX/dt = -X div(F_rho) cadvf`, i.e.
`dX = X |dlnMdot/dcell| dt/t_cross` with `t_cross = dr/v`. MEASURED on that
state at the runaway cell 140 (r = 1.0515): `dlnMdot/dcell = -2.903e-02`,
`dr = 1.037e+06 cm`, `v = 2.008e+02 cm/s`, `D_eff = 1.910e+09 cm^2/s`,
`K_zz = 1e+09`, so `t_diff = 370 s`, `t_cross = 5165 s`, `dt_code = 370 s`
and the predicted `dX/grow = 2.08e-03`. The measured 1.5e-3 to 2.3e-3 is
that number.

The consequence is the one in the verdict: to bring 5.2e-3 (the first LHS
refusal) below 1e-12 takes `log2(5.2e9) = 33` halvings against a budget of
20, and the limit of the halving chain is no advance at all.

## 4. What the fixed point is, measured

`F_He`, `F_H` and `Mdot_face` are the operator's own diagnostic columns,
written from the same face-averaged `rho_f`, `v_f` and the `X` and `J` of
the step (`write_element_flux_profile`, lines 3786-3806, READ), so their
ratio is an internally consistent measurement.

### T3. LHS 1140 b, the state pass 1 hands to the relaxation (refused)

| cell | r [R_p] | X_face | F_He [g/s] | F_H [g/s] | Mdot_face [g/s] | F_He/Mdot |
|---|---|---|---|---|---|---|
| 100 | 1.0252 | 0.845002 | 1.10024e+08 | 2.19552e+07 | 1.31979e+08 | 0.833647 |
| 150 | 1.0621 | 0.848046 | 1.01334e+08 | 5.26116e+07 | 1.53946e+08 | 0.658246 |
| 200 | 1.1502 | 0.901418 | 9.39149e+07 | 2.81770e+07 | 1.22092e+08 | 0.769215 |
| 250 | 1.3609 | 0.930568 | 8.79799e+07 | 1.51968e+07 | 1.03177e+08 | 0.852711 |
| 300 | 1.8646 | 0.934666 | 8.38790e+07 | 9.33777e+06 | 9.32168e+07 | 0.899827 |
| 350 | 3.0688 | 0.954084 | 8.27863e+07 | 4.32962e+06 | 8.71159e+07 | 0.950300 |
| 400 | 5.9475 | 0.988929 | 8.21051e+07 | 9.44535e+05 | 8.30496e+07 | 0.988627 |
| 450 | 12.829 | 0.999938 | 8.13610e+07 | 6.08894e+03 | 8.13671e+07 | 0.999925 |
| **455** | **13.906** | **0.9999999** | 8.12851e+07 | **-2.50741e-02** | 8.12851e+07 | **1.000000** |
| 470 | 17.764 | 0.999474 | 8.10560e+07 | 3.99753e+04 | 8.10960e+07 | 0.999507 |
| 495 | 26.919 | 0.996642 | 8.06691e+07 | 2.66105e+05 | 8.09352e+07 | 0.996712 |

`X_face` and `F_He/Mdot` agree to seven digits, which is the check that the
reading "the fixed point is `X = F_He/F_rho`" is what the operator is doing
(the diffusive term contributes the difference and is 1e-3 of the advective
one out there). Between 1.36 and 26.9 R_p `F_He` falls 9.4 percent and
`Mdot_face` 27 percent; the ratio has nowhere to go but up, and it reaches
1 where the hydrogen elemental flux reaches zero.

### T4. HD 209458 b, the same measurement on the certified route

| cell | r [R_p] | X_face | F_He [g/s] | F_H [g/s] | Mdot_face [g/s] | F_He/Mdot |
|---|---|---|---|---|---|---|
| 50 | 1.0099 | 0.245857 | 9.20017e+09 | 3.08309e+10 | 4.00311e+10 | 0.229826 |
| 150 | 1.0467 | 0.243330 | 4.84668e+09 | 1.75615e+10 | 2.24082e+10 | 0.216290 |
| 250 | 1.1639 | 0.208116 | 2.75459e+09 | 1.10575e+10 | 1.38121e+10 | 0.199433 |
| 300 | 1.2978 | 0.199767 | 2.73938e+09 | 1.10154e+10 | 1.37548e+10 | 0.199159 |
| 400 | 1.9623 | 0.199179 | 2.73673e+09 | 1.10057e+10 | 1.37424e+10 | 0.199145 |
| 495 | 3.9023 | 0.199145 | 2.73674e+09 | 1.10058e+10 | 1.37425e+10 | 0.199144 |

Above 1.30 R_p both fluxes are constant to 2e-5 and the
partition is 0.1991 flat. Below that the ratio FALLS, 0.2298 to 0.1991,
which is helium being left behind, the direction settling gives. On
LHS 1140 b the ratio RISES over the same kind of interval, which is the
signature that the state's own mass flux is decaying faster than the
element row's.

## 5. The experiments

All from seed S1c of L5a (`<scratch>/seeds/S1c`), the input file of
`<scratch>/L5a_S1c` unchanged except where the row says, `Load IC? True`,
`Restart intent: stationary`, `Solver: Newton`, `EXHALE_PTC_DTAU0=1.0`.

### T5. Pass 1 of each variation

| run | what changed | hydro | mass | momentum | energy | element relaxation |
|---|---|---|---|---|---|---|
| L5a control (READ) | - | info=2 | 1.08E+00 | 2.64E-02 | 1.84E+00 | 24 steps, then outcome 3 at grow 1.605E-02 |
| `L5b_diag_S1c` | `EXHALE_DIFFUSION_CHECK=1` | info=2 | 1.08E+00 | 2.64E-02 | 1.84E+00 | 24 steps, then outcome 3 at grow 1.605E-02 |
| `L5b_omega_S1c` | `EXHALE_DIFF_OMEGA=0.125` | info=2 | 1.08E+00 | 2.64E-02 | 1.84E+00 | 24 steps, then outcome 3 at grow 1.605E-02 |
| `L5b_hydroonly_S1c` | `EXHALE_OUTER_PASSES=1` | info=2 | 1.08E+00 | 2.64E-02 | 1.84E+00 | not entered: no update is taken on the last pass |
| `L5b_maxit1_S1c` | `EXHALE_JFNK_MAXIT=1` | info=1 | 1.96E+00 | 4.78E-02 | 1.54E+00 | 15 steps, then outcome 3 at grow 4.176E-04 |
| `L5b_frozen_S1c` | `EXHALE_JFNK_MAXIT=1`, `EXHALE_PTC_DTAU0=1.0e-12` | info=1 | 1.63E+00 | 3.92E-02 | 1.63E+00 | **fixed point in 36 steps**, drift 1.35E-01 |
| `L5b_equil_S1c` | `Restart intent: stationary equilibrate` | info=2 | 1.78E+00 | 6.82E-03 | 1.47E+00 | **fixed point in 30 steps**, drift 3.98E-01 |

The diagnostic run reproduces the L5a control row for row, so the
diagnostic does not change the arithmetic. The omega row is identical to it
in every printed digit, which is the measurement behind the code reading
that omega cannot act on a refused pass.

The hydro-only row is the answer to the item's last question. Its refusal is
the hydrodynamics':

```
worst refusing entry: hydrodynamic mass row, measure  3.236E-03
                      against  3.018E-12 at cell 224, r =  1.2269
hydrodynamic mass row      max= 1.075E+00  vol= 2.690E-02  cell=69   ABOVE
hydrodynamic momentum row  max= 2.643E-02  vol= 3.092E-05  cell=500  ABOVE
hydrodynamic energy row    max= 1.842E+00  vol= 3.609E-01  cell=5    ABOVE
```

### T6. The two runs that get past pass 1

`L5b_frozen_S1c` (the hydrodynamic solve made nearly the identity, so the
relaxation sees the loaded wind) and `L5b_equil_S1c`:

| run | pass | hydro | mass | energy | element relaxation | drift |
|---|---|---|---|---|---|---|
| frozen | 1 | info=1 | 1.63E+00 | 1.63E+00 | fixed point, 36 steps | 1.35E-01 |
| frozen | 2 | info=1 | 1.84E+00 | 1.89E+00 | fixed point, 32 steps | 5.55E-01 |
| frozen | 3 | info=1 | 1.53E+00 | 1.96E+00 | 24 steps, outcome 3 at grow 1.605E-02 | refused |
| equilibrate | 1 | info=2 | 1.78E+00 | 1.47E+00 | fixed point, 30 steps | 3.98E-01 |
| equilibrate | 2 | info=2 | 1.27E+00 | 8.64E-01 | fixed point, 31 steps | 1.94E-01 |
| equilibrate | 3 | info=2 | 1.27E+00 | 8.64E-01 | fixed point, 30 steps | 9.83E-02 |
| equilibrate | 4 | info=2 | 1.27E+00 | 8.64E-01 | fixed point, 32 steps | 7.41E-02 |
| equilibrate | 5 | info=2 | 1.27E+00 | 8.64E-01 | no update: this pass ends the iteration | - |

The frozen run is the positive control for the verdict: on the loaded wind,
whose `rho v r^2` is uniform to 1.2e-16 per cell (MEASURED on
`seeds/S1c/Hydro_ioniz_IC.txt`), pass 1's relaxation converges with its
maximum pinned at `X_base = 0.865918` and only the minimum moving. The same
operator, the same planet, the same composition: what changes between an
admissible relaxation and a refused one is the wind. The frozen run's
pass-3 refusal is that state after two composition updates, whose mass flux
is no longer uniform (`dlnMdot/dcell = -4.25e-03` at cell 140, MEASURED).

The equilibrate run shows the element drift contracting, 3.98e-1, 1.94e-1,
9.83e-2, 7.41e-2, which is the HD 209458 b control's behaviour and not the
LHS one. Its hydrodynamic rows stop falling at mass 1.27, energy 0.864, and
the run ends at pass 5 on `outer_no_progress`:

```
outer pass 5: REFUSED -- the joint distance of the state (its largest entry
over its own tolerance) has not fallen in 3 consecutive passes; the
alternation is not approaching a joint fixed point at these step lengths.
   it stands at  1.896E-04 against  1.000E-05
   (elemental transport He/H partition) at cell 217
```

## 6. What would change, and where

No source change was made. These are proposals with their locations and the
reproduction each would be judged on.

1. **The refusal message names neither the cell nor the bound.**
   `binary_element_diffusion.f90` lines 1236-1244 print the step count, the
   outcome integer, the step length and the mass-closure departure;
   `element_step_outcome_text` (line 931) says "helium fraction left [0,1]"
   without saying which end. L5a could not tell from the log whether the
   excursion was at 0 or at 1, and this item had to re-measure it with
   `EXHALE_DIFFUSION_CHECK`. Adding the cell, its radius, which bound and the
   excursion would have made that immediate. Reproduction: any run of
   section 5; the entry text is the message as it stands.

2. **The retry budget is not sized against the tolerance it has to clear.**
   `he_relax_retry_max = 20` (line 585) buys a factor 1e6 of excursion
   against `element_fraction_bound_tol = 1e-12` (line 452). Where the
   excursion falls in proportion to the step, the budget needed is
   `log2(Xover/tol)`, about 33 here; where it does not fall in proportion,
   the fixed point is on the bound and no budget helps. The measurement that
   separates the two cases is already taken every retry (`Xover` against the
   halved `grow`), so the loop could stop at the second refusal whose
   excursion did not halve with the step, and say that the operator's fixed
   point is outside the admissible set rather than that no admissible step
   was found. This changes nothing on HD 209458 b, which refuses no step
   (MEASURED: 30, 30 steps at passes 1 and 2 of `L5b_ctrl_hd209`, no retry).

3. **A comment's boundedness argument does not cover the case measured
   here.** Lines 1049-1053 argue that the conservative form is safe because
   "an element inherits the mass row's own imbalance instead of having it
   subtracted out, and a cell can only lose the fraction of its mass the
   mass row loses". That covers `div(F_rho) > 0`. The case here is the other
   sign: a cell whose discrete mass row GAINS drives the element mass
   fraction UP, the drift's `X(1-X)` shutoff does not act on the advective
   term, and X passes 1. The statement should say that the bound holds only
   where the state's mass row is itself converged, which is the condition the
   HD 209458 b route meets (mass row 1.18e-9) and this one does not (1.08).

4. **Outside the item, reported and not fixed:** line 861 of the same file
   uses the word "pipeline" in a comment ("nothing else in the pipeline can
   see it"), which `~/.claude/CLAUDE.md` forbids as a name for a physical or
   numerical object. "nothing else in the run can see it" carries the same
   meaning.

**What this hands back to L5.** The element relaxation is a consequence of
the wind it is given and not an obstruction of its own: on a wind whose mass
flux is uniform it converges (frozen run, pass 1), on the wind the solve
produces it cannot. With the relaxation removed entirely the state is
refused by the hydrodynamic mass and energy rows anyway. The route's entry
should be `Restart intent: stationary equilibrate`, which is measured to
carry the relaxation through every pass of a five-pass run; the obstruction
that remains is the hydrodynamic one, L5c, and with that entry the route
says so itself instead of naming helium.

## 7. Method and reproduction

`EX` is the repository, `S` the session scratch directory
`.../scratchpad/L`. Binary `$S/L5a_bin/EXHALE.x` throughout, which is newer
than every `src/*.f90` of this tree (checked with `find -newermt`).

```
setup() {                       # $1 tag, $2 seed
  mkdir -p $S/$1/output
  \cp -f $S/L5a_S1c/input.inp $S/$1/input.inp
  \cp -f $S/seeds/$2/*.txt $S/$1/output/
}

# the diagnostic run: the X range of each step on stderr, and the last
# candidate's X profile in output/element_flux_profile.txt
setup L5b_diag_S1c S1c
cd $S/L5b_diag_S1c && OMP_NUM_THREADS=16 EXHALE_PTC_DTAU0=1.0 \
   EXHALE_DIFFUSION_CHECK=1 $S/L5a_bin/EXHALE.x > run.log 2> diag.log

# hydrodynamics alone, composition held at the loaded value
setup L5b_hydroonly_S1c S1c
cd $S/L5b_hydroonly_S1c && OMP_NUM_THREADS=16 EXHALE_PTC_DTAU0=1.0 \
   EXHALE_OUTER_PASSES=1 $S/L5a_bin/EXHALE.x > run.log 2>&1

# the relaxation on (nearly) the loaded wind
setup L5b_frozen_S1c S1c
cd $S/L5b_frozen_S1c && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0e-12 \
   EXHALE_JFNK_MAXIT=1 EXHALE_DIFFUSION_CHECK=1 \
   $S/L5a_bin/EXHALE.x > run.log 2> diag.log

# omega, and the equilibrate entry
setup L5b_omega_S1c S1c
cd $S/L5b_omega_S1c && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 \
   EXHALE_DIFF_OMEGA=0.125 $S/L5a_bin/EXHALE.x > run.log 2>&1
mkdir -p $S/L5b_equil_S1c/output
sed 's/^Restart intent: stationary$/Restart intent: stationary equilibrate/' \
    $S/L5a_S1c/input.inp > $S/L5b_equil_S1c/input.inp
\cp -f $S/seeds/S1c/*.txt $S/L5b_equil_S1c/output/
cd $S/L5b_equil_S1c && OMP_NUM_THREADS=16 EXHALE_PTC_DTAU0=1.0 \
   EXHALE_DIFFUSION_CHECK=1 $S/L5a_bin/EXHALE.x > run.log 2> diag.log

# the HD 209458 b control with the same diagnostic, two passes
mkdir -p $S/L5b_ctrl_hd209/output
\cp -f $S/ctrl_hd209_elem/input.inp $S/ctrl_hd209_elem/metals.inp $S/L5b_ctrl_hd209/
\cp -f $S/ctrl_hd209_elem/output/*_IC.txt $S/L5b_ctrl_hd209/output/
cd $S/L5b_ctrl_hd209 && OMP_NUM_THREADS=16 EXHALE_PTC_DTAU0=1.0 \
   EXHALE_OUTER_PASSES=2 EXHALE_DIFFUSION_CHECK=1 \
   $S/L5a_bin/EXHALE.x > run.log 2> diag.log
```

`EXHALE_DIFFUSION_CHECK=1` writes three lines per accepted or
bound-refused operator call to stderr (`report_step`, line 3929) and
replaces `output/element_flux_profile.txt` on every one, so after a run that
file holds the candidate of the LAST call. A call whose Newton row did not
solve returns before both and leaves no trace, which is why the diagnostic
line count is short of the accepted-plus-refused total.

`$S/l5b_parse.py <diag.log>` replays the acceptance sequence to reconstruct
`grow` (start 1, x1.5 on accept, x0.5 on refuse) and prints tables T1 and
T2; `$S/l5b_flux.py <tag> <Hydro_ioniz.txt> ...` prints the mass-flux decay
per cell of a state. Both are measurements of this item and are kept in the
scratch directory, not in the repository.

Code units, as in L1, L3 and L5a: `R_0 = 1.1273716e9 cm`,
`n_0 = 3.204863e13 cm^-3`, `v_0 = 1.365459e5 cm/s`. `Mdot(wind) = 6.3487e7 g/s`.

Run directories: `$S/L5b_{diag,hydroonly,maxit1,frozen,omega,equil}_S1c`,
`$S/L5b_ctrl_hd209`.
