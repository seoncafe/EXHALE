# L26: continuity of the wind-window base boundary, measured

Item L26 of `docs/PLAN_20260913_lhs_stationary.md`, specified in
`docs/PLAN_20260916_rev3.md` section 6. **Probes only. No line of
`src/modules/states/base_boundary.f90` was changed, and no repair is
proposed here: the choice among the three repair directions of section 6 is
the user's.**

Binary and objects: a private build of the tree of 2026-09-16
(`make OBJDIR=build_L26 EXE=EXHALE_L26.x`), deleted after the measurements.
Everything below is MEASURED with
`src/tests/grid_and_gates/base_boundary_continuity_probe.f90`, which calls
`input_read` and `init` exactly as `EXHALE_main.f90` does and then calls the
production entry point `base_boundary_states` on copies of the loaded state.
The probe re-implements no part of the boundary condition. Single-threaded
(`OMP_NUM_THREADS=1`).

## 1. Verdict

**The base face density is discontinuous at a zero window flux, and the
discontinuity is directional, exactly as the second Codex review derived
(R26).** On the two LHS 1140 b wind states the gap between the two
directional limits is **37.4 per cent of the face density**; on the cold-start
initial condition it is **0.41 per cent**. The gap is not an artifact of the
special point: it equals `|w_i - 1/2| |rho_res - rho_rev|` to the last printed
digit on all three states, so it is the branch weight jumping from the value
the uniform direction gives it (1/2) to the value the zero window gives it
(`w_i`), multiplied by the distance between the two isentropes the boundary
mixes. The two isentropes do differ where the gap is large: `rho_rev/rho_res`
is 0.4555 on the certified state and 0.4548 on the transient, and 0.9918 on
the cold start, which is why the same weight jump is worth 37 per cent on the
first two and 0.4 per cent on the third.

**A second, independent defect is measured: the `max`/`min` spread has a
kink.** When the cell carrying the window maximum changes, the two one-sided
derivatives of `rho_b` with respect to that cell's relative flux converge at
three step sizes to **-242.2 and +0.0128** on the certified state (-232.0 and
+0.0130 on the transient): the value is continuous and the slope jumps by a
factor of 1.9e4. This has nothing to do with the zero window and survives any
repair that keeps `(max - min)/|mean|`.

**At a fully static state the value is continuous, as R27 says, and the
neighbourhood is differentiable in the direction tested.** With `v = 0` in
every cell, `w_i = 1/2` and `rho_b` agrees with the zero-window limit of the
uniform path to 6.5e-15 relative; the derivative along a uniform window flux
is the same from both sides and converged at three step sizes.

**The four identities of R48 pass.** Ten repeated evaluations, an evaluation
of another state in between, a complete reconstruction and Riemann evaluation
in between, and a changed call order all return the same bits. No stale flux
is a hidden argument of the boundary today.

## 2. What was probed, and on which states

Three states, all read-only copies taken into a scratch directory; the probe
writes nothing into the catalog.

| state | source | what it is |
|---|---|---|
| certified 0.03 XUV | `LHS1140b/models/.L14/x003_HeH2.13/output/` (`Hydro_ioniz.txt`, `Ion_species.txt`, its `input.inp`) | header `certified=T cert_reason=certified_in_wind` |
| stalled 0.02 transient | `LHS1140b/models/.L14/x002_HeH2.13/output/` | header `certified=F cert_reason=no_stationary_claim` |
| cold-start initial condition of the 0.02 case | the same `input.inp` with `Load IC? False` (and the `Restart intent:` line removed, which `input_read` refuses without a loaded state), so `set_IC` builds it | the code's own isothermal hydrostatic start |

Note on the state named in the plan: the plan asks for the `x0.02` transient
of `x002_HeH9.7`. Neither `.L14/x002_HeH9.7` nor `.L4h/x002_HeH9.7` holds a
written state, only the `_IC` seed pair it was started from; `.L14/x002_HeH2.13`
is the only `x0.02` case in the catalog with a complete written
`Hydro_ioniz.txt` and `Ion_species.txt`, and it is the one used. Both restart
pairs were placed as `output/*_IC.txt` in the scratch copy, which is how the
code reads a restart (`molecular_seed_state_file`).

All three carry the same grid and the same window: `j_flux = 217`, 284 window
cells, `r(j_flux) = 1.2006911`. The boundary options are the ones each case
runs with (`Base BC: pressure 1.0`, `Well balanced: True`,
`base_face_mach_blend = 1e-8`, `base_wind_window_spread = 1e-2`).

The paths are those of section 6. In each, only the VELOCITIES are changed;
density and pressure are untouched, so the composition, the reservoir
isentrope and the interior isentrope are those of the loaded state throughout
and only the discriminant of the branch moves.

## 3. The state as it stands

| quantity | certified 0.03 | transient 0.02 | cold start |
|---|---|---|---|
| window flux `F0` [code units] | 1.673449989E-08 | 1.106486318E-08 | 8.448286912E-06 |
| `du_window` | 3.580686E-03 | 3.771033E-03 | 1.398195E+01 |
| `M_i` | -8.492845E-07 | -8.498428E-07 | +1.270140E-04 |
| `M_wind` | +1.132468E-08 | +7.493160E-09 | +3.748064E-06 |
| `w_i` | 1.0000000 | 1.0000000 | 0.0000000 |
| `s_wind` | 1.0000000 | 1.0000000 | 0.0000000 |
| `w_rev` | 0.0000000 | 0.0431934 | 0.0000000 |
| `rho_res` | 2.966075941E+00 | 2.966075941E+00 | 2.966075941E+00 |
| `rho_rev` | 1.351126958E+00 | 1.349226309E+00 | 2.941819408E+00 |
| `rho_b` | 2.966075941E+00 | 2.896238651E+00 | 2.966075941E+00 |

## 4. Path (i): a uniform window flux going to zero, and the two limits

The window velocities are set so that `rho v r^2 = eps F0` in every window
cell, for `eps = +/-1e-2, 1e-4, ... 1e-14` and then exactly zero. The
resulting `du_window` is at the rounding floor (5e-16 to 7e-16), so
`s_wind = 1` and the window decides, until `eps = 0` makes the window mean
exactly zero, `have_F` false and the first interior cell decide.

Certified 0.03 state, `rho_b` against `eps` (both signs give the same limit):

| `eps` | `rho_b` (eps > 0) | `w_rev` | `rho_b` (eps < 0) | `w_rev` |
|---|---|---|---|---|
| 1e-02 | 2.172317451106499 | 0.4915069 | 2.144885448187360 | 0.5084931 |
| 1e-04 | 2.158738615524708 | 0.4999151 | 2.158464283769151 | 0.5000849 |
| 1e-06 | 2.158602821305713 | 0.4999992 | 2.158600077988146 | 0.5000008 |
| 1e-08 | 2.158601463363517 | 0.5000000 | 2.158601435930341 | 0.5000000 |
| 1e-10 | 2.158601449784095 | 0.5000000 | 2.158601449509763 | 0.5000000 |
| 1e-12 | 2.158601449648301 | 0.5000000 | 2.158601449645558 | 0.5000000 |
| 1e-14 | 2.158601449646943 | 0.5000000 | 2.158601449646915 | 0.5000000 |
| **0** | **1.351126958008022** | **1.0000000** | | |

The gap, on all three states:

| state | `rho_b` (eps -> 0) | `rho_b` (eps = 0) | gap | per cent of `rho_b(eps -> 0)` | `\|w_i - 1/2\| \|rho_res - rho_rev\|` |
|---|---|---|---|---|---|
| certified 0.03 | 2.158601449646943 | 1.351126958008022 | +8.074744916E-01 | 37.407 | 8.074744916E-01 |
| transient 0.02 | 2.157651125119984 | 1.349226308954115 | +8.084248162E-01 | 37.468 | 8.084248162E-01 |
| cold start | 2.953947674438118 | 2.966075941285836 | -1.212826685E-02 | 0.411 | 1.212826685E-02 |

The bound and the gap agree exactly, on every state.

The one-sided derivative of `rho_b` at `eps = 0` along this path, by forward
differences at `h = 1e-8, 1e-10, 1e-12`:

| state | direction | `D(1e-8)` | `D(1e-10)` | `D(1e-12)` | ratio per decade pair |
|---|---|---|---|---|---|
| certified 0.03 | + | 8.074745E+07 | 8.074745E+09 | 8.074745E+11 | 100.000, 100.000 |
| certified 0.03 | - | -8.074745E+07 | -8.074745E+09 | -8.074745E+11 | 100.000, 100.000 |
| transient 0.02 | + | 8.084248E+07 | 8.084248E+09 | 8.084248E+11 | 100.000, 100.000 |
| cold start | + | -1.212820E+06 | -1.212827E+08 | -1.212827E+10 | 100.001, 100.000 |

Each estimate is one hundred times the last when the step is a hundred times
smaller, which is `gap/h`: the difference quotient diverges, which is what a
jump does and what a kink or a smooth point does not.

## 5. Path (ii): nonuniform shapes, and the direction dependence

Four shapes, each scaled by the same `eps`: `1 + 0.5 sin` (spread 1.000),
alternating 1 and 2 (spread 0.6667), a linear ramp 0.5 to 1.5 (spread 1.000),
and `1 + 0.002 sin` (spread 4.000E-03, below the 1e-2 threshold).

On the three shapes whose spread exceeds twice the threshold, `s_wind = 0` at
every `eps` down to 1e-14 AND at `eps = 0`: `rho_b` is the SAME number at
every point of the path, the derivative is exactly zero at all three step
sizes, and the limit is `w_i`'s value.

| state | limit along the spread-above-threshold shapes | limit along the uniform shape | difference, per cent |
|---|---|---|---|
| certified 0.03 | 1.351126958008022 | 2.158601449646943 | 37.407 |
| transient 0.02 | 1.349226308954115 | 2.157651125119984 | 37.468 |
| cold start | 2.966075941285836 | 2.953947674438118 | 0.411 |

**Two paths reach the same zero window with different limits.** The fourth
shape, whose spread is 4.0e-3 and therefore below the threshold, behaves
exactly like the uniform path: same values at every `eps`, same jump at zero,
same 1/h derivative. So the split is not between "uniform" and "nonuniform"
but between "spread below the threshold" and "spread above twice it", which
is the calibrated gate itself.

## 6. Path (iii): zero mean with nonzero variance

Alternating `+eps F0` and `-eps F0` over the window, the last cell zeroed
where the cell count is odd, so the mean is zero by construction. The mean is
not exactly zero in floating point at every `eps` (at `eps = 1e-4` on the
certified state it rounds to exactly zero and `have_F` is false; at the other
`eps` it is a rounding residue and `have_F` is true), but `du_window` is then
1.8e17 to 4.6e17, far above twice the threshold, so `s_wind = 0` either way.

`rho_b` is the `w_i` value at every point of this path and at zero, on all
three states, with a zero derivative. **The zero-mean direction is continuous
and agrees with the nonuniform limit, not with the uniform one.** Whether
`have_F` is true or false makes no difference to `rho_b` here, because both
routes end at `w_i`; the discontinuity of section 4 is not the `have_F` test
by itself but the pair (`have_F` false, or spread above the gate) against
(spread below the gate).

## 7. Path (iv): the first cell and the window disagreeing in sign

The window is held at a uniform `+F0` (gas enters) or `-F0` (gas leaves) and
`v(1)` is swept through zero from `-|v1|` to `+|v1|`, so the two discriminants
disagree over half of each sweep.

Certified 0.03, window at `+F0`:

| `t` (multiplier on `\|v(1)\|`) | `w_i` | `s_wind` | `w_rev` | `rho_b` |
|---|---|---|---|---|
| -1.0 | 1.0000000 | 1.0000000 | 0.0000000 | 2.966075941285836 |
| -0.5 | 1.0000000 | 1.0000000 | 0.0000000 | 2.966075941285836 |
|  0.0 | 0.5000000 | 1.0000000 | 0.0000000 | 2.966075941285836 |
| +0.5 | 0.0000000 | 1.0000000 | 0.0000000 | 2.966075941285836 |
| +1.0 | 0.0000000 | 1.0000000 | 0.0000000 | 2.966075941285836 |

`rho_b` does not move at all: with `s_wind = 1` the local weight has no weight,
and the one-sided derivatives are exactly zero at all three step sizes in both
directions. With the window at `-F0` the same sweep holds `rho_b` at
1.351126958008022 (the certified state) and 1.419063599238400 (the transient,
whose `|M_wind|` is inside the handover so `w_rev` is 0.9569 and not 1). On
the cold start `s_wind = 0` and the derivative is also exactly zero, because
`M_i` is 1.27e4 blend widths and `w_i` is saturated.

## 8. Path (v): the extremal cell of the window changes

Two window cells (`j = 227` and `j = 237`) are bumped above an otherwise
uniform window, the first by a fixed `alpha = 1.5e-2` and the second by `t`.
The maximum is carried by the first while `t < alpha` and by the second while
`t > alpha`, and `alpha` was chosen so that `du_window` sits inside the
handover band (threshold to twice it), where `s_wind` is strictly between 0
and 1.

Certified 0.03 state:

| `t` | `du_window` | `s_wind` | `w_rev` | `rho_b` |
|---|---|---|---|---|
| 1.000E-02 | 1.499868E-02 | 0.5001980 | 0.4998020 | 2.158921283739681 |
| 1.200E-02 | 1.499857E-02 | 0.5002139 | 0.4997861 | 2.158946868033456 |
| 1.400E-02 | 1.499847E-02 | 0.5002297 | 0.4997703 | 2.158972451966618 |
| 1.500E-02 | 1.499842E-02 | 0.5002377 | 0.4997623 | 2.158985243797964 |
| 1.600E-02 | 1.599825E-02 | 0.3522515 | 0.6477485 | 1.919995133621702 |
| 1.700E-02 | 1.699808E-02 | 0.2162414 | 0.7837586 | 1.700345736898498 |
| 1.800E-02 | 1.799791E-02 | 0.1042008 | 0.8957992 | 1.519406005259923 |
| 1.900E-02 | 1.899773E-02 | 0.0281229 | 0.9718771 | 1.396544072023231 |
| 2.000E-02 | 1.999754E-02 | 0.0000002 | 0.9999998 | 1.351127252220484 |

One-sided derivatives at `t = alpha`:

| state | direction | `D(1e-6)` | `D(1e-7)` | `D(1e-8)` |
|---|---|---|---|---|
| certified 0.03 | + | -2.422040E+02 | -2.422039E+02 | -2.422039E+02 |
| certified 0.03 | - | +1.279179E-02 | +1.279179E-02 | +1.279177E-02 |
| transient 0.02 | + | -2.320207E+02 | -2.320207E+02 | -2.320207E+02 |
| transient 0.02 | - | +1.295564E-02 | +1.295532E-02 | +1.295222E-02 |
| cold start | + | 0.000000E+00 | 0.000000E+00 | 0.000000E+00 |
| cold start | - | -0.000000E+00 | -0.000000E+00 | -0.000000E+00 |

Both sides converge, and to different numbers: the value is continuous, the
slope is not. On the two wind states the slope jumps by 242.2 and 232.0 in
`rho_b` per unit relative bump of ONE window cell, i.e. by 1.9e4 times the
slope on the other side. Over the sampled interval `rho_b` falls from 2.1590
to 1.3511 (37 per cent) while one cell's flux changes by 5e-3 of the window
mean. On the cold start both sides are zero, because `w_i = 0` there and the
window weight is also 0, so the two branches `s_wind` mixes are the same
number and moving `s_wind` moves nothing.

## 9. Path (vi): a base flux the window does not carry

Cells 1 to 6 scaled by 0.5, 1, 2, 0 and -1, the window left exactly as the
state carries it.

| state | `t = 0.5` | `t = 1` | `t = 2` | `t = 0` | `t = -1` |
|---|---|---|---|---|---|
| certified 0.03, `rho_b` | 2.966075941285836 | 2.966075941285836 | 2.966075941285836 | 2.966075941285836 | 2.966075941285836 |
| certified 0.03, `w_i` | 1.0000000 | 1.0000000 | 1.0000000 | 0.5000000 | 0.0000000 |

`rho_b` is one number across the whole range, including the sign reversal of
the base velocity, and the one-sided derivatives at `t = 1` are exactly zero.
**On these states the boundary's entropy source is decided entirely by cells
at `r >= 1.2007` and not at all by the six cells above it.** That is the
stated design (the source's own note at `base_branch_on_wind_flux`), measured;
it is also the stated LIMITATION of that design, that a transient whose base
carries a flux the wind does not is bounded by the wind's direction.

## 10. Path (vii): the static states

With `v = 0` in every cell, the window mean is exactly zero, `have_F` is
false, and `w_i = 1/2` because `M_i = 0`:

| state | `rho_b` at `v = 0` | `rho_b` (eps -> 0 of path (i)) | relative difference |
|---|---|---|---|
| certified 0.03 | 2.158601449646929 | 2.158601449646943 | 6.5E-15 |
| transient 0.02 | 2.157651125119975 | 2.157651125119984 | 4.2E-15 |
| cold start | 2.953947674438050 | 2.953947674438118 | 2.3E-14 |

The two limits coincide at a fully static base, as R27 states. The
neighbourhood was also tested: from `v = 0` everywhere, a uniform window flux
`t F0` was switched on and the one-sided derivatives taken at `t = 0`.

| state | `D(1e-8)` | `D(1e-10)` | `D(1e-12)` | direction |
|---|---|---|---|---|
| certified 0.03 | 1.371659E+00 | 1.371658E+00 | 1.371347E+00 | + |
| certified 0.03 | 1.371659E+00 | 1.371663E+00 | 1.371792E+00 | - |
| transient 0.02 | 9.086485E-01 | 9.086509E-01 | 9.086065E-01 | + |
| transient 0.02 | 9.086485E-01 | 9.086465E-01 | 9.086065E-01 | - |
| cold start | 6.818629E+00 | 6.818626E+00 | 6.818546E+00 | + |
| cold start | 6.818629E+00 | 6.818635E+00 | 6.818546E+00 | - |

Converged and equal from both sides to five digits, the residual movement
being the cancellation of two nearly equal doubles at the smallest step. In
this direction the static state is differentiable.

The closed-gate direction out of the same static state was measured too: from
`v = 0` everywhere, a `1 + 0.5 sin` window of amplitude `t F0` (spread 1.000,
above twice the gate, so `s_wind = 0`). `rho_b` is then bitwise the static
value at every `t` from 1e-2 down to 1e-14 and at `t = 0`, on all three
states, with a derivative of exactly zero at all three step sizes. So both
directions out of the static base agree with each other and with the static
value: no jump appears at a static base in either direction tested, which is
R27 measured rather than argued.

With `v = 0` only in cells 1 to 6 and the window as the state carries it,
`rho_b` is the state's own value (2.966075941285836 on the certified state,
2.896238651001553 on the transient), which is section 9 again.

## 11. R48: the boundary is a function of its argument

Four assertions, all PASS on all three states, printed by the probe as
`PASS|FAIL` rows.

| assertion | certified 0.03 | transient 0.02 | cold start |
|---|---|---|---|
| ten repeated evaluations of one state give one `rho_b` | PASS, 0 of 10 differ | PASS | PASS |
| an evaluation of ANOTHER state in between changes nothing | PASS, bitwise | PASS | PASS |
| a complete `W_to_U`, `Reconstruct`, `RK_rhs` evaluation in between changes nothing | PASS, bitwise | PASS | PASS |
| boundary, flux, boundary, flux, boundary gives the same `rho_b` as boundary alone | PASS, bitwise | PASS | PASS |

These are invariants and not diagnostics, so they are the probe's only
verdict rows.

## 12. What the probes say about the three repair directions

Neutral, and in the order section 6 lists them. **No direction is chosen
here.**

**A local face-flux reading (the P44 rule: the base FACE flux, not the cell
velocity).** Section 9 measures that the six base cells have no influence on
`rho_b` at all while the window gate is open, and that their velocity may
reverse without moving the face density. So a local reading would change the
answer on exactly these states, and by the full 37 per cent if it returned the
first cell's own sign: `w_i` is 1 on both wind states while the window says
inflow. Whether the base FACE flux agrees with the window is not measured
here; it is READ from item L21, which reports the Riemann face mass flux at
+1.00000 `F_wind` at every face of the certified state including the base face,
against a cell-centred product of -2.00 `F_wind` at cell 1. The constraint
R48 puts on this direction, that the local construction be a function of the
current state and the reservoir alone and not of a stored flux, is now
testable: section 11 is that test, and it passes today, so a change that
introduced a stored flux would be caught by it.

**A nonlocal diagnostic kept with an amplitude weight.** The numbers such a
weight has to remove are in sections 4 and 5: a gap of 37.4 per cent of
`rho_b` at the zero window on both wind states and 0.41 per cent on the cold
start, reached along the uniform direction and along any shape whose spread is
below the 1e-2 gate (the 4.0e-3 shape of section 5 behaves exactly like the
uniform one), and absent along every shape whose spread exceeds twice the
gate. Section 8 adds a requirement the amplitude weight does not address: the
`max`/`min` functional has a slope jump of 242 at the argmax switch, in the
band where the gate is half open, and that jump is present on both wind states
and independent of the window's amplitude. R28 already refuted the RMS
substitute on three separate grounds.

**A base closure from the reservoir and the local characteristics.** Section
10 measures that at a fully static base the present boundary is already
continuous and differentiable, so a local closure would agree with it there;
sections 4, 5 and 8 are all properties of the window diagnostic, so all three
disappear with it. What the probes do NOT say is what a local closure gives on
the two wind states: with `w_i = 1` on both, a closure that read the first
cell's own sign would put `rho_b` at `rho_rev`, a factor 0.4555 below the
present value, and item L21 measured that the resulting state has a mass
certification row of 1.31 instead of 2.17e-09. A local closure that is a face
construction and not the cell-centred product is a different object from
`w_i`, and nothing here measures it.

## 13. Noticed outside the item, reported and not fixed

**The 0.02 transient's base operating point sits INSIDE the reversal
handover.** `M_wind = +7.493160E-09` against `base_face_mach_blend = 1e-8`,
i.e. 0.749 of a blend width, so `w_wind = 0.0431934` and the face carries
4.32 per cent of the interior isentrope although the window reads a clean
inflow: `rho_b = 2.896238651` against `rho_res = 2.966075941`, 2.35 per cent
low. The certified 0.03 state sits at `M_wind = +1.132468E-08`, 1.13 widths,
so it clears the handover by 13 per cent and `w_rev` is exactly zero. This is
the same condition item L21 measured at a blend width of 1e-6 and repaired by
moving the width to 1e-8, recurring at 1e-8 one XUV rung lower; the memo's own
requirement is that the width stand "far BELOW the physical operating point".
It is a width, not an L26 question, and nothing was changed.

**`docs/PLAN_20260916_rev3.md` section 6 names `x002_HeH9.7` as the transient
state.** That directory holds no written state. The plan row should name
`x002_HeH2.13`, which is the only `x0.02` case in the catalog with a complete
written pair; the row was not edited, since the item's own row in
`docs/PLAN_20260913_lhs_stationary.md` is the one this item updates.

## 14. Reproduce

```bash
make OBJDIR=build_L26 EXE=EXHALE_L26.x
# a scratch copy of a run directory, with the written state placed as the
# restart pair the code reads:
#   <case>/input.inp
#   <case>/output/Hydro_ioniz_IC.txt   <- the run's Hydro_ioniz.txt
#   <case>/output/Ion_species_IC.txt   <- the run's Ion_species.txt
EXHALE_OBJDIR=$PWD/build_L26 EXHALE_EXE=$PWD/EXHALE_L26.x \
EXHALE_L26_STATE=<the scratch copy> \
   src/tests/grid_and_gates/run.sh base_continuity
```

The probe writes `EXHALE_setup.out` and an `./output` directory into the
directory it runs in, which is why it is given a scratch copy and never a
catalog case. Without `EXHALE_L26_STATE` the suite prints one DIAGNOSTIC line
and skips it.

---

# Repair (2026-09-17)

Item L26 of `docs/PLAN_20260913_lhs_stationary.md`, the repair the probes above
were taken before, approved by the user on 2026-09-17. Everything in this
section is MEASURED with a private build of the tree of 2026-09-17
(`make OBJDIR=build_L26r EXE=EXHALE_L26r.x`) against a control built out of
tree from the same snapshot with `src/modules/states/base_boundary.f90` and the
two test programs at their entry text. Single-threaded
(`OMP_NUM_THREADS=1`) for every frozen-state measurement.

## R1. Verdict

**The closure adopted is the nonlocal one the plan's section 6 allows as its
second direction: the wind window is kept, and its say is made an AMPLITUDE
times a SHAPE,**

```
s_wind = A(M_wind) C(d_window) ,
d_window = sqrt(<F^2> - <F>^2)/|<F>| over r >= r_flux ,   F = rho v r^2 .
```

`A` is the cubic smoothstep in `|M_wind|`, zero at a zero window flux and one
at `|M_wind| >= 1e-12`; `C` is the cubic smoothstep in `d_window`, one at or
below `2e-3` and zero at or above twice it. The two measured defects are
removed by one term each: `A` vanishes quadratically with the window flux, so
`s_wind` and its derivative go to zero along EVERY direction and the zero
window is an ordinary point of the closure; `d_window` is a second moment and
has no argmax, so the slope jump at the extremal-cell switch has nowhere to
come from.

**The three local directions the brief named were measured first and all
three fail, on a property that is not a matter of tolerance.** On the four
LHS 1140 b wind states every reading of the base mass flux that is a function
of the base cells alone is NEGATIVE, while the branch those states carry is
inflow; on the two states whose base face the Riemann solve reads as an inflow
the local readings disagree with it in sign AND by one to three orders of
magnitude. Section R3 has the table. The base velocity of these winds is 1e-6
of the sound speed and the base numerical artifact is twenty to three hundred
times the wind's own flux, so the direction of the flow at the base face is
not resolved by the discretization there at all: there is no local reading to
choose, and this is a statement about the states and not about the candidates.

**Nothing the code meets moves.** On all five probed states the branch weight
`w_rev`, the face density `rho_b` and both ghost densities are bit for bit
what the entry text returns; the two certified LHS 1140 b states re-evaluated
under `Restart intent: stationary evaluate` write outputs that differ from the
control's in the provenance line alone; the fiducial re-solved from its own
certified state under the new boundary returns `info = 0`, CERTIFIED, and
reproduces the catalog profile to 2e-12 in velocity and 2e-16 in density.

## R2. What the defect was, in one line each

- **The zero-window jump.** `C` alone is scale free, so it has no limit at the
  zero window: along a uniform `eps F` it is 1 at every `eps` and along a
  shape whose spread is above the gate it is 0 at every `eps`. MEASURED gap
  between the two limits, entry text: 37.407 per cent of `rho_b` on the
  certified 0.03 state, 37.468 on the 0.02 transient, 4.188 on the certified
  0.10 state, 3.794 on the fiducial, 0.411 on the cold start.
- **The argmax kink.** `(max - min)/|mean|` keeps its value and changes its
  slope where the cell carrying an extremum changes. MEASURED one-sided
  derivatives of `rho_b`, entry text, converged at three step sizes:

  | state | from above | from below | ratio |
  |---|---|---|---|
  | certified 0.03 | -2.422039E+02 | +1.279177E-02 | 1.89E+04 |
  | transient 0.02 | -2.320207E+02 | +1.295222E-02 | 1.79E+04 |
  | certified 0.10 | -3.576421E+01 | +1.888800E-03 | 1.89E+04 |
  | fiducial | -3.251745E+01 | +1.717382E-03 | 1.89E+04 |
  | cold start | 0 | 0 | the two branches coincide there |

## R3. The local readings, measured, and why each is refused

`src/tests/grid_and_gates/base_boundary_continuity_probe.f90` gained the group
`local readings of the base mass flux`, which prints every candidate beside the
mass flux the Riemann solve puts through the base face on the same state. All
in units of the window's own mean flux `F0`.

| reading | fiducial | 0.10 certified | 0.03 certified | 0.02 transient | cold start |
|---|---|---|---|---|---|
| Riemann flux at the base face | **+32.38** | **+336.55** | **-404.32** | **-584.45** | -4735.68 |
| cell 1 alone | -18.74 | -229.56 | -74.99 | -113.42 | +33.89 |
| mean over cells 1-2 | -14.65 | -176.52 | -80.27 | -121.42 | +50.16 |
| mean over cells 1-4 | -6.83 | -87.98 | -38.87 | -59.05 | +81.39 |
| mean over cells 1-6 | -4.20 | -58.14 | -25.38 | -38.73 | +110.96 |
| least squares, cells 1-2 to the face | -22.82 | -282.59 | -69.72 | -105.41 | +17.62 |
| least squares, cells 1-4 to the face | -20.65 | -247.60 | -99.34 | -150.13 | +18.92 |
| least squares, cells 1-6 to the face | -15.51 | -187.69 | -78.79 | -119.18 | +21.00 |
| the face state's own `rho_b v_b r_b^2` | -249.90 | -2780.77 | -3254.55 | -4797.28 | -6368.82 |
| the branch these states carry | inflow | inflow | inflow | inflow | inflow |

**Candidate (1), the mass flux through the base face reconstructed from cells
1..k without the Riemann solve, fails "a converged state must not change
branch".** Every mean and every extrapolation is negative on all four wind
states, so each of them puts the boundary on the reversal branch, which is the
`w_i = 1` state item L21 measured at a mass certification row of 1.31 against
2.17e-09. The two-cell mean does not annihilate the mode as the brief's
reasoning supposed: the base artifact of these states is not a Nyquist
component but a two-cell boundary layer, cells 1 and 2 carrying -18.7 and
-10.6 `F0` on the fiducial and -75.0 and -85.5 on the 0.03 state, both of one
sign, so averaging them deepens it instead of cancelling it. The bias test the
brief asks for cannot even be taken: the reference it names, the face flux,
is +32.4 and +336.6 `F0` on two of the states and -404.3 and -584.5 on the
other two, so the local readings are wrong in sign against the reference on
two states, wrong in magnitude by one to three orders on all four, and the
reference itself does not agree with the wind on two of them.

**Candidate (2), the same with a smooth Mach weight on that reconstructed
flux, fails for the same reason and adds nothing.** A weight cannot change the
sign of what it weighs; it can only hand the branch back to the same negative
reading with a smooth transition.

**Candidate (3), the window kept only as a printed diagnostic, is candidate
(1) with `k = 1`.** With the window out of the closure the branch is `w_i`,
and `w_i = 1` on all four wind states (`M_i` is -8.5e-07 to -5.6e-06 against a
blend of 1e-8, tens to hundreds of widths), which is the pre-L21 boundary and
its 418 K base against `base.inp`'s 226 K.

**Why no local reading can succeed on these states, stated as physics.** In a
steady escaping atmosphere the base face carries the wind's mass flux, and at
the base that flux corresponds to a face Mach number of order 1e-6. The
numerical artifact of the base layer carries twenty to three hundred times
that flux with the opposite sign. A discriminant is being asked for the sign
of a quantity that the discretization does not resolve at that location; the
window is the nearest place where the same conserved flux IS resolved, which
is the design note at `base_branch_on_wind_flux` and remains the reason the
window is read. The plan's preferred direction is therefore not available at
this resolution, and the plan's second direction, with the regularization it
names, is the one taken.

One further measurement, reported and not acted on: **the base face flux the
Riemann solve produces is not the wind's flux on any of these four states**,
by factors of +32, +337, -404 and -584, and the faces above it read -18.6,
-203.7, -153.9 and -232.5 `F0` before settling near +1 at face 4 to 8. These
states are certified IN THE WIND, `r >= 1.2`, and the mass row below that is
gated only against the rounding anchor of cells whose density is 1e14, so a
flux imbalance of hundreds of times the wind's own flux passes it. That is a
statement about what "certified" covers at the base, not about this item's
closure, and it belongs with L25's base odd-even mode and P44.

## R4. The threshold, recalibrated on the new functional

The range and the standard deviation are different functionals and the
calibrated threshold does not carry over, which is what review 1 section 3.3
says. MEASURED 2026-09-17 on `d_window`, over the same window
(`r >= r_flux = 1.2`), on the physical cells of every written state that
carries a claim:

| population | `(max-min)/|mean|` | `sqrt(<F^2>-<F>^2)/|<F>|` |
|---|---|---|
| 302 certified states of the two LHS 1140 b catalogs | 1.51e-04 to 3.58e-03, median 3.83e-04 | 4.06e-05 to 1.82e-04, median 6.09e-05 |
| the stalled 0.02 transient | 3.77e-03 | 1.90e-04 |
| carrier_elem_newton, its solved state | 5.04e-02 | 1.39e-02 |
| carrier_model_a_newton, its IC | 7.50e-02 | 1.51e-02 |
| wasp_he23off | 1.45e-01 | 4.43e-02 |
| wasp_full | 1.46e-01 | 4.45e-02 |
| hp_front | 8.13e-01 | 1.63e-01 |
| mol_base_handoff | 9.13e-01 | 1.83e-01 |
| carrier_elem_newton, its IC | 1.02e+00 | 2.24e-01 |
| lower_profile | 2.32e+00 | 7.87e-01 |
| atomic_elem_newton | 3.51e+00 | 8.22e-01 |
| hydrostatic_column | 1.26e+01 | 2.40e+00 |

The gap on the new functional runs from 1.90e-04 to 1.39e-02, a factor 73, and
`sqrt(1.90e-04 * 1.39e-02) = 1.6e-03`, which is **2e-3 to one digit**: ten
times above the loosest state that is a wind, and the tightest state that is
not a wind sits 3.5 times above TWICE it, where the window falls silent. Both
populations are therefore in the saturated parts of the weight and no state
sits in the handover, which is the same margin the range threshold had.

## R5. The amplitude scale, and the sensitivity to it

`base_wind_window_amplitude = 1e-12` in face Mach number. It is a
regularization and not a gate: `A` is exactly 1 for `|M_wind| >= 1e-12`, and
the MEASURED window Mach numbers are 2.99e-07 (fiducial), 2.68e-08 (0.10
certified), 1.13e-08 (0.03 certified), 7.49e-09 (0.02 transient) and 3.75e-06
(the cold start), with 3.7e-09 READ from item L21 as the smallest any state has
been measured at, at the cold start of `wasp_full`. 1e-12 stands 3700 times
below that.

The sensitivity is therefore a step function and was measured as one: the
probe's `state` row is bit for bit the entry text's on all five states, so the
answer on every state is independent of this scale over the whole range in
which the scale stays below the states' own `M_wind`. `A` varies only over
`|M_wind| < 1e-12`, which on the 0.03 state is a coherent change of 1e-4 of
the window's flux; a Jacobian difference quotient perturbs the state by about
1e-8 relative, eight decades short of reaching it, so the curvature `A`
carries is not differenced by the solver either.
`EXHALE_BASE_WIND_AMPLITUDE=<value>` overrides it.

## R6. The probe tables after the repair

The probe is now the suite: 32 `PASS|FAIL` rows, all PASS on all five states
(the three of section 2, plus the certified 0.10 control
`atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13` and the fiducial
`atomic_scalar_gj1132_kzz1e9/HeH2.13`). The rows are

- `base_continuity_*_limit`, seven of them: the face density at the smallest
  sampled amplitude of a path against its value at amplitude zero, relative to
  itself, tolerance 1e-12. **MEASURED 0.000000E+00 exactly, on every path and
  every state.** The entry text reads 8.07e-01 on the uniform path of the 0.03
  state, which is 37.4 per cent of `rho_b` and 3.7e11 times this tolerance.
- `base_continuity_*_derivative`, sixteen of them: the difference quotient must
  not GROW as the step shrinks, which is the one thing a jump does. Bound
  1.001. The entry text's uniform path reads 100.000 per step of the sweep,
  `gap/h` exactly.
- `base_continuity_extremum_no_slope_jump`: the two one-sided derivatives at
  the extremal-cell switch, relative, tolerance 1e-3. MEASURED 4.2e-08
  (certified 0.03), 1.6e-07 (0.02 transient), 1.2e-08 (certified 0.10),
  9.2e-08 (fiducial), 0 (cold start), against the 1.9e4 of the entry text.
  The bump of that path is now set at run time to
  `1.5 * threshold * sqrt(nw/2)` so that the sweep sits INSIDE the handover
  band whatever the window's size is; on these states it is 3.57e-02 and
  `s_wind` runs from 0.657 to 0.349 across the sweep, so the switch happens
  where it can move `rho_b` and the row is not passing vacuously.
- the four R48 identities, unchanged and PASS.

The state rows, entry text against repair, on all five states:

| state | `rho_b`, entry text | `rho_b`, repair | `w_rev` both | `s_wind` both | `d_window` (was `(max-min)/|mean|`) |
|---|---|---|---|---|---|
| fiducial | 2.966075941285836 | 2.966075941285836 | 0 | 1 | 5.585392E-05 (2.457839E-04) |
| certified 0.10 | 2.966075941285836 | 2.966075941285836 | 0 | 1 | 1.520590E-04 (2.939874E-03) |
| certified 0.03 | 2.966075941285836 | 2.966075941285836 | 0 | 1 | 1.817283E-04 (3.580686E-03) |
| transient 0.02 | 2.896238651001553 | 2.896238651001553 | 0.0431934 | 1 | 1.904913E-04 (3.771033E-03) |
| cold start | 2.966075941285836 | 2.966075941285836 | 0 | 0 | 2.474035E+00 (1.398195E+01) |

Both ghost densities agree bitwise as well. Only the printed spread changes,
which is the functional and not the answer.

## R7. Impact

**The two certified LHS 1140 b states, re-evaluated.** `Restart intent:
stationary evaluate` on `atomic_scalar_gj1132_kzz1e9/HeH2.13` and
`atomic_scalar_gj1132x0.10_kzz1e9/HeH2.13`, from their own certified pair,
control against repair. Both are CERTIFIED under both binaries with the same
row values (`hydrodynamic energy row max = 1.700E-08 at cell 45` on the
fiducial, `3.886E-07 at cell 1` on the 0.10 case) and the same gated element
row (4.051E-01 and 2.262E-06 at cell 2). `Hydro_ioniz.txt`,
`Ion_species.txt` and `Hydro_ioniz_adv.txt` differ in the `# provenance` line
alone, which names the control's out-of-tree build. **The base rows do not
move at all.**

**The fiducial re-solved under the new boundary**, its own recipe
(`OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 EXHALE_OUTER_PASSES=40`), from its own
certified state, in `LHS1140b/models/.L26/fid_resolve/`: `info = 0`,
CERTIFIED, rows mass 1.734E-09 at cell 1 against a tolerance of 3.8E-11
volume-weighted 4.298E-11, momentum 3.495E-12 at cell 484, energy 1.700E-08 at
cell 45. Against the catalog state, over all 504 rows:

| quantity | largest relative difference |
|---|---|
| rho | 1.83e-16 |
| v | 2.05e-12 |
| p | 6.01e-14 |
| T | 6.01e-14 |

The base rows `r = 1.000193` to `1.000967` are the catalog's to every printed
digit: T 238.003, 257.716, 278.100, 297.682, 316.272 K, v -0.597, -0.386,
+0.0752, +0.00845, +0.0617 cm/s. The base odd-even amplitude is therefore the
catalog's, cells 1 and 2 at -18.7 and -10.6 of the wind flux. `rho v r^2` at
cell 401 is 2759848546092.128 in both, so the mass-loss rate is unchanged to
the last digit. The transit equivalent widths were NOT re-synthesized: they
are computed from a profile that is the catalog's to 2e-12 relative, so a
synthesis would measure the transit tool's own rounding and not this change.

**The regression bases** `wasp_full`, `wasp_he23off`, `mol_base_handoff` and
`lower_profile`, on scratch copies, single-threaded, with their own `maxsteps`
where they have one, control against repair:

| case | steps, control | steps, repair | outputs |
|---|---|---|---|
| `wasp_full` | 12251 | 12251 | identical but for the provenance lines |
| `wasp_he23off` | 12236 | 12236 | identical but for the provenance lines |
| `mol_base_handoff` | 12000 | 12000 | identical but for the provenance lines |
| `lower_profile` | 12000 | 12000 | identical but for the provenance lines |

`Hydro_ioniz.txt`, `Ion_species.txt` and both `_adv` products agree on every
row; the differing lines are `# provenance: git=... run=<timestamp>` and
`# source git=...`, which name the build and the wall clock. The step counts
reproduce those of the goldens' own logs. **No golden is touched**, and the
question the paragraph above left open is answered: neither trajectory ever put
a state in the band where the two spread functionals fall on opposite sides of
their thresholds.

**The 0.02 transient re-solved**, `Numerical flux: HLLC` and `ROE`, 40 passes
each, in `LHS1140b/models/.L26/x002_hllc` and `.L26/x002_roe`, from the stalled
written state, under the repaired boundary (which is bitwise the entry text's
on that state). **Neither certifies, `info = 1` in both.**

| | HLLC | ROE |
|---|---|---|
| mass row | 2.906E-01 at cell 1, tolerance 1.2E-06 | 5.947E-06 at cell 382, tolerance 3.0E-12 |
| momentum row | 1.441E-08 at cell 2, tolerance 1.0E-08 | 6.076E-07 at cell 500, tolerance 1.0E-08 |
| energy row | 6.299E-01 at cell 1, tolerance 1.0E-06 | 1.773E+00 at cell 6, tolerance 1.0E-06 |

Under HLLC the base rows that held this transient do NOT close: mass and energy
still refuse it at cell 1, five and six decades outside their tolerances. That
is the expected answer, since the closure is bitwise inert on the state and the
obstruction is the base layer's own residual. Under ROE the base mass row DOES
close, the worst mass row moving to cell 382 and falling from 2.9E-01 to
5.9E-06 and the momentum row moving to the outer boundary; what remains is the
energy row at cell 6, and the base sits at a different thermal structure, the
ghost reaching 226.02 K, `base.inp`'s `T_base` to the digit it is given in,
with cells 1 to 5 at 220.4 to 203.5 K against 470 to 526 K under HLLC. That
contrast is the flux family of item L25 step 3 seen on this transient, and it
is that item's question.

**What a golden refresh will move.** Nothing, unless a marching trajectory
passes through a state whose two spread functionals fall on opposite sides of
their thresholds. Both are saturated on every written state measured in
section R4, and the ratio of the two functionals over the 302 certified states
runs from 3.28 to 19.7, so a state can only straddle if its range sits in
[2e-2, ...] while its standard deviation sits in [2e-3, 4e-3], i.e. if its
ratio exceeds 5 to 10 at a spread ten times the certified band. That is a
transient property and it cannot be decided from written states; the
regression measurement above is what decides it.

## R8. The option keys

| key | kept or removed |
|---|---|
| `EXHALE_BASE_BRANCH_ON_CELL1` | KEPT. It restores the local discriminant and is the control for every measurement of what the window is worth. |
| `EXHALE_BASE_MACH_BLEND` | KEPT. It is the width of the reversal handover, untouched by this item, and its own operating-point measurement is still live (section 13 above). |
| `EXHALE_BASE_WIND_SPREAD` | KEPT, and it now overrides the recalibrated threshold `2e-3` of the new functional. |
| `EXHALE_BASE_WIND_SAY_ON_MACH` | REMOVED, with `base_wind_says_on_its_mach`. It keyed the window's say on a Mach number at `base_face_mach_blend`, which is exactly the statement the amplitude weight now carries at its own scale; keeping both would give the same quantity two widths again, which is the defect item L21 removed. |
| `EXHALE_BASE_WIND_AMPLITUDE` | NEW. It overrides the regularization scale, which is how section R5's sensitivity is measured. |

No `input.inp` key changes, so `docs/input_schema.md` is unchanged: these are
environment controls for measurements and that file documents no environment
control.

## R9. Noticed outside the item, reported and not fixed

**The base face flux the Riemann solve produces is not the wind's flux on any
of the four LHS 1140 b wind states**, MEASURED in section R3: +32.4, +336.6,
-404.3 and -584.5 of the window's mean flux at the base face, with the faces
above it at -18.6, -203.7, -153.9 and -232.5 before settling near +1 at faces
4 to 8. Mass is conserved, so a stationary state carries one flux through every
face and the base layer of these states does not. They are certified IN THE
WIND, `r >= 1.2`, and the mass row below that is read against the rounding
anchor of cells whose density is 1e14, so an imbalance of hundreds of times the
wind's own flux passes it. That belongs with the base odd-even mode of
`docs/lhs1140b_stationary_L25_20260916.md` section 3 and with
`docs/p44_base_sawtooth.md`, not with this item's closure.

The same measurement **refutes, on the catalog states of 2026-09-16, the
sentence item L21 wrote about the fiducial**: that the Riemann face mass flux
is +1.00000 `F_wind` at every face of the grid including the base face
(`docs/lhs1140b_stationary_L21_20260915.md` section 3). It was true of the
state L21 measured on 2026-09-15; it is not true of the state
`atomic_scalar_gj1132_kzz1e9/HeH2.13` now carries, whose base face reads
+32.4. That memo was not edited, since it records a measurement made on a state
that has since been re-solved; the correction belongs here, with the state it
was re-measured on.

**`flux_spread_of_state` still reads `(max - min)/|mean|`** and so still has
the argmax kink this item removed from the boundary. That is admissible where
it is: it is a convergence gate, evaluated once on a finished state and never
differenced. The two functionals are no longer the same one, and the design
note at `base_wind_window_spread` says so. Not changed: another file, another
item.
