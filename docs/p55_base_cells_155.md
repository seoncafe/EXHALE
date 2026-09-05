# Section 155 draft: the base cells are reproducible, and the residual is not a function of the state alone

Measured 2026-09-03 on two builds:

* **Phase C** -- `.../scratchpad/phaseC/src` copied read-only into
  `.../p55base/pc155` and rebuilt; the clean build reproduces Phase C's own
  binary bit for bit (`583aaa98`), so the boundary being measured is theirs.
  Roots: `phaseC/hotU_root` and `phaseC/wasp_newton`, restarted from their own
  output.
* **Pre-Phase-C** -- the tree as of section 145, copied into `.../base123`.
  Roots: the two accepted states of section 144.5 (`s11/hotU_cellmax`,
  `s11/WASP_cellmax`).

Probes (scratch only, both env-gated, both off by default):
`EXHALE_SWEEP_CONSERVE` reports what one sweep does to `rho`, `n_H` and `n_He`;
`EXHALE_RESID_REPRO=<n>` repeats `eval_residual` on one accepted `Y` (negative
`n` restores the composition before each repetition). Patch:
`p155_selfconsistent.diff`.

---

## 1. The premise of item (4) does not reproduce

Section 11.6 recorded that re-running the ionization sweep on an accepted state
moves the mass residual of cells 1 and 2 from `9e-6` to `2.7e-1` and `8.6e-1`,
and attributed it to `ioniz_eq` rewriting `rho` through `calc_rho`.
**Neither half of that survives measurement.**

**`calc_rho` does not move `rho`.** Over every sweep of both roots on both
builds, the density the sweep writes back and the density it was handed differ
by at most:

| build | state | max abs `drho/rho` | max abs `dn_H/n_H` | max abs `dn_He/n_He` | sweeps |
|---|---|---|---|---|---|
| Phase C  | hot Uranus | 8.13e-16 | 4.67e-16 | 4.32e-16 | 57 |
| Phase C  | WASP-121b  | 1.71e-15 | 2.22e-16 | 3.04e-16 | 86 |
| pre-C    | hot Uranus | 9.68e-16 | 3.64e-16 | 4.42e-16 | 346 |

This is roundoff, and it is roundoff for a reason: the sweep's unknowns are
fractions of the element totals `n_H`, `n_He` and `n_M,tot` that it computes
from the state it is handed, and every write-back is one of those totals times a
fraction. Because the code's masses are the nucleus counts exactly
(H2 = 2 m_H with 2 H nuclei, H3+ = 3 with 3, HeH+ = 5 = 4 + 1), conserving
nuclei conserves mass identically. The `max(...,0.0d0)` clamps on `nhi` and
`nhei` in the molecular write-back CAN break this -- they are the one place the
element budget is not closed by construction -- but on these states they do not
bite.

**The base cells' residuals are not O(1).** Repeating `eval_residual` on one
accepted `Y`, cells 1-3, first evaluation against the twentieth:

| build | state | cell | mass row, eval 1 | eval 20 | ratio |
|---|---|---|---|---|---|
| Phase C | hot Uranus | 1 | 7.904e-05 | 7.765e-05 | 0.982 |
| Phase C | hot Uranus | 2 | 2.521e-05 | 2.470e-05 | 0.980 |
| Phase C | WASP-121b  | 1 | 1.1725e-07 | 1.1725e-07 | **1.000** |
| Phase C | WASP-121b  | 2 | 1.6359e-07 | 1.6359e-07 | **1.000** |
| pre-C   | hot Uranus | 1 | 2.553e-07 | 2.491e-07 | 0.976 |
| pre-C   | hot Uranus | 2 | 8.154e-08 | 7.762e-08 | 0.952 |
| pre-C   | WASP-121b  | 1 | 9.476e-09 | 5.163e-09 | **0.545** |
| pre-C   | WASP-121b  | 2 | 2.177e-08 | 2.010e-08 | 0.923 |

The largest excursion anywhere is a factor two, at WASP-121b cell 1 before
Phase C. **Under Phase C that same cell is reproducible to four significant
figures.** The base cells are not the problem, and to the extent they were one,
Phase C's boundary closed it: it is the prescription that removed the base
cell's dependence on the previous sweep's composition (no `n_part_cell1`, no
valve).

---

## 2. What is real: the residual is a function of (state, composition)

`f_sp` is `intent(inout)` in `ioniz_eq` and the residual mutates it. It is not
part of the Newton unknown `Y`. So `eval_residual(Y, f_sp)` is a function of
two arguments, only one of which the state carries.

**It is deterministic in both.** With the composition restored before every
repetition (`EXHALE_RESID_REPRO=-4`), all four states give bitwise identical
residuals on every repetition. Nothing else -- not the frozen chemical
background, not the excited-hydrogen populations, not the stored WENO
smoothness factors -- contributes anything at all. That was worth establishing
before looking further.

**But the number reported at acceptance is evaluated at the WRONG composition.**
In the line search the accepted trial is evaluated as
`eval_residual(Ytry, f_sp_j)` with `f_sp_j` copied from the PREVIOUS iterate, and
`rnorm` is taken from that `F`. The state handed back and written to disk carries
the composition the accepted trial LEFT. Evaluating the same `Y` at its own
composition, and iterating the sweep at fixed `Y` until it stops moving:

| build | state | reported `\|\|R\|\|` | own composition, 1 sweep | at the composition fixed point | sweeps | factor |
|---|---|---|---|---|---|---|
| Phase C | hot Uranus | 2.211e-06 | 1.283e-05 | 1.257e-05 | 3 | **5.7** |
| Phase C | WASP-121b  | 6.287e-06 | 3.865e-05 | 3.593e-05 | 9 | **5.7** |
| pre-C   | hot Uranus | 4.141e-06 | 1.762e-05 | 1.722e-05 | ~10 | **4.2** |
| pre-C   | WASP-121b  | 6.264e-06 | 1.737e-04 | 3.539e-04 | ~14 | **56** |
| pre-C   | `wasp_full_newton` | 8.873e-06 | 1.851e-04 | 3.412e-04 | ~10 | **38** |

The composition at fixed `Y` has its own fixed point and reaches it in 3 to 14
sweeps; one evaluation is still short of it by up to a factor two. **Phase C cuts
the WASP-121b factor from 56 to 5.7 and leaves the hot Uranus at about 5.7.**

### What this does and does not change about acceptance

It does not, today, change which states are accepted. Only 7 of the 63
regression cases set `Resid tol` at all; without it `resid_th = -1`, the
residual gate is off (`||R|| < -1.000E+00` in the acceptance line), and every one
of the roots measured here was accepted on the FLUX GATE alone. The flux gate
reads the conserved face mass flux of `u` and is composition-independent, so it
is untouched by all of this.

What it changes is what the reported number means. `||R||` is quoted in run
logs, in `docs/`, and in item (AD)'s own measurements, and on these states it
understates the residual of the state that was written by 4 to 56 times. Every
one of the five states above is above `1e-5` at its own composition, so a run
that DID set `Resid tol: 1.0e-5` would be accepting a state that does not meet
it.

---

## 3. Prescription

**(a) The gate and the report read the residual of the state that is written, at
the composition that state carries.** `p155_selfconsistent.diff` adds
`resid_at_own_composition` (default `.false.`, `n_selfconsistent_max = 25`): at
hand-back, re-evaluate at `(Y, f_sp)` and iterate the sweep at fixed `Y` until
`||R||` changes by less than `1e-3` relative. Measured cost: 3 sweeps on the hot
Uranus, 9 on WASP-121b, one-off. The coupled route already does a single
re-evaluation for its carrier half (`nvar_jac >= 4`); this is the same idea taken
to the fixed point and applied to all routes.

Default off because turning it on makes every root accepted under `Resid tol`
report a residual 4 to 56 times larger, which is a physics decision about what
"converged" means, not a refactor. The recommendation is to turn it on: a steady
state is a joint fixed point of the hydrodynamic balance and the chemical
equilibrium, and a measure that reads the hydro balance under last iterate's
chemistry is not measuring that.

**(b) `rho` should be an input to the sweep, not an output.** The measurement
says this changes no number today -- the write-back is element-conserving on
every state tested, to 1e-15. It is still the right shape: the sweep determines
the composition at a given `(rho, T)` and the particle count that follows; mass
is the hydrodynamics' variable. Making `n_io` `intent(in)` and normalizing
`f_sp` by the input density would make that structural instead of accidental,
and it would turn the `max(...,0.0d0)` clamps in the molecular write-back --
the one place the element budget is NOT closed by construction -- into a
detectable condition rather than a silent mass adjustment. It is a
golden-refreshing change at the 1e-15 level, so it belongs in a change series
that is refreshing them anyway.

**Not adopted here:** anything that would make the composition part of `Y`. That
is the coupled route's design and it exists (`nvar_jac >= 4`, the carrier
unknown); extending it to the full composition is a much larger question than
this measurement supports.

---

## 4. Contact with Phase C

The patch touches `src/modules/time_step/steady_newton.f90` (the hand-back
block, after `unpack_carrier`/`Apply_BC` and before `steady_gates_met`) and adds
two parameters in `src/modules/init/parameters.f90`. **No boundary code is
touched**: `Apply_BC.f90`, the new `base_boundary.f90`, `init.f90` and
`utilities.f90` are exactly Phase C's. The one place the two meet is the
hand-back `Apply_BC(u)` call, which Phase C rewrites and this patch only calls;
if Phase C changes the call signature there, this patch follows it and nothing
else.

The probes are separable and belong in scratch: `EXHALE_SWEEP_CONSERVE` in
`ionization_equilibrium.f90` and `EXHALE_RESID_REPRO` in `steady_newton.f90`.
The conservation probe is worth keeping in some form -- as a cheap assertion
that the sweep's write-back closed the element budget, reported once per run
rather than printed per sweep.

---

## 5. What is not established

* One accepted state per planet per build, all with the residual gate off. A
  case that actually sets `Resid tol` and converges under it has not been
  re-measured at its own composition.
* Whether the composition fixed point at fixed `Y` always exists and is unique.
  It was reached in 3 to 14 sweeps on all five states, monotonically on the two
  that start below it and with one overshoot on the three that do not
  (WASP-121b under Phase C: 3.86, 4.36, 3.79, ... 3.59e-05), but nothing here
  proves existence or uniqueness.
* Whether accepting on the self-consistent residual changes the mass-loss rate.
  It cannot change the state that is written -- the patch only measures -- but a
  run that keeps solving instead of stopping will end somewhere else.
* The `max(...,0.0d0)` clamps: shown not to bite on these five states, not shown
  never to bite. The conservation probe is what would say.
