# Ionization-equilibrium roots validated against the physical simplex

*2026-08-11. Diagnosis raw material: the session scratchpad
`/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/31dad923-90fb-485a-bbf3-1cf82bd83337/scratchpad/negden_diag/`
-- negative-density scans of the `wasp_full` golden, a
secondary-ionization-off control run, and the heating breakdown across the
affected band. That directory is transient; the numbers it produced are
reproduced below.*

## Symptom

`Ion_species.txt` carries negative number densities in a thin band just above
the base:

| run | cells | radii [R_p] | columns |
|---|---|---|---|
| `wasp_full` golden | 17 (i = 55-71) | 1.0106-1.0140 | H II, O II, O III, Fe I |
| WASP-121b `output/` | 15 (i = 55-69) | 1.0106-1.0136 | H II, O II, O III, Fe I |
| HD189733b `output/` | 5 | 1.0054-1.0064 | H II, O II, O III, Fe I |

Fe I is the worst: in WASP-121b its most negative value is 146 times the
largest positive Fe I density anywhere in the domain, so the neutral-iron
profile in that band is meaningless rather than merely noisy. Re-running
HD189733b single-threaded from that same stored state reproduces the band at
10 cells rather than 5, so the count is trajectory-dependent; the mechanism is
not.

## Diagnosis

The equilibrium systems are polynomial in the stage fractions and possess
roots outside the physical simplex. At the step where secondary ionization is
switched on (step 7186 of `wasp_full`), MINPACK `hybrd1` returned `info = 1`
on such a root: x(H II) = -5.2e-6, and the Fe stage fractions summed to
1.2209, i.e. an Fe I density 22% below zero. Nothing tested the returned root,
so it was accepted, and the next step's warm start re-seeded the solve from
those same values. The state therefore reproduced itself for the remaining
6302 steps.

No branch of `ioniz_eq` tested the physics of the root it accepted. The
merged He 2^3S + metals branch was the most exposed -- it called `hybrd1`
directly and did not even read `info` -- but the metals-without-triplet branch
turns out to produce inadmissible roots as well (35 of them over a
`wasp_he23off` run, see the gates below); there they healed before the final
snapshot, which is why that golden looks clean. The molecular branch was the
only one with any defense: it retries from a chemical-equilibrium guess and
clips the extraction.

The consequence is not cosmetic, and it reaches the heating through two
separate routes.

* The metal photoheating integrand is built from the neutral metal densities,
  so a negative Fe I turns it negative: in the `wasp_full` band the metal
  channel came out at -6.6e-7 erg cm^-3 s^-1. A negative photoheating rate is
  impossible.
* The ionized fraction of the H + He nuclei that drives the Spitzer & Scott
  (1985) secondary-ionization partition is clamped into [0,1], so a negative
  H II sends it to exactly zero. The heating fraction is
  f_heat(x) = 0.9971 (1 - (1-x^0.2663)^1.3163), and f_heat(0) = 0, so *every*
  photoheating channel above the photoelectron threshold is switched off, not
  only the metal one.

Together they left the total heating 93% below its neighbors on both sides of
the band (7.2e-6 outside, 5.4e-7 inside). The second route is why clamping the
offending fractions to zero is not a repair: the clamp is exactly what
produces f_heat = 0. The cure has to be a physical root, not a clipped one.

A control run with `Secondary ionization: False` produced no negative
densities at all, which locates the trigger but not the cause: the cause is
that an inadmissible root is accepted whenever the solver lands on one.

## Fix

All in `src/modules/radiation/ionization_equilibrium.f90` unless noted.

1. **`ionization_fractions_physical`** -- tests a root against the simplex
   every unknown lives in: each stage fraction >= 0, and, for each element
   separately, the tracked ionized stages summing to <= 1 (the neutral stage
   is the remainder). Comparisons carry a 1e-10 tolerance so a root sitting on
   a face of the simplex is not rejected for round-off; the violations seen
   here are 1e-6 to 2e-1.

2. **`ionization_balance_at_fixed_ne`** -- each element in its own
   photoionization + electron-impact ionization against radiative
   recombination balance at the electron density of the incoming state, the
   couplings between elements dropped. With u_k = gamma_k + beta_k n_e out of
   stage k and d_k = alpha_k n_e back into stage k-1, the stage populations of
   a three-stage element are proportional to (d_1 d_2, u_0 d_2, u_0 u_1),
   which is non-negative and normalizes to one: the result is inside the
   simplex by construction, whatever state it replaces. The He 2^3S fraction
   is the steady state of the metastable level at those populations, capped by
   the neutral He fraction. Rates are read from the same cell state the
   residuals read (`ieq_cell`, and the `met_*` coefficients), so this is the
   cell's own physics rather than a generic guess.

3. **Root validation with physically seeded retries.** Every atomic cell solve
   now tries up to three starting points -- the previous state, the ionization
   balance of (2), and the optically thick limit with all nuclei neutral --
   and keeps the best admissible root, ranked (physical and converged) >
   (physical) > none. Attempt 1 is accepted as it stands whenever it converges
   to a physical root, so a healthy cell takes exactly the same solver path as
   before. If no attempt produces a physical root the cell falls back to the
   ionization balance of (2), which is admissible, and is counted.

4. **Warm-start self-sticking blocked.** If the stored previous state is not
   physical it no longer seeds the solve; the ionization balance is used
   instead. This is what breaks the loop that kept the bad state alive for
   6302 steps.

5. **`solve_ieq`** (`nonlinear_system_solver/newton_solver.f90`) gained an
   optional `converged` argument so the caller can tell an accepted solution
   from a solver that gave up. Existing callers are unaffected.

6. **Absorbed energy sign** (`radiation/util_ion_eq.f90`). `q_abs` is an
   absorbed energy rate and cannot be negative; the old guard
   `Hea_1/max(q_abs, 1e-99)` only protected against zero and would have turned
   a negative q_abs into a huge efficiency. The efficiency is now zero unless
   q_abs > 0. Where q_abs = 0 the heating integral vanishes with it, so the
   value is unchanged for every physical state.

Reporting is in two places. A line per equilibrium sweep names the step and
counts the cells whenever a state left the simplex -- a stored state rejected,
a first root outside the simplex, or a cell with no admissible root at all.
Restarting a cell merely because the solver did not reach its tolerance is
routine and stays silent there; it is counted in a run-wide total printed at
the end of the run next to the Newton/hybrd1 usage line
(`src/EXHALE_main.f90`). Both are silent for a run that never leaves the
simplex.

## Scope

The molecular network keeps its own chemical-equilibrium retry and is left
alone; the H-only system (no helium) and the post-processing advection solves
in `post_process_adv.f90` are not covered by this validation.

That last exclusion is visible in the production outputs. Measured 2026-08-11
after the re-convergences below, the equilibrium `Ion_species.txt` of all four
paper planets is free of negative densities, but the advected
`Ion_species_adv.txt` of HD209458b still carries 42 negative entries in the
first 14 cells above the base (r = 1.0000-1.0033 R_p): H II down to
-3.7e7 cm^-3 against a base H density of ~1e14, plus He III and He 2^3S at
-7.6e-6 and -7.3e-7 cm^-3. HD189733b, WASP-52b and WASP-121b are clean in both
files. Those entries come from the advection correction, not from the
equilibrium solve, so they are outside what this change addresses.

## Gates

All runs single-threaded (`OMP_NUM_THREADS=1`), so every comparison below is
deterministic. The tree carried unrelated uncommitted work at the time; a
baseline `make check` on it passed all three cases byte-identically, so the
differences reported here come from this change alone.

### Regression matrix

| case | verdict | Mdot [log g/s] | steps | negative densities |
|---|---|---|---|---|
| `wasp_full` | FAIL (intended) | 13.22 -> 13.22 | 13488 -> 13483 | 68 entries / 17 cells -> 0 |
| `wasp_he23off` | FAIL | 13.22 -> 13.22 | 13482 -> 13478 | 0 -> 0 |
| `mol_base_handoff` | PASS, byte-identical | - | 12000 | 0 -> 0 |

Goldens were **not** re-snapshotted.

`wasp_full` is the case the defect lives in. Away from the affected band the
solution barely moves: the largest relative differences against the golden are
0.89% in T, 0.77% in rho, 1.3% in v and 0.13% in p. Inside the band the
heating changes by up to 93%, which is the point of the exercise. The metal
photoheating channel is positive everywhere in the new run (domain minimum
+1.03e-7, against -8.9e-7 in the golden), and the total heating across the
band is continuous with its neighbors (7.06e-6 in the band, against 5.36e-7
in the golden and ~7.1e-6 on either side).

`wasp_he23off` was expected to be untouched, since its golden shows no
negative densities. It is not: over its 13478 steps the run reports 35 roots
outside the physical simplex and 9 stored states rejected, with 9 cells that
no starting point could resolve (left on the ionization balance). Those
violations were previously accepted in silence and healed before the final
snapshot, which is why the golden looks clean. The steady state is
insensitive to them -- Mdot unchanged at 13.22, no negative densities either
way, and maximum relative differences of 0.28% in T, 0.086% in rho, 0.46% in
v and 0.039% in p. The run also restarted 12675 cell solves (0.19% of ~6.8
million) that had merely failed to reach the solver tolerance; those restarts
are kept only when they produce a converged, physical root, so they can only
improve the root, but they do perturb the trajectory and contribute to the
golden difference.

`mol_base_handoff` exercises the molecular branch, which this change does not
touch, and reproduces its golden bit for bit.

### Where the validation fires

In `wasp_full` the sweep summary fires at exactly two places in a 13483-step
run, and the run totals are 30 cell solves restarted, 21 roots outside the
simplex, 5 stored states rejected, 4 cells left on the ionization balance:

* steps 585-589 -- one cell per step, which no starting point could resolve
  for four of those steps; it clears on its own.
* step 7183 -- 17 roots outside the simplex, every one of them resolved by a
  restart, none left unresolved. This is the secondary-ionization activation
  step (7186 in the golden trajectory, which had already drifted by three
  steps), and 17 is exactly the number of cells that used to end up with
  negative densities.

Nothing fires afterwards: the state never leaves the simplex again, so the
self-sticking is gone.

The contrast with `wasp_he23off`, which restarts 12675 cell solves against
`wasp_full`'s 30, is in the solver rather than the physics: with the He
triplet off, the metals system goes through `solve_ieq`, and a cell counts as
unconverged only if both the analytic-Jacobian Newton and the `hybrd1`
fallback miss the tolerance, which happens for roughly one cell per step
there. Those restarts are accepted only when they yield a converged physical
root.

### Production runs

Both were re-run from their own `input.inp`/`metals.inp` (and, for HD189733b,
its stored IC), once with the pre-change binary and once with the fix, so the
comparison isolates this change.

| | WASP-121b | HD189733b |
|---|---|---|
| Mdot [log g/s] | 13.17 -> 13.17 | 9.05 -> 9.05 |
| steps | 8228 -> 8228 | 2002 -> 2002 |
| negative densities | 60 entries / 15 cells -> 0 | 40 entries / 10 cells -> 0 |
| min heat_metals | -9.64e-7 -> +7.98e-8 | -2.24e-7 -> +3.54e-36 |
| max rel diff T, rho, p | 8.0%, 10.1%, 7.0% | 11.6%, 9.8%, 2.6% |

Neither run's mass-loss rate moves, and in both the metal photoheating is
non-negative everywhere afterwards. In WASP-121b the total heating across the
band recovers from 1.60e-6 to 1.03e-5 erg cm^-3 s^-1, continuous with the
1.03e-5 of the neighboring cells. The profile differences are concentrated
in and just above the affected band.

HD189733b is the clearest demonstration of the self-sticking: its stored
initial condition is one of the poisoned outputs, and at step 1 the run
rejects 5 stored cell states -- exactly the 5 cells that carried negative
densities -- and resolves all of them. One further restart at step 2, then
silence.

The remaining negative entries in HD189733b `Hydro_ioniz.txt` are negative
velocities at the base (24 cells before, 15 after). They are the known base
inflow of that configuration, not a product of the ionization solve.
