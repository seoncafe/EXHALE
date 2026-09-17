# L4b: why the flux-closure rung's wind is not certified, and what a recipe would take

Item L4b of `docs/PLAN_20260913_lhs_stationary.md`. The rung measured is
`LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH9`, iteration 0
(`k00`), solved by the recipe of `LHS1140b/MODELS.md` section 6. Every
number below is MEASURED on this tree unless it is marked READ. No source
file was changed; the measurements are configuration changes and
environment variables only.

## 1. Verdict

**The rung is refused by the outer loop, and what refuses it is the
hydrodynamic energy row of the second and every later outer pass. The
cause is not the profile lower boundary and not the helium-rich
composition: it is the C, N, O reservoirs. With them on, the radiative
cooling between 1.15 and 1.5 R_p is 79 to 94 percent neutral-carbon line
emission and stands at 0.53 to 0.66 of the heating, against 0.04 in the
certified metals-off case at the same cell; the elemental relaxation of
one outer pass moves X_He there by 0.4 to 2.6 percent and leaves an energy
row of 3.7e-2, and the stationary solve of the next pass cannot take it
back down.**

**It cannot take it down for a reason that is in the step control, not in
the state.** The solve is entered from a state whose mass and momentum rows
are already at their floors (2.11e-9 against an anchored 8.0e-9, and
1.31e-14 against 1e-8) and whose energy row alone is off. The ledger the
stagnation detector ranks iterates on is the distance from certification,
and a Newton iterate raises the continuity row from 2.11e-9 to 6.5e-4
against a 3e-12 tolerance, so the entry state is the best iterate BY
CONSTRUCTION and no trial can improve on it. After `n_stall_best` = 20
iterations twice the solve stops as STAGNATED and hands the entry state
back unchanged (`||R||` identical to seven digits), while the merit it
actually descends on has fallen by a factor 4.8 in those same iterations.
Three such passes exhaust `outer_no_fall_max` = 3 and the run ends
`info = 1`.

**No configuration reaches a certified rung.** Smaller composition steps
(`EXHALE_DIFF_OMEGA` = 0.125) scale the insult down linearly and are
refused identically; the outer-iteration cap was never the binding limit
(the solve stopped at iteration 40 of 3000); a hundredfold stronger
pseudo-transient shift damps the continuity excursion by a factor 4.5, three
decades short of what it would take; the trust-region restart, the one
route that resets the stall counter, is unreachable because the partitioned
solve registers no species row; and the coupled route, which would remove
the alternation altogether, does not move this planet's base layer at all.
The parameters that would let the solve go on -- `n_stall_best`,
`outer_no_fall_max`, the floor of `omega` and the Krylov forcing term --
are all compile-time constants. `run_closure.sh`, `closure.json` and
`make_models.py` are therefore left as they are; the minimal change is
stated in section 7 and belongs to `src/`, which this item does not touch.

## 2. What the rung does, pass by pass

READ from `LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH9/k00/run.log`
(the run of 2026-09-13 14:24 to 14:33, 16 threads).

| pass | hydro | worst gated species row | mass | momentum | energy (cell, r) | omega |
|---|---|---|---|---|---|---|
| 1 | `info = 0` | 5.55e-4 at 255 | 2.11e-9 | 1.31e-14 | 3.95e-8 (40, 1.008) | 0.500 |
| 2 | `info = 2` | 2.64e-4 at 260 | 2.11e-9 | 1.31e-14 | 3.69e-2 (233, 1.266) | 0.250 |
| 3 | `info = 2` | 1.95e-4 at 262 | 2.11e-9 | 1.31e-14 | 4.60e-2 (234, 1.270) | 0.125 |
| 4 | `info = 2` | 1.70e-4 at 262 | 2.11e-9 | 1.31e-14 | 4.93e-2 (235, 1.275) | 0.125 |

The elemental side is converging: the gated He/H transport row falls
5.55e-4, 2.64e-4, 1.95e-4, 1.70e-4, and its relaxation reaches the fixed
point of the element operator in 30 to 32 steps in every pass. It slows
only because the progress control halves `omega` -- and it halves it
because the joint distance did not fall, which is the ENERGY row's doing,
not the element row's. Four more passes at the pass-1 rate would have
carried the element row to its 1e-5 tolerance.

## 3. Where the energy row is, and what it is made of

### 3.1 The row over the column

MEASURED with `EXHALE_RESIDUAL=1` on the state `k00` wrote
(`scratchpad/L4b/resid_k00`); it reproduces the run's own certification
(4.931e-2 at cell 235) to every digit.

| cell | r/R_p | R_energy/scale | cool/heat | T [K] | v [cm/s] |
|---|---|---|---|---|---|
| 10 | 1.0019 | +1.55e-3 | 0.012 | 759 | 0.054 |
| 40 | 1.0077 | +2.77e-3 | 0.034 | 1229 | 0.14 |
| 100 | 1.0250 | +8.29e-3 | 0.446 | 1869 | 0.52 |
| 150 | 1.0615 | -7.29e-3 | 0.562 | 2355 | 2.4 |
| 200 | 1.1489 | -2.97e-2 | 0.649 | 3869 | 26 |
| 235 | 1.2751 | **-4.93e-2** | 0.647 | 5336 | 139 |
| 267 | 1.4816 | -2.43e-2 | 0.497 | 5873 | 532 |
| 300 | 1.8571 | +1.11e-2 | 0.266 | 5423 | 1.8e3 |
| 450 | 12.726 | +3.22e-2 | 0.007 | 849 | 6.7e4 |

The row is not a defect at one cell. It is a smooth, single-signed excess
over the whole subsonic column, negative (net heating in excess of the
enthalpy-flux divergence) from 1.06 to 1.85 R_p and positive on either
side. The mass and momentum rows of the same state are 1.49e-9 and
2.35e-14, both inside their tolerances.

The signed terms at the binding cell, in code units:

```
cell 235, r = 1.27514
R = -3.292944E-06   flux divergence = 2.030526E-05
heat = 6.677633E-05  cool = 4.317812E-05  heat - cool = 2.359821E-05
scale = heat = 6.677633E-05
```

so the row is a 3.7 to 4.9 percent imbalance of two terms that themselves
cancel to a third of the heating.

### 3.2 What the cooling is

MEASURED from `output/Cooling_breakdown.txt` written by a post-processing
pass on the same state (`scratchpad/L4b/pp_k00`).

| cell | r/R_p | T [K] | C I | He I excitation | recombination | N I | C I / C |
|---|---|---|---|---|---|---|---|
| 100 | 1.025 | 1869 | 93.6% | - | 4.7% | - | 0.989 |
| 200 | 1.149 | 3869 | 93.7% | 0.9% | 3.1% | 1.0% | 0.981 |
| 233 | 1.266 | 5273 | 87.0% | 5.3% | 3.6% | 2.5% | 0.969 |
| 260 | 1.426 | 5834 | 79.0% | 10.7% | 4.8% | 3.0% | 0.944 |
| 300 | 1.857 | 5423 | 60.3% | 24.2% | 9.0% | 2.0% | 0.855 |

Neutral carbon carries the cooling of the whole band the energy row
refuses on, at C/H = 2.778e-4 (the reservoir the Photochem column hands
over) and with carbon 94 to 99 percent neutral there. The photon grid does
reach below the H I edge on this run -- `EXHALE_setup.out` states "Below
13.6 eV (down to 4.77 eV): the same SED file as every other band", the
floor being the He 2^3S threshold -- so the C I ionization threshold at
11.26 eV is inside the field and the neutral fraction is the code's own
balance, not a missing band. Whether that balance is right against an
independent calculation was NOT checked here and is flagged in section 8.

For comparison, the certified prescribed case
`atomic_scalar_gj1132_kzz1e9/HeH1.60` (metals off) has cool/heat = 0.042 at
cell 233 and 0.053 at cell 260, and its cooling has no metal channel at
all.

### 3.3 The composition change and the row are in the same band

MEASURED, He mass fraction of the seed against the state `k00` returned:

| cell | r/R_p | X_He seed | X_He pass 4 | dX/X | R_energy/scale |
|---|---|---|---|---|---|
| 1 | 1.0002 | 0.974902 | 0.973002 | -0.20% | small |
| 64 | 1.0127 | 0.974651 | 0.973448 | -0.12% | +4.5e-3 |
| 150 | 1.0615 | 0.970320 | 0.970929 | +0.06% | -7.3e-3 |
| 200 | 1.1489 | 0.957834 | 0.961763 | +0.41% | -3.0e-2 |
| 233 | 1.2657 | 0.938006 | 0.948275 | +1.09% | -4.9e-2 |
| 267 | 1.4816 | 0.899357 | 0.922333 | +2.55% | -2.4e-2 |
| 300 | 1.8571 | 0.851214 | 0.884835 | +3.95% | +1.1e-2 |

The band where the elemental relaxation first moves the composition by a
percent, 1.15 to 1.5 R_p, is the band where the energy row binds, and it is
the band where C I carries the cooling. Below 1.03 R_p, where the
composition moves by 0.1 to 0.2 percent, the row is thirty times smaller.
`docs/figures/lhs1140b_L4b_energy_row.pdf` puts the three panels on one
radius axis.

## 4. Why the stationary solve returns `info = 2`

READ from the pass-2 block of `k00/run.log` and from
`src/modules/time_step/steady_newton.f90`.

The step control descends on the merit `||F/Drow||_2`; the ledger that
decides which iterate is "best", and the stagnation detector that reads it,
rank iterates on `distance_from_certification` -- each row over the
tolerance it carries at that cell, continuity against 3e-12 where the
cell's arithmetic resolves it (steady_newton.f90,
`hydrodynamic_distance_from_certification_by_cell`, and the comment there
says the two functionals are different). The pass-2 trajectory:

| point | judged distance | which row sets it | `\|\|R\|\|` | merit `\|\|Fs\|\|_2` |
|---|---|---|---|---|
| entry | 3.685e+04 | energy 3.685e-2 / 1e-6 | 3.685e-2 | 1.50e-2 |
| iteration 1 | 2.181e+08 | continuity ~6.5e-4 / 3e-12 | 3.031e-2 | 1.42e-2 |
| iteration 20 | (restart to best) | | 3.685e-2 | 7.00e-3 |
| iteration 40 | 1.354e+08 | continuity | 3.685e-2 (best) | 3.15e-3 |
| returned | 3.685e+04 | the entry state | 3.685e-2 | |

and the last iterate's own rows over the subsonic column read
`mass = 3.366e-02`, `mom/grav = 5.517e-08`, `energy = 1.356e-01`, the
largest of them at cell 1 -- the continuity row has risen seven decades from
the 2.11e-9 the entry state carried.

That is the whole mechanism. The base cell of this planet runs at Mach
3e-7; its continuity row is a difference of two nearly equal face fluxes
measured against a scale that bounds the row's largest term, so holding it
at 1e-9 is a statement about the base velocity at the 1e-9 level, and any
Newton transient breaks it long before it is restored. Pass 1 could climb
that hill because its entry distance was 3.317e+11 and every iterate
improved on it. Pass 2 cannot, because its entry distance is 3.685e+04 --
four decades BELOW the plateau any Newton transient sits on -- so
`n_since_best` never resets, `n_stall_best` = 20 fires twice, and
`solve_steady_jfnk` returns the state it was given.

The solve is not out of iterations: the cap was 3000 and it stopped at 40.
It is not out of descent either: the merit fell by 4.8 over those 40
iterations and was still falling.

**Why a Newton step puts its error in the continuity row.** The linear
solve of the partitioned route is inexact by construction: `solve_steady_jfnk`
calls `pgmres` with the relative tolerance written in as the literal
`1.0d-1` (steady_newton.f90, the two calls in the `.not. use_tr` branch of
the step block), and the log confirms every one of the 40 solves reached
that tolerance. GMRES minimizes the 2-norm of the SCALED residual, and in
that norm the entry state's energy rows stand at 3.7e-2 while its
continuity rows stand at 2.1e-9, so the ten percent the linear solve is
allowed to leave costs nothing where it is left in the continuity rows and
everything where the certification then measures them, against 3e-12. The
two functionals of the code's own comment -- the merit on the Newton's
column scales, and the certification's row maxima against row tolerances --
are as far apart here as they can be made.

The judged distance is not decomposed in the log; the arithmetic names the
row. 2.181e+08 x 3e-12 = 6.5e-4, while the energy row of that same iterate
is 3.031e-2/1e-6 = 3.0e+04, three decades below.

## 5. The experiments

Each on a scratch copy of `k00` (the same `input.inp`, `base.inp`,
`lower_atmosphere_profile.dat` and the same mapped seed as the rung), 8
threads, `EXHALE_PTC_DTAU0=1.0`. The scratch copy reproduces the rung's
pass 1 to every digit (`info = 0`, `||R|| = 3.952e-08`).

| # | what | result | verdict |
|---|---|---|---|
| a | `EXHALE_DIFF_OMEGA=0.125` (a quarter of the composition step) | pass 1 `info = 0` unchanged; pass-2 entry energy row 8.506e-3 at cell 232 against 3.685e-2 at omega 0.5, i.e. linear in omega; pass 2 STAGNATED at 40 iterations, best judged distance 8.506e+03 = the entry, entry state returned | NO HELP |
| b | `EXHALE_JFNK_MAXIT` (more iterations) | the cap in force was already 3000 and pass 2 stopped at iteration 40 on `n_stall_best`, a `parameter` in `steady_newton.f90` with no environment override, so the cap cannot be what stops the solve. MEASURED with the variable set to 400 in run (b'): pass 1 ended at the same `info = 0` and `\|\|R\|\| = 3.952e-08`, and pass 2 followed the same trajectory as the control | NOT THE BINDING LIMIT |
| b' | `EXHALE_TR_RESTART_STALL=10`, `EXHALE_TR_RESTART_WHAT=rtca` (the only configuration route that resets `n_since_best` without returning to the best iterate, so the solve keeps descending on the merit), `EXHALE_JFNK_MAXIT=400`, `EXHALE_OUTER_PASSES=2` | no restart fired in 24 outer iterations of pass 2 and the trajectory matched the control value for value (distance 2.181e+08, 1.887e+08, 4.124e+08, 1.711e+08, ... at iterations 1, 2, 3, 4). The block that takes the restart is guarded by `use_tr`, and `use_tr = (nspec_row .gt. 0)`: the partitioned route registers no species row, so it is never entered | INERT ON THIS ROUTE |
| c | `EXHALE_OUTER_PASSES=1` (the composition held at the seed's) | hydro `info = 0`; mass 2.11e-9, momentum 1.31e-14, energy 3.95e-8, all inside; the ONLY refusing entry is the elemental transport He/H partition, 5.55e-4 against 1e-5 at cell 255 | the hydrodynamics of this rung is solvable; what is not is the coupled fixed point |
| d | the `HeH2.09` rung of the same group (the composition of the stored `lower_profile`) | pass 1 `info = 0` (mass 2.41e-9, momentum 4.95e-14, energy 2.71e-8, element row 1.11e-3 at cell 255, relaxation drift 1.68e-1 in 34 steps); pass 2 entered at `||R|| = 8.309e-2`, STAGNATED at 40 iterations at a judged distance 6.757e+07, best 8.309e+04 = the entry, and returned the entry state (`done info = 2`) | THE SAME FAILURE AT He/H = 2.09 |
| e | `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13`: the certified scalar family's own configuration with the C, N, O reservoirs added and NOTHING else changed (scalar base, He/H = 2.13, no profile) | pass 1 `info = 0` (mass 6.02e-9, momentum 4.19e-13, energy 1.69e-7); pass 2 STAGNATED at judged distance 1.777e+08, best 8.643e+04 = the entry, `||R|| = 8.643e-2`, energy row at cell 243, r = 1.3165 | THE ISOLATING CONTROL: the same failure, with no profile and at He/H = 2.13 |
| e' | `atomic_scalar_gj1132_kzz1e9/HeH9.7` (metals off, helium-rich) | the campaign had not reached this case; NOT MEASURED | open |

The `HeH2.09` rung (d) was run to its own ending and reaches the same one as
`HeH9`, four passes and `info = 1`:

| pass | hydro | worst gated species row | energy (cell) | omega |
|---|---|---|---|---|
| 1 | `info = 0` | 1.11e-3 at 255 | 2.71e-8 | 0.500 |
| 2 | `info = 2` | 5.95e-4 at 264 | 8.31e-2 | 0.250 |
| 3 | `info = 2` | 4.61e-4 at 267 | 1.01e-1 | 0.125 |
| 4 | `info = 2` | 4.09e-4 at 268 | 1.08e-1 | 0.125 |

ending "REFUSED -- the joint distance of the state ... has not fallen in 3
consecutive passes", NOT CERTIFIED on the same two entries. Its element row
falls 1.11e-3, 5.95e-4, 4.61e-4, 4.09e-4, the same shape as `HeH9`'s. So
the refusal does not depend on how helium-rich the reservoir is: it is
there at He/H = 2.09, 2.13 and 9.01 alike.

Experiment (e) is what makes the verdict of section 1 a measurement rather
than an inference: `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13` differs from the
certified `atomic_scalar_gj1132_kzz1e9` family only by the three element
ratios in its `base.inp`, and it fails at the same outer pass, in the same
row, in the same band of radii, with the same ledger signature.

The pseudo-transient shift was measured on the state the rung returned
(`||R|| = 4.931e-2`, energy row 4.93e-2 at cell 235, continuity row
1.49e-9), reloaded and re-solved with `EXHALE_OUTER_PASSES=1`. That state
has exactly the structure every pass after the first is handed, and outer
pass 4 of the rung itself entered from it with the same numbers, so the
rung's own pass 4 is the control and the comparison is exact. Row (g) is a
different experiment: it starts from the mapped seed, like the rung's own
pass 1, because the coupled route is a different partition of the unknowns
and not a re-solve of one state.

| # | what | first Newton step's judged distance | verdict |
|---|---|---|---|
| f (control) | `EXHALE_PTC_DTAU0=1.0` -- outer pass 4 of the rung | entry 4.931e+04 -> **2.853e+08** | stagnated after 40, entry returned |
| f' | `EXHALE_PTC_DTAU0=0.01` -- a hundredfold stronger pseudo-transient shift | entry 4.931e+04 -> **6.272e+07** | the same refusal: the entry is still the best judged iterate at iteration 20, STAGNATED at 40, `done info = 2`, entry state returned |
| g | `Coupled carrier solve: True` -- the He/H element row as a Newton unknown instead of the alternation, with the same well-balanced option | the residual norm is 1.996 at entry and 1.664 to 1.977 at iterations 9 to 12, judged distance pinned at 3.31e+11 throughout | the base layer does not move at all (the element row does fall, 5.53e-3 to 1.53e-3); stopped at iteration 12. This is the coupled-route behaviour `input_read` already records as MEASURED on 2026-09-11, and the well-balanced option does not change it |

The pseudo-transient shift damps the continuity excursion by a factor 4.5,
and would have to damp it by a further 1300 to keep the judged distance
inside the entry value.

## 6. Recipe

**No configuration reaches a certified rung, so `run_closure.sh`,
`closure.json` and `make_models.py` are left exactly as they are.** Writing
a setting into them that does not certify would put a false recipe in the
tree.

What was tried, and why each fails:

- a smaller composition step (a) scales the insult linearly and the refusal
  does not depend on its size;
- more outer iterations (b) is not the binding limit, and the limit that is
  -- `n_stall_best` -- has no environment override;
- the trust-region restart (b'), the one configuration route that resets the
  stall counter, is unreachable here: `use_tr = (nspec_row .gt. 0)` and the
  partitioned route registers no species row;
- a stronger pseudo-transient shift (f') moves the continuity excursion by
  a factor 4.5, three decades short;
- the coupled route (g), which removes the alternation altogether, does not
  solve this planet's base layer at all, with the well-balanced option or
  without it.

What the campaign should do until the solver item is taken: the metals-free
`atomic_scalar` groups are unaffected and their cases certify as
`MODELS.md` section 6 reports. The 15 closure rungs and the
`atomic_scalarCNO` group carry this refusal, and the right record for them
is the one `MODELS.md` already has -- the route, the adapter and the seed
all work over the profile boundary; a converged rung waits on the solver.
The closure driver's behaviour is right as it stands: it stops at the first
iteration whose wind does not report `done info = 0`, and it should.

## 7. The minimal change, and where it belongs

Not made here: the item's rules confine it to `LHS1140b/models/*` and
`src/utils/element_flux_closure.py`, and the change is in
`src/modules/time_step/steady_newton.f90`.

The stagnation detector (`solve_steady_jfnk`, the block at the
`n_since_best .ge. n_stall_best` test) counts iterations since the best
JUDGED DISTANCE improved, and stops the solve at twice that count. That
rule is sound when a solve starts far from stationarity, where the judged
distance and the merit fall together. It cannot work when a solve starts
from a state whose continuity and momentum rows are already at their
floors and whose energy row alone is off, which is exactly what the
partitioned outer loop hands it at every pass after the first: the entry
state is then the best judged iterate by construction, because a Newton
iterate raises the continuity row by five to seven decades against a
3e-12 tolerance before it comes back down.

The smallest repair that keeps the detector's purpose is to reset
`n_since_best` when the MERIT improves on its own best, and to stop only
when neither functional has improved for the count. The evidence that this
is the distinction that matters is in the pass itself: over the 40
iterations the solve was stopped after, the judged distance never improved
and the merit fell from 1.50e-2 to 3.15e-3.

**The second half of the same repair is the forcing term.** The linear
tolerance of the partitioned route is the literal `1.0d-1` at both `pgmres`
calls of the `.not. use_tr` branch. Ten percent of the scaled residual left
over is ten percent of a norm the energy rows own, and the leftover lands
in the rows the certification measures against 3e-12. A forcing term that
tightens as the residual falls, which is what an inexact Newton method
normally carries, would make the continuity excursion second-order in the
step instead of first, and it is the quantity the whole refusal turns on.
It is a literal in one branch of one routine.

Two smaller observations that belong with it:

- the progress control of the outer loop responds to a joint distance that
  did not fall by halving `omega`, the element relaxation's step. In this
  rung the distance is held by the hydrodynamic energy row in every pass,
  and `omega` damps the one row that was converging. The response and the
  refusing row are not connected.
- `n_stall_best` (20), `outer_no_fall_max` (3) and the floor of `omega`
  (0.125) are `parameter`s. A measurement of any of them needs a rebuild,
  which is why (b) above could only be refuted by reading.

## 8. Noticed, not fixed

- **The neutral-carbon fraction.** Carbon is 94 to 99 percent neutral
  between 1.0 and 1.5 R_p at T = 1900 to 5900 K with the hydrogen ionized
  fraction at 1 to 4 percent, and its line cooling is then 79 to 94 percent
  of all the cooling there. C I has an ionization threshold of 11.26 eV,
  below the H I edge, so its photoionizing photons are not attenuated by
  hydrogen at all; a neutral fraction above that of hydrogen is worth a
  check against an independent calculation before the campaign's metal-bearing
  results are quoted. The sub-Lyman band is present in the field (section
  3.2), so this is a question about the rates, not about the grid. NOT
  CHECKED here.
- The failure is not confined to the closure rungs. Every metals-bearing
  case of the campaign that uses the element diffusion alternation carries
  it, `atomic_scalarCNO_gj1132_kzz1e9/HeH2.13` included.

## 9. Reproduction

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
D=$EX/LHS1140b/models/atomic_photochem_gj1132_kzzprofile/HeH9

# the scratch copy: the rung's own input, base state, profile and seed
mkdir -p /tmp/L4b/run/output
cp $D/k00/input.inp $D/k00/base.inp $D/k00/lower_atmosphere_profile.dat /tmp/L4b/run/
cp $D/k00/output/Hydro_ioniz_IC.txt $D/k00/output/Ion_species_IC.txt /tmp/L4b/run/output/
# the spectrum path in input.inp is relative to the case directory
sed -i "s|^Spectrum file: .*|Spectrum file: $EX/LHS1140b/sed/lhs1140_sed_gj1132_at_b.txt|" \
    /tmp/L4b/run/input.inp

# (c) the composition held at the seed's
cd /tmp/L4b/run && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 EXHALE_OUTER_PASSES=1 $EX/EXHALE.x

# (a) a quarter of the composition step
cd /tmp/L4b/run && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 EXHALE_DIFF_OMEGA=0.125 $EX/EXHALE.x

# the residual profile of the state the rung returned
cp $D/k00/output/Hydro_ioniz.txt /tmp/L4b/run/output/Hydro_ioniz_IC.txt
cp $D/k00/output/Ion_species.txt /tmp/L4b/run/output/Ion_species_IC.txt
cd /tmp/L4b/run && OMP_NUM_THREADS=2 EXHALE_RESIDUAL=1 $EX/EXHALE.x
#   -> output/residual_profile.txt

# the cooling channels of the same state: a post-processing pass
sed -i -e '/^Restart intent:/d' -e '/^Solver:/d' -e 's/^Do only PP: .*/Do only PP: True/' \
    /tmp/L4b/run/input.inp
cd /tmp/L4b/run && OMP_NUM_THREADS=2 $EX/EXHALE.x
#   -> output/Cooling_breakdown.txt, output/Heating_breakdown.txt
```

`docs/figures/lhs1140b_L4b_energy_row.pdf` is three panels against r/R_p,
built from the files those commands write: the scaled energy row from
`output/residual_profile.txt` with its 1e-6 tolerance drawn; cool/heat of
this rung and of `atomic_scalar_gj1132_kzz1e9/HeH1.60` from the two
`Hydro_ioniz.txt` files, with the C I share of the cooling from
`output/Cooling_breakdown.txt` beside them; and the relative change in the
helium mass fraction from the rung's `Ion_species_IC.txt` against its
`Ion_species.txt`. The band 1.15 to 1.50 R_p is shaded in all three.
