# D9fix: the molecular seed of a certified atomic state, and the x0.01 XUV refusals (2026-09-19)

Two defects reported by the D9 step 3 catalog refresh
(`docs/lhs1140b_catalog_refresh_20260919.md`). Every number is MEASURED in
this item unless labeled READ. Control binary: the tree's `EXHALE.x` of the
entry text, md5 `7670f31031fb4db91d27b44cb0da6f70` (copied to scratch before
any edit). Measured binary: a private build of the delivered text
(`make OBJDIR=build_gc EXE=EXHALE_gc.x`).

## 1. Verdicts

**Defect 1 is two code defects, both fixed.**

1. The ghost H2 closure of D5b-2 is a fixed-point substitution
   x -> g(x) = x2 (1 - x_ion(x)) whose slope on the LHS 1140 b ghosts is
   positive and close to one: 0.74 at He/H = 2.13 and 0.90 at He/H = 9.7
   (ratio of consecutive moves at pass 30). The sequence is a **monotone
   linear contraction**, not an oscillation or a cycle. 30 passes reach
   1.19e-6 and 2.9e-4 against 1e-10. Extrapolated from those ratios,
   substitution would need about 60 and 180 passes. The closure now takes a
   secant step on the scalar equation F(x) = g(x) - x and closes in 10 and 11
   passes, with the residual at 5.6e-16 and 2.2e-16.
2. Behind the closure stop was a second defect in the conversion
   (`molecular_seed_from_atomic_state.f90`, `transfer_h2`). load_IC does not
   read the lower ghost rows of the atomic pair. It fills them with the
   molecular run's reservoir, which already carries H2 (f_H2 = 5.29e-2
   against f_HI = 1.78e-6 on the HeH2.13 ghost). The transfer partitioned
   H I alone, capped the transfer at f_HI and overwrote f_H2, which destroyed
   the reservoir's H2 nuclei. The conversion's own check then refused the
   seed: the H nuclei changed by 0.99998 relative and the mass by 0.106. On
   the entry text this same defect stops every conversion of the
   `src/tests/molecular_seed` suite (H nuclei 0.985), so that suite gave
   5 PASS / 14 FAIL on the control.

Both molecular seeds now build. H nuclei are conserved to 4.0e-16 and 3.7e-16,
He nuclei to 0, and the equation-of-state closure is 3.4e-16 and 4.3e-16.

**Defect 2 is physical, and nothing was changed.** D2b alone moved the
three x0.01 XUV states: with `base_boundary.f90` and `viscous_conduction.f90`
reverted to their pre-D2b text, all three certify on the evaluate route.
The direction rule is not at fault. On all three, the wind window has
standing (d_window = 2.03e-4 to 2.05e-4, below 2e-3; M_wind = +1.8e-9 to
+2.4e-9, above 1e-12). It reads outflow, gives w_rev = 0 (the reservoir is
upwind), and the certification report states that this agrees with the base
face flux. What changed is the face density. Before D2b the continuous
handover took w_rev = 0.32 to 0.37, an interpolation between the reservoir's
and the interior's isentropes, which the D2b decision rejects as having no
physics. The parent states were solved on that blended face. Their interior
at the face sits 2.1 to 2.5 K warmer than the reservoir's isentrope, and the
acoustic matching reads a pressure mismatch of v_b - v_i = -7.7 to -13.3 cm/s.
Against winds of M_wind ~ 2e-9 (a face velocity of about 2e-4 cm/s), a state
with that mismatch is not stationary under the correct law. The old binary
certified a base that the accepted boundary law does not have.

## 2. Defect 1: the measurement

Reproduced in scratch through the runner's route (`run_case.sh`, the
molecular seed block, READ): the certified atomic generation copied as the
`_IC` pair into a source directory, the molecular case's `input.inp` and
`base.inp` (the SED path made absolute), and the binary run with
`EXHALE_MOLECULAR_SEED=<source>` and `EXHALE_MOLECULAR_SEED_X2=local`. The
control reproduces the catalog's stop exactly (ghost cell 0, 30 passes, last
move 1.18696e-6, residual 8.77853e-7). The run takes 0.3 s.

The sequence of the imposed partition, ghost cell 0, from a temporary print
(removed) in a private build:

| pass | HeH2.13 move | HeH9.7 move |
|---|---|---|
| 2 | 1.46e-1 | 3.48e-1 |
| 5 | 3.45e-3 | 8.61e-3 |
| 10 | 5.35e-4 | 2.98e-3 |
| 20 | 2.43e-5 | 8.25e-4 |
| 29 | 1.60e-6 | 3.20e-4 |
| 30 | 1.19e-6 | 2.89e-4 |

The move ratio settles to 0.74 and 0.90. Every move has the same sign
(x_H2 rises monotonically: 0.813, 0.959, 0.971, ... toward 0.98799). The fixed
point exists and is attractive; substitution is simply slow on it.

With the secant step, ghost cell 0 at HeH2.13 goes through the moves
1.46e-1, 1.36e-2, 1.02e-2, 4.19e-3, 1.14e-3, 1.08e-4, 2.29e-6, 4.25e-9,
1.60e-13, and the residual after the last pass is 5.6e-16. Cell -1 follows
the same course, and HeH9.7 closes at pass 11 with a residual of 2.2e-16.

**The tolerance keeps its value and meaning (1e-10 absolute in x_H2).** The
secant takes the residual of x_H2 - x2 (1 - x_ion) to 2e-16 to 6e-16, so the
ionization solve resolves the equation to rounding and 1e-10 asks nothing of
it that it cannot deliver. The closure test now reads that residual at the
returned state instead of the size of the last move. A failure is still an
error stop that prints the ghost cell, the passes, the last move and the
residual.

A first version guarded the secant point with the room check
`x + x_ion(previous pass) <= 1`, and it rejected every secant point until
pass 25. That check uses the ionization of the wrong partition. At the
secant point, the linear model of g gives x_ion = 1 - x/x2, which always
fits. The delivered guards are: both evaluations in the same radiation field
(past the self-field passes), F decreasing through the pair, and
0 <= x <= x2. Otherwise the pass substitutes.

## 3. Diff by file

| file | what changed |
|---|---|
| `src/modules/radiation/ionization_equilibrium.f90` | the ghost closure imposes the secant root of F(x) = x2 (1 - x_ion(x)) - x once two evaluations in the same field exist; the closure test is the residual at the returned state, at most `base_ghost_closure_tol`; the infeasible-handoff check reads the prescription x2 (1 - x_ion) and not the secant trial; six new private locals in the OpenMP private list; the tolerance comment restated with the measured resolution |
| `src/modules/init/molecular_seed_from_atomic_state.f90` | `transfer_h2` partitions the neutral hydrogen pool f_HI + 2 f_H2 (exactly f_HI in an atomic cell) and no longer zeroes H2+, H3+ and HeH+ (zero in every atomic cell, part of the reservoir in a ghost); unused `isp_H2p`, `isp_H3p`, `isp_HeHp` imports dropped |
| `src/tests/boundary_state/molecular_seed_of_an_atomic_state.sh` | new, two assertions (section 4) |
| `src/tests/boundary_state/run.sh`, `README.md` | the `molecular_seed` test and its statement |
| `docs/lhs1140b_stationary_D5b2_20260918.md` | the tolerance row of section 2.3 says what the criterion now is |

Nothing was changed for defect 2.

## 4. RED and GREEN

`src/tests/boundary_state/run.sh`, whole suite, on the measured build:
**37 of 37 PASS** (the 35 existing rows plus 2 new ones). The D5b-2 items
hold: the ghost rows of a restart are not an input (three rows the same), the
ghost H2 partition closure residual is 4.65605e-12 (the D5b-2 value, READ from
its memo section 6), the reservoir and the solved count are two quantities,
the face-cache reads are all current, the boundary model travels with the
state, and both evaluation sites derive the boundary after the sweep. The
22 contact-law rows all pass.

The new rows, on the control: **2 of 2 FAIL**.

```
FAIL ghost_h2_closure_of_the_molecular_seed measured=failed reference=closed tol=0
FAIL hydrogen_nuclei_of_the_molecular_seed measured=missing reference=at_most tol=1.0e-12
```

On the measured build: **2 of 2 PASS** (`closed`; `3.9996E-16`). The second
row is RED on its own as well. A private build with the closure fixed and the
transfer as it entered stops in the conversion with H nuclei changed by
9.9998e-1 (section 1).

Other suites whose fixtures run the changed code (a molecular ghost with a
stated handoff, or the seed conversion), on the measured build, with the
control run wherever mine did not pass everything:

| suite | measured build | control | reading |
|---|---|---|---|
| `molecular_seed` | 23 PASS / 2 FAIL | 5 PASS / 14 FAIL | the 12 rows the transfer defect broke now pass; the 2 left are section 6, item 1 |
| `run_mode` | 31 / 0 | | |
| `coupled_block_jacobian` | 5 / 2 | 5 / 2 | the same two rows fail on both builds |
| `coupled_source_step` | 30 / 0 | | |
| `species_face_flux` | 1 / 0 | | |
| `steady_species_rows` | 197 / 0 | | |
| `residual_determinism` | 5 / 1 | 5 / 1 | the driver reads `EXHALE_RESID_EXE`, so both runs used the control; atomic case, unreachable by this diff |
| `attempted_step` | 70 / 0 | | |
| `grid_and_gates` | 299 / 1 | `carrier_movement_bound_hold` 5 / 0 | section 6, item 2 |

## 5. Movement on the regression cases

Scratch copies of `backup/regression/mol_base_handoff` and
`backup/regression/mol_carrier`, `EXHALE_MAXSTEPS=12000` (their pinned step
count), `OMP_NUM_THREADS=1`, control against measured build. Data rows were
compared with the `#` lines dropped; the numbers are the largest relative
difference over every column.

| case | Hydro_ioniz.txt | Ion_species.txt | log10 Mdot | wall |
|---|---|---|---|---|
| `mol_base_handoff` | not byte-identical: 504 of 504 rows differ, max 1.39e-9 (row 83, column 6) | 504 of 504, max 2.57e-9 (row 83, column 37) | 10.57 on both | 1664 s / 1662 s |
| `mol_carrier` | 504 of 504, max 9.04e-12 (row 38, column 2) | 504 of 504, max 1.26e-12 (row 8, column 1) | 10.57 on both | 1662 s / 1655 s |

The branch that fires is the ghost closure. It now stops on the residual at
the returned state, and it steps by secant where it used to substitute, so
the ghost composition these snapshots carry differs at the 1e-10 level of the
closure tolerance. The whole column follows that difference through 12000
relaxation steps. Every difference is below the 1e-3 regression tolerance by
at least six decades. No golden was refreshed.

## 6. Defect 2: the evidence

Evaluate route (a scratch copy of `run_case.sh`'s evaluate block,
`OMP_NUM_THREADS=8`, lart4) on each case's certified parent generation g0001,
which the parent binary `EXHALE_3146d11b.x` wrote:

| case | parent binary 3146d11b | pre-D2b text (control with D2b reverted) | control 7670f310 |
|---|---|---|---|
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` | CERTIFIED | CERTIFIED | REFUSED |
| `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` | CERTIFIED | CERTIFIED | REFUSED |
| `atomic_photochem_gj1132x0.01_kzzprofile/HeH9.7` | CERTIFIED | CERTIFIED | REFUSED |

The pre-D2b build is the entry text with `src/modules/states/base_boundary.f90`
from the D5b-3 scratch tree (after D5b-2 and D5b-3, before D2b) and
`src/modules/time_step/viscous_conduction.f90` from the D4b tree (the D2b
memo section 4 lists every difference of that file as D2b's, READ). Since
reverting D2b alone restores all three, D5b-2 and D5b-3 were not reverted.

The base face, `EXHALE_RESIDUAL=1`:

| | HeH2.13 control | HeH2.13 pre-D2b | HeH9.7 control | photochem HeH9.7 control |
|---|---|---|---|---|
| v_i [cm/s] | -0.333 | -0.333 | -0.451 | -0.576 |
| M_wind | 2.426e-9 | 2.426e-9 | 2.150e-9 | 1.799e-9 |
| d_window | 2.048e-4 | 2.048e-4 | 2.043e-4 | 2.031e-4 |
| window standing | 1 | s = 1 | 1 | 1 |
| v_b [cm/s] | -8.04 | -8.04 | -10.89 | -13.90 |
| w_rev | 0 | 0.322 | 0 (pre-D2b 0.341) | 0 (pre-D2b 0.367) |
| T_res / T_i at the face [K] | 223.20 / 225.28 | same | 222.57 / 225.08 | 182.10 / 184.57 |
| rho_rev / rho_res | 0.9908 | same | 0.9889 | 0.9867 |
| base face flux / window mean | 38.4 | 0.99998 (the parent binary: the same) | 82.2 (pre-D2b 0.99998) | 177.8 (pre-D2b 0.99997) |
| certificate: direction | read from the window, agrees with the face flux | | same | same |

The brief's hypothesis, a silent window that leaves v_b's sign to decide, does
not hold: the window has standing with margins of 10 in d_window and 2e3 in
|M_wind|, and it decides outflow, which is right. The branch is right too. The
v_b the matching produces, -8 cm/s (inward), is a pressure mismatch of the
parent state against the reservoir, not a direction. The pre-D2b handover
read v_b/c against a width of 1e-8 and gave w_rev between 0 and 1. That blend
made the face flux match the wind, and it is the "average of the two
isentropes" that D2b section 3 removes as a law with no physics.

The re-solves of D9 step 3 did not certify within 30 minutes (hydro info=2
from pass 2). That is a solver and budget question on these very weak winds,
not a boundary defect, and it is left to the catalog work.

The three D2b catalog states (atomic HeH2.13, molecular HeH2.13 and HeH9.7)
are untouched: `base_boundary.f90` was not edited.

## 7. Noticed outside this item, reported and not fixed

1. `src/tests/molecular_seed/run.sh`: `molecular_seed_local_root_governs_somewhere`
   and `molecular_seed_local_crossover_is_recorded` fail on the measured
   build. `carrier_h2_chemical_root` returns 0 in every cell of this fixture
   because `bg_ready` is set only when `transported_rows_exist()`, and the
   hot-Uranus seed fixture transports no carrier. The same suite passes 25 of
   25 on `EXHALE_3146d11b.x`, where the root is 0.971 at cell 1. Neither line
   is in this item's diff, and the LHS 1140 b seeds (carrier transport on) do
   carry the root (287 and 259 cells). The regression lies in
   `diffusive_photochemistry.f90` or `ionization_equilibrium.f90`'s
   `bg_ready` statement between the two binaries and has not been localized.
2. `src/tests/grid_and_gates/carrier_movement_bound_hold.sh`,
   `carrier_bound_default_halves`: the row needs the control run of the
   fixture to halve its bound within 14 passes. On the entry text the halving
   at pass 14 is triggered by the mass row rising from 3.13e-11 to 7.16e-11,
   which is rounding. On the measured build the same passes print the same
   species rows to three digits and the hydrodynamic rows at 2e-11 to 8e-11,
   and pass 14 goes from 4.06e-11 to 4.03e-11, so no halving. The row's
   outcome is set by rounding noise; its pass count (or its fixture) should
   be chosen so that the halving follows from a stall of the species rows.
3. `coupled_block_jacobian` (two rows) and `residual_determinism` (one row)
   fail on the control as well.

## 8. How to repeat it

Scratch: `/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/59a81fec-9ef0-4992-9a57-6494f229e6bb/scratchpad/gc/`.
`evalgen.sh <case> <generation dir> <binary> <tag>` is the evaluate route;
`mol_HeH2.13/` and `mol_HeH9.7/` are the seed reproductions (see section 2
for the environment); `src/tests/boundary_state/run.sh molecular_seed` is the
row.
