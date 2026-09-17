# L6c: what moved the regression goldens of the consolidated tree

Item L6c of `docs/PLAN_20260913_lhs_stationary.md`. The consolidated private
binary `EXHALE_C.x` (md5 `f3f7fccdcbd3`) was run over the regression matrix
with `REGRESSION_EXE`, and the log `scratchpad/consol/check_C.log` reports
seven of the sixteen cases outside the 1e-3 relative tolerance. This item
attributes that movement to the source changes that entered the tree on
2026-09-13 and 2026-09-14, describes it physically case by case, and decides
which goldens are stale.

Every number below is MEASURED here unless it is marked READ. No golden was
refreshed, and no golden refresh was run.

---

## 1. Verdict

**All of the movement in all seven cases is the C/N/O line cooling change of
items L6 and L6b, and nothing else in today's tree contributes anything a
golden can see.** A control binary built from the consolidated tree with only
`Cool_coeff.f90` and `util_ion_eq.f90` reverted to HEAD reproduces the golden
of six of the seven cases BITWISE, including the marching step at which each
run stops, and reproduces the seventh, `wasp_full_newton`, to a maximum
relative difference of 2.2e-8. The difference between the measured and the
control binary is, file by file, row by row and value by value, the same
difference the log reports against the golden.

**The seven failing cases are exactly the seven cases that carry a C/N/O
elemental reservoir.** The nine passing cases carry helium only. The
correspondence is one to one over the whole matrix (table 4), and it is the
correspondence the change predicts: L6 and L6b touch the line cooling of six
ions, C I, C II, N I, N II, O I and O II, and no other coolant in the file.
Channel-resolved measurement confirms it directly: with the gas state held,
the column-integrated cooling of those six ions changes by factors 0.001 to
0.5 while Mg I, Mg II, Ca II, Na I and Fe II change by at most 1.5 percent
(table 5).

**The changed cooling is the cause and the changed profiles are the effect,
not the other way round.** With one equilibrium solve from a single held
state, the two builds differ in temperature by at most 3.7e-3 while the
cooling differs by factors 0.02 to 0.78 (table 6). The size of each case's
golden movement is ordered by how much of the local heating that cooling
change is worth (table 7): `lower_profile` 51 percent of the heating at the
median cell, the three WASP-121 b cases 17 to 30 percent, the two molecular
gates 0.3 percent, `oxygen_chemistry` under 0.01 percent.

**Judgement on the goldens: the seven are stale and should be refreshed once,
at the end of the change series, and the nine others should not be touched.**
The cause is a physics correction, not a numerical accident: the coronal,
low-density-limit evaluation of six C/N/O coolants was being applied one to
three decades above the critical densities of the metastable levels that
carry 99.7 percent of their emission, which overran the LTE ceiling on what
those ions can radiate (READ, `docs/lhs1140b_stationary_L6_20260913.md`
section 1 and `..._L6b_20260913.md` section 1). Physical correctness is the
acceptance criterion, so the references built on the old formulation are the
things to renew. The command is written in section 8 and was NOT run.

**One caveat on what a refreshed golden will pin.** Only `wasp_full_newton`
converges: it reaches `info = 0` in both builds, so its movement is a real
change of the steady solution, log10 Mdot 13.30 to 13.33 and a 5.06 percent
rise in the mass flux. The other six stop on a `du` threshold or at a fixed
step count and are relaxation snapshots, so their new goldens will pin a
transient exactly as their old ones do.

---

## 2. The two builds, and what separates them

Both trees are isolated copies, not the working tree, because other workers
were editing files in it throughout this item. The snapshot was taken once, at
01:38 on 2026-09-14, and both trees were made from that one snapshot, so the
file that item L4g is still editing, `steady_newton.f90`, is the SAME copy in
both. (That file's modification time had moved again by the end of this item,
to 05:51, which the fixed snapshot is deliberately blind to.)

| tree | contents |
|---|---|
| control | HEAD `43bc28cef587` plus every working-tree change of the snapshot, EXCEPT `src/modules/radiation/Cool_coeff.f90` and `src/modules/radiation/util_ion_eq.f90`, which are HEAD's |
| measured | HEAD `43bc28cef587` plus every working-tree change of the snapshot |

MEASURED: `diff -rq` over the two `src/` trees reports exactly two differing
files, the two named above. Both were built with the conda-forge gfortran
16.2.0 on the default PATH, the compiler the goldens were taken with; a first
pair built with `/usr/bin/gfortran` 13.1.0 was discarded, because every object
differed in size from the delivered `build_C/` objects and `EXHALE_C.x` cannot
have come from that compiler.

### The measured binary is not md5-equal to `EXHALE_C.x`, and why

| binary | md5 (first 12) | size |
|---|---|---|
| measured (this item) | `b261c3287e64` | 3622336 |
| control (this item) | `f1c5dbffcd8e` | 3622248 |
| `EXHALE_C.x` (delivered) | `f3f7fccdcbd3` | 3622248 |

MEASURED, object by object against the delivered `build_C/`: of the 122
objects, exactly ONE differs, `steady_newton.o` (989704 bytes against 988464);
the other 121 are byte-identical. `src/modules/time_step/steady_newton.f90`
carried the modification time 2026-09-14 01:37, after the 23:12 build that
produced `EXHALE_C.x`, so the one difference is item L4g's edit landing
between that build and this item's snapshot. Nothing else in the tree moved in
between.

**The consequence is bounded by measurement, not argued.** `wasp_full_newton`
is the only case in the matrix that reaches the steady solver at all; there
the measured binary and `EXHALE_C.x` agree to 7.6e-9 (`Hydro_ioniz.txt`) and
1.3e-8 (`Ion_species.txt`). On all six other cases the two binaries are
BITWISE identical in all four compared files (table 3, last column). So the
measured build is the consolidated build for every purpose this item serves.

---

## 3. Attribution

Each case was run on a scratch copy with each binary, single threaded
(`OMP_NUM_THREADS=1`) and with the case's own `maxsteps` exported as
`EXHALE_MAXSTEPS`, which is what `backup/regression/run_check.sh` does. The
comparison is `backup/regression/compare_within_tolerance.py` at 1e-3, the
same one the harness uses.

### Table 1. Where each run stops

| case | golden | control | measured |
|---|---|---|---|
| `wasp_full` | count 17635, du 9.9551e-04 | **17635, 9.9551e-04** | 12251, 9.7613e-04 |
| `wasp_he23off` | count 17623, du 9.8323e-04 | **17623, 9.8323e-04** | 12236, 9.7535e-04 |
| `wasp_full_newton` | count 4700, du 9.9386e-03 | **4700, 9.9386e-03** | 4709, 9.9192e-03 |
| `mol_metals` | count 12000, du 7.9064e-01 | **12000, 7.9064e-01** | 12000, 7.9165e-01 |
| `mol_ir_bands` | count 12000, du 7.9062e-01 | **12000, 7.9062e-01** | 12000, 7.9161e-01 |
| `lower_profile` | count 12000, du 4.1694e+00 | **12000, 4.1694e+00** | 12000, 3.1666e+00 |
| `oxygen_chemistry` | count 1000, du 4.5171e+00 | **1000, 4.5171e+00** | 1000, 4.4960e+00 |

The control reproduces the golden's stopping point to every printed digit in
all seven cases, including the two that stop on a `du` threshold after more
than seventeen thousand steps.

### Table 2. The attribution itself

"control vs golden" asks whether reverting the two cooling files restores the
reference. "measured vs control" asks what those two files alone produce.

| case | file | control vs golden | measured vs control |
|---|---|---|---|
| `wasp_full` | `Hydro_ioniz.txt` | IDENTICAL | 1.976e+00 at row 209 col 3 |
| | `Ion_species.txt` | IDENTICAL | 8.168e-01 at row 299 col 26 |
| | `Hydro_ioniz_adv.txt` | IDENTICAL | 1.976e+00 at row 209 col 3 |
| | `Ion_species_adv.txt` | IDENTICAL | 1.000e+00 at row 290 col 35 |
| `wasp_he23off` | `Hydro_ioniz.txt` | IDENTICAL | 1.992e+00 at row 211 col 3 |
| | `Ion_species.txt` | IDENTICAL | 8.188e-01 at row 299 col 26 |
| | `Hydro_ioniz_adv.txt` | IDENTICAL | 1.992e+00 at row 211 col 3 |
| | `Ion_species_adv.txt` | IDENTICAL | 1.000e+00 at row 292 col 35 |
| `wasp_full_newton` | `Hydro_ioniz.txt` | 1.150e-08 | 6.532e-01 at row 1 col 7 |
| | `Ion_species.txt` | 2.171e-08 | 6.187e-01 at row 482 col 23 |
| | `Hydro_ioniz_adv.txt` | 1.000e+00 at row 201 col 10 (see section 7) | 1.000e+00 at row 201 col 10 |
| | `Ion_species_adv.txt` | 2.171e-08 | 6.187e-01 at row 482 col 23 |
| `mol_metals` | `Hydro_ioniz.txt` | 3.562e-08 | 7.315e-01 at row 32 col 7 |
| | `Ion_species.txt` | 5.148e-09 | 2.609e-02 at row 491 col 7 |
| | `Hydro_ioniz_adv.txt` | 2.374e-06 | 7.561e-01 at row 33 col 7 |
| | `Ion_species_adv.txt` | 2.114e-08 | 2.609e-02 at row 491 col 7 |
| `mol_ir_bands` | `Hydro_ioniz.txt` | 3.125e-08 | 3.381e-01 at row 32 col 7 |
| | `Ion_species.txt` | 3.351e-09 | 2.605e-02 at row 491 col 7 |
| | `Hydro_ioniz_adv.txt` | 2.470e-06 | 9.117e-01 at row 436 col 10 |
| | `Ion_species_adv.txt` | 3.599e-08 | 2.605e-02 at row 491 col 7 |
| `lower_profile` | all four | IDENTICAL | 1.269e+00 at row 474 col 3, and the three others |
| `oxygen_chemistry` | all four | IDENTICAL | 1.947e+00 at row 15 col 6, and the three others |

In every case the "measured vs control" column names the same file, row,
column and pair of values as the "measured vs golden" comparison, which is
the line `check_C.log` reports. The residual of the control against the
golden, where it is not zero, is between 3.4e-9 and 2.5e-6: that is the
combined footprint of every other change in today's tree, three to six
decades below the acceptance tolerance.

### Table 3. The measured build against the delivered `EXHALE_C.x`

| case | `Hydro_ioniz.txt` | `Ion_species.txt` | `Hydro_ioniz_adv.txt` | `Ion_species_adv.txt` |
|---|---|---|---|---|
| `wasp_full` | IDENTICAL | IDENTICAL | IDENTICAL | IDENTICAL |
| `wasp_he23off` | IDENTICAL | IDENTICAL | IDENTICAL | IDENTICAL |
| `wasp_full_newton` | 7.604e-09 | 1.277e-08 | zero against 1.4e-14 (section 7) | 1.277e-08 |
| `mol_metals` | IDENTICAL | IDENTICAL | IDENTICAL | IDENTICAL |
| `mol_ir_bands` | IDENTICAL | IDENTICAL | IDENTICAL | IDENTICAL |
| `lower_profile` | IDENTICAL | IDENTICAL | IDENTICAL | IDENTICAL |
| `oxygen_chemistry` | IDENTICAL | IDENTICAL | IDENTICAL | IDENTICAL |

### What this leaves for the other items of the day

The changes that were in the tree beside L6 and L6b, and what the control
measurement says about each:

| item | files | what a golden sees |
|---|---|---|
| L7b | `diffusive_photochemistry.f90`, `certification.f90`, `energy_semi_implicit.f90`, `ionization_equilibrium.f90`, `EXHALE_main.f90` | the molecular gates move by 3.4e-9 to 2.5e-6; `lower_profile` and `oxygen_chemistry` do not move at all |
| L4d, L4e, L4g | `steady_newton.f90` | `wasp_full_newton` reaches the same stopping step and the same state to 2.2e-8, from `info = 0` at both `||R||` 6.2e-9 (control) and 3.1e-9 (measured) |
| L10 | `post_process_adv.f90` | the two `*_adv.txt` files are bitwise identical to the golden in six of the seven cases |
| L12-A | the `System_*` family, `constrained_chemical_equilibrium.f90`, `ion_residual_core.f90` | nothing separable from the above |
| L7 seed mode | `load_IC.f90`, `write_output.f90`, `molecular_seed_from_atomic_state.f90` | nothing; every branch is behind `molecular_seed_on()`, which no regression case sets |

### Table 4. Why these seven and not the others

| case | elemental reservoir written in the golden header (READ) | verdict in `check_C.log` |
|---|---|---|
| `wasp_full` | He, C, O, N, Mg, Ca, Na, Fe | FAIL |
| `wasp_he23off` | He, C, O, N, Mg, Ca, Na, Fe | FAIL |
| `wasp_full_newton` | He, C, O, N, Mg, Ca, Na, Fe | FAIL |
| `mol_metals` | He, C, O, N, Mg, Ca, Na, Fe | FAIL |
| `mol_ir_bands` | He, C, O, N, Mg, Ca, Na, Fe | FAIL |
| `lower_profile` | He, C, O, N (carried by the profile, no `metals.inp`) | FAIL |
| `oxygen_chemistry` | He, C, O, N | FAIL |
| `mol_base_handoff` | He only | PASS |
| `mol_lyman_werner` | He only | PASS |
| `mol_diffusion` | He only | PASS |
| `mol_sec_ion` | He only | PASS |
| `mol_carrier` | He only | PASS |
| `hydrostatic_column` | He only | PASS |
| `hp_zero_seed` | He only | PASS |
| `hp_trace_seed` | He only | PASS |
| `hp_front` | He only | PASS |

---

## 4. What changed in the cooling, measured at a held state

`Do only PP: True` performs one equilibrium solve at a loaded state and writes
the heating and cooling of the composition it returns. Running it with both
binaries from the SAME state, the control run's own final state, separates the
change in the cooling coefficient from the change in the gas the coefficient
is evaluated on.

### Table 5. Column-integrated cooling by channel, control to measured

Only channels carrying more than 0.1 percent of the control total are listed.

| case | channel | control | measured | ratio | share of the control total |
|---|---|---|---|---|---|
| `wasp_full` | total | 9.8474e-03 | 6.3734e-03 | 0.647 | |
| | Mg II | 1.3501e-03 | 1.3642e-03 | 1.010 | 13.7 % |
| | O I | 1.1326e-03 | 3.7809e-06 | 0.003 | 11.5 % |
| | O II | 1.1292e-03 | 6.1967e-07 | 0.001 | 11.5 % |
| | Fe II | 7.9578e-04 | 8.0542e-04 | 1.012 | 8.1 % |
| | N II | 5.0471e-04 | 4.6588e-06 | 0.009 | 5.1 % |
| | C I | 4.3366e-04 | 3.9073e-06 | 0.009 | 4.4 % |
| | Ca II | 3.3008e-04 | 3.3260e-04 | 1.008 | 3.4 % |
| | Na I | 3.0572e-04 | 3.0573e-04 | 1.000 | 3.1 % |
| | C II | 2.8053e-04 | 1.3846e-04 | 0.494 | 2.8 % |
| | N I | 2.4431e-04 | 2.6227e-07 | 0.001 | 2.5 % |
| `lower_profile` | total | 1.4318e-05 | 1.0512e-06 | 0.073 | |
| | C I | 1.2233e-05 | 2.7798e-08 | 0.002 | 85.4 % |
| | O I | 7.1735e-07 | 1.0297e-07 | 0.144 | 5.0 % |
| | N II | 2.3164e-07 | 7.0218e-10 | 0.003 | 1.6 % |
| | N I | 1.5063e-07 | 8.8148e-10 | 0.006 | 1.1 % |
| `mol_metals` | total | 2.5232e-06 | 2.0707e-06 | 0.821 | |
| | C I | 4.5185e-07 | 3.7294e-10 | 0.001 | 17.9 % |
| | O I | 3.0567e-07 | 3.0515e-07 | 0.998 | 12.1 % |
| | Fe II | 2.3218e-07 | 2.3216e-07 | 1.000 | 9.2 % |

The six ions L6 and L6b rewrote are the only ones that move. Mg I, Mg II,
Ca II, Na I and Fe II stay within 1.5 percent, and what moves them at all is
the electron density of the re-solved composition, not their own coefficient.
The one entry above that is not a fall, O I in `mol_metals` at ratio 0.998, is
the same statement: at 900 to 2600 K that ion radiates through its ground-term
fine structure, which the change does not touch.

### Table 6. The cooling ratio and the temperature it was evaluated at

| case | largest relative change in T over the one solve | cooling ratio measured/control, over cells above 1 percent of the peak | column-integrated ratio |
|---|---|---|---|
| `wasp_full` | 3.58e-03 | 0.350 to 0.779 | 0.647 |
| `wasp_he23off` | 3.66e-03 | 0.350 to 0.777 | 0.642 |
| `wasp_full_newton` | 3.71e-03 | 0.349 to 0.779 | 0.663 |
| `mol_metals` | 1.49e-01 in exactly 2 of 504 rows, 1.1e-04 elsewhere | 0.702 to 1.257 | 0.821 |
| `mol_ir_bands` | 1.07e-06 | 0.421 to 0.999 | 0.845 |
| `lower_profile` | 1.05e-04 | 0.020 to 0.465 | 0.073 |
| `oxygen_chemistry` | 1.88e-06 | 0.090 to 1.712 | 1.001 |

The `mol_metals` entry is the one that needs a word: the two rows that move
are 198 and 199, at r = 1.0871 and 1.0883, which simply exchange values
(875.1 and 1033.8 K become 1028.8 and 879.8 K). That is the two-cell
alternation this gate is known to carry at the H2 front, re-picked, not a
temperature change. All 502 other rows hold to better than 1.1e-04.

### Table 7. The size of the perturbation, as a fraction of the local heating

`|cool(measured) - cool(control)| / heat`, from the same one-solve pair, over
cells carrying more than 1 percent of the peak cooling.

| case | median | maximum | cool/heat of the control over those cells |
|---|---|---|---|
| `lower_profile` | 0.514 | 3.330 at r = 1.0031 | 0.155 to 3.396 |
| `wasp_he23off` | 0.304 | 0.410 at r = 1.4241 | 0.194 to 0.909 |
| `wasp_full` | 0.303 | 0.398 at r = 1.0334 | 0.195 to 0.909 |
| `wasp_full_newton` | 0.174 | 0.393 at r = 1.3928 | 0.177 to 0.910 |
| `mol_ir_bands` | 0.0031 | 0.071 at r = 1.0827 | 0.010 to 0.262 |
| `mol_metals` | 0.0029 | 0.056 at r = 1.0059 | 0.000 to 0.189 |
| `oxygen_chemistry` | 0.0000 | 0.002 at r = 1.0061 | 0.000 to 0.919 |

This column is what orders the case descriptions below. Where the C/N/O lines
were a large part of the energy budget the solution moves a great deal; where
another coolant carries the layer the solution barely moves, and what the
golden then reports is the position of a transient.

---

## 5. Case by case

### 5.1 `wasp_full` and `wasp_he23off`, WASP-121 b, marching to a `du` gate

Both stop when the fractional radial spread of the mass flux falls below
1e-3, and the corrected cooling moves that stopping step from 17635 to 12251
(and 17623 to 12236). The two states compared are therefore the same wind at
two different moments of its approach to the gate, and the regression notes
already put the mass-loss spread of a `du`-threshold stop at several percent
(READ, `CLAUDE.md`). log10 Mdot is 13.35 against 13.36 for both.

What the cooling does to the structure, control against measured at the
stopping state:

| r/Rp | cool, control | cool, measured | T, control | T, measured | rho, control | rho, measured |
|---|---|---|---|---|---|---|
| 1.0169 | 5.0130e-06 | 1.0313e-06 (-79 %) | 3072 K | 2743 K | | |
| 1.0500 | 3.1760e-06 | 7.7368e-07 (-76 %) | | | | |
| 1.1147 | | | | | 1.0122e+11 | 7.4224e+10 (-27 %) |
| 1.1230 | | | 5255 K | 6326 K (+20 %) | | |
| 1.2000 | 3.8211e-05 | 4.3995e-05 (+15 %) | 8240 K | 8805 K (+6.9 %) | | |
| 1.5006 | | | 10871 K | 12411 K (+14.2 %) | | |

The chain is: less cooling below 1.1 R_p, so the base runs cooler in the
lowest cells but loses less of the deposited energy, so the wind arrives
hotter and denser higher up, where the surviving metal channels then radiate
20 percent more.

The entry the log names for both cases, `Hydro_ioniz.txt` row 209 (211)
column 3, is the VELOCITY, and it is at r = 1.064 (1.065), inside the base
region where this configuration breathes. MEASURED over rows 196 to 224: the
golden state has v between +2730 and +6694 cm/s there while the measured
state has v between -5770 and -726 cm/s. Both are far below the local sound
speed and both are the base oscillation caught at a different phase, which is
also why the relative measure reports a number near 2: it divides a difference
of two comparable numbers of opposite sign by the larger. It is not a reversal
of the wind. The two states are 5384 marching steps apart.

`Ion_species.txt` row 299 column 26 is Na I at r = 1.133, 906.9 against
4950.1. Sodium has the lowest ionization potential of the trace metals in this
run, so its neutral fraction is the most sensitive quantity in the file to the
temperature rise of the same cell, and MEASURED it is the largest relative
move of any species column here. The whole metal block moves with it: at their
own worst rows, Ca I -68 %, Fe III +174 %, He III +146 %, O III +119 %,
Mg I -35 %, C I -25 %.

`Ion_species_adv.txt` column 35 is `adv_T_status`, a validity flag of the
advection-corrected profile, 1 against 0 at one row: the post-process refuses
or accepts the correction at a different cell because the profile moved.

### 5.2 `wasp_full_newton`, the only converged case

This case marches to a loose gate and then runs the steady solver, and both
builds reach `info = 0`: the control at `||R||` 6.195e-09 with a flux spread
of 2.172e-13, the measured at 3.110e-09 and 1.075e-12. So its movement is NOT
path dependence. It is the change of the steady solution itself.

| r/Rp | T control | T measured | cool control | cool measured | rho control | rho measured |
|---|---|---|---|---|---|---|
| 1.0049 | 2523.9 | 2563.0 (+1.5 %) | 3.0852e-06 | 1.2355e-06 (-60 %) | 2.9147e+12 | 2.8749e+12 |
| 1.0199 | 2815.6 | 2938.7 (+4.4 %) | 2.6606e-06 | 1.3110e-06 (-51 %) | 1.4880e+12 | 1.4525e+12 |
| 1.1004 | 4796.6 | 4776.9 | 5.7873e-06 | 2.4845e-06 (-57 %) | 1.0293e+11 | 1.1434e+11 (+11 %) |
| 1.2000 | 8825.9 | 8947.6 (+1.4 %) | 6.5022e-05 | 5.6505e-05 | 1.9465e+10 | 2.2228e+10 (+14 %) |
| 1.4011 | 10607.4 | 11900.2 (+12.2 %) | 2.0373e-05 | 1.8308e-05 | | |
| 1.5589 | 11312.9 | 12768.5 (+12.9 %) | 6.8066e-06 | 5.6167e-06 | | |

Median `rho v r^2` over r > 1.3 R_p: 7.68156e+15 against 8.07035e+15,
**+5.06 percent**. log10 Mdot 13.30 against 13.33.

The log names `Hydro_ioniz.txt` row 1 column 7, which is the cooling of the
INNER GHOST cell at r = 0.9998, 1.03e-06 against 2.98e-06. That is the same
minus-65-percent base cooling as row 11 of `wasp_full`; it is the largest
RELATIVE move in the file only because the ghost sits deepest, where the
fraction of the cooling carried by C I and O I is largest.
`Ion_species.txt` row 482 column 23 is Ca I at r = 1.494, 0.1476 against
0.3870: neutral calcium at 12500 K instead of 11000 K, a factor 2.6 fall, the
strongest single response in the file to the outer temperature rise.

### 5.3 `lower_profile`, the largest movement, and why

Here the C/N/O lines are not a correction to the energy budget, they are most
of it. MEASURED at the held state: C I alone carries 85.4 percent of the
column-integrated cooling of the control, and it falls by a factor 440. The
resulting change in the NET heating is 51 percent of the heating at the median
cell and 3.3 times the heating at r = 1.0031 (table 7), because there the
control's `cool/heat` stands at 3.4, that is, the cell was being cooled harder
than it was heated.

The case is a fixed 12000-step snapshot of a wind that has not arrived: `du`
is 4.17 in the control and 3.17 in the measured, so the flux spread is
several hundred percent and the outer half of the grid is still the initial
column being overrun by a heated front. With the base no longer over-cooled,
that front runs further in the same 12000 steps.

MEASURED, the outermost cell still above 500 K at step 12000:

| | row | r/Rp | T |
|---|---|---|---|
| control (bitwise the golden) | 472 | 3.1590 | 764 K |
| measured | 482 | 3.4247 | 511 K |

Ten cells, 0.27 R_p. Everything the log reports for this case is that
displacement read off row by row: `Hydro_ioniz.txt` row 474 column 3 is the
velocity at r = 3.2097, +518501 cm/s in the measured state where the front has
passed against -139261 cm/s in the golden where it has not; row 475 column 11
is O I, 1.160 against 0.00675, the same cell before and after the front. The
comparison aligns the two states by row, so a front at two positions reports
as an order-unity difference at every row between them.

### 5.4 `mol_metals` and `mol_ir_bands`, the hot-Uranus gate with metals

These two move the least of the seven, and the reason is physical. The layer
is at 900 to 2600 K, where the C/N/O channels are down in the exponential
tail; what cools it is the H3+ infrared emission and the ground-term fine
structure, neither of which this change touches. MEASURED at the held state,
the change in the net heating is 0.3 percent of the heating at the median cell
and at most 7 percent anywhere (table 7), against `cool/heat` of 0.19 at its
largest.

At the 12000-step snapshot the interior is accordingly unmoved:

| r/Rp | cool control | cool measured | ratio | T control | T measured |
|---|---|---|---|---|---|
| 1.0061 | 1.7536e-08 | 4.7276e-09 | 0.270 | 1247.9 | 1248.0 (+0.01 %) |
| 1.0201 | 1.0104e-08 | 3.2794e-09 | 0.325 | 1227.5 | 1227.6 (+0.01 %) |
| 1.0500 | 2.7708e-09 | 1.6605e-09 | 0.599 | 1170.9 | 1171.0 |
| 1.2005 | 5.7386e-10 | 5.6285e-10 | 0.981 | 1147.0 | 1146.9 |
| 3.0053 | 1.3631e-10 | 1.3600e-10 | 0.998 | 2026.0 | 2026.7 |
| 4.1718 | 1.9972e-11 | 2.0680e-11 | 1.035 | 1566.8 | 1585.2 (+1.17 %) |
| 4.2103 | 1.2928e-11 | 1.3367e-11 | 1.034 | 1357.8 | 1373.9 (+1.19 %) |

The cooling of the shielded base falls by 40 to 73 percent and the temperature
there does not notice, because that cooling was 8 percent of the heating. The
only place anything exceeds 1 percent is the outermost two cells, r = 4.17 and
4.21, where the wind front of this snapshot is still advancing and the same
small absolute change is the largest relative one. That is where the log's
`Ion_species.txt` row 491 column 7 sits: He 2^3S at r = 4.1718, 544.5 against
559.1, a 2.6 percent move, and every other species column at that same row
moves by the same 1 to 2 percent.

The log's other entry, `Hydro_ioniz.txt` row 32 column 7, is the cooling at
r = 1.0059 itself: 4.75e-09 against 1.77e-08 for `mol_metals` and 2.53e-08
against 3.82e-08 for `mol_ir_bands`. That is the change in the coefficient,
reported directly, on a quantity that is 8 percent of the local heating.

`mol_ir_bands` additionally reports `Hydro_ioniz_adv.txt` row 436 column 10,
3.93e-05 against 4.45e-04. That column is `adv_mass_row`, the discrete mass
row the post-process writes beside its two validity flags. Unlike the
`wasp_full_newton` entry of section 7 both values here are far from zero, so
this one is a real factor 11 in a residual measure at r = 2.633, in the outer
region the case has not relaxed; it is not the exact-zero artefact.

### 5.5 `oxygen_chemistry`, and the row 173 item

This case receives the SMALLEST perturbation of the seven: the change in the
net heating is below 0.01 percent of the heating at the median cell and 0.2
percent at its largest, confined to a thin shell between about 1.005 and 1.012
R_p. MEASURED at the held state, the twelve cells carrying more than 10
percent of the peak cooling have ratio 1.000 to four figures, and the column
integral is 1.001. The layer is cooled by H3+, H2, H2O and OH.

It nevertheless produces order-unity entries in the log, and the reason is the
state, not the change. The case is a 1000-step cap on a relaxation that is
nowhere near stationary: `du` is 4.50, and the base already carries a
cell-to-cell alternation in both states. MEASURED: the GOLDEN state itself has
34 cells with a negative `heat` column and five isolated temperature spikes
more than 15 percent above both neighbours; the measured state has 45 and
five. A negative entry in that column is ordinary: `heat` is the sum of every
heating channel including the signed collisional oxygen chemistry
(`molecular_reaction_heat::oxygen_chemical_heating`, whose reverse rows O1r
and O2r carry negative reaction energies), so a cell in which the reverse
channels outrun the photon deposition is endothermic on balance. The log's
`Hydro_ioniz.txt` row 15 column 6, -3.24e-03 against +3.42e-03, is one such
cell at r = 1.0025 changing membership of that set.

**The row 173 column 37 item is H3+, and it is one transient cell that moved
nine cells outward.** MEASURED, with the states printed at their own
positions:

| state | row | r/Rp | T | H2 | H3+ | H2O |
|---|---|---|---|---|---|---|
| golden (= control) | 163 | 1.0537 | 1596.0 | 4.915e+08 | 3.636e-01 | 7.227e+03 |
| golden (= control) | **164** | **1.0544** | **1862.1** | **1.082e+11** | **4.779e+04** | **4.841e+07** |
| golden (= control) | 165 | 1.0551 | 1248.1 | 2.438e+08 | 7.787e-02 | 2.964e+03 |
| measured | 172 | 1.0603 | 1575.0 | 1.444e+08 | 6.378e-02 | 8.551e+02 |
| measured | **173** | **1.0611** | **1835.2** | **7.489e+10** | **3.085e+04** | **3.351e+07** |
| measured | 174 | 1.0619 | 1242.6 | 7.049e+07 | 1.579e-02 | 3.347e+02 |

The excursion is a single cell whose H2 stands more than two decades above
both of its neighbours, and BOTH states carry exactly one of them: the golden
at r = 1.0544 (a factor 220 above its neighbours), the measured at r = 1.0611
(a factor 520). At the peak the two agree to within a factor 1.4 in H2, H3+
and H2O and to 1.5 percent in temperature. What differs is the position, nine
cells. The reported pair 30845.7 against 0.019181 is the measured peak read
against the golden's off-peak value at the same row; the mirror image holds at
the golden's own row, 164, where the values are 47787.3 (golden) against
0.112352 (measured).

So the item is not a species changing by six decades. It is a single-cell
excursion of a violently transient snapshot standing nine cells further out,
which a row-aligned comparison must report at full magnitude. The same reading
applies to `Hydro_ioniz_adv.txt` row 431 column 3 (the velocity at r = 2.341,
+293.6 against -1428.4, both small) and to `Ion_species_adv.txt` row 3 column
42, which is `adv_T_status`, a flag.

**This deserves recording beyond this item.** A perturbation worth at most 0.2
percent of the local heating, in a handful of cells near the base, displaces
this snapshot's structure by nine cells at step 1000. That golden is not a
stable reference against any change that touches the base energy budget at
all, and it will keep failing for reasons that are not physics defects. Whether
the case should be capped where its solution is less violently transient is a
question for the plan, not for this item.

---

## 6. What was NOT measured

- Whether any of the seven refreshed states would be a better or worse match
  to a published solution. This item compares the code against itself.
- The nine passing cases were not rerun here. The consolidated log reports
  them PASS at 1.1e-07 to 1.2e-04, and none of them carries a C/N/O
  reservoir, so the change this item attributes cannot reach them.
- The correctness of the L6 and L6b cooling itself. That is those items' own
  measurement (READ, their sections 1 and 6); this item takes it as given and
  asks only what it moved.
- Any case outside `DEFAULT_CASES`.

---

## 7. Noticed outside the scope of this item, and not fixed

**`Hydro_ioniz_adv.txt` column 10, `adv_mass_row`, cannot be compared by a
relative measure.** It is the discrete mass-row residual the post-process
writes, and MEASURED on `wasp_full_newton` it is exactly zero in 1 of 504 rows
of the golden and in 12 rows of the control, with a median nonzero magnitude
near 1e-14 and a maximum of 7.6e-03. `compare_within_tolerance.py` divides by
`max(|a|,|b|)`, so an exact zero against a 1.7e-14 roundoff reports a relative
difference of 1.000 and fails the file. This is why table 2 shows
`wasp_full_newton` failing `Hydro_ioniz_adv.txt` in the control column while
its two physics files agree to 2.2e-08; the largest ABSOLUTE difference in
that column between control and golden is 2.6e-12. Any byte-identity claim on
that file is hostage to it. An absolute floor in the comparison for that
column would settle it. Not fixed here: `compare_within_tolerance.py` is the
shared acceptance interface and other workers were running it during this
item.

**The regression case directories that use "arm" as a noun** (`armA_LW`,
`armD_D2`, `armHeH_*` and the rest) are already recorded by
`docs/lhs1140b_stationary_L6_20260913.md` section 8 and
`..._L6b_20260913.md` section 9. Renaming them moves stored data, so they wait
for the user. Recorded again only to say the earlier entries still stand.
**DONE 2026-09-16** (PLAN_20260916_rev3 section 10): all sixteen were
renamed; the old-to-new mapping is `docs/named_case_audit.md` section 6.

---

## 8. Judgement on the goldens, and the command

**Refresh these seven, once, at the end of the change series:**

```
wasp_full  wasp_he23off  wasp_full_newton  mol_metals  mol_ir_bands
lower_profile  oxygen_chemistry
```

**Do not touch the other nine.** They pass at 1.1e-07 to 1.2e-04 and carry no
C/N/O reservoir.

The grounds: the movement is entirely a correction to the C/N/O line cooling,
which was being evaluated in the coronal limit above the critical densities of
the levels carrying its emission and above the LTE ceiling those levels can
radiate to. The old references encode the superseded formulation. The standing
rule that a sub-0.1-percent movement never spends a golden refresh does not
apply: the movement here is 2.6e-02 to 2.0e+00 in six cases, and in the one
converged case it is a 5.06 percent change of the mass flux and 12 to 13
percent in the outer temperature.

The procedure the log of record uses (READ,
`docs/Update_EXHALE_stage2.md`, the 2026-09-08 19:05 entry) is: archive the
present set first, verify the archive against `golden/`, then `check`, then
`golden`, then `check` again, because `golden` mode snapshots whatever sits in
the case directory without re-running.

```bash
cd backup/regression
cp -a golden golden_pre_L6b_20260914          # archive first, then verify
diff -r golden golden_pre_L6b_20260914
./run_check.sh check  wasp_full wasp_he23off wasp_full_newton \
                      mol_metals mol_ir_bands lower_profile oxygen_chemistry
./run_check.sh golden wasp_full wasp_he23off wasp_full_newton \
                      mol_metals mol_ir_bands lower_profile oxygen_chemistry
./run_check.sh check                           # the whole matrix, confirming
```

**This was NOT run, and this item refreshed nothing.** The refresh belongs at
the close of the series, on the delivered binary, not on this item's private
build.

Two things to carry into that refresh:

1. `wasp_full_newton` is the only one of the seven that converges. Its new
   golden pins a steady solution. The other six pin a transient, as they do
   now.
2. `oxygen_chemistry` will move again under any later change that reaches the
   base energy budget, for the reason in section 5.5. A refreshed golden for
   it buys less than the others do.

---

## 9. Reproduction

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
S=$EX/scratchpad/L6c

# the snapshot: HEAD plus the working tree, taken ONCE so that the file another
# worker is editing is the same copy in both builds
mkdir -p $S/snap && cd $EX && git archive HEAD | tar -x -C $S/snap
git status --porcelain -- src Makefile        # overlay each M and ?? entry

\cp -rf $S/snap $S/tree_M; \cp -rf $S/snap $S/tree_C
git show HEAD:src/modules/radiation/Cool_coeff.f90  > $S/tree_C/src/modules/radiation/Cool_coeff.f90
git show HEAD:src/modules/radiation/util_ion_eq.f90 > $S/tree_C/src/modules/radiation/util_ion_eq.f90
diff -rq $S/tree_M/src $S/tree_C/src          # must name exactly those two files

# build with the gfortran of the default PATH (conda-forge 16.2.0), the one the
# goldens were taken with; /usr/bin/gfortran 13.1.0 gives different objects
( cd $S/tree_M && make OBJDIR=build_L6c_x EXE=EXHALE_L6c_M.x -j8 )
( cd $S/tree_C && make OBJDIR=build_L6c_x EXE=EXHALE_L6c_C.x -j8 )

# each case on a scratch copy, single threaded, with its own maxsteps, at most
# six at a time; then compare against backup/regression/golden/
$S/run_case.sh <case> <binary> <tag>          # the runner used here
$S/attrib.sh   <case>                         # the four comparisons of table 2

# the held-state pair of section 4: one equilibrium solve from the control
# state, both binaries, "Do only PP: True"
$S/pp_case.sh  <case> <binary> <tag>
```

A note for whoever repeats this. The state written by an older build cannot be
reloaded by today's: the restart metadata check refuses the golden pair
because its option string predates the `wellbal` token. The held-state
comparison of section 4 therefore starts from the control run's own final
state, which both builds accept.

The private builds (`tree_M/`, `tree_C/`, their `build_L6c_x/` and the two
executables) and the isolated trees were deleted at the end of the item. No
process of this item was left running.
