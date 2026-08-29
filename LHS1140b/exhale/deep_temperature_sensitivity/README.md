# What a 2.1e-6 K difference in the deep boundary temperature does to the handoff column

Measurement record, 2026-08-29, this directory: README, six scripts, three
tables, thirteen adapter runs under `runs/`. LHS 1140 b, He/H = 9.0, the
configuration of `../nh_refusal_diagnosis/scan_reservoir_pc090.sh` (deep
boundary 20 bar, 60 climate layers, Zahnle H/He/N/O/C, `K_zz = 1e9`, the
GJ 1132 SED at the b orbit, trial fluxes of the stored crossing). Interpreter:
the repository's own `env/photochem/bin/python` (photochem 0.9.0). Every run
here is chemistry plus climate only -- no wind step, no closure iteration.

The observation this record answers is
`../crossings_pc090/results.txt` section 7: two builds put the deep boundary
2.1e-6 K apart and the handoff column came back differing by up to 2.4 percent
in pressure, 15 percent in one trace species and 0.5 percent in the elemental
O/H, with the mechanism left open.

## What was measured

1. `climate_guess_ladder.py` -> `climate_guess_ladder_heh9.txt`
   The climate solve alone, from ten initial guesses for the deep boundary
   temperature spanning 300-450 K. The solve is the same physical problem in
   every row; only the starting point of the root solve moves.

2. `run_guess.sh` -> `runs/g4*`
   The full adapter at seven of those guesses, one handoff column each.
   `g400a` and `g400b` are the same command run twice.

3. `run_tightened.py`, `run_tight.sh` -> `runs/t1em4_*`
   The same four of those runs with the photochemical stopping test tightened:
   `conv_longdy` 1e-2 -> 1e-4, `conv_longdydt` 1e-6 -> 1e-8, `equilibrium_time`
   1e17 -> 1e22 s. The adapter itself is unchanged; the wrapper only replaces
   `EvoAtmosphereGasGiant` with a subclass that states those three numbers.

4. `compare_columns.py`, `spread_table.py` -> `spread_table.txt`,
   `spread_table_tightened.txt`
   The columns compared three ways: at fixed level index, at fixed pressure,
   and at the matching level.

5. `runs/fix_g400`, `runs/fix_g400_1em6`
   Two of the same runs repeated after `insert_level` in
   `src/utils/lower_profile_schema.py` was corrected to recompute `n_tot` and
   `rho` at the inserted node from `p` and `T` instead of interpolating them.

## What the measurements say

**The 2.4 percent is the extent of Photochem's altitude grid at the moment the
run is declared steady. It is not an amplification of 2.1e-6 K, and it is not
the chemistry's convergence tolerance.**

- The climate solve reproduces its root to 2.121e-06 K over the whole guess
  ladder (424.95340719 - 424.95340876 K), its tropopause to 2.5e-8 relative
  and the water above it to 6.6e-8. Run twice in separate processes the ladder
  is bit-identical, and `g400a` / `g400b` write byte-identical profiles: the
  code is deterministic and none of this is run-to-run scatter.
- Every quantity of the written column is **two-valued**, and the label is the
  model top pressure: 1.0934e-08 bar or 1.1207e-08 bar, 2.5 percent apart.
  Within a group the columns agree to 1e-6 or better; between groups the
  pressure at fixed level index differs by 2.470e-02, the elemental O/H handed
  over at the matching level by 5.2e-03, and `q_OH` at fixed index by 0.58.
- The difference is a grid, not a state. `vertical_grid`
  (`photochem/src/photochem_eqns.f90`) lays the levels out uniformly in
  ALTITUDE between a fixed bottom and `top_atmos`, so a different `top_atmos`
  stretches every cell: measured, the pressure difference is zero at the pinned
  bottom (16.73098 bar in both) and grows linearly in level index to 2.5
  percent at the top.
- `top_atmos` is re-pinned to the requested top pressure only every
  `freq_update_TOA = 1000` internal steps, and the exit test accepts any top
  pressure within a factor 3 of it (`TOA_in_range` in
  `photochem/extensions/gasgiants.py`). Which side of that cycle a run exits on
  is what the 2.1e-6 K decides.
- **The chemistry tolerance is excluded by measurement.** With `conv_longdy`
  100x tighter, `conv_longdydt` 100x tighter and `equilibrium_time` 1e5x
  longer, both groups sit exactly where they were (`spread_table_tightened.txt`:
  pressure 2.466e-02, O/H at the match 5.199e-03) and each tightened run
  reproduces its loose-tolerance twin to 4e-7 in O/H.
- **The matching pressure does not move.** It is an exact inserted node and
  reads `1.000000e-06` bar in every profile written here. The 2.4 percent is a
  comparison at fixed level index; EXHALE reads the file as a table over
  pressure, and at fixed pressure the spread is 5e-4 in T, 1e-3 in `n_tot` and
  1.3e-2 in the oxygen near the cold trap.

## What it does to the base state EXHALE reads

At the matching level, over the seven runs: `T` 2.7e-09, `r` 4.8e-09, `q_H2`
5.7e-08, He/H 6.5e-08, C/H 3.1e-07, N/H 3.2e-07, O/H 5.2e-03. `n_tot` and
`rho` are carried by the file and explicitly not imposed
(`input_read.f90`). The only place the 2.4 percent reaches the wind is
`lap_r_top_RJ`, the upper edge of the elemental flux window, which moves by
8.2e-06 R_J = 5.2e-05 R_p.

## A defect found on the way, and corrected

`insert_level` interpolated every column linearly in log p, including `n_tot`
and `rho`, which are not free: `n = p/(k_B T)` holds at every node the solution
was computed on. Measured on the columns here, the identity was violated at
**exactly one level of 102 -- the matching level**, by 7.7e-04 and 1.24e-03,
and nowhere else. `insert_level` now recomputes `n_tot` from the inserted
pressure and temperature and `rho` from the interpolated mean molecular weight.
Verified by rerunning `g400a` (`runs/fix_g400`): two cells of the file change
and nothing else, the header and `solution_id` are unchanged, the identity
holds to 2.2e-16 across the whole table, and the run-to-run spread of `n_tot`
at the matching level falls from 4.65e-04 to 4.68e-10.
