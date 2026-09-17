# L9: a no-step evaluation that also produces the post-processed profiles

Item L9 of `docs/PLAN_20260913_lhs_stationary.md`, planned in
`docs/PLAN_20260916_rev3.md` section 9 with rows R19, R42, R49 and R50 of its
section 0 tables, and shaped by `docs/PLAN_20260916_review2.md` section 9.

Binaries, both built from the working tree of 2026-09-16 (git `43bc28c`) and
differing only in `src/EXHALE_main.f90` and
`src/modules/files_IO/write_output.f90`:

| | md5 |
|---|---|
| control, the entry text of those two files | `4197cca55a555815d6084e3878e6aed6` (`EXHALE_L9ctl.x`) |
| measured | `66500a8d7ffb6bdacaa2ecf99fef6892` (`EXHALE_L9.x`) |
| the tree binary the catalog was produced with | `c2e9c9990b9f14f1be8cd77abca68945` (`EXHALE.x`) |

The first two md5 are of the builds every number below was measured with. The
working tree is being edited by other items at the same time, so a rebuild of
the same two files a few minutes later links a different binary; what fixes
the comparison is that the control and the measured build differ in these two
files and in nothing else.

## 1. Verdict

**A state can now be measured and have every product of a run made ON that
measurement, with no step taken anywhere.** The `CFL 1e-12` step of the
runner's post-processing pass is gone, and with it three things it was
carrying that nobody had asked for: it was taken with the reconstruction
`input.inp` names rather than the one the solution was reached with, it
overwrote the solved state's files with a relaxation snapshot written
`certified=F cert_reason=no_stationary_claim`, and the mass row the
post-process weighs every row of its own product by came out of that wrong
operator, four to nine decades above the state's own (MEASURED, section 5).

Three states are named on the route and every product says which one it
describes. The work state is thermodynamically defined for the first time:
the conserved density, momentum and total energy are held and the pressure
and the temperature follow the refreshed composition through the caloric
equation of state. On an atomic gas that is the same map as the rule it
replaces, bit for bit; where H2 is present it is not (section 4).

## 2. The three states, and their contracts

| state | what is held | what it is | which product describes it |
|---|---|---|---|
| **loaded** | everything | the conserved variables and the composition the two `_IC` files carry, copied at entry and not written to again | the file's own `certified=`/`cert_reason=` pair, and the first of the two answers of section 3 |
| **work** | the conserved density, momentum and total energy of the loaded state | the composition one equilibrium sweep returns from the loaded one, with p and T derived from THAT composition at THAT conserved energy | `Hydro_ioniz.txt`, `Ion_species.txt`, `Cooling_breakdown.txt`, `Heating_breakdown.txt`, the residual breakdown, the certification, the mass-loss line |
| **advection-derived** | nothing of the above: it is another composition | what `post_process_adv` builds from the work state, row by row, where its validity conditions allow | `Hydro_ioniz_adv.txt`, `Ion_species_adv.txt` |

**Why the work state holds u and not p.** The evaluation is an evaluation of
the state the file carries, so the conserved variables are the ones that may
not move. The pressure of a molecular cell is the inverse of the caloric
equation of state and therefore a function of the composition at a given
thermal energy: holding p across a change of composition assigns that
composition a thermal energy which is NOT the file's, which is a change of
the conserved state and not a measurement of it. This is the contract
`pressure_and_temperature_at_fixed_conserved_state`
(`diffusive_photochemistry.f90`) states for the molecular relaxation, and the
evaluation now uses the same one. The entry text held p and recomputed T from
it alone.

**"Refreshed" is not "closed".** One sweep is a single Picard step of a
nonlocal coupling. The work state reports its own closure defect, in two
parts, both of them quantities the route already forms:

- **chemical**: the largest normalized reaction residual the sweep accepted a
  cell state at (`acc_resmax` of the sweep ledger, the same record the
  certification reads its rootless-cell count from), so it is the residual on
  the composition the sweep RETURNED and not on the one it was given;
- **thermal**: `max |T_work - T_sweep|/T_sweep`, the distance between the
  temperature each cell was solved at and the temperature the returned
  composition has at the held conserved energy.

MEASURED on the certified atomic fiducial
`LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH2.13`: chemical 4.408e-17 with
0 cells without a chemical root, thermal 5.959e-14. On the molecular
`backup/regression/golden/mol_base_handoff` state: chemical 2.195e-17,
thermal 2.828e-10.

## 3. Two answers, kept apart

The route prints, side by side:

```
 (EXHALE_main) the two answers of this evaluation:
   original claim of the file: certified=T cert_reason=certified_in_wind
   original claim: REPRODUCED -- the evaluation of the state the file carries certifies it
   work state verdict: CERTIFIED -- every active equation of the inventory is within its tolerance
   the pair written into the state: certified=T cert_reason=certified_in_wind
```

The first two lines answer a question about the FILE: the pair its
`# coupling:` header states (`ic_certified`, `ic_cert_reason`, item L23) is a
stationary claim about the state it carries, and re-measuring that state
either confirms or refuses it. A file that claimed nothing has nothing to
reproduce and the line says so. The last two answer a question about THIS
run: what the certification measured, and what pair is written into the state
this run emits. They are printed as two statements because the written reason
is not always the verdict's: a run measuring a state whose file claimed
nothing refuses nothing, so it writes `no_stationary_claim` however the
entries came out, and the entries that refuse the work state are named by the
certification report above it. MEASURED on the molecular state, which claims
nothing and whose three hydrodynamic rows refuse it:

```
   original claim: NONE MADE -- the file states no stationary claim, so nothing is reproduced or refused
   work state verdict: NOT CERTIFIED -- the certification report above names the entries that refuse it
   the pair written into the state: certified=F cert_reason=no_stationary_claim
```

## 4. The p-versus-u contract, measured

The two rules, on the work state's own composition: `p(keep u)` is the
pressure the caloric map assigns the refreshed composition at the loaded
conserved energy, `p(keep p)` is the loaded state's pressure; the temperature
follows each of them at the same refreshed particle count.

**The certified atomic fiducial** (`atomic_scalar_gj1132_kzz1e9/HeH2.13`),
cells 1 to 6, code units:

| cell | p(keep u) | p(keep p) | dp/p | T(keep u) | T(keep p) | dT/T |
|---|---|---|---|---|---|---|
| 1 | 9.4139417847e-01 | 9.4139417847e-01 | 0 | 1.05311242 | 1.05311242 | 0 |
| 2 | 8.8950914988e-01 | 8.8950914988e-01 | 0 | 1.14033784 | 1.14033784 | 0 |
| 3 | 8.4406653042e-01 | 8.4406653042e-01 | 0 | 1.23052907 | 1.23052907 | 0 |
| 4 | 8.0390324087e-01 | 8.0390324087e-01 | 0 | 1.31717906 | 1.31717906 | 0 |
| 5 | 7.6800163653e-01 | 7.6800163653e-01 | 0 | 1.39943453 | 1.39943453 | 0 |
| 6 | 7.3559797511e-01 | 7.3559797511e-01 | 0 | 1.47764116 | 1.47764116 | 0 |

and 0 over the whole column. **That is not a null result, it is the
statement the contract makes**: with no H2 anywhere, `caloric_mixture_active`
is false and the map is `(gamma_ad - 1) rho e`, which carries no composition,
so the two rules are one map and the fiducial's numbers are unchanged to the
last bit. The rule bites only where a composition change moves the heat
capacity.

**Where it bites** (`backup/regression/golden/mol_base_handoff`, an 80 per
cent H2 layer), same columns:

| cell | p(keep u) | p(keep p) | dp/p | dT/T |
|---|---|---|---|---|
| 1 | 5.3934703855e-01 | 5.3934703855e-01 | -4.879e-14 | -4.881e-14 |
| 2 | 5.3523592991e-01 | 5.3523592991e-01 | -4.895e-14 | -4.894e-14 |
| 3 | 5.3114769681e-01 | 5.3114769681e-01 | -4.996e-14 | -4.995e-14 |
| 4 | 5.2708352235e-01 | 5.2708352235e-01 | -5.076e-14 | -5.084e-14 |
| 5 | 5.2304350217e-01 | 5.2304350217e-01 | -5.137e-14 | -5.129e-14 |
| 6 | 5.1902765686e-01 | 5.1902765686e-01 | -5.262e-14 | -5.264e-14 |

with a column maximum of **5.929e-10 at cell 198**, four decades above the
base cells: the difference follows where the sweep moves x(H2), not where H2
is most abundant. It is small on this state because the state is near its own
chemical fixed point; it is not a rounding artefact, and on a state further
from it the sweep moves the heat capacity by as much as it moves the
composition.

## 5. The two routes, on the certified atomic fiducial

Scratch copies, `OMP_NUM_THREADS=1`, the same solved state
(`output/Hydro_ioniz_IC.txt`, `Ion_species_IC.txt`) handed to both:

- the present route: `Load IC? True`, `Do only PP: True`, `Restart intent`,
  `Solver` and `du_th` deleted, `CFL: 1.0e-12` appended, tree binary
  `c2e9c9990b9f`;
- the new route: `Load IC? True`, `Restart intent: stationary evaluate`,
  nothing else changed, `EXHALE_L9.x`.

They are two different operations -- one step against no step, and PLM
against WENO3 -- so a difference is expected; none of it is judged by a
certification tolerance.

**The state files written back over the solved state.** MEASURED, largest
relative difference over the 500 physical cells:

| file, column | max rel. difference |
|---|---|
| `Hydro_ioniz.txt` rho | 1.8e-16 |
| `Hydro_ioniz.txt` v | 2.0e-12 |
| `Hydro_ioniz.txt` p | 6.0e-14 |
| `Hydro_ioniz.txt` T | 1.3e-14 |
| `Hydro_ioniz.txt` heat, cool | 2.1e-13, 2.8e-13 |
| `Ion_species.txt` H I, H II, He I, He II, He III, He 2^3S | 3.5e-13 or below |

and the headers differ in what they say about the state:

| | present route | new route |
|---|---|---|
| `Hydro_ioniz.txt` | `recon=PLM certified=F cert_reason=no_stationary_claim` | `recon=WENO3 certified=T cert_reason=certified_in_wind` |
| `Hydro_ioniz_adv.txt` | `# coupling: ... certified=F cert_reason=no_stationary_claim` | `# derived_from: certified=T cert_reason=certified_in_wind`, and no coupling line |

The first row is the finding that matters for the catalog: **the solved state
of every case in `LHS1140b/models/` has been overwritten by the
post-processing pass with a snapshot that claims nothing**, and the
certification of the solve survives only in the `_IC` copy beside it. The
same file re-read is certified again on the new route.

**The advection-corrected profiles.** The measure the post-process weighs
every row by, `adv_mass_row`:

| cell | present route | new route |
|---|---|---|
| 1 | 1.5731e+00 | 1.7338e-09 |
| 5 | 1.3542e-01 | 7.7533e-10 |
| 26 | 1.4464e-05 | 9.5855e-11 |
| 200 | 4.1550e-06 | 9.0005e-13 |
| 500 | 3.0840e-05 | 3.9957e-14 |

The present route's column is the mass row of a state assembled with PLM,
which the solution is not a solution of; the new route's is the mass row of
the state itself, and it is the same number the certification judged
(1.733811e-09, the k=1 row of the evaluation). Nine decades at the base, four
in the wind.

That changes which base rows are corrected and by how much:

| quantity | max rel. difference, physical cells | where |
|---|---|---|
| `Hydro_ioniz_adv` p, T | 5.80e-01 | cell 5, r = 1.00097 |
| `Hydro_ioniz_adv` cool | 2.64e-01 | cell 5 |
| `Hydro_ioniz_adv` heat | 3.19e-03 | cell 5 |
| `Ion_species_adv` He 2^3S | 2.38e-02 | cell 26 |
| `Ion_species_adv` H II | 7.52e-03 | cell 26 |
| `Ion_species_adv` H I | 4.01e-07 | cell 500 |

**and it is confined to the base.** The 42 cells whose temperature differs by
more than one per cent lie between r = 1.00097 and r = 1.00889; above
r = 1.05 the largest difference in p or T is 4.4e-04, above r = 1.10 it is
6.3e-05, and above r = 1.20 it is 1.3e-05. The `adv_T_status` column changes
on three rows (cells 4, 5, 6): the present route retained them, because their
PLM mass row stood above the 1e-2 conditional tolerance, and the new route
corrects two and reports one as a failed cell solve. The composition status
column is identical on every row.

**The mass-loss line** is `log10 Mdot = 7.87` on both, to the two decimals
the line carries (the `rho v r^2` the two states carry at cell N-20 agree to
1e-12).

**The transit spectrum**, `EXHALE_transit.py` with the WINERED HIRES-Y
kernel, `MPLBACKEND=Agg`, on each of the two `_adv` pairs:

| | present route | new route | relative |
|---|---|---|---|
| He I 10830 red-pair equivalent width [%A] | 1.374922 | 1.374915 | 5.1e-06 |
| red-pair depth [%] | 4.957682 | 4.957660 | 4.4e-06 |

The line forms far above the cells that moved, so the spectrum does not see
the change.

**A finding this item does not settle.** The corrected temperature of cell 5
is 132.94 K against the 316.27 K the state carries there, and the base of
this model is at 226 K. The advection correction is a one-way steady
correction on a fixed (rho, v) field, and now that its own validity measure at
the base is 1e-9 instead of 1e-1 it corrects rows it used to retain. Whether
the correction is meaningful in the first ten cells of a subsonic base is a
question about `post_process_adv`, not about which state it is handed, and it
is left open (section 8).

## 6. What changed in the source

| file | what |
|---|---|
| `src/EXHALE_main.f90` | `stationary_state_of_the_loaded_restart`: the loaded state copied and named; the work state derived at fixed u through `pressure_from_energy_density`; the closure defect reported; `stationary_claim_and_work_state_verdict` (the two answers); the evaluate branch calls `write_excited_H`, `post_process_adv` and the mass-loss line; `evaluation_state_dump`, the in-memory diagnostic keyed by `EXHALE_EVAL_STATE_DUMP`. The mass-loss block of the marching route is now the contained `steady_mass_loss_rate`, called by both routes, so one text states that number |
| `src/modules/files_IO/write_output.f90` | the `'ad'` write states `# derived_from: certified=<T\|F> [cert_reason=...]` with the sentence that fixes it as provenance, in place of the `# coupling:` header it used to copy from the run state (`write_derived_provenance_header`) |
| `src/tests/grid_and_gates/stationary_evaluate_products.sh` | new, 20 assertions (section 7) |
| `src/tests/grid_and_gates/run.sh` | the `evaluate_products` row |
| `LHS1140b/models/run_case.sh` | the post-processing pass takes the evaluate route; `CFL 1.0e-12`, `Do only PP` and the deletion of `Solver`/`du_th` are gone; the two answers are read out of `pp.log` and reported |
| `LHS1140b/models/write_reproduce.py` | the two answers in the outcome list; the "Re-measuring" section describes the pass and the three states |

**What the `_adv` header no longer says, stated because it is a loss.** The
coupling line it used to copy carried `sec_ion`, `sec_ion_step`, `recon` and
`mode` beside the certification pair, and the `# derived_from:` line carries
the pair alone. Those four are properties of the run that produced the state,
and the state file stands in the same directory and states all of them; the
`_adv` files are never read as a restart, so nothing reads them from there.
Restoring them would mean a coupling header without its certification field,
which `write_coupling_state_header` does not offer and which is not this
item's file.

## 7. The tests

`src/tests/grid_and_gates/stationary_evaluate_products.sh`, run through the
suite as `evaluate_products`. The fixture is a copy of
`backup/regression/wasp_full_newton/IC/`, the certified pair the restart rows
of the same suite already read; the regression directory is never written to.
The output directory is created with a **stale sentinel** in each file a
producer must write, so a skipped producer is caught by its sentinel and not
only by a missing file.

The in-memory comparison is made through `EXHALE_EVAL_STATE_DUMP=1`, which
writes the conserved state, the composition and the live
composition-derived module state at four named points. An unchanged
`Hydro_ioniz_IC.txt` proves nothing, because the loader never writes one
back.

MEASURED, control binary (RED) and measured binary (GREEN), same script:

| assertion | control | measured |
|---|---|---|
| `produced_Hydro_ioniz`, `produced_Ion_species` | PASS | PASS |
| `produced_Hydro_ioniz_adv` | FAIL, `stale_sentinel` | PASS |
| `produced_Ion_species_adv` | FAIL, `stale_sentinel` | PASS |
| `produced_Cooling_breakdown`, `produced_Heating_breakdown` | PASS | PASS |
| `mass_loss_line_written` | FAIL, absent | PASS |
| `certification_metadata_written` | PASS | PASS |
| `adv_header_states_derived_from` | FAIL, 0 of 1 | PASS |
| `adv_header_makes_no_coupling_claim` | PASS on the sentinel | PASS |
| `adv_certified_only_as_provenance` | PASS on the sentinel | PASS |
| `in_memory_dump_has_four_points` | FAIL, no dump written | PASS |
| `loaded_state_untouched_in_memory` | (no dump) | PASS |
| `loaded_block_and_work_block_are_two_states` | (no dump) | PASS |
| `products_do_not_move_the_work_state` | (no dump) | PASS |
| `products_do_not_move_the_module_state` | (no dump) | PASS |
| `changed_option_original_claim_not_reproduced` | FAIL | PASS |
| `changed_option_work_state_verdict_reported` | FAIL | PASS |
| `changed_option_written_pair_reported` | FAIL | PASS |
| `changed_option_products_written` | FAIL, missing | PASS |

7 PASS and 9 FAIL with the control, **20 PASS and 0 FAIL with the measured
binary**.

Two of the rows pass with the control for a reason that is not the code's:
the negative assertions on the derived header find the stale sentinel, which
carries no `certified=` either. Their RED is the file the present route
actually writes, MEASURED on the fiducial:
`# coupling: sec_ion=T sec_ion_step=0 recon=PLM certified=F cert_reason=no_stationary_claim mode=init`
at the head of `Hydro_ioniz_adv.txt` -- a coupling header of the run state,
standing over rows that are another composition.

The last row of the table is the separation of the two answers: the state
claims `certified=T` and is evaluated under a stellar EUV luminosity twice
the one it was solved at, so the claim cannot reproduce; the run says
`original claim: NOT REPRODUCED` and states the work state's own verdict
beside it.

**Suites rerun with the measured binary** (`EXHALE_OBJDIR=build_L9`), every
suite whose driver links `write_output.o` or reads an `_adv` file:

| suite | measured binary |
|---|---|
| `grid_and_gates` (the whole suite, including the new row) | 240 PASS, 1 FAIL |
| `adv_static_limit` | 53 PASS, 0 FAIL |
| `transit_census` | 21 PASS, 0 FAIL |
| `certification` | 84 PASS, 0 FAIL |
| `run_mode` | 31 PASS, 0 FAIL |

The one `grid_and_gates` failure is
`outer_iteration_ending_is_the_stagnation_one` of
`output_state_consistency.sh`, which fails identically with the control
binary and belongs to whatever moved the outer iteration; it is the same row
`docs/lhs1140b_stationary_L23_20260916.md` section 4 reports.

## 8. What remains

- **The advection correction at a subsonic base.** Section 5: with the mass
  row now measured by the operator the state solves, the post-process
  corrects rows at r < 1.01 that it used to retain, and cell 5 comes back at
  132.94 K against the state's 316.27 K over a 226 K base. Whether the
  one-way correction is meaningful there is a question about
  `post_process_adv`'s closure and not about which state it is handed.
- **The catalog's state files.** Every `output/Hydro_ioniz.txt` in
  `LHS1140b/models/` was written by the old post-processing pass and states
  `certified=F cert_reason=no_stationary_claim recon=PLM`; the certified
  state is the `_IC` copy beside it. Nothing was rewritten by this item. A
  re-run of a case now leaves the certification in both.
- **The closure defect is reported and not gated.** Nothing refuses a state
  because its work state is far from closed; the numbers are printed and the
  certification judges what it judged before.

## 9. Reproduce

```bash
EX=/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00
cd $EX && make OBJDIR=build_L9 EXE=EXHALE_L9.x

# the two routes, on a scratch copy of the certified atomic fiducial
S=$EX/LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH2.13
for d in pp_old eval_new; do
   mkdir -p $d/output
   cp $S/input.inp $d/
   cp $S/output/Hydro_ioniz_IC.txt $S/output/Ion_species_IC.txt $d/output/
done
( cd pp_old && sed -i -e 's/^Load IC?.*/Load IC? True/' \
     -e 's/^Do only PP:.*/Do only PP: True/' -e '/^Restart intent:/d' \
     -e '/^Solver:/d' -e '/^du_th /d' input.inp
  echo 'CFL: 1.0e-12' >> input.inp
  OMP_NUM_THREADS=1 $EX/EXHALE.x > pp.log 2>&1 )
( cd eval_new && sed -i -e 's/^Load IC?.*/Load IC? True/' \
     -e 's/^Restart intent:.*/Restart intent: stationary evaluate/' input.inp
  OMP_NUM_THREADS=1 EXHALE_EVAL_STATE_DUMP=1 $EX/EXHALE_L9.x > run.log 2>&1 )
# (the Spectrum file key of the case is relative to the case directory; make
#  it absolute in the copies)

# the transit spectrum of each
. $EX/LHS1140b/winered_hires_y.sh
MPLBACKEND=Agg PYTHONPATH=$EX python3 $EX/EXHALE_transit.py

# the p-versus-u contract and the in-memory comparison, from the dump
#   eval_new/output/eval_state_dump.txt, four blocks:
#   loaded, work_before_products, work_after_products, loaded_kept

# the tests
EXHALE_EXE=$EX/EXHALE_L9.x bash $EX/src/tests/grid_and_gates/stationary_evaluate_products.sh
```
