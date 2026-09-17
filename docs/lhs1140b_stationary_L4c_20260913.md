# L4c: the run route of the molecular LHS 1140 b cases

Item L4c of `docs/PLAN_20260913_lhs_stationary.md`. What is asked is a route
that takes a `molecular_*` case of `LHS1140b/MODELS.md` section 3 to a
certified stationary state, as the atomic prescribed cases reached one on
2026-09-13 (`docs/lhs1140b_stationary_L5c_20260913.md`, MODELS.md section 6).

Every number is MEASURED on this tree unless it says READ. No source file was
touched.

## 1. Verdict

**The molecular cases are still blocked, and the reason L4 recorded is not the
reason.** Three statements, in the order they were measured:

1. **The base H2 fraction standing at its ceiling does not cause the cold
   march to stop.** Loosening it by a factor 12 to 19 (f = 0.9998 in place of
   the photochemical column's 0.99998) moves the stop from step 1199 to step
   1396 and changes nothing else. The stopping cell is a ghost cell below the
   base face whose element budget drifts monotonically out of chemical
   equilibrium from step 984 on. The column's own f is therefore kept.
2. **A second, independent stop sits at the He-rich end of the ladder**:
   the well-mixed He/H = 0.083 case never reaches the chemical stop, because
   at step 980 the energy update of a cell at 7.1e4 K that still counts as
   molecular has no bracket below the 5e4 K top of the H2 rovibrational
   table, and the run exits without writing a state.
3. **The stationary route from a clean snapshot fails on the hydrodynamic
   rows, not on the chemistry.** A step-900 snapshot (below both stops) has a
   mass-flux spread `du` = 9.58, and the partitioned route from it ends
   `info = 1` in all four configurations tried. The obstruction is mass 1.0
   to 1.4, momentum 1.7 to 2.0 and energy 1.5 to 2.0 at cells 330 to 420 --
   the same obstruction the atomic COLD-march snapshots show, which the
   atomic cases escape only by starting from an archived solved state. The
   molecular groups have no archived state and, by the restart contract, can
   never be given one built from an atomic run.

So the molecular blocker is the absence of a near-stationary molecular seed,
plus two domain limits of the molecular closures at the base and at the
He-rich end. None of them is closed by an input change.

## 2. What was run

Binary `EXHALE.x`, md5 `97e10317a710b9ccc63addbedde3586a`, the build of
2026-09-13 03:20 the atomic campaign is running on. 8 threads a run, at most
two at a time, alongside that campaign's 64 threads (system load average 50
to 85 on 72 cores).

Scratch tree `LHS1140b/models/.L4c/`. The case directories under
`LHS1140b/models/molecular_*` carry only input files.

## 3. The base H2 fraction is not what stops the cold march

The hypothesis L4 left was that `q_H2_base` standing 2.0e-5 below the ceiling
`0.5/(0.5 + He/H)` leaves the base element budget no room. To test it,
`make_models.py` was set to hold the hydrogen-nucleus fraction bound into H2
at f = 0.9998, which puts `q_H2_base` 2.1e-4 (He/H = 9.7) to 3.7e-4
(He/H = 0.083) below the ceiling, and the `molecular_scalar_*` `base.inp`
files were rewritten from it.

MEASURED, `molecular_scalar_gj1132_kzz1e9/HeH2.13`, cold two-stage march,
`Well balanced: True`, `EXHALE_MAXSTEPS=3000`:

| f | relative room below the ceiling | stop | cell | element violation | tolerance |
|---|---|---|---|---|---|
| 0.99998 (READ, L4) | 2.0e-5 | step 1199 | the base cell | 3.4e-4 | 1e-6 |
| 0.9998 | 2.4e-4 | step 1396 | ghost cell -1, r = 0.999807 R_p | 2.395e-4 | 1e-6 |

Both are `ERROR STOP ioniz_eq: persistent non-root chemical equilibrium`
after 1000 consecutive non-root sweeps of one cell
(`ionization_equilibrium.f90`, `nonroot_streak_update`). Twelve times the
room buys 197 steps, so the proximity to the ceiling is not the cause and the
constant was put back to the column's value.

What the run does show is the two cells below the base face leaving chemical
equilibrium and never returning. MEASURED at f = 0.9998, the element
violation of ghost cell -1:

| step | 985 | 1006 | 1046 | 1086 | 1126 | 1166 | 1396 |
|---|---|---|---|---|---|---|---|
| violation | 6.5e-6 | 2.7e-5 | 6.6e-5 | 1.0e-4 | 1.4e-4 | 1.6e-4 | 2.4e-4 |

The first sweep above the 1e-6 tolerance is at step 984 (cell 0) and 985
(cell -1); no sweep of either cell is a root afterwards (933 reported
acceptances at cell -1, 927 at cell 0). At the stop the cell is at
T = 271.06 K, n_tot = 4.798e9 cm^-3, n_e = 2.722e6 cm^-3, with stage
fractions (the `System_HeH_mol` layout)

    x(H+)/n_H 4.29e-4   2n(H2)/n_H 0.99956   H in H2+ 4.0e-7
    H in H3+ 8.5e-6     H in HeH+ 1.3e-6
    n(He+)/n_He 2.33e-3 n(He++)/n_He 1.4e-6  n(He 2^3S)/n_He 4.8e-6

so 6.6e-7 of the hydrogen is left atomic and the violation stands four orders
above it. This is a base-boundary statement about the molecular network at
270 K and 5e9 cm^-3, not about `q_H2_base`.

## 4. A second stop, at the He-rich end of the ladder

MEASURED, `molecular_scalar_gj1132_wellmixed/HeH0.083` (no element
diffusion), the same cold march: no chemical-equilibrium stop at all, and
instead

```
ATTEMPTED STEP EXHAUSTED at step 980 after 9 attempts.
  the operation that refused it: the energy update
  energy update FAILURE, step 980, 4 cell(s), first cell 228 (T 7.06E+04 K)
  reason: no bracket below the equation-of-state ceiling,
          |R|/scale 2.9179E-01, last iterate 5.0000E+04 K,
          missing source -9.6781E-03 erg cm^-3 s^-1
```

with exit status 2 and no state written. The ceiling is
`T_ceiling_mol_K = 5.0d4` of `energy_semi_implicit.f90`, the top of the H2
rovibrational table, applied to every cell `caloric_eos` marks molecular; the
atomic ceiling beside it is 1e7 K. At He/H = 0.083 the wind heats a cell that
still counts as molecular to 7.1e4 K, above the table, and the update then
has no bracket. Nine halvings to dt/2^8 give the same reason, so it is not a
step-size failure. Reported as found; the fix is a source change and outside
this item.

## 5. The stationary route from a marched snapshot

Both stops arrive after step 980, so a snapshot below them is clean.
MEASURED with `EXHALE_MAXSTEPS=900`: both cases finish exit 0 and write a
full state, with mass-flux spread `du` = 9.576 (kzz1e9, He/H = 2.13) and
9.585 (well-mixed, He/H = 0.083) -- the flux `rho v r^2` flat to no better
than a factor 10, which is the marching behaviour item L2 measured for the
atomic wind of this planet.

Each snapshot was then entered into the partitioned stationary route with the
keys the atomic recipe uses (`Reconstruction scheme: PLM`, `Load IC? True`,
`Solver: Newton`, `Restart intent: stationary`,
`Secondary_ionization: Immediate`, `Well balanced: True`, no `du_th`),
`EXHALE_PTC_DTAU0=1.0`, 8 threads, 40 min cap.

### The pass table

| case, configuration | passes | hydro info | mass | momentum | energy | worst species row | verdict |
|---|---|---|---|---|---|---|---|
| kzz1e9 HeH2.13, carriers transported | 1 | 2 | 1.03E+00 | 1.85E+00 | 1.65E+00 | 9.91E-01 at cell 335, elemental transport He/H | REFUSED, `info = 1` |
| wellmixed HeH0.083, carriers transported | 4 | 2 | 1.42E+00 | 1.99E+00 | 1.49E+00 | 2.33E-02 at cell 464, carrier balance H2 | REFUSED, `info = 1` |
| kzz1e9 HeH2.13, `Molecular carrier transport: False` | 1 | 2 | 1.03E+00 | 1.70E+00 | 1.99E+00 | 9.93E-01 at cell 338, elemental transport He/H | REFUSED, `info = 1` |
| kzz1e9 HeH2.13, `Coupled carrier solve: True` | 0 | -- | -- | -- | -- | -- | 40 min cap reached inside the first JFNK |

against the gates 3.0e-12 (mass), 1.0e-08 (momentum), 1.0e-06 (energy) and
1.0e-05 (the gated species rows at r >= 1.2). The refusal texts are the
solver's own:

- kzz1e9, both carrier settings: *"the element composition relaxation found
  no admissible advance and its entry composition was restored, so no further
  pass could differ from this one"*, after the hydrodynamic solve reported
  *"no descent along the Newton or the damped Gauss-Newton direction"* at
  iteration 29 with `||grad merit||` = 7.4e6;
- well-mixed: *"the joint distance of the state (its largest entry over its
  own tolerance) has not fallen in 3 consecutive passes; the alternation is
  not approaching a joint fixed point at these step lengths"*. Its four
  passes are identical to three digits in every row while the trust radius
  halves from 1.0e-2 to 1.3e-3;
- coupled carriers: 39 JFNK iterations in 40 min with `||R||` 2.000 ->
  1.882, the worst row the momentum of cells 390 to 419 (r = 5.1 to 7.8).

The first certification, of the loaded snapshot, already names the same rows
(kzz1e9: mass 1.587 at cell 333, momentum 1.989 at cell 394, energy 1.548 at
cell 330, carrier H2 0.174, elemental transport 0.994), so the route neither
improves nor damages the state: it finds no descent at all.

### The variants, and what they cost to reach

The two carrier settings are not input switches on an existing state.
`carrier` and `carrier_newton` are layout tokens of the restart contract
(`load_IC.f90`): a state's option block must match the run's, and these two
may never be named on a `Restart option change:` line, because they decide
how many unknowns the state has. MEASURED: reloading the transported-carrier
snapshot under `Coupled carrier solve: True` is refused at startup with
*"carrier_newton: F -> T ... The state in the file is a state of another
equation set"*. Each variant therefore needed its own 900-step cold march.
Both marched cleanly to step 900. The march itself is insensitive to either
key: `du` at step 900 is 9.5763566287991679 with `Coupled carrier solve:
True`, digit for digit the reference value, and 9.5763531226384213 with
`Molecular carrier transport: False`.

For the same reason the route the brief listed as (b), a longer march to
10000 steps, cannot be taken: the march ends at step 1396 (kzz1e9) or 980
(well-mixed) with the current base, and only a base far off the ceiling
(`q_H2_base` = 0.15 at He/H = 2.13, 21 percent below it, READ from L4) marches
past 10000 -- with `du` 9.5 to 10.1 and a lower boundary that is no longer the
photochemical column's.

## 6. What blocks it, and the smallest change that would not

The blocking rows are the hydrodynamic rows of the wind between r = 2.5 and
7.8, at O(1) of their own scale, on a state whose mass flux is not flat to a
factor 10. That is not a molecular defect: the atomic cases show the same
thing from a cold-march snapshot (MODELS.md section 6, "the snapshots are not
near stationarity anywhere") and reach `info = 0` only from an archived
solved state mapped onto the current grid. The molecular groups have no such
state, and the restart contract forbids making one from an atomic solution.

The smallest change that would open the route is therefore a molecular
initializer inside the code: given a certified atomic solution of the same
He/H and the same boundary, write the molecular state it implies (the
carriers at the base fraction, decaying with the atomic hydrogen profile) so
that `Load IC? True` starts the molecular stationary route from an almost
stationary wind. This respects the contract, because the state would then be
written by the code that owns the layout, rather than a file assembled
outside it. It is a source change and is NOT made here.

Two further items, also source changes, are recorded and not made:

- the molecular equilibrium of the two cells below the base face at
  T ~ 270 K, n ~ 5e9 cm^-3 has no root the sweep can find once the wind has
  developed (section 3);
- the 5e4 K ceiling of the molecular caloric equation of state is applied to
  cells the wind has heated to 7e4 K, with no path for such a cell to leave
  the molecular classification (section 4).

## 7. What was changed in the tree

- `LHS1140b/models/make_models.py`: the comment at `H2_NUCLEUS_FRACTION` now
  records that the ceiling proximity was tested and is not the cause. The
  value is unchanged at 0.9999831600354949, and the seven
  `molecular_scalar_*/HeH*/base.inp` files it writes are byte-identical to
  what they held before (verified by rewriting them from the restored
  constant).
- `LHS1140b/MODELS.md`: section 4 records the f = 0.9998 test and its
  outcome; section 6's molecular block is rewritten around the measurements
  above.
- `LHS1140b/models/run_case.sh.molecular.patch` (new): the refusal of a
  molecular case, as a diff against `run_case.sh`. It is held as a patch and
  not applied, because the atomic campaign re-reads `run_case.sh` for every
  case and an in-place rewrite of it breaks the cases mid-read.

## 8. How to reproduce

From `LHS1140b/models/`, with `EX` the repository root:

```
# the cold march of a molecular case, to the stop
d=.L4c/kzz1e9_HeH2.13; mkdir -p $d/output
cp molecular_scalar_gj1132_kzz1e9/HeH2.13/{input.inp,base.inp} $d/
( cd $d && OMP_NUM_THREADS=8 EXHALE_MAXSTEPS=3000 $EX/EXHALE.x > march.log 2>&1 )

# the same, capped below the stop, to produce a snapshot
( cd $d && OMP_NUM_THREADS=8 EXHALE_MAXSTEPS=900 $EX/EXHALE.x > march.log 2>&1 )

# the partitioned stationary route from that snapshot
s=.L4c/kzz1e9_HeH2.13_stat; mkdir -p $s/output; cp $d/base.inp $s/
for f in Hydro_ioniz Ion_species; do cp $d/output/$f.txt $s/output/${f}_IC.txt; done
sed -e 's/^Reconstruction scheme:.*/Reconstruction scheme: PLM/' \
    -e 's/^Load IC?.*/Load IC? True/' -e '/^du_th /d' $d/input.inp > $s/input.inp
printf 'Restart intent: stationary\nSecondary_ionization: Immediate\n' >> $s/input.inp
( cd $s && OMP_NUM_THREADS=8 EXHALE_PTC_DTAU0=1.0 $EX/EXHALE.x > run.log 2>&1 )
```

A variant is the same with `Coupled carrier solve: True` or
`Molecular carrier transport: False` appended to BOTH the march input and the
stationary input, since the reload refuses the change otherwise.
