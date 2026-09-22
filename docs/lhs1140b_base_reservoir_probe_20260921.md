# LHS 1140 b: the base-discriminating probe

2026-09-21. The measurement `PLAN_20260920_rev9.md` sections 11.3 and 17.2
name as the one that would decide the open physical question: is the base
state of these refusing columns WRONG, or merely FAR FROM ITS FIXED POINT?
Two probes, both Mode R.

## 1. The answer, first

**Reading 2 of section 17.2.** The base reservoir statement is not what these
states fail on. A scan of the base reservoir particle count `ntot_bc` over
the whole family of reservoir states consistent with the prescribed base
pressure and temperature moves the base face mass flux monotonically from
-5998 wind means to +5638 wind means, through +1, and **the mass and energy
rows at cells 1 and 2 never fall below 1.0E-01 anywhere in that family, and
are of order unity at the crossing**. The lowest value the continuity row of cell 1
reaches anywhere in that family is 5.28E-01, against a tolerance of 3E-12 or
the cell's own floor, and the lowest the energy row of cell 1 reaches is
4.75E-01 against 1E-06. So no reservoir state at this base level closes these
rows: the base is far from a fixed point that the reservoir does not
determine, and a continuation is the instrument, not a boundary repair.

**What the reservoir DOES determine, measured exactly.** Forcing `ntot_bc` to
the count of the same composition with no H nuclei bound into H2 reproduces
the atomic reservoir density to 0.03 per cent and the atomic face velocity to
3.0 per cent, and the characteristic face condition follows it:

| | atomic configuration | molecular, as it stands | molecular, `ntot_bc` forced atomic |
|---|---|---|---|
| `ntot_bc` | (atomic) | 9.501493355011453E-01 | 1.000036022632164E+00 |
| derived `n0` [cm^-3] | 3.9061E+13 | 4.1112E+13 | 3.9061E+13 |
| reservoir density at the face [mH/cm3] | 1.397205E+14 | 1.468023E+14 | 1.396814E+14 |
| face velocity the condition returns [cm/s] | -1.126628E+01 | -1.381929E+02 | -1.092319E+01 |
| base face mass flux, in wind means | +1.000000 | -5998.470 | +19.799 |
| mass row, cell 1 / cell 2 | 8.87E-10 / 3.00E-09 | 1.0006 / 7.03E-01 | 1.0174 / 1.3444 |
| energy row, cell 1 / cell 2 | 4.14E-09 / 1.16E-08 | 9.84E-01 / 1.0581 | 1.2004 / 9.60E-01 |

The base face mass flux does not follow: it changes sign and falls by a factor
303, and stops at twenty times the wind flux, because that flux is a strongly
cancelling quantity at this base. And no setting brings a row into tolerance.

**The term export says the two checkpoints fail the same way.** At cell 2 of
both the atomic and the molecular checkpoint the mass row and the energy row
are carried by the LOWER FACE flux, with the upper face carrying 1.9 and 4.1
per cent of it and the radiative source two decades below the residual. The same term dominates on both sides, so this is one phenomenon
with two signs: the atomic checkpoint piles mass into cells 1 and 2 from
below, the molecular one drains it downward. Neither row can be closed by any
plausible change in a source term: `heat` would have to move by 2.0E+04 per
cent (atomic) and 9.0E+03 per cent (molecular) to close the energy row at
cell 2, while the lower energy face would have to move by 97.7 and 102.2 per
cent.

**Found on the way, and reported as a separate item.** In the molecular
configuration the collisional-network chemical heating raises the energy row
of the base region by two decades over its atomic value on a column with no
H2 in it: with `Molecular reaction heat: False` the energy row at cells 3, 4
and 5 falls from 9.91E-01, 9.85E-01, 9.77E-01 to 7.31E-03, 7.60E-03,
7.88E-03, and the energy scale at cell 3 falls from 4.638E-02 to 4.172E-04,
which is the atomic configuration's own 4.189E-04. It does not move the base
cells: the mass rows and the base face flux are unchanged to every digit
printed. This is a separate defect candidate and it is not the base
question.

## 2. Classification

Every number below is **Mode R** (rev9 section 3): the current operator
evaluating an immutable stored state, one stationary assembly per run, no
physical step taken, nothing published, nothing written into
`LHS1140b/models/`. Each run was made on a copy of the state in its own
scratch directory with its own `output/`, and carries the sidecar section 4.2
requires. Values read from a stored certificate or a source file are READ or
INSPECTED with the file and line; values this work produced are MEASURED.

Two of the six evaluations are the current operator on a state written under
another boundary model, `A0_atomic` and `B1_atomic_ckpt`: both those states
carry `..._contact_upwind_v2` (READ, their headers) and this binary solves
`..._ghost_fixed_point_seed_reservoir_row_v3`, so the boundary is rebuilt from
the physical column and this run's reservoir and the residual measured is not
the residual the state was written with (INSPECTED, `load_IC.f90:1278-1291`,
and the run log says so). The molecular checkpoint and the converted column
carry the `_v3` model the binary solves.

## 3. Identities

### 3.1 The operator

| | |
|---|---|
| path | `<scratch>/EXHALE_probe.x`, a private build (`make OBJDIR=<scratch>/obj EXE=<scratch>/EXHALE_probe.x -j8`) |
| md5 MEASURED | `39604c4237fe13eb9818676e2a977f5a` |
| the tree's delivered binary, for comparison | `EXHALE.x` md5 `89ec67149aa452eea3deeb19052f83d1` (MEASURED); it does contain the phase-5 export and does NOT contain the diagnostic of section 4 |
| repository HEAD | `93eed8667483`, working tree dirty (MEASURED, 60 entries) |
| source digest | `295cfd910c35bed8980cac269de8bde3`, the md5 of the sorted md5 list of the 245 `.f90`/`.inc` files of `src/` and `Makefile` (MEASURED 2026-09-21, WITH the diagnostic of section 4 applied) |
| the tree rebuilds to it | MEASURED: a second build from an empty object directory, of the source exactly as this work leaves it, gives the same md5 `39604c4237fe13eb9818676e2a977f5a`, so the delivered source text is the source behind every number below |
| compiler strings | `GCC: (conda-forge gcc 16.2.0-5) 16.2.0`; `GCC: (GNU) 8.5.0 20210514 (Red Hat 8.5.0-20)` (MEASURED, `strings -a`) |
| assembly selector | `EXHALE_RESID_QUAD` unset in every run, so `ROWS_PRODUCTION` in double |

Every evaluation reports, in its own export header (MEASURED, all six):

```
# assembly rows_kind=production_double ieq_sweep_state_kind=marching reconstruction=WENO3 momentum_branch=weno3_well_balanced
# flags well_balanced=T transport_active=F use_plm=F use_weno3=T recon_lambda_on=F
```

so the kind-generic assembly never fires on this route, the state is not
inside a reconstruction continuation, and the transport sources `Smom` and
`Sene` are exactly zero in every cell of every file (MEASURED).

### 3.2 The states

| run | state | `Hydro_ioniz` md5 MEASURED | `Ion_species` md5 MEASURED |
|---|---|---|---|
| `A0_atomic` | `atomic_photochem_gj1132_kzzprofile/HeH9/k04/output`, the certified atomic wind, in its OWN atomic configuration | `f55a75c2b1347177b217bcd8ffff1362` | `b6dd7275ed06caea52217baec63f8987` |
| `A1_mol_asis` | the `EXHALE_MOLECULAR_SEED_X2=0` conversion of that column (physical cells equal to it to 1.3E-15, H2 zero everywhere), in the MOLECULAR configuration | `97eae5c9e1204d9b4637aeff9c36bcf5` | `031cd879416bb1a05329fa72a20918f3` |
| `A2_mol_forced` | the same pair and the same configuration, with `EXHALE_BASE_PARTICLE_COUNT=atomic` | the same | the same |
| `A3_mol_forced_noreacheat` | the same, with `Molecular reaction heat: False` and `Restart option change: mol_heat` | the same | the same |
| `B1_atomic_ckpt` | `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7`, generation `g0002_20260919T004843Z_f7485b14` | `25e11562ec3cd1f6d204bf36f8e503eb` | `78ac49d7e1ce90a3a4f271785400efda` |
| `B2_mol_ckpt` | `molecular_photochem_gj1132_kzzprofile/HeH9`, generation `g0005_20260920T053700Z_f59de41d` | `348a6ab88a84b1b967d7582f14f9a723` | `e39c8734708d77894f578fe1df25283c` |

The two Probe B md5 pairs equal the ones the generation manifests record and
the ones `docs/lhs1140b_identity_records_20260921.md` MEASURED. The A1 pair is
the one `docs/lhs1140b_m2_conversion_audit_20260921.md` section 11 records for
`zero_p`, and the A0 pair is the atomic source that audit converted.

**The operator reproduces both certificates to every printed digit**
(MEASURED, against the values READ from each `certification.txt`):

| checkpoint | row | READ from the certificate, with its tolerance | MEASURED here |
|---|---|---|---|
| atomic `g0002` | mass | 9.806E-01 above 1.9E-08 at cell 2 | 9.806145E-01 at cell 2 |
| | momentum | 8.724E-06 above 1.0E-08 at cell 1 | 8.723527E-06 at cell 1 |
| | energy | 1.071E+00 above 1.0E-06 at cell 1 | 1.071147E+00 at cell 1 |
| molecular `g0005` | mass | 1.041E+00 above 3.4E-10 at cell 2 | 1.040949E+00 at cell 2 |
| | momentum | 1.450E-03 above 1.0E-08 at cell 1 | 1.449595E-03 at cell 1 |
| | energy | 1.021E+00 above 1.0E-06 at cell 1 | 1.020731E+00 at cell 1 |

The consumer of section 5 rebuilds all six of those measures from the
exported terms alone, using the code's own three row scales, and reproduces
every one of them to every printed digit. That is the check that its reading
of the scales is the assembly's.

### 3.3 How each evaluation was made, and the section 4.2 gate

```
env -i PATH=/usr/bin:/bin HOME=$HOME OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 EXHALE_RESIDUAL=1 EXHALE_CONSERVATION_BUDGET=1 \
    EXHALE_BOUNDARY_TRACE=1 [EXHALE_BASE_PARTICLE_COUNT=<value>] ./EXHALE_probe.x
```

in a directory holding a copy of the case `input.inp` with absolute paths and
`Restart intent: stationary evaluate`, the state pair under the names
`load_IC` reads, a copy of the executable and an empty `output/`.
`EXHALE_RESIDUAL=1` assembles the stationary residual of the loaded state once
and stops (`src/EXHALE_main.f90:1263-1341`), so exactly one assembly is
exported.

Each run directory carries `SIDECAR.txt` with the complete environment set
including the variables set to nothing, the mode, the input and state
identities, the executable md5, the evaluation stage, the exit class, and the
export's existence, byte count, row count and data md5. **All six exited 0 and
all six wrote `output/conservation_budget_0001.txt` with 503 rows.** The data
md5 of each export (the header carries a wall-clock line and cannot be
compared by md5):

| run | export data md5 MEASURED |
|---|---|
| `A0_atomic` | `37704a9d1829a5394940023f73447962` |
| `A1_mol_asis` | `a108978a99247aa20e1537e911885ada` |
| `A2_mol_forced` | `bebe04687da9f7a5bb6ec0cb0fe44bfe` |
| `A3_mol_forced_noreacheat` | `3d48f4ec228d327be578663665e579d2` |
| `B1_atomic_ckpt` | `bf951b72d1a8c59e0cabdbee07fae144` |
| `B2_mol_ckpt` | `d7b5cc7dd8a40c6c72086a27392c6c6b` |

## 4. The diagnostic that makes Probe A possible

### 4.1 Nothing already reaching `ntot_bc` would do

INSPECTED. `ntot_bc` has exactly one assignment,
`input_read.f90:1941`, `ntot_bc = comp_ntot_bc()`, and `comp_ntot_bc`
(`composition.f90:265-277`) is the single source of the base composition
policy: it returns 1, adds the metal nuclei when they are in the equation of
state budget, and subtracts `h2_bound_fraction()` when `molecular_base`. No
input key, no `base.inp` field and no environment variable reaches it
independently of the base composition:

- `Molecular base: False` cannot be stated, because `input_read.f90:1869-1873`
  sets `molecular_base = .true.` whenever the molecular chemistry is on and
  says so;
- `q_H2_base` changes the ghost composition and the particle count together,
  which is precisely the pair the probe has to separate;
- `Base BC: pressure` and the profile's matching level fix `p_base`, not the
  count.

So the probe needs a diagnostic, and it is added.

### 4.2 What was added

`EXHALE_BASE_PARTICLE_COUNT`, read once at `input_read.f90` immediately after
the base composition scalars are resolved and the molecular-base consistency
refusal has run, and documented there in full. It forces `ntot_bc`, the free
particles the base reservoir holds per (H+He) nucleus. Its value is either a
positive real, or `atomic`, which is `ntot_bc + h2_bound_fraction()`, the
count of the SAME composition with no H nuclei bound into H2, taken from the
one module that defines that fraction. When it is set the run writes

```
 (input_read) EXHALE_BASE_PARTICLE_COUNT: the base reservoir particle count is FORCED; this run is not the configuration its input.inp describes.
   ntot_bc from the base composition =  9.501493355011453E-001
   ntot_bc forced                    =  1.000036022632164E+000
```

so both counts are in the record at seventeen digits. Unset, empty,
unreadable or not a positive number leaves `ntot_bc` exactly as
`comp_ntot_bc` returned it, writes no line and touches nothing; a value that
is neither a positive number nor `atomic` stops the run rather than being
silently ignored.

**What it changes, stated exactly.** `ntot_bc` enters the run at two places
and the diagnostic sits before both, so both follow it:

1. the base level, `n0 = p_base/(k_B T0 ntot_bc)` (`input_read.f90:2084`,
   `:2150`): at the prescribed base pressure and temperature a smaller count
   means more nuclei, so the reservoir is denser;
2. the reservoir itself,
   `set_base_reservoir(ntot_bc + dp_bc, 1, (ntot_bc + dp_bc)/rho_bc, 1)`
   (`input_read.f90:2189`), and the pre-solve cell-1 count `n_part_cell1`.

Because `load_IC` reads the stored state in dimensional units and divides by
`n0` (`load_IC.f90:1154-1180`), and `R0`, `T0`, `v0 = sqrt(k_B T0/mu)` and
`b0` do not contain `n0`, **the physical interior column is held fixed under
this change**: what moves is the base reservoir's physical density at the same
`(p_base, T0)`, which is the one quantity the probe is asked to vary. That is
not only an argument. MEASURED, with `atomic` forced on the molecular
configuration against the atomic control: the interior velocity at the base
face is bitwise -1.484002E+00 cm/s in both; the volume-divided upper-face mass
contribution of cell 2 is +2.9899455293662942E-03 against
+2.9899455293662986E-03, the same to fifteen significant digits; and the area-weighted mass flux at the
certification gate and at the outer boundary is +5.7845966674E-07 and
+5.7845966675E-07 in both, to the ten digits printed. The derived `n0` is
3.9061E+13 cm^-3 in both, to the four figures the log prints, and the
reservoir density at the face is 1.396814E+14 against 1.397205E+14 mH/cm3,
0.03 per cent apart.

It does NOT change the ghost composition: the lower ghost still carries the
handoff H2 partition. That is deliberate, and it is what makes the probe a
one-variable experiment.

### 4.3 The regression proof, with the key unset

A scratch copy of the harness (`run_check.sh`, `compare_within_tolerance.py`
and the two case directories, copied out of the repository; the repository's
own `backup/regression/` was neither written to nor run):

```
REGRESSION_EXE=<scratch>/EXHALE_probe.x \
REGRESSION_GOLDEN_DIR=<repo>/backup/regression/golden \
    ./run_check.sh check hp_front mol_base_handoff
```

`mol_base_handoff` is the case that exercises the touched path with a
non-unit count (`Molecular base` on, `ntot_bc` below 1); `hp_front` is the
atomic control. MEASURED, the complete output:

```
[build] skipped: REGRESSION_EXE=<scratch>/EXHALE_probe.x (39604c4237fe)
[hp_front] running (OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=100) ...
            final: count=100  du= 6.6225E-01  dtu= 3.9604E-03
       PASS Hydro_ioniz.txt (data identical)
       PASS Ion_species.txt (data identical)
       PASS Hydro_ioniz_adv.txt (data identical)
       PASS Ion_species_adv.txt (data identical)
[mol_base_handoff] running (OMP_NUM_THREADS=1 EXHALE_MAXSTEPS=12000) ...
            final: count=12000  du= 7.9925E-01  dtu= 6.0928E-03
       PASS Hydro_ioniz.txt (data identical)
       PASS Ion_species.txt (data identical)
       PASS Hydro_ioniz_adv.txt (data identical)
       PASS Ion_species_adv.txt (data identical)
==> REGRESSION PASS (all cases byte-identical)
```

The verdict line, verbatim:

```
==> REGRESSION PASS (all cases byte-identical)
```

Eight files over two cases, byte-identical to their references with the key
unset. No golden was read for writing and none was refreshed.

## 5. The consumer

`base_rows.py` in the scratch directory, Python 3 standard library only. It
reads the export and forms the two rows again from the face fluxes, the face
areas, the cell volume, the potential and the sources, and never sums the
exported residual and calls that a rebuild. It applies geometry exactly once,
following the export's own `# kinds:` line: the `face_*` columns are FLUX
DENSITIES and are multiplied by that face's area, the `dF_*`, `S_*`, `R_*`,
`heat`, `cool`, `Smom`, `Sene` and `grav_work_over_volume` columns are
CONTRIBUTIONS ALREADY DIVIDED BY THE CELL VOLUME and are not touched, and the
common solid-angle factor 4 pi is omitted from the areas and the volume as
the assembly omits it.

MEASURED, over all six exports: the mass row and the energy row rebuilt from
the exported terms differ from the exported `R_mass` and `R_energy` by
between exactly zero and 1.4E-17, against rows whose own terms reach 1E+01.
Nothing is missing from the file and nothing in it is counted twice.

The three row scales are the code's own, INSPECTED at
`steady_residual.f90:392-424`, `:522-553`, `:713-785` and `:789-811`:

| row | scale |
|---|---|
| mass | `max(|A_hi F_hi|, |A_lo F_lo|)/V` |
| momentum | `max(|ram|, |pressure gradient|, |gravity|, |Smom|)` |
| energy | `max(|dF_3|, |S_3|, heat, cool, |Sene|)` |

and the tolerances are the certification's: momentum 1E-08, energy 1E-06 and,
for the continuity row, `max(3E-12, the cell's own rounding floor)`. That
floor needs the caloric sound speed and is not rebuilt here, so a continuity
refusal printed against 3E-12 is a LOWER BOUND; the two certificates READ
above give the binary's own value at cell 2 as 1.9E-08 and 3.4E-10, and the
measures stand seven and nine decades above them either way.

## 6. Probe A: the boundary alone, on the molecular case

### 6.1 The three states, held interior, varied reservoir

The interior is the same stored column in A0, A1 and A2, and the interior
velocity at the base face is -1.484002E+00 cm/s in all three, to the seven
digits the report prints (MEASURED, the `[base face]` block of each run).

| MEASURED | A0, atomic configuration | A1, molecular as it stands | A2, molecular with `ntot_bc` forced atomic |
|---|---|---|---|
| `ntot_bc` | (the atomic count) | 9.501493355011453E-01 | 1.000036022632164E+00 |
| derived `n0` [cm^-3] | 3.9061E+13 | 4.1112E+13 | 3.9061E+13 |
| lower ghost after `Apply_BC`, rho [n0 mu] | +3.6752678242E+00 | +3.6752976426E+00 | +3.6752805903E+00 |
| the same, in nuclei [mH/cm3] | 1.435601E+14 | 1.510988E+14 | 1.435606E+14 |
| lower ghost momentum density [n0 mu v0] | -3.2589459849E-04 | -3.9905375421E-03 | -3.1588188587E-04 |
| lower ghost energy density [p0] | +1.5003518494E+00 | +1.4721771635E+00 | +1.5490224903E+00 |
| lower ghost velocity [cm/s] | -10.967 | -134.289 | -10.630 |
| reservoir density at the face [mH/cm3] | 1.397205E+14 | 1.468023E+14 | 1.396814E+14 |
| reservoir temperature [K] | 182.104 | 181.984 | 182.156 |
| contact upwind `w_rev` | 0.0 | 0.0 | 0.0 |
| face velocity the condition returns [cm/s] | -1.126628E+01 | -1.381929E+02 | -1.092319E+01 |
| base face mass flux `A F_m` [code units] | +5.7845966877E-07 | -3.2967788400E-03 | +1.1452926825E-05 |
| the same, in wind means | +1.000000 | -5998.470 | +19.799 |
| gate face `A F_m` at r = 1.198946 | +5.7845966674E-07 | +5.4960326981E-07 | +5.7845966674E-07 |
| outer face `A F_m` | +5.7845966675E-07 | +5.4960326982E-07 | +5.7845966675E-07 |

The -5998.470 of A1 is the -5997 of `lhs1140b_m2_conversion_audit_20260921.md`
section 8, measured again here with a different instrument (the phase-5
export rather than the boundary report) and a different binary, and agreeing
to 0.02 per cent.

Two things the table says on its own. The reservoir statement IS the whole of
the difference in the boundary's own quantities: with the atomic count the
reservoir density comes back to 0.03 per cent of the atomic control's and the
face velocity to 3.0 per cent. And it is NOT the whole of the difference in
the flux the assembly then forms at that face: the base face mass flux changes
sign and falls by a factor 303, and stops at twenty wind means rather than
one. The face flux is a strongly cancelling quantity at this base, which is
why a 3 per cent change in the face velocity moves it by a factor 300.

### 6.2 The rows at cells 1 to 5

MEASURED, against the code's own scales; `|R|/s` is the certification's own
measure.

**A0, the atomic configuration on the atomic column (the control).**

| cell | r | mass `|R|/s` | momentum `|R|/s` | energy `|R|/s` |
|---|---|---|---|---|
| 1 | 1.0001933 | 8.873400E-10 | 3.550649E-14 | 4.137880E-09 |
| 2 | 1.0003866 | 3.004380E-09 | 1.518055E-14 | 1.159687E-08 |
| 3 | 1.0005800 | 1.823988E-10 | 1.256937E-15 | 4.734168E-09 |
| 4 | 1.0007733 | 4.668500E-10 | 4.363346E-17 | 4.897086E-09 |
| 5 | 1.0009666 | 2.216685E-10 | 1.939688E-15 | 2.431289E-09 |

**A1, the molecular configuration on the same column.**

| cell | r | mass `|R|/s` | momentum `|R|/s` | energy `|R|/s` |
|---|---|---|---|---|
| 1 | 1.0001933 | 1.000562E+00 | 3.123368E-02 | 9.837073E-01 |
| 2 | 1.0003866 | 7.033530E-01 | 1.787643E-05 | 1.058074E+00 |
| 3 | 1.0005800 | 7.641860E-04 | 7.243071E-09 | 9.910120E-01 |
| 4 | 1.0007733 | 4.126450E-10 | 1.139480E-14 | 9.852031E-01 |
| 5 | 1.0009666 | 5.064178E-10 | 3.977811E-15 | 9.774301E-01 |

**A2, the same with the atomic reservoir count forced.**

| cell | r | mass `|R|/s` | momentum `|R|/s` | energy `|R|/s` |
|---|---|---|---|---|
| 1 | 1.0001933 | 1.017394E+00 | 8.351394E-05 | 1.200432E+00 |
| 2 | 1.0003866 | 1.344389E+00 | 9.752029E-06 | 9.597193E-01 |
| 3 | 1.0005800 | 1.824006E-10 | 2.781692E-15 | 9.909660E-01 |
| 4 | 1.0007733 | 1.280866E-10 | 1.748731E-15 | 9.852031E-01 |
| 5 | 1.0009666 | 2.356138E-10 | 3.738174E-15 | 9.774301E-01 |

The momentum row is the one the reservoir does largely fix: 3.12E-02 to
8.35E-05 at cell 1, a factor 374, although still four decades above its 1E-08
tolerance. The mass row of cell 1 does not improve, and the mass row of cell 2
gets worse by its own measure while the integrated numerator below r = 1.03
falls from 3.300E-03 to 1.243E-05, a factor 265 (READ from each run's own
window block): the absolute imbalance falls and the scale it is measured
against falls faster, because the scale is the base face flux itself.

The mass row above cell 2 is at the arithmetic floor in A2 exactly as in A0,
and the upper face of cell 2 carries +2.9899455293662942E-03 against the
atomic control's +2.9899455293662986E-03, the same number to fifteen significant digits.
So with the atomic count forced, everything above cell 2 IS the atomic wind
and the disturbance is the two base cells alone.

### 6.3 The scan over the reservoir count, and what it decides

The family the plan asks about is "any reservoir state consistent with the
base pressure and temperature". That family is exactly one parameter, and
that parameter is `ntot_bc`. INSPECTED, `input_read.f90:2189`: the reservoir
is stored as `(p, T, nhat)` with `p = ntot_bc + dp_bc`, `T = 1` and
`nhat = (ntot_bc + dp_bc)/rho_bc`, so in physical units its pressure is
`p_base` and its temperature `T0` at every value of the count, while its
density is `rho_bc n0 mu` and `n0 = p_base/(k_B T0 ntot_bc)`. Varying the
count at fixed `rho_bc` is varying the particles a unit mass of the base gas
carries, which is exactly what the H2 binding does and exactly what a
reservoir at fixed `(p_base, T0)` has left to choose. MEASURED, fourteen values of it, the same column and configuration
throughout, `|R|/s` against the code's own scales:

| `ntot_bc` | base face, wind means | mass c1 | mass c2 | momentum c1 | energy c1 | energy c2 |
|---|---|---|---|---|---|---|
| 0.9501493355011453 (its own) | -5.998470E+03 | 1.0006E+00 | 7.0335E-01 | 3.1234E-02 | 9.8371E-01 | 1.0581E+00 |
| 0.96 | -4.786663E+03 | 1.0026E+00 | 9.2002E-01 | 2.4618E-02 | 9.7971E-01 | 1.2993E+00 |
| 0.97 | -3.584573E+03 | 1.0070E+00 | 9.6024E-01 | 1.8082E-02 | 9.7316E-01 | 1.6319E+00 |
| 0.98 | -2.396067E+03 | 1.0174E+00 | 9.7595E-01 | 1.1809E-02 | 9.6048E-01 | 1.9411E+00 |
| 0.99 | -1.186934E+03 | 1.0457E+00 | 9.8156E-01 | 5.9432E-03 | 9.2305E-01 | 1.7184E+00 |
| 0.9995 | -4.194052E+01 | 1.1709E+00 | 8.6050E-01 | 2.8156E-04 | 4.7489E-01 | 1.1564E+00 |
| 0.9998 | -8.066431E+00 | 1.3759E+00 | 6.7017E-01 | 7.8368E-05 | 8.7874E-01 | 1.0482E+00 |
| 0.99985 | -2.437026E+00 | 1.9540E+00 | 5.6986E-01 | 4.4332E-05 | 9.4610E-01 | 1.0297E+00 |
| 0.9999 | +3.417872E+00 | 5.2807E-01 | 3.8004E-01 | 9.4147E-06 | 1.0141E+00 | 1.0110E+00 |
| 0.99995 | +9.443581E+00 | 9.0507E-01 | 1.0352E-01 | 2.4716E-05 | 1.0826E+00 | 9.9227E-01 |
| 1.000036022632164 (atomic) | +1.979901E+01 | 1.0174E+00 | 1.3444E+00 | 8.3514E-05 | 1.2004E+00 | 9.5972E-01 |
| 1.01 | +1.221905E+03 | 1.1038E+00 | 1.0079E+00 | 6.6028E-03 | 1.0725E+00 | 7.0385E-01 |
| 1.02 | +2.467207E+03 | 1.0570E+00 | 1.0071E+00 | 1.1915E-02 | 1.0375E+00 | 7.3262E-01 |
| 1.05 | +5.637947E+03 | 1.0178E+00 | 1.0100E+00 | 2.8439E-02 | 1.0167E+00 | 6.2675E-01 |

The base face flux is monotone in the count and crosses +1 wind mean between
0.99985 and 0.9999, a hundredth of a per cent below the atomic count. At the
crossing the mass row of cell 1 reads 5.3E-01 to 2.0E+00 and the energy row
of cells 1 and 2 reads 1.01 and 1.01. **The rows never close.** Their minima
over the whole family are 5.28E-01 (mass, cell 1, at 0.9999), 1.04E-01 (mass,
cell 2, at 0.99995), 4.75E-01 (energy, cell 1, at 0.9995) and 6.27E-01
(energy, cell 2, at 1.05), reached at four different points, so not one
reservoir state is even close on all four. Every one of them stands at least
five decades above the tolerance of its row (energy 1E-06; continuity
`max(3E-12, the cell's own floor)`, whose value at cell 2 the two certificates
READ as 1.9E-08 and 3.4E-10). The momentum row is the one with a genuine
minimum at the crossing, 9.41E-06, three decades above its own 1E-08.

The scan varies the count and holds the ghost composition, which is the
molecular base's other statement; section 8 says what that leaves open.

### 6.4 What remains, measured

With the atomic count forced, the energy row is still of order unity not only
at the base cells but at cells 3, 4 and 5, where the mass row is at the
arithmetic floor. That is a source imbalance, and it is identified. In the
molecular configuration `heat` carries the collisional-network chemical
heating of the H2/He reactions (INSPECTED, `util_ion_eq.f90:2118-2130`,
`heat = heat + heat_chan(:,15)` under `with_molecules .and.
mol_reaction_heat`, default `.true.`, `parameters.f90:646`). Turning it off
with the existing production key, one variable, everything else held
(`A3_mol_forced_noreacheat`, MEASURED):

| cell | energy `|R|/s`, A2 | energy `|R|/s`, A3 | energy scale, A2 | energy scale, A3 | the atomic control's scale |
|---|---|---|---|---|---|
| 1 | 1.200432E+00 | 1.012561E+00 | 2.019739E-01 | 4.048365E-02 | 5.135638E-04 |
| 2 | 9.597193E-01 | 8.765127E-01 | 9.170036E-02 | 3.692041E-03 | 4.608341E-04 |
| 3 | 9.909660E-01 | 7.309266E-03 | 4.637561E-02 | 4.171902E-04 | 4.189493E-04 |
| 4 | 9.852031E-01 | 7.602696E-03 | 2.619805E-02 | 3.858351E-04 | 3.876435E-04 |
| 5 | 9.774301E-01 | 7.881074E-03 | 1.611160E-02 | 3.617766E-04 | 3.636292E-04 |

Above the base the energy row falls by two decades and the energy scale comes
back to the atomic control's own value to half a per cent. At the base cells
it does not: cell 1 stays at 1.013 and cell 2 at 0.877. And the change is
confined to the energy row: MEASURED over all 500 physical cells of the two
exports, `face_mass_lo`, `face_mass_hi`, `R_mass`, `R_momentum`, `rho`,
`momentum_density` and `energy_density` are BITWISE identical between A2 and
A3, and only `heat` and `R_energy` differ, by at most 2.015E-01. The base
face mass flux is +1.1452926825E-05 in both and the mass row measures are
1.017394E+00 and 1.344389E+00 in both.

So the base cells' failure survives both the reservoir count and the chemical
heating, and the large energy row of the layer ABOVE the base is a separate
matter: on a column with no H2 in it, the molecular configuration deposits a
chemical heat two decades above anything the atomic configuration has there,
with nothing balancing it. That is reported here and not pursued.

## 7. Probe B: the term export at cells 1 and 2 of both checkpoints

### 7.1 The convention, stated before any term

INSPECTED, `steady_residual.f90:296` and `:304`, and `RK_rhs.f90:273-275`:

```
R_mass   = dF_mass   - S_mass
R_energy = dF_energy - S_energy - (heat - cool) - Sene
dF3p     = dA_p F_m,p (Gphi_i(j)   - Gphi_c(j))
         - dA_m F_m,m (Gphi_i(j-1) - Gphi_c(j))
dF_energy = (dA_p F_E,p - dA_m F_E,m + dF3p)/dV
```

so **the gravitational work rides on the MASS face flux INSIDE `dF_energy`
with a positive sign**, and `S_energy` is zero in every branch
(`Source.f90:44-50`). No second gravitational work term is to be subtracted
anywhere. The energy zero is the caloric mixture's `1.5 T + x2 u_rv/T0`
(INSPECTED, `caloric_eos.f90:486`), translational energy plus the H2
rotational and vibrational excitation, not a binding energy measured from
separated atoms; the dissociation energy is in the chemical source
accounting. Transport was NOT active in either evaluation (`visc=F cond=F`,
`Smom` and `Sene` exactly zero in every cell), so the viscous and conductive
sources are absent rather than omitted.

Units: the mass row in `n0 mu / t_s`, the energy row in `p0 / t_s`, with the
five normalization constants in each export header.

### 7.2 The atomic checkpoint, `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` `g0002`

Base face `A F_m` = +4.2300976616E-07, gate face +5.1472049116E-09, outer
face +5.1472049117E-09, wind-window mean over 284 faces +5.1472049117E-09, so
the base face carries **+82.18242** wind means and the same over the gate
(MEASURED). That is the 82.18 of
`lhs1140b_conservation_audit_20260921.md` section 7.6 and the 82.2 of
`lhs1140b_stationary_D9fix_20260919.md` section 6, now a third time, with a
third binary. The three flux values above equal the conservation audit's to
every digit it printed.

The rows at the first five cells, MEASURED against the code's own scales:

| cell | r | mass `|R|/s` | momentum `|R|/s` | energy `|R|/s` |
|---|---|---|---|---|
| 1 | 1.0001933 | 3.723110E-01 | 8.723527E-06 | 1.071147E+00 |
| 2 | 1.0003866 | 9.806145E-01 | 2.680323E-06 | 1.004968E+00 |
| 3 | 1.0005800 | 4.176999E-09 | 1.246033E-15 | 1.847701E-08 |
| 4 | 1.0007733 | 1.936872E-08 | 2.293790E-15 | 8.819762E-08 |
| 5 | 1.0009666 | 2.786895E-08 | 6.785531E-16 | 1.523548E-07 |

Two cells, and then nothing: from cell 3 upward the rows are seven to nine
decades smaller, which is the localization
`lhs1140b_conservation_audit_20260921.md` section 7.6 reached from the
telescoping sums.

| term | cell 1 | cell 2 |
|---|---|---|
| **MASS ROW** | | |
| `-A_lo F_m,lo / V` | -2.1873005070828411E-03 | -1.3724139698733325E-03 |
| `+A_hi F_m,hi / V` | +1.3729445453627333E-03 | +2.6604901236554130E-05 |
| `-S_mass` | -0.0 | -0.0 |
| `= R_mass` | -8.1435596172010764E-04 | -1.3458090686367784E-03 |
| largest term | the lower face, -2.187301E-03 | the lower face, -1.372414E-03 |
| to close the row, the lower face would have to change by | +8.143560E-04, 37.23 % of itself | +1.345809E-03, **98.06 %** of itself |
| the same for the upper face | +8.143560E-04, 59.31 % | +1.345809E-03, 5059 % |
| **ENERGY ROW** | | |
| `-A_lo F_E,lo / V` | -1.0326812944999951E-03 | -9.3461574735278076E-04 |
| `+A_hi F_E,hi / V` | +9.3497706997003685E-04 | +1.1554636844312712E-05 |
| `+grav work dF3p/V` | +3.6537930161049454E-05 | +1.4353258979878131E-05 |
| `-S_energy` | -0.0 | -0.0 |
| `-heat` | -4.3669000095215123E-06 | -4.5304524549202689E-06 |
| `+cool` | +1.5102487891922102E-08 | +1.5678242805776501E-08 |
| `-Sene` | -0.0 | -0.0 |
| `= R_energy` | -6.5518091890538451E-05 | -9.1322262574070447E-04 |
| largest term | the lower face, -1.032681E-03 | the lower face, -9.346157E-04 |
| to close the row, the lower face would have to change by | +6.551809E-05, 6.344 % of itself | +9.132226E-04, **97.71 %** of itself |
| the same for `heat` | +6.551809E-05, 1500 % | +9.132226E-04, 2.016E+04 % |
| the same for the gravitational work | +6.551809E-05, 179.3 % | +9.132226E-04, 6362 % |

The reading: **cell 2 keeps essentially everything that enters it from
below.** The upper mass face carries 1.9 per cent of the lower one and the
upper energy face 1.2 per cent. The radiative source is 0.5 per cent of the
energy residual and the gravitational work 1.6 per cent, so neither can be
the explanation and neither is missing.

### 7.3 The molecular checkpoint, `molecular_photochem_gj1132_kzzprofile/HeH9` `g0005`

Base face `A F_m` = -4.3358126362E-06, wind-window mean +5.8873638938E-07,
base face **-7.364608** wind means (MEASURED), which is the -7.364 of
`lhs1140b_m2_conversion_audit_20260921.md` section 8. This state's wind is
not flat to the last bit like the atomic one's: gate +5.9528429054E-07 and
outer +5.8700177230E-07, 1.4 per cent apart.

The rows at the first five cells, MEASURED against the code's own scales:

| cell | r | mass `|R|/s` | momentum `|R|/s` | energy `|R|/s` |
|---|---|---|---|---|
| 1 | 1.0001933 | 5.251252E-01 | 1.449595E-03 | 1.020731E+00 |
| 2 | 1.0003866 | 1.040949E+00 | 3.172691E-04 | 9.890065E-01 |
| 3 | 1.0005800 | 2.068871E-01 | 7.537021E-07 | 1.187836E-01 |
| 4 | 1.0007733 | 1.783302E-01 | 8.251526E-07 | 1.776548E-01 |
| 5 | 1.0009666 | 1.273595E-01 | 6.847826E-08 | 1.046377E-01 |

Here the drop above cell 2 is a factor 5, not seven decades: this state's
base disturbance extends into the cells above it, which the atomic one's does
not.

| term | cell 1 | cell 2 |
|---|---|---|
| **MASS ROW** | | |
| `-A_lo F_m,lo / V` | +2.2419636463560915E-02 | +4.7193431549301383E-02 |
| `+A_hi F_m,hi / V` | -4.7211676538488637E-02 | +1.9325247158529530E-03 |
| `-S_mass` | -0.0 | -0.0 |
| `= R_mass` | -2.4792040074927725E-02 | +4.9125956265154341E-02 |
| largest term | the upper face, -4.721168E-02 | the lower face, +4.719343E-02 |
| to close the row, the lower face would have to change by | -2.479204E-02, 110.6 % of itself | -4.912596E-02, **104.1 %** of itself |
| the same for the upper face | -2.479204E-02, 52.51 % | -4.912596E-02, 2542 % |
| **ENERGY ROW** | | |
| `-A_lo F_E,lo / V` | -3.9500844199135951E-02 | +8.1253631719582331E-02 |
| `+A_hi F_E,hi / V` | -8.1285044388326289E-02 | +3.2598096465453782E-03 |
| `+grav work dF3p/V` | -8.4692702418722853E-04 | -5.5037320065529874E-04 |
| `-S_energy` | -0.0 | -0.0 |
| `-heat` | -8.8665763220114203E-04 | -9.2543601469430639E-04 |
| `+cool` | +2.8548223442427206E-06 | +2.3873147629308765E-06 |
| `-Sene` | -0.0 | -0.0 |
| `= R_energy` | -4.3514930023234463E-02 | +8.3040019465541029E-02 |
| largest term | the upper face, -8.128504E-02 | the lower face, +8.125363E-02 |
| to close the row, the lower face would have to change by | +4.351493E-02, 110.2 % of itself | -8.304002E-02, **102.2 %** of itself |
| the same for `heat` | +4.351493E-02, 4908 % | -8.304002E-02, 8973 % |
| the same for the gravitational work | +4.351493E-02, 5138 % | -8.304002E-02, 1.509E+04 % |

### 7.4 One phenomenon or two

**One.** At cell 2 of both checkpoints the mass row and the energy row are
carried by the same term, the LOWER face flux, which the upper face does not
carry on:

| | atomic `g0002` | molecular `g0005` |
|---|---|---|
| the term that carries the mass row of cell 2 | the lower face, 98.06 % of itself | the lower face, 104.1 % of itself |
| the upper mass face, as a fraction of the lower | 1.94 % | 4.09 % |
| the term that carries the energy row of cell 2 | the lower face, 97.71 % of itself | the lower face, 102.2 % of itself |
| the upper energy face, as a fraction of the lower | 1.24 % | 4.01 % |
| `heat` as a fraction of `|R_energy|` at cell 2 | 0.50 % | 1.11 % |
| the gravitational work as a fraction of `|R_energy|` at cell 2 | 1.57 % | 0.66 % |
| the sign of the base face flux | + (into the column) | - (out of it) |

The only difference is the sign: the atomic checkpoint accumulates mass and
energy in the two base cells, the molecular one drains them downward. In both,
the row is a flux difference at the base that nothing balances, the radiative
source is one to two decades too small to matter and the gravitational work
smaller still, and no source term change of any plausible size would close
either row.

At cell 1 the largest term is the lower face on the atomic checkpoint and the
upper face on the molecular one, which is the same statement from the other
side: on the molecular checkpoint the face between cells 1 and 2 carries twice
what the base face carries, so the pair drains in sequence.

## 8. The decision, in the plan's own terms

Section 17.2 gives two readings and asks which one the measurements pick.

**Reading 1, that the boundary law or its reservoir states something these
winds cannot satisfy, is not what the measurement shows.** The reservoir is
not making an unsatisfiable statement about the base face: there is a
reservoir count in the admissible family for which the base face mass flux is
exactly the wind's, and it is within a hundredth of a per cent of the atomic
count. If the reservoir statement were the defect, that point would close the
base rows. It does not.

**Reading 2 is what the measurement supports.** The rows do not close for ANY
reservoir state consistent with the prescribed base pressure and temperature.
Over a fourteen-point scan spanning a factor 1.105 in the base particle
count, a sign change of the base face flux and magnitudes of it from 2.44 to
5998 wind means, neither the mass row nor the energy row of cells 1 and 2
comes within five decades of its tolerance anywhere. The smallest values
reached anywhere in the family are 5.28E-01 (mass, cell 1, at 0.9999),
1.04E-01 (mass, cell 2, at 0.99995), 4.75E-01 (energy, cell 1, at 0.9995) and
6.27E-01 (energy, cell 2, at 1.05), at four different points, so no single
reservoir state comes close even on one pair of them. The base of these
columns is far from a fixed point that the reservoir does not determine, and
the instrument is a continuation, not a boundary repair.

**What the measurement decides beyond that.** It removes the molecular base
particle count as the cause of the M2 refusal, which section 11.3 left open:
the count is what sets the base face velocity, exactly as that audit measured,
and setting it to the atomic value returns the face velocity and the reservoir
density to the atomic control's to 3.0 and 0.03 per cent, without bringing a
single row into tolerance. And it makes the atomic and the molecular failures
one phenomenon rather than two, measured term by term at the cells where they
live.

**What it does not decide.**

- Whether a fixed point exists at this base for these columns. The scan says
  no reservoir count closes the rows of THIS state; it says nothing about
  whether some other state closes them at some reservoir count. That is the
  existence question rev9 section 16 already lists as open, and only a Mode C
  continuation with a declared budget can address it.
- The ghost COMPOSITION, the molecular base's other statement, was held
  throughout: the lower ghost carries the handoff H2 in every point of the
  scan while cell 1 has none. A probe that varies the ghost composition at
  fixed particle count is a different experiment and was not made. The M2
  audit measured that the composition step across the base face is worth 0.6
  per cent of the base face mass flux (READ,
  `lhs1140b_m2_conversion_audit_20260921.md` section 8), which is three
  decades below what the count is worth, but that is a measurement of a
  different pair of states and it is inherited, not repeated here.
- Whether the discrete base boundary law is the right law. Both probes read
  the assembly's own terms; two implementations can share one formula error,
  and neither probe is an independent statement of the physics of the
  boundary.
- The carrier and elemental rows of the molecular checkpoint, which are the
  two its certificate actually refuses beside the hydrodynamic ones, are not
  addressed. Neither is a hydrodynamic row.

## 9. Things noticed, reported and not fixed

1. **The molecular configuration's chemical heating at and above the base, on
   a column with no H2.** Section 6.4. The energy row of cells 3 to 5 is
   9.8E-01 of its own scale with the term on and 7.5E-03 with it off, and the
   scale itself falls by two decades to the atomic control's value. The
   loaded column carries zero H2 in its physical cells, so the H2 the
   collisional network acts on there is H2 the composition sweep makes
   itself; that step was NOT verified here, and the composition after the
   sweep was not written out on this route. This is a candidate defect in its
   own right, it is not the base question, and it was not pursued.
2. **`EXHALE_resolved.out` is not written on the `EXHALE_RESIDUAL` route**
   (MEASURED: `EXHALE_setup.out` is created empty and no resolved file
   appears). Rev9 section 3 names that file as one of the three a result's
   effective configuration is assembled from, so a diagnostic evaluation made
   this way cannot produce the record the plan asks for, and the derived route
   keys have to be taken from the case's own stored copy instead.
3. **The base particle count and the base level are one quantity in the
   code.** Because `n0 = p_base/(k_B T0 ntot_bc)`, a change in the base
   composition's particle count is a change in the base mass density at fixed
   pressure. That is physically right and it is worth saying plainly in the
   user manual, where "Molecular base: True" is described as an
   equation-of-state correction: it moves the base density by 5.3 per cent at
   this `q_H2`, and the base face velocity by a factor 12.
4. **`EXHALE_BASE_PARTICLE_COUNT` is not in `README.md`.** The other
   diagnostic keys are listed there, but `README.md` is being edited by
   another worker in this tree and was left alone; the key is documented in
   full at its reader.

## 10. Artifacts

Everything is in the worker scratch, `<scratch>` below:

```
<scratch> = /tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/
            59a81fec-9ef0-4992-9a57-6494f229e6bb/scratchpad/basereservoir
```

| path | what it is |
|---|---|
| `<scratch>/EXHALE_probe.x` | the private build, md5 `39604c4237fe13eb9818676e2a977f5a` |
| `<scratch>/build.log`, `build2.log`, `build3.log` | its builds; `build3.log` is the rebuild from an empty object directory that gives the same md5 |
| `<scratch>/runs/<run>/` | the six evaluations, each with `input.inp`, `output/`, `run.log` and `SIDECAR.txt` |
| `<scratch>/scan/n<value>/` | the fourteen points of the reservoir-count scan |
| `<scratch>/base_rows.py` | the consumer of section 5 |
| `<scratch>/probeA_rows.txt`, `probeB_rows.txt` | its output for the six evaluations |
| `<scratch>/scan_table.txt` | the table of section 6.3 |
| `<scratch>/reg/regression/` | the scratch copy of the harness |
| `<scratch>/reg/check_unset.log` | the regression of section 4.3 |

Nothing was written into `LHS1140b/models/`, nothing into
`backup/regression/`, no golden was refreshed, and the only file added to the
repository is this one. The one source change is the diagnostic of section 4,
in `src/modules/files_IO/input_read.f90`, uncommitted.
