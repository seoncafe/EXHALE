# Stage F: the continuation in XUV of the LHS 1140 b prescribed-composition branch

`PLAN_20260920_rev9.md` section 14, authorized by the user on 2026-09-21.
Mode C throughout (section 3: a retained trajectory on the current branch,
with a declared budget and a recorded ending).

Every number below is labeled MEASURED (this work ran it), READ (from a file
of this tree) or INSPECTED (source, with `file:line`).

## What it found, in five lines

- **He/H = 9.7** (`atomic_scalar_gj1132x<XUV>_kzz1e9/HeH9.7`): the continuation from the
  certified x0.10 state reaches **0.07 and stops there**. Three step sizes
  below it (0.146, 0.067 and 0.032 dex out of the same certified rung) refuse
  in the same shape, so halving the step is not the instrument and a turning
  point is suspected; this work stops at it, as section 14 directs.
- **He/H = 2.13** (the same group, `HeH2.13`): the SAME ladder, same
  binary, grid, base law, flux, spectra and controls, **descends the whole
  decade to 0.01 with all six rungs certified**, none refused, no step halved
  and the under-relaxation never shortened.
- What refuses, whenever anything does, is one equation: the **elemental
  transport He/H partition row in the wind**. The three hydrodynamic rows are
  within tolerance at every one of the 158 outer passes of this work (49 on
  the He/H = 9.7 ladder, 109 on the He/H = 2.13 one).
- So the obstruction is not the spectrum, the grid, the boundary law, the
  flux, the pseudo-time start or the budget, all of which the two ladders
  share. It follows the composition, and the evidence points at the
  alternation rather than at the absence of a solution (section 5, reading 2).
- The result language is section 7, and it is the plan's: a lowest CERTIFIED
  XUV on a tested branch with this grid, boundary law and numerical method,
  never a statement that no stationary atmosphere exists below it.

---

## 1. Identity

### 1.1 The executable

| field | value | source |
|---|---|---|
| binary | `/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/EXHALE.x` | the tree's binary of record for 2026-09-21 |
| md5 | `89ec67149aa452eea3deeb19052f83d1` | MEASURED, `md5sum` |
| compiler | conda-forge gcc 16.2.0 (the `GCC:` notes in the image also carry the system 8.5.0 of the C runtime) | MEASURED, `strings EXHALE.x \| grep GCC:` |
| OpenMP runtime | `/opt/miniconda3/lib/libgomp.so.1` | MEASURED, `ldd` |
| BLAS/LAPACK | `/opt/miniconda3/lib/libopenblas.so.0` | MEASURED, `ldd` |
| other linked | `libgfortran.so.5`, `libquadmath.so.0`, same prefix | MEASURED, `ldd` |
| git | `93eed86`, working tree dirty (67 entries) | MEASURED, `git rev-parse --short HEAD`, `git status --porcelain \| wc -l` |
| host, threads | `lart4` (72 cores), `OMP_NUM_THREADS=8` | MEASURED |

No build was performed by this work.

### 1.2 The effective parsed configuration of every rung

READ from `EXHALE_resolved.out` of the rungs (identical in all of them):

```
well_balanced             T
carrier_transport         F
carrier_in_newton         F
carrier_newton_on_stall   F
ionization_transport      F
oxygen_chemistry          F
```

READ from `EXHALE_setup.out` (the two the resolved summary does not carry):
`Numerical flux: ROE`, `Reconstruction method: PLM` (the stationary route
reconstructs with WENO3; the certification block of every rung states
`operator: WENO3, numerical flux: ROE, well balanced: T`).
`EXHALE_RESID_QUAD` was NOT set by this work, in any rung.

`EXHALE_CARRIER_DEBUG=1` was in force in every rung of this work. It writes
nothing here, and the reason is recorded rather than passed over: the case is
atomic, `carrier_transport` is `F`, and every pass prints `carrier residual of
the returned composition 0.00E+00 ... [the run carries no relaxed carrier]`
(READ, `run.log`). **The ledger of the passes, for this family, is the
`(EXHALE_main) outer pass N:` table of `run.log`**, which is what section 3
of this document tabulates.

### 1.3 The seed of the He/H = 9.7 ladder

(the second ladder's seed is in section 6.)

| field | value | source |
|---|---|---|
| case | `LHS1140b/models/atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7` | the brief |
| generation | `g0004_20260919T212852Z_de543034` | READ, `state_index.json`: `latest_certified` AND `latest_complete` |
| certification | `CERTIFIED`, `state_claim_certified: true` | READ, the same index entry |
| `Hydro_ioniz.txt` md5 | `ee5bd9ed574518dce590c4cd626f682c` | MEASURED, equal to the index's field |
| `Ion_species.txt` md5 | `4b20e8f310ab268f18087816ee01b28c` | MEASURED, equal to the index's field |
| `certification.txt` md5 | `a434c5ab48d14f1dab846540d93b333d` | MEASURED |
| boundary model | `characteristic_face_ps_reservoir_C_minus_contact_upwind_ghost_fixed_point_seed_reservoir_row_v3` | READ, the state header |
| reservoir | `He/H 9.6999999999999993E+00` | READ, the state header |
| grid | `N 500 R0[cm] 1.1273716464000001E+09 r_min[Rp] 1.0001933187778382E+00 r_max[Rp] 2.9031224061195065E+01 mode Mixed`, bitwise the `# grid` line of `models/current_grid_Hydro_ioniz.txt` | MEASURED, both lines compared |
| the binary that wrote it | `75d55d9d4fd0e748cd01d6e35713e34c` | READ, the index entry |

The generation the ladder starts from was written by a binary that is not the
one this work runs. That is the ordinary situation of a continuation: the seed
is a predictor, and what the rung certifies is measured by the current binary
on the state the current binary writes.

### 1.4 The case the rungs are built from

Every rung directory holds the input of
`LHS1140b/models/atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` with **one line
changed**, `Spectrum file:`, made absolute. That is the rule of the earlier
continuation (`models/.L25/ladder/run_ladder.sh`), and it is what keeps the
grid, the base law, the flux, the tolerance and the composition fixed along
the ladder.

MEASURED, `diff` of the two catalog inputs `x0.10` and `x0.01`: they differ in
three lines only, the leading comment, `Planet name:` and `Spectrum file:`.
So building the rungs from the target input or from the seed's input is the
same physics; the rungs inherit the target's `Planet name` label, which is
inert (INSPECTED: `p_name` is read at `input_read.f90:240` and used only by
`write_setup_report.f90:68,704`).

---

## 2. The spectrum scaling, as a checked fact

Section 14 requires the band definition to be fixed along the ladder. It is,
and this is what was checked rather than assumed.

**Where the files come from.** `LHS1140b/sed/make_sed.py` builds the two
reference SEDs at the planet orbit (`lhs1140_sed_gj1132_at_b.txt` and the
GJ 699 alternate) and writes NO scaled file: INSPECTED, the script's only
outputs are `lhs1140_sed_{tag}_at_b.txt`. The scaled grid was written by the
rule stated in its own file header and in the update log
(`docs/Update_EXHALE_stage2.md` line 9927, `docs/lhs1140b_stationary_L14_20260914.md`
section 3.3): **the whole flux column of the reference file times the factor,
on the same wavelength grid.** The header line of each scaled file says it,
e.g. `# flux column scaled by 0.01 relative to lhs1140_sed_gj1132_at_b.txt`
(READ).

**The rule is now a script, and it reproduces the whole grid.** The scaling
had never been written down as anything but prose, so this work put it in
`LHS1140b/sed/scale_sed.py` beside `make_sed.py` and MEASURED that it
reproduces **every scaled spectrum on the tree, byte for byte, all 30003
lines each** (`scale_sed.py --check`: 0.33, 0.30, 0.28, 0.26, 0.25, 0.20,
0.18, 0.16, 0.15, 0.10, 0.07, 0.05, 0.03, 0.02, 0.015, 0.01, zero differing
lines in each, and the two this work added afterwards are in the check list
too). The two spectra the halved steps of section 4.3 needed were written
with it. The script refuses to overwrite an existing spectrum.

**The check.** MEASURED here (`audit/measure_sed.py`, `audit/sed_bands.txt`):

| factor | file md5 | rows | same wavelengths | F_XUV(10-1300 A) | F_X(0.25-2 keV) | F_EUV(100-911 A) | F(911-2583 A) | F_XUV / reference | max rel. deviation from column x factor |
|---|---|---|---|---|---|---|---|---|---|
| 1 (reference) | `ac73a75c6396793435ed72e1cf34c427` | 29992 | - | 3.344870E+01 | 4.407194E+00 | 1.005124E+01 | 2.634986E+01 | 1.000000 | - |
| 0.10 | `3c46ba24ea8165b5071cee50e4deac3d` | 29992 | yes | 3.344870E+00 | 4.407194E-01 | 1.005124E+00 | 2.634986E+00 | 0.100000 | 2.2E-16 |
| 0.07 | `43751af0e1f45cfa968e54fe549cba6a` | 29992 | yes | 2.341409E+00 | 3.085036E-01 | 7.035871E-01 | 1.844490E+00 | 0.070000 | 5.0E-07 |
| 0.05 | `bc2090e7dd7324cfdf8bd7a10d5e22d7` | 29992 | yes | 1.672435E+00 | 2.203597E-01 | 5.025623E-01 | 1.317493E+00 | 0.050000 | 5.0E-07 |
| 0.03 | `0860cfad7704c37b984238c47c02d7e6` | 29992 | yes | 1.003461E+00 | 1.322158E-01 | 3.015373E-01 | 7.904957E-01 | 0.030000 | 5.0E-07 |
| 0.02 | `219ab99e40dcf7d1a51ce45b49bb5549` | 29992 | yes | 6.689739E-01 | 8.814388E-02 | 2.010249E-01 | 5.269971E-01 | 0.020000 | 4.0E-07 |
| 0.015 | `a082358b94d1f420f00d5ef81bc25be1` | 29992 | yes | 5.017304E-01 | 6.610790E-02 | 1.507687E-01 | 3.952478E-01 | 0.015000 | 5.0E-07 |
| 0.01 | `7da5fa41a0706feacac74b7d245fbec3` | 29992 | yes | 3.344870E-01 | 4.407194E-02 | 1.005124E-01 | 2.634986E-01 | 0.010000 | 2.2E-16 |

Band fluxes in erg cm^-2 s^-1, trapezoidal on the file's own abscissa.

Three statements follow, and they are the ones section 14 needs:

1. **No band is redefined.** Every file carries the same 29992 wavelengths as
   the reference (MEASURED, exact array equality), so the boundaries 10, 100,
   911, 1300 A and the X-ray band fall on the same bins in every rung.
2. **The scaling is achromatic.** Every band ratio equals the nominal factor
   to six digits, and the whole column equals the reference column times the
   factor to 5E-07 (the six-digit `%.6e` text rounding of the written file) or
   to machine precision where the factor is a power of ten.
3. **The spectrum sets the field, and nothing else in the input does.**
   INSPECTED, `sed_read.f90:472-496`: when the spectrum is loaded from a file,
   `LX` and `LEUV` are OVERWRITTEN by the file's own band integrals, so the
   `Log10 of X-ray luminosity:` and `Log10 of EUV luminosity:` lines of
   `input.inp` (identical at every rung, and stale) do not enter. The XUV of a
   rung is the file's.

The ladder is therefore a one-parameter family in a single achromatic scale
factor on the irradiation, and the ratios quoted in the tables below are
MEASURED integrals, not nominal labels.

---

## 3. How a rung is run, and what its budget is

A rung is one directory holding the target case's `input.inp` with the
spectrum line replaced, run by the tree's own `LHS1140b/models/run_case.sh`,
so the four passes, the route, the restart intent, the seed mapping and the
publishing contract are the catalog's and not this work's. The driver is
`ladder.sh` in the scratch directory named in section 9.

Nothing was written under `LHS1140b/models/`. The runner was given the rung
directory as a path out of `models/` into the scratch tree, so every file it
wrote (`output/`, `states/`, `state_index.json`, `runs/`, `REPRODUCE.md`,
`ENDING`) is inside that rung directory. Publication into the catalog is a
separate step and is not part of this work.

**The controls, fixed for every rung:**

| control | value | why this value |
|---|---|---|
| `OMP_NUM_THREADS` | 8 | the catalog's, and the brief's |
| `EXHALE_OUTER_PASSES` | 40 | the budget of class `atomic_prescribed` (READ, `models/run_campaign.sh:175`), which is also `run_case.sh`'s own default |
| JFNK iteration cap | 3000 per solve, from the caller | READ, `run.log`: `(JFNK) outer iteration cap 3000, from the caller` |
| `EXHALE_PTC_DTAU0` | 1.0e8 | section 3.1 |
| `SEED_ATTEMPTS` | 1 | the seed IS the continuation predictor; no walk down `pick_seed.py`'s list is admissible in a continuation |
| `FORCE` | 1 | re-run a rung directory that already carries products |
| `EXHALE_CARRIER_DEBUG` | 1 | asked for by the brief; section 1.2 records that it writes nothing in an atomic case |
| wall | 60 min per rung, SAFETY ceiling only | the catalog campaign's 30 min is a BUDGET it stops cases at, and a budget in wall clock would make the ending depend on how many other workers share the host; here the budget is the 40 passes, and the wall exists only so that a run cannot hang forever. An external termination is recorded as one and never converted into a measured ending |

**The budget is 40 outer passes and 3000 Newton iterations per solve.** The
wall ceiling is not a budget, and no rung of either ladder reached it.

What the whole item spent, MEASURED: **158 outer passes and 946 JFNK
iterations** over the ten rungs, plus 50 JFNK iterations in the pseudo-time
control experiment of section 3.1, at `OMP_NUM_THREADS=8`; 5489 s of runner
wall clock, of which 2761 s is the one publication stall of section 4.5. The
outer-pass budget was reached by no rung: the most any rung spent is 34 of
40 (rung 1 of the He/H = 9.7 ladder), and the He/H = 2.13 rungs spent 17 to
19 each. No single solve came near the 3000-iteration Newton cap; the whole
34-pass rung 1, every solve of it together, took 220 iterations.

### 3.1 The pseudo-time start, and the one control experiment on it

The catalog recipe for a prescribed-composition case is
`EXHALE_PTC_DTAU0=1.0` (READ, `MODELS.md` section 6, the four-pass table).
`run_case.sh`'s own header states the exception, and it is this family at a
neighbouring composition: "On the He/H = 2.13 column of LHS 1140 b below 0.10
of the GJ 1132 spectrum the dtau0 = 1 solve falls into the ramp collapse of
plan item L4h ... Where that is the case, state EXHALE_PTC_DTAU0=1.0e8 in the
environment and the first solve is taken there" (READ,
`run_case.sh:1105-1112`). The header names the He/H = 2.13 column; the case
here is He/H = 9.7, and the control experiment below MEASURES the same
collapse on it, so the exception is established for this column and not
carried over by assumption. The continuation of 2026-09-17 that produced the
currently certified states of this family used 1.0e8 at every rung (READ,
`models/.L25/ladder/run_ladder.sh`).

This work MEASURED the alternative before fixing the control, because a
continuation whose control is known to fail at its first step measures the
control and not the branch.

**The control experiment** (scratch `probe_dtau0_default/x0p07_dtau0_1p0/`):
the x0.07 rung, the same seed, the same input, `EXHALE_PTC_DTAU0=1.0`. The
first stationary solve enters the ramp collapse immediately and never leaves
it: from `||R|| = 3.509E-01` the step length stays at `lam = 9.54E-07` while
`dtau` walks down 2.50E-01, 6.25E-02, ..., 8.03E-05 over seven iterations with
`||R||` unmoved to four digits and 500 of 500 energy cells outside tolerance;
the damped Gauss-Newton escape at iteration 8 raises `||R||` to 1.073E+00, and
by iteration 41 it stands at 1.929E+00 with all 500 mass cells and 466
momentum cells outside tolerance and the worst row at cell 1. MEASURED, the
archived `run.log`. **The run was then stopped by me** and the ending is
recorded as what it is: `stopped_by_signal ... 331 s after it started, not at
a wall ceiling the runner was told of, in the wind pass; that pass had written
no complete state of its own, so nothing new is published` (READ, its
`ENDING`). It is not an ending the solver reached, and no claim is made here
that the solve would never have recovered.

So the ladder's fixed control is `EXHALE_PTC_DTAU0=1.0e8`, and this is part of
"the numerical method" the closing statement of section 7 is conditioned on.

### 3.2 The control the catalog itself provides

The most recent catalog attempt at the target case is the SAME seed taken in
ONE step: `atomic_scalar_gj1132x0.01_kzz1e9/HeH9.7` seeded from
`g0004_20260919T212852Z_de543034` of the x0.10 case (READ, its
`not_solved.md`), `EXHALE_PTC_DTAU0=1.0`, binary `75d55d9d`. It completed 4 of
its 40 outer passes in 30 minutes and was stopped by the campaign's wall
ceiling, refusing on the hydrodynamic energy row 9.978E-01 at cell 248 and the
He/H partition row 9.768E-04 at cell 354. That is the one-step control against
which the ladder below is the graded descent: same seed, same endpoint,
different number of steps and a different pseudo-time start.

## 4. The rungs of the He/H = 9.7 ladder

The family is `atomic_scalar_gj1132x<XUV>_kzz1e9/HeH9.7`: H, He, He(2^3S) and
electrons above the prescribed scalar base, `He_diffusion: True` with
`He_Kzz: 1.0e9`, `Numerical flux: ROE`, `Well balanced: True`, 500 cells,
reservoir He/H = 9.7.

### 4.1 The ladder table

Every row is MEASURED by this work. `dlog10` is the descent step from the rung
above. F_XUV is the MEASURED 10-1300 A integral of the rung's own spectrum.
"gated row" is the elemental transport He/H partition row, which is the entry
that decides every rung of this ladder.

| # | XUV | F_XUV | dlog10 from | seeded from | passes used | ending | info | gated row at the end | cell (r/R_p) | wall |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.07 | 2.341 | 0.1549 (0.10) | `g0004_20260919T212852Z_de543034` of the catalog x0.10 case | 34 of 40 | **solved, CERTIFIED** | 0 | 9.924E-06 within 1.0E-05 | 351 (3.087) | 578 s |
| 2a | 0.05 | 1.672 | 0.1461 (0.07) | rung 1, `g0002_20260921T095220Z_fbeeef3d` | 5 of 40 | composition_refusal | 1 | 4.100E-02 above 1.0E-05 | 352 (3.124) | 2903 s |
| 2b | 0.06 | 2.007 | 0.0669 (0.07) | rung 1, the same state | 5 of 40 | composition_refusal | 1 | 3.231E-02 above 1.0E-05 | 350 (3.051) | 144 s |
| 2c | 0.065 | 2.174 | 0.0322 (0.07) | rung 1, the same state | 5 of 40 | composition_refusal | 1 | 2.031E-02 above 1.0E-05 | 350 (3.051) | 137 s |

**The lowest certified XUV of this ladder is 0.07 of the fiducial GJ 1132
spectrum**, F_XUV(10-1300 A) = 2.341 erg cm^-2 s^-1 at the orbit, log10 Mdot
= 6.64 g/s (READ, `pp.log` of rung 1: `Log10 of steady-state Mdot = 6.64`).
The smallest step that succeeded is the one that reached it, 0.1549 dex; the
LARGEST step that failed is 0.1461 dex and the SMALLEST step that failed is
0.0322 dex, all three from the same certified rung 1 state.

### 4.2 Rung 1, the one that certified

Seeded from the catalog x0.10 certified generation, mapped onto
`models/current_grid_Hydro_ioniz.txt` (an identity of the grid line, section
1.3) with no reservoir rescale, the composition being the same.

The alternation: the gated elemental He/H row rises for two passes, the
progress control halves `omega` once at pass 3, the row turns at pass 4 and
falls monotonically to acceptance at pass 34.

| pass | gated row | omega | hydro info |
|---|---|---|---|
| 1 | 4.77E-03 | 0.500 | 0 |
| 2 | 8.86E-03 | 0.500 | 0 |
| 3 | 9.25E-03 | 0.250 | 0 |
| 4 | 8.67E-03 | 0.250 | 0 |
| 8 | 5.24E-03 | 0.250 | 0 |
| 16 | 8.90E-04 | 0.250 | 0 |
| 24 | 9.95E-05 | 0.250 | 0 |
| 33 | 1.27E-05 | 0.250 | 0 |
| 34 | 9.92E-06 | 0.250 | 0, ACCEPTED |

The full table is `audit/rungs/x0p07.md`. Every pass returned `hydro info=0`:
**the hydrodynamics was never the obstruction at any rung of this ladder.**

The certificate of the published generation `g0002_20260921T095220Z_fbeeef3d`,
every active equation (READ, its `certification.txt`; 6 of the 30 in the
inventory are active in this configuration):

```
   hydrodynamic mass row          max= 4.686E-08  cell=2    tol= 1.4E-07  within
   hydrodynamic momentum row      max= 3.016E-14  cell=492  tol= 1.0E-08  within
   hydrodynamic energy row        max= 9.931E-08  cell=2    tol= 1.0E-06  within
   elemental transport He/H       max= 9.924E-06  cell=351  tol= 1.0E-05  within
   level balance He 2^3S          max= 2.993E-20  cell=481  tol= 1.0E-06  within
   eliminated-species closure     max= 2.974E-17  cell=486  tol= 1.0E-06  within
   CERTIFIED: every active equation was evaluated and is within its tolerance
```

The elemental row certifies at 0.992 of its tolerance. The rung is certified
and it is certified with no margin on the one row that decides it.

The four catalog passes all ran, so the rung also carries its
advection-corrected profiles and its line: log10 Mdot = 6.64 g/s, and the
He I 10830 red pair through the WINERED HIRES-Y kernel at EW = 0.0445 %A,
red depth 0.1406 per cent (READ, the runner's own `DONE` line). For scale,
the catalog's certified x0.10 case of this same series stands at log10 Mdot
6.81 and EW 0.2353 (READ, `MODELS.md` section 7), so a factor 1.43 down in
XUV takes the line down by a factor 5.3.

### 4.3 The three refusals below it, and why the step was halved twice

The step 0.07 -> 0.05 was refused. Following the procedure, no control was
changed and the step was halved: 0.07 -> 0.06 (0.458 of the refused step).
That was refused in the same way, so it was halved again: 0.07 -> 0.065
(0.48 of the previous, 0.22 of the original). Refused again. The two
spectra this needed did not exist and were written by the rule of section 2
(`scale_sed.py 0.06 0p06 "..."`, `scale_sed.py 0.065 0p065 "..."`), their md5
`4279de66b46464f6caaadc925db10740` and `036e3b99cf82bf042059dde893e6cdac`.

**The three refusals are one phenomenon.** In all three the wind is solved and
the composition is not:

| pass | 0.05 | 0.06 | 0.065 | omega | what the control said |
|---|---|---|---|---|---|
| 1 | 1.25E-02 | 5.83E-03 | 2.81E-03 | 0.500 | the coupled residual fell |
| 2 | 3.14E-02 | 2.25E-02 | 1.40E-02 | 0.500 | the composition residual, not yet measured twice |
| 3 | 3.88E-02 | 2.98E-02 | 1.89E-02 | 0.250 | neither the coupled nor the composition residual fell |
| 4 | 4.04E-02 | 3.17E-02 | 2.00E-02 | 0.125 | neither fell |
| 5 | 4.10E-02 | 3.23E-02 | 2.03E-02 | 0.125 | neither fell, third strike, REFUSED |

(the gated elemental He/H row; MEASURED, the three `run.log`s. `hydro info=0`
at every one of these fifteen passes.)

The ending is the outer loop's own progress control, not the budget: 5 of the
40 passes were spent. INSPECTED, `EXHALE_main.f90:6885`,
`outer_no_fall_max = 3`, and the rule at `:6804-6819`: passes are counted, and
`omega` halved to its floor of 0.125, only where the gated row AND the
composition's own distance to its fixed point both failed to fall. In these
three rungs both stand still from pass 3 on.

The certificates of the three refused states are identical in shape (READ,
their `certification.txt`): the three hydrodynamic rows within tolerance, the
He 2^3S level balance and the eliminated-species closure within by thirteen
decades, and the elemental transport He/H partition alone above, by a factor
2.0E+03 (0.065), 3.2E+03 (0.06) and 4.1E+03 (0.05), at cell 350 to 352,
r = 3.05 to 3.12 R_p.

### 4.4 The branch marker: the response of the partition row per unit step

The gated row after the FIRST pass, divided by the step in dex, is a measure
of how far the He/H partition's fixed point moves per unit log XUV. MEASURED:

| step | dlog10 | gated row after pass 1 | per dex |
|---|---|---|---|
| 0.10 -> 0.07 | 0.1549 | 4.770E-03 | 3.079E-02 |
| 0.07 -> 0.065 | 0.0322 | 2.810E-03 | 8.731E-02 |
| 0.07 -> 0.06 | 0.0669 | 5.830E-03 | 8.708E-02 |
| 0.07 -> 0.05 | 0.1461 | 1.250E-02 | 8.554E-02 |

Two things follow, and they are the quantitative content of this ladder:

1. **The response is linear in the step.** The three steps out of rung 1 give
   8.73, 8.71 and 8.55 E-02 per dex over a factor 4.5 in step size. The
   predictor is first-order accurate and the refusal is NOT a step-size
   effect: a step four and a half times smaller leaves the same residual per
   unit step and ends the same way.
2. **The response nearly triples between 0.10 and 0.07**, from 3.08E-02 to
   8.7E-02 per dex. The partition the wind holds at r about 3 R_p is a
   rapidly steepening function of the irradiation as the XUV falls through
   this decade. That, and not the arithmetic of the step, is what the
   continuation ran into.

### 4.5 One thing the runner did that is not a result

The 0.05 rung's wall clock, 2903 s, is not solver time. MEASURED from the file
timestamps: the binary's solve ended at 18:55:02 (`run.log` and
`output/Hydro_ioniz.txt`) and `run_case.sh` returned at 19:41:03, with
`states/`, `state_index.json`, `not_solved.md` and `REPRODUCE.md` all stamped
19:41:03. So 46 minutes were spent inside the publication of the refused
generation, on a host carrying four other workers. The same publication took
under two seconds on the rung before it and the two after it. The cause was
NOT established here and it did not touch the physics: the state, the log and
the certificate are those the binary wrote at 18:55:02. It is recorded because
a wall clock in a budget table has to mean what it says.

## 5. The path of the He/H = 9.7 ladder

What that continuation traced, in the terms section 14 asks for. The second
composition, which behaves differently and decides how this section is to be
read, is section 6.

**Which row refuses first as the XUV falls, and where.** One row, at every
rung: the **elemental transport He/H partition** in the wind, at cell 350 to
352, r = 3.05 to 3.12 R_p, gated at r >= 1.2 R_p against 1.0E-05. The three
hydrodynamic rows never refuse: `hydro info=0` at all 49 outer passes of all
four rungs, and in the three refused states the mass and energy rows stand
inside their tolerances by factors of 3 to 6 and the momentum row by six
decades. **The obstruction of this
continuation is the composition, in the wind, and it is the same obstruction
the plan's phase 9 found on the x0.01 atomic checkpoint from its own stored
state** (section 17.2: "the elemental transport He/H partition row in the
wind, 7.218E-05 against 1.0E-05 at cell 295"), and the same one phase 6 found
on the molecular side. Three different states, one row.

**The smallest step that succeeded and the largest that failed.** Succeeded:
0.1549 dex, 0.10 -> 0.07, and it is also the only step that succeeded.
Failed: 0.1461, 0.0669 and 0.0322 dex, all out of the certified 0.07. The
largest failure is smaller than the success. **The ladder is not ordered by
step size**, and that is the central negative result: between 0.10 and 0.07 a
step of 0.155 dex is taken, and out of 0.07 a step of 0.032 dex is not.

**The residual structure at each ending.** Identical in the three refusals
(section 4.3): hydrodynamics solved, the He 2^3S level balance and the
eliminated-species closure satisfied by eleven to fourteen decades, the
partition row alone above tolerance by 2E+03 to
4E+03, and the row RISING over the five passes rather than stalling at a
level. The refusing cells are the same three cells in all three.

**Why each correction ended.**

- rung 1 (0.07): acceptance. Every active equation within tolerance at outer
  pass 34 of 40, `info = 0`, `||R|| = 9.931E-08`.
- rungs 2a, 2b, 2c (0.05, 0.06, 0.065): the outer loop's progress control.
  Three consecutive passes in which neither the joint distance nor the
  composition's own distance to its fixed point fell, with `omega` already at
  its floor of 0.125 (INSPECTED, `EXHALE_main.f90:6885` and the rule at
  `:6804-6819`). 5 of the 40 passes were spent. Neither the pass budget nor
  the Newton iteration cap nor the wall ceiling was reached in any rung.
- the control experiment on the pseudo-time start (section 3.1): stopped by
  me, recorded as `stopped_by_signal`.

**A turning point is suspected here, and this work STOPS at it.** The
condition of section 14 is met: no step size below the last succeeds
(0.1461, 0.0669, 0.0322 dex, a factor 4.5 apart, all refused in the same
shape and at the same cells) while the rung above still certifies. Section 14
says an augmented continuation equation is the instrument for that, and this
item does not carry one, so nothing further was tried and no control was
changed to force a descent.

**Section 6 then measured the same ladder at He/H = 2.13 and it descends the
whole decade**, which does not resolve the suspicion but bounds it: whatever
stops the helium-rich column is not shared with the helium-poor one, and
every part of the method is.

**What the suspicion is a suspicion OF, stated apart from the evidence.** The
measurements above are consistent with at least three different things, and
this ladder does not separate them:

1. a fold of the certified branch in this parameter, somewhere below 0.07;
2. a fixed point that still exists but that the ALTERNATION cannot reach,
   the wind solve and the element relaxation each converging while the
   composite map does not contract. The progress control's own text says this
   is what it measures, not the absence of a solution, and the control ends
   the loop after 5 of 40 passes at `omega = 0.125`;
3. a discretization limit of this 500-cell grid at r about 3 R_p, where the
   refusing cells sit and where the element flux is the quantity the row
   measures.

Reading 2 has the most direct support, because the response measure of
section 4.4 is clean and linear while the outcome is not: a map whose
first-order response is well behaved and whose iteration nevertheless stops
moving is a statement about the iteration. It is NOT established here.

## 6. The second composition, He/H = 2.13

The same ladder on `atomic_scalar_gj1132x<XUV>_kzz1e9/HeH2.13`, everything
fixed except the reservoir: the input of the catalog x0.01 case of that
composition, with only `Spectrum file:` changed, seeded from its certified
x0.10 generation `g0005_20260919T212847Z_80bee6f1` (READ, that case's
`state_index.json`; MEASURED md5 `53e80a06549d375559c319bacd2daa05` and
`185040241b2bb3bf7fde80f7cc2ca25c`, equal to the index's fields). Same
controls as section 3.

**One difference has to be stated.** The catalog's x0.10 case of this
composition is solved with `Numerical flux: HLLC` and its x0.01 case states
`ROE` (READ, the two `input.inp`), so rung 1 of this ladder takes an
HLLC-solved state as its predictor and solves with ROE, as the catalog's own
recipe for the low-XUV members of this composition requires. Along the
ladder, from rung 1 on, only the spectrum moves.

### 6.1 The rungs

**All six rungs certified. The ladder reached 0.01 of the fiducial spectrum,
the bottom of the existing grid, without one refusal and without one halved
step.** Every ending is `solved`, `info = 0`, and the certificate is the
current binary's on the bytes the current binary wrote.

| # | XUV | F_XUV | dlog10 | seeded from | passes | ending | \|\|R\|\| | gated row at acceptance | cell (r/R_p) | log10 Mdot | wall |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.07 | 2.341 | 0.1549 | the catalog x0.10 case, `g0005_20260919T212847Z_80bee6f1` | 19 | solved | 1.976E-07 | 6.414E-06 | 254 (1.384) | 6.66 | 292 s |
| 2 | 0.05 | 1.672 | 0.1461 | rung 1, `g0002_20260921T105424Z_44ee896b` | 18 | solved | 2.969E-07 | 7.732E-06 | 250 (1.358) | 6.51 | 247 s |
| 3 | 0.03 | 1.003 | 0.2218 | rung 2, `g0002_20260921T105830Z_2d661cc0` | 19 | solved | 1.169E-07 | 6.723E-06 | 247 (1.340) | 6.28 | 239 s |
| 4 | 0.02 | 0.6690 | 0.1761 | rung 3, `g0002_20260921T110229Z_434efc7f` | 18 | solved | 1.096E-07 | 8.226E-06 | 245 (1.328) | 6.10 | 266 s |
| 5 | 0.015 | 0.5017 | 0.1249 | rung 4, `g0002_20260921T110656Z_ba1a27aa` | 17 | solved | 9.233E-07 | 9.233E-06 | 243 (1.316) | 5.97 | 173 s |
| 6 | 0.01 | 0.3345 | 0.1761 | rung 5, `g0002_20260921T110949Z_6252a706` | 18 | solved | 3.154E-07 | 8.011E-06 | 244 (1.322) | 5.79 | 178 s |

The published generation of rung 6 is
`g0004_20260921T111247Z_0aebdf49` (corrected 2026-09-22: this line first read
`g0002`, the name the generation carried while the continuation ran. Rungs 1
to 5 opened new cases and kept their `g0002`; rung 6 was published into
`atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13`, which already held `g0001` to
`g0003`, so it was renumbered on publication. Timestamp and hash are
unchanged, and `state_index.json` of that case gives
`latest_certified = g0004_20260921T111247Z_0aebdf49`, MEASURED).
Tolerances as in section 4.2; every rung's
certificate has the same six active equations all within, and the elemental
row lands at 0.64 to 0.92 of its gate, which is where the alternation stops
by construction.

**A cross-check that was not arranged.** The catalog's own
`atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` entry, which is an UNCERTIFIED
evaluation of a state a 30-minute wall ceiling interrupted, records log10 Mdot
5.790 (READ, `MODELS.md` section 7). Rung 6 here, reached by six certified
steps from a different state and certified, gives 5.79. The two agree to the
printed digits.

### 6.2 The contrast with He/H = 9.7, which is the result of this section

| | He/H = 2.13 | He/H = 9.7 |
|---|---|---|
| lowest certified XUV | **0.01**, the bottom of the grid | **0.07** |
| steps taken | 6, all certified | 1 certified, 3 refused |
| step halvings needed | none | two, and both refused |
| `omega` | 0.500 at every pass of every rung, never halved | halved once on the certified rung, twice on each refused one, to its floor |
| gated-row trajectory | falls monotonically from pass 1, by a factor 0.572 to 0.809 per pass on rung 1, the later passes settling at 0.572 | rises for two to five passes; on the certified rung it turns at pass 4, on the refused ones it never turns |
| passes per rung | 17 to 19 | 34 on the certified rung, 5 before the control stops the refused ones |
| response of the row per dex of step | 0.115 at 0.10, rising to 0.188 at 0.015 | 0.031 at 0.10, 0.087 at 0.07 |

The response measure of section 4.4 on this composition, MEASURED the same
way (the gated row after pass 1, divided by the step in dex):

| step | dlog10 | gated row after pass 1 | per dex | passes | smallest omega |
|---|---|---|---|---|---|
| 0.10 -> 0.07 | 0.1549 | 1.780E-02 | 1.149E-01 | 19 | 0.500 |
| 0.07 -> 0.05 | 0.1461 | 1.960E-02 | 1.341E-01 | 18 | 0.500 |
| 0.05 -> 0.03 | 0.2218 | 3.150E-02 | 1.420E-01 | 19 | 0.500 |
| 0.03 -> 0.02 | 0.1761 | 2.960E-02 | 1.681E-01 | 18 | 0.500 |
| 0.02 -> 0.015 | 0.1249 | 2.340E-02 | 1.873E-01 | 17 | 0.500 |
| 0.015 -> 0.01 | 0.1761 | 3.310E-02 | 1.880E-01 | 18 | 0.500 |

**The column with the LARGER response is the one that descends.** The
helium-poor column's partition row moves 0.115 to 0.188 per dex of log XUV,
1.6 times the helium-rich column's at the one step both took out of 0.07 and
3.7 times it at the step out of 0.10, and it converges at every rung with
no under-relaxation halving at all; the helium-rich column moves less per dex
and stops. So the size of the residual a step produces is NOT what decides
whether the rung certifies, and neither is the step: what decides is whether
the alternation contracts, which is reading 2 of section 5.

The ladder also says what the obstruction is not. It is not the spectrum
scaling, the grid, the base law, the boundary model, the flux, the pseudo-time
start, the pass budget or the binary, because all of those are shared by the
two ladders. It is the composition: at He/H = 9.7 the helium column the wind
must carry through r about 3 R_p, where the refusing cells sit, is four and a
half times the one at He/H = 2.13, and it is there that the element relaxation
and the wind solve stop agreeing.

### 6.3 What these six states are

Six certified states of this composition at six XUV normalizations, produced
by the current binary on the current branch, sitting in scratch. Only one of
the six, `0.01`, is a case of the catalog (the catalog's XUV grid for this
composition is 0.01, 0.10, 0.15, 0.20, 0.25, 0.30, 0.33); the other five are
the continuation's own intermediate points and belong to no catalog case.
**The 0.01 one matters**: the catalog's
`atomic_scalar_gj1132x0.01_kzz1e9/HeH2.13` carries no certified generation
today (READ, its `state_index.json`: `latest_certified` is null), and this is
a certified state of exactly that configuration under the current binary.

They are NOT published, by the rule of this item. They sit in the publisher's
own layout, each with its `states/<generation>/` pair, `certification.txt`,
`manifest.json` and seed identity, and the decision to publish is the
advisor's.

## 7. What this result is, in the plan's words

Two branches were tested and they answer differently, so the sentence is
given twice.

**On the He/H = 9.7 branch, the lowest certified XUV on the tested branch,
with this grid, boundary law and numerical method, is 0.07 of the fiducial
GJ 1132 spectrum**, F_XUV = 2.341 erg cm^-2 s^-1 over 10-1300 A at the orbit
of LHS 1140 b.

**On the He/H = 2.13 branch, the lowest certified XUV on the tested branch,
with the same grid, boundary law and numerical method, is 0.01**, F_XUV
= 0.3345, which is the bottom of the spectrum grid and not a limit the ladder
found: the continuation was still taking 0.12 to 0.22 dex steps at 17 to 19
passes each when it ran out of spectra, so on that branch this work measured
no lowest XUV at all, only that 0.01 is reached.

The method both sentences are conditioned on, for
`atomic_scalar_gj1132x<XUV>_kzz1e9`: 500 cells on the Mixed grid of
`models/current_grid_Hydro_ioniz.txt`, the characteristic-face reservoir base
of boundary model `_v3`, the Roe flux with the well-balanced option, WENO3 in
the stationary route, the partitioned stationary solve alternating a JFNK
hydrodynamic solve with the element relaxation, `EXHALE_PTC_DTAU0 = 1.0e8`,
40 outer passes, and the certification tolerances of the inventory as the
binary `89ec67149aa452eea3deeb19052f83d1` applies them.

**The first sentence is never a statement that no stationary atmosphere
exists below 0.07 at He/H = 9.7.** A fold of the branch, an unresolved branch,
a discretization artifact of this grid and a failure of this numerical method
all produce the picture measured here, and this ladder does not separate them;
section 5 lists the three that are consistent with the numbers and says which
has the most direct support and that none is established. The He/H = 2.13
ladder is direct evidence for the caution: the same code, grid, base law,
flux, spectra and controls descend the whole decade at a neighbouring
composition, which is what a numerical or discretization limit would be
expected to do and what the absence of a solution would not. Time dependence cannot be read off it either:
that would take an accurately integrated physical time calculation, and none
was run.

Two further limits on the sentence, stated so that it is not read wider than
it is:

- It is a statement about a CONTINUATION FROM A NEARBY CERTIFIED BRANCH. The
  x0.01 cases of this family are not thereby shown to be unreachable: they may
  be reachable from their own stored states, by a longer alternation, or by an
  instrument this item did not use. A separate run of this campaign is
  advancing exactly that experiment from the stored x0.01 states, and the two
  are complementary (section 8).
- The earlier continuation of 2026-09-17 (`models/.L25/ladder/`, L25 step 4)
  reached 0.01 on this same case with a different binary, a different boundary
  model and a different seed, and its rungs are recorded as certified. This
  work does not contradict that record and does not reproduce it: what it
  measures is the current operator on the current branch, which is the only
  thing a Mode C result can be about (plan section 3).

## 8. What this is not a control for

Another run of this campaign is advancing the three x0.01 cases as a longer
alternation of twelve passes each, from THEIR OWN stored states. The x0.01
rung of this ladder, had it been reached, would have been the same case from a
DIFFERENT seed, the certified x0.015 state. **The two experiments are
complementary and neither is the other's control**: they share the case, the
binary and the grid, and they differ in the initial state, which is exactly
the variable a continuation exists to vary. A disagreement between them would
be information about the basin, not about either run's correctness.

## 9. The artifacts

All under the scratch directory
`/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/59a81fec-9ef0-4992-9a57-6494f229e6bb/scratchpad/stageF/`
(its `README.md` gives the layout, and the paths of every rung).

| what | where |
|---|---|
| the driver | `stageF/ladder.sh` |
| the He/H = 9.7 ladder | `stageF/L_scalar_HeH9.7/`, one directory per rung, with `ladder.log` |
| the refused rungs, kept apart by their seed | `x0p05_seeded_from_x0p07/`, `x0p06_seeded_from_x0p07/`, `x0p065_seeded_from_x0p07/` |
| the He/H = 2.13 ladder, all six rungs certified | `stageF/L_scalar_HeH2.13/x0p07`, `x0p05`, `x0p03`, `x0p02`, `x0p015`, `x0p01`, with `ladder.log` |
| the one-line state of every rung of both ladders | `stageF/audit/summarize.sh` |
| the pseudo-time control experiment | `stageF/probe_dtau0_default/x0p07_dtau0_1p0/` |
| the record of each rung, as this document quotes it | `stageF/audit/rungs/*.md` |
| the spectrum measurement | `stageF/audit/measure_sed.py`, `sed_bands.txt`, `sed_md5.txt` |

Every rung directory carries `SIDECAR.txt` (the gate of plan section 4.2: the
environment in force, the input state and its md5, the executable and its md5,
the mode, the stage, the exit class), and the state pair of every rung that
certified is a published generation in the contract of `MODELS.md` section 9,
with its `certification.txt`, its `manifest.json` and its seed identity,
inside the rung directory.

**Nothing was published into the catalog.** No file under `LHS1140b/models/`
was created or modified by this work. Three files were added to
`LHS1140b/sed/`: the two spectra the halved steps needed
(`..._xuv0p06.txt`, `..._xuv0p065.txt`) and `scale_sed.py`, the script that
writes them and that reproduces the existing grid byte for byte. Its
`README.md` was corrected where it described a state the directory had left
(section 10).

## 10. Two things repaired outside the item, and reported

Both are in `LHS1140b/sed/`, the directory this ladder's one-parameter family
is built from, and both were found by checking the scaling clause of section
14 rather than by looking for them.

1. **The scaling rule existed only as prose.** Sixteen scaled spectra were on
   the tree and no script wrote them, so the rule the model catalog's XUV axis
   rests on could not be re-run or checked. It is now
   `LHS1140b/sed/scale_sed.py`, and `--check` MEASURED that it reproduces all
   sixteen byte for byte. Nothing on the tree was rewritten: the script
   refuses an existing path.
2. **`LHS1140b/sed/README.md` described a normalization the directory no
   longer uses.** Its closing paragraph said "The shipped products are now
   integral-matched (`A_eff` in the headers)", while the shipped
   `lhs1140_sed_gj1132_at_b.txt` carries `A = 0.5900 (catalog X-ray flux
   ratio)` and `make_sed.py` sets exactly that, the integral-matched switch
   having been reverted on 2026-08-24 for the reason its own docstring gives.
   Its product table also quoted the superseded `adapt-var-res` generation
   (F_XUV 39.0) rather than the file on disk. Both were corrected with values
   MEASURED on the files as they stand (F_X 4.407, F_XUV 33.449, F_EUV 10.051
   for the GJ 1132 product; 3.492, 169.931, 70.181 for the GJ 699 one), and a
   section on the scaled grid was added. This matters to more than the
   README: the fiducial normalization is the zero point of every XUV factor in
   the model catalog, and the file said it was one thing while being another.
3. **Three smaller corrections in the same README, made while editing it.**
   Its "Inputs kept" paragraph named only the `adapt-var-res` FITS files while
   `make_sed.py` reads the `const-res` pair, so both are now named with which
   one the builder reads; the sentence "That catalog-ratio C is what
   `make_sed.py` uses (0.5926 ..., 0.0559 ...)" quoted the unrounded ratios
   while the script uses the paper's rounded 0.59 and 0.056, so it now says
   which value the files carry; and the forbidden wording was removed, a
   software-engineering term in the X-ray band paragraph replaced by what it
   meant there, "any model built on the proxy SED", together with the file's
   em-dashes.
