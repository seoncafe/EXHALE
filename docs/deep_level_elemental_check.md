# The deepest-level elemental check: what it can see, and why it refused reservoir He/H = 8.5-8.7 on LHS 1140 b

Written 2026-08-29. Primary measurement record:
`LHS1140b/exhale/nh_refusal_diagnosis/` (README, four scripts, five tables,
38 adapter runs under `scan/`). Every number quoted here was read back from
that record or recomputed from the profiles it wrote; no EXHALE run was made
for this document.

> **Superseded in part, the same day.** Sections 1--7 stand: they diagnose the
> photochem 0.8.4 solver, and that diagnosis is what led to the correction.
> Sections 8, 9 and 10 have been rewritten against the corrected Photochem
> build (`photochem 0.9.0` with the Equilibrate, gasgiants and clima patches of
> `docs/photochem_solver_modification_implementation.md`). The tolerance
> decision of section 8 was made for the old solver and no longer holds; the
> `--abundance-tol` default is now `1e-10`, not `1e-3`.

---

## The judgement first

The refusal is **the residual of Photochem's chemical-equilibrium solver, not
a failure of the photochemistry and not a leak of nitrogen**. The level the
check is stated at is Photochem's Dirichlet lower boundary, so the check
cannot see a photochemical or transport error at all; what it reads is the
equilibrium initialization, and the equilibrium solver is used deliberately
without enforcing convergence. The residual is a factor 4072 more visible in
N than in He, C or O, for a stoichiometric reason that has nothing to do with
nitrogen chemistry being harder.

The tolerance was therefore first raised from `1.0e-4` to `1.0e-3`
(`--abundance-tol`, `src/utils/lower_profile_schema.py`). That was a judgement
about the old solver. It was then superseded: with the solver corrected the
residual it was accommodating is gone, and the default was tightened to
`1.0e-10`, the value the check was designed at. Section 8.

---

## 1. The symptom

In the LHS 1140 b closure ladder
(`LHS1140b/exhale/crossings_gm25/`, `results.txt`), reservoir values
He/H = 8.5, 8.55, 8.6, 8.65 and 8.7 could not be run. The adapter
(`src/utils/photochem_to_lower_profile.py`, through
`lower_profile_schema.write_handoff`) refused the column before any wind step:

```
REFUSED: the N/H ratio at the deepest level, 8.185823e-05, departs from the
input 8.184650e-05 by 1.43e-04 (tolerance 1.0e-04): the element columns do
not sum to the input abundances
```

The measured departures across the band, from the closure logs:

| reservoir He/H | N/H departure |
|---|---|
| 8.50 | 1.02e-4 |
| 8.55 | 1.22e-4 |
| 8.60 | 1.43e-4 |
| 8.65 | 1.67e-4 |
| 8.70 | 1.94e-4 |

**Independent of the trial flux.** Two closure lineages with different trial
escape fluxes reach exactly the same numbers: `fc_*` at
`--trial-flux-H 5477925.0 --trial-flux-He 19195010.0` and `fcA_*` at
`--trial-flux-H 4768206.9882 --trial-flux-He 20433648.55`, agreeing to the
three digits the message prints at every one of the five reservoirs. The
diagnosis scan of section 5, run standalone with no wind and no closure
iteration at the `fcA_` flux, reproduces the same values a third time.
(`results.txt` states this as three different trial fluxes; what is verifiable
from the surviving logs is two distinct flux pairs plus the standalone scan.)

---

## 2. The level the check reads is a Dirichlet boundary

`gasgiants._initialize_atmosphere`
(`photochem/extensions/gasgiants.py`, line 371 of the installed copy) sets, for
every gas species of the mixture handed over,

```python
self.set_lower_bc(sp, bc_type='press', press=Pi)
```

so the bottom cell is held at a fixed partial pressure per species. It carries
the equilibrium initialization unchanged, and the photochemical solve cannot
move it.

Measured, at He/H = 8.6
(`probe_initial_vs_steady.py` -> `probe_heh8p6_summary.txt`): the deepest-level
El/H ratios of (a) chemical equilibrium plus the quench overwrite as handed
over on the climate grid, (b) the photochemical model's initial state, and
(c) the converged steady state read back from the written profile, agree to
**every one of the 13 digits printed**:

```
     He  8.599999697179e+00   dev -3.5212e-08
     C   2.775879902313e-04   dev -3.5191e-08
     N   8.185823483630e-05   dev +1.4338e-04
     O   6.061779786554e-04   dev -3.5212e-08
```

identical in all three. The deep `q_NH3` is `8.995632484e-06` converged
against `8.995632e-06` initial. The stoichiometry of the (a) -> (b) change is
at the level of double-precision rounding (relative `1.1e-16` on H and He).

Two consequences, both of them measured rather than argued:

- **A photochemical or transport leak never reaches this check.** Whatever the
  network does above the bottom cell, the bottom cell is pinned.
- **The refusal cannot depend on the trial escape flux**, which enters only as
  an upper boundary condition. Section 1 is the direct confirmation.

---

## 3. The equilibrium solver is used unconverged, and still reports convergence

`composition_at_metallicity` (`photochem/extensions/gasgiants.py`, around line
920) loops over levels, retries `equilibrate.ChemEquiAnalysis.solve` at five
temperature perturbations, takes the first that reports convergence, and then
keeps the result either way. Its own comment, line 921:

```python
        # Do not enforce convergence.
```

Called directly at the deepest level of the LHS 1140 b climate column
(`equilibrium_closure.py` -> `equilibrium_closure_heh.txt`, at
P = 1.673098e7 dyn/cm2, T = 428.10 K), the solver returns `converged = True`
at **every** reservoir value swept, while its elemental closure on N misses
the abundances it was handed by up to **1.72e-4** (at He/H = 9.85; 6.9e-5 at
8.50, 8.1e-5 at 8.70). The same table's last column shows the departure is
present in the solver's own `molfracs_atoms_gas`, not introduced by the
adapter's carrier summation: `dev_N` and `dev_N_atomsgas` are equal to the
printed precision in every row.

So the check, at the level where it is stated, is (i) a guard on the adapter's
own carrier summation, which is what it was written for, and (ii),
unintentionally, a convergence test of a solver that is used without enforcing
convergence.

---

## 4. Why nitrogen alone: NH3 stoichiometry, a factor 4072

In all 38 runs of the scan the ratio of the N departure to the He departure is
the same number:

```
dev_N / dev_He = 4072   (32 of 38 runs; 4069-4086 over all 38,
                         the outliers being runs whose residual is
                         so small that the printed digits round)
```

against

```
1 / (3 N/H) = 1 / (3 x 8.18465e-5) = 4073
```

The reading: at the deep level essentially all nitrogen is NH3 (measured on
the He/H = 8.6 profile, the N2 nuclei contribution is 1.0e-5 of `X_N` and HCN
is 2.2e-25 of it), so **one excess N rides on three excess H**. He, C and O
carry only the common error in the hydrogen denominator; N carries that plus
the relative error of the equilibrium NH3 mole fraction, and the second term
is larger by 1/(3 N/H).

**The check is therefore not element-neutral.** At this level it tests the
equilibrium NH3 mole fraction 4072 times more tightly than it tests anything
else, which is not a property anybody chose.

---

## 5. The composition is smooth; only the residual jumps

`equilibrium_deep_composition.py` -> `equilibrium_deep_composition.txt` reports
the deep mixing ratios themselves, to separate "the chemistry changed branch"
from "the solver residual grew". Across 8.45 -> 8.75, H2, H2O, CH4, NH3, N2
and He all vary **monotonically and by less than 1 %**, with no jump anywhere
(NH3 9.14461e-06 -> 8.84864e-06, N2 3.25565e-10 -> 3.36526e-10, He
0.944168 -> 0.945979).

What jumps is the residual. From `scan_table.txt`:

| reservoir He/H | N/H departure | jump |
|---|---|---|
| 8.71 | 1.995e-4 | |
| 8.72 | 5.227e-6 | factor 38 across a 0.1 % change in composition |
| 9.15 | 6.570e-5 | |
| 9.20 | 4.500e-8 | factor 1500 |
| 10.25 | 2.018e-4 | |
| 10.30 | 3.478e-8 | factor 5800 |

A smooth composition with a residual that changes by three orders of magnitude
between neighboring reservoir values is what a solver stopping criterion looks
like, not what chemistry looks like. The branch-change hypothesis is refuted by
the composition table.

---

## 6. Not a special band, and not a special planet

At the old tolerance of 1e-4 the same scan also refuses:

- **He/H = 6.0**, departure 1.31e-4;
- **He/H = 10.20 and 10.25**, departures 1.48e-4 and 2.02e-4
  (10.1 sits at 6.7e-5, just under the old threshold).

And at **solar composition**, He/H = 0.0969, the standalone equilibrium solver
misses N/H by **1.50e-4 at 800 K** and **1.14e-4 at 1600 K**
(`equilibrium_closure_temperature.txt`) -- both would have been refused.

Elsewhere in the LHS 1140 b tree, the same refusal is already on record at
He/H = 14.643 (1.23e-4) and 15.090 (**2.12e-4**) with `K_zz = 1e8`, and at
He/H = 25.0 (1.28e-4) in the climate probe of the same scan
(`LHS1140b/exhale/kzz_profile_scan/`). 2.12e-4 is the largest departure
recorded anywhere in this tree.

The two HD 209458 b departures that fixed the original 1e-4 default --
**1.6e-7** with the Zahnle H/He/N/O/C set and **2.7e-5** through the VULCAN
NCHO network (`docs/Update_EXHALE.md` section 77) -- were a fortunate sample,
not a floor.

---

## 7. There is no accounting error

At the deepest level the profile's diagnostic columns account for **100.00 %**
of the `X_N` column: recomputed from
`scan/heh8p6/lower_atmosphere_profile.dat`,

```
(q_NH3 + 2 q_N2 + q_HCN) / (X_N * n_H)  =  1.0000   at the deepest level
                                        =  0.9391   at the 1 microbar match
```

where `n_H` is the H-nuclei mixing ratio summed over its own carriers.

**This distinction matters and has been confused before.** A reading that
"the nitrogen carriers sum to only about 96 % of `X_N`" is a statement about
the **matching level** and about the `q_*` **diagnostic** columns, whose set is
deliberately short (the remainder there is N, NH, NH2 and the like). The check
does not read the matching level and does not read the `q_*` columns: it reads
`X_N` at the deepest level, which is summed over every carrier the network
carries. The shortfall in the diagnostic columns is a property of that column
set, not a miscount, and it is not what the refusal is about.

---

## 8. What was changed, and what replaced it

**The default is `1.0e-10`.** The route there had two steps, and only the
second one is current.

### 8.1 The first step, against the old solver: `1.0e-4` -> `1.0e-3`

`src/utils/lower_profile_schema.py`: the `--abundance-tol` default was raised
to `1.0e-3`, and both the option help and the comment at the check were made to
state what the check can and cannot see. The grounds were:

1. **Headroom.** The worst departure recorded anywhere in this tree was
   2.12e-4, so 1e-3 left a factor 4.7.
2. **Almost no detection power given up.** At the deepest level each element
   sits in one dominant carrier, so a genuine miscount -- a carrier omitted
   from the sum -- is an O(1) error and is caught by any tolerance at all. The
   smallest conceivable miscount is dropping N2, and on the He/H = 8.6 profile
   that is **1.0e-5** of `X_N`: it was already invisible at the old 1e-4. Where
   N2 is a large enough share to matter it stays visible at 1e-3 as well --
   measured over the 38 profiles of the scan, the N2 nuclei share of `X_N`
   runs from 4.2e-6 (He/H = 12, deep T = 379 K) to 4.0e-3 (He/H = 1, deep
   T = 588 K), rising with the deep temperature, so exactly where dropping N2
   would be a real loss it is above 1e-3.
3. **The check was measuring the wrong thing anyway.** Sections 2 and 3: it
   cannot test conservation through the photochemistry, and what it was
   actually rejecting on was the stopping point of a solver that is used
   unconverged on purpose.

Grounds 2 and 3 still hold. Ground 1 does not: it was a statement about how far
the 0.8.4 equilibrium solver missed, and that is what got fixed.

### 8.2 The current default: `1.0e-10`

Photochem was corrected instead
(`docs/photochem_solver_modification_implementation.md`): Equilibrate now tests
each element against its own abundance rather than scaling every elemental
residual by the largest one, and `composition_at_metallicity` refuses a
condensate-free level that does not close. The residual section 3 measured is
therefore no longer there to accommodate, and the default was set to `1.0e-10`,
the value the check was designed at.

The grounds are measurement on the corrected build:

- Over the **68 handoff writes** made by the corrected build under
  `LHS1140b/exhale/` (43 in `nh_refusal_diagnosis/scan_pc090/`, 25 in
  `crossings_pc090/`), the largest deepest-level departure is **3.28e-13** and
  the median **1.73e-14**. `1e-10` stands a factor 305 above the worst of
  those. Re-solved standalone, the 38 saved deep states of section 5 give a
  largest residual of **3.539391003e-13** and a median of **1.915134717e-14**,
  against **2.018197773e-4** for the same measure on 0.8.4.
- Detection power is unchanged in kind: ground 2 above is about carrier
  bookkeeping, not about the solver, and `1e-10` is orders of magnitude below
  the 1.0e-5 to 4.0e-3 range a dropped N2 would produce.
- **A 0.8.4 handoff is refused at this default, and that is the intended
  signal.** Over the 43 stored 0.8.4 handoffs of `crossings_gm25/` the
  departures run 3.48e-8 (best) to 8.50e-5 (worst), median 6.21e-6: even the
  best 0.8.4 write misses `1e-10` by a factor 348. Reproducing a stored 0.8.4
  result therefore needs an explicit `--abundance-tol` on the command line,
  which is how a deliberate reproduction states that it accepts that residual.
  The scans under `nh_refusal_diagnosis/` already pass `--abundance-tol 1.0`
  and are unaffected.

Still not changed: the check itself, its location, or the refusal policy.

## 9. What this does to the LHS 1140 b result: nothing measurable

The flux-closure crossing was re-determined end to end on the corrected build,
11 arms, same EXHALE binary, same scripts, same parent seed, same tolerances --
only the interpreter differs. Record:
`LHS1140b/exhale/crossings_pc090/results.txt`.

| quantity | 0.8.4, 16 arms | corrected build, 11 arms |
|---|---:|---:|
| crossing, chord | 9.0484 | 9.0484 |
| crossing, quadratic | 9.0461 | 9.0461 |
| 1 sigma band | 8.4923 -- 10.7716 | 8.4649 -- 10.7357 |

- **The crossing does not move**, to every digit either solve prints. On the
  seven arms seeded identically in both builds the He I 10830 equivalent width
  moves by at most `1.2e-6` in relative terms and `log10 Mdot` is unchanged in
  its fourth decimal. That is the scientific result of the correction: the
  equilibrium residual was real and is now gone, and it was never reaching the
  observable.
- **The handoff profile does move.** At He/H = 9 the deepest-level N/H
  departure goes from 2.1e-5 to 1.2e-14, and trace species such as `q_H2O`,
  `q_CO` and `q_HCN` move by up to 1e-2 relative. None of that reaches He I
  10830, which is set by H, He and the temperature.
- **The 1 sigma band moves, but not because of the solver.** Its low end used
  to be interpolated across the un-runnable band He/H = 8.5--8.7; that band now
  runs (departures 1.4e-14 to 3.2e-14) and real arms at 8.5 and 8.6 carry the
  end, moving it from 8.4923 to 8.4649. Its high end moves from 10.7716 to
  10.7357 because the 10.6 arm was not seeded the same way in the two ladders
  -- a 1.1e-3 seed difference in equivalent width, the size of the seed spread
  the stored ladder already measured.
- **The refusals at other reservoir values are gone too.** He/H = 6.0, 10.20
  and 10.25, and the band 8.5--8.71, all pass at 1e-14 on the corrected build.
  Solar composition, He/H = 0.0969, is runnable end to end for the first time:
  it was previously blocked not by the tolerance but by the climate step
  (section 10).

---

## 10. The climate root solve: three compositions recovered, two others refused

Section 10 previously recorded three unexplained failures of `clima`'s
`surface_temperature_bg_gas`, reached through
`src/utils/radiative_convective_column.py:140`. The corrected build changes the
picture in both directions.

**Recovered.** He/H = 0.0969 (solar), 8.1 and 15 now return
radiative-convective solutions, with deep boundary temperatures 720.5, 432.4
and 393.1 K. He/H = 8.1 lies between its 8.0 and 8.2 neighbors as it should.
The bounded, scaled solves of
`docs/photochem_solver_modification_implementation.md` section 2.3 are what
removed them.

**Refused instead: He/H = 9.4 and 9.5, and this is not a regression.** Both
solved on 0.8.4 and now raise

```text
_clima.ClimaException: hybrd1 root solve failed: Could not bracket background
pressure in make_profile_bg_gas.
```

This section first read that as an isolated failure of the new bracketing scan,
with a net of three failures removed and one introduced. **The measurement
refutes both statements.** It is recorded in
`LHS1140b/exhale/clima_bracket_diagnosis/` and written up in
`docs/Update_EXHALE.md` section 91, which is the single place the diagnosis
lives. In short: the bracketing scan is not the cause -- the root is inside the
interval it scans and the residual crosses zero once, monotonically, in the
last scan step -- and the failure comes from the outer solve's
forward-difference Jacobian, whose `hybrd1` step of `sqrt(epsmch)*|x|`
(`2.4e-5` K at 400 K) is below the flux residual's own round-off jitter. The
bad step it produces sends the trial surface temperature to `6.5e8` K, where
every scan point is thermodynamically invalid.

The same defect is in 0.8.4, at the same rate and at other compositions:

| grid | corrected build | 0.8.4 |
|---|---|---|
| He/H 9.30--9.60, step 0.01 (31 points) | 9.40, 9.50 | 9.37, 9.45 |
| He/H 8.00--11.00, step 0.05 (61 points) | 9.40, 9.50, 10.45 | 8.10, 9.45, 9.85, 10.15 |

The earlier table in this section listed 9.45 as `not attempted` on 0.8.4; it
is a **failure** there. With the two builds measured on the same grid the "three
removed, one introduced" accounting does not stand, and it is withdrawn.
Neither value is a bracket arm of the crossing, and the solution is smooth
across the failures -- the deep temperature falls monotonically from 422.72470
K at He/H = 9.30 to 420.61128 K at 9.60.

**What is open.** `epsfcn = 1e-4`, which requires `hybrd` in place of `hybrd1`,
takes the 61-point grid from 54/61 to 61/61 in replication and does the same on
0.8.4's parameterization; applying it to the clima patch and re-measuring on a
rebuilt build has not been done. Until then, a composition that refuses is
recovered by passing `--climate-t-deep-guess` anything other than its default
400.0. What produces the `~4e-7` relative jitter in `ISR - OLR`, and whether
the same trap is reachable at other planets or other `K_zz`, is still unknown.
Item (M) of `TO_BE_DONE.md`.

---

## 11. Files

Measurement record, `LHS1140b/exhale/nh_refusal_diagnosis/`:

| file | what it holds |
|---|---|
| `README.md` | the record's own summary |
| `scan_reservoir.sh`, `scan/`, `scan_table.txt` | the adapter at 41 reservoir values, check disarmed on the command line, 38 columns produced |
| `probe_initial_vs_steady.py`, `probe/`, `probe_heh8p6_summary.txt` | equilibrium / photochemical-initial / converged at He/H = 8.6 |
| `equilibrium_closure.py`, `equilibrium_closure_heh.txt`, `equilibrium_closure_temperature.txt` | the equilibrium solver alone, against the abundances it was handed |
| `equilibrium_deep_composition.py`, `equilibrium_deep_composition.txt` | the deep mixing ratios through the refused band |
| `read_deepest.py` | the same departure, read from any written profile |
| `b1.log` .. `b5.log` | the scan batches |

The corrected build adds `scan_pc090/` beside `scan/` in the same record: the
same reservoir scan re-run on photochem 0.9.0, three points wider (9.45, 9.55,
9.6), from which the recovery of section 10 is read. The climate solve's own
diagnosis is a separate record, `LHS1140b/exhale/clima_bracket_diagnosis/`,
written up as `docs/Update_EXHALE.md` section 91.

Code: `src/utils/lower_profile_schema.py` (the check, the default, the two
explanatory texts). Ladder and its crossings:
`LHS1140b/exhale/crossings_gm25/results.txt` (photochem 0.8.4) and
`LHS1140b/exhale/crossings_pc090/results.txt` (the corrected build). Earlier
instances of the same refusal: `LHS1140b/exhale/kzz_profile_scan/`. The origin
of the 1e-4 default: `docs/Update_EXHALE.md` section 77. The correction itself:
`docs/photochem_solver_modification_investigation.md` (diagnosis) and
`docs/photochem_solver_modification_implementation.md` (what was built and
measured). The check's place in the handoff design:
`docs/phase_e_flux_closure_design.md` section 4.2.
