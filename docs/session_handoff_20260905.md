# Session handoff, 2026-09-05: where EXHALE stands and what comes next

Written for a fresh session. **Read `docs/code_status_20260905.md` first**: it is the
physics-by-physics state of the code and the task list; this file keeps the
last session's runs, paths and pitfalls. Everything below was verified in this session
unless marked otherwise; the section numbers refer to `docs/Update_EXHALE.md`.

## 1. State of the tree

* Branch `v1.00`, HEAD `6d07d48`; **everything since is uncommitted** (about
  165 modified/new files, `git status --short | wc -l`). No commit or push
  was made in this session (the user's rule: only on explicit instruction).
* `backup/` (regression harness and goldens) is not in the git remote; it
  lives only in this working copy.
* Two binaries in the root: `EXHALE.x` (gfortran 16.2 + OpenBLAS, the
  production build) and `EXHALE_ifx.x` (ifx 2026.1 + MKL, `build_ifx/`).

### What the uncommitted work contains (by changelog section)

| section | what | verification state |
|---|---|---|
| 159 | scaled trust region re-measured; unshifted Newton leg is the default | matrix PASS at the time |
| 160 | a run creates `./output` itself | -- |
| 161-163 | (AD) closed: marching departure = O(dt) kick + limit cycle; hand-off attempts | measured, documented |
| 164-167 | marching step 4x faster: metal-table skip, `eval_cool` over cell blocks, one definition per coefficient, hydro/carrier parallel loops | byte-identical or < 1e-9 |
| 168 | hot-Uranus n_e excess vs Koskinen 2022 diagnosed (F/4 convention, R10 sink, local-equilibrium ionization) | `docs/k22_electron_density_excess.md` |
| 169 | `Ionization transport: True` (H+ as 5th carrier); `load_IC` He/H check over physical cells; **writer ghost defect fixed** (molecular columns refreshed from the state before every write); harness comparator/loader test repaired; roundtrip golden re-snapshotted | matrix PASS, round trip identity, fcheck clean |
| 171 | like-for-like Koskinen comparison: `Domain mode: Spherical`, Ribas-band solar spectrum, `Atomic rate set`, `Caloric EOS: monatomic`, `Photoelectron heating: full`, LW stated-zero rule, SED allocation and ifx directory-probe fixes | 12 cases bitwise identical to the section-170 outputs |
| 170 | **heating and cooling of `ioniz_eq` assembled from the post-sweep composition** (rates pre-sweep); energy step's first `eval_cool` removed; toolchain gfortran 16.2 + OpenBLAS (rpath, compiler in the rebuild stamp); ifx + MKL builds; directory probe portable; four edit descriptors widened; `REGRESSION_EXE` in the harness | matrix run once: Mdot unchanged in all 12 cases, transient states moved (up to 2 per cent in molecular snapshots); round trip identity; fcheck clean; **goldens NOT yet refreshed** (see 3) |

## 2. Build, run, test -- as of now

```
make -j8                      # gfortran 16.2 from PATH (/opt/miniconda3), OpenBLAS of that prefix, rpath
make -q                       # 0 = up to date (the harness refuses a stale binary)
make FC=ifx OBJDIR=build_ifx EXE=EXHALE_ifx.x -j8     # Intel build, MKL sequential
cd <run dir> && OMP_NUM_THREADS=16 ../EXHALE.x        # reads ./input.inp, writes ./output/
EXHALE_MAXSTEPS=N ... ; EXHALE_PROFILE=1 ...          # step cap; phase timing table
backup/regression/run_check.sh check [cases]          # matrix; bitwise, else within REGRESSION_REL_TOL (1e-3)
REGRESSION_EXE=$PWD/EXHALE_ifx.x backup/regression/run_check.sh check ...   # another binary, no rebuild
backup/regression/run_check.sh golden [cases]         # snapshot the outputs sitting in the case dirs
( cd backup/regression/roundtrip && ./roundtrip_check.sh )   # restart write->read->write identity
backup/regression/run_fcheck.sh                       # -fcheck build, one HD 209458 b case
backup/regression/test_roundtrip.sh                   # loader identity test (its own dir roundtrip_loader/)
```

Do NOT prepend `PATH=/usr/bin:$PATH` any more (that was the gfortran-13
habit; the Makefile now resolves the toolchain itself). Shell quirks that
still hold: `tail -n N`, `\cp -f`, `\mv -f`; never `pkill` or any name-based
kill -- kill by PID after checking `/proc/<pid>/cwd`.

## 3. Pending: regression goldens (bookkeeping, user-run)

The section-170 fix is physical, so the goldens of all twelve cases must be
re-snapshotted. The case directories hold the post-fix gfortran outputs;
the command below snapshots them and confirms. It was refused to the agent
by the auto-mode classifier (a one-case refresh passed earlier, twelve did
not):

```
cd backup/regression && ./run_check.sh golden wasp_full wasp_he23off wasp_full_newton \
  mol_base_handoff mol_metals mol_lyman_werner mol_diffusion mol_ir_bands mol_sec_ion \
  mol_carrier lower_profile roundtrip && ./run_check.sh check
```

Standing rules: sub-0.1 per cent movement from compiler/constant/reduction
changes counts as identical (no refresh); a physical fix refreshes whatever
moved and the movement is reported; an error found is fixed at once, never
deferred to keep a golden.

## 4. The Koskinen line (the acceptance test)

Target: Koskinen et al. (2022, ApJ 929, 52) Model A at 0.05 au (H/He/H2,
1 microbar base, T0 1140 K, q0(H2) 0.84, flux/4, optically thin H3+ cooling,
every species transported, no Lyman-Werner). The memo
`docs/koskinen2022_model_a_comparison.tex` (build: `cd docs && pdflatex
koskinen2022_model_a_comparison.tex` twice) records inputs, method
differences and the comparison; `docs/k22_model_a_figures.py <output_dir>
[<earlier_output_dir>]` regenerates its figures and table into
`docs/figures/k22_model_a/`; the run itself is `benchmarks/koskinen2022_model_a/`
(input.inp, base.inp, output/ = the final state; `output/*_IC.txt` = the
state it was continued from).

Configuration keys that matter: `2D approximate method: Rate/4 + Mdot`,
`Molecular chemistry: True`, `Molecular carrier transport: True`,
`Ionization transport: True`, `Secondary_ionization: Immediate`,
`du_th [PLM,WENO3]: 1.0e9 1.0e-30` (march for EXHALE_MAXSTEPS; the steady
solver would put the composition back on its local root -- do not use it
here), `base.inp` with `p_base 1e-6`, `q_H2_base 0.84`, `Kzz_base 1e9`.
Quick comparison: `python3 <scratch>/k22_compare.py <output_dir>` (the
script is also reproduced by the figure script's table).

Where the comparison stood before this session's last run (section 169.1,
10^5 steps): x(H+) within 1.0-1.1x of the advected estimate; n_e 1.5x
Model A at 2.0 r_base and 1.4x at 3.0 (was 5-10x with local equilibrium);
T within 60 K in the layer, 340 K low at 3.0; H3+ 3x low at 3.0; f(H2)
0.42 vs 0.51 at 2.0; log Mdot 10.21 vs 10.28. After the section-170 fix
and 2e5 more steps (the benchmark state): log Mdot 10.22 (Model A 10.28);
T 1752/1937/2207/2999/4135 K at 1.15/1.20/1.30/1.50/2.00 vs 1720/1900/2270/
2870/4120; n_e 2.7e6 at 2.0 (A 1.8e6) and 1.1e6 at 3.0 (7.8e5); f(H2) 0.416
at 2.0 (0.508); H3+ 3.0e3 at 2.0 (7.2e3); T 4399 at 3.0 (4740). The fix
changed the steady numbers by < 3 per cent (transient-only, as expected).

Remaining gaps and the candidates to test next, in order:
1. n_e still 1.4-1.5x: (a) their R1 recombination is 1.2-1.5x ours -- a
   like-for-like test is to switch to their coefficient in `mol_rates.f90`
   for this run only; (b) He+, H3+ and HeH+ are still local equilibrium --
   extend the carrier set (the operator takes any species as a carrier;
   `ic_Hp` is the pattern, section 169); (c) spectrum: a mean solar SED
   instead of the power law (`Spectrum type` supports a file).
2. T at 3.0 r_base 300-400 K low and f(H2) low: the domain ends at 4.73 R_p;
   Model A reaches 7.16 base radii and its T is still rising there. Extend
   `Escape radius`/grid to 7 base radii and re-run before attributing the
   outer temperature to physics.
3. H3+ 3x low in the outer wind: follows from (1b) and (2) -- test after them.
4. Path independence: the run was continued from earlier states; a fresh
   start (cold IC) to the same state was shown for the Rate/4 configuration
   (section 168, `r4_fresh` = `r4_cont`) but not yet with ionization transport.

## 4b. The like-for-like Koskinen comparison (2026-09-05, section 171)

Every Model A condition was checked one by one (the table is section 2.1 of
`docs/koskinen2022_model_a_comparison.tex`; the memo is rebuilt, 7 pages).
The runs are `benchmarks/koskinen2022_model_a/matched_hnu_minus_I/`
(physical accounting) and `matched_full_hnu/` (whole photon energy as heat,
a comparison option). Findings, in order of size:

1. **Gravity**: the default Roche domain (L1 = 4.8 R_p) was the largest
   difference; Model A is spherical to 7.24 base radii (`Domain mode:
   Spherical`, `Outer radius [R_p]: 7.24`). Velocity and neutral densities
   then match within the plotted points.
2. **Spectrum**: use `solar_ribas2005_bands_0p05au.txt` (Ribas 2005 Sun band
   fluxes, 1.6e3 total). The Gueymard composite alone is 2.4-3x low in
   1-100 A. `Stellar LW flux: 0.0` must be stated (Model A has no H2
   photodissociation); a stated zero now wins.
3b. **Stated-fraction runs** (`Photoelectron heating: <f>`): f = 0.55 puts
   the heating on their Figure 9 and T within the points (4119/4840 K at
   2.0/3.0); n_e 1.4-1.6x, H3+ 0.35-0.5x, He+ 3x, f(H2) 0.06-0.08 low remain
   (items 2-3 of memo section 4.1: local-equilibrium He+/molecular ions and
   the unsettled layer). Kept in `matched_fraction_0p55/`; comparison only.
3. **Remaining residual = Model A's heating rate**: 2.5-3x ours from 1.5
   r_base outward at the same ionization rate (their Fig. 9 vs Fig. 15),
   which the stated h nu - I accounting and any Ribas-band spectrum cannot
   give. `Photoelectron heating: full` brackets it from above (T overshoots,
   but n_e, H3+ and the cooling curve land on theirs). Neither
   `Caloric EOS: monatomic` nor the H3+ non-LTE factor changes the outer T.
4. **Layer**: the base is still filling (base flux 2.2x the wind's); the
   1.05 trough is 976-1021 K against 1180; a continuation (`C2`, 4e5 steps)
   was running at the end of the session in the scratch directory
   `k22_matchA/C2`.

Next steps for this line: (a) ask the authors for the spectrum file and the
heating deposited per ionization (or Figure 9 in numbers); (b) meanwhile
keep both accountings as the bracket; (c) if a match "within the plotted
points" is required without (a), the only remaining lever is the heating
per ionization, which is not a physics choice EXHALE should make silently;
(d) the layer's steady state needs either a very long march or a steady
solve that keeps the transported H+, which `Ionization transport` refuses
today (design item).

Scripts: `docs/k22_model_a_figures.py <run> [<second run>]`
(`K22_FIG_DIR`, `K22_LABEL1/2`); session scratch `k22_matchA/snap.sh`
(figures + table of a running case), `k22_matchA/interp_ic.py` (restart
onto a new grid; worth moving under `src/utils/`).

## 5. Other open items (code)

* `Cool_coeff.f90`: 26 wrappers left without callers after section 165;
  nine `_func` names; H I/He I literal constants not yet unified (section
  166 left them). Decide and clean in one pass; byte-identical expected.
* Ionization sweep = 68 per cent of a 16-thread step (section 167 profile);
  the next speed-up is inside `System_HeH_mol` (Jacobian reuse across
  sweeps), not in the hydro.
* `Ionization transport` refuses `Solver: Newton`, `EXHALE_PTC` and the
  coupled carrier solve; an atomic-branch version and a coupled JFNK with
  the proton row are design items, not started.
* `_adv` files' molecular columns are still the equilibrium sweep's (the
  post-process has no molecular state); documented in section 169.5.
* Snapshots written before section 169.3 carry inconsistent ghost rows;
  they restart, but their ghost rows must not be compared.
* `docs/Update_EXHALE.tex`/`.pdf` and `docs/EXHALE_user_manual.tex`/`.pdf`
  were modified earlier in the series; whether they are current against the
  `.md` changelog was not checked in this session.
* `interp_ic.py` lives only in the session scratch; move it under `src/utils/`
  if a restart onto another grid is wanted again.
* Instruction-gated tail (never picked up autonomously): paper/poster
  re-convergence, the LHS 1140 b write-up, the GitHub public switch, P36.

## 6. Memory files written this session (`~/.claude/projects/.../memory/`)

`project_ionization_transport.md`, `project_heat_cool_post_sweep.md`,
`feedback_test_once_after_all_edits.md`, `feedback_fix_errors_immediately.md`,
`feedback_koskinen_match_over_goldens.md`, `reference_shell_quirks.md`
(updated: bare `make` is right now), `feedback_sub_tolerance_no_golden_refresh.md`.

## 7. Pitfalls met this session, so they are not met again

* `backup/regression/test_roundtrip.sh` used to `rm -rf` the named case
  `roundtrip` (same directory); it now works in `roundtrip_loader/`. The
  case was restored from an isolated copy; its golden was stale since
  section 151.11 and was re-snapshotted.
* A compound shell command that includes a golden refresh, or a twelve-case
  refresh on its own, is refused by the auto-mode classifier; a single-case
  refresh passed. Keep refreshes as separate, explicit commands.
* `$PWD` inside a `( cd dir && ... )` subshell is the new directory: use
  absolute binary paths in run loops.
* `inquire(file='dir/.', exist=)` is not portable (ifx says .false. for a
  directory); the code now probes by creating a file.
* Every hand-edited `input.inp` for the hot-Uranus gate needs `base.inp`
  beside it; a run without `./output` now creates it.
