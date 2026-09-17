# Session handoff, 2026-09-15 (evening): the LHS 1140 b stationary series, where it stands and what comes next

Written 2026-09-15 at the user's request, before the campaign is re-run.
Every number is READ from the memos and logs named beside it; nothing here
was re-measured for this document. The previous handoffs are
`session_handoff_20260912.md` and earlier.

## 1. Where things stand, in five sentences

The LHS 1140 b campaign (`LHS1140b/models/`, 95 cases) was solved on
2026-09-14 on the well-balanced partitioned route and 83 cases certified;
the memo `docs/lhs1140b_exhale_vs_pwinds.pdf` and `Update_EXHALE_stage2.md`
section 11 were written from those results. The three cases at 0.01 of the
fiducial XUV do not solve (item L17, closed as diagnosed) and the nine
molecular cases were held by the H2 carrier row (L7e, closed) and then by the
base face (L7f, L7g). While closing those, item L21 found that the base
boundary condition took its entropy branch from the cell-1 velocity, which
on this planet is the odd-even artifact of P44, so **every LHS solution of
the campaign was solved with the interior isentrope at the base instead of
the stated 226 K reservoir**; the corrected boundary moves Mdot of the
fiducial case by -1.8 percent (log10 7.8756 to 7.8679, -0.0077 dex) and the
base temperature from 418 K to 238 K.
**The campaign therefore has to be re-run**; the user approved the re-run at
20:35 and it is running (section 11).

## 2. The campaign as it is now (pre-L21 base boundary)

| group (see `LHS1140b/MODELS.md` section 3) | cases | state |
|---|---|---|
| atomic prescribed, GJ 1132 and GJ 699, K_zz ladder 0 to 1e11, well-mixed, metals | 56 | certified (pre-L21) |
| XUV grid, scalar base, x0.10 to x0.33 | 14 | certified (pre-L21) |
| XUV grid, scalar base, x0.01 (HeH 2.13, 9.7) | 2 | NOT SOLVED (L17) |
| flux-closure rungs, fiducial column, He/H 2.09 to 12 | 9 | converged, k = 4 to 5 (pre-L21) |
| XUV grid above the fiducial column, x0.10 to x0.33 | 6 | certified (pre-L21) |
| XUV grid above the fiducial column, x0.01 | 1 | NOT SOLVED (L17) |
| molecular (three groups) | 9 | not run; gated on L7g and on the re-run |

Results of record (pre-L21): `MODELS.md` sections 7 and 8 (written by
`models/status.py --write`), every case's `REPRODUCE.md`, the memo tex and
its 20 figures (`LHS1140b/make_memo_figures.py`, `docs/figures/lhs1140b_*`).
Crossings He/H against the measured 1.108 %A: well-mixed 0.401, GJ 699 0.053,
K_zz 0 to 1e11: 2.72, 2.71, 2.63, 2.48, 1.93, 1.49, 1.14, 0.88; the closure
ladder has no crossing (EW 1.50 to 2.27 %A, all above the measurement); the
closed XUV grid crosses the 0.6 percent depth limit at about 0.08 of the
fiducial, the scalar one at 0.28. All of these are to be re-measured.

## 3. The items of `docs/PLAN_20260913_lhs_stationary.md`

Memos: `docs/lhs1140b_stationary_L<item>_2026091{3,4,5}.md` (32 files).

| item | one line | status |
|---|---|---|
| L1, L2, L3 | archived states never stationary below 1.1 Rp; no march settles; the base state exists | closed 09-13 |
| L4, L4b, L4d | recipe (`Well balanced: True`, `Secondary_ionization: Immediate`, `Restart intent: stationary`, dtau0); C I cooling band; stagnation detector | closed 09-13 |
| L4e, L4f, L4g | row scaling refuted; outer-cell energy row; pseudo-time ramp (dtau0 = 1) | closed 09-14 |
| L4h | ramp collapse on damped steps: x0.1 clip, counter only at dtau <= 10 dt_CFL, return to the best iterate; two designs refuted by measurement; citations verified against Mulder & van Leer 1985, Kelley & Keyes 1998, Coffey et al. 2003, Gropp et al. 1998 (the lam-keyed cut has no published counterpart) | closed 09-15 |
| L5a, L5b, L5c | seeds; element relaxation; the O(1) base energy residual is the HLLC contact-speed truncation, removed by the well-balanced option | closed 09-13 |
| L6, L6b, L6c | C I balance right; six C/N/O coolants moved to density-resolved CHIANTI tables; seven goldens attributed | closed 09-13/14 |
| L7, L7b, L7c, L7d | molecular initializer from a certified atomic state and its corrections | closed 09-13/14 |
| L7e | H2 carrier row 1.000 at 1.6 Rp: four solver-control defects fixed; species rows certify; H2 in the hot wind is made in place by He(2^3S)+H -> HeH+ -> H2+ -> H2 | closed 09-15 |
| L7f | base H2 chemistry vs the handoff: the 7-21x was the seed's root formula (production depends on n(H2)); the two agree to 1.7 percent; R20 retired (Schauer 1989 bounds the HeH+ channel at <= 1e-14), R17 corrected from Johnsen 1980 and Boehringer & Arnold 1986, H2+ + He -> HeH+ + H added from Black 1978 | closing (regression table pending) |
| L7g | molecular base 1023 K vs 418 K: the reaction-heat ledger does NOT double count (cycles close on the ionization potential once); what remains is the vibrational heat fraction (~1) against the H2 band cooling (0.2 percent of the deposit); on the L21 boundary the fiducial molecular case certifies at pass 11 (ghost 226 K, cell 1 808 K, mass row 4.9e-9) | closed 09-15 as diagnosis; the vibrational heat fraction named, not opened |
| L8 | outer boundary: the wind is transonic, critical point at 40 Rp outside the 30 Rp domain; Mdot insensitive, EW +-1.7 percent | closed 09-13 |
| L9 | read-only diagnostic and output path | open |
| L10, L11, L12 | `_adv` profiles; composition ladder seeding; ionization closure above 4 Rp (stage A done, B-F open) | closed / decision recorded |
| L13 | closure seed contract: every element reservoir and the grid R0 carried onto the iteration (`map_state_to_grid.py --reservoir El/H v`, `src/utils/profile_match_level.py`) | closed 09-14 |
| L14 | low-XUV progress control: a pass counts as progress when the composition's own distance falls and the hydro rows do not bind (`EXHALE_OUTER_PASSES` 40, `EXHALE_PTC_DTAU0` a default) | closed 09-15 |
| L15 | the JFNK solve is not reproducible at eight OpenMP threads (one thread is bitwise) | open |
| L16 | L4g section 7.3 numbers not reproduced by the current source | open |
| L17 | the 0.02-XUV wall: the base is subsonic to Mach 1e-8 with alternating-sign velocities and the HLLC contact selection puts every base cell on a kink of the residual; not solvable without a low-Mach flux | closed 09-15 (diagnosed, not solved) |
| L18 | a certified state did not reproduce its certification from its files: the loader rebuilt the density from the species instead of reading the `rho` column; fixed (conserved quantity is the authority); 120/120 states re-evaluate exit 0 | closed 09-15 |
| L19 | mass closure now reported in the certification and kept at 1 ulp by projecting each sweep's composition onto the conserved density (default on) | closed 09-15 |
| L20 | the `cool` column: negative as posed (the 0.35-0.76 ratio was a 2026-09-08 fixture IC); the identity is now asserted for `cool` as for `heat` | closed 09-15 |
| L21 | base boundary branch keyed on the cell-1 velocity (the P44 artifact) and the window 1e-6 above this planet's base Mach 4.6e-7: fixed (branch on the wind-window mass flux, window 1e-8); other regression bases byte-identical; **LHS campaign must be re-run** | closed 09-15 |

## 4. Source and binaries

- HEAD is still `43bc28c` ("minor bugs fixed"); everything since is
  uncommitted (73 files, about 9.6k insertions at the last count; the
  `LHS1140b/models/` tree, the memos and the plan are untracked).
- The tree binary `EXHALE.x` was md5 `a11038c24505` (built 08:23, L14 only)
  until 20:38, when it was rebuilt from the current source: md5
  `db87b88d1ce5`, carrying L4h, L7e, L7f, L18, L19, L21 (section 11). It is
  the campaign binary and must not be rebuilt before the campaign ends.
  Build with bare `make`
  (conda-forge gfortran 16.2 + its OpenBLAS); never `PATH=/usr/bin:$PATH
  make` (gfortran 13, wrong LAPACK; the 2026-09-14 incident).
- Private verification builds still on disk, each from the source as it
  stood for that item: `EXHALE_L4h.x` (`aaf762175e54`, comment-only
  rebuilds after the citation checks), `EXHALE_L7e.x` (`e4ebc751bcf1`),
  `EXHALE_L7f.x` = `EXHALE_L18.x` (`db87b88d1ce5`, carries L21),
  `EXHALE_L14.x` (`a11038c24505`). Object directories `build_L4h/`,
  `build_l14/`, `build_L7e/`, `build_L7f/`, `build_L18/`, `build_L17/`.
  Delete them at the end; the record is the memos.
- Files changed by the post-campaign items (this session): `steady_newton.f90`
  (L4h), `EXHALE_main.f90` (L14, L7e, L21 diagnostics), `diffusive_photochemistry.f90`,
  `System_HeH_mol.f90`, `certification.f90`, `molecular_seed_from_atomic_state.f90`
  (L7e), `mol_rates.f90`, `write_output.f90` (L7f), `load_IC.f90` (L18),
  `ionization_equilibrium.f90` (L19), `base_boundary.f90` (L21),
  `binary_element_diffusion.f90` (a "pipeline" word), tests in `src/tests/`
  (`carrier_retry`, `molecular_seed`, `run_mode`, `grid_and_gates`,
  `physics_probe`), `src/utils/map_state_to_grid.py`, `element_flux_closure.py`,
  `profile_match_level.py` (new), `LHS1140b/models/*.sh|*.py`.
- Regression matrix on `a11038c24505` (08:24-10:42, one thread): 7 pass,
  9 fail, all nine attributed (L6b cooling tables: the seven C/N/O cases;
  L7b/L7d: `mol_sec_ion`, `hp_trace_seed`). L19 adds `roundtrip` (0.39
  percent, a 1e-16 perturbation amplified by that fixture). L21 moves none
  of the non-LHS bases. The goldens have NOT been refreshed; the
  `wasp_full_newton/IC/` fixture is stale (written 2026-09-08 by an
  earlier code) and must be refreshed together with the goldens
  (`backup/regression/README.md` rule; command in the L19/L20 memo section 2.5).

## 5. Documents produced this series

- `docs/PLAN_20260913_lhs_stationary.md` (plan of record, items L1-L21),
  the 32 item memos, `docs/well_balanced_and_hydrostatic_20260913.{tex,pdf}`.
- `LHS1140b/MODELS.md` (catalog; sections 7-8 by `status.py --write`),
  `LHS1140b/models/README.md`, every case's `REPRODUCE.md`
  (`write_reproduce.py` now adds a "re-measuring this state" section, L18).
- `docs/lhs1140b_exhale_vs_pwinds.{tex,pdf}` (73 pp, results of 2026-09-14,
  humanized 2026-09-14; backup `.tex.before_humanize_20260914`);
  `LHS1140b/make_memo_figures.py` ported to `models/`.
- `docs/Update_EXHALE_stage2.md` section 11 (858 lines, 31 items;
  the golden paragraph is a placeholder `<!-- GOLDEN: ... -->`), and its
  generated `Update_EXHALE_stage2.{tex,pdf}` (191 pp).
- `docs/TO_BE_DONE.md`: a dated section for L4h, L7e, L4e proposals, L9,
  the `output_state_consistency` test row, the L6b coolant notes.
- `docs/restart_contract_design_20260909.md` section 6 and
  `docs/certification_tolerance_anchoring_20260910.md` anchor 7 (L18, L19).
- Papers obtained by the user this session, in `references/`: Kappeli 2014,
  2016; Schauer 1989; Johnsen 1980; Boehringer & Arnold 1986 (84, 1459;
  84, 2097 is unrelated); Black 1978; Kelley & Keyes 1998; Coffey et al.
  2003; Mulder & van Leer 1985; and Gropp et al. NASA CR-1998-208435 fetched
  from NTRS. Every citation these items make was checked against them.

## 6. Physics findings of the series (what a reader should take away)

1. The LHS wind certifies only with `Well balanced: True`: the O(1) base
   energy residual was the WENO3 truncation of the near-hydrostatic column
   entering the HLLC contact speed at Mach 1e-6 (L5c; the tex memo).
2. With the C/N/O coolants at the local electron density (L6b) the base
   metals are no longer a thermostat; the flux-closed ladder sits above the
   measured line at every He/H, the scalar branch crosses at He/H 1.49
   (K_zz 1e9).
3. H2 in the hot wind is made in place: He(2^3S) + H -> HeH+ + e,
   HeH+ + H -> H2+ + He, H2+ + H -> H2 + H+ (about half the metastable-He
   associative ionizations end as H2); it is not carried across the front
   (L7e section 10).
4. The wind network's base H2 chemistry and the lower-atmosphere handoff
   agree to 1.7 percent at 418 K (L7f); the He+ + H2 channels are now as the
   primary measurements have them (R20 retired, R17 two-mechanism form).
5. The molecular base at 1023 K is the ionization energy of H+ and He+
   deposited as heat when the molecular channels neutralize them; the ledger
   closes on the ionization potential once (L7g). Open: the vibrational
   heat fraction of R15.
6. The base boundary must key its entropy branch on a flux, not on the
   cell-1 velocity (L21); the LHS results of 2026-09-14 carry the interior
   isentrope at the base and are superseded.
7. The 0.02-XUV wall is the HLLC contact selection on a base subsonic to
   Mach 1e-8 (L17); the line is gone by 0.05 anyway.

## 7. Running at the time of writing (HISTORICAL, 20:30; section 11 is current)

- Worker on L7f/L7g: the `fixed5` continuation on the L21 boundary
  (`LHS1140b/models/.L7f/`), an old-boundary control, the five-case
  molecular regression (control and measured), a `mol_metals` control.
- Worker on L18/L19/L20/L21: the tail of the L19 regression matrices
  (`.L18/regm`, `.L18/regn`: `mol_ir_bands`, `mol_lyman_werner`,
  `mol_metals`, `mol_sec_ion`, `wasp_he23off`).
- No campaign case is running. `.L14/`, `.L17/`, `.L4h/`, `.L7e/`, `.L7f/`,
  `.L18/`, `.L21/` under `LHS1140b/models/` hold the experiments (the
  certified 0.03 ladder states are in `.L14/`).

## 8. The plan from here (in order; steps 1-5 were taken by 20:53, see section 11)

1. ~~User decision pending~~ APPROVED 20:35 ("재실행 승인"): the atomic
   part runs now; the nine molecular cases run alongside it on the same
   binary, from the certified HeH2.13 state (L7g closed, section 11).
2. Rebuild the tree binary with bare `make`; record the md5 in
   `MODELS.md` section 6 and in every new `REPRODUCE.md`.
3. `make check`; refresh the goldens once (the attributed movers of L6b,
   L7b/L7d, L19's `roundtrip`, and whatever L7f moves in the `mol_*`
   cases), refresh `backup/regression/wasp_full_newton/IC/` with them, run
   `make check` again, fill the `<!-- GOLDEN -->` paragraph of
   `Update_EXHALE_stage2.md` section 11 with the measured movements.
4. Preserve the current tree: `mv LHS1140b/models LHS1140b/models_20260914_preL21`
   (with its `REPRODUCE.md`s) and regenerate `models/` with
   `python3 make_models.py`; point `pick_seed.py`'s tiers at the preserved
   tree so every case continues from its own pre-L21 certified state (only
   the base changes, so a few passes each).
5. `./run_campaign.sh 8 8 'atomic_scalar*'`, then the nine closure rungs
   (`run_rungs.sh`), then `make_models.py --only atomic_photochem_gj1132x`
   and those seven; the three 0.01 cases stay "not solved (L17)".
6. After L7g: the nine molecular cases by continuation (HeH2.13 from the
   `.L7f` state, the others from it through the mapper), on the same binary.
7. `status.py --write`; regenerate the memo figures and the tex tables
   (the pending/"not solved" wording, Table 2 archived-vs-new becomes a
   three-column table: 2026-08-30, 2026-09-14, 2026-09-15); rebuild both
   PDFs; update `MODELS.md` sections 5-6 with the L21 reason for the
   re-run; `Update_EXHALE_stage2.md` section 11 gets the L15-L21 items and
   the re-run; memory files; `docs/code_status_20260910.md` pointer; README
   "Last updated" only on a push.
8. Delete the private binaries and object directories, and
   `EXHALE_L14.x`, `.L*` scratch trees that the memos do not name as records.

## 9. Open decisions and items for the user

- L7g's remainder: is the near-unity vibrational heat fraction of R15
  right at 3e13 cm^-3, and is the H2 band opacity right? (a physics item;
  the worker stops at the measurement).
- L15 (eight-thread reproducibility) is a prerequisite for calling the
  eight-thread campaign one reproducible solution set: a bounded study
  (a representative subset at 8 and 1 threads on frozen inputs and the
  frozen binary, comparing arrays, residuals, certification, Mdot and the
  synthetic line) is scheduled after the campaign; L9 (a read-only
  evaluation path) serves it. L16, the L4e proposals, L12 stages B-F: open.
- `backup/regression/arm*` case directories carry the banned noun; the
  goldens are keyed to the names (`docs/named_case_audit.md`).
- Whether `Johnsen`-based R17 should carry its three-body term at higher
  densities (five decades down at 1000 K; recorded at the code site).

## 10. Rules this session added to the standing ones

- Build only with bare `make`; a `PATH=/usr/bin` build ran one case for ten
  minutes on the wrong toolchain (2026-09-14).
- Kill by PID chosen by `cwd`, never `pkill -f` with a string that is in
  the shell's own command line (exit 144 three times).
- Replace running scripts atomically (`\mv -f` of a new file).
- Control/measured comparisons at one thread (L15).
- A worker's "verdict" is checked against the log before it is acted on;
  three of this session's briefs carried a wrong premise that the worker
  refuted by measurement (L14's energy row, L7e's implicit chemistry,
  L4h's earned floor), and the record keeps both the premise and the
  refutation.

## 11. Update, 2026-09-15 20:53

- **Tree binary rebuilt** at 20:38 with bare `make` from the current source:
  `EXHALE.x` md5 `db87b88d1ce53facf1d61084fa535ca5` (gfortran 16.2 conda,
  `make -q` up to date), identical byte for byte to the L7f/L18 verification
  builds. It carries L4h, L7e, L7f (R20 retired, R17 two-mechanism, Black
  channel), L18, L19, L21. **This is the campaign binary**; no `make` may run
  until the campaign ends (the source is still being edited by the L7g
  worker for diagnostics; a rebuild would change the campaign md5).
- **Regression matrix** running on a frozen copy of that binary
  (`REGRESSION_EXE=<scratchpad>/EXHALE_check_db87.x`, log
  `<scratchpad>/check2.log`, one thread); the goldens are refreshed once
  after it, together with `backup/regression/wasp_full_newton/IC/`.
- **L7g closed as diagnosis**: the reaction-heat ledger does not double
  count (every reaction heat is a difference of one formation energy per
  species; photoionization deposits hv - IP; the He+, H+ and Lyman-Werner
  cycles close on the ionization potential once; Schauer's -6.51/-9.16 eV
  reproduced as +6.5109/+9.1615 without fitting). The hot molecular base is
  the ionization energy of H+ and He+ deposited when molecular channels
  neutralize them. On the L21 boundary the molecular fiducial case
  CERTIFIES at outer pass 11 (ghost 226.02 K, cell 1 808 K, x2 0.355,
  mass row 4.9e-9, energy 1.0e-8); on the old boundary it never did.
  Named, not opened: the near-unity vibrational heat fraction of R15
  against the H2 band cooling (0.2 percent of the deposit; both IR channels
  on move the base by 0.2 K). `photochem` prescribes T and cannot be
  compared on this point.
- **Campaign re-run started 20:50:53**: the 2026-09-14 results were moved
  to `LHS1140b/models_20260914_preL21/` (inputs left in place),
  `pick_seed.py` got a tier 0 (the same case's own pre-L21 certified
  state), `CAMPAIGN_EXCLUDE='x0\.01'` drops the two scalar 0.01 cases,
  68 atomic prescribed cases at 6 x 8 threads (`models/campaign_20260915.log`),
  then the 9 closure rungs, then the 6 XUV photochem cases (the fiducial
  9.7 rung's new profile). Worker writes `docs/lhs1140b_rerun_20260915.md`
  (pre-L21 vs re-run table) and runs `status.py --write` at the end.
- **Molecular campaign**: `molecular_scalar_gj1132_kzz1e9/HeH2.13` is
  being re-solved in its catalog directory on the tree binary (outer pass
  1: mass 7.2e-9, energy 1.4e-8, tracking the scratch run accepted at pass
  11); the other eight launch automatically two at a time at 8 threads with
  `SEED=` that state (mapper continuation), `run_case.sh`'s molecular branch
  now takes `SEED=`.
- **L19 regression tail** (`.L18/regm`, `.L18/regn`): `mol_lyman_werner`,
  `wasp_he23off` still marching; `lower_profile` 3.3e-10, `mol_carrier`
  2.7e-11 already in the L19/L20 memo.
- Load: 72 cores, about 20 solver processes (campaign 6 x 8, molecular
  1-2 x 8, regressions at one thread).

## 12. Codex review of this handoff (2026-09-15 21:13), what was taken

`docs/session_handoff_20260915_review.md`. Each finding was checked
against the code before acting:

| finding | check | action |
|---|---|---|
| R23 (radiative charge transfer) deposits its whole 9.16 eV as gas heat | confirmed: `molecular_reaction_heat.f90:454` adds `k23 n(H2) n(He+) q(R23)` to `gamma_chem` with no photon term | reopened under L7g: photon leaves, product excitation thermalizes per the vibrational rule; every radiative channel of the network audited; molecular cases re-solved after the fix (their present runs serve as seeds); an atomic case checked byte-identical |
| L21 switch between the two discriminants is discontinuous; a distant-window flux is a stationary approximation used as a physical-time boundary law | confirmed: `base_boundary.f90` ~507-512 selects `M_wind` only above the window, else `M_i` | C1 blend of the two discriminants; wind-window keying restricted to the stationary route, marching keeps the local face; four tests; LHS branch unchanged (window far exceeded), verified bitwise |
| L18 density authority decided by the size of the discrepancy | confirmed: `load_IC.f90` ~884-912 | authority by intent: composition only where the loader itself transformed it (flags), else conserved density with a stated rounding allowance and refusal above it |
| L19 projection unbounded and unreported | confirmed: `ionization_equilibrium.f90` ~2870-2910 | maximum correction measured and reported; warning above rounding, refusal (sweep failure) above 1e-8; element budgets reported |
| rate commentary overstates (Black "upper bound"; R17 extrapolation "established"; three-body omission by a representative density) | confirmed in `mol_rates.f90` | wording to the stated assumptions; three-body guard on the local n and T |
| R15 heat fraction lacks a formation-state model | agreed; already named under L7g | stays open as the L7g remainder |
| pseudo-time CFL comment stated as general | confirmed at `steady_newton.f90` ~1609-1613 | reworded as the measured control it is; L17 stated as failure of the routes tried |
| L15 labeled nonblocking | agreed | relabeled a prerequisite; bounded study scheduled (section 9) |
| handoff: case rows summed to 101 | confirmed (the atomic row said 60; it is 56; the catalog has 95) | fixed in section 2 |
| "log Mdot -1.8 percent" | confirmed wrong wording | fixed: Mdot -1.8 percent, -0.0077 dex |
| sections 7-8 read as live instructions | agreed | marked historical |
| golden refresh should not follow attribution alone | noted; the standing rule of this project is one deliberate refresh at the end of a series with the movement reported, and the user decides | unchanged; the refresh waits for the R23/L21/L18/L19 fixes above so that one snapshot follows one source |
| a frozen source manifest should accompany the binary md5 | agreed | `LHS1140b/models/BINARY_MANIFEST_db87b88d1ce5.txt` written (md5 of every `src/` file at 21:20, with the list of files edited after the 20:38 build named as such) |
