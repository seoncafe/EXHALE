# Session handoff, 2026-09-15, revision 1 (22:00): the LHS 1140 b stationary series after the Codex review

Revision note. This is the 20:30 handoff (`session_handoff_20260915.md`)
carried forward with: the campaign re-run and the molecular campaign
running (section 11); the Codex review of 21:13 and what was taken from it
(section 12); the L7f/L7g close-out numbers (section 3); the corrected
case count and Mdot wording (sections 1-2); the source manifest of the
campaign binary; and one more paper to obtain (section 9). Sections 7-8 are
kept as the historical plan and marked so. Nothing here was re-measured
for the document.


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
| L7f | base H2 chemistry vs the handoff: the 7-21x was the seed's root formula (production depends on n(H2)); the two agree to 1.7 percent; R20 retired (Schauer 1989 bounds the HeH+ channel at <= 1e-14), R17 corrected from Johnsen 1980 and Boehringer & Arnold 1986, H2+ + He -> HeH+ + H added from Black 1978; golden movement attributed (HeH+ to nothing, He+/He2+ x2.2, hydro 1-3 percent; the `mol_metals` `cool` failure predates the item) | closed 09-15 |
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
  densities (five decades down at 1000 K; a guard on the local n, T is
  being added after the review).
- Paper to obtain: Moses & Bass 2000, JGR 105, 7013 (the source of R17's
  Arrhenius branch above ~700 K, where nothing is measured).
- The Codex review's structural points on L21 (a local face-consistent
  entropy condition rather than a distant-window flux) and on R15 (a
  formation-state model of the association energy) are recorded as the
  direction of the next items, not opened.

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

## 13. State at 21:19

- Atomic re-run: 6 certified, 0 failed so far of 68 (log
  `LHS1140b/models/campaign_20260915.log`); the closure rungs and the XUV
  photochem cases follow automatically.
- Molecular campaign: the fiducial case certified (log Mdot 7.91, EW 1.576
  %A, depth 5.63 percent, pass 11); the other eight run two at a time from
  its state; all nine are to be re-solved from their own states once the R23
  heat recipient is fixed (section 12).
- Review fixes in progress on private builds: R23 and the radiative
  channels (L7g reopened), the base-branch C1 blend and stationary gating
  (L21), the loader intent and the projection bounds (L18/L19), the
  pseudo-time comment (L4h). The tree binary stays `db87b88d1ce5` until the
  atomic campaign ends; the manifest is
  `LHS1140b/models/BINARY_MANIFEST_db87b88d1ce5.txt`.
- Regression matrix on the frozen copy of that binary still running; the
  golden refresh follows the review fixes, so that one snapshot follows one
  source.

## 14. Decision 2026-09-15 21:35 ("그렇게 진행")

The running atomic re-run is not stopped (no review fix reaches the atomic
fixed points). When the four review fixes have landed and the tree compiles
as one source: rebuild once, write its manifest, run the regression matrix
on a frozen copy, refresh the goldens once (with `wasp_full_newton/IC/`),
then **re-solve every catalog case from its own certified state on that
final binary** (atomic, closure rungs, XUV photochem, and the nine
molecular cases after the R23 fix) so the whole catalog stands on one
binary and one manifest; the present re-run supplies the seeds and the
pre-L21 comparison table. Then `status.py --write`, the memo figures and
tables, section 11 of the update log, the memories.
- 23:44: the regression matrix on the frozen campaign binary
  (`db87b88d1ce5`, one thread, five cases in parallel; log kept as
  `LHS1140b/models/.stopped/check_db87b88d1ce5_20260915.log`) fails every
  case against the goldens of 2026-09-10, `hydrostatic_column`'s hydro
  files excepted (9.6e-8). The largest movers are the ones the item memos
  attribute: the HeH+ column to nothing in every molecular case (R20
  retired, L7f), the C/N/O-cooled cases (L6b), the velocity sign flips at
  fronts (`wasp_full` row 209, `lower_profile` row 474, `mol_sec_ion` row
  164), and "1.0" entries that are zero-against-nonzero flags of the
  comparator. A whole-matrix attribution on the FINAL binary (control =
  each item's own measurement) is part of the golden refresh; no golden
  was refreshed.

## 15. Update, 2026-09-16 04:30

- Re-run on `db87b88d1ce5`: 64 of 83 atomic cases certified (55 prescribed,
  9 closure rungs, 0 failed), 13 low-XUV cases (x0.10-0.33 scalar and the
  six photochem-x) were grinding at dtau0 = 1 on the outer-cell energy row
  (L4f) and are restarted at `EXHALE_PTC_DTAU0=1.0e8` (the L14 recipe;
  approved 04:30); the 6 GJ 699 cases follow. Worker records:
  `docs/lhs1140b_rerun_20260915.md` (pre-L21 vs re-run, third column ready),
  `models/compare_trees.py`, `models/archive_results.sh` (written, not run),
  `pick_seed.py` tier 0 generalized to "the case's own most recent certified
  state", `write_reproduce.py` records the hostname.
- What moved (64 cases): cell-1 T -3 to -58 percent (median -45), density
  x1.03 to x2.31 (median x1.82), Mdot -0.2 to -2.9 percent (median -2.0),
  red EW median -1.9 percent (helium-poor rungs far more, e.g. kzz0/HeH0.55
  0.168 -> 0.031 %A). Crossings He/H: wellmixed 0.401 -> 0.404, kzz1e11
  0.878 -> 0.893, kzz1e10 1.141 -> 1.158, kzz1e9 1.492 -> 1.516, kzz1e8
  1.933 -> 1.973, kzz1e7 2.484 -> 2.537, kzz1e6 2.633 -> 2.893, kzz1e5
  2.707 -> 3.071, kzz0 2.717 -> 3.101; the closed ladder still has no bracket.
- A runner defect found and fixed on the way: `run_case.sh` carried only
  He/H onto a profile-case seed while `load_IC` compares every element of
  the `# reservoir` line; the rule now lives once in
  `src/utils/profile_match_level.py` and both the runner and the closure
  driver call it (recorded under L13).
- Molecular: kzz1e9/HeH2.13 and HeH0.55 certified with transit, HeH0.083
  and HeH9.7 (retry) solving, wellmixed 0.083/0.55 solving on lart3, the
  three remaining wait for their seeds. R23 fix moves the base by 1.9e-5.
- Machines: lart4 (this one) load ~49, lart3 ~48.

## 16. What is not going well (2026-09-16 07:10), in detail

Everything below is READ from the run logs named; nothing is a guess.

### 16.1 Four of the nine molecular cases do not certify on the present recipe

All four run the molecular recipe (`run_case.sh` molecular branch: seed,
`Well balanced: True`, `Secondary_ionization: Immediate`, stationary
route, `EXHALE_OUTER_PASSES=40`, movement bound on the particle count,
the L21 base boundary) on binary `db87b88d1ce5` or the L7f build
`42ec3a6c` (same source but for the molecular heat terms). The two
fiducial-group cases that certify (`kzz1e9/HeH2.13` at pass 11,
`kzz1e9/HeH0.55` with transit) show the recipe is sound near the
fiducial composition; the four that do not are at the composition
extremes or in the well-mixed group.

| case | seed | where it stands | reading of the trajectory |
|---|---|---|---|
| `molecular_scalar_gj1132_kzz1e9/HeH9.7` | certified molecular `HeH2.13`, mapper `--reservoir He/H 9.7` (a 0.66 dex step) | first attempt (lart4, `db87`): `info = 1`, `||R|| = 1.925`; second (lart3, `42ec3a6c`): outer pass 1 REFUSED, "the element composition relaxation found no admissible advance and its entry composition was restored, so no fixed point was approached"; a third attempt is in its first JFNK (`||R||` 0.91 to 1.48 over iterations 5-6, worst row the energy of cell 1) | the L5b shape: the element relaxation cannot take a first admissible step from this seed. A 0.66 dex composition step is beyond the 0.3 dex continuation rule the atomic ladder measured; the atomic `kzz1e9/HeH9.7` of the re-run (certified) is the seed to try through the atomic-to-molecular branch, once, before recording the case as not solved |
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | atomic `wellmixed/HeH0.083` (own pair, atomic-to-molecular branch) | pass 1: carrier row 4.79e-2 at cell 342, hydro within (mass 6e-6 falling to 4e-9 by pass 10); pass 10: 3.16e-2; pass 15-20: 2.5e-2 to 2.4e-2 **at cell 500** (the outer boundary cell); pass 24: REFUSED by the stagnation rule ("the joint distance has not fallen in 3 consecutive passes") | the carrier row moved from the interior (cells 342-348) to the outer boundary cell at pass 15 and stopped falling there (2.42e-2, 2.43e-2, 2.45e-2, 2.46e-2 at passes 21-24, slightly rising). Its index locates a maximum, not a cause: the outer face carries zero diffusive flux by construction but advective escape is allowed, and a last-cell residual can come from the advective flux difference, the diffusive-flux transition, chemistry, shielding, or an under-relaxed coupled update; that decomposition has not been measured (Codex review 1 of 09-16, R5) |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | atomic `wellmixed/HeH0.55` (own pair) | pass 1: 1.27e-2 at cell 233; pass 5: 7.59e-2 at 248 (rose); then 3.9e-2, 3.1e-2, 2.7e-2, 2.3e-2, 2.27e-2 at passes 10-30, the binding cell walking outward 258 -> 293; hydro within throughout | the index of the largest row moves outward (294 -> 299 over passes 31-35) while its value is flat then rising (0.0227 -> 0.0232); that is not evidence of a front passing, only of the largest row changing identity among several structures (Codex review 1 of 09-16, section 6.2). At 2.3e-2 against 1e-5 it will not certify inside the cap; a longer identical run is not justified without the row decomposition |
| `molecular_scalar_gj1132_kzz1e9/HeH0.083` | certified molecular `HeH2.13`, mapper `--reservoir He/H 0.083` (a 1.4 dex step) | pass 1: carrier row 7.8e-4 at cell 500; from pass 5 the binding row is the He/H elemental transport at cell 346: 1.18e-3, 1.30e-3, 1.13e-3, 9.87e-4, 8.94e-4, 8.35e-4 at passes 5-30; hydro within (mass 1.5e-9) | a slow fall of the elemental row at one cell that is not globally monotone (1.18e-3 at pass 5 rose to 1.30e-3 at pass 10 before falling); a constant 2 percent per pass from 8.35e-4 would need about 220 further passes, an extrapolation and not a forecast. The 1.4 dex composition step of the seed is a plausible contributor, not an established cause (the ladders were built with steps of at most 0.3 dex, an empirical continuation choice) |

Common threads: (i) the molecular recipe was proven on one composition and
is now applied with composition steps of 0.66 and 1.4 dex, against a
measured continuation limit of 0.3 dex; (ii) the two well-mixed cases bind
at or walk toward the outer boundary cell (500), where the carrier's outer
boundary has not been examined (L8 examined the wind's); (iii) every
"ionization roots" line of these runs reports "1 root(s) outside the
simplex" at step 1, which is the molecular clamp L7e section 5.1 records
as admissible and is not the cause. Decision in force (07:00): no blind
retries; each case gets one classified action (a 40-pass continuation only
where the binding row is still falling; the atomic-to-molecular seed once
for HeH9.7), and otherwise "not solved: <reason>" in its `REPRODUCE.md`
and `MODELS.md`. The catalog definitions are not loosened.

### 16.2 Thirteen low-XUV atomic cases needed the dtau0 = 1e8 start

On the corrected base boundary the first solve at `EXHALE_PTC_DTAU0=1.0`
of the x0.10-0.33 scalar-base cases and the six photochem-x cases burned
30-90 minutes per outer pass on the outer-cell energy row (cell 500, item
L4f) with the species row already at 1e-12; pre-L21 the same cases took 90
minutes end to end. Restarted at dtau0 = 1e8 (the L14 recipe, which
`run_case.sh` documents): `x0.20/HeH2.13` accepted at pass 34, photochem
x0.33 at pass 4. That the restart converges is evidence about these runs,
not proof that the explicit-stable interval caused the stalls; a
controlled comparison (same seed, binary, thread count, boundary; only
dtau0 varied; the shift against the local Jacobian scale recorded) has not
been made. The first attempts are kept in
`models/.stopped/<case>_dtau0-1_<ts>/`. Open: why the dtau0 = 1 stage,
which certified these cases on 2026-09-14, no longer converges on the
corrected base (the base density doubled and the cell-1 temperature
halved, so the explicit-stable interval and the base rows changed; not
measured beyond that).

### 16.3 The regression matrix on the campaign binary fails every case against the goldens of 2026-09-10

Section 13 (23:44 entry). Expected in direction (six items moved the
physics: L6b cooling tables, L7b/L7d molecular corrections, L7f rate
repairs incl. HeH+ to nothing, L19 mass projection, L21 base boundary),
but the whole-matrix attribution exists only per item (each worker's
control-versus-measured), never as one table on one binary, and two
signatures are not yet explained by any item: `hydrostatic_column`'s
`Ion_species` row 1 column 6 (a base ghost cell, zero against nonzero)
and the "1.0" flags of the comparator on zero-valued columns. The golden
refresh waits for that attribution on the final binary.

### 16.4 Eight-thread runs are not reproducible (L15)

**NOTE ADDED 2026-09-17: closed, no thread dependence found.** The traced
trajectory of the element reload, 105 outer iterations, 152 Krylov products and
923 recorded quantities, is bitwise identical at 1, 8 and 16 threads, and so are
both molecular configurations and the closing run on the catalog binary. The
2026-09-15 divergence is not reproducible on the present tree and is attributed
to nothing. `docs/lhs1140b_stationary_L15_20260916.md`; the account of the
series that closed it is `docs/session_handoff_20260917.md`. The section below
is left as written on 2026-09-16.


The JFNK stationary solve at eight OpenMP threads diverges from itself
after about 20 iterations (L4h memo section 5.0); one thread is bitwise.
The campaign runs at eight threads, so "reproducible" for a campaign case
means: the certified state re-evaluates within its certification (L18,
120/120 at one thread), not that the solve path repeats. A bounded study
(a representative subset at 8 and 1 threads, comparing arrays, residuals,
certification, Mdot and the line) is scheduled after the campaign. The cause is not located and
"an order-dependent reduction" is only one of the possibilities (a race, a
mutable cache, statefulness, or rounding amplification are the others);
the first test is repeated residual evaluations of one frozen state at one
and at eight threads, comparing the returned arrays before any iteration
history.

### 16.5 The three x0.01 cases (L17) and the closed ladder's lack of a bracket

Unchanged: at 0.02 of the fiducial XUV the base is subsonic to Mach 1e-8
with alternating-sign velocities and the HLLC contact selection puts the
base cells on a kink of the residual; the ladder reaches 0.03 on all three
branches. The Codex review's point stands that this is the failure of the
routes tried, not a proof. The flux-closed ladder sits above the measured
line at every He/H on the corrected base too (re-run rungs EW 2.20-2.22 %A
against 1.108 measured).

### 16.6 Housekeeping that is still open

- `backup/regression/wasp_full_newton/IC/` and `atomic_elem_newton/IC/`
  are pinned fixtures written on 2026-09-08 by an earlier code (mass
  closure 1.3e-3); three `grid_and_gates` assertions read them and fail on
  every binary; they are refreshed with the goldens.
- `backup/regression/arm*` case directories carry the banned noun; the
  goldens are keyed to the names.
- Private builds and object directories (`EXHALE_L4h.x`, `EXHALE_L7e.x`,
  `EXHALE_L7f.x`, `EXHALE_L18.x`, `EXHALE_L21b.x`, `EXHALE_L14.x`,
  `build_L*/`) and the scratch trees `LHS1140b/models/.L*` are to be
  deleted at the end except where a memo names one as a record.
- The tree binary has been rebuilt once during the series without a
  manifest of every private build's source; from `db87b88d1ce5` on there
  is one (`LHS1140b/models/BINARY_MANIFEST_db87b88d1ce5.txt`).

Addendum to 16.1 (07:20, from the molecular worker's own classification,
which agrees with the table above): `kzz1e9/HeH0.55` certified (log Mdot
7.64, EW 0.159 %A, depth 0.63 percent, 5 passes); `kzz1e9/HeH2.13`
re-solved on the R23-corrected build reproduces the db87 run to the fourth
digit (EW 1.5758 vs 1.5759). The stagnant well-mixed cases bind in the far
outer wind where H2 is a trace of a cold thin gas (cell 500 is 29.0 Rp at
285 K with x2 = 1.7e-2) and the worst cell of `wellmixed/HeH0.55` walks
outward pass by pass (233 -> 298) with the trust halved three times: which the first reading took for a feature drifting through the carrier's outer boundary; the second Codex review (section 6.2) is right that the logs show only the largest row changing identity, and the reading is withdrawn until the row is decomposed. `kzz1e9/HeH9.7`'s first-pass refusal came with hydrodynamic rows
of order one (mass 4.6e-2, momentum 0.78, energy 1.68): the seed carried
from HeH2.13 is nowhere near a solution; its one authorized attempt from
the re-certified atomic `kzz1e9/HeH9.7` is running. `kzz1e9/HeH0.083` gets
one 40-pass continuation (its row falls 2 percent per pass).

## 17. Second Codex review (2026-09-16 07:09, `session_handoff_20260915_review1.md`), what was taken

Each finding checked against the source before acting:

| finding | check | action |
|---|---|---|
| R1: the stationary and the physical-time routes use different base-boundary operators, so a certified state is certified against R_stat, not R_time | confirmed: that gating was my own decision after review 1 (to keep the marching goldens bitwise) | reopened: measure both constructions on one saved certified state (ghost rho, p, T, entropy source, base mass and energy flux, first residual rows); then one boundary for both routes, the flux-keyed branch wherever the wind window exists and the local face only where it does not (the guard already in the code), regression re-measured; the certification label stays "of the stationary operator" until the two agree |
| R2: `composition_closing` reads a bounded displacement (the carrier drift is the accepted change, limited by the movement bound), so a tighter bound can look like contraction | confirmed for the carrier half (`comp_drift` <- carrier drift after the bound); the element half already reads the undamped distance | the carrier's own fixed-point defect (the unbounded relaxation residual, or the residual change at fixed component scales before and after the coupled update) replaces the displacement in the progress test; logged beside it |
| R3: one global particle-count budget for the whole grid can let one region limit another; `worst drift cell`, `bound cell` and `worst residual cell` are three different indices | confirmed in `carrier_particle_count_change` and the trial loop | the three indices logged separately every pass; regional relaxation only if the bound is measured to bind at an unrelated location |
| R4: the carrier matrix omits the copied-ghost derivative at j = N when Frho(N) < 0 and at the inner ghost when no base composition is imposed | confirmed in the matrix assembly and `carrier_face_mass_fraction` | boundary derivatives derived from the same imposed-versus-copied ghost rule as the residual; directional-derivative probes for both signs and both branches; not claimed as the cause of any failure |
| R5: the outer-boundary reading of the well-mixed failures exceeds the evidence | agreed | wording withdrawn in 16.1 (above); the row decomposition at cell 500, its neighbor, the bound cell and the drift cell (both face fluxes, production and loss, the row scale, the accepted interval, T, radiation, boundary flux sign) is the next measurement, before any domain or boundary change |
| R6: R23's post-photon remainder is deposited at once; a validity criterion is needed | agreed; already named under L7g | stays open with the R15 formation-state question; an effective thermalization fraction only with its timescale conditions stated |
| R7 and section 9: "order-dependent reduction" too specific; test frozen-state residuals first | agreed | wording corrected in 16.4; the frozen-state residual test at 1 and 8 threads is the first step of L15 |
| section 6.4: "monotone 2 percent per pass" is not monotone and the extrapolation is 220 passes | agreed | corrected in 16.1 |
| section 8: the dtau0 = 1e8 recovery is evidence about those runs, not a cause | agreed | corrected in 16.2 |
| section 9: a zero-versus-nonzero ghost value needs its boundary and loader semantics checked, not a tolerance | agreed | part of the golden attribution on the final binary |

The review's section 10 order of experiments (same saved state under both
boundaries; row decomposition of a frozen well-mixed state; frozen-state
residuals at 1 and 8 threads; boundary directional derivatives; the
composition-bound comparison; own-abundance versus mapped seeds; domain
sensitivity last) is adopted as the order of the next work, split between
the two workers that own those files. No pass budget is raised, no
tolerance loosened, no boundary row suppressed, no rate changed to obtain
convergence.

## 18. Update, 2026-09-16 10:25

- Review-2 items closed by measurement: R1 (one boundary for both routes;
  the window's standing read from its own flux spread, threshold 1e-2 set
  by the measured gap between 147 certified winds (<= 2.95e-3) and every
  non-wind fixture (>= 7.5e-2); `wasp_full`, `hydrostatic_column`,
  `hp_front` identical to the pre-L21 boundary; the certified fiducial and
  the L14 fixture bitwise; `base_branch` 6/6), R2 (carrier progress read
  from its own fixed-point defect; three false stagnation counts removed;
  verdicts unchanged), R3 (drift, bound and residual cells logged
  separately; measured distinct), R4 (two boundary derivatives and a factor
  m_c in every advective carrier-Jacobian entry fixed; new suite
  `carrier_boundary_jacobian` 10/0), L15 first step (frozen-state residual
  evaluation md5-identical at 1 and 8 threads: the thread dependence is not
  in the residual). The certification label stays "against the stationary
  operator" (`restart_contract_design_20260909.md` section 7) although the
  boundary is now one.
- Campaign: 79 of 83 atomic cases certified on `db87b88d1ce5`; still
  solving `x0.10/HeH9.7`, photochem x0.15 and x0.20; photochem x0.25 ended
  info = 1 twice; the four hold the cell-1 mass row at ||R|| about 1.5 on
  the colder, denser base (a base obstruction of the L17 kind at higher
  XUV than before; not yet measured). Pre-L21 versus re-run table in
  `docs/lhs1140b_rerun_20260915.md`.
- Molecular: 2 certified, HeH0.083 in its one continuation, the
  well-mixed row decomposition pending, three cases waiting for seeds.
- Next: the molecular worker's last diagnostic edit, then the final
  rebuild, its manifest, the matrix with per-item attribution, the golden
  refresh with the two pinned IC fixtures, `archive_results.sh
  models_20260915_db87`, the final pass of every case on the final binary,
  `status.py --write`, the memo tables and figures, `Update_EXHALE_stage2.md`
  section 11, the memories, and the cleanup of private builds.

## 19. Final binary, 2026-09-16 11:05

`EXHALE.x` rebuilt from the finished source: md5
`c2e9c9990b9f14f1be8cd77abca68945` (identical to the molecular worker's
last build), manifest `LHS1140b/models/BINARY_MANIFEST_c2e9c9990b9f.txt`
(256 source files at build time). The regression matrix runs on a frozen
copy of it; its per-item attribution and the golden refresh follow.
Decisions taken at the same time: (i) the two `molecular_photochem`
cases are redefined on the re-run closure rungs' own converged columns
(He/H 2.0924 and 9.0103, the rungs 2.09 and 9) instead of the stored
2026-08-30 profiles, so that each sits on the same column as its atomic
counterpart and can be seeded from its certified state; recorded in
`MODELS.md`; (ii) the well-mixed molecular failures are the one scalar
movement bound over a column that holds a slowly moving H2 front and a
far wind that could converge at once (measured: every relaxation of both
runs ended on the bound, the binding cell at the front, the refusing cell
in the far wind); a regional or coupled relaxation is a new solver item
(L22), proposed and not opened; the three well-mixed molecular cases are
recorded as not solved with that reason.

## 20. Final pass done, 2026-09-16 12:10

- Atomic: **79 / 79 certified on `c2e9c9990b9f`** (67 prescribed in 5 min at
  24 s each, 9 closure rungs k = 4 in 40 min on lart3, 3 XUV photochem in
  43 s), every one in a single outer pass from its db87 state; the db87
  results are preserved in `LHS1140b/models_20260915_db87/` (83 cases,
  413 MB) beside `models_20260914_preL21/`. db87 against final over the
  79: red-pair depth bit-identical 79/79, log Mdot to 8.7 digits, base T
  and rho to 8.3, EW to 7.6 digits (the largest departures are closure
  rungs, where the alternation stops at a marginally different point of
  one fixed point). `docs/lhs1140b_rerun_20260915.md` carries the
  three-column table and an `agreement` block from `compare_trees.py`.
- Not solved (7): x0.01 x3 (L17) and `x0.10/HeH9.7`, photochem
  x0.15/0.20/0.25 (hydrodynamic rows under a still-travelling composition,
  the same at dtau0 = 1e8 as at 1).
- Molecular: certified on the final binary `kzz1e9/HeH2.13`, `HeH9.7`
  (`HeH0.55` in its final pass); `photochem/HeH2.09` and `HeH9` running on
  their own re-converged rung columns; `kzz1e9/HeH0.083` continuation
  rising toward its cap; the three well-mixed recorded not solved (one
  scalar movement bound; L22 proposed).
- In progress: the regression matrix on the final binary (then per-item
  attribution and the golden refresh), the memo update to the final
  results, section 11 of the update log through L22 and the re-run.

## 21. Goldens refreshed and the tree cleaned, 2026-09-16 14:40

- Regression: the 2026-09-10 goldens preserved as
  `backup/regression/golden_blockJ_20260910/`; all 16 default cases
  re-snapshotted from the final binary `c2e9c9990b9f`; `run_check.sh check`
  at one thread: **PASS, 64/64 files byte-identical**. Attribution table,
  the two signature verdicts (He III at 1e-180 cm^-3 in the base ghost is
  L19's projection, proved bit for bit with `EXHALE_MASS_PROJECTION=0`;
  the "1.0" columns are the integer `adv_*_status` columns following the
  physics) and the timeline point (the goldens predate L6b, so the metals
  cases move most) are in `Update_EXHALE_stage2.md` section 11's golden
  block. `wasp_full_newton/IC/` refreshed (the three `stationary_evaluate_*`
  rows now pass; `grid_and_gates` 205/1, the one being the pre-existing
  one-bit stagnation row); `roundtrip` (outside the default cases) is being
  refreshed on the same decision; `atomic_elem_newton/IC/` awaits its
  70k-step re-run.
- Tree: the sixteen private binaries (md5 record in
  `LHS1140b/models/.stopped/private_builds_deleted_20260916.txt`) and the
  eleven private object directories deleted; only `EXHALE.x` and `build/`
  remain, `make -q` up to date. The scratch trees `LHS1140b/models/.L*`
  stay where a memo names them as records.
- Still running: the two `molecular_photochem` cases on lart3 (final
  binary), the `atomic_elem_newton` fixture re-run.
- 15:50: `roundtrip` golden refreshed (4/4 identical; its old `_adv`
  goldens had 7 and 38 columns against the 10 and 40 the code writes).
  Two findings from that refresh, recorded as plan items: **L23**, the
  loader does not restore `cert_reason` (one-line fix, deferred so the
  catalog and goldens keep one binary); **L24**, the `atomic_elem_newton`
  fixture no longer stands still (N7's obstruction is gone) but ramp-collapses
  at iteration 60, a case L4h did not measure; its 2026-09-09 pin is kept.
  No compute process of the golden worker remains.
