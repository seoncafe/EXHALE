# EXHALE: state of the code, what physics is in and verified, where the problems begin, and what comes next (2026-09-05)

Written to restart from in a fresh session. Read this first; then
`docs/session_handoff_20260905.md` (the last session's runs and scratch
paths), `docs/Update_EXHALE_stage1.md` sections 159-171 (the dated record), and
`docs/koskinen2022_model_a_comparison.tex` (the current validation target).
Every statement below was verified in the sessions of 2026-09-02 to 09-06
unless marked "per the earlier record" (memory notes of earlier sessions).

## 0. The one-paragraph verdict

**The atomic line -- H, He (with the 2^3S metastable) and the ten trace
metals, with Roche or spherical gravity, the two-stage march and the Newton
finish, and the transit post-processing -- is implemented and validated**
(Huang et al. 2023 WASP-121b reproduced in mass-loss rate, metal cooling,
n=2 physics and transit depths; Newton-converged on four planets), with
four qualifications listed in section 3.1. **The molecular line -- H2,
H2+, H3+, HeH+, the oxygen cycle, the molecular carriers and the transported
ionization state -- is implemented, regression-guarded, and only partly
validated**: against Koskinen et al. (2022) Model A it reproduces the
velocity, the neutral densities and the layer temperature within the
plotted points, but the outer wind's temperature and ionization hinge on
a heating rate the published paper does not let us reproduce, and its
lower layer has no reachable steady state in the present code. The
problems begin exactly there (section 4).

## 1. The tree

* Branch `v1.00`, HEAD `6d07d48` (2026-09-02). **175 modified or new files
  are uncommitted**; nothing since HEAD has been committed or pushed.
  `backup/` (regression harness, goldens) is not in the git remote.
* Toolchain: plain `make` builds with the conda-forge gfortran 16.2 on PATH
  and the OpenBLAS of the same prefix (rpath'd); the compiler's path and
  version are part of the rebuild stamp. `make FC=ifx OBJDIR=build_ifx
  EXE=EXHALE_ifx.x` builds with ifx 2026.1 + MKL (verified: the converged
  atomic cases within 3e-7 of the gfortran outputs, unconverged molecular
  snapshots up to 2 per cent in trace species, mass-loss rates identical).
  Do not prepend `/usr/bin` to PATH any more.
* Regression: `backup/regression/run_check.sh check` (12 named cases: 11
  defaults plus `roundtrip`); bitwise first, else within
  `REGRESSION_REL_TOL` (1e-3); `REGRESSION_EXE=<binary>` runs another
  binary. **The goldens are the pre-section-170 snapshots**; section 170 is a
  physical fix that moved the molecular snapshots by up to 2 per cent, so
  all twelve goldens must be re-snapshotted (the case directories hold the
  post-170 outputs): `cd backup/regression && ./run_check.sh golden <12
  cases> && ./run_check.sh check`. The agent's attempt was refused by the
  auto-mode classifier; the user runs it. Everything after 170 (section 171)
  left the default path bitwise identical (verified in scratch copies).
* Documentation: `docs/Update_EXHALE_stage1.md` (26,000 lines, sections 1-171),
  104 memos in `docs/`, `docs/input_schema.md` (every input key),
  `README.md`, `README_HOWTO.md`, `docs/EXHALE_user_manual.tex`. Whether the
  manual and `Update_EXHALE_stage1.tex/.pdf` are current against the `.md`
  changelog was not checked this session.

## 2. Physics inventory

Legend: **V** verified against a published solution or an analytic limit;
**I** implemented and regression-guarded, not validated against an
external result; **L** implemented with a known limitation stated in the
row; **X** not implemented.

### 2.1 Hydrodynamics and the steady state

| piece | status | evidence / limitation |
|---|---|---|
| Godunov (HLLC) + RK3, PLM and WENO3 reconstruction, Mixed grid, runtime cell count | V | ATES heritage; the regression matrix |
| gravity: Roche potential truncated at L1 (default) or spherical to a stated radius (`Domain mode: Spherical`, `Outer radius [R_p]`) | V | Huang 2023 Case B (Roche, per the earlier record); Koskinen 2022 Model A velocity and neutral densities (spherical, this session) |
| lower boundary: fixed rho or p at the base, isothermal ghost, base velocity from the mass flux; characteristic face condition (section 152) | L | the molecular base "breathes": a marching limit cycle at the default CFL (section 161: absent below CFL 0.15-0.3). Harmless to the steady state, but it keeps `du` from converging in molecular runs and distorts fixed-step snapshots |
| outer boundary: outflow | I | -- |
| viscosity, conduction (Watson et al. 1981 kappa) | I | keys `Viscosity`, `Conduction`; conduction negligible on Model A, as they state |
| caloric EOS: H2 rovibrational ladder, atoms/ions/electrons 3/2 k; `Caloric EOS: monatomic` (comparison option) | I | P53; the Koskinen gate's outer temperature is insensitive to it (20 K) |
| convergence: two-stage `du` (flux spread) hand-off, JFNK Newton finish with the scaled trust region (sections 153/159), PTC | L | atomic winds: Newton info = 0 on four planets (per the earlier record); `wasp_full_newton` stops at `||R||` ~5e-4 against a 1e-5 target (base energy row). A `du` stop alone leaves Mdot path-dependent at the percent level: **always Newton-finish an atomic run before quoting Mdot.** Molecular winds with transported H+ cannot be Newton-finished (2.4) |
| restart: write -> read -> write is the identity (`roundtrip` case); restart onto another grid by `interp_ic.py` (session scratch) | I | section 169.3 fixed the ghost-row inconsistency of the writer; files written before it carry it in their two lower ghost rows |
| OpenMP over cells (ionization sweep, cooling, hydro stages, carrier Jacobian) | V | 16 threads: 4x a 2026-09-04 step; serial-vs-parallel bitwise for the atomic cases |

### 2.2 Radiation and heating

| piece | status | evidence / limitation |
|---|---|---|
| spectrum: power law, monochromatic, or a loaded SED (`Spectrum type: Load` + file: two columns, A and erg cm^-2 s^-1 A^-1 at the planet, no header) | V | the mean-solar file of the Koskinen comparison; a molecular run with a loaded SED crashed at startup until 171.2 (an unallocated H2 channel vector) |
| dayside average `2D approximate method: Rate/4 + Mdot` (rates and heating / 4, full-sphere Mdot) | V | the same as Model A's "flux divided by 4" (section 168) |
| photoionization H, He, He+, He 2^3S, metals (Verner et al. 1996 and Verner & Yakovlev 1995; opacity tables optional) | V | ATES/AIOLOS heritage; Huang 2023. The seventeen metal stages take each fit only where it is valid and carry the SHELL SUM, outer plus subshells, reproducing Verner's own `phfit2` to 1e-10 at 17 ions x 12 energies (MEASURED). Huang et al. 2023 (ApJ 951, 123, section 2.4.1) also combine outer- and inner-shell cross sections, which supports checking a shell-complete opacity and does not by itself validate the effective charge-stage and heating approximation of the next row, which keeps its applicability statement and its own charge/energy audit (PLAN_20260909_review, section 10) |
| metal inner-shell absorption: the ion advances one stage, one electron of h nu - E_th(outer) to the cascade | L | the true Auger event leaves the atom two or more stages up, which the three-stage metal ladder cannot hold; the one electron delivers 0.95-1.09x the heat and 1.10-1.37x the secondary H I ionizations of the two the real event has (MEASURED, C I 300 eV, N I 530 eV, O I 550 eV). Fluorescence neglected, and quantified: with e_top = 1240 eV only C, N, O and Na I can have a K hole and their omega_K are 0.28, 0.52, 0.83 and 2.3 per cent, while every reachable L-shell event is radiationless to better than 0.7 per cent (READ, Krause 1979 Tables 3 and 5). What the published photochemistry models do (Cecchi-Pestellini et al. 2009; Locci et al. 2022); decision item 7, decided (b) 2026-09-07 |
| H2 photoabsorption: four channels (H2+; H+ + H; 2H+; neutral dissociation window) sharing one cross section (section 151) | I | moved molecular snapshots 2-24 per cent (151.11); no external check |
| photoelectron energy: h nu - I heats (default); secondary ionization (Shull & van Steenberg, staged after first convergence, or `Immediate`, or off) | V/L | standard for the atomic line; the Koskinen comparison shows Model A deposits 2.5-3x h nu - I at the same ionization rate -- a property of their code. `Photoelectron heating: full` or `<fraction>` are comparison-only options (171.5) |
| He recombination photons ionizing H (`He_rec_coupling`, default on) | I | off for the Model A runs (they lose I to recombination) |
| H(n=2): Balmer photoionization, Ly-alpha pumping, in-situ Ly-alpha emissivity from LaRT | V | Huang 2023 phases 3a/3b (per the earlier record) |
| Lyman-Werner H2 photodissociation with a Cloudy-derived self-shielding table, fluorescence heating (sections 122/150) | I | default off unless a flux is stated or a spectrum is loaded; a stated zero switches it off (171.2) |

### 2.3 Ionization, chemistry, cooling

| piece | status | evidence / limitation |
|---|---|---|
| coupled H/He (+ He 2^3S, + metals, + molecules) local-equilibrium solve per cell (`System_HeH_*`, analytic-Jacobian Newton, hybrd1 second attempt) | V | Huang 2023 phases 2-3 |
| rates: Badnell case B recombination, Voronov collisional ionization; `Atomic rate set: Koskinen2022` (their Table 1 R1-R4) | V | the switch moved the Koskinen gate by < 5 per cent |
| metals C, N, O, Mg, Si, Ca, Na, K, S, Fe; closed-form CHIANTI cooling; `eos_metals` in the mass/electron/particle budget | V | Huang 2023 WASP-121b Fe II / Mg II cooling and Mg II, Na I, Ca II transit depths (per the earlier record) |
| He 2^3S: Penning ionization (Taylor 2025 fit), triplet cooling in the energy solver, charge exchange, He 10830 transit | V | Falorca 2026 review items closed (per the earlier record) |
| molecular network H2, H2+, H3+, HeH+ (Koskinen 2022 Table 1 R1-R23); collisional reaction heat (P53, default on); H2 caloric EOS | I | the Koskinen gate; the H2 loss rate of one molecule at 1.5 r_base agrees with Model A (1.0e-5 vs 1.2e-5 s^-1) |
| oxygen cycle OH, H2O, CO, four FUV photolysis bands, CO reservoir, `Oxygen_chemistry.txt` | I | no external check |
| CO destruction: `He+ + CO` (UMIST RATE22 4068, with its He+ sink in the molecular systems' He+ row since B3b-CO2) and shielded photodissociation on the Lyman-Werner beam (Visser et al. 2009 Table 6), with a domain record evaluated in each cell (2026-09-06, B3b-CO) | I | the shielding table's own spread between its 5 K and 50 K excitation temperatures is a factor 22 at the worst grid point, and the layer runs at 800-3000 K. The thermodynamic CO ceiling that stood above the layer was deleted on 2026-09-06 (item CEILING-DEL); MEASURED at the deletion, it never fired in the current tree, so every output column was byte-identical across it |
| cooling: recombination, collisional excitation (Ly-alpha and He), bremsstrahlung, H3+ (Miller 2013 LTE x Table 6 non-LTE), H2/H2O/CO infrared bands, the infrared field of the atmosphere below | I/L | H3+ per molecule is ~3x what Model A's Figure 9 implies at 3 r_base (their Koskinen 2009 detailed-balance table is not reproduced here) |
| the sweep returns the heating and cooling OF THE COMPOSITION IT RETURNS (section 170) | V | fixed 2026-09-05; before, both belonged to the pre-sweep composition |

### 2.4 Transport of the composition

| piece | status | evidence / limitation |
|---|---|---|
| ATES-style advection correction of the ionization state in post-processing (`_adv` profiles, atoms and metals; feeds the transit tool) | V | Huang 2023 transit depths (per the earlier record); this is how the atomic line treats advected ionization |
| carrier transport operator in the hydro: implicit advection-diffusion of H2 (and OH, H2O, CO), 2nd-order limited, base inflow, Blanc's-law molecular diffusion, eddy `Kzz` (sections 137/158) | I | Picard alternation with the wind |
| `Ionization transport: True`: H+ as a fifth carrier, the sweep handed the transported fraction (section 169) | V | reproduces the advected estimate where P r/v < 1; the Koskinen gate's n_e excess fell from 5-10x to 1.5x. **Molecular branch only** (refuses a run without molecular chemistry) |
| He/H binary element diffusion (`He_diffusion`) | I | Phase 1 (per the earlier record); off in the Koskinen runs (their He/H is not depleted) |
| **He+, H2+, H3+, HeH+ (and every metal stage) are local equilibrium in the hydro** | X | where the Koskinen comparison's remaining ionization excess sits: He+ 3x theirs above 1.5 r_base -- a local root overstates an ion where P r/v < 1, as it did for H+ |
| **a steady solve that keeps a transported H+** | X | the proton row is not a Newton unknown; `Ionization transport` refuses `Solver: Newton`, `EXHALE_PTC` and `Coupled carrier solve`. Consequence: the molecular layer's steady state (time scale 1e7-1e8 s) is unreachable by marching and cannot be closed by Newton either |

### 2.5 Lower atmosphere and hand-off

| piece | status |
|---|---|
| Koskinen 2022 analytic column (Tier 1); `base.inp` scalar handoff (p_base, T_base, HeH_base, Kzz_base, q_H2_base); lower-atmosphere PROFILE handoff from VULCAN/photochem (Tier 3) | I (per the earlier record); the Koskinen runs use the scalar handoff at 1 microbar |

### 2.6 Post-processing

| piece | status |
|---|---|
| transit spectra (He 10830, Ly-alpha, H-alpha/beta, Mg II, Ca II, Na I), Roche geometry, LaRT coupling (`EXHALE_transit.py`) | V (Huang 2023 phase 5, per the earlier record). Every transit depth computed before 2026-09-03 was 0.2-3.8 per cent too deep (an outer ghost row was counted); the paper/poster tables are stale for that and for the defaults changed since |
| `_adv` profiles: molecular columns | L | the equilibrium sweep's values (no molecular post-process state) |

## 3. Validation status, by target

### 3.1 The atomic + metals line: implemented and validated, with four qualifications

Reproduced (per the earlier record): Huang et al. (2023) WASP-121b -- Mdot
(0.056 / 0.337 for Cases A / B against 0.052 / 0.32), the Fe II / Mg II
cooling structure, the n=2 proton source and Balmer heating, the Mg II,
Na I, Ca II transit depths; Newton-converged (info = 0) winds for HD 209458b,
HD 189733b, WASP-52b, WASP-121b. The two physical fixes of this session
(169.3, 170) moved the converged atomic states by at most 2e-4.

1. **Convergence depth.** `wasp_full_newton` reaches `||R||` ~5e-4, not the
   1e-5 target (base energy row). Answers are stable; they are not
   "precisely converged". A `du`-only stop is path-dependent at the
   percent level in Mdot.
2. **Advected ionization is post-processed, not in the hydro**, for the
   atomic line (ATES method). Where P r/v >> 1 in the wind -- the hot
   Jupiters validated -- local equilibrium is right; at low flux or large
   radius the same limitation the molecular line met applies, and the
   in-hydro `Ionization transport` exists only with molecular chemistry on.
3. **The paper and poster tables** predate several defaults (staged
   secondary ionization, He recombination coupling, triplet cooling,
   section 170) and the ghost-row transit correction; re-convergence is
   instruction-gated.
4. The photoelectron accounting (h nu - I) is the standard one and is what
   the validation used; the "heat per ionization" question of the Koskinen
   comparison is about their code, not a defect of this line.

### 3.2 The molecular line: implemented, partly validated

Target: Koskinen et al. (2022) Model A (every condition checked one by
one, memo section 2.1). Within the plotted points: velocity, H2/H/He
densities, the layer temperature (1.05-1.30, within 100 K after 4e5
steps). Decided by the heat deposited per ionization: the outer
temperature, velocity and ionization -- their Figure 9 has 2.5-3x the
h nu - I rate at the same ionization rate; with a stated 0.55 of the
photon energy (comparison only) the heating and temperature match, and
n_e stays 1.4-1.6x, H3+ 0.35-0.5x, He+ 3x, f(H2) 0.06-0.08 low. Mdot
0.66-0.78 of theirs, from the unsettled layer. Memo section 4.1 lists all
ten differences with sizes.

## 4. Where the problems begin (ordered by size and depth)

1. **No steady state for a molecular layer with transported H+** (design).
   Marching cannot reach it (1e7-1e8 s at 2 s per step); Newton cannot
   hold it (the proton row is not an unknown, so the finish returns the
   composition to its local root, which section 168 showed is wrong where
   P r/v < 1). Blocks: a converged Model A comparison, a converged
   hot-Uranus Mdot, any molecular case that is not a fixed-step snapshot.
2. **Local-equilibrium ions in a transported wind** (He+, the molecular
   ions; for the atomic line every stage): the same error the H+ had
   before section 169; fixable with the existing carrier operator.
3. **Heating per ionization, external.** Model A's heating profile is not
   what its stated accounting and spectrum produce. EXHALE keeps h nu - I;
   the comparison is a bracket until the authors' spectrum and heating
   implementation are known. Do not tune the fraction silently.
4. **Base breathing / limit cycle** at the default CFL: numerics, not
   physics; it spoils fixed-step molecular snapshots and `du` convergence.
5. **H3+ non-LTE cooling** differs from the published detailed-balance
   treatment by ~3x per molecule at wind densities.
6. Smaller: the spectrum file has nothing below 5 A; `_adv` molecular
   columns; pre-169 snapshots' ghost rows; the manual and the `.tex`
   changelog may lag the `.md`; 26 callerless `Cool_coeff` wrappers and
   nine `_func` names; the H I / He I literal constants not unified.

## 5. Tasks, in priority order, with acceptance criteria

1. **Goldens (user).** Refresh the twelve goldens for section 170 and
   re-check. Acceptance: `run_check.sh check` all identical.
2. **Proton (then carrier) rows as Newton unknowns.** Extend the coupled
   steady system (`Coupled carrier solve` already carries n(H2)) with the
   H+ continuity row -- advection plus the same `mol_heh_rows` source the
   carrier operator uses -- so `Ionization transport` and the Newton
   finish coexist. Acceptance: (a) a Newton-finished Koskinen gate whose
   x(H+) equals the marching transported value where P r/v < 1 and the
   local root where P r/v >> 1 (limit test on `examples/18_oxygen_chemistry`);
   (b) the layer closed to `||R||` 1e-5 with the base mass flux equal to
   the wind's; (c) atomic cases bitwise unchanged with the key off.
3. **Transport He+ (then H2+, H3+, HeH+) as carriers**; `ic_Hp` is the
   template (section 169). Acceptance: He+ on the Koskinen gate within the
   plotted points (theirs 1e4, flat); element census closed to 1e-10;
   key-off bitwise. Then decide whether the atomic line gets in-hydro
   ionization transport too (today it refuses without molecular chemistry).
4. **Ask the authors** (user decision) for the mean solar spectrum file
   and the heating deposited per ionization, or Figure 9 in numbers.
5. **H3+ non-LTE factor**: the Koskinen et al. (2009) Table 2 ratios as an
   alternative table (a key selects); measure the outer temperature
   response. Acceptance: cooling of each molecule at 3 r_base within 1.5x of
   what their Figures 8 and 9 imply.
6. **Base limit cycle**: local time-step control in the base cells (CFL
   below 0.15-0.3 there) or the characteristic face condition revisited.
   Acceptance: `mol_sec_ion` base velocity free of the odd-even pattern;
   goldens refreshed once.
7. **Convergence depth of the atomic Newton finish**: the base energy row
   that holds `wasp_full_newton` at `||R||` 5e-4. Acceptance: 1e-5 reached
   on the four planets; Mdot movement reported.
8. **Housekeeping**: `interp_ic.py` under `src/utils/`; the `Cool_coeff`
   wrappers and `_func` names; the literal constants; regenerate
   `Update_EXHALE_stage1.tex/.pdf` and the manual from the `.md`; a spectrum file
   with the 1-5 A part.
9. **Instruction-gated (never autonomous)**: paper/poster re-convergence
   with the current defaults, the LHS 1140 b write-up, the GitHub public
   switch, P36.

## 6. Conventions and pitfalls worth remembering

* A comparison run states every condition (global CLAUDE.md rule of
  2026-09-05) and is judged on plotted profiles, not on two numbers.
* Physics defects are fixed when found; goldens follow (refresh, report the
  movement); sub-0.1 per cent technical movement needs no refresh.
* Comparison-only options (`Photoelectron heating: full|<f>`, `Caloric EOS:
  monatomic`, `Atomic rate set`) are labeled so in the setup report and in
  `EXHALE_resolved.out`; a run with one must say so in any report.
* A stated `Stellar LW flux: 0.0` turns the band off even with a spectrum
  file; an absent key with a spectrum file integrates the band.
* `Ionization transport` runs are marching only; `du_th` may stay at its
  normal values (a `du` stop only ends the march); do not add `Solver:
  Newton` to them.
* Kill by PID after checking `/proc/<pid>/cwd`; never `pkill`. `\cp -f`,
  `tail -n N`; `$PWD` inside `( cd dir && ... )` is `dir`; the interactive
  `cp` alias hangs a script.
* The auto-mode classifier refuses multi-case golden refreshes and compound
  commands that include one.

## 7. Where things are

* Koskinen comparison runs: `benchmarks/koskinen2022_model_a/`
  (`matched_hnu_minus_I/`, `matched_full_hnu/`, `matched_fraction_0p55/`,
  the spectrum files with READMEs, `model_a_fig9_digitized.txt`).
* Figures and table: `docs/k22_model_a_figures.py` (up to three runs;
  `K22_LABEL1/2/3`, `K22_FIG_DIR`).
* Digitized published profiles: `docs/p23_published_profiles.md`.
* Session scratch (volatile): `/tmp/claude-1000/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/4b5720ff-a4b5-4ae7-b461-c2d91dc0a3af/scratchpad/k22_matchA/`
  (`interp_ic.py`, `snap.sh`, the A-G run directories).
* Memory index: `~/.claude/projects/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/memory/MEMORY.md`.
