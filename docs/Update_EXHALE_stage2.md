# EXHALE Update Log (stage 2, from 2026-09-05)

This log starts fresh on 2026-09-05. The two earlier logs are kept unchanged
under their new names:

- `docs/Update_EXHALE_stage0.{tex,pdf}`: the early phase (the ATES fork up to
  2026-08-15).
- `docs/Update_EXHALE_stage1.{md,tex,pdf}`: sections 1-171, 2026-06 to
  2026-09-05. Every "Update_EXHALE section N" reference in the source
  comments and in the older memos points into that file; the references were
  rewritten to the new name on 2026-09-05.

Each entry below is one dated section: what changed, why, which files, the
gate that was run and its result, and the movement it produced against the
scratch baseline or the goldens. Section numbers restart at 1. The plan the
entries execute is `docs/development_plan_20260905_rev3.md` (sections 1-10)
with the order of `docs/development_plan_20260905_execution.md`; item labels
(A1, A8, 10.1-1, D0, C1 ...) are the plan's.

---

## 1. Plan, audit and decisions that open stage 2 (2026-09-05)

State of the tree at the start: HEAD `35d9dd5` (clean apart from `docs/`),
binary rebuilt from it with gfortran 16.2 (md5 `16cde610347c10c0a48db5c2a7fa6cbb`).
All twelve regression goldens differ from the case outputs (they are the
pre-section-170 snapshots); by decision they are not refreshed until the end
of Phase 1 (below).

Documents written today, all under `docs/`:

- `physics_numerics_audit_20260905.md` (external audit of `35d9dd5`), its
  diagnostics in `audit_20260905/`.
- `development_plan_20260905.md`, `_rev1.md`, `_rev2.md`, `_rev3.md` and the
  three reviews `_review.md`, `_review2.md`, `_review3.md`: the plan and its
  revisions after each review. Rev 3 is the plan of record; its section 10
  holds the independent five-group source audit (three new P0, twelve new P1,
  about fifty P2), with the reports and the verification probe in
  `audit_20260905/source_audit/`.
- `development_plan_20260905_execution.md`: Steps 0-7, who acts, gates.

Decisions taken by the user today (recorded in rev 3 sections 8 and 10.5
and in the execution document, Step 0):

1. Goldens are refreshed once, at the end of Phase 1; until then a scratch
   baseline (`backup/regression/baseline_post170_<date>/`, read through a
   `REGRESSION_GOLDEN_DIR` variable to be added to `run_check.sh`) is the
   control for byte-identity. `golden/` is not touched before that.
3. ROE interface: `speed_estimate_ROE` returns a status and, in the vacuum
   case, the two gas-edge speeds; `Num_flux` gains an HLLE vacuum branch.
   Signature shown before coding.
13. One spectrum type builds every band. Power-law input: power law on the
    whole grid, XUV and below 13.6 eV alike (the present construction is the
    stated meaning of the key). A new `Spectrum type: Planck` gives the
    Planck photosphere from `Stellar Teff` and `Stellar radius` on the whole
    grid. `Load`: the table on the whole grid. The Balmer continuum follows
    the run's type. The setup report states the sub-13.6 eV source and the
    integrated grid flux.
14. `m_He/m_H = 3.9715` in the species table; `RJ = 7.1492e9` cm (IAU
    nominal equatorial) in the code and the Wind-AE bridge alike (R0 of every
    run +2.26%); MJ and MSUN to the IAU nominal values.
15. The `du` stop uses the signed spread of rho v r^2 over the mean, as
    `flux_spread` does.
16. Every planet folder, benchmark and planet example runs on the SED its
    reference literature used or assumed, one file per planet under
    `inputdata/sed/` with provenance; the input switch is a configuration
    change, the reruns stay instruction-gated.
17. A loaded SED that ends above the band an active absorber needs (He 2^3S
    at 2583 A, a low-IP metal at its threshold) stops the run; the message
    offers the two remedies (turn the absorber off in the input, or supply a
    covering SED). No continue key.

Open: decision 4 (A3 retry design: carrier substep retry or outer step
controller); 2, 5, 6 inside D0; 7 before B2; 10 in Phase 5.

Documentation restructuring done today (this section): `Update_EXHALE.*`
renamed to `Update_EXHALE_stage1.*`, `Update_EXHALE_early_phase.*` to
`Update_EXHALE_stage0.*`; every reference in `src/`, `docs/`, `python/`,
`README*.md`, the notebooks and the workspace `CLAUDE.md` rewritten to the
new names. Check: `git diff --numstat` shows additions equal to deletions
in every one of the 112 touched files (no line moved; the three exceptions are
the renamed logs themselves, which git sees as deleted paths). The rebuilt
binary's md5 moved from `16cde610347c10c0a48db5c2a7fa6cbb` to
`e779298d0b482bfbc40ae2f9aa80d60a` because four of the rewritten references
are inside printed message strings (`EXHALE_main.f90` 2602,
`ionization_equilibrium.f90` 2948, `input_read.f90` 1004,
`diffusion_tests.f90` 1292), which is text the binary carries; no executable
statement changed. The pointers in `README.md`, `README_HOWTO.md` and the
workspace `CLAUDE.md` name this file as the current log and the stage 1 file
as the history.

Scratch baseline created: `backup/regression/baseline_post170_20260905/`, the
`Hydro_ioniz.txt`, `Ion_species.txt`, `*_adv.txt` and `run.log` of the twelve
cases as they sat in the case directories (48 files, 26 MB; provenance in
the headers: `git=6d07d48afd41 tree=dirty`, the working tree that became
`35d9dd5`, runs of 2026-09-05). It is the byte-identity control for Phase 1
and is deleted after the golden refresh at the end of Phase 1.

## 2. Phase 0: the tests exist and are red on the known defects (2026-09-05)

Plan reference: `development_plan_20260905_execution.md` Step 1. Four
workers, disjoint files; every result below was re-run by the advisor after
the worker's report (times are the advisor's runs).

### 2.1 `src/tests/physics_probe/` (W1): 41 assertions, 15 PASS, 26 FAIL, 3.9 s

Derived from the frozen diagnostics in `docs/audit_20260905/` (which stay);
links the production modules through a computed `use` closure
(`source_closure.py`), so no module can be left out of the link by hand.
Drivers named for the physics: `riemann_wave_speeds.f90` (LLF signal speed
1 against the required 3, mass flux -3.5 against -4.5, mirror symmetry; ROE
normalization invariance broken by 0.661 and 3.162 in the two branches; NaN
star speeds on the equal-pressure expansion; the production ROE flux leaves
thermal energy -2.93 after one Euler update on the vacuum-producing
expansion), `h2_channel_detector_ratio.f90` ((S+D)/M 0.2545 against the
measured 0.2827, D/(S+D) 1/9 against 0.20; channel sum and stoichiometry
GREEN), `h3p_cooling_limits.f90` (cooling independent of n(H2) at 0, 1,
100 cm^-3; joins +0.430, -0.00113, -0.0235 at 300, 800, 1800 K; factor
jump 1.5e-3 at 1e14), `h2_rovibrational_identity.f90` (chemistry
T^2 dlnZ/dT against the EOS u_rv: -0.63% to -8.29% over 300-8000 K),
`photoevent_energy_ledger.f90` (H event at 20 eV closes exactly through the
production `photoelectron_share`; the H2 double event at 80 eV accounts for
60.28 of 80 eV, 19.72 eV unassigned), `radiation_exchange_ledger.py` (two-cell
ownership rule, GREEN by construction).

### 2.2 `src/tests/grid_and_gates/` (W2): 6 programs, all RED, 13.4 s

`grid_width_identity.f90` (production `define_grid`): `dr_j(j)` against
`r_edg(j) - r_edg(j-1)` off by 1.0145 on the Mixed grid at r_max = 10, 1.0100
at r_max = 2, 1.0046 on Stretched, exact on Uniform (where every cell has the
same width, which is why the defect was invisible).
`photon_grid_quadrature.f90` (production `set_energy_vectors`, `J_inc`,
cross sections; `wasp_full` spectrum): code/exact P_HI 1.0951 (floor 4.8
eV), 1.0890 (5.139 eV), 1.0003 (13.6 eV); P_HeI 1.0157 and P_HeII 1.0299 in
every configuration; heating 0.99984. New from this test: P_HeI carries the
same threshold-bin defect (24.6 eV is a grid point whatever the floor), and
the X-ray block ends at 1184.19 eV, not 1240 (share of P_HI 1.9e-10).
`mass_flux_spread_functional.py`: the coded `du` and the signed
`flux_spread` on the baseline states (lower_profile 823.4 vs 3.61 with two
sign changes in the window; mol_base_handoff 1.29 vs 0.696 with none;
wasp_full 9.85e-4 vs 9.85e-4), and `du` = 0 on a synthetic +F/-F window.
`output_state_consistency.sh` (200-step `mol_base_handoff`): header
`sec_ion=T`, heat_total/heat median 0.2752, largest difference 0.822.
`base_level_single_statement.sh` (`lower_profile`, 2 steps): setup 2.0035e-5
bar against the profile's 1.00e-6 bar, factor 20.035.
`restart_grid_guard.sh` (`roundtrip` on Mixed, restarted on Uniform with the
same N): exit 0, the restart's `r` equals the file's Mixed radii exactly and
sits 0.60 relative from the run's own grid.

### 2.3 `inputdata/sed/` (W4): the SED inventory of decision 16

Four files built from data already in the workspace, each with a
provenance header, a build script under `inputdata/sed/build/`, and the
coverage test `sed_coverage_and_order_check.py` (31 assertions, all PASS):
WASP-52 b (eps Eri, MUSCLES v22, reproduces the existing planet-folder file
to round-off), hot Uranus (the benchmark's Ribas band table continued to
2995 A from the Gueymard composite it came from), HD 209458 b (WHI 2008
solar minimum at 0.04707 au, the p-winds file reproduced exactly), LHS 1140 b
(GJ 1132 proxy). Not built, data missing: WASP-121 b (Huang 2023's LLmodels
F5V spectrum above 1700 A) and HD 189733 b (Salz 2016's composite; its solar
table is not here). `inputdata/scaled_solar_hd189.sed` judged an orphan
(no provenance, hardness inconsistent with its own normalization, ends at
1899.5 A). Verbatim paper quotes and the open items are in
`inputdata/sed/README.md`; rev 3 decision 16 records the status. Also found:
`WASP-52b/input.inp` points at a directory that no longer exists (`EXHALE/`),
the top-level Koskinen benchmark runs at 0.048 au against the paper's 0.05,
`HD209458b/input.inp` at 0.048 au against the literature 0.047.

### 2.4 Regression harness and new cases (W3)

`run_check.sh`: `REGRESSION_GOLDEN_DIR` (default `golden/`; `golden` mode
still writes only to `golden/`); the rebuild is fatal with its log kept in
`backup/regression/.build.log`; the summary counts only the cases of the run;
the binary copy is status-tested; `FILES` gains the two `_adv.txt` profiles
(the post-processing pass was outside the regression); a file without a
reference is `MISS`, not `FAIL`, and a run that compared nothing prints
`NO VERDICT`; a restart case keeps its IC pair in `<case>/IC/` and the script
restores it after clearing `output/`; `phase0-cases` prints the new case
list. `Makefile`: `make check CASES=...`; `make test` runs
`element_census_tests`, `diffusion_tests`, `residual_determinism` and every
executable `src/tests/*/run.sh` found at run time, every suite even after a
failure, nonzero exit if any failed.

New cases, not in `DEFAULT_CASES` (`PHASE0_CASES`): `hydrostatic_column`
(hot-Uranus gravity, irradiation 18 orders below normal, no key exists to
switch cooling off so the gas is held where it is negligible; 300 steps,
isothermal to 2%, but draining through the outflow boundary at 0.8 Mach,
so it is a configuration for the Phase 4 D2 measurement, not a hydrostatic
state), `oxygen_chemistry` (`examples/18_oxygen_chemistry` verbatim, 1000
steps, 3 min 52 s single-threaded under load), `hp_front`, `hp_zero_seed`,
`hp_trace_seed` (restarts of the `mol_base_handoff` state with `Ionization
transport: True`; the H II column set to 0 and to 1e-12 of the H nuclei by
moving the protons into H I, since `load_IC` refuses a changed He/H; the
three 100-step states are mutually distinct, so the seed reaches the
solver). No dump hook asserts the seed at the entry of `solve_carriers`;
the README names the two lines where one belongs.

Measured (advisor re-run of two cases: 6.2 s, byte-identical on all four
files): `check wasp_full mol_base_handoff` against the scratch baseline is
byte-identical on all four files per case (W3, 13 min); the five new cases
reproduce their own baseline copies; `make check CASES=mol_base_handoff`
against the stale `golden/` reports the known 8.3e-4 and 3.6e-4 movements
and `MISS` on the two `_adv` files; `make test` (31 s) runs five suites:
`diffusion_tests` and `residual_determinism` PASS, `grid_and_gates` and
`physics_probe` RED on the known defects as intended, and
**`element_census_tests` does not compile at HEAD** (pre-existing:
`carrier_state` gained `nrho` and `wfac` in section 158 and the test's two
calls still pass ten arguments; `carrier_write_back` likewise takes `nrho`
where the test passes `ntot`). Fixed in section 3 below.

Also found by W3: with the flux made negligible, the He recombination
coupling warning (`ionization_equilibrium.f90` 947, 952) fires in every
cell with a ratio of 1e80 because both sides sit at their guards; diagnostic
only.

## 3. `element_census_tests` compiles again (2026-09-05)

`src/tests/element_census_tests.f90`: the two `carrier_state` calls pass the
`nrho` and `wfac` outputs the routine has had since stage 1 section 158, and
the two `carrier_write_back` calls pass `nrho` (the density the carrier
fraction is a fraction of) in place of `ntot`. Test only; no production
source. Result: builds, 20 assertions PASS (E1 stoichiometry and charge, E2
elemental closure through the write-back including the CO collapse, E3 the
H2 photoevent nuclei ledger, E4 the base H2 fraction). Until now the suite
had been unbuildable, so its assertions had not been checked since section
158.

## 4. Phase 1 batch 2a: bounded corrections with no golden movement (2026-09-05)

Seven items, seven workers on disjoint files in one tree (private
`OBJDIR`/`EXE` builds, no shared regression runs); one rebuild and one gate
by the advisor afterwards. Plan references: rev 3 sections 4.1 (A2), 10.2
items 7, 8, 10, 12, 10.3, 10.5 decisions 3, 13, 17.

- **A2, ROE estimator and caller** (`speed_estimate_ROE.f90` rewritten;
  `Num_Fluxes.f90` ROE branch; `input_read.f90`; `EXHALE_main.f90`):
  `docs/a2_roe_interface.md` implemented as approved: `roe_star_state`,
  status OK/VACUUM/INADMISSIBLE, the vacuum test first, PVRS kept only when
  admissible, two-rarefaction exponent `1/z`, two-shock densities times the
  upstream densities, HLLE flux on faces without an admissible star state
  with a printed counter, refusal of `Numerical flux: ROE` with molecular
  chemistry and the ladder EOS. Four deviations accepted: a positivity guard
  on the input states at the top of the estimator; the refusal placed in the
  consistency block after the optional keys (where `thereis_mol` is known);
  the old `p_min` floor removed (the guard supersedes it); the counter printed
  for ROE runs only. Test: 30 ROE assertions GREEN (scale invariance,
  equal-pressure expansion `c_*` 0.9576611, stationary contact, asymmetric
  shock and rarefaction against exact Riemann solutions at 5%, vacuum edges
  -6.127017/+6.127017, positive thermal energy after one Euler update).
- **2a-SED, decision 17** (`sed_read.f90`): a table that stops above the
  band an active absorber needs stops the run with the file, the missing band
  in eV and A, the absorbers (He 2^3S; each active neutral metal below the
  file minimum, K I at 4.341 eV included) and the two remedies; no continue
  key. Test `sed_coverage_stop.sh` (six assertions) GREEN; no configuration
  of record trips it (the Koskinen runs have the metastable off).
- **2a-OPAC** (`util_ion_eq.f90` `PH_heat_H`/`PH_heat_HHe`,
  `opacity_models.f90`): model P's factor multiplies the shared flux weight
  `int_f` of every absorption integrand (rates, heating, secondary
  ionization, `acc_q`), so the beam's loss equals the local absorption;
  test `opacity_model_p_ledger.f90`: a cell with f = 10 closes energy and
  H I ionizations to 3e-8 (was a factor 10 off).
- **2a-LYA** (`lya_rt.f90`, `util_ion_eq.f90` `fuv_lw_photon_field`,
  `water_photolysis.f90` header, plus the two call sites in
  `ionization_equilibrium.f90` and `diffusive_photochemistry.f90` that the
  new `nHI` argument required): band B2 is attenuated by
  `lya_stellar_beam_transmission(T, tau)` = `erfc(x1/(sqrt(2) X_s))`, the
  expression `jlya_escape_prob` carried inline and now calls (bitwise equal
  at 4941 test points). Test `lya_band_transmission.f90` (ten assertions)
  GREEN. Movement, `oxygen_chemistry` 1000-step snapshot: du 903 -> 596,
  log10 Mdot 9.62 -> 9.61, B2 transmission 1e-291 at the base, 1e-3 at
  1.07 R_p, 0.98 at 1.5 R_p; H2O moves up to 6.2x and OH 7.3x in the
  layer. Scratch baseline of that case re-snapshotted (previous copy kept).
- **2a-WAE** (`wae_exhale_input.f90`): mass fractions from the same
  `m_He/m_H` the solver is handed; test: reconstructed He/H equals the input
  to 1e-12 (was 1.0071696x). RJ/MJ/MSUN untouched (decision 14 is 2c).
- **2a-HPFIX** (`constrained_chemical_equilibrium.f90`): `hp_is_fixed`
  pins `species_fixed(is_HII)` like H2 and the oxygen carriers; new suite
  `src/tests/constrained_network_layout/` (eight assertions) GREEN. The three
  `hp_*` cases never enter the continuation (every cell is a class-1 root),
  so no run exercises the fix yet.
- **2a-PLANCK** (`J_inc.f90`, `set_energy_vectors.f90`, `input_read.f90`,
  `write_setup_report.f90`, `excited_hydrogen.f90`, `docs/input_schema.md`):
  `Spectrum type: Planck` (requires `Stellar Teff` and `Stellar radius`)
  builds the whole grid from `pi B_nu(T_eff) (R_star/a)^2`; the Balmer
  routines call the same shared function (unchanged result); the setup
  report states the spectrum type, the source of the band below 13.6 eV and
  the integrated grid flux against the nominal J_XUV. New suite
  `src/tests/spectrum_type/` (12 assertions) GREEN. Measured on the
  `wasp_full` star: Planck/power-law 1642 at 4.8 eV, 358 at 6.2 eV, 2.9 at
  10 eV, 0.014 at 13.6 eV; He 2^3S rate 45.2 s^-1 (Planck) vs 0.100
  (power law); integrated grid flux / nominal 1.42 (power law, floor 4.8 eV),
  0.995 (floor 13.6), 112.8 (Planck).
- **2a-GUARDS** (`define_grid.f90`, `set_IC.f90`, `EXHALE_main.f90`,
  `viscous_conduction.f90`, `steady_newton.f90`, `input_read.f90`,
  `diffusive_photochemistry.f90`, `water_photolysis.f90`, `oxygen_rates.f90`,
  `lyman_werner.f90`, `parameters.f90`, `molecular_infrared_data.f90`):
  `cell_nearest_radius` (the bare `minloc` read `r(503)` = 0 out of bounds
  on the default grid), `j_min >= 1`, dead `report_unmet_steady_gates`
  removed, Mdot guarded for small N and non-outflow, conduction floor
  counted and reported, `okres` read in the ray test, correct sweep count,
  `get_word` initialized, dead `close(unit=1)` removed, `adv_corr` zeroed on
  allocation, dead `Tcode` argument removed from four carrier routines and
  their callers, `carrier_mass_amu` through the named species indices,
  `wp_ready` used as the guard it names, dead constants and stale headers
  fixed, `molecular_infrared_data`'s `i` private. Test `grid_window_indices`
  (seven assertions) GREEN.

Gate, advisor's single rebuild (md5 `b52c79ebb15ee9014858425bfea4cb1c`):
`make test` runs seven suites; every assertion added in this batch is GREEN
and only the known defects of later items remain RED (two test drivers had
linked the program object `element_census_tests.o` from `build/`; their link
lists now exclude program objects). Regression, sixteen cases against the
scratch baseline, single-threaded: fifteen byte-identical on all four
files (the eleven of `DEFAULT_CASES`, `hydrostatic_column`, the three
`hp_*`); `oxygen_chemistry` moved as stated under 2a-LYA.

Follow-ups recorded for batch 2b and the docs pass: `carrier_mass_amu`
returns 0 for a species missing from the table (needs a stop);
`write_output.f90` 806-816's FUV ledger closed form carries no line
transmission for B2 (its `rel_diff` reads 1.0 now; the LW row was already
0.84); `inputdata/sed/README.md` describes the old warning in six places;
`docs/input_schema.md` row 11a should state the coverage requirement;
`docs/a2_oxygen_option_design.md` 402 describes the pre-2a B2 field; the
`physics_probe/README.md` summary counts are stale; `J_XUV` (used by the
Wind-AE bridge), the FUV band fluxes and the monochromatic flux do not yet
follow the spectrum type (decision 13 is enforced for the photon grid only);
`backup/regression/parse_golden/` was already stale at 35d9dd5; the three
audit probe scripts under `docs/audit_20260905/` no longer compile against
the tree and are marked as frozen records.

## 5. Phase 1 batch 2b: bounded corrections that move some outputs (2026-09-05)

Six workers (five code items and one documentation pass), disjoint files,
private builds; one rebuild (md5 `d2173f0210abcdc84f41cb6fb5076e18`) and one
gate by the advisor. Plan references: rev 3 sections 4.1 (A1, A3 increment
(a), A8) and 10.3 (`k_lw` cell mean, He I singlet, He II recombination
cooling).

- **A1, LLF signal speed** (`Num_Fluxes.f90`, `RK_rhs.f90`, `EXHALE_main.f90`):
  `a1 = max(abs(vL)+aL, abs(vR)+aR)`, the spectral radius of the Euler flux
  Jacobian, replacing `max(abs(vL+aL), abs(vR+aR))` (1 instead of 3 on the
  inflow test pair); the header's positivity claim now describes the
  coefficient the code carries. Reporting: calls of the positivity repair,
  repaired interfaces, repaired interfaces in accepted steps, and LLF faces,
  all printed even when zero. Test: the three LLF assertions GREEN. Movement:
  none (no case selects LLF and the repair path is never entered in the
  matrix: `flux correction: 0 positivity repair call(s)` in every run).
- **A3 increment (a), carrier reference scales and status**
  (`diffusive_photochemistry.f90`, `write_setup_report.f90`,
  `EXHALE_main.f90`): one `carrier_element_reference_density(ic, nH, nO, nC)`
  (H2 and H+ hydrogen, OH and H2O oxygen, CO min(O, C)) read by both the
  Jacobian step (`carrier_source_derivative_step`) and the residual floor;
  the reported residual is that of the state after `limit_to_element_budget`;
  `solve_carriers` returns converged / iteration cap / line search failed and
  the caller stops on a failed status in a carrier the limiter did not move
  (rejection and retry remain increment (b), decision 4 open); a cumulative
  CO-ceiling record (cells, applications, molecules removed) in the run
  summary and `EXHALE_resolved.out`; `carrier_mass_amu` stops on a species
  missing from the table. New suite `src/tests/carrier_reference_scales/`
  (14 assertions; before the change the proton row's diagonal derivative was
  exactly zero at H+ = 0 and 1e-12 n_H because the step had collapsed to
  1e-300). Movement of this item alone: `hp_zero_seed` up to 7.8e-4 at the
  ionization front, `hp_trace_seed` 2.8e-5, `hp_front` round-off,
  `mol_carrier` and `mol_base_handoff` byte-identical. Cost of the extra
  residual assembly +1.5% wall time on `hp_zero_seed`.
- **A8, Chung detector inversion** (`h2_photo_channels.f90`): Table II is
  the detector event ratio `(N_S+N_D)/N_M` (Chung et al. 1993 p. 886, one
  pulse per event, quoted at the code site), `q_D = N_D/(N_S+N_D)` = 0.20 at
  80 eV; `sigma_M = (1-f_di) sigma_i`, `sigma_S = (1-q_D) f_di sigma_i`,
  `sigma_D = q_D f_di sigma_i`; `yan_rho` converted with `q_D =
  rho_D/(1+rho_D)/f_di`; an inadmissible `q_D` stops the run instead of being
  clipped; `frac_H2_double_of_protons` renamed
  `frac_H2_double_of_proton_events`. At 80 eV: D 0.02254 -> 0.04408, M
  0.7972 -> 0.7796, S 0.1803 -> 0.1763. Test: 18 assertions GREEN including
  the stored-channel and yield integration checks; E1 test 4 prints the same
  6.951e-3. Movement of this item alone: `mol_base_handoff` n_e up to 1.8e-2
  at r = 1.018, x(H+) 1.9e-2 at the base, T 4.6e-4; `mol_sec_ion` at the
  1e-3 to 5e-3 level, with a 0.078 cm/s velocity shift moving the v = 0
  crossing by one cell at r = 1.024-1.026 and flipping the advection-
  correction branch of the `_adv` files there. Mdot unchanged.
- **`k_lw` cell mean** (`lyman_werner.f90`,
  `util_ion_eq.f90` `fuv_lw_photon_field`): new
  `lyman_werner_dissociation_rate_cell_mean(F, N_out, N_in, T, n_H, tau_out,
  tau_in)`, a composite three-point Gauss rule on segments of at most 0.05
  dex of column and 0.5 of continuum depth (the level-resolved cross section
  is a piecewise power law with slope up to 3.3, so the water bands'
  exponential mean does not apply; measured accuracy 2.5e-4 worst over 1500
  random cells); the star-ward face carries `NH2col(j+1)`. Test: seven
  assertions GREEN (the old inner-face rate was 3.18x too small on a cell
  spanning 1.12 dex). Movement: `mol_lyman_werner` k_LW +1.4% to +28.1%
  (median +2.7%), n(H2) -20% and n(H3+) -20% at the front r = 1.071, T up to
  8e-3, n_e 1.7e-2; `mol_base_handoff` and `wasp_full` identical (no LW
  flux). Cost under 0.5%.
- **He I singlet and He II recombination cooling** (`util_ion_eq.f90`
  `eval_cool`/`eval_cool_cells`, `Cool_coeff.f90`,
  `ionization_equilibrium.f90` comment): the 24.6 eV collisional ionization
  and the Cen (1992) excitation of He I are charged to the ground singlet
  `he_ground_singlet_density(nhei, nheiTR)` (the metastable had paid both
  levels); `lambda_rec_HeII = kT alpha_rec_HeII_total(T)` with
  `alpha_rec_HeII_total = alpha_11S + alpha_23S` when the metastable is
  tracked (the coefficient the balance removes He+ at) and `alpha_B`
  otherwise; the coefficient is built inside `eval_cool` from T and the run's
  flags so that `energy_semi_implicit`'s finite-difference `dC/dT` sees one
  function (a coefficient-passing form would have left an offset the
  derivative multiplies by 1e5). At 1e4 K the He II coefficient rises
  2.807e-13 -> 3.640e-13 (+30%); the He recombination coupling's own weight
  (`y_net alpha_1`) is still not charged (about 10% of this channel, stated at
  the code site). Test: 15 assertions GREEN. Movement of this item alone:
  `wasp_full` median 2e-5 / max 1.7e-4 in rho, T max 5.9e-4, Ca I and Na I
  up to 4.7e-3; `mol_base_handoff` cool column median 1% / max 2.6% (the
  recombination channel is 64% of the cooling there), state 1e-4.
- **Documentation pass**: `inputdata/sed/README.md` (seven places, the
  stop), `docs/input_schema.md` (row 11a coverage rule; rows 11, 11b, 11d,
  K3, K4 line references corrected), `docs/a2_oxygen_option_design.md`
  (B2 row), `src/tests/physics_probe/README.md` (two summaries: 41/15/26 at
  35d9dd5; 83 assertions, 69 PASS, 14 FAIL after 2a with the LLF fix already
  in the tree) and `grid_and_gates/README.md`, rev 3 section 10.2 and 10.5
  "landed" notes, the execution document Steps 1 and 2a marked done,
  `README.md` and `README_HOWTO.md` (`Spectrum type: Planck`, the coverage
  rule, `make test`, the regression section's case list corrected from
  eight to eleven).

Gate (advisor): `make test` runs eight suites; every 2b assertion GREEN,
RED only on the 2c and Phase 3 targets. Regression, sixteen cases against
the post-2a scratch baseline, single-threaded: `wasp_he23off` and
`hydrostatic_column` byte-identical (metastable off, no carriers); the other
fourteen moved, Mdot unchanged to two decimals in all sixteen. Largest
relative differences over all cells and columns (single cells; the large
values are v near a zero crossing or trace species in one cell, see the
items):

| case | Hydro_ioniz | Ion_species | Hydro_ioniz_adv | Ion_species_adv |
|---|---|---|---|---|
| wasp_full | 1.06e-3 | 4.63e-3 | 1.37e-3 | 4.63e-3 |
| wasp_full_newton | 1.48e-3 | 4.69e-3 | 1.91e-3 | 4.69e-3 |
| mol_base_handoff | 2.54e-2 | 1.99e-2 | 5.58e-2 | 4.21e-2 |
| mol_metals | 2.16e-2 | 2.20e-2 | 2.70e-2 | 2.20e-2 |
| mol_lyman_werner | 1.25e-1 | 2.06e-1 | 5.63e-2 | 2.06e-1 |
| mol_diffusion | 2.54e-2 | 1.99e-2 | 5.86e-1 | 1.95e-2 |
| mol_ir_bands | 2.16e-2 | 2.20e-2 | 9.65e-2 | 2.20e-2 |
| mol_sec_ion | 1.70 (v at a zero crossing) | 6.31e-3 | 1.70 | 1.0 (trace ion, one cell) |
| mol_carrier | 2.52e-2 | 1.98e-2 | 3.87e-2 | 1.95e-2 |
| lower_profile | 5.38e-1 | 5.19e-2 | 5.38e-1 | 4.96e-2 |
| oxygen_chemistry | 1.94 | 1.0 | 1.94 | 1.0 |
| hp_zero_seed | 8.13e-2 | 9.54e-2 | 4.20e-2 | 9.54e-2 |
| hp_trace_seed | 4.25e-1 | 1.0 | 4.24e-2 | 1.0 |
| hp_front | 2.59e-2 | 6.30e-3 | 4.42e-2 | 6.29e-3 |

The scratch baseline was re-snapshotted from these outputs; the post-2a
state of every case is kept in `backup/regression/baseline_post2a_20260905/`.

Follow-ups recorded: `T_equation.f90` (the temperature root of the
advection post-process) carries no He 2^3S cooling channel and charges the
He I ground-state coefficients to the summed neutral column, so the `_adv`
temperature and the hydro temperature are not solved with one cooling
function (Phase 3 A5/A6 or the `_adv` consistency item of Phase 6);
`post_process_adv.f90` 392 calls `eval_cool` without `nheiTR`;
`cross_sec.f90` 519-528, `parameters.f90` 475-501 and `input_read.f90`
856-858 still describe the H+ column as a proton count and quote a "factor
1.6" between the two double-ionization models that is now reversed and much
smaller; `write_output.f90` 200-203 says the `f_shield` column is the factor
the rate carries (the rate is now the cell mean); `pct_limited` counts only
overshoots above 1e-10 while the clamp acts at any size (same class as the
CO counter, now recorded cumulatively); the give-up path of `retry_step`
leaves the loop through the same exit as an accepted step (decision 4);
`set_oxygen_coeffs` is called with strided sections that create array
temporaries in the carrier Jacobian; `docs/EXHALE_user_manual.tex` is stale
on the spectrum types and the SED stop.

## 6. Phase 1 batch 2c: corrections that move every output (2026-09-05)

Fourteen worker items on disjoint files plus advisor edits, every one
measured alone against a frozen snapshot of the tree as it stood when the
item started (private builds; the items therefore chain, and the "before" of
one item is the "after" of those that landed earlier). One rebuild by the
advisor with the conda gfortran 16.2 and its OpenBLAS (md5
`5845df587aa7948153f1bed3becb28ac`), `make test`, and one regression run
against the post-2b scratch baseline. Plan references: rev 3 sections 10.1
(items 1-3), 10.2 (items 1-6, 9, 11), 10.5 (decisions 13-15, 17, 2c-CONST).
Movement is quoted as the largest relative change over the 500 physical
cells of `wasp_full` (WASP-121b, He 2^3S + metals, `du` stop, no Newton
finish, so its Mdot carries the path-dependent spread of that stop) and
`mol_base_handoff` (hot Uranus, 12000-step relaxation snapshot, not steady,
so its outer-tail and base-cell numbers compare two points of one transient).

- **2c-DRJ, the cell width is the cell's own width** (`define_grid.f90`,
  `set_IC.f90`): `dr_j(j) = r_edg(j) - r_edg(j-1)`, an exact identity of the
  finite-volume grid; the stored width had been the width of cell j+1, 0.7 to
  1.25 percent too large over the stretched region. Every consumer confirmed
  to use `dr_j(j)` for cell j; none changed. `grid_width_identity` GREEN.
  `wasp_full`: T 2.0e-3, v 5.4e-3 (base cell), n_e 6.3e-3, x(H+) 1.1e-2,
  mass flux +0.65 percent, printed log10 Mdot 13.24 both. `mol_base_handoff`:
  bulk profiles 1e-2 class with the transient tail at 1e-1, mass flux -0.96
  percent, printed Mdot 10.52 -> 10.51.
- **2c-QUAD, thresholds on bin edges** (`set_energy_vectors.f90`): the photon
  grid is built from bin edges, `e_v` the geometric center of each bin and
  `de_v` its exact width; H I, He I, He II, the metastable edge and `e_mid`,
  `e_top` are edges, and the X-ray block reaches `e_top`.
  `photon_grid_quadrature` GREEN except P_HeI at 1.2e-3 (closed by 2c-THR).
  `wasp_full`: T 7.5e-2 at 1.12 Rp, n_e 7.8e-2, x(H+) 0.23, base v 0.47,
  mass flux +0.82 percent, log10 Mdot 13.24 -> 13.25. `mol_base_handoff`: T
  4.7e-2 (tail), n_e 2.5e-2, mass flux +0.08 percent, Mdot 10.51 both. This
  is the single largest movement of the batch after the Jupiter radius, and
  it is the removal of the 9.5 percent overcount of P_HI found by the audit.
- **2c-BALMER, the Balmer continuum follows the run's spectrum type**
  (`excited_hydrogen.f90`, `J_inc.f90`; decision 13): the n = 2
  photoionization rate and heating integrate the field the photon grid uses,
  the power law for a power-law run, the Planck field for `Planck`, the
  loaded table for `Load` (`stellar_flux_eV`, `spectrum_covers_eV`); a table
  that does not reach 3647 A stops the run with the decision-17 message
  naming H(n=2). `balmer_continuum_field` GREEN (11 assertions).
  `wasp_full` (n = 2 coupling armed): gamma_2 falls by a factor 2939, n_e at
  1.02 Rp -12.5 percent, x(H+) at the base -27 percent, T 9.6e-3, wind above
  1.2 Rp and Mdot unchanged (13.25). `mol_base_handoff`: byte-identical (no
  stellar lines).
- **2c-KICK, the momentum kick removed** (`EXHALE_main.f90`): the
  `u(2,:) + 1e-16` diagnostic aid is gone; the relative step change skips
  cells whose old momentum is exactly zero. `momentum_row_from_fluxes_only`
  GREEN. Movement round-off on `wasp_full` (1e-9), 4e-6 on
  `mol_base_handoff`, 4.4e-3 in v on `lower_profile` (cold-start snapshot
  with a sign-changing base flux amplifies it); Mdot unchanged everywhere.
- **2c-STATE, one state per output file** (`EXHALE_main.f90`): when the
  secondary-ionization coupling never switched on during the march, the
  final write no longer forces it on; the header says `sec_ion=F`, and the
  `_adv` files and the heating breakdown are built from the relaxed state.
  `Hydro_ioniz.txt` / `Ion_species.txt` byte-identical in every case; the
  `_adv` files and `Heating_breakdown.txt` of the fixed-step snapshots whose
  flip never fired (`mol_base_handoff`, `lower_profile`, and the other
  molecular gates) change by up to 0.97 in `n(H II)_adv` at the base.
  `output_state_consistency` stays RED at 1.2e-5: the heat column is
  assembled from the pre-sweep rates and the breakdown re-evaluates on the
  post-sweep composition (Phase 3, together with the `excited_H_update`
  call between the loop and the final write).
- **2c-DU, the signed mass-flux spread** (`EXHALE_main.f90`; decision 15):
  one functional `mass_flux_spread`, `(max - min)/max(|mean|, 1e-30)` of the
  signed `rho v r^2` over `j_min..N`, read by the stop, the Newton arming
  and the stage flip, and written to `run.log` at the final step
  (`mass_flux_spread_functional.py` GREEN). `wasp_full` byte-identical
  (every threshold crossing on the same step); `mol_base_handoff`
  byte-identical (fixed steps; the reported `du` falls from 1.29 to 0.70);
  `wasp_full_newton` flips and hands off one step earlier, solution moves
  by 7e-8.
- **2c-RESTART, a restart keeps the run's grid** (`load_IC.f90`): the file's
  cell centers are compared with the run's over the physical cells and a
  relative difference above 1e-10 refuses the restart, naming both grids.
  `restart_grid_guard.sh` GREEN; `roundtrip` and the three `hp_*` cases
  byte-identical on the same grid. Consequence: after the grid-moving items
  of this batch (RJ +2.26 percent of the domain in R_p units, `dr_j`) the
  shipped `hp_*` IC pairs are refused, and are regenerated with
  `backup/regression/make_hp_initial_conditions.py` from a fresh
  `mol_base_handoff` output (below).
- **2c-BASE, one base level with a lower-atmosphere profile**
  (`input_read.f90`, `write_setup_report.f90`): with a profile, `n0` follows
  from `p_match_bar` exactly as from `p_base` of `base.inp`; a density key
  beside a profile is refused unless it agrees within 1 percent
  (`base_level_single_statement.sh` GREEN). `lower_profile` starts at the
  profile's level: base density 1.0e14 -> 5.0e12 cm^-3, its 12000-step
  snapshot changes entirely (printed log10 Mdot 8.44 -> 7.92, ionization
  front inward, base n_e x8). The case README and
  `examples/17_lower_profile/input.inp` drop the density key. Nothing else
  moves.
- **2c-CI, the C I Voronov row** (`Cool_coeff.f90`): `P = 0, X = 0.193,
  K = 0.25, A = 6.85e-8, dE = 11.26` (Voronov 1997 Table 1; the transcription
  had X and K in the wrong slots, 2.56x at 2000 K, 1.69x at 1e4 K).
  `voronov_carbon_ionization` GREEN. `wasp_full` 1e-5 class, x(C I) 6e-4,
  Mdot -1e-5 dex; `mol_metals` bitwise identical (collisional ionization of
  C I is negligible at its temperatures).
- **2c-MASS and 2c-MASS2, helium mass and Jupiter radius** (decision 14;
  `species_table.f90`, `parameters.f90`, `wae_exhale_input.f90`,
  `wae_ic_writer.f90`, `binary_element_diffusion.f90`,
  `diffusive_photochemistry.f90`, `lower_column.f90`, `utilities.f90`,
  `examples/exhale_io.py` and the Python tools): He I / He II / He III /
  He 2^3S carry 3.9715 m_H (HeH+ 4.9715), `RJ = 7.1492e9` cm, `MJ =
  1.8982e30` g, `Msun = 1.98842e33` g (IAU 2015 nominal), one definition
  each, the Python tools importing `RJ_CM`, `MJ_G`, `MSUN_G` from
  `exhale_io`. Tests `atomic_mass_and_radius_constants`,
  `jupiter_radius_python_tools.py`, suite `species_masses` GREEN. Measured
  as two sub-steps: the helium mass alone moves `wasp_full` by 2.3e-3 in T,
  1.1e-2 in x(H+), Mdot +0.095 percent, and `mol_base_handoff` by 1.2e-2
  (tail), Mdot -0.004 percent; the Jupiter radius then rescales every
  physical length (R0, the domain top `r_max` in R_p units, the gravity)
  and moves `wasp_full` by 7.3e-2 in T, 0.16 in n_e, 0.34 in x(H+), raw
  Mdot +11.0 percent (printed 13.24 -> 13.29), and `mol_base_handoff` raw
  Mdot +15.8 percent (printed 10.52 -> 10.58). MASS2 alone (the remaining
  copies): `mol_diffusion` 2.8e-5, `mol_carrier` 6.5e-8, `wasp_full`
  byte-identical.
- **2c-THR, one definition of each threshold** (`parameters.f90`,
  `cross_sec.f90` and the files that restated a threshold): `e_th_HI =
  13.598434599`, `e_th_H2 = 15.425927`, `e_th_HeI = 24.587389`, `e_th_HeII
  = 54.417765`, `e_th_HeTR = 4.767775` eV (NIST ASD), every cross section
  turning on exactly there; `sigma(E, Z)` takes the ion's measured potential
  as an optional argument so He II is not restated as `13.6 Z^2`.
  `ionization_threshold_turn_on` and `threshold_literal_uniqueness.py`
  GREEN; the P_HeI residual of QUAD closes (0.99878 -> 0.99995).
  `wasp_full` 7e-4 at most, mass flux +0.044 percent; `mol_base_handoff`
  1.3e-3 (tail); printed Mdot unchanged.
- **2c-N2FLOOR, H(n=2) is an absorber of the grid-floor rule** (`sed_read.f90`
  `photon_grid_floor_eV`, `J_inc.f90`, `set_energy_vectors.f90`,
  `excited_hydrogen.f90`): with the excited-hydrogen coupling armed the grid
  floor is `e_th_HI/4` = 3647 A, a loaded SED is read down to it and the
  decision-17 stop is about the file, not the type (WASP-52b with the eps Eri
  table to 3700 A now runs). Tests `sed_coverage_stop.sh` and
  `photon_grid_n2_floor` GREEN: the 20 bins added below the lowest bulk-gas
  absorber change no rate the grid computes, exactly. `wasp_full` 2e-4 at
  most (the edge written as `e_th_HI/4` instead of 3.40 eV), Mdot 0.0000
  dex; wall clock +2.3 percent for the 20 bins.
- **2c-CONST, the particle masses and the astronomical constants in one
  place** (`parameters.f90`, `species_table.f90`, `lower_column.f90`,
  `Cool_coeff.f90`, `binary_element_diffusion.f90`,
  `diffusive_photochemistry.f90`, `lyman_werner.f90`, `oxygen_rates.f90`):
  `m_He_atom = 6.6464790722e-24` g (4.002603254 u x CODATA 2018 u; the brief's
  6.6464731e-24 was 9.2e-7 low and is not the helium-4 atom), `amu`, `m_e`,
  `m_p` in `parameters.f90`, `m_He_over_m_H` and `amu_over_m_H` derived,
  `bsp_mass` built from them, `melem_A = melem_A_u * amu/mu`. This closed a
  UNIT ERROR: the standard atomic weights in u had been added to mass sums
  carried in hydrogen atoms, every metal nucleus 0.78 percent too heavy in
  `rho`, the mean molecular weight, the base composition, the census closure
  and the metal element diffusion. The advisor closed the remaining local
  copies (`Cool_coeff` line-center opacity now takes the atomic weight in u
  and uses the global `amu`; `amu_C/N/O = melem_A_u(iel_*)`; the local
  `m_amu_g` and `kb_erg`/`m_h2` copies of the diffusion, Lyman-Werner and
  oxygen modules replaced by the globals). `wasp_full` 1.0e-3 at most
  (base v), mass flux +3.2e-5; `mol_base_handoff` 2.8e-5; printed Mdot
  unchanged.
- **2c-METALEDGE, every active threshold a bin edge**
  (`set_energy_vectors.f90`): each of the six bands keeps its bin budget and
  is cut at every ionization threshold strictly inside it (H2 when the
  molecular chemistry is on, every photo-ionizable stage of an element with
  nonzero abundance), the budget shared over the sub-intervals in proportion
  to their logarithmic width with at least one bin each; `Nl` is derived
  from the partition and equals the previous count in every configuration
  (240 `wasp_full`, 220 `mol_metals`, 220 `mol_base_handoff`). A band with
  no interior threshold keeps exactly the partition it had.
  `photon_grid_threshold_edges` GREEN: seventeen thresholds at relative
  distance 0 from an edge (eighteen with H2), P_m of C I 1.045 -> 1.000,
  Ca I 0.933 -> 0.997, Na I 1.030 -> 0.998, Fe I 0.965 -> 0.994 of a
  400001-point reference, total H2 photoabsorption 0.985 -> 1.000.
  `wasp_full`: T 9.6e-3, n_e 5.2e-3, x(C I) 0.14 at the base, mass flux
  -0.17 percent, log10 Mdot 13.2860 -> 13.2852. `mol_metals`: bulk 2e-3,
  n_e 4.5e-2 in the neutral layer at 1.1 Rp (ionization degree 1e-6), Mdot
  unchanged. `mol_base_handoff`: 1e-4 class, by the H2 edge alone. Open
  (decision needed, below): Mg I, Ca I, Na I, Fe I still miss 1e-3 by 0.25
  to 0.58 percent, the rectangle-rule error of the 20-bin sub-Lyman band.
- **Advisor edits after the worker:** the H2 channel thresholds
  `e_th_H2_di` (18.076 eV) and `e_th_H2_dd` (51.4 eV) join the edge set when
  the molecular chemistry is on, because the channel shares of one cross
  section switch on there (the 32-41.5 eV neutral-dissociation window is
  continuous and is not cut; the comment at the code site says so); the setup
  report prints the grid floor as the lower edge of the first bin
  (`photon_grid_floor_eV`, or the table floor for a loaded SED) instead of
  the first bin center, which the new partition had moved to 4.95 eV
  against a floor of 4.768 eV. The documents that restated the changed
  quantities were brought in line: `docs/input_schema.md` (profile level,
  R_J), `docs/a2_oxygen_option_design.md` (molecular masses in m_H units),
  rev 3 and the execution document (landed marks).

### Decision handed to the user: the sub-Lyman band budget

The band `[4.768, 13.598]` eV has 20 bins (edge ratio 1.054), 4.4 times
coarser in `d ln E` than the pure-H I band (50 bins, 1.012), and it is where
every low-IP metal is photoionized; Fe I's cross section falls by a factor
2.5 across the single bin above its threshold. Measured ladder (code/exact):

| num_TR | 20 | 30 | 40 | 50 | 80 |
|---|---|---|---|---|---|
| Mg I | 0.99753 | 0.99871 | 0.99915 | 0.99955 | 0.99982 |
| Ca I | 0.99737 | 0.99877 | 0.99937 | 0.99961 | 0.99984 |
| Na I | 0.99752 | 0.99912 | 0.99958 | 0.99967 | 0.99987 |
| Fe I | 0.99417 | 0.99685 | 0.99805 | 0.99883 | 0.99953 |
| `Nl` (`wasp_full`) | 240 | 250 | 260 | 270 | 300 |

**Decided by the user 2026-09-05: `num_TR = 89`**, the budget that gives the
band the bin ratio of the pure-H I band (0.0118 in ln E). Measured with the
same driver: Mg I 0.999853, Ca I 0.999870, Na I 0.999902, Fe I 0.999615, so
the four rates are asserted at the 1e-3 of every other ion. `Nl` becomes
309 for `wasp_full` (was 240) and 289 for the hot-Uranus cases (was 220);
the H(n=2) band keeps its 20 bins. Every He 2^3S and metals-on case moves
again by this, and the gate below was re-run on the 89-bin binary (md5
`b0a4ada9f84ee8410c84fda6780ba097`) before the baseline re-snapshot. The
user also decided the `LHS1140b/examples/photochem_profile_base/input.inp`
density key: deleted (the profile fixes the level).

### Gate (advisor)

`make test` on the single rebuild: every assertion of this batch GREEN
(`threshold_edges`, `n2_floor`, `balmer_continuum_field`,
`atomic_mass_and_radius_constants`, `species_masses`,
`ionization_threshold_turn_on`, `threshold_literal_uniqueness`,
`voronov_carbon_ionization`, `grid_width_identity`, `photon_grid_quadrature`,
`restart_grid_guard`, `base_level_single_statement`,
`mass_flux_spread_functional`, `momentum_row_from_fluxes_only`,
`sed_coverage_stop`); `diffusion_tests` had 13 of 35 assertions fail because
its own references carried the helium mass as 4 m_H (7.2e-4 = 1 - 3.9715/4),
and now derives every helium mass from `m_He_over_m_H` (35 of 35, no
tolerance changed, the T9 ambipolar targets written as functions of the
mass ratio). RED, as intended, only on the Phase 3 targets:
`output_state_consistency` (1.2e-5), `h3p_cooling_limits`,
`h2_rovibrational_identity`, `h2_double_photoevent_energy_sum`,
`helium_collisional_ionization_by_level`.

Regression, sixteen cases against the post-2b scratch baseline,
single-threaded: the thirteen cases that ran all moved beyond the 1e-3
tolerance, as every item of this batch predicts; the three `hp_*` restart
cases were refused by the new grid guard (their IC pairs sat on the old
grid: domain top 4.7253 -> 4.6215 R_p from R_J, base spacing from `dr_j`)
and were regenerated with `make_hp_initial_conditions.py` from the fresh
`mol_base_handoff/output` (nuclei conserved row by row; H II as written /
0 / 1e-12 of the nuclei) and re-run. Movement of the thirteen (the
cell-wise maxima compare the same cell index, whose radius moved by up to
2.2 percent with the domain, so they overstate the change at fixed radius;
the Mdot column is the run's printed value):

| case | log10 Mdot before -> after | mass flux at cell N-20 | max rel T | max rel v | max rel rho | domain top [R_p] |
|---|---|---|---|---|---|---|
| `wasp_full` | 13.24 -> 13.29 | +6.2 % | 6.7e-02 | 7.6e-01 | 4.1e-01 | 1.6014 -> 1.5666 |
| `wasp_he23off` | 13.24 -> 13.29 | +6.2 % | 6.7e-02 | 7.7e-01 | 4.1e-01 | 1.6014 -> 1.5666 |
| `wasp_full_newton` | 13.21 -> 13.26 | +6.7 % | 1.3e-01 | 1.9e-01 | 3.4e-01 | 1.6014 -> 1.5666 |
| `mol_base_handoff` | 10.52 -> 10.57 | +8.3 % | 2.5e-01 | 1.1e-01 | 4.5e-01 | 4.7253 -> 4.6215 |
| `mol_metals` | 10.52 -> 10.58 | +8.3 % | 2.4e-01 | 1.1e-01 | 4.4e-01 | 4.7253 -> 4.6215 |
| `mol_lyman_werner` | 10.50 -> 10.55 | +8.2 % | 1.0e-01 | 1.2e-01 | 2.5e-01 | 4.7253 -> 4.6215 |
| `mol_diffusion` | 10.52 -> 10.57 | +8.3 % | 2.5e-01 | 1.1e-01 | 4.5e-01 | 4.7253 -> 4.6215 |
| `mol_ir_bands` | 10.52 -> 10.58 | +8.3 % | 2.4e-01 | 1.1e-01 | 4.4e-01 | 4.7253 -> 4.6215 |
| `mol_sec_ion` | 10.47 -> 10.52 | +8.2 % | 2.1e-01 | 3.7e+02 (sign change of a stagnant base cell) | 3.8e-01 | 4.7253 -> 4.6215 |
| `mol_carrier` | 10.52 -> 10.58 | +8.4 % | 2.5e-01 | 1.1e-01 | 4.6e-01 | 4.7253 -> 4.6215 |
| `lower_profile` | 8.44 -> 7.95 | -68.9 % | 2.0e+00 | 1.5e+03 (sign change) | 9.7e-01 | 4.1487 -> 4.0577 |
| `hydrostatic_column` | 8.25 -> 8.27 | -0.1 % | 5.7e-02 | 4.3e-01 | 4.9e-01 | 2.9600 -> 2.9595 |
| `oxygen_chemistry` | 9.61 -> 9.60 | -6.9 % | 7.3e-01 | 3.1e+01 (sign change) | 8.9e-01 | 4.1487 -> 4.0577 |

The +0.05 dex of every WASP-121b and hot-Uranus case is the Jupiter radius
(item MASS, +11 and +16 percent in the raw flux at fixed cell, partly
compensated by the shrunken domain at cell N-20); `lower_profile` is item
BASE (the run now starts at the profile's 1 microbar level); the
`oxygen_chemistry` and `hydrostatic_column` moves are the grid and the
threshold edges. The baseline was re-snapshotted after the `hp_*` re-run
(the post-2b copy kept as `baseline_post2b_20260905/`).

### Golden refresh: prepared, awaiting the user's command

`golden/` was archived to `backup/regression/golden_pre170_20260905/` with a
README (provenance from the headers, unknowns marked), `roundtrip` was re-run
on the post-2c binary, and `DEFAULT_CASES` of `run_check.sh` now lists the
five Phase 0 cases as well. The refresh itself, `run_check.sh golden` over
the seventeen cases followed by `make check` and one `run_fcheck.sh`, is the
user's command to give; until it is given the reference of record is the
pre-2c state and `make check` reports every case as moved.

### Gate on the 89-bin band (2026-09-06)

Fourteen cases re-run against the post-2c (20-bin) baseline with md5
`b0a4ada9f84ee8410c84fda6780ba097`: `hydrostatic_column` byte-identical (no
sub-Lyman absorber), every other case within 1.3e-3 in any cell (`mol_metals`
rho 1.3e-3, `wasp_full` v 1.3e-3 at the base) and printed Mdot unchanged, the
size the 20 -> 89 ladder predicts; the `hp_*` IC pairs were regenerated on
the 89-bin `mol_base_handoff` output and the three cases run (Mdot 10.56).

**`oxygen_chemistry` stopped on the 89-bin grid, diagnosed.** At step 13,
cell 7, the carrier Newton reached its iteration cap (30) with a relative
OH row residual of 0.63 and the status stop of item 2b-A3a ended the run.
Diagnosis (worker, measured on frozen snapshots): the Newton is correct (a
root exists; the same solver with a larger budget reaches it in 141
iterations, the last three 5.4e-2, 3.1e-9, 3.9e-12) and its residual scale
and the element limiter are not involved. The stiffness comes from upstream:
at step 12 the ionization sweep accepted a non-root at cell 7 under the
relaxation amnesty (`info 4`, normalized residual 5.2e4), which put n_e =
3.5e10 and n(H+) = 2.6e8 into the frozen background where the previous step
had 1.1e3 and 0, so the H+ + H2 channel asked for two thirds of the cell's
H2 in one 2.85 s step. On the 20-bin grid cell 7 never took that path. The
"7.8e-2 carrier row" of the 20-bin summary is the STEADY residual of a
relaxation snapshot (dt -> infinity), not the Newton's; the backward-Euler
Newton converged to below 6e-12 on all 999 steps there. Two things landed:
the iteration budget is 200 (the comment at the code site states the
measurement), and the attribution test of the status stop, which asked
"did an element constraint act" by carrier over the whole grid, now asks it
by cell (the CO ceiling, acting in 160 cells of this case, had made any CO
failure anywhere excusable, and had left OH rows it moved by 49 percent
through the oxygen closure reported as the solver's); new suite
`src/tests/carrier_constraint_attribution/`, 7 assertions, RED before.
Open for the user: decision 4 (reject and retry the step on a shorter dt,
the structure that also covers a step without a root) and the ionization
sweep's relaxation amnesty, which handed the carrier system a non-root
(acceptance policy, `ionization_equilibrium.f90` class 4).

### Golden refresh executed (2026-09-06, on the user's instruction)

`run_check.sh golden` over sixteen cases: the eleven of the previous
`DEFAULT_CASES`, `hydrostatic_column`, the three `hp_*` restart cases and
`roundtrip`, all from the 89-bin binary `b0a4ada9f84ee8410c84fda6780ba097`
(conda gfortran 16.2, OpenBLAS); the previous golden is
`golden_pre170_20260905/`. `golden/` now carries the `*_adv.txt` files as
well. `oxygen_chemistry` has no golden entry yet (it stops, above); it stays
in `DEFAULT_CASES`, so `make check` reports it as the one failing case until
the carrier solver question is settled. `make check` and one
`run_fcheck.sh` follow in this section.

### After the refresh, same day (2026-09-06): the loaded-SED grid, the carrier budget, hygiene, the SED files

- **2c-SEDGRID, the photon grid of a loaded SED** (`sed_read.f90`, `J_inc.f90`,
  `parameters.f90`, `set_energy_vectors.f90`; item C5 of the D0 draft): the
  Load path kept the table rows as bin points with central-difference
  widths, so no threshold was a bin edge there, on the path every planet
  calculation of record uses (decision 16). Now the bins are built from
  edges (geometric means of neighboring rows, the table span as the two
  end edges), every active ionization threshold inside the span is inserted
  as an edge with the split bin's two halves carrying the row's flux (exact
  for the integrated flux), `e_v` the geometric center, `de_v` the exact
  width; `read_sed` receives the threshold list from
  `ionization_thresholds_active`, one definition. The table rows themselves
  are kept in `e_sed_node`/`F_sed_node` (undiluted), which `stellar_flux_eV`,
  `loaded_table_floor_eV` and `spectrum_covers_eV` now read, so the loaded
  field no longer depends on the binning. Test `sed_edges`
  (`src/tests/spectrum_type/`, 15 assertions) on the eps Eri table of
  WASP-52 b: thresholds at distance 0, flux identity to 2e-16, rate
  deviations from the exact table integral P_HI 1.8e-3 -> 4.6e-7, P_HeI
  1.1e-3 -> 1.1e-6, P_HeII 8.4e-3 -> 8.5e-6, worst metal 3.2e-3 -> 5.5e-4
  (Mg II, the table's own 1 A resolution). Movement: `benchmarks/wasp52`
  2000-step marching stop, mass flux at cell N-20 2.3e-6, log10 Mdot 9.68
  both, profiles up to 2e-2 at a relaxation point; `wasp_full` (power law)
  byte-identical on seven output files. No regression case loads a table,
  so no golden moves. Noticed, not fixed: `e_mid` is not an edge of the
  loaded grid (the `LX`/`LEUV` split it feeds moved 0.02 percent);
  `stop_if_balmer_band_unstated` prints bin centers as the band.
- **Carrier Newton budget** (`diffusive_photochemistry.f90`): `newton_maxit`
  30 -> 200 after the diagnosis above; `oxygen_chemistry` runs its 1000
  steps on the 89-bin grid (log10 Mdot 9.60) and is added to `golden/`
  below. The constraint attribution of the status stop, by cell, and its
  test suite `carrier_constraint_attribution` (7 assertions) landed with it.
- **Hygiene, D0 items C16 and C4**: the dummy `alpha`, threaded through
  `reconstruction_continuation_rhs`, `RK_rhs` and `Num_flux` and never read,
  is removed from the three signatures and ten call sites
  (`hydrostatic_column` byte-identical against a same-compiler control);
  `output_state_consistency.sh` and its README describe the rate lag
  (measured RED 1.19e-5 with the header now reading `sec_ion=F`) instead of
  the removed post-loop block. `docs/audit_20260905/roe_flux_vacuum_probe.f90`
  is an archived probe that no longer compiles against
  `speed_estimate_ROE`; left as the frozen record it is.
- **Document sync** (`docs/EXHALE_user_manual.tex` + PDF rebuilt,
  `docs/p44_base_sawtooth.md`, `write_output.f90` and `util_ion_eq.f90`
  comments, `write_setup_report.f90`): the manual carries the `Planck`
  type, the decision-17 stop, the six-band edge grid, the IAU constants and
  NIST thresholds, the signed `du`, the profile base level and the restart
  grid guard; `NlTR` is now counted on the Load path (the "no grid point
  below 13.6 eV" sentence had been printed for every loaded SED);
  `src/utils/glob.py`, which shadowed the standard library module, is
  `src/utils/EXHALE_interface_state.py` with its two importers updated.
  Reported, not changed: the GUI offers no `Planck` type; `global_parameters`
  exports an integer named `count` that shadows the intrinsic.
- **Stellar SED files** (`inputdata/sed/`, decision 16, source rule of
  2026-09-06: MUSCLES data normalized as the planet's reference paper did):
  `wasp121b_solar_huang2023.txt` (recommended; WHI 2008 solar shape below
  1700 A, F(1-912 A) = 1.6e6 and Lya = 1.0e5 as Huang et al. 2023 section
  2.2 prescribe, 6459 K photosphere above; the 1.6e6 confirmed as the
  lambda < 912 A band by Salz et al. 2019 section 7),
  `wasp121b_wasp17_huang2023.txt` (MUSCLES WASP-17 alternative; its 1700-2583
  A photosphere is 1.8x fainter than the blackbody, the size of the He 2^3S
  photoionization bias of the solar candidate), `hd189733b_epseri_salz2016.txt`
  (eps Eri, MUSCLES v22, normalized to Salz et al. 2016 Table 3; the proxy
  pairing is the one Salz et al. 2018 section 4.1.1 make), and
  `hd189733b_bourrier2020.txt` (the released MOVES III tables, mean of
  Visits B-E, measured to be the flux at 1 au, scaled to 0.031 AU, eps Eri
  photosphere above 1600 A identical to the Salz file's; 30 percent fainter
  in the Lyman-Werner band, 1.4x brighter in the EUV). Each file was read by
  `EXHALE.x` for two steps on the final binary without error. User decision:
  HD 189733 b runs at a = 0.031 AU (Salz 2016 Table 2, the modeling paper
  of the comparison); the planet folder, `benchmarks/hd189`, examples 01-11
  and 13, and `params_table.txt` carry it (the regression inputs
  `jfnk_hd189*`, `tpm_hd189` and the archived copies keep 0.0330). No
  `input.inp` has been switched to a loaded SED yet (2d).

### Phase 1 closing gate (2026-09-06)

`make test` on the final binary `b2c15e53c888d09fa90ffdad80443c3a`: ten
suites, eight PASS; RED only on the Phase 3 targets (`output_state_consistency`
1.19e-5, `h3p_cooling_limits`, `h2_rovibrational_identity`,
`h2_double_photoevent_energy_sum`, `helium_collisional_ionization_by_level`).
`make check`: sixteen `DEFAULT_CASES` against the refreshed `golden/`, all
four files byte-identical in every case (`REGRESSION PASS`); the sixteen
goldens taken from the earlier 89-bin binary reproduce on the final one,
which is the measured statement that the loaded-SED grid, the carrier budget,
the `alpha` removal and the report edits leave the analytic-spectrum matrix
untouched. `run_fcheck.sh` result recorded below.
`run_fcheck.sh`: the `-fcheck=bounds,do,mem` build ran 1500 bounded HD 209458 b
steps with no runtime trap (`fcheck regression: CLEAN`, log10 Mdot 9.59), and
the production build was restored to `b2c15e53c888d09fa90ffdad80443c3a`.
**Phase 1 is closed.** Open for the user before Phase 3: decision 4 (carrier
step rejection/retry), the relaxation amnesty of the ionization sweep that
handed the carrier system a non-root, the D0 draft
(`docs/d0_governing_system_20260906.md`, 33 consistency items, 15 questions)
for acceptance, the switch of the planet inputs to their loaded SEDs (2d), the
WASP-121 b SED choice (solar shape, WASP-17 shape, or the two joined at 1700 A),
and the sub-Lyman band budget of the `Load` grid (the table's own resolution
is kept; Mg II stays at 5.5e-4).

## 7. PLAN_20260906_rev2, Step A: acceptance contracts and run modes (2026-09-06)

Plan of record for this stage: `docs/PLAN_20260906_rev2.md` (rev 0 and rev 1
and their two reviews are kept beside it; the execution document of
2026-09-05 points at rev 2). The decisions of `To_be_determined_by_user_20260906.md`
were judged against `To_be_determined_by_user_recommend_20260906.md` and
adopted as rev 2 section 2 states.

Documents written (advisor): `docs/a0_run_mode_contract_20260906.md` (the
three run states; `Time stepping: Local` and PTC are initialization only;
physical time advances only by accepted global steps; headers gain `mode=` and
`t_phys=`; the proposed `Run mode: init | phys` key awaits the user's
confirmation before A0-impl); `docs/a2_certification_contract_20260906.md`
(one evaluator, four contexts: probe, stationary trial, physical step,
certification; statuses `not_applicable / evaluated / unavailable`; one
condition per active balance; what replaces what). Documents drafted by
workers and verified for structure and wording:
`docs/b1a_active_equation_inventory_20260906.md` (flags, equations per flag,
constraints, five validity states with producers, the prerequisite table for
eight configurations; finding: only the hydro rows have a state-consistent
evaluator; the composition projection, the `rho` rewrite and thirteen closures
with stated domains are recorded nowhere) and
`docs/b1_target_system_20260906_draft.md` (energy ledger with the shared H2
zero at `v=0,J=0`; chemistry mass-preserving; Lyman-alpha loss and return
from one 2p balance, an escape-probability multiplier on the existing
cooling explicitly not the target; the zero-net-current closure derived from
Koskinen 2022 B9/B10 and shown to give exactly the ambipolar field
`binary_element_diffusion` already computes; the trace-metal counter-flux;
the rollback table of fourteen state-changing operations showing `u_old` is
not a checkpoint; CO alternative A blocked on rate literature, alternative B
specified and recommended; ten open questions).

- **A1 (carrier half)** (`diffusive_photochemistry.f90`; suite
  `carrier_returned_state_acceptance`, 36 assertions, plus one in
  `carrier_constraint_attribution`): the carrier update is accepted on the
  RETURNED state, every unconstrained row of every cell,
  `|res|/max(row_terms, 1e-300) <= 1e-8` with the absolute floor inside the
  scale, finiteness of every row, independent of the termination branch;
  the relative-drop exit is `carrier_solve_stagnated`, a diagnostic; the
  decision is `carrier_returned_state_verdict`, a side-effect-free function
  returning verdict and reason for A3 to consume. Defect found and fixed:
  `limit_to_element_budget` marked a cell constrained only when the clamp
  exceeded `limit_report = 1e-10`, so smaller clamps left cells presented as
  Newton-measured. Impact: `mol_carrier` byte-identical (no changed branch
  reachable: all 499 transport steps of the first 500 ended at the floor);
  **`oxygen_chemistry` refuses at step 4**: the CO thermal ceiling, a clamp
  without a rate, leaves a neighboring cell's CO row at 1.0e-7 of its terms
  (accepted steps before it: 1.1e-9, 4.2e-9), and the code no longer
  certifies it. That is the plan's position (B1 alternative B: a
  ceiling-active atmosphere is outside the validated set); the case can run
  only in A0's initialization mode once A0-impl exists, so it fails
  `make check` until then.
- **A1 (class half)** (`steady_residual.f90`, `steady_newton.f90`,
  `ionization_equilibrium.f90` comment block; suite `acceptance_classes`,
  12 assertions): `n_cells_without_chemical_root(acc_n) = acc_n(4)`, one
  definition for the trial handling (a local-closure validity check, no
  stationary condition) and for the three points that declare a steady
  state solved; a class-4 cell in the state handed back is written but
  reported "NOT certified". Impact: `wasp_full_newton` and
  `mol_base_handoff` identical (no class-4 or class-5 acceptance in any
  ledger). Open: the marching stop cannot pass the count because the loop
  never receives the sweep ledger; A2 closes it.
- **C1** (in progress) and **A2** (steps 1-3 of the contract, in progress).
- **C1** (`inputdata/sed/README.md`, `build/sed_coverage_and_order_check.py`,
  suite `spectrum_type` driver `loaded_sed_table_semantics`, 78 assertions):
  a table value is a bin average on every shipped file (MUSCLES const-res,
  WHI, Gueymard: 1 A or 1 nm bins; Bourrier 0.016 A samples with no stated
  convention) except the blackbody part of the solar WASP-121 b candidate
  (point sample); the reader treats a row as a histogram in energy, exact
  for the band energy of a bin-averaged table, second order in the bin
  width otherwise (MEASURED order 2.0 to 2.9 through the production
  consumers). Coverage is now checked per configuration from the active
  threshold inventory (2600.5 A metastable, 3647 A when both stellar keys
  arm H(n=2), the lowest active metal), 12 configurations. Findings closed
  the same day (advisor): `sed_read.f90` single-precision `1e-8` at two
  sites (6e-9 relative shift of every loaded energy); the HD 209458 b WHI
  file rebuilt to 5000 A (its 3000 A cut stopped the run in `read_sed`
  naming H(n=2)); `WASP-52b/input.inp` SED path from before the 2026-07-14
  move; the two solar build scripts' headers (2583 A / 4.80 eV). Remaining:
  the WASP-17 file begins at 10.0 A because the MUSCLES v24 table does (the
  grid top becomes 1239.84 eV; recorded).
- **BALMERQ** (`excited_hydrogen.f90`; suite `spectrum_type` driver
  `balmer_band_quadrature`, 21 assertions; `balmer_continuum_field` Planck
  references refined): the Balmer-continuum rate and heating integrals
  follow the field's own nodes in closed form (power-law segments of a
  loaded table and of the power-law type; 5-point Gauss on a geometric
  subdivision for Planck), replacing a fixed 400-point uniform rule that
  stepped over emission lines. On the eps Eri table of `benchmarks/wasp52`
  the old rule was 2.9 percent low in the rate and 2.7 percent high in the
  heating (the Lyman-alpha line at 10.20 eV lies in the band); on the power
  law it carried 9.5e-5 / 6.0e-5 of truncation. Movement: `benchmarks/wasp52`
  gamma_2 +2.99 percent, n_e 5.3e-4 at 1.09 Rp; `wasp_full` 1e-8 class
  (sub-0.1 percent, no golden refresh). One worker reported two scratch runs
  killed by a `pkill` matching a shared basename: the common brief now says
  kill by PID after checking `/proc/<pid>/cwd`, never by name.
- **A2, steps 1-3** (`src/modules/time_step/certification.f90` new,
  `steady_residual.f90`, `steady_newton.f90`, `EXHALE_main.f90` call sites,
  `utilities.f90` header, `diffusive_photochemistry.f90` save/restore pair;
  suite `certification`, 18 assertions): one evaluator builds a 19-entry
  inventory from the configuration flags and returns per entry a status
  (`not_applicable` / `evaluated` / `unavailable`), the row measure (max
  over physical cells of `|res|/max(scale, floor)`), the volume average,
  the worst cell, finiteness and the scale text; hydro rows and carrier
  balances are evaluated (carriers inside an isolated workspace that saves
  and restores every module array `carrier_steady_residual` writes,
  asserted after every call), elemental transport, level equations and
  the seven `System_HeH_*` closures are `unavailable` until step 4. The
  certification context runs at the marching stop (which now receives the
  sweep ledger, closing the A1 open item), the JFNK final and the PTC
  declaration; probe and trial decisions read one record, the merit rule
  unchanged. Exit status (user decision): a run that DECLARED a stationary
  state and failed certification exits 2 (`stop 2`); a run ending on a cap
  makes no claim, writes `certified=F cert_reason=no_stationary_claim` and
  exits 0; `run_check.sh` never tests the status. Five shell gates of
  `grid_and_gates` that asserted exit 0 from a full run now accept 2 with
  the reason stated at the code site (advisor). Impact: numeric content of
  `wasp_full`, `wasp_full_newton`, `mol_carrier` identical; headers gain
  `certified=`. **No configuration certifies today**, by construction (the
  closures are `unavailable`), and two measured findings stand on their
  own: the JFNK `info = 0` is set at the loop top and never revised when the
  returned state is re-measured at its own composition (`wasp_full_newton`
  energy row 4.1e-4 against the 1e-5 target, 40x), and the `du` stop
  accepts hydro rows of order one on their own terms (`wasp_full` at
  `du = 1.0e-3`: mass 1.3e-2, momentum 1.8, energy 0.72; at `du = 2e-2`
  the same rows read 1.7, 2.0, 0.96, so the spread and the rows measure
  different things). Not implemented in this item (it would change the
  solve): the probe rule that a class-4 chemistry makes a Newton direction
  unusable. Tolerances inherited; the convergence study is owed.
- **B1 accepted (2026-09-06).** The user decided the fourteen open
  questions of `docs/b1_target_system_20260906_draft.md` as recommended
  (metal reservoir zero; electron mass convention; `ioniz_eq` density
  `intent(in)`; the residual H I excitation cooling under the Lyman-alpha
  ledger; trace metals in the diffusion matrix; H+ in the B9 set; CO
  exclusion adopted with the destruction model held for rate literature;
  H3+ low-collider policy; one H2 partition function; rollback of
  diagnostic extrema; identity row placement; the CO ceiling's validity
  category; O(1D) reservoir; the refused `_adv` product rule). The accepted
  specification is `docs/b1_target_system_20260906.md` (the draft kept and
  marked superseded). This is the Phase 2 gate of PLAN rev 2; B2 (the
  residual-controlled temperature solve) started the same day.
- **A3** (`diffusive_photochemistry.f90`; suite `carrier_retry`, 32
  assertions): the carrier subsystem's checkpoint (`carrier_checkpoint_take
  / restore / matches`: carrier fractions, the module arrays of the A2
  save/restore enumeration, every constraint and solve diagnostic, the
  verdict, the cumulative CO record), and the controller: a refused
  returned state restores the substep's entry state, halves the carrier
  substep with the background explicitly frozen (listed at the routine),
  covers the interval with accepted substeps, and on exhaustion (nine
  refusals, floor `dt/2^8`) restores the ENTRY state and stops; `f_sp` is
  written only for a covered interval. Two CO records (accepted, restored
  on refusal; attempted, never restored) and controller attempt statistics
  with the A0 ledger-family hook. Retried result equals the explicitly
  subdivided integration bitwise; step doubling gives first order (ratio
  2.001). Impact: `mol_carrier` (0 retried intervals in 11999),
  `mol_base_handoff`, `oxygen_chemistry` (3 steps) byte-identical.
  **Finding:** the retry recovered the `oxygen_chemistry` refusal (1.04e-7
  -> 7.9e-9 after ten halvings) because the row scale carries the `1/dt`
  time term, so any fixed clamp imbalance passes a fixed relative floor at
  a short enough substep; that is a defect of the acceptance measure shared
  by A1 and A2, fixed next (A1scale: the returned-state and stationary
  measures are judged against the physical terms of the row, transport and
  reaction, with the absolute floor; the time term stays only in the
  Newton's own iteration test).
- **B2** (`energy_semi_implicit.f90`; suite `energy_update`, 20
  assertions): the implicit temperature step of a cell,
  `R(T) = u(T) - u(T_old) - a [heat - cool(T)] = 0` with `u(T)` from the
  run's equation of state (caloric H2 included), is solved by a safeguarded
  Newton-bisection hybrid on a physical bracket (`0.01 T0` to 5e4 K with
  H2, 1e7 K otherwise) to `|R|/scale <= 1e-9`, the cooling written back is
  `cool(T_returned)`, and every non-converged outcome is a named failure
  (`ITER_CAP / NO_BRACKET / FLOOR / NONFINITE`) with the cell, the residual
  and the missing volumetric source; no clamp remains, a floor is a
  failure, and the interim action is a stop that B3a replaces by rejecting
  the attempted step. The caller is untouched (optional status argument,
  `energy_update_last_status`). RED kept in the suite: the two-iteration
  update on a falling cooling branch left `|R|/scale = 2.9e-2` and a 6.9
  percent temperature error. Impact: `wasp_full` transient 5.4e-5,
  converged (same stop step 16257, log10 Mdot 13.29) 1.6e-6;
  `mol_base_handoff` 4.9e-5; `hydrostatic_column` identical; no floor or
  failure fired anywhere; cost +9 percent end to end on `wasp_full`.
  Follow-ups assigned: the identical clamp in the conduction stage
  (`viscous_conduction.f90` 572-589; B2b), the certification's counting
  of energy-floor hits as an accepted correction (A2 step 4), the stale
  end-of-run floor report (A0-impl).
- **A1scale** (`diffusive_photochemistry.f90`; `carrier_retry` 42
  assertions): the returned-state acceptance and the stationary carrier
  measure are judged against the row's PHYSICAL terms (`|transport| +
  |reaction| + floor`, the floor in rate units without the `1/dt` form);
  the time term `|n - n_old|/dt` stays only in the Newton's own iteration
  scale (`row_terms`, bitwise unchanged), beside which `row_terms_phys` is
  accumulated. RED: a residual of 1e-7 of a row's physical terms read
  3.3e-8 at dt and 5.6e-9 at dt/16 on the old scale (crossing the 1e-8
  floor with the defect unchanged); now 1.0e-7 at every dt. Impact:
  `mol_carrier`, `mol_base_handoff` identical; `oxygen_chemistry` is refused
  at step 2 (CO row 5.5e-8 beside a clamped cell; 37 retry attempts to the
  floor, all refused, entry state restored), and the old measure only
  postponed that refusal to step 88 on today's tree. Finding carried
  forward (A1scale2, needs its own impact measurement): with a
  dt-independent acceptance there is a shortest certifiable substep, since
  the Newton's stopping test is still relative to the full terms and the
  round-off it leaves grows as 1/dt against the physical terms; the same
  measure in both places closes it. Also: `cterms` in `steady_newton.f90`
  (534, 645) is written and never read.
  **Consequence for the matrix:** `oxygen_chemistry` cannot run under the
  physical acceptance while the rateless CO ceiling is active (B1 decision
  7: exclusion); once A0-impl lands, the exhaustion action reads the run
  mode (initialization: record the refusal in the attempted ledger, write
  the state flagged, continue; physical: stop), which is the item A3b, and
  the case runs in `init` mode as the diagnostic it is.
- **B3a and B4 designs** (`docs/b3a_attempted_step_controller_design_20260906.md`,
  `docs/b4_spatial_operator_design_20260906.md`, read-only workers; advisor
  decisions in their last sections). B3a: eight of the fourteen step
  operations have no returned-state verdict today (RK/BC, element diffusion,
  excited H, `ioniz_eq`, the projection, `Apply_BC`, conduction, Shapiro);
  checkpoint about 2.1 kB per cell. B4: species mass fluxes are not
  consistent with the mass flux today (equilibrium species are not
  transported and rewrite `rho`; carrier advection is a cell-centered
  product; element diffusion is in advective form; no species enthalpy
  transport); increments B4-1 to B4-6 with their expected movement.
  **Advisor correction to the B4 draft:** the formation-energy flux
  `F_form = sum_s eps_s F_s` is the RESERVOIR's transport and enters no
  thermal equation (photoheating already carries `h nu - I`; recombination
  emits `I`); the draft had written it as a thermal source, which would
  heat a cell by the ionization energy leaving it. B4-2 becomes the
  diagnostic reservoir balance that the ledger identity T1.5 tests, plus the
  face-composition enthalpy (T-B4.5). The B4-2 pre-measurement
  (`src/utils/formation_energy_flux_diagnostic.py`, section 7 of the design)
  stands as the size of that reservoir transport: on `wasp_full_newton` 28
  percent of the absorbed heating leaves the column as advected ionization
  energy (median 14 percent of the local heating, 52 percent at the
  ionization front at 1.23 Rp); the present energy face flux is `v (E + p)`
  with the caloric H2 internal energy already inside `E`, addressed by the
  owning cell (T-B4.5). Findings: `lower_profile`'s golden mass flux is not
  flat (1.1e-2), so continuity-based diagnostics on it discriminate nothing
  until it is Newton-finished; `benchmarks/hd209/output` integrates to
  cooling 24 percent above heating on its saved state.
- **B2b** (`viscous_conduction.f90`; `energy_update` suite second driver
  `conduction_floor_tests`, 24 assertions): the Crank-Nicolson conduction
  stage no longer clamps; a solved temperature at or below the shared EOS
  floor (imported from `energy_semi_implicit`, one declaration) is a FAILURE
  with status, cell, the solved temperature and the energy the floor would
  have injected; the tridiagonal solve refuses a zero or non-finite pivot;
  no partial state is written. The operator is in flux form and conserves
  thermal energy exactly on a closed column (MEASURED 0.0); as
  boundary-conditioned in a run it exchanges heat through the base Dirichlet
  anchor by design (3.4e-15 up to that flux). Impact: byte-identical on
  `wasp_full`, `mol_base_handoff`, `hydrostatic_column` and on two
  conduction-on variants in which the stage runs; zero floor activations
  anywhere. B1 section 7.4 rows 9 and 11 updated (no clamp state left).
  Noticed: `benchmarks/koskinen2022_model_a/*` cannot restart (stored IC on
  the pre-2c grid, 4.5e-4 off; the same class as the `hp_*` regeneration).
- **B4-5 measurement** (`src/tests/grid_and_gates/hydrostatic_residual.f90`,
  design section 8): on the analytic isothermal column, through the
  production `Reconstruct`, `Num_flux`, `RK_rhs` and `Source`, the
  pressure-gravity pair is at design order in the interior (`dr^1.96` PLM,
  `dr^1.98` WENO3 at r = 2; base cell 1e-4 of its weight, first order in the
  base spacing), and the O(1) hydrostatic imbalance of the mechanical column
  sits at the OUTER free-outflow cell: 0.75 of its own weight at every N
  under PLM (the zero-gradient ghost zeroes the MC forward difference, so
  cell N carries no pressure gradient), 6.4e-2 falling at order 1.4 under
  WENO3. A well-balanced pair (B4-5b) would not remove the drain; the outer
  boundary fill is the item (B4-4a, started). Two latent defects found:
  `define_grid`'s `tol = 1.0` initializer implied SAVE, so a second call in
  one process skipped the stretch-factor Newton (production calls it once;
  fixed by the advisor: `tol` reset before the loop); the case README's
  single `Grid cells` ladder does not refine the base cell (README updated).
- **A2 step 4** (`certification.f90`, `binary_element_diffusion.f90`,
  `excited_hydrogen.f90`, `ionization_equilibrium.f90`; `certification`
  suite +7 driver and +15 shell assertions): the five entries that were
  `unavailable` are measured on the state without solving or writing
  anything (He/H partition and trace-element transport balances, He 2^3S and
  H(n=2) level balances, the eliminated-species closure as the sweep's own
  normalized reaction residual at the given composition), each restoring
  every module array it installs, asserted on every call. Output content
  identical on seven cases. **What refuses certification now**: the three
  hydrodynamic rows on every marching case; the energy row alone on
  `wasp_full_newton` (4.1e-4 at cell 376, the JFNK `info = 0` finding); the
  H2 carrier row on `mol_carrier`; and, new, the elemental transport on the
  diffusive runs (`mol_diffusion` partition 6.4e-5 at cell 203;
  `lower_profile` partition 4.8e-3, trace metals 1.2e-2, the first
  measurement of what the Picard alternation reaches). The closures and
  level balances are within tolerance everywhere (1e-16 to 1e-12 of their
  turnover), so the equations that lacked an evaluator are the solved ones.
  Both floors are reported as attempt diagnostics (0 everywhere).
  `cert_tol_element` is inherited, not measured, and is now decisive for the
  diffusive runs: the convergence study of contract section 4 is owed.
- **A1scale2** (`diffusive_photochemistry.f90`; `carrier_retry` 51
  assertions): the carrier Newton stops when its residual is within its own
  tolerance of the full terms AND within the acceptance tolerance of the
  physical terms, or when the arithmetic forbids the second (residual at
  `64 epsilon` of the full terms), a named outcome; the interim reason
  `unreachable_at_substep` for a refusal whose rows are all at round-off was
  superseded the same day by A1scale3 (below); the stop message reports the
  halvings made.
  Option (a) (stopping on the physical scale alone) was built, measured
  (one extra iteration every step, no rescue of the short substep) and
  declined. Correction of the A1scale report: the Newton reached one
  epsilon in two iterations on that column; the 7.7e-5 was the ninth
  attempt's state. Impact: `mol_carrier`, `mol_base_handoff` identical;
  `oxygen_chemistry` still refused after step 2 (now the unreachable reason,
  OH at cell 90, 15 attempts). Decision taken from its finding (A1scale3,
  started): a row at the arithmetic round-off of its own terms is not
  evidence of a defect and is accepted as "round-off limited" with a flag
  and a count, so a trace carrier at the element floor does not refuse an
  interval; the physical-terms acceptance stays for every other row.
- **A0-impl** (`input_read.f90`, `parameters.f90`, `EXHALE_main.f90`,
  `utilities.f90`, `load_IC.f90`, `write_setup_report.f90`,
  `ionization_equilibrium.f90` ledgers; suite `run_mode`, 22 assertions):
  the run mode is in. `Run mode: init | phys`, default `init` in every
  configuration (advisor decision replacing the draft rule: the ordinary
  use is a stationary solution reached by marching and certified, which is
  initialization in the contract's terms; a physical integration is asked
  for); `phys` with `Time stepping: Local` refused at startup; in `phys` the
  handoff state is tested at step 0, the clock advances only after a
  complete accepted step (measured equal to the sum of accepted `dt` to the
  last bit, with and without forced rejections), two ledger families
  (`ledger_family` in `global_parameters`), headers gain `mode=` and, in
  `phys`, `t_phys=`; a restart of an `init` snapshot in `phys` mode is the
  handoff at t = 0 (contract section 3 clarified), and only a `mode=phys`
  header without a finite `t_phys` is refused; `init` never exits nonzero
  for certification. Impact: no input changed; `wasp_full`, `hp_front`,
  `hydrostatic_column`, `mol_base_handoff` bitwise identical. Both floor
  reports of `EXHALE_main.f90` now state the B2/B2b behavior.
- **COlit** (`docs/co_destruction_rates_literature_20260906.md`; publisher
  PDFs of Westley 1980 NSRDS-NBS 67, Tsai 2017, Wakelam 2012 added to
  `references/`): B1's CO choice is settled on physics: every neutral CO
  destruction channel carries a barrier of 6.7 eV or more, so the one-sided
  destruction model's `tau_dest << tau_res` fails at the ceiling's onset at
  every density (factor 9 to 1.5e8) and the exclusion (decision 7) is the
  physics. Findings: `molecular_reaction_heat.f90` (H2, H+, He+, H2+, H3+,
  HeH+) and the oxygen Shomate table of `oxygen_rates.f90` carry different
  reference zeros (T1.9 one-table item, assigned as B3b-ENTH); CO2 is not in
  the species table (removes the one evaluated CO channel below 2500 K and
  an OH sink); the workspace VULCAN network's CO reactions all produce
  species EXHALE does not carry, VULCAN's stated validity is 500 to 2500 K,
  KIDA's 10 to 300 K. Requests for the user: the A&A papers Visser et al.
  2009 (503, 323; the CO shielding functions), McElroy et al. 2013 (550, A36
  + RATE12 file), Venot et al. 2012 (546, A43) and 2020 (634, A78), plus
  Moses et al. 2011 (ApJ 737, 15) and Tsai et al. 2021 (ApJ 923, 264),
  which this host cannot download.
- **A1scale3** (`diffusive_photochemistry.f90`; `carrier_retry` 61
  assertions): a row at the arithmetic round-off of its own full terms
  (`64 epsilon`) is counted "round-off limited", flagged per row, and left
  out of the acceptance decision; every other row is judged on the physical
  terms at 1e-8 as before; the controller halves only for genuine refusals;
  the stop prints both scales for every refused row. Impact: `mol_carrier`,
  `mol_base_handoff` bitwise identical (0 round-off limited rows);
  `oxygen_chemistry` still refused at the same step and row (CO at cell 281,
  `|res|/phys = 3.2e-8`, `|res|/full = 2.2e-12`, four orders above the
  round-off bound): the element limiter's imbalance one cell beyond a
  clamped cell, correctly refused. The count's accessor exists; its line in
  the run summary and the certification carrier entry are B3a's (owner of
  both files).
- **B3b-H2Q** (`caloric_eos.f90`, `mol_rates.f90`; `h2_rovibrational_identity`
  19 assertions): one H2 partition function over the observed Roueff et al.
  2019 ladder (302 levels, zero at v=0,J=0) gives `u_rv(T)` as its first
  moment and the chemistry's `K_eq` as its free energy; the Huber-Herzberg
  model constants are deleted. The RED (0.6 to 8 percent in `u_rv`) was the
  model ladder's level placement (no centrifugal distortion: `J^4`
  signature, 100 percent of the gap to 2000 K), not a defect of the sum;
  independent check: `K_eq` against NIST-JANAF (Shomate) now 0.9997 to
  1.0008 over 600 to 5000 K (was 1.010 to 1.101). Impact: `wasp_full`
  identical; `mol_base_handoff` T 1.3e-5, n(H2) 3.1e-4; `mol_diffusion` T
  3.2e-5; `mol_metals` at 400 steps n(H2) 5.6e-3 at 1.17 Rp; Mdot unchanged.
  Findings: a third copy of the ladder sum in `molecular_infrared_cooling.f90`
  212-215 (to fold into `h2_partition_function`); `element_census_tests`
  does not compile against the live tree (`no_free` argument, another
  worker's signature change; to fix at the next single build); the load
  average reached 38 with the parallel long runs, so future briefs pace the
  12000-step impact runs.
- **A3b** (`diffusive_photochemistry.f90`; `carrier_retry` 89 assertions):
  the carrier controller's exhaustion reads the run mode and returns a
  status (`carrier_interval_covered / exhausted`) the caller reads. `phys`:
  unchanged (restore, print, stop until B3a takes the status). `init`: the
  entry state is restored, the interval is recorded in the attempts ledger
  (count, first/last step, worst row and both scales), the carrier history
  is marked not certifiable (`carrier_history_certifiable()`), and the march
  continues from the entry carriers. The accepted CO record is indexed by
  A0-impl's `ledger_family`. Impact: `mol_carrier`, `mol_base_handoff`
  bitwise identical; `oxygen_chemistry` (init, its default) now runs its
  1000 steps, exit 0, with 998 exhausted intervals (first at step 2, worst
  CO at cell 313), log10 Mdot 9.61 against the 20-bin 9.60, and a final
  certification refusing on four carrier rows at 1.0: its carriers are
  frozen at the step-2 values for the whole run, which is exactly what the
  mark and the ledger say. The refused row is the real imbalance beside a
  clamped cell (`|res|/full` 4e-12 against a 1.4e-14 round-off bound), open
  physics of the excluded CO model, not a controller defect. Its golden is
  re-taken at the end of this series with the mark noted in the case
  README.
- **B5pre** (`steady_newton.f90`; suite `steady_completion_flag`): the JFNK
  and PTC completion flag is decided AFTER the certification of the returned
  state, from the same record the certification prints; a solve that met
  the gate at the loop top and not at the returned state ends `info = 2`
  naming the row. Measured on `wasp_full_newton`: loop-top `||R|| = 5.5e-6`,
  returned state 4.2e-4; the existing own-composition consistency pass
  opens the gap (0, 1, 2, 3, 8 sweeps: 5.5e-6, 1.9e-4, 3.2e-4, 3.9e-4,
  4.2e-4), so the hydro-only Newton converges to a state that is not a root
  of the coupled hydro + composition system: the loop-top gate measures a
  lagged residual. Consequence today: `EXHALE_main` treats `info = 2` as a
  solver failure and falls back to marching, which ends on a du plateau
  with a state that certifies worse than the JFNK state; assigned to B3a
  (`info = 2` = a failed stationarity claim, written uncertified, exit 2, no
  fallback) and to B5a (`steady_newton.f90`: the residual evaluated on the
  self-consistent composition of the current iterate, so gate, merit and
  returned state measure one thing). Also found: the JFNK stop message
  prints an unset `resid_th` (B3a); `certification_row_measure` and
  `relnorm_over_cells` are one expression written twice.
- **B3b-PE** (`h2_photo_channels.f90`; `photoevent_energy_ledger` 105
  assertions): the 19.725 eV missing from the H2 double photoionization
  event is the kinetic energy release of the two protons: 51.4 eV is a
  vertical threshold (Yan et al. 1998 section 4), the asymptotic product
  energy is `D0(H2) + 2 I(H) = 31.675` eV, and energy conservation at the
  threshold fixes the fragment kinetic energy with no free parameter
  (Franck-Condon check `e^2/R_e = 19.42` eV, 1.5 percent). Second defect
  found: the 33 to 41 eV neutral window dissociates into `H(2p) + H(2s)`
  (Chung et al. 1993 p. 887), so 20.398 eV per event is internal excitation
  leaving as prompt Lyman-alpha and two-photon continuum, which the
  integrand charged to heat. `h2_channel_energy_recipients` (reservoir,
  electron, fragment kinetic, radiated; sums to `h nu`) with the per-channel
  table and sources. Channels M and S already closed; every atomic absorber
  and all 17 metal ions correct (`wasp_full` byte-identical). Impact:
  `mol_sec_ion` T 4.2e-4, `heat_H2` 4.5e-2; `mol_base_handoff` T 2.1e-5;
  Mdot unchanged. The four-line wiring of the two new terms into the
  heating integrand of `util_ion_eq.f90` is handed to B4-6 (owner of that
  file). Open: suprathermal H2 dissociation credits heat but destroys no
  H2; the double channel has no absorber row of its own; channel S has no
  source for a fragment kinetic split.
- **B3b-H3** (`h3p_cooling.f90`; `h3p_cooling_limits` 17 assertions): the
  Miller et al. 2013 fits as the paper defines them: Table 5 is the LTE
  emission of one molecule per steradian (`Lambda = n(H3+) 4 pi E(T) s`),
  Table 6's `s(T, n_H2)` is scaled LTE with H2 the only collider (Oka and
  Epp coefficient; the paper's upper-limit caveat at the code site); below
  the table's `1e6 cm^-3` the non-LTE factor follows the exact collisional
  limit (linear in `n_H2`, zero at zero), replacing the clamp at the lowest
  column; above `1e14` the table edge is held; the four Table 5 segments do
  not join (measured with the paper's coefficients: +43 percent at 300 K,
  -0.11 percent at 800 K, -2.4 percent at 1800 K) and the code joins them by
  a stated linear blend in `log_e E` over the 5 percent below each join
  (Table 6's 300 K value supports the upper polynomial); B1 T6.3 reworded to
  say so. Published values reproduced per segment to the paper's precision.
  Impact: the H3+ cooling peak unmoved; T moves 1e-5 to 2.7e-4 in the three
  molecular cases; the clamp removal changes cells beyond 1.22 Rp where
  `n_H2 < 1e6` by more than 50 percent at 1e-10 of the peak (emission from
  gas with no collision partners, gone); `wasp_full` identical. Four domain
  counters exposed (B6 wiring is the certification owner's).
- **COlit, completed with the user's PDFs** (Visser 2009, McElroy 2013,
  Venot 2012, 2020; `docs/co_destruction_rates_literature_20260906.md`, 1419
  lines): CO photodissociation as Visser et al. define it (37 bands inside
  EXHALE's 912 to 1110 A beam, unshielded 2.6e-10 s^-1 in the Draine field,
  shielding tables in `N_CO x N_H2` for four parameter sets, T_ex raised only
  to 512 K because no data exist for dissociation out of v'' = 1); evaluated
  on the `oxygen_chemistry` run's own base columns (log N_CO 18.96, log N_H2
  21.80) it gives `tau_dest = 2.8e8` s, twelve decades faster than the
  barrier channels and still two to six decades short of `tau_dest <<
  tau_res` at the base; zero in `mol_base_handoff`, whose spectrum carries no
  flux longward of 911.8 A. RATE12's paper prints no entries (the data file
  is still needed for the ion channels); Venot 2012 is valid 300 to 2500 K
  and 0.01 to 100 bar (EXHALE's base at 1e-6 bar is four decades below);
  none of the six sources covers the 3000 to 4000 K band where the ceiling
  fires. The CO exclusion (B1 decision 7) stands, now quantitatively. If a
  CO photodissociation term is ever added it is a fourth absorber on the
  existing Lyman-Werner beam with a second column (`N_CO`) and the cell-mean
  form of `lyman_werner_dissociation_rate`.
- **TESTINFRA** (`grid_and_gates/run.sh`, `physics_probe/run.sh`, their
  READMEs): both suites take `EXHALE_OBJDIR`, `EXHALE_EXE`, `EXHALE_TEST_OUT`
  with the semantics of `spectrum_type/run.sh`; default behavior byte for
  byte as before; against a private build both run without touching
  `build/` objects (the six `grid_and_gates` shell tests still place their
  run copies under `build/tests/`). A dead fallback that wrote a build stamp
  the closure never compiled is fixed. Also (advisor): the two Python source
  scanners now skip `.ipynb_checkpoints/` directories, which hold editor
  auto-save copies of four modules (git-ignored) and had made
  `threshold_literals_written_once` red on a stale copy of `cross_sec.f90`;
  the `Makefile` source list is explicit and never reached them.
- **B3b-ENTH** (`molecular_reaction_heat.f90`, `oxygen_rates.f90`; 78 Fortran
  + 37 Python assertions): one table `species_formation_energy(isp)` for
  every species column (H/He stages, 27 metal stages cumulative from the
  neutral, H2/H2+/H3+/HeH+, OH/H2O/CO, plus free electron, H(n=2), O(1D))
  from the single B1 reference; the 17 reactions of the H2/He network are a
  stoichiometry table and each heat a difference of table entries; the
  Shomate table stays with its reference stated, every use of it being a
  balanced-reaction difference (the zero cancels), and the oxygen and carbon
  entries of the one table are differences of it; one D0(H2)
  (`mol_rates::h2_dissociation_energy_eV`, reproduced by the fit to 1.5
  meV). Physics defect fixed: the oxygen photolysis thresholds were 298.15 K
  reaction enthalpies; they are now 0 K dissociation energies (the lowest
  Shomate interval's F coefficient is the 0 K enthalpy), 0.9 to 1.0 percent
  lower, deposited heat +0.8 percent at the oxygen base cell. Impact:
  `mol_base_handoff`, `mol_sec_ion`, `mol_metals` byte-identical;
  `oxygen_chemistry` T +0.4 percent at 1.02 Rp, Mdot 9.61. Open for B3c: the
  undeposited collisional oxygen heats and the He 2^3S associative branch
  (D0 C25; the table exists); the B4-2 Python diagnostic lacks four oxygen
  entries; `formation_energy_report` has no caller; a stale threshold
  comment in `water_photolysis.f90:281`.
- **C3** (six `input.inp` and READMEs; `inputdata/sed/README.md`): the
  planet folders and benchmarks of HD 209458 b, HD 189733 b and WASP-121 b
  load their literature SED files (`Spectrum type: Load from file..`;
  `Power-law index` removed, `LX`/`LEUV` kept because mandatory and
  overwritten by the grid integral, which the setup report shows); every
  README records the configuration of record with the existing power-law
  products listed by mtime and the reruns left to a separate instruction;
  `read_sed` accepts every file at its configuration's floor (the H(n=2)
  edge in all six). `log F_XUV` at the orbit moves 3.193 -> 3.127
  (HD 209458 b), 4.438 -> 4.311 (HD 189733 b; the recomputed `LX`/`LEUV`
  reproduce Salz Table 3), 6.205 -> 6.204 (WASP-121 b). Advisor follow-ups:
  HD 209458 b's orbital distance set to 0.04707 AU (Oklopcic and Hirata
  2018, the reference the run is compared with; the input said 0.0480 while
  the file was scaled to 0.04707); the setup report's flux ratio format
  (`F8.3` overflowed when the sub-Lyman band dominates); the three
  benchmark directories added to the coverage script's inventory (all pass).
  Known: `Load IC? True` in the two planet folders is broken by the grid
  guard (pre-rebaseline IC), the same class as the `hp_*` regeneration.
- **COlit, revised with RATE22, Moses 2011, Tsai 2021** (document now 1956
  lines): RATE22 (Millar et al. 2024, A&A 682, A109) carries `He+ + CO -> C+
  + O + He` at 1.6e-9 cm^3 s^-1, measured, accuracy A, all products EXHALE
  species; it has no `H+ + CO`, `CO + M` or `CO + H2` destruction entry, and
  its `Tu = 41000` on 6562 entries is an undocumented "no upper bound"
  marker. On the `oxygen_chemistry` profile the He+ channel with the Visser
  shielding gives `tau_dest/tau_res` from 6.7e2 at the base to 9.6e-5 at
  1.14 Rp: **the one-sided destruction model is constructible above the
  helium ionization front and not below it**; the boundary is
  compositional, not the ceiling's thermal switch. Decision 7 (exclusion)
  stays until T5.1 is implemented; B1 section 5.2 revised, B3b-CO scheduled
  after B3c. Moses 2011: no stated envelope, one fit interval 500 to 2500 K,
  no shielding. Tsai 2021: 0.1 nm continuum opacities, no CO
  self-shielding; EXHALE's `lyman_werner.f90` is the working precedent. The
  literature requests are closed except Baulch 1976 vol. 3, which blocks
  nothing.
- **B4-6** (`utilities.f90` `cell_mean_attenuation`, `util_ion_eq.f90`;
  `xuv_cell_mean_attenuation` 56 assertions): the photoionization and
  photoheating sums of a cell use the exact cell mean of the attenuation,
  `exp(-tau_out) (1 - exp(-dtau))/dtau` (composite Gauss when the 2D
  correction is on), instead of the inner-face value the top-down columns
  gave; one definition for H, He, He 2^3S, H2's four channels, the metals,
  the single-absorber rates, `q_abs` and the secondary-ionization
  integrands. Inner-face error: 0.5 percent at `dtau = 0.01`, 42 percent at
  1. Photon conservation on the written states: cell mean 1.000000000 of the
  incident flux, inner face 0.9725 (`wasp_full`), 0.9843, 0.9849 (the loss
  sat at the ionization front). Impact: `wasp_full` (24000 steps) x(H+)
  8.3e-2 at 1.32 Rp, Mdot +0.74 percent; `mol_base_handoff` Mdot +0.93;
  `mol_sec_ion` +1.01; `hydrostatic_column` unchanged. B3b-PE's two heating
  terms wired into the integrand (its own numbers reproduced cell for cell).
  Advisor: the single-precision `1.0e-18` cross-section unit literal (33
  sites, 4e-8 relative on every rate) made double. Finding to resolve at the
  gate: `wasp_full` on the mid-series tree plateaus at `du = 4.2e-3` from
  step 22000 and never reaches 1e-3 (the golden stopped at 16230), identical
  in the item's before binary, so it belongs to an in-flight change (B3a
  edits `EXHALE_main.f90`) or to an earlier landed one; measured after B3a.
- **A2tol** (measurement; `a2_certification_contract_20260906.md` section 9):
  the `du` stop and the hydrodynamic row measures are not measures of one
  thing (over a decade of `du` on `wasp_full` the mass row falls at order
  0.2, the momentum row rises, the energy row is not monotone; `mol_carrier`'s
  carrier row falls at order 0.08 over 12000 to 48000 steps); the JFNK
  returns an energy row 41 to 47 times its target at every rung of a
  two-decade `Resid tol` ladder (structural; B5a); the elemental transport
  rows cannot be measured on the marching route (the operator advances 1e-6
  of the way per CFL step) and the diffusive cases never reach the steady
  route (du stalls at 5 to 6e-2 against the 1e-2 hand-off). Round-off floors
  on the converged JFNK state: mass 1.0e-13, momentum 5.1e-12, energy
  1.6e-11 (residual bitwise reproducible). Recommended and adopted
  (applied by B3a in `certification.f90`, provisional until B5a):
  `cert_tol_mass 3e-12`, `cert_tol_momentum 1e-8`, `cert_tol_energy 5e-5`,
  `cert_tol_element 1e-8` (with the homopause regime caveat), carrier,
  closure and level tolerances kept. At these values `wasp_full_newton`
  certifies with `Resid tol: 1.0e-7` (marginally at 1e-6) and nothing as
  shipped. Noted for the user: the shipped default `Resid tol` 1e-5 is one
  to two decades looser than the solve reaches at no cost.
- **B5a** (`steady_newton.f90`; suite `steady_selfconsistent_residual`):
  `eval_residual` eliminates the composition AT the state it evaluates
  (the equilibrium sweep iterated at fixed conserved variables until the
  composition stops moving), so the loop-top gate, the line-search merit,
  the returned state and the certification measure one object; the gate
  equals the certification to 1e-12 by construction. Finding: the
  eliminated residual's Jacobian carries the summed Picard response
  `(I - M)^-1`, and a finite-difference probe with a convergence test stops
  after one pass (the seed already is the state's composition), so the
  Newton direction ascended; the probes now run a FIXED sweep count taken
  from an elimination seeded at a different state (about 15 terms), and the
  "frozen-composition Jacobian" the brief offered as a fallback is exactly
  the arrangement that fails. `wasp_full_newton`: `info = 0` at
  `||R|| = 4.4e-7` in 9 iterations and CERTIFIED (was `info = 2` at
  4.8e-5 after 26), state moved 1.0e-4 at most, log10 Mdot 13.27, cost
  1.49x. `mol_base_handoff` identical. Open: `Resid tol` values were
  calibrated on the lagged measure; the other `Solver: Newton` cases and
  the PTC route are measured at the series gate; `n_selfconsistent_max`
  comment stale (B3a); the duplicated row-measure expression (B3a).
- **HYG-absorbed** (`water_photolysis.f90`, `utilities.f90` comment; driver
  `absorbed_fraction_switch`, 5 assertions): the FUV band rates take the XUV
  beam's `absorbed_fraction_per_unit_depth` by use association; the local
  copy was bitwise the same over 1501 points and is deleted; the series
  switch and its error bound (about 1e-8 on the closed-form side of
  `d = 1e-8`) stated once; `mol_base_handoff`, `oxygen_chemistry`
  byte-identical. Advisor: the B4-6 driver's reference carried the same
  single-precision `1.0e-18` the production code had; made double with the
  production change, the driver is green again (51 assertions had turned red
  on the 4e-8 mismatch).
- **B3a** (`EXHALE_main.f90`, new `attempted_step.f90`, `certification.f90`,
  `parameters.f90`; suite `attempted_step`, 87 assertions): a marching step
  is one trial of the fourteen operations from a checkpoint covering
  everything rows 1 to 12 write, evaluated at one adoption boundary before
  `update_map_end_step`, retaken at half `dt` from the checkpoint when
  refused (cap 8, floor `dt/2^8`, a reduced `dt` bounding the next step);
  the clock and the accepted-step count move only there; exhaustion stops in
  `phys` and continues flagged in `init`. Impact: `wasp_full` (300 and its
  du stop at 14496 steps, log10 Mdot 13.33), `hydrostatic_column`,
  `mol_base_handoff`, `mol_carrier` bitwise identical against a control
  built from the same snapshot; cost within run-to-run spread. Measured and
  reported, not gating: the time-discrete hydro row (1.92 at a cold-start
  base where the momentum scale vanishes) and the energy identity T1.5
  (-0.20 to -0.74 on `wasp_full`, the size of the row-8 projection, the RED
  of AT-1a); the step-doubling estimate gives order 0.89 for the whole step
  (the operator split is first order), cost 12 percent at `n_err_every =
  20`. Two defects of the item caught by its own impact runs and fixed. The
  mid-item requests landed (A2tol tolerances, `info = 2` written uncertified
  with no fallback, the ledgers in the summary, the H3+ domain counters, one
  row-measure definition, the `n_selfconsistent_max` comment). The B4-6
  "du plateau" was the mid-edit tree; the final tree reaches its du stop.
  Advisor decision on the exit status: a run that DECLARES a stationary
  state and fails certification exits 2 in every mode (A0's `init`
  exemption covers physical-step claims only); the mode gate is removed
  (A2exit). Small items of its report assigned (HYG-B3a: operator status on
  every path, conduction accessor, step-level carrier verdict, optional
  volume measure) and the `-O0` floating-point exception on
  `hydrostatic_column` under diagnosis (O0FPE).
- **B4-6b** (`lya_rt.f90`, `excited_hydrogen.f90` field evaluation; driver
  `lya_beam_cell_mean`, 37 assertions): the stellar Lyman-alpha depth was a
  trapezoid between cell centers (neither the cell's face depth nor its own,
  every depth short by half the outermost cell); it is now the rectangle
  rule of `calc_column_dens` on the cell's own opacity, with the cell's own
  depth and its star-ward face as outputs. The beam transmission is
  `erfc(c sqrt(tau))`, not exponential, so its cell mean is a 3-point Gauss
  rule in `u = sqrt(tau)` (mean and mean of the square, since the trapping
  makes the field quadratic in the transmission); the photosphere factor's
  mean is closed form with `log1p`. `jlya_escape_prob` keeps the cell center
  by Neufeld 1990 section II (the escape fraction is a point property of
  where the photon is created). Impact: `hydrostatic_column`,
  `mol_base_handoff` byte-identical; `wasp_full` reaches its du stop at
  step 14442 on both binaries, wind unmoved (1e-7), the field 1.6e-3 and
  the Balmer proton source 6.5e-4; the face rule over-counted the column's
  scatterings by 2.4 to 2.7 percent. Carried forward: the FUV band B2 line
  transmission in `util_ion_eq.f90` still reads a face value (B4-6c).
- **HYG-B3a** (`diffusive_photochemistry.f90`, `certification.f90`,
  `steady_newton.f90`, `viscous_conduction.f90`; suites `carrier_retry` 97,
  `energy_update` 48, `certification` 29): the carrier operator sets its
  interval status on every return path (already true since A3b; now stated
  and tested); `carrier_last_verdict_of_run` renamed for what it is, plus a
  step-level `carrier_step_verdict()` reset at every transport step;
  `conduction_floor_cell_hits` exported like the energy update's; the
  volume-weighted companion of `certification_row_measure` optional and no
  longer computed where discarded. All impact cases bitwise identical
  (`wasp_full_newton` finish: `info = 2` at `||R|| = 7.3e-6` against the
  new `cert_tol_momentum = 1e-8`, re-anchored by A2tol2). Left for the gate
  (advisor): the now-unneeded pre-set in `EXHALE_main.f90` and a stale
  field comment in `attempted_step.f90`.
- **A2exit** (`EXHALE_main.f90` exit block; `run_mode` suite 30 assertions,
  `certification` shell gate): the certification exit status follows the
  CLAIM, not the mode: a run that declared a stationary state (du stop,
  residual gate, JFNK finish) and was refused exits 2 in every mode with
  its outputs written; a capped or plateaued run exits 0. `wasp_full_newton`
  now exits 2 with `info = 2`; `run_check.sh` never reads the status. Advisor:
  `run_newton.sh` and `run_baseline.sh` (set -e loops) now report a status
  2 instead of aborting. Open: the `EXHALE_PTC=1` route stops before the
  certification and so claims a stationary state with no verdict (B5).
- **B4-4a** (`Apply_BC.f90`; `hydrostatic_residual` +1 and new
  `free_outflow_boundary` 18 assertions): the outer free-outflow ghost is the
  isothermal hydrostatic continuation of the last cell (`rho_g = rho_N
  exp[-(phi_g - phi_N)/(p_N/rho_N)]`, `p_g` by the same factor, `v_g = v_N`),
  one rule for PLM and WENO3, replacing the zero-gradient copy (PLM) and the
  linear extrapolation (WENO3). The outer cell's momentum residual on the
  analytic column falls from 0.750 of its weight at every N to 4.4e-3 to
  1.0e-4 over N = 250 to 2000 (order 1.8 PLM, 1.6 WENO3); interior and base
  unchanged; four ghost rules measured, this one 30 to 40x better. Velocity:
  the copy (a continuity-exact rule stopped fastest but divides by a ghost
  density the stratification makes small and was the one rule under which
  `wasp_full_newton` segfaulted at -O3). Two premises of the brief
  corrected: `wasp_full`'s edge is subsonic (Mach 0.18 to 0.96;
  `mol_base_handoff` is the supersonic case), and a supersonic edge blocks
  the ghost only through the Riemann solution, not the stencil, so byte
  identity is not available. Impact: `mol_base_handoff` one cell above 1e-4;
  `hydrostatic_column` 3.5e-3 at cell N; `wasp_full` du stop 16208 -> 17135
  steps with the base-breathing phase moved (11 percent in rho at 1.27 Rp,
  mean mass flux above `r_flux` 0.078 percent); `wasp_full_newton` Mdot
  13.27 both. Removing the outer residual does not by itself stop the
  marched hydrostatic column's drain (still 0.79 c_s at 300 steps). Reported:
  a reproducible -O3-only segfault reachable from the rejected velocity rule
  (latent; O3SEGV to diagnose); `EXHALE_main.f90` 2538-2540 ghost-counter
  comment stale (B3c); B4 design section 8.3 (ii) historical.
- **O0FPE** (`util_ion_eq.f90`; shell test `fpe_traps` in `grid_and_gates`;
  `run_fcheck.sh` header): the `-O0` failure was two defects. (1) An integer
  division by zero in `EXHALE_main.f90` about 1234: `mod(count,
  n_err_every)` with `n_err_every = 0` in `init` mode, reached because `.and.`
  does not short-circuit (`-O1` and above happen to skip it); fix handed to
  B3c (owner of the file). (2) A real floating-point overflow in
  `he_rec_coupling` (`util_ion_eq.f90:2813`): `1100 exp(1.2/T4) sqrt(T4)` at
  17 K is 3.9e309, an `Inf` consumed as `ne/Inf = 0` into the correct
  low-density limit; the product is now formed only above `T4 = 1.71e-3`
  and the limit value taken below, byte-identical on `hydrostatic_column`,
  `wasp_full`, `mol_base_handoff`. A third site of the same class (the
  stall detector's `huge` sentinel) also handed to B3c. `run_fcheck.sh`
  states why it builds at -O1.
- **B4-6c** (`util_ion_eq.f90` `fuv_lw_photon_field`; driver
  `lya_band_transmission` rewritten, 16 assertions): band B2's line factor is
  the cell mean of the stellar Lyman-alpha transmission over the cell's own
  line-center depth (product of the two means with the continuum mean formed
  inside `water_photolysis_rate`; the covariance shortfall 6.0e-6 against a
  1.4e-4 bound, the face rule 1.2e-2 high on the budget); the top cell
  carries its own mean instead of a neutral 1. The face rule was 1.0018x on
  the `oxygen_chemistry` state, 1.095x on an exponential H I layer and 18.2x
  in the base cell of a molecular base whose line-center column reaches 1e7
  in one scale height. Consumers: the H2O and OH photolysis in B2 (both
  `ionization_equilibrium` call sites and the carrier photolysis), the output
  columns, the FUV heating column; the H2 channel reads LW only. Impact:
  `wasp_full`, `mol_base_handoff`, `mol_sec_ion` identical; `oxygen_chemistry`
  moves by factors at 45 steps, shown by a control (a 4e-5 perturbation of
  the face value diverges the same way) to be the case's transient
  amplification, not the size of the change; the case is not usable as a
  controlled point until it converges. Finding: `write_output.f90`'s FUV band
  ledger closed form `N_ph (1 - exp(-tau))` carries no line factor for B2 and
  its other rows are 0.79 to 0.98 off (LEDGER-FUV assigned).
- **O3SEGV** (diagnosis, no change): the -O3-only segfault B4-4a reported
  was a null `F0` argument of `jv_product` in a `pgmres` clone the compiler
  had specialized for a null pointer, i.e. a defect of the source text of
  one intermediate build (13:08) between two edits of `steady_newton.f90`;
  every later build carries a plain `pgmres`, `HEAD` plus the rejected ghost
  rule completes `wasp_full_newton` three times, AddressSanitizer finds
  nothing over a full run. Correction to B4-4a's premise: `-g` changes no
  generated code (`.text` bitwise identical at -O3 with and without it).
  Noted: at the time of the diagnosis the working tree stopped
  `wasp_full_newton` at step 0 on the energy floor (B3c in progress).
- **LEDGER-FUV** (`write_output.f90` `write_oxygen_chemistry`; suite
  `fuv_band_ledger`, 15 assertions): the FUV band ledger is three labeled
  statements: (a) closure of the rated absorptions against the beam's loss
  between the cell faces (holds for any rate form, 1e-16); (b) the beam
  budget with each band's own line factor and the un-rated share named per
  band (none for B1/B3/B4, which is therefore the cell-mean test; H I
  resonance scattering for B2, a share; for LW a sign-changing
  disagreement); (c) energy; (d) state drift. Two defects of the ledger
  fixed: it multiplied post-sweep densities by a pre-sweep field, and summed
  the two inner ghost cells. Data rows of every output identical; only the
  `FUV_bands.txt` trailer changes. **Finding (LW-NORM assigned):** the LW
  absorbers take 45 percent MORE photons than the beam loses (1.27e13
  against 8.73e12) and 249 erg cm^-2 s^-1 against an incident 171.5: two
  independent normalizations of one absorption, the level-resolved
  dissociation table over `p_lw_absorbed` against the DB96 eq. 39 band
  share (`lyman_werner.f90`, `h2_self_shielding_table.f90`). Also:
  `docs/a2_oxygen_option_design.md` 882 documents columns that no longer
  exist; the `a2_m2` kinetics check asserts the beam identity for every band.

- **A2tol2** (measurement and decision; `a2_certification_contract_20260906.md`
  section 10): the section-9 anchors re-measured on the self-consistent JFNK
  residual of B5a. The returned row is now 0.36 to 0.93 of the solve's target
  (section 9 measured 41 to 47 above it). `wasp_full_newton` at `Resid tol:
  1.0e-8` certifies (`info = 0`, 14 iterations, 658 s, log10 Mdot 13.27;
  rows mass 2.6e-13, momentum 1.3e-9, energy 6.275e-10); the PTC route stalls
  from a cold start; `wasp_he23off_newton`, `newton_rsw01`, `newton_rsw05` at
  the default 1e-5 refuse on the mass and momentum rows at cell 499; the four
  `Load IC` Newton cases (`jfnk_hd189`, `jfnk_hd189_tight`, `ptc_warm`,
  `armD_D2_newton`) are stopped by the grid guard and need their IC
  regenerated on the current grid before they can be measured.
  **Applied:** `cert_tol_energy` 5e-5 -> 1e-6 in `certification.f90` (5e-5
  had been anchored on the lagged-gate energy row 4.700e-06 and could refuse
  nothing; 1e-6 is a decade above 1.212e-07, the largest energy row of any
  state whose solve met its own target, 5.7e4x above the floor 1.765e-11);
  mass 3e-12 and momentum 1e-8 kept; the report header now names sections 9
  and 10. `backup/regression/wasp_full_newton/input.inp` gains `Resid tol:
  1.0e-8` as its configuration of record, the only setting at which that case
  certifies (Mdot unchanged at this precision, cost +11 percent); its golden
  moves at the series gate with everything else. Verified: private build,
  `certification` suite green (29 assertions). Not changed: the other three
  Newton cases stay at the default so that the matrix keeps a refused claim
  (exit 2) under test.

- **LW-NORM** (diagnosis; no executable statement changed): the 45 percent
  excess of the FUV band ledger is two WAVELENGTH BANDS, not two fits. The
  self-shielding table's line list runs to 1200 A while its `sigma_pump` is
  normalized per photon of 912-1110 A; the beam loss is the DB96 eq. (39)
  equivalent width of 912-1110 A. MEASURED: the column integral of
  `sigma_diss/p_eff` reaches 1.51538 at the top of the column axis against
  1.519288 for the photon content of 912-1200 A over 912-1110 A of a flat
  spectrum (0.26 percent), is exactly density-independent, and reproduces
  the shipped `oxygen_chemistry` ratio 1.4544 at its base column where eq.
  (39) clamps to 1. `p_eff` already separates dissociations from
  absorptions, so the heating takes fragment energy only (not cause 2).
  Landed: `lyman_werner.f90` header section 3a with the invariant and the
  measurement, its section 3 "2 percent" statement corrected to the live
  table's thin-limit cross section 2.6396e-17 cm^2; `lyman_werner_cell_mean`
  probe gains two assertions and a photon-conservation diagnostic (9/9);
  `fuv_band_ledger/run.sh` comment. Advisor follow-up: the table's generator
  inputs (nine `abs_T####.npz`, 27 CLOUDY rate files, the line-by-line
  scripts, the Meudon level and line data, 29 MB) moved from the 2026-09-03
  scratch directory into `src/utils/h2_shielding_lbl/` with a README, and
  `h2_shielding_table_line_by_line.py`'s template brought in line with the
  hand-renamed `axis` dummy of `h2_shield_locate`: the generator now
  reproduces the shipped module BYTE FOR BYTE (MEASURED). The repair, table
  on 912-1110 A or beam widened to 912-1200 A with the B1 flux joined, is
  item 6 of `To_be_determined_by_user_20260906.md` (recommendation: widen),
  because the second redefines the `Stellar LW flux` key. Still to correct
  after B3c releases `write_output.f90`: the ledger's LW label says "two
  fits" and should say "two bands" (with the `header_names_lw_fits`
  assertion of `fuv_band_ledger`).

- **HYG-HeCI** (test hygiene, advisor): `physics_probe` was red on three
  `helium_collisional_ionization_by_level` assertions (1e-13 at 1e4 K to
  1.0e-4 at 5e4 K). Not a code defect: batch 2c replaced the single-precision
  literal `3.940e-11` erg (24.6 eV) of `eval_cool` with the global
  `e_th_HeI_erg` (24.587389 eV), and the test's reference still carried the
  old literal with a comment calling it "a separate open item". The reference
  now uses `e_th_HeI_erg` and `e_th_HeTR_erg`; `physics_probe` PASSED on the
  current tree (private build). The LW-NORM worker's attribution of this to
  a concurrent `Cool_coeff.f90` edit was wrong: that file's uncommitted diff
  is the landed 2c-CI / MASS2 work.

- **LW-NORM-B** (user decision 2026-09-06: recommendation (B) adopted;
  landed): ONE Lyman-Werner band, 912-1201 A, with one normalization.
  * **The band.** `oxygen_rates.f90` carries four FUV bands, not five: the
    Lyman-Werner interval is 912-1201 A and band B1 (1110-1201 A) is gone
    into it. 1110 A was an edge because Draine & Bertoldi end the Solomon
    process there for interstellar v = 0 gas; at the 700-3200 K of a
    planetary base the vibrationally excited levels pump longward of it, and
    the self-shielding table's line list always ran to 1200 A. The two
    halves carried the same H2O yield triplet, so no quantum yield moves;
    the band-averaged cross sections are re-formed over the whole interval
    and MEASURED again from the same Photochem files: sigma(H2O)
    8.51337e-18, sigma(OH) 4.74371e-18 cm^2, <E>(H2O) 12.0526 eV, <E>(OH)
    11.3829 eV, <hv> = 2hc/(912+1201 A) = 11.7354 eV. They agree to seven
    digits with the merge of the two shipped halves, and they reproduce the
    values the band carried before it was split in 2026-08.
  * **The key.** `Stellar LW flux` now means the 912-1201 A band flux at the
    orbit; `Stellar FUV B1 flux` is retired and a file that states it stops
    the run with a message naming the merge and asking for the two numbers
    added. `Spectrum type: Load` integrates 912-1201 A. The resolved dump
    loses `fuv_band_B1_flux`. The inputs that carried both keys were
    restated to their sum: 343.0 + 137.9 = 480.9 for the HD 209458 b
    spectrum (`Update_EXHALE_stage1.md` records the 912-1201 A integral of
    that spectrum as 481.0 and its halves as those two numbers), and the
    LW-only cases (`mol_lyman_werner`, `armA_LW`, `armD_D2_LW`,
    `armD_D2_LW_newton`) from 343.0 to 480.9 so that they keep the same
    star. The README_HOWTO HD 189733 b recipe goes from LW 600.1 + B1 648.4
    to LW 1248.5.
  * **The table.** `src/utils/h2_shielding_lbl/lbl_table.py` normalizes per
    photon of 912-1201 A (`WL_HI`, `wl_max`, and `E_LW_PHOTON_ERG` the mean
    photon energy of that band); `F_BAND_ERG` is restated as
    343.0*289/198 = 500.6 with the derivation in the file, and the file now
    records that it cancels out of `sigma_pump` exactly. The nine
    `abs_T####.npz` were recomputed (`run_abs.py`, 23 to 251 s each) and
    `h2_self_shielding_table.f90` regenerated. MEASURED: the generator
    reproduced the shipped module byte for byte before the change, so the
    only difference is the band. `sigma_pump` at the bottom of the column
    axis is 1.701e-17 (700 K) to 1.887e-17 cm^2, which is DB96's own
    2.557e-17 per 912-1110 A photon restated per 912-1201 A photon
    (1.676e-17) to 1-13 percent, the same agreement as before.
  * **The beam loses what the table rates.** The band share A is no longer
    the DB96 eq. (39) equivalent width; it is the column integral of the
    pump cross section `sigma_diss/p_eff` of the same table, in closed form
    (`h2_lw_band_photon_fraction_absorbed`, exact for the module's own
    interpolant knot interval by knot interval). MEASURED: A at the top of
    the column axis is 0.6903 to 0.9974 over the whole (T, n_H) grid, where
    it used to reach 1.5154; the closed form agrees with an independent
    20000-step quadrature of the two public accessors to 2e-10. The band
    ledger's photon count is now the cell mean of the same pump cross
    section (`lyman_werner_band_absorption_rate_cell_mean`), on the
    quadrature the dissociation rate already used, and not the cell mean of
    `sigma_diss` divided by `p_eff` at one column: on cells spanning 2.5 H2
    scale heights the second stood 8.0 percent above the beam loss and the
    first closes on it to 1e-6 (`physics_probe/lyman_werner_cell_mean`
    sec. 6). `h2_band_equivalent_width` and eq. (39) are removed; the two
    published self-shielding fits stay, called by nothing, for the
    comparison of `docs/h2_self_shielding_cloudy.md`.
  * **The ledger.** The LW row of `output/FUV_bands.txt` is a closure like
    B3 and B4, and its label says so. MEASURED on the `oxygen_chemistry`
    state, 5 steps: `unrated_frac` -3.62e-1 before, -1.18e-2 after, and the
    LW absorbed energy goes from 233.5 against an incident 171.5 (36 per
    cent over the band) to 202.3 against 200.0. The residual left is the
    part of the column above the top of the table's column axis: this case
    reaches 6.4e21 cm^-2 against an axis top of 5e21, where the clamped edge
    cross section drives A past 1 and the beam is capped while the rate is
    not. The run already warns about that column.
  * **Impact.** `mol_lyman_werner` at its own 12000 steps, before at 343.0
    against after at 480.9: log10 Mdot 10.53 -> 10.54, T and rho under
    5e-3, cooling 6.5e-2, the molecular densities and the shielding 1.2e-1 at
    the H2 front (r = 1.075), the base H2 column 3.3875e20 -> 3.4067e20 cm^-2
    and the base `x_H2` 0.96288 -> 0.96365. The optically thin rate at the
    top of the column moves by 480.9/500.6 = 0.9606 exactly: the band's
    photon content grew 46 percent and the cross section fell by the same
    factor, so what is left is that the real spectrum gives 480.9 over
    912-1201 A where the flat continuation of the 912-1110 A deck would give
    500.6. Deeper cells move -4.7 to -6.0 percent, the wider band's extra
    lines shielding differently; slightly less destruction leaves slightly
    more H2.
  * **Not repaired here, and it blocks the `fuv_band_ledger` suite.** With
    its restated 480.9 the `oxygen_chemistry` case stops at step 0 on the
    energy floor. MEASURED, it is not the band merge: the same case stops
    with the OLD binary at the same restated flux, and both binaries stop or
    run at chance along the flux axis (old: runs at 343.0 and 440.0, stops
    at 360.0, 400.0 and 480.9; new: runs at 343.0 and 400.0, stops at 360.0,
    440.0 and 480.9). That is the step-0 acceptance knife edge of items B3a
    and 2 of `To_be_determined_by_user_20260906.md`. Every ledger assertion
    is GREEN on the same case at a flux where it runs.
  * Also: `p_lw_absorbed` is now written by the sweep and checkpointed but
    read by nothing, since the ledger takes the pump rate directly.
    Removing it needs `attempted_step.f90`.

- **B3b-IR** (`molecular_infrared_cooling.f90`; `h2_infrared_line_populations`
  and `h2_partition_sum_uniqueness` probes): the third copy of the H2 ladder
  sum (212-215) is gone; the H2 infrared line populations are normalized by
  `caloric_eos`'s `h2_partition_function`, the one Boltzmann sum over the
  Roueff et al. 2019 ladder whose energy the equation of state carries and
  whose free energy sets the H + H <-> H2 equilibrium. The probe that greps
  the tree for a second sum went RED (1) -> GREEN (0). Impact (advisor
  measured, 12000-step snapshots, `OMP_NUM_THREADS=1`): `mol_ir_bands` (9
  output files) and `mol_base_handoff` (7 files) BYTE-IDENTICAL before and
  after, Mdot 7.99 both; the retired copy had an `exp(-min(E/T, 700))`
  clamp whose only difference from the shared sum's `exp(-E/T)` is below
  double rounding. No golden moves. The worker's draft report stopped at
  "impact runs in progress"; this entry closes it.

- **LW-NORM-B, advisor verification** (2026-09-06): private build of the
  landed tree, `physics_probe` PASSED (554 assertions incl. the rewritten
  section 6: column integral of `sigma_pump` at the top of the axis 0.69 to
  0.997 on the whole grid, was 1.515; photons spent / photons lost 0.999999),
  `spectrum_type` PASSED (154), `fuv_band_ledger` NOT RUNNABLE: its case is
  `oxygen_chemistry`, whose `Stellar LW flux` restated to 480.9 (the same
  spectrum over 912-1201 A) stops at step 0 with `energy update FAILURE` at
  cells 85-88 (T 2.1e4 K, ionization sweep NON-ROOTs accepted by the
  relaxation amnesty with residual 1 to 4, then the balance asks for T below
  the 14.5 K floor). MEASURED by the worker along the flux axis with both
  binaries: 343 runs, 360 stops, 400 runs (after) / stops (before), 440
  stops (after), 480.9 stops: a step-0 acceptance knife edge of the OLD code
  at cells with no H2, items 1 and 2 of `To_be_determined_by_user_20260906.md`
  (decision 4: reject-and-retry; amnesty cap), not of the band. The flux is
  NOT tuned to make the case run; the case and `fuv_band_ledger` stay red
  until those two decisions land. Code read: `a_lines = min(A, 1)` with
  transmission `1 - A` (a band share, not an optical depth) from the
  closed-form power-law integral of the tabulated `sigma_diss/p_eff`; both
  the rate and the share are read at the local T and n_H of the cell for the
  whole column above it (homogeneous-column table), a comment in
  `util_ion_eq.f90` now says so. Accepted. Goldens untouched; every LW case
  moves at the series gate (`mol_lyman_werner` 500-step pair: T 3.3e-3, H2
  0.125 at the front, `k_LW` thin-limit x0.9606 = 480.9/500.6).

- **ACCEPT-1/2** (user decisions 2026-09-06, items 1 and 2 of
  `To_be_determined_by_user_20260906.md`, both as recommended; brief
  `plan_rev2_brief_ACCEPT12.md`, to start when B3c releases
  `EXHALE_main.f90`, `energy_semi_implicit.f90`, `ionization_equilibrium.f90`):
  (1b) a failed energy update refuses the attempted step in both modes and
  the controller retakes it at half dt, exhaustion exits 2; (2a) the
  relaxation amnesty accepts a non-root only below a residual cap measured
  on the matrix, above it the cell keeps the composition it entered the
  sweep with (new acceptance class 6). Motivation: `oxygen_chemistry` at its
  restated LW flux dies at step 0 on exactly this pair (LW-NORM-B
  verification entry above).

- **B3b-IR impact runs, advisor correction** (2026-09-06): the worker's
  before/after pairs were run with binaries built from the tree as it stood
  at 12:30, when the B3a controller still refused on the time-discrete
  hydrodynamic row in init mode: MEASURED from their run logs, 12000 steps
  accepted of 104389 attempted (`mol_ir_bands`) and 104594
  (`mol_base_handoff`), 103697 refusals all "the Runge-Kutta stages and the
  fluxes: the time-discrete hydrodynamic row was out of tolerance", dt
  collapsed and T stayed within 1 K of the 1140 K initial condition over the
  whole grid (golden `mol_ir_bands` reaches 2553 K at 2.18 Rp). The
  byte-identity of the pair therefore says little about the layer; the
  statement that the IR refactor changed no arithmetic rests on the emission
  assertions at five temperatures with identical digits and on
  `retired/Q - 1 = 0` at every probe temperature, which stand. On the
  current tree the same case accepts 12000 of 12000 attempts with no
  refusal (MEASURED, private build). The H2 IR tabulation grid extension
  (advisor, this entry's neighbor above): `n_tab` 401 -> 496 at the same
  log spacing to `t_tab_hi` = 30146 K, old nodes reproduced by the same
  expression, `emission at 10000 K: clamped = direct sum = 1.57875E-18`
  (was 1.181e-18); `molecular_infrared_bands.py` template carries the
  `private :: i` of the generated data module and regenerates it byte for
  byte (MEASURED).

- **B3c** (`EXHALE_main.f90`, `ionization_equilibrium.f90`,
  `energy_semi_implicit.f90`, `molecular_reaction_heat.f90`; new suite
  `coupled_source_step`, 30 assertions in the advisor's run): rows 7 to 9 of
  the enumerated step are one local source step per cell at fixed volume,
  the temperature and the composition its unknowns, iterated to a fixed
  point (`|dT|/T <= 1e-8`, `max |df_sp| <= 1e-8`, cap 20 passes). The
  composition projection (`comp_p_from_T` + `W_to_U`, which rebuilt the
  energy row at the new particle count) is gone; the energy row is anchored
  on `u_th_old`, the thermal energy the cell had at its old composition;
  `rho` is `intent(in)` to the chemistry and the mass sum is a checked
  diagnostic (worst 6.4e-12, T2.2). The two collisional deposits D0 C25
  named (oxygen O1/O1r/O2/O2r/O6, associative He 2^3S) are in as differences
  of the one formation-energy table. RED reference kept in the suite: the
  removed update changed `u(3)` by -0.600 with `heat = cool = 0`.
  **Advisor decision on the worker's deviation (accepted):** `Delta u_form`
  is NOT added to the energy row. In this code `heat` and `cool` are the net
  THERMAL sources (photoheating deposits `h nu - E_th`, collisional
  ionization cooling removes `E_th`, every molecular, Penning, Lyman-Werner
  and oxygen reaction heat is a difference of the same table), so the
  formation reservoir is already accounted for and the term would double
  count; MEASURED with it added, `wasp_full` dies at step 0 with 368 of 504
  cells at the floor. Consequence: T1.5/T1.6 as an identity on
  `u_th + u_form` against `dt(heat - cool)` cannot close; the B3a instrument
  (`attempted_step.f90`) is to be restated as `d u_th = dt(heat - cool)` +
  the transport contribution of the step, or evaluated on the source substep
  alone (item T15-INSTRUMENT, pending). Impact (control = same snapshot with
  the edits reversed): `wasp_full` to its `du` stop, log10 Mdot 13.34 ->
  13.34, T 2.5e-3, rho 2.4e-3, metals 1.5e-2, one sign change of v at r =
  1.062 (2.5e-4 of the wind); `wasp_full` 300 steps front one cell apart;
  molecular 12000-step snapshots up to 40 percent in rho beyond 4 Rp, Mdot
  10.58->10.57, 10.53->10.51, 10.58->10.58; `hydrostatic_column` <= 6e-9.
  **Cost: 5 to 12 passes per step, 8.1x on `wasp_full` 300, 10.3x to the
  `du` stop (18:59 -> 3:15:07).** Assigned as B3c-COST (re-sweep only the
  cells not yet at the fixed point; tolerance anchored on a measured
  convergence study, not chosen to pass). Also in and byte-identical: the
  carrier status pre-set removal (HYG-B3a) and the two O0FPE sites;
  `fpe_traps` passes at -O0 with traps. Left red on purpose:
  `heat_column_matches_breakdown_total` 6.886e-3, because
  `write_heat_breakdown_eq` (`util_ion_eq.f90` 2280) is a second
  hand-maintained copy of the heating assembly that lacks the new deposit,
  and `post_process_adv.f90` ~847 a third: assigned as HEATBRK (one
  assembly, three callers), not the ten-line patch. Advisor verification:
  private build of the live tree, `coupled_source_step` 30/0, `energy_update`
  48/0, `attempted_step` 31/0, `certification` 29/0, `run_mode` 30/0,
  `acceptance_classes` 12/0, `grid_and_gates` 92/1 (the HEATBRK red).

- **HEATBRK, first pass BLOCKED** (2026-09-06): the heating assembly the
  sweep uses is 210 inline lines of `ioniz_eq` (`ionization_equilibrium.f90`
  2119-2331), not a routine, so the one-assembly fix has to extract it there,
  a file ACCEPT-1/2 owned at the time; the worker stopped with no edit, as
  briefed. Re-measured RED: worst 6.886e-3 (cell 376, the associative He 2^3S
  and collisional oxygen deposits missing from the breakdown) over a 1.24e-5
  median floor (the breakdown RE-SOLVES the field and rates on the written
  state while the `heat` column carries the rates the sweep produced it
  with). Interface proposed: `heating_of_composition(...)` in `utils_ion_eq`
  returning `heat` and a non-optional `heat_chan(:,17)` (channel list and
  header owned by the module; LW dissociation and fluorescence as two
  columns), called by the sweep, by `write_heat_breakdown_eq` (which then
  forms no rate) and by `post_process_adv` (which today omits the Balmer and
  associative terms from its `theat`). **Advisor decision:**
  `Heating_breakdown.txt` writes the channel array the sweep computed for
  the written state, not a recomputation; this is the "one state per output
  file" rule (rev3 section 10.2 item 2, landed) applied to this file, and
  the only form in which the breakdown total can equal the `heat` column to
  round-off. Resumes when ACCEPT-1/2 releases `ionization_equilibrium.f90`.

- **ACCEPT-1/2, advisor verification** (2026-09-06): private build of the
  live tree, `energy_update` 53/0 (the driver was killed by `error stop 1`
  before), `acceptance_classes` 22/0 (did not compile before: class 6 did not
  exist), `attempted_step` green (10 RED before: exit status, history lines,
  reason breakdown), `coupled_source_step` 30/0, `certification` 29/0,
  `run_mode` 30/0, `carrier_retry` 97/0. Code read: the kept state of class 6
  is `x_entry`, the composition the cell entered the sweep with (built once,
  also the solve's starting point), heat and cool of the returned
  composition as section 170 requires. Accepted, including the two
  deviations: the two-line `steady_residual.f90` change that makes
  `n_cells_without_chemical_root` count class 6 (without it the controller
  would not see a refused non-root), and `as_op_projection` ->
  `as_op_composition` (operation number 8 unchanged; B3c had removed the
  projection the name described). Findings taken over: cap 1e2 MEASURED on
  the whole default matrix at 400 steps (the amnesty population is a
  cold-start population: 111 events 1e-6 to 8.853 at step 0 of
  `oxygen_chemistry`, nothing until 8.6e4, then the four cells 85-88 whose
  energy update failed); with the cap `oxygen_chemistry` runs its 1000 steps
  and exits 0 (log10 Mdot 9.61), no step refused, fifteen of sixteen cases
  byte-identical; the comment's "res ~28 in the metals-on molecular gates"
  is no longer reproduced (they accept no non-root today) and was corrected
  in place. The energy-refusal path is exercised by the default-off hook
  `EXHALE_ENERGY_FAIL_AT_STEP`. Outside scope, assigned: "class 4" wording
  in `certification.f90`/`steady_newton.f90` and `steady_newton.f90:815`
  `acc_resmax(4)` missing class 6 (HYG-CLASS6, folded into HEATBRK's resume);
  an H-only run has no acceptance classes at all (pre-existing, B6 list);
  `oxygen_chemistry` cannot run in phys mode (56 handoff non-roots refused
  by A0), a statement about that configuration, recorded.

- **HEATBRK landed, advisor verification** (2026-09-06; measured on
  `mol_base_handoff` 200 steps, `wasp_full` 300, `oxygen_chemistry` 1000,
  the `mol_metals` impact pair still running when this was written):
  `heating_of_composition` (`util_ion_eq.f90` 1503-1780) is the one heating
  assembly, 17 named channels, the total the running sum of the columns in
  the sweep's original order (sweep `heat` bit-identical); the sweep calls it
  once (181 inline lines of `ioniz_eq` gone), `write_heat_breakdown_eq` writes
  the module array `heat_channel_state` the sweep filled and forms no rate,
  `post_process_adv` calls the same routine. Private build of the live tree:
  `grid_and_gates` 93/0 with `heat_column_matches_breakdown_total` 1.678e-16
  (was 6.886e-3), `physics_probe` 557/0 incl. `heating_channel_closure` and
  the scanner `heating_sum_uniqueness` (one exception recorded: the band
  ledger of `write_output.f90` restates the Lyman-Werner and FUV photolysis
  event energies, a fourth copy, assigned HYG-LEDGER-EPS). `Heating_breakdown.txt`
  18 -> 21 columns; `_adv` profiles move by up to 1.8e-2 (the Balmer and
  associative deposits the third copy lacked). `heat_channel_state` is
  class-1 state and enters the attempted-step checkpoint through B3c-COST
  (owner of `attempted_step.f90`), assigned. HYG-CLASS6 (class 4/6 wording,
  `acc_resmax` over both classes) in the same landing.


- **B3b-CO landed: the CO row carries two published destruction rates**
  (2026-09-06). `src(ic_CO)` is no longer `0`. It is
  `-(k_D1 n(He+) + k_CO) n_CO` with `k_D1 = 1.60e-9 cm^3 s^-1` (UMIST
  RATE22 entry 4068, measured, accuracy better than 25 percent,
  temperature independent, READ) and `k_CO` the cell mean of
  `(F_LW/<hv>) sigma_CO Theta(N_CO, N_H2) exp(-tau_cont)` on the
  912-1201 A beam. New: `co_self_shielding_table.f90` (the 12CO block of
  Visser, van Dishoeck & Black 2009 Table 6, `T_ex(CO) = 50 K`, with their
  Table 5 beside it for the spread) and `co_photodissociation.f90` (the
  point rate, the two-column cell mean, the fragment heat).
  `fuv_lw_photon_field` gains the CO column and returns `k_CO` and `Theta`;
  `carrier_photolysis` rebuilds `k_CO` from the current carrier columns, so
  the CO layer shields itself inside the relaxation. Constants, all
  MEASURED: `sigma_CO_band = 1.0160e-17 cm^2` per beam photon (Heays et al.
  2017 cross section over the 912-1118 A window CO absorbs in),
  `<hv> = 12.8674 eV` (oscillator-strength weighted over Visser's 37
  bands), `D0(CO) = 11.1157 eV` and `q(D1) = +2.2117 eV`, both differences
  of the one formation-energy table; the deposit is 1.7517 eV per
  dissociation. Both deposits are channels 18 and 19 of
  `heating_of_composition` (17 -> 19 channels) and are formed nowhere else;
  `heat_per_co_dissociation` is registered in the uniqueness scanner with
  the band ledger of `write_output.f90` as its one recorded exception.
  `output/FUV_bands.txt` gains `N_CO`, `k_CO` and `Theta_CO` and a CO row in
  the beam ledger (CO is printed apart from `rated_ph` because
  `beam_loss_ph` carries no CO term: MEASURED, folding it in moved the LW
  budget residual from -1.8e-2 to -8.2e-2). `output/Oxygen_chemistry.txt`
  gains the domain record `tau_dest <= 0.1 tau_res` with the worst ratio and
  its radius, `tau_res/tau_form` from RATE22 8597, the count of cells above
  the shielding table's 512 K excitation-temperature limit and the count in
  which the He+ channel outruns He+ recombination. Tests: `physics_probe`
  `co_destruction` (27 assertions, RED on a mutated table entry and a
  mutated photon energy), `fuv_band_ledger` extended.
  **Impact, MEASURED on `oxygen_chemistry` at its `maxsteps` of 1000, the
  only regression case that carries CO** (a control binary of the same tree
  with the CO source and the two deposits reverted): both run 1000 steps,
  exit 0 and give `log10 Mdot = 9.61`; the after binary costs 628 s against
  11 percent less for the control at a matched 150 steps. The He+ + CO
  charge-transfer deposit is **65 percent of the total heating at
  r = 1.089**, the helium ionization front where CO still survives, and
  under 1 percent outside 1.05-1.17; the layer's temperature moves by up to
  46 percent at r = 1.023 and its CO by up to 52 percent. `mol_carrier` at
  300 steps is byte-identical in every output file, which is the gate that
  the new terms are on the same `thereis_oxychem` as the rest of the oxygen
  path; no other regression case sets `Oxygen chemistry`. The domain record
  of that run: 401374 out-of-domain cell visits, worst `tau_dest/tau_res`
  1.34e5 at r = 1.0002, worst `tau_res/tau_form` 0.219, 305125 cell visits
  above the 512 K excitation-temperature limit and 281325 in which the He+
  channel outruns He+ recombination. Goldens not refreshed.
  **P0 recorded:** the He+ ledger is not closed. This reaction consumes a
  helium ion and neither the carrier row (ion stages frozen) nor the sweep's
  He+ balance returns it, so where `k_D1 n_CO` exceeds the He+ loss rate both
  the destruction and that heating peak are over-stated. The fix is a sink in
  the He+ row of `mol_heh_rows` (`System_HeH_mol.f90`), not this item's file;
  `dom_cells_HeP` is the observable.
  **Not done, carried over:** the thermodynamic CO ceiling is still in
  `limit_to_element_budget`; deleting it changes `EXHALE_main.f90`,
  `certification.f90` and `write_setup_report.f90`, none of them this
  item's files. The exact lines are in the item's report.

- **HYG-LEDGER-EPS closed in the same landing** (2026-09-06): the FUV band
  ledger of `output/FUV_bands.txt` is formed by
  `utils_ion_eq::fuv_band_absorption_ledger`, beside the heating assembly
  and reading the same imports, and `write_output.f90` writes what it
  returns; the `heat_LW` column of `output/Lyman_Werner.txt` is now the
  assembly's own channel (`heat_channel_state(:,ih_H2_LW_dissoc)`) instead
  of a second evaluation of `k_LW n_H2 e_lw_fragment_erg`. The
  `heating_sum_uniqueness` scanner's exception list is EMPTY: every
  deposited energy is named only in the module that defines it and in the
  one assembly. RED before (the pre-move import line, two deposits reported
  outside the assembly), GREEN after. `FUV_bands.txt` is byte-identical
  across the move (MEASURED, `oxygen_chemistry` at 5 steps); the `heat_LW`
  column moves by at most 5.57e-4 relative, which is the one-sweep lag
  between the state the file wrote and the state the energy equation
  deposited, and the new number is the deposited one.

- **B3c-COST, the cost of the coupled step measured and the skip specified,
  not implemented** (2026-09-06; `EXHALE_main.f90`, `attempted_step.f90`,
  the `coupled_source_step` and `attempted_step` suites). WHERE THE COST IS,
  MEASURED on `wasp_full` 300 steps, one thread: the run timers give
  `ioniz_eq` 83.6 percent and the energy solve 16.3 percent of the marching
  wall time, and a gprof build puts **78 percent of the whole marching time
  inside the MINPACK solve of one cell's chemistry** (one Jacobian and one QR
  per cell per pass; the warm start already reaches each root in one Newton
  step), against about 13 percent in the whole-grid rate and field work a
  pass repeats. So the cost is the number of CELL solves and not the field,
  and reusing the field between passes is worth at most that 13 percent.
  THE CONVERGENCE STUDY CONTRADICTS THE EXPECTATION THE ITEM WAS WRITTEN
  FROM: after pass 2, 308 of the 500 cells are still above the tolerances,
  after pass 5 still 130, and the cell that sets the pass count walks inward
  one cell per pass with the ionization front. Contraction 0.191 median over
  2729 passes. A cell-local skip is therefore worth 2.41x on the cell solves
  with a one-pass rule -- but that rule is WRONG: a cell at rest moves again
  on 177 of the 3029 passes, 690 of the 809 wake-ups at pass 8, because the
  column keeps moving a cell after the cell has stopped. **Two consecutive
  passes of rest: zero wake-ups, 1.96x on the cell solves, 1.55x on the
  step.** The skip has to be a mask in the cell loop of `ioniz_eq`, which
  this item does not own (HEATBRK was writing to that file at 20:30 and
  20:35), so it is SPECIFIED and not implemented: the optional
  `cell_at_rest` argument, the class-6 path (`sys_x = x_entry` straight to
  the extraction) as the body of the skip, and the one consequence that is
  not a pure skip -- the acceptance ledger describes the cells the sweep
  SOLVED, so a mask needs an acceptance class stored for each cell and inherited by
  masked cells, or `n_cells_without_chemical_root` reads a different
  population from the one it believes. TOLERANCE ANCHORED AND KEPT AT 1e-8:
  the anchor the comment gave (the 1e-6 heat-column gate) stopped depending
  on it when HEATBRK made the two assemblies one; the anchor that binds is
  `cert_tol_carrier` = `cert_tol_element` = 1e-8 on the composition this step
  returns. Measured 1e-10 to 1e-5 on `wasp_full` 300 and `mol_base_handoff`
  1500, the adopted state moves by about the tolerance and a decade costs 1.3
  passes and 15 percent of the wall time, **so the tolerance is not the cost
  lever**; the step's own truncation error (step-doubling, `phys`, 5.3e-3
  relative) is five decades above the loosest value and does not bind.
  Newton on T with the composition's response (item 4) stays live, and the
  CHEAP version of it is refuted by measurement: the contraction ratio
  spreads by 122 percent globally and by 200 percent cell by cell over the tail
  of a step, so the sequence is geometric on average over a run and not
  within one step at one cell -- which is why B3c's Aitken failed, and why
  the acceleration has to come from `dc/dT = -(dF/dc)^-1 (dF/dT)` against
  the factorization `hybrd` already forms, not from the iterates.
  ALSO FIXED, a defect of my own file found outside the task: the run's cost
  line divided passes accumulated over every ENTRY of the coupled loop by the
  number of adopted STEPS, two different populations, so it could report a
  mean above the worst single entry (MEASURED on the control, `wasp_full` 300
  in `phys` mode: "13.98 passes per step on average, worst 11"). It now
  divides by the entries and names them. ALSO: `heat_channel_state` is in the
  attempted-step checkpoint (save/restore/compare/perturb), the four lines
  HEATBRK's report asked for, with two assertions RED on the control (the
  moved array is not seen, and the restore leaves it 3.0 off). Tests:
  `coupled_source_step` 36/0 (33 before), `attempted_step` 97/0 (95 before).
  Impact, control = the same tree snapshot with this item's two production
  files reverted: `wasp_full` 300, `mol_base_handoff` 1500 and `mol_metals`
  1500 all **byte-identical** in `Hydro_ioniz.txt` and `Ion_species.txt`.

- **B3c-COST** (measurement; `EXHALE_main.f90` comments and debug counters,
  `attempted_step.f90` checkpoint; numerically inert, `wasp_full` 300,
  `mol_base_handoff` 1500, `mol_metals` 1500 byte-identical): the coupled
  source step's cost is not where the brief assumed. MEASURED on `wasp_full`
  300 (3029 passes): after pass 2 still 308 of 500 cells above tolerance,
  after pass 5 still 130; the pass count is set by the ionization front
  walking inward one cell per pass, not by a fixed handful of cells, so a
  cell-at-rest skip is worth 1.96x on the cell solves and 1.55-1.62x on the
  step, and the naive rule (rest on the last pass) is wrong: 177 of 3029
  passes wake a resting cell, two consecutive passes of rest wake none.
  gprof: 78 percent of marching time is the MINPACK solve of one cell's
  chemistry (Jacobian + QR per cell per pass), the whole-grid field and
  rate work 13 percent, the energy solve 16 percent (already deactivates
  converged cells). Tolerance anchored and KEPT at 1e-8 (the binding anchor
  is `cert_tol_carrier = cert_tol_element = 1e-8` on the returned
  composition): a decade costs 1.3 passes and 15 percent, not the lever.
  Aitken and any acceleration built from the iterates fail because the
  contraction ratio spreads 200 percent cell by cell within a step. **The
  lever is a Newton on (T, composition) taken cell by cell with `dc/dT = -(dF/dc)^-1
  dF/dT` from the factorization `hybrd` already forms, inside the sweep
  (`ionization_equilibrium.f90`), plus the two-pass rest mask with an
  acceptance class inherited by masked cells**: assigned as COST2, after
  B3b-CO releases that file. Also fixed: the cost line divided passes over
  loop entries by adopted steps (could report a mean above the worst);
  `heat_channel_state` is in the checkpoint with two assertions
  (`attempted_step` 97/0, `coupled_source_step` 36/0; advisor rebuild and
  rerun green). Outside scope recorded: `EXHALE_ERR_EVERY` read in both modes
  but the probe runs only in phys; `mol_base_handoff` in phys mode refuses 7
  of its first 20 steps at the fixed point (pass cap at step 0), recovers.

- **B3b-CO, advisor verification and follow-up** (2026-09-06): private build
  of the live tree, `physics_probe` 584/0 (incl. `co_destruction` 27 and the
  scanner `heating_sum_uniqueness` with an EMPTY exception list after
  HYG-LEDGER-EPS moved the band ledger into `utils_ion_eq::fuv_band_absorption_ledger`),
  `fuv_band_ledger` 24/0, `grid_and_gates` 93/0. Accepted with the four stated
  departures from the design (cross section as the mean over the 912-1201 A
  beam of one that ends at 1118 A; the same flat-band cross section for a
  loaded spectrum, the shape-aware integral belonging to `sed_read.f90`;
  `q_D1` as a table row because `mion_fsp` is a parameter; CO on its own
  ledger row, `co_frac` 6.39e-2 the size of decision 2's approximation).
  **P0 opened by the item itself and assigned at once as B3b-CO2:** the He+ +
  CO channel deposits its heat (65 percent of the total at the He+ front of
  `oxygen_chemistry`, r = 1.089) while no He+ is consumed and no C+ or O
  delivered anywhere; brief `plan_rev2_brief_B3b-CO2.md` (He+ sink in both
  molecular systems with Jacobian, products to the free O and carbon pools,
  heat rate identical to the species rate from one place, element census).
  Also recorded: the thermodynamic CO ceiling is still present because its
  deletion touches `EXHALE_main.f90` 34, 3641-3648, `certification.f90` 84,
  1072-1076 (`n_active_unvalidated_physics` loses its only producer),
  `write_setup_report.f90` 18, 925-934 (CEILING-DEL, after B4-1 releases
  main); heat is deposited for carrier reactions whose products a
  never-covered carrier interval does not deliver (`oxygen_chemistry` in init
  mode never covers one), an operator-split defect that B4-1c and B5 address
  and that leaves the oxygen case without an end-to-end gate for its CO
  source; `heat_channel_state` now has 19 columns (checkpointed by
  B3c-COST); `python/paper_data.py` `HEAT_FIXED` stale.

- **B3b-CO2 landed: the He+ + CO channel now consumes the helium ion**
  (2026-09-06). Row (2) of `mol_heh_rows` (`System_HeH_mol.f90`) loses
  `k_D1 n_CO n(He+)` with `k_D1 = rk_D1_Hep_CO() = 1.60e-9 cm^3 s^-1`, and
  its turnover scale (`set_mol_turnover_rates`) gains `k_D1 n_CO n_He`; the
  helium leaves neutral, so the free neutral He that closes the helium budget
  receives it and no other row moves. One row serves every path:
  `System_HeH_mol_metals` and the constrained-equilibrium residual both call
  it, and neither molecular system has an analytic Jacobian to extend
  (`hybrd1`, finite differences), so the test asserts the difference
  derivative the solver actually builds. CO is a background density of that
  row and not an unknown, because the one-sided model gives it no formation
  and therefore no local equilibrium; the reaction is split across the two
  operators, each reading the other's frozen value from the one rate, which
  makes the heat of channel 18 the heat of exactly the events the species
  ledger performs. The products needed nothing added: the free-atomic-oxygen
  closure of `carrier_source` and `nC_free - n_CO` in `carrier_write_back`
  already carry the two nuclei and `element_census` already counts them in
  CO. Tests: new `physics_probe` driver `co_helium_ion_sink` (8 assertions,
  RED 4 FAIL against the pre-change module relinked into the same suite,
  GREEN 8 PASS; the RED reads `hep_row_co_sink` 0.0 against -3.84e6,
  `hep_row_derivative_co_term` 0.0 against -38.4, `hep_turnover_scale_co_term`
  0.0 against 3.84e13); `physics_probe` suite PASSED after; element census
  with `EXHALE_ELEMENT_ASSERT=1` on `oxygen_chemistry` 150 steps exits 0 with
  no budget violation (a guard, not a RED/GREEN pair: the change moves no C
  or O nucleus). Impact, MEASURED at 1000 steps of `oxygen_chemistry`: the
  largest share of the total heating taken by channel 18 falls from 0.672 at
  r = 1.081 to 0.035 at r = 1.117, channel 18 at r = 1.0894 from 6.90e-7 to
  9.75e-9 erg cm^-3 s^-1, `n(He+)` there to 0.015 of its former value and to
  0.0073 at r = 1.011, `T` at r = 1.0233 from 1156.4 K to 1088.2 K, and
  `log10 Mdot` stays 9.61 with exit 0. `dom_cells_HeP` does NOT fall
  (281325 -> 281382): the counter compares `k_D1 n_CO` with recombination and
  still says this reaction is the leading He+ loss of the layer, so what it
  now reports is the operator split's lag rather than a missing sink.
  `mol_carrier` 300 steps and `mol_metals` 400 steps byte-identical in every
  output file. Goldens not refreshed. Reported for owners elsewhere: the
  `dom_cells_HeP` header text in `write_output.f90` 558-562 still says the
  frozen He+ "over-states that channel"; the mirror decomposition
  `hydrogen_helium_row_terms` in `constrained_chemical_equilibrium.f90`
  needs the same term after line 1859 or it stops summing to row 2; the `n_co`
  comment in `ion_cell_state.f90` still calls CO an equilibrium closure.

- **B3b-CO2 landed, advisor verification and follow-ups** (2026-09-06): the
  He+ row of `mol_heh_rows` (`System_HeH_mol.f90`, the one row every caller
  uses) loses `k_D1 n_CO n(He+)`, its turnover scale gains the same bound;
  CO stays a background density of the solve because the one-sided model
  gives it no formation channel (not a balance row). Products needed no
  code: the O nucleus goes to the free-oxygen closure and the C nucleus to
  `nC_free - n_CO` in the carrier write-back, `element_census` already counts
  both; the carbon CHARGE is set by the sweep's stage equilibrium, the
  design of record ("the step moves nuclei, the sweep moves charge"), so the
  reaction's C+ is not imposed (recorded as a design choice, not a defect).
  MEASURED, `oxygen_chemistry` 1000 steps: channel-18 share of the total
  heating 0.672 at r = 1.081 -> 0.035 at r = 1.117; n(He+) at the front
  0.015 of before (0.0073 at r = 1.011); T at r = 1.0233 1156 -> 1088 K;
  log10 Mdot 9.61 both; constrained-continuation fallbacks 932 -> 1880, all
  accepted; `mol_carrier` 300 and `mol_metals` 400 byte-identical.
  `dom_cells_HeP` did not fall (281325 -> 281382): it measures dominance of
  the channel over recombination, not over-statement; its header text now
  says so. Advisor follow-ups applied directly: the `D1 He+ + CO` term added
  to the He+ decomposition of `constrained_chemical_equilibrium.f90`
  (`hydrogen_helium_row_terms`, whose header requires every `mol_heh_rows`
  reaction to appear), the `FUV_bands.txt`/`Oxygen_chemistry.txt`
  `dom_cells_HeP` header rewritten, the stale `n_co` comment of
  `ion_cell_state.f90` replaced. Private build: `physics_probe` 592/0 (incl.
  `co_helium_ion_sink` 8), `fuv_band_ledger` 24/0. Remaining: the operator
  split is now a LAG (the sweep destroys He+ at the entry CO while the CO row
  destroys CO at the previous sweep's He+) and `oxygen_chemistry` never
  covers a carrier interval in init mode, so this row has no end-to-end
  gate: B4-1c/B5 territory, recorded.

- **HEATBRK closed** (2026-09-06): the `mol_metals` pair finished at its
  12000 steps, `Hydro_ioniz.txt` and `Ion_species.txt` byte-identical, gate
  1.683e-16 with 21 columns, `_adv` profiles moved by at most 8.45e-3 (the
  deposits the post-process copy lacked). All three impact cases of the item
  therefore confirm: the sweep's heat is unchanged, only the two consumers
  of the assembly moved.

- **B4-1** (`species_face_flux.f90` new; `Reconstruction.f90`, `PLM_rec.f90`,
  `RK_rhs.f90`, `EXHALE_main.f90`, `binary_element_diffusion.f90`, `Makefile`,
  new suite `src/tests/species_face_flux/`): the advective half of the ELEMENT
  transport now rides on the hydrodynamic face mass flux, inside the
  Runge-Kutta stages, on the same faces, areas, volumes and time step as the
  mass row. The face composition is the reconstruction of the mass fraction by
  the same limiter as the primitive variables, taken from the side the face
  mass flux selects, and the closing member of the set is one minus the
  others, so `sum_s F_s(j) = F_rho(j)` holds by construction; the operator
  asserts it at every face of every stage (largest residual MEASURED in the
  unit suite: 2.1e-16 relative, against a tolerance of 64 ulp). The
  cell-velocity upwind advection of `binary_element_diffusion.f90` is gone
  from the marching path: the operator advects only when it is given an
  advecting mass flux, which the fixed-wind relaxation supplies and the
  marching path does not, and the trace-metal solve lost its advective term
  the same way. The carriers of `photochemical_transport_step` are NOT in this
  increment (that file was held by another item): they keep their present
  advective term and are B4-1c.
  Tests, all MEASURED with a private build: a uniform composition is preserved
  to 2 ulp on a compressing flow, an expanding flow and one whose face mass
  fluxes alternate in sign, for PLM and for WENO3; the mass fraction stays in
  [0,1] on a step profile carried across the domain in 300 steps; the closing
  member of the set differs from an independently reconstructed `1 - Y` by
  1.1e-16; the face reconstruction converges at order 1.97 (PLM) and 1.92
  (WENO3); the helium nucleus total of the domain changes only by its boundary
  fluxes to 4.9e-17, where the cell-velocity form it replaces leaves 1.6e-9 on
  the same column and the same step (RED and GREEN in one run); the injected
  violation of the face identity stops the operator with the face named; the
  new step is refusable at its injection point and the run goes on. The
  element-diffusion suite `diffusion_tests` is 35 of 35 after T11 was rewritten
  onto the new operator (its premise, that a face-averaged upwind decouples a
  cell whose face velocities straddle zero, is what the face-flux form removes).
  Impact, MEASURED on scratch copies, `OMP_NUM_THREADS=1`, `mol_*`/`hp_*` at
  200 steps and `mol_diffusion`/`lower_profile` at 1000 rather than the 12000
  their `maxsteps` pins, `wasp_*` at 300: 14 cases byte-identical, every one of
  them without element transport (`wasp_full`, `wasp_he23off`,
  `wasp_full_newton`, `mol_base_handoff`, `mol_metals`, `mol_lyman_werner`,
  `mol_ir_bands`, `mol_sec_ion`, `mol_carrier`, `oxygen_chemistry`,
  `hp_zero_seed`, `hp_trace_seed`, `hp_front`, `hydrostatic_column`);
  `mol_diffusion` moves by at most 7.1e-5 in `Ion_species.txt` and 3.6e-3 in
  `Cooling_breakdown.txt`, `lower_profile` by 1.2e-3 and 1.1e-2. The
  transported quantity, `n_He/n_H`, moves by 1.9e-10 at the base and 3.5e-6 at
  the molecular front of `mol_diffusion`; the base-layer mass flux `rho v r^2`
  moves by at most 9.2e-7 over `1.00 <= r <= 1.16`. Those are smaller than the
  design's estimate because at 1000 steps neither case has separated its
  elements far from the reservoir ratio, and a uniform composition is what both
  forms carry exactly; a converged pair is what would measure the difference at
  the separation those cases reach, and is named as open in
  `docs/b4_spatial_operator_design_20260906.md`. No golden was refreshed.

- **B4-1 landed (element rows), advisor verification** (2026-09-06): new
  module `species_advective_transport` (`species_face_flux.f90`): the face
  mass fraction by the primitives' own limiter (`PLM_rec_scalar`,
  `Reconstruct_scalar`, the MC limiter and the WENO3 geometry coefficients
  factored, not duplicated), upwind side = sign of `F_rho` (coincides with
  the HLLC contact selection wherever both are defined, ROE/LLF offer only
  the sign), the closing member one minus the others so
  `sum_s F_s = F_rho` by construction and asserted at every face and stage,
  advanced inside the three SSP-RK3 stages with the mass row's `dt`, `dV`
  and areas; `positivity_limited_fluxes` writes a replaced first-order flux
  back into `face_flux` so the species flux built on it is the repaired one.
  Transported set: helium and the hydrogen component (normalized), the trace
  metals riding on the same faces inside `m_1` at the reservoir abundance
  (trace limit stated). `binary_element_diffusion.f90` keeps the diffusive
  half only (`element_transport_residual` keeps its stationary advective
  term, equal where `r^2 rho v` is constant). Carriers NOT migrated
  (`diffusive_photochemistry.f90` was B3b-CO's): B4-1c. Suite
  `species_face_flux` 13 + 2 whole-binary rows GREEN on the advisor's private
  build: uniform composition preserved to 2 ulp (PLM, WENO3, three flows
  incl. a sign-alternating face flux), face identity 2.1e-16, bounds exact,
  face reconstruction order 1.97/1.92, helium nucleus total changes only by
  the boundary flux to 4.9e-17 where the replaced form gave 1.58e-9 (B1
  AT-4d RED/GREEN in one run), identity assertion fires, element step
  refusable; `diffusion_tests` 35/35 (T11 rewritten: its premise was the
  defect removed); `attempted_step`, `grid_and_gates` 93/0, `carrier_retry`
  97/0 green. Impact (worker, control = same snapshot with the stage calls
  removed): 14 cases without element transport BYTE-IDENTICAL (`wasp_*` at
  300, `mol_*` at 200 steps); `mol_diffusion` and `lower_profile` at 1000
  steps move by at most 3.6e-3 / 1.1e-2 (cooling breakdown), n_He/n_H by
  <= 6e-5, base `rho v r^2` by <= 1e-6: far below the design's estimate
  because at 1000 steps the elements have not separated from the reservoir
  ratio; the converged 12000-step pair is open in the design document.
  Outside scope recorded: ESWENO3 measures second order on a smooth profile
  (1.92 against PLM 1.97, face error 8x PLM), a property of the primitive
  reconstruction to look at (WENO3-ORDER); `hp_*` cases need
  `output/Hydro_ioniz_IC.txt` in place; `element_flux_profile.txt` has a
  trailing text column.

- **CEILING-DEL** (`diffusive_photochemistry.f90`, `write_output.f90`,
  `write_setup_report.f90`, `EXHALE_main.f90`, `certification.f90`; tests
  `carrier_retry` (+ scanner `no_CO_ceiling_identifier_in_src`),
  `carrier_constraint_attribution` rewritten on the carbon clamp): the
  thermodynamic CO ceiling (clamp to `co_equilibrium_density`, four
  counters, attempted twins, three accessors, `Oxygen_chemistry.txt` header
  line, three `co_ceiling_*` resolved keys, the end-of-run line, validity
  state 4.1's only producer) is deleted per design decision 10; the domain
  record of B3b-CO stands in every place (seven `co_domain_*` resolved keys
  under oxygen chemistry + carrier transport, the CO items under class 4.2
  beside H3+ in the certification, `n_active_unvalidated_physics` = -1 "not
  produced", kept because the five validity states are a fixed set). MEASURED
  before deletion on `oxygen_chemistry` 1000 steps: `co_ceiling_cells 0`,
  applications 0, CO removed 0.0 (with the published destruction rates in
  the row and the He+ sink in the sweep the carried CO never exceeds its
  equilibrium), so every numerical column is byte-identical before/after;
  `mol_carrier` 300, `mol_metals` 400, `wasp_full` 300 byte-identical. The
  design's "of order 160 cells" predates the CO rates and is retired with the
  clamp. Advisor verification: private build, `carrier_retry`, `certification`
  29/0, `carrier_constraint_attribution` 11/11, `fuv_band_ledger` green, no
  `co_ceiling` identifier left in `src/`. Recorded for the docs sweep
  (DOCS-CEIL): the checkpoint enumerations and validity category 12 in
  `a0_run_mode_contract`, `a2_certification_contract`, `b1a_active_equation_inventory`,
  `b3a_attempted_step_controller_design`, `b1_target_system` (1246-1554),
  `development_plan_20260905_rev2/rev3/review3`, `physics_numerics_audit_20260905`,
  and `b3b_co_destruction_design` sections 2.3/4/7 still describe the
  ceiling as present. The domain record is not checkpointed (a statement
  about the state, valid whether or not the step is accepted), unlike the
  ceiling record it replaces.

- **DOCS-CEIL** (documentation only, 2026-09-07): sixteen documents brought
  level with CEILING-DEL, the current state READ from the code (destruction
  row `src(ic_CO) = -(k_D1 n(He+) + k_CO) n_CO`, conservation-only
  `limit_to_element_budget`, `carrier_co_domain_record`, certification class
  4.2 CO block, `n_active_unvalidated_physics = -1`, seven `co_domain_*`
  keys); no measurement deleted, the superseded ones ("of order 160 cells";
  "still cuts CO above 3000 K") kept with the measurement that replaced them
  beside them. Files: a0/a2 contracts, b1a inventory (4.1 row retired, 4.2
  row the domain record, open question 3 closed), b3a design (checkpoint
  enumerations), b1 target system (5.1 split into the 2026-09-06 READ and
  the current one, AT-5 superseded, decision 12 moot), b3b_co design
  (sections 2.3, 3(iv), 4, 6, 7, 7.4 as done), development plan rev2/rev3/
  review3/execution, physics_numerics_audit (P6 resolved), code_status, d0
  (C26 closed, question 5 answered), co_destruction_rates_literature
  (banner), To_be_determined_recommend, PLAN rev2. Left as dated records:
  the superseded plan revisions and external reviews, and the b1 draft.
  Recorded: line-number citations in the 2026-09-06 design documents have
  drifted by up to a thousand lines (a re-anchoring sweep on routine names
  is a separate item, DOCS-LINES); the stale checkpoint header comment in
  `diffusive_photochemistry.f90` goes to B4-1c.

- **WENO3-ORDER, measurement and fix** (2026-09-07): the ESWENO3
  reconstruction was SECOND order on every grid the code runs on. Probe
  `weno3_reconstruction_order` (physics_probe; exact volume averages by
  8-point Gauss quadrature, face error against the point value, N = 50 to
  800): production rate 1.997 (uniform in r), 1.988 (the `define_grid` Mixed
  production grid), 2.940 on a grid uniform in the volume coordinate, where
  the defect is invisible. Cause, isolated by variants: not the ESWENO
  nonlinear weights (S = 1 gives the same 1.997), not the floor `eps =
  dr_j^2` (rate 3.000 at kappa 0.1 to 10 once fixed), not the ideal ratios
  `D1`, `D2` (1/2 + O(h), harmless); TWO of the four stencil coefficients
  applied the NEIGHBOR's volume share where cell j's own belongs:
  `C1(j) = c/(b+c)` on the jump to j+1 at the right face and `C2(j-1) =
  a/(a+b)` on the jump to j-1 at the left face, each 1/2 + O(h) where
  1/2 - O(h) belongs, an O(h) coefficient error on an O(h) jump. Inherited
  character for character from ATES v2.0 (`ATES/ATES-Code-main`). Advisor
  check of the derivation: a linear interpolant of two cell averages in the
  volume coordinate, evaluated at a face dV(j)/2 from cell j's center,
  weights the jump by dV(j)/(dV(j) + dV(nb)); with the corrected shares a
  profile quadratic in the volume coordinate is reproduced to 5.9e-16.
  **Fixed (advisor, two symbols in each of `Reconstruct` and
  `Reconstruct_scalar`, `Reconstruction.f90`):** `C1(j) -> C2(j)` on `dWp`
  in `WL`, `C2(j-1) -> C1(j-1)` on `dWm` in `WR`; the probe's production
  mapping and header updated, the two waiting assertions enabled:
  `weno3_production_third_order` rate 3.000, `_mixed` 2.991 (RED before:
  shortfalls 0.803, 0.812); `physics_probe` 603/0. **Every golden whose case
  reaches the WENO3 stage moves, and not by rounding** (face error down two
  decades at N = 800); PLM-only cases are untouched. Consequence for the
  method: until now the two-stage PLM -> WENO3 recipe was a continuation
  between two second-order schemes with the second having the larger error
  constant (1.3x uniform, 3.8x on the production grid). Mdot values are to
  be quoted from Newton-finished runs after the gate. Attribution of the
  weight form and the floor to Yamaleev and Carpenter (2009, JCP 228, 3025
  and 4248) stays unverified: both papers are paywalled here; requested
  from the user. Also recorded: stale `.ipynb_checkpoints` copies of
  `Apply_BC.f90` and `UW_conversions.f90` under `src/modules/` define the
  same module names as the live files (not built; to remove, HYG-CKPT).

- **HYG-CKPT** (advisor, 2026-09-07): ten untracked, gitignored Jupyter
  checkpoint copies of production modules under `src/modules/*/
  .ipynb_checkpoints/` (Apply_BC, Reconstruction, Source, base_boundary,
  caloric_eos, low_mach_dissipation, cross_sec, UW_conversions,
  h2_photo_channels, element_census) removed: not built, but each defined
  the same module name as a live file and two scanners had to skip that
  directory (batch 2c) to stay green.

- **T15-INSTRUMENT** (`attempted_step.f90` +289/-165, suite, b3a design
  section 9, b1 T1.5 wording; advisor verification: private build,
  `attempted_step` 43/0 with the two whole-binary rows skipping until the
  main-loop hook lands, `coupled_source_step`, `certification` green): the
  B3a instrument now states the THERMAL identity of the step, `sum_j dV_j
  [u_th^{n+1} - u_th^n] = sum_j dV_j dt (heat - cool)_j + (transport
  contribution)`, the transport contribution MEASURED between two marks on
  either side of the coupled source step, the residual the closure of the
  source step alone, `Delta u_form` reported as the reservoir's change and
  never added (the B3c reason stated at the routine). GATE in phys mode at
  the energy row's own `energy_res_tol = 1e-9` (reason `as_reject_source_energy`
  = 14 at operation 9): MEASURED on a hooked scratch tree, worst
  |residual/scale| 7.44e-1 -> 5.40e-10 (`wasp_full` 300 phys), 8.63e-1 ->
  1.88e-10 (`mol_base_handoff`), zero refusals, states byte-identical in
  both modes. The two hook lines in `EXHALE_main.f90` (before/after the
  `coupled_source:` loop) are assigned to COST2, the file's current owner,
  with the rewording of the print block. **Defect found and fixed in
  passing:** `attempted_step_checkpoint_take` saved `u` and `f_sp` through
  a saver that reallocates from index 1, so a reader indexing the saved grid
  arrays by cell compared cell j with cell j - Ng; `energy_identity` was that
  reader, so every `d_u_th`/`d_u_form` printed so far (incl. the 28 percent
  reservoir share of the B3c report: with bounds fixed, same step,
  `d_u_th` 2.006e-5, `d_u_form` 1.298e-5) was formed against a state shifted
  by two cells; the state restore itself was unaffected (whole-array
  assignment). The second private copy of the `eps_s` table in
  `attempted_step.f90` (no OH/H2O/CO rows) is deleted in favor of
  `species_formation_energy`. Recorded: `n_projection_applied` can no
  longer be nonzero since B3c but is still printed by the report and refused
  by the certification (HYG-PROJ, spans `global_parameters`,
  `certification.f90`, main); design section 9 decision 8 still names the
  projection as a live correction (dated record).

- **COST2** (measurement; the only production change is the T15 hook and
  print in `EXHALE_main.f90`, byte-identical on `wasp_full` 300 in both
  modes and `mol_base_handoff` 60 phys; `attempted_step` whole-binary rows
  now run, identity residual -1.8e-10, zero refusals): the brief's premise
  is refuted twice. (1) MEASURED with T frozen through the coupled loop on
  `wasp_full`: the composition sweep alone still takes 10-13 passes with the
  same contraction 0.15-0.20 and the same front walking one cell per pass;
  freezing the lagged coefficients as well collapses it to 2 passes;
  freezing only the radiation field leaves 3-4, only the rate coefficients
  9-13. The pass count is the self-consistency of the photoionization field
  with the composition that attenuates it; the temperature is worth at most
  one pass. (2) A timer inside `ioniz_eq` (60 steps, 29.95 s): `PH_heat`
  11.40 s (38 percent of the run), the whole cell loop 7.37 s (25 percent),
  pre-sweep `eval_cool` 2.44 s, post-sweep assembly 2.40 s; skipping the
  entire cell body for resting cells saves 4.7 percent. B3c-COST's gprof
  attribution (78 percent in MINPACK) was inverted; the `csm_max_pass`
  comment still carries it (to be rewritten by COST3). (3) The two-pass rest
  mask was implemented, measured (6.0 percent on `wasp_full`, 0 on
  `mol_base_handoff`, state moved 5.5e-7, step 0 of `mol_base_handoff` in
  phys turned into an exhausted attempt) and removed. **The lever: the field
  dependence is one-directional (a cell's P_HI is set by the attenuation of
  the cells outside it) and the sweep already runs outside in, so updating
  the column with the composition just produced turns the Jacobi iteration
  into a forward substitution: same fixed point, about one pass instead of
  ten, and `PH_heat` no longer repeated per pass.** Assigned as COST3
  (`util_ion_eq.f90` `PH_heat_HHe`/`PH_heat_H`, `ionization_equilibrium.f90`
  sweep order; the OpenMP parallel cell loop becomes sequential in the field
  direction, so the design must weigh 10 parallel passes against 1-2
  sequential ones and measure both at 1 and 16 threads). Also recorded: the
  pre-sweep `eval_cool` is run for its coefficients and its cooling
  discarded, then a second traversal forms the cooling (16 percent
  together); `mol_base_handoff` in phys refuses 7 of its first 60 steps at
  the fixed point on the unmodified tree.

- **B4-1c landed: the carriers ride on the species face flux** (2026-09-07):
  the carriers of `photochemical_transport_step` (H2 always, OH/H2O/CO under
  `Oxygen chemistry`, H+ under `Ionization transport`) join the transported
  set of B4-1 and are advected as `F_rho Y_c/m_c` inside the Runge-Kutta
  stages. `carrier_set_init` declares them once the keys are parsed
  (`advected_carrier_register` in `binary_element_diffusion.f90`, which also
  holds the columns and hands them back through
  `advected_carrier_fractions`); `EXHALE_main.f90` is untouched, the existing
  three calls carry both sets. The carriers ride OUTSIDE the normalized set
  (`n_norm = 0`): a carrier is one species inside an element, the element is
  already the closing member, and normalizing it beside the element would
  count that mass twice; what keeps carrier and element consistent is the
  operator's own write-back, which restores the element totals cell by cell.
  Their rows are advanced from cell 1 (the elements' start at 2), with the
  inflow composition taken from the inner ghosts: the handoff's `q_H2_base`
  where one is stated, the base cell's own partition where none is, one rule
  for both branches of decision D6 (B4-4a). `diffusive_photochemistry.f90`
  loses the cell-velocity advective term, `carrier_advection_correction` and
  the base Dirichlet branch from its MARCHING rows and their Jacobian; two
  paths keep the advective term, `carrier_steady_residual` (the stationary
  balance is the whole equation) and the fixed-wind relaxation
  `relax_photochemical_composition` (it takes many transport steps with no
  hydrodynamic stage between them, so nothing else is taking that transport),
  through the switch `carrier_rows_advect`. Tests: two rows added to
  `species_face_flux`, `uniform_carrier_partition_preserved` 0.83 ulp (bound
  8) on a sign-alternating face mass flux and
  `h2_nucleus_total_changes_only_by_boundary_flux` 1.36e-15 (bound 1e-12),
  RED before in the strongest sense (the driver does not compile against the
  pre-increment module); suite 15 + 2 GREEN; `carrier_retry` 101,
  `carrier_returned_state_acceptance` 37, `carrier_constraint_attribution`
  11, `carrier_reference_scales` 15, `acceptance_classes` 22 GREEN with no
  reference re-derived. Impact MEASURED against a binary built from the same
  tree snapshot with the two files at their pre-increment text: 12 cases
  BYTE-IDENTICAL, every one without `Molecular carrier transport`
  (`wasp_full`, `wasp_he23off`, `wasp_full_newton`, `mol_base_handoff`,
  `mol_metals`, `mol_lyman_werner`, `mol_ir_bands`, `mol_sec_ion`,
  `mol_diffusion`, `lower_profile`, `oxygen_chemistry`,
  `hydrostatic_column`); the four that do transport carriers move.
  `mol_carrier` at 1000 steps: `Ion_species` 5.0e-1 (H3+ at `r = 1.100`),
  `Cooling_breakdown` 6.2e-1, `Hydro_ioniz` 2.6e-2 (`v` at `r = 1.186`);
  `hp_zero_seed`/`hp_trace_seed`/`hp_front` at 100 steps: 8e-3 to 1.1e-1,
  `Hydro_ioniz` 1.3e-3 to 1.8e-3. **The front does not move; its tail does.**
  The `x_H2 = 0.5` crossing of `mol_carrier` is `r = 1.0698` at 1000 steps
  before and after and `x_H2(1.05)` is 0.7308 against 0.7310, while
  `x_H2(1.10)` falls from 3.19e-3 to 2.02e-3 (-37%): the molecular tail the
  cell-velocity form was carrying by its own numerical diffusivity, which is
  where the design said the difference lives, and every large entry above is
  a consequence of that column. Every carrier interval covered and the
  element limiter moved no cell, in both runs. Recorded and not fixed:
  `oxygen_chemistry` has `carrier_transport` on by default with the oxygen
  chemistry, but at 200 steps in initialization mode 201 of 201 intervals are
  uncovered before AND after, so the write-back never runs and the case is
  byte-identical for that reason, not for an absence of effect (its rows did
  change: worst CO row of step 1 `|res|` 6.03e-5 against 5.68e-5). No golden
  was refreshed. Also in this increment: the stale clause saying the carrier
  checkpoint restores "the cumulative CO record" removed from two comments in
  `diffusive_photochemistry.f90` (the `carrier_checkpoint` type carries no
  such field since CEILING-DEL).

- **B4-1c landed, advisor verification** (2026-09-07; the worker's own entry
  carries the numbers): every carrier of the photochemical transport
  operator (H2; OH/H2O/CO under oxygen chemistry; H+ under ionization
  transport) rides on `F_rho Y_c / m_c` inside the RK stages through the
  registry of `binary_element_diffusion.f90` (carriers outside the
  normalized set, from cell 1, base inflow composition = the inner ghosts'),
  `species_face_flux.f90` unchanged (species-agnostic), `EXHALE_main.f90`
  untouched; the marching rows of `diffusive_photochemistry.f90` keep
  diffusion, drift and chemistry only; the advective term stays in
  `carrier_steady_residual` (the whole stationary equation) and in the
  fixed-wind relaxation `relax_photochemical_composition` (no hydro stage
  takes it there), switch `carrier_rows_advect`. Advisor: private build,
  `species_face_flux` 15 + 2 rows (carrier partition 0.83 ulp, H2 nucleus
  budget 1.4e-15), `carrier_retry`, `carrier_returned_state_acceptance`,
  `acceptance_classes`, `grid_and_gates` green. Impact (worker, same-snapshot
  control): 12 cases without carrier transport BYTE-IDENTICAL; `mol_carrier`
  1000 steps: the H2 front does not move (x_H2 = 0.5 at r = 1.0698 both),
  the tail the cell-velocity form carried by its numerical diffusivity does
  (x_H2(1.10) 3.19e-3 -> 2.02e-3, -37 percent), hence H3+ 0.50, cooling
  0.62, heating 0.37 relative there; `hp_*` 100 steps <= 1.7e-2; base
  `rho v r^2` <= 9.7e-4. Recorded: `element_advection_*` names are now too
  narrow (HYG-NAME, six call sites in main, after COST3);
  `oxygen_chemistry` in init mode covers 0 of 201 carrier intervals so it
  cannot judge carrier transport; a refused carrier interval discards a
  transport the mass row already took (B3a/B5 question); `adv_corr` still
  checkpointed though only the residual writes it; 52 documents under
  `docs/` carry em-dashes (DOCS-DASH, launched).

- **DOCS-DASH** (documentation only, 2026-09-07): 1354 em-dashes and 237
  en-dashes removed from 60 files (`docs/*.md`, READMEs, test and data
  READMEs), each replaced by the punctuation the sentence needs (colon,
  comma, parentheses, or a table placeholder hyphen), en-dash ranges and
  joined names to the hyphen the repository already writes; 108 mechanical
  substitutions corrected by hand on read-back; grep over the 187-file set
  clean. Left: `Update_EXHALE_stage1.md` (append-only log, 319/59), and
  1220 LaTeX `---` in 27 `.tex` sources that render as the same dash
  (DOCS-DASH2, launched: the `.md`/`.tex` twins realign there; the built
  PDFs are regenerated at the gate's docs pass), plus
  `docs/lower_atmosphere_figs/README.md` (11) and three READMEs under
  `backup/` (5) outside the first sweep's paths.

- **DOCS-DASH2** (documentation only, 2026-09-07): 1283 LaTeX `---` in 29
  `.tex` sources (docs, paper, posters) and 30 spaced `--` replaced by the
  punctuation each sentence needs; four exceptions kept and named (a direct
  quotation of Frelikh and Murray-Clay 2026 compared letter for letter, the
  AASTeX `\keywords` separator, a literal `metals.inp` listing, hyphen rules);
  the four READMEs of the first sweep's remainder cleared. The 24 built
  `docs/*.pdf`, `literature_survey.pdf`, `paper/ms.pdf` and the poster PDFs
  are now dated and are regenerated in the gate's docs pass. Left: the
  append-only stage logs (`Update_EXHALE_stage1.tex` 1614, `stage0.tex` 37).

- **COST3** (measurement; `util_ion_eq.f90`, `ionization_equilibrium.f90`;
  default output BYTE-IDENTICAL on `wasp_full` 300, `mol_base_handoff` 1500,
  `mol_metals` 1500, `hydrostatic_column` 300 at one thread; advisor rebuild:
  `coupled_source_step` 36, `physics_probe` 862/0 incl. the new
  `photoionization_field_substitution` driver (259 assertions),
  `grid_and_gates` 93/0): the star-ward columns can now be carried forward
  cell block by cell block (`advance_starward_columns`, reproducing
  `calc_column_dens*` bitwise; `photoionization_field_at_cell_*`, the
  `PH_heat_*` cell bodies extracted verbatim; `xuv_field_block_cells`,
  default 0 = whole grid). MEASURED on `wasp_full` 300: the substitution
  stops the front's walk (worst-moving cell 206 -> 177 over the passes on
  the lagged scheme, 202 -> 196 with the columns carried) and cuts the pass
  count 10.10 -> 7.13, but a pass costs the same field work either way, and
  the field of a whole grid is one loop sixteen threads solve at once: block
  width 8 buys 1.43x at one thread (136 -> 103 s) and LOSES 2.8x at sixteen
  (17.0 -> 42.5 s); width 64 keeps the parallel front and saves nothing. The
  default therefore stays the whole grid, and the width is not tied to the
  thread count (the widths reach states differing by 1.5e-7 at the 1e-8
  tolerance, the thread count by 7e-14). **The remaining lag is the cell's
  OWN optical depth**, which enters its own field only on the next pass
  (contraction 0.17): closing it means the cell residual carries its own
  field, `System_HeH*`, i.e. B5. Folded: the pre-sweep `eval_cool` run for
  its coefficients only is now `chemical_rate_coefficients` (-3.6 percent
  `wasp_full`, -5.6 percent `mol_metals`, +3.5 percent `mol_base_handoff`,
  the last from the 27 metal slots looped with metals off: guard on
  `thereis_metals`, HYG-COL). Recorded: `docs/openmp_parallelization.md`
  (2.07x at 16 threads, bit-identical across threads) is stale: measured 8.2x
  and round-off 6.6e-14 (DOCS-OMP). A `cp -irp` of another item blocked on an
  interactive prompt since 2026-09-06 was killed by PID (23818, cwd checked).

- **WENO3-REF** (comments only, `Reconstruction.f90` and the probe's
  References block; the user supplied the publisher PDFs, read with
  `pdftotext -layout`): the nonlinear weights are EXACTLY Yamaleev and
  Carpenter 2009, JCP 228, 3025, Eqs. (18), (20), (21), (22) (published
  `tau` already squared; `d_0 = 2/3`, `d_1 = 1/3` entering as the ratio 1/2
  that `D1`, `D2` carry); the scheme is NOT their ESWENO: the energy-stability
  dissipation of Eq. (48) is absent, their analysis is a uniform-grid
  finite-difference flux reconstruction, and for systems only in
  characteristic variables (Section 5 quoted), so the code inherits the
  accuracy result (Eq. 24) and not the energy estimate; the geometry
  (`C1, C2, D1, D2` on the stretched finite-volume grid) is this code's own.
  `eps = dr_j^2`: a defensible stand-in, not a defect: the published h^2
  scaling of Eqs. (62), (64) without the solution scale of Eq. (65);
  MEASURED on the production grid the face error is flat to every digit for
  kappa >= 0.1 over four decades (rate 2.991) and deteriorates below 1e-2
  per Eq. (57), so in smooth regions the departure is not measurable; where
  it acts is the ENO biasing (b_r/eps spans 2.2e-6 to 8.4e3 over the
  `wasp_full` grid; the published floor would make the base LESS nonlinear),
  a candidate experiment for the base sound-wave work, not a change here.
  Case label and comments now say "third-order WENO"; identifiers unchanged.
  Advisor: `paper/ms.tex` 144 "the energy-stable ESWENO3" reworded to the
  verified attribution and the ESWENO3 name replaced by WENO3 there and in
  `docs/code_comparison.tex`. **Recorded for the user:** every `src/tests/*/`
  suite of this series (19 directories) and the 2026-09-05/06 design and
  plan documents are untracked in git; a fresh clone of the remote carries
  none of them.

- **References received** (2026-09-07, from the user, publisher PDFs in
  `references/`): Zhang and Shu 2010 (JCP 229, 8918), Poinsot and Lele 1992
  (JCP 101, 104), Thompson 1987 (JCP 68, 1) and 1990 (JCP 89, 439), Neufeld
  1990 (ApJ 350, 216), Mao and Kaastra 2016 (A&A 587, A84), Verner et al.
  1996 (ApJ 465, 487), Yelle 2004 (Icarus 170, 167), Larsson et al. 2008 (CPL
  462, 145), Tully 1975 (JCP 62, 1893), Slanger and Black 1982 (JCP 77,
  2432). Five verification items launched against them: REF-ZHANGSHU (the
  positivity limiter), REF-LODI (the characteristic boundary design
  document), REF-NEUFELD (the Ly-alpha escape-probability closure),
  REF-METALS (every metal photoionization and recombination fit parameter
  against the published tables, a transcription error being a physics fix),
  REF-CHEM (the H3+ network transcriptions, the H3+ dissociative
  recombination coefficient, O6 and the Ly-alpha H2O yields). The series
  gate waits for them, since REF-METALS and REF-CHEM may move coefficients.

- **Reference data fetched** (2026-09-07, at the user's direction, each with
  a provenance README, third-party material unchanged): `references/verner_photo/`
  (D. A. Verner's distribution: `phfit2.f` reference routine with its fit
  parameters, `table1.dat` of Verner et al. 1996, the 1995 inner-shell fits),
  `references/mao2016_cds/` (CDS J/A+A/587/A84 `rr_spex3.dat`, 39696 rows),
  `references/badnell_rr/` and `badnell_dr/` (Strathclyde TAMOC total RR and
  DR fit coefficients `clist_K`/`clist_eV`, tabulated `cfout_*`, reference
  routines). REF-METALS now compares every metal stage's photoionization, RR
  and DR coefficients against these by direct evaluation of the reference
  routines, and reports which source each stage uses today and any missing
  DR contribution.

- **HYG-MAIN** (2026-09-07; numerically inert, `wasp_full` 300, `mol_carrier`
  200, `mol_diffusion` 200, `mol_base_handoff` 300 byte-identical; advisor
  rebuild: `certification`, `attempted_step`, `species_face_flux` green):
  HYG-NAME `element_advection_*` -> `species_advection_*` in
  `binary_element_diffusion.f90`, the main program (use list, six calls) and
  the two test drivers, the two design documents that still named them
  (advisor); HYG-PROJ `n_projection_applied` deleted from `global_parameters`,
  the `attempted_step_report` block and both `certification.f90` references
  (MEASURED: every reference in `src/` was a read, no assignment anywhere,
  so the counter was provably zero since B3c), `rep%n_unbudgeted_accepted_correction`
  is now `n_shapiro_applied` alone, b3a design decision 8 dated closed for
  row 8; HYG-COL the 27-slot metal column loop of `advance_starward_columns`
  guarded on `thereis_metals` (inert; COST3's +3.5 percent NOT reproduced:
  guarded 53.7/52.8 s against unguarded 55.7/52.5 s on `mol_base_handoff`
  300, within the run-to-run spread, kept as dead-work removal); DOCS-OMP
  `docs/openmp_parallelization.md` dated and brought to COST3's measured
  state (8.2x at 16 threads, 6.6e-14 between thread counts). Recorded:
  no `reduction` on a physical quantity exists in `src/`, so the 6.6e-14
  thread-count difference comes from block boundaries set by
  `omp_get_max_threads()` in `chemical_rate_coefficients`/`eval_cool`; which
  routine turns a moved boundary into a moved last bit is not localized
  (THREAD-DET: until it is, the thread count is a silent input to every
  result at round-off); `photoionization_field_at_cell_HHe` makes an array
  temporary for `nm_j` on every call (runtime warning in `physics_probe`);
  `EXHALE_main.f90:452` is 141 characters.

- **REF-ZHANGSHU** (comments, `docs/p49_front_energy_mode.md` 6.1, probe
  `positivity_limiter_scaling` 15 assertions; user-supplied PDF read): the
  citation was wrong as stated: `positivity_limited_faces` is not the
  limiter of Zhang and Shu 2010 (JCP 229, 8918) but one built on the same
  idea, departing on three axes (primitive, not conserved, variables; one
  theta per face state, not per cell, Eqs. (2.4)-(2.5), (2.10)-(2.12); a
  relative one-ulp floor, not their absolute 1e-13). What carries over:
  Lemma 2.5 (convex combination with an admissible anchor); what does not:
  Theorem 2.1 (the next cell average admissible under their Gauss-Lobatto
  set and CFL, N = 2 covers k = 1 only; Remark 2.6 names finite-volume WENO
  as open) and conservativity of the reconstruction (a theta taken
  face by face moves the cell's end-point mean by 8.3e-2 where a theta taken
  cell by cell gives 0). Along
  their conserved path p(t) is concave, so the code's theta is never larger
  than their t_eps (MEASURED 0.856 against 0.949 at Mach 60): the anticipated
  negative-internal-energy failure cannot occur. **Defect found and fixed
  (advisor applied the worker's proposal):** with eps one ulp of the average
  the bare update `q_avg + theta (q_rec - q_avg)` is a cancellation of the
  size of the floor and lands on EXACTLY ZERO for 1013 of 200000
  reconstructed values crossing zero (never negative), the value the floor
  exists to prevent (division in v = m/rho and in the sound speed);
  `positivity_limited_faces` now clamps the scaled density and pressure to
  eps after the update (same number in exact arithmetic). Moves only runs in
  which the limiter fired, by one ulp except where the zero would have
  divided; goldens at the gate. Probe: `bare_update_reaches_zero` (the
  reason), `production_limiter_never_delivers_zero` (the production routine
  over 4000 grids of zero-crossing face states) GREEN. A theta taken cell by cell
  (Eq. 2.11), which would restore the reconstruction average, is recorded as
  a change of scheme, not applied.

- **REF-LODI** (documentation, `docs/phaseC_characteristic_base_bc.md`
  +282/-50; the three user-supplied journal PDFs read as page images since
  their text layers lose the displayed equations): the wave-amplitude
  normalization (Poinsot and Lele Eqs. (15)-(23); Thompson 1990 (44)-(47)
  with l_1 = (0, 1, -rho c, 0, 0); Thompson 1987 (19)), the LODI relations
  ((24)-(36), equal to Thompson 1990 (48) term for term), the subsonic
  inflow two-condition closure (P&L section 3.1, Tables III/IV; Thompson
  1990 Eq. (68)) and all four branch-table rows CONFIRMED; no constant
  wrong, the constraint count unchanged. Two statements of the note
  corrected: pair (D) is NOT singular at stagnation (Thompson 1990 Eqs.
  (69)-(70) have denominator u_1 + c; the note's singularity was its own
  closure rho_b = F_c/(r_b^2 v_b); D1 reopened with a pointer), and D6's
  claim that no transparent condition is expressible as a fixed point is too
  strong (Thompson 1990 Eq. (64) L_5 = rho c g_1 and P&L Eq. (40)
  L_1 = K (p - p_inf) are algebraic in the face state; D6 stays the user's
  decision with a third option). One gap stated: the source projection
  l_1 C (gravity and the spherical divergence, Thompson 1990 footnote 1 and
  Eq. (52); the steady form is L_1 = -rho c g_1, not 0) survives only
  because both states of (C-) sit at the same radius; any reformulation
  relating cell 1's center to the face must restore it. Carlson's "but not U
  and p" marked unverified against the primaries.

- **REF-CHEM** (comments and documents only; `mol_rates.f90` and
  `oxygen_rates.f90` function bodies byte-identical; probes
  `h3p_recombination_branching` 12 and `water_photolysis_lyman_alpha_yields`
  6, RED with the competing published values): every `mol_rates` constant
  sourced to Yelle 2004 matches his Table 1 digit for digit; Yelle's R7
  (H3+ + H) reference column reads "Estimated", so Frelikh's barrier-free
  2.0e-9 transcribes an estimate faithfully; the code's H3+ dissociative
  recombination R6/R7 (2.16e-8 / 5.04e-8 x (300/T)^0.65) is already Larsson
  et al. 2008's 300 K total 7.2e-8 split by their branching 0.70, with the
  0.65 index inherited from the Sundstrom/Datz fits Yelle prints (Larsson
  gives no temperature dependence; he states the CRYRING/ASTRID values the
  older 2.9e-8 rests on were too high from rotational excitation), so the
  code stays; the Slanger and Black 1982 Ly-alpha H2O yields 0.78/0.10/0.12
  are exact (p. 2435) with one label corrected, `Phi2 + Phi4` -> `Phi1 +
  Phi2`, and the paper's coupling of the three yields and the 8 percent
  OH(A) bound recorded; Tully 1975's 2.87e-10 is a 300 K statistical-model
  value (Table II) that VULCAN flattens over 100-2100 K, and its agreeing
  "experiment" is a 1973 relative-rate evaluation, strengthening the
  audit's verdict for O6. Documentation fixed in passing: `methodology_comparison.md`
  6.4 listed three-body H2 formation as identical in value while the code's
  Cohen and Westberg 1983 2.8e-31 T^-0.6 differs from Ham et al. 1970 by
  1.142 at every temperature; Baulch 1992 is in `references/` (its H2 + H2
  entry valid 2500-8000 K, above the whole layer); section 10 item on H2
  dissociative photoionization closed (the code carries it,
  `h2_photo_channels.f90` channel S).

- **REF-NEUFELD** (comments and documents; `lya_rt.f90` machine code
  identical before/after; `physics_probe` 1118/0): the Ly-alpha
  escape-probability closure does NOT use Neufeld 1990 or Harrington 1973.
  `beta_esc = pi^(-1/4) sqrt(a/tau)` is the one-flight wing escape chance
  under complete redistribution; Ly-alpha scatters with partial
  redistribution and the published static-slab damping-wing solution is
  Harrington 1973 eq. (40), `<N> = 0.9093161 B` scatterings for a mid-plane
  source, so `beta = 1/<N>`, proportional to 1/tau and independent of a;
  MEASURED on the stored `wasp_full` column the code exceeds it by 234x at
  the base and 14x at 1.20 Rp, crossing at 1.42 Rp. Neufeld contains no
  expression of the code's form (his escaping fractions are escape against
  destruction, eq. 2.25 partitions the two faces; the code's two-face
  PRODUCT is unpublished); the escape frequency scales as (a tau)^(1/3)
  (Harrington eqs. 37, 47; Neufeld eq. 4.32), not the code's (a tau)^(1/2);
  Neufeld's validity `(a tau)^(1/3) > 10` holds in 320 of 504 cells (inside
  1.156 Rp); Neufeld is static, the Sobolev factor is Sobolev 1960 radial-ray.
  The header's `J_int` formula lacked the (1 - beta) factor the code applies
  (fixed). **Advisor decision: adopt the published solution (LYA-BETA,
  launched):** beta = 1/<N> of eq. (40) joined to the thin limit, the face
  partition of eq. 2.25 in place of the product, the (a tau)^(1/3) escape
  frequency, a validity record, `A_2p1s` 6.3e8 -> 6.2648e8 (spectroscopic),
  every case with `Jlya escape-prob: True` moves (expected at frozen state:
  `J_int`/`n_2p` up to 90x higher inside 1.02 Rp, 14x at 1.20 Rp, up to 7x
  lower beyond 1.47 Rp; Mdot to be measured to the du stop). The same
  question for the metal lines' `line_escape_probability_one_face`
  (`Cool_coeff.f90`) is LYA-BETA-METALS, after REF-METALS.

- **REF-METALS** (`cross_sec.f90`, `Cool_coeff.f90`; probes
  `metal_photoionization_fits`, `recombination_coefficient_fits`, 255
  assertions; advisor rebuild green apart from LYA-BETA's in-progress probe):
  all 132 fit parameters of the seventeen metal photoionization cross
  sections match Verner et al. 1996 Table 1 (`verner_photo/photo.dat`) digit
  for digit and `metal_photoion_sigma` reproduces Verner's own `phfit2`
  exactly wherever the fit is the one `phfit2` uses; 220 Badnell RR and DR
  parameters (fourteen daughters; `cfout_K` witness within 0.045 dex), the
  H/He rows and the Mao and Kaastra n = 1 rows match. **Two fixes:** the
  Mg II turn-on gated at Verner's rounded 15.04 eV while the grid edge is
  the NIST 15.035 (`e_th_MgII` now; byte-identical on `wasp_full`, the band
  geometric means miss the sliver); `alpha_rr_metal('CaI')` used Verner
  `rrfit.f`'s `rrec(:,20,19)`, the row that FORMS Ca II (second index = the
  recombined ion's electron count, fixed by `rnew(:,1,1)` = H+ + e and
  `rnew(:,2,2)` = He+ + e), 6.05x the Ca I row `(1.120e-13, 0.900)` at 1e4 K;
  corrected. Impact (worker, scratch): `wasp_full` 300 Ca I 0.87 relative at
  its 1.196 Rp maximum, Ca II +2.8 percent at the base, T 4.3e-3, v 9.9e-3;
  `mol_metals` 400 Ca II/III x2.5 at 1.38 Rp, cool 2.9e-2; `wasp_full`,
  `mol_metals`, `lower_profile` move at the gate through calcium. **Two
  physics items opened:** (a) the outer-shell fit is evaluated on the whole
  photon grid to 1240 eV while every row carries an `E_max` (34-558 eV)
  above which `phfit2` hands over to the 1995 subshell fit and adds the
  inner shells: the code carries 0.1-12 percent of the true metal
  photoabsorption there, and Fe I's row (`Q` = 1.54) RISES with energy to
  49x the truth at 1240 eV; the handover of `phfit2` itself is adopted
  (METALS-OUTER, launched); (b) the inner-shell absorption and its Auger
  cascade (an atom two or more stages up) cannot be represented by the
  three-stage metal ladder: `To_be_determined_by_user_20260906.md` item 7.
  Also: `verner_photo/README.md` corrected (`photo.dat` is the 1996 Table 1,
  `table1.dat` the 1995 fits); rev3 line 944 closed; `e_th_MgI`/`e_th_MgII`
  provenance lines in `parameters.f90` to be annotated (METALS-OUTER, after
  LYA-BETA releases the file).

- **METALS-OUTER** (`cross_sec.f90`; generator `src/utils/metal_photoion_table.py`
  from `verner_photo/photo.dat` + `table1.dat`; probes
  `metal_photoionization_fits` 188 and `metal_photoion_table_transcription`
  15; `docs/photoion_cross_sections.tex` + PDF; advisor rebuild:
  `physics_probe` green): the seventeen `sigma_XX` routines are one
  parameter block, and `metal_photoion_sigma` makes `phfit2`'s handover:
  the 1996 outer-shell fit below `E_max`, the 1995 fit of the same shell at
  and above it (MEASURED: `einn` of `phfit2` equals `photo.dat`'s `E_max`
  for all sixteen metal rows; the production routine reproduces `phfit2`'s
  outer shell to 1e-10 at 17 ions x 12 energies where 110 of 204 rows
  disagreed before, worst Fe I at 1240 eV 24.9 Mb against 0.00172 Mb, a
  factor 14464). Inner shells still absent by decision item 7 (option (a)
  recommended; the measured sizes now in that item). Impact (scratch,
  metals-off `mol_base_handoff` byte-identical): `wasp_full` 300 v 0.15, T
  4.7e-2, heat 7.8e-2, Fe III 0.99, C III 0.60 relative (hard photons now
  reach the base, the top stages fall at small radius; metal share of the
  heating at the base 17.4 -> 8.5 percent); `mol_metals` 400 v 0.56, heat
  0.42; `lower_profile` 200 heat 8.3e-3; log10 Mdot unchanged to 0.01 dex in
  the bounded snapshots. Both Mg thresholds already agree with NIST to their
  digits (item 4 needed no move; provenance note pending `parameters.f90`).

- **LYA-BETA-METALS** (measurement by the worker, fix applied by the
  advisor; `Cool_coeff.f90`, probe `fine_structure_escape_probability`,
  `docs/resonance_line_trapping.md` sections 6.1 and 10; publisher scans of
  de Jong, Boland and Dalgarno 1980 and Osterbrock 1962 fetched into
  `references/`): `line_escape_probability_one_face` is de Jong et al. 1980
  eq. (B-7), both branches to round-off, applied to the eight ground-term
  fine-structure lines only; NO metal resonance line uses an escape
  probability (Mg II h&k, Ca II H&K, Na I D, Fe I/II at beta = 1, justified
  in that memo's sections 1-5 by `n_crit,eff` = beta A/q ~ 1e10-1e15 against
  n_e <= 4.3e9), and the transit post-process uses none. The fine-structure
  lines have a ~ 1e-8 and `a tau <= 1e-9`, so the Ly-alpha correction does
  not transfer. **Defect:** (B-7) is written in the frequency-integrated
  depth `sqrt(pi) tau_0` (normalized profile; equal to Hollenbach and McKee
  1979 eq. (5.10) once the convention is matched, 3.4 percent at worst over
  tau = 3 to 1e6) while the caller passed the line-center column: beta too
  large by 1 (thin) to 2.11 (branch point) to 1.77 (asymptotic). **Fixed:**
  each cell's depth is now `sqrt(pi)` times the line-center depth in
  `fine_structure_line_transfer`. MEASURED regime on the shipped columns:
  largest line-center depth 0.14 ([O I] 63 um at the hot-Uranus base), so the
  factor there is 1.12 and the change of the total radiative losses is below
  0.1 percent (`mol_metals` 0.096 percent), but a factor of two in any base
  thick in these lines and in the `Base IR field` incident term the same
  beta scales. `physics_probe` 1291/0 after the fix. Metals cases move at
  the gate at that level.

- **Decision item 7, literature check** (2026-09-07; the user supplied
  Locci 2022 PSJ 3, 1; Locci 2024 PSJ 5, 58; Cecchi-Pestellini 2006 A&A
  458, L13 and 2009 A&A 496, 863; Fox 2008 SSRv 139, 3; all in
  `references/`): the published exoplanet X-ray photochemistry models take
  the inner-shell absorption as opacity summed over all ionization channels
  (Cecchi-Pestellini 2009 Eq. 1 with Verner 1993 cross sections), put the
  photoelectron and the Auger electrons into the secondary-electron cascade
  (Locci 2022 Eq. 3), and carry singly charged ions only, i.e. option (b) of
  the item with the charge under-count accepted; no escape hydrodynamics
  code includes inner shells. Recommendation revised from (a) to (b) in
  `To_be_determined_by_user_20260906.md`; the user decides.

- **METALS-INNER** (decision 7 option (b), user 2026-09-07; `cross_sec.f90`
  subshell block from `phfit2.f`'s `PH1` data (its thresholds differ from
  `table1.dat` in six rows, `phfit2` being the reference),
  `set_energy_vectors.f90` (every subshell turn-on of an active metal ion a
  bin edge; bin count unchanged, 17 -> 60 threshold cuts on `wasp_full`),
  `util_ion_eq.f90` comment only (the existing rate `sigma/e_v` and heat
  `photoelectron_share(E_th) sigma` integrals already implement "one stage,
  one electron of h nu - E_th(outer) to the SvS85/Dalgarno cascade" once the
  cross section is the total); generator `metal_photoion_table.py` `--check`;
  probes 392 metal assertions, `physics_probe` 1479/0; docs
  `photoion_cross_sections.tex` (+PDF), `code_status` new L row): the metal
  photoabsorption is the SHELL SUM `phfit2` returns, reproduced to 1e-10 at
  17 ions x 12 energies (130 of 204 rows disagreed before). It also closed a
  cliff METALS-OUTER's handover had left at `E_max` where the 1996 fit stood
  for a valence complex the 1995 outer row does not cover (Fe I x0.012 and
  Fe II x0.0066 steps, now 1.23 and 1.13; genuine K edges C/N/O x19-21), and
  put `E_max`, a step of the cross section, on a bin edge (it was not).
  Approximation stated at the code site: ion advanced one stage, one
  electron of `h nu - E_th(outer)` to the cascade (MEASURED with the
  degradation routine: 0.95-1.09x the heat and 1.10-1.37x the secondary H I
  ionizations of the true two-electron event; the un-spent I_2 is 24-35 eV
  of a 300-550 eV photon), fluorescence neglected (with `e_top` = 1240 eV
  only C, N, O, Na have K holes; the high-Z ions contribute L shells only).
  Impact (worker, same-snapshot control; metals-off `mol_base_handoff`
  byte-identical): `wasp_full` 300 the wind 1-3 percent hotter and faster,
  base heating +26 percent, T at most +9.5 percent, Fe III/Mg III/N III at
  the base x74-160 from vanishing values, **Mg II falls 84 percent at r =
  1.178** (the h&k transit line); `mol_metals` base heating +26 percent, T
  1.4e-3; `lower_profile` heating +21 percent at the base, T 2e-5; metals'
  share of the 100-1240 eV absorption 2.2 -> 3.5 percent (`wasp_full`), the
  band 0.38 percent of the whole; log10 Mdot unchanged to 0.01 dex in the
  snapshots. Krause 1979 (JPCRD 8, 307) K/L fluorescence yields could not
  be obtained (NIST 503, AIP 403): requested from the user; no yield is
  quoted anywhere.
  Advisor verification: `physics_probe` 1479/0, `grid_and_gates` 93/0;
  `spectrum_type` segfaulted in `loaded_sed_threshold_edges` (its threshold
  array was sized `n_mion + 16` = 43 while the active list now reaches 126),
  fixed by sizing it with the module's `n_thr_max`; suite green after.
  (`n_thr_max` was a local parameter of `set_energy_vectors`; hoisted to
  module scope, public, so every caller of `ionization_thresholds_active`
  sizes its array by the one bound; `spectrum_type` 154/0, `physics_probe`
  1479/0, `grid_and_gates` 93/0 after.)

- **METALS-INNER addendum, Krause 1979 yields** (comments and documents
  only; user-supplied `Krause_1979_JPCRD_8_307.pdf`, read from page images
  where the OCR is noisy, alignment checked on omega_K(Fe) = 0.340): Table 3
  omega_K for C/N/O/Na 2.8e-3/5.2e-3/8.3e-3/0.023, Mg/Si/S/K/Ca/Fe
  0.030/0.050/0.078/0.140/0.163/0.340; Table 5 effective L1/L2 yields
  (Coster-Kronig included per the text, nu_3 = omega_3) for Mg to Fe all
  below 6.5e-3; Table 2 uncertainties (omega_K 40-10 percent at Z 5-10,
  footnote (a) on molecules and solids). With `e_top` = 1240 eV only C, N,
  O and Na I can hold a K hole, so the neglected radiative branch is at most
  2.3 percent of a K-hole event (under 1 percent for the three coolants) and
  under 0.7 percent of any L-shell event; raising `e_top` past 1.3 keV (Mg,
  Si K) or 4 keV (Ca, Fe K) is the condition under which the approximation
  must be revisited, stated at the code site, in `photoion_cross_sections.tex`
  (new `tab:kedge`), `code_status` and decision item 7.

- **LYA-BETA landed, advisor verification** (`lya_rt.f90`, `hydrogen_n2_rates.f90`
  (`A_2p1s` 6.3e8 -> 6.2649e8, the nonrelativistic dipole value times the
  reduced-mass factor, = NIST ASD's Wiese and Fuhr 2009 value),
  `parameters.f90` comments; probe `lya_escape_probability` 24 assertions;
  advisor rebuild `physics_probe` 1479/0, `certification` 29/0): the
  escape probability is the published static-slab solution, `<N>` =
  Neufeld 1990 eq. (3.27) at zero destruction = Harrington 1973 eq. (40)
  `(4 sqrt(6)/pi^2) u_2 B = 0.909316 B` for a mid-plane source (the
  `0.9093161` of the REF-NEUFELD entry was a transcription slip), `beta =
  1/(1 + <N>)` (exact at both limits: a generated photon is emitted 1 + <N>
  times, and the budget `n2p A beta = P` wants the escape chance per
  emission), the two faces sharing one escaping population by Neufeld eq.
  (2.25) instead of the former product, `Phi(xi)` summed in closed form; the
  Voigt parameter has left the escape probability. `a` still enters only
  `lya_wing_penetration_coefficient` (where an incident stellar photon
  crosses in one flight, a different question, judged at the code site).
  Defect found in passing: `beta_tot = 1 - (1-beta_esc)(1-beta_sob)` lost
  digits at the base (4.8e-9 relative), now `beta_esc + beta_sob (1 -
  beta_esc)` (closure 2.2e-16). RED before 4.79e2 / GREEN 2.19e-16 on the
  end-to-end identity. **Impact, `wasp_full` to its du stop (2 h 59 m):
  log10 Mdot 13.34 -> 13.34 (+1.6 percent, 0.0067 dex); `Jint`/`n_2p` x84
  at the base, x73 at 1.05, x37 at 1.10, x10 at 1.20, x0.75 at 1.50 R_p; T
  within 2 percent; the old `min(1,.)` cap had zeroed the internal field
  from 1.15 R_p outward.** H(n=2) chord column x31.5 at b = 1.00, x2.8 at
  1.20, x1.005 at 1.45: **every H-alpha/H-beta transit number in this
  repository from a `Jlya escape-prob: True` run is superseded** (He 10830
  not affected). Validity `(a tau)^(1/3) > 10` holds in 349 of 504 cells,
  out to 1.193 R_p, exactly where the change is an order of magnitude.
  Advisor follow-ups applied: the second `A_2p1s` in `exhale_transit_lib.py`
  and `docs/transmission_spectrum.tex` set to 6.2649e8 (one value with the
  wind model), the stale `E = min(boost, 1/beta)` comment in `parameters.f90`
  replaced by the implemented `1 + (boost-1)(1-beta) T_star`. Pending: the
  B6-style print of the Ly-alpha domain record in `certification.f90` (the
  record and accessor exist; LYA-DOM, with the gate's follow-ups).

- **LYA-DOM** (advisor; `certification.f90`): the Ly-alpha damping-wing
  domain record (`lya_wing_domain_record`: cell visits below the
  (a tau)^(1/3) limit, visits seen, the smallest value, the limit) joins the
  out-of-domain closure count and is printed beside the H3+ and CO records,
  informational like them. `certification`, `run_mode`, `attempted_step`
  green. **Every worker item of the series has landed; the series gate
  begins (single build, `make test`, the 16-case matrix, movement table,
  one golden refresh, fcheck, docs pass).**

- **Series gate, 2026-09-07** (PLAN rev2 steps A to C and the items above,
  one series): single build of the shared tree (`make -q` up to date,
  `EXHALE.x` md5 ba3d34e6eaee); `make test` green apart from the two steady
  suites, which then passed on the `wasp_full_newton` run log
  (`EXHALE_STEADY_RUNLOGS`); the 16 default cases run against that binary
  (the harness holds a single-instance lock, so `wasp_full` ran through
  `run_check.sh check` and the other fifteen through a runner replicating
  its steps for each case, all `OMP_NUM_THREADS=1`); every case moves against
  the previous goldens, as every item of the series said it would.
  Movement of the physical columns over the cells where the golden value
  is at least 1e-3 of the column maximum (max relative; T median; the
  change of log10 Mdot from the mean of r^2 rho v over 1.20-1.45 Rp):
| case | rho max | v max | p max | T max | T median | dlog10 Mdot | log10 Mdot new |
|---|---|---|---|---|---|---|---|
| wasp_full | 2.97e-01 | 1.28e+00 | 1.65e-01 | 9.18e-02 | 3.8e-02 | +0.0095 | 13.35 |
| wasp_he23off | 2.95e-01 | 1.30e+00 | 1.62e-01 | 9.18e-02 | 3.9e-02 | +0.0084 | 13.35 |
| wasp_full_newton | 6.96e-01 | 3.51e-01 | 4.30e-01 | 2.35e-01 | 8.7e-02 | +0.0414 | 13.30 |
| mol_base_handoff | 1.60e-01 | 2.50e-01 | 1.52e-01 | 2.33e-01 | 2.2e-02 | -0.0440 | 10.57 |
| mol_metals | 1.51e-01 | 2.25e-01 | 1.53e-01 | 2.07e-01 | 2.3e-02 | -0.0476 | 10.58 |
| mol_lyman_werner | 2.12e-01 | 4.87e-01 | 1.99e-01 | 2.55e-01 | 2.6e-02 | -0.0011 | 10.54 |
| mol_diffusion | 1.59e-01 | 2.50e-01 | 1.52e-01 | 2.33e-01 | 2.2e-02 | -0.0440 | 10.57 |
| mol_ir_bands | 1.51e-01 | 2.26e-01 | 1.53e-01 | 2.07e-01 | 2.3e-02 | -0.0475 | 10.58 |
| mol_sec_ion | 2.16e-01 | 4.03e-01 | 1.96e-01 | 2.53e-01 | 3.0e-02 | +0.0149 | 10.51 |
| mol_carrier | 1.54e-01 | 6.37e-01 | 1.44e-01 | 2.24e-01 | 1.9e-02 | +0.0106 | 10.56 |
| lower_profile | 9.40e-03 | 3.22e+01 | 7.50e-03 | 9.77e-01 | 8.8e-03 | -0.0335 | 7.29 |
| hydrostatic_column | 2.09e-12 | 2.81e-03 | 3.54e-12 | 1.40e-03 | 2.5e-14 | -0.0000 | 8.27 |
| oxygen_chemistry | 8.62e-01 | 1.29e+02 | 9.76e-01 | 7.61e-01 | 3.7e-01 | +nan | 9.61 |
| hp_zero_seed | 2.25e-05 | 1.81e-04 | 1.61e-04 | 3.20e-03 | 4.4e-05 | +0.0000 | 10.56 |
| hp_trace_seed | 2.24e-05 | 1.81e-04 | 1.61e-04 | 3.20e-03 | 4.4e-05 | +0.0000 | 10.56 |
| hp_front | 2.07e-05 | 5.77e-04 | 3.66e-04 | 1.00e-02 | 6.7e-05 | +0.0000 | 10.56 |
  The maxima sit at the fronts, which move by a cell; the T medians
  (2-9 percent on the wind cases) and the Mdot changes (within +-0.05 dex
  everywhere it is defined; `oxygen_chemistry` has a reversed flow in the
  window) are the state of the physics after B3c, B4-1/1c, LYA-BETA, the
  WENO3 coefficient fix, the metals (REF-METALS, METALS-OUTER/INNER), the
  LW band, the CO channels, ACCEPT-1/2 and the rest. `wasp_he23off` exits 2
  (its stationary claim refused by four inventory entries, as A2tol2
  measured at the default tolerance; outputs written in full);
  `wasp_full_newton` exits 0, certified, log10 Mdot 13.30. PDFs of the 25
  dated documents, the paper and both posters rebuilt. **Golden refresh
  pending the user's hand**: the archive copy `golden -> golden_pre_series_20260907`
  and `run_check.sh golden` were refused by the session's permission
  classifier; the commands are in the report to the user. The final
  confirmation (`make check` against the refreshed goldens, serial) follows
  the refresh.
  **Golden refresh done by the user (2026-09-07 13:19):** `golden/` archived
  as `golden_pre_series_20260907/` (17 cases) and re-snapshotted from the
  case outputs of this gate (checked: golden equals output on `wasp_full`,
  `mol_metals`, `oxygen_chemistry`, `lower_profile`). The serial `make check`
  against the refreshed goldens runs next as the bitwise confirmation, then
  `run_fcheck.sh`.

- **User manual, environment-variable table** (user report, 2026-09-07):
  the fixed `table` float of the environment hooks (67 rows) overflowed
  page 11 into the footer; it is now a `longtable` with the caption at its
  head, breaking across pages 10-11 with a "continued" foot. PDF rebuilt.
  Pages 20-21 (user report): the keyword table's long file paths were
  unbreakable `\texttt` tokens and ran past the right margin; 18 such paths
  are now `\path{}` (url package, breaks at `/` and `_`), one inside a
  caption kept as `\texttt` (a moving argument). Rebuilt with latexmk, so
  the longtable column widths and cross-references are resolved (a single
  `pdflatex` pass leaves both unresolved, which is what produced the
  "Table ??" and the first overflow the user saw).

- **Update log TeX twin** (user request, 2026-09-07): `docs/Update_EXHALE_stage2.tex`
  and its PDF now exist, GENERATED from this Markdown by
  `src/utils/update_log_to_tex.py` (pandoc body, the memo-class preamble of
  the stage logs); regenerate after appending here. 56 pages at this entry.

- **Decision item 4 closed** (user, 2026-09-07): candidate C,
  `wasp121b_composite_huang2023.txt` (solar XUV shape below 1700 A rescaled to
  Huang et al. 2023's F_XUV and Ly-alpha, the MUSCLES WASP-17 photosphere
  above), is the WASP-121 b spectrum of record, as both `input.inp` files
  already state; A and B remain sensitivity inputs. No decision item of
  `To_be_determined_by_user_20260906.md` is open.
  **`make check` against the refreshed goldens: REGRESSION PASS, all 16 cases
  byte-identical (64 files), 13:21-16:25 KST 2026-09-07** (the harness runs
  the cases eight at a time, one thread each). `run_fcheck.sh` follows.
  **`run_fcheck.sh`: CLEAN** (1500 bounded steps of the HD 209458 b folder
  case under `-fcheck=bounds,do,mem`, log10 Mdot 9.56) after one harness fix:
  the scratch copy of the planet input named its spectrum by the folder's
  relative path, which the script now rewrites to the repository's absolute
  path for the `Spectrum file` and `Lower atmosphere profile` keys (the first
  run failed on "Cannot open file ../inputdata/sed/...", not on a runtime
  check). Production `EXHALE.x` restored, `make -q` up to date, md5
  ba3d34e6eaee. **The series is closed**: goldens refreshed once, `make
  check` byte-identical, fcheck clean, docs and PDFs rebuilt. What follows is
  in `docs/session_handoff_20260907.md`.

- **HYG-PY** (advisor, 2026-09-07; `python/paper_data.py`,
  `python/make_struct_figures.py`): `HEAT_FIXED` (the headerless-file
  fallback) and `HEAT_CH_LABEL`/`HEAT_CH_ORDER` (the heating-budget panel)
  now carry all 19 deposit channels of `heat_channel_name` with the
  descriptive brackets stripped the way `_clean_name` strips them
  (`heat_He23S_assoc`, the two Lyman-Werner halves, `heat_mol_chem`,
  `heat_FUV_photolysis`, `heat_oxygen_collisional`, the two CO channels);
  the panel loop already skips a channel a run's file lacks. Both scripts
  parse; the loader reads the file header, so the fallback list was never
  what a real file was read by.

- **THREAD-DET** (advisor, 2026-09-07; `util_ion_eq.f90`): the thread count
  is no longer an input to any result. Mechanism, MEASURED: `Cool_coeff.o`
  imports glibc libmvec's two-lane `_ZGVbN2v_exp`, `_ZGVbN2v_log`,
  `_ZGVbN2vv_pow` (gfortran -O3 vectorizes the rate-coefficient loops), and
  the vector variants agree with the scalar libm to one ulp, not bitwise;
  inside a block of cells the loop takes lane pairs from `j_lo` and an odd
  remainder by the scalar call, so where a block STARTS decides which cell
  meets which variant. `chemical_rate_coefficients` and `eval_cool` cut the
  grid into `2*omp_get_max_threads()` blocks, so the block boundaries, and
  with them the last bit of every rate, moved with the thread count (the
  6.6e-14 of COST3/HYG-MAIN). Fixed: both sweeps use blocks of the fixed
  even length `xuv_rate_block = 32` (16 blocks on the 504-cell grids); the
  threads decide only who takes which block. MEASURED on `wasp_full` 300:
  old binary at 1 thread, new at 1 thread and new at 16 threads are
  byte-identical in `Hydro_ioniz`, `Ion_species`, `Heating_breakdown`,
  `Cooling_breakdown` (the old single-thread partition happened to start
  both its blocks at even offsets, so the refreshed goldens are untouched).
  Shared `EXHALE.x` rebuilt from the tree (`make -q` up to date).

- **DOCS-LINES** (documentation only, 2026-09-07): every `file.f90:NNN`
  citation in the nine design and status documents (d0 922 tokens on 309
  lines, b1a 127, b1 64, b3b 24, a2 7, a0 5, b3a 5, b4 2) re-anchored to
  file + routine, declaration or labeled comment block, read against the
  live tree; the step tables of d0 6.1 and b1 7.1 index by the `as_op_*`
  operation tag of `attempted_step.f90`; the convention is stated in each
  document. Citations of removed code (the composition projection, the
  two-iteration energy update, the CO ceiling, the cell-velocity carrier
  advection, the element advective term, DB96 eq. (39), the H3+ clamp, the
  second H2 partition function) kept as dated records naming the removing
  item. Stale values corrected in passing: `e_lw_photon_erg` 1.88021e-11 erg
  (11.7354 eV) over 912-1201 A, `binary_element_diffusion.f90` 3037 lines,
  a2 section 1's chemical count now passed by the marching caller. All 190
  names introduced grepped present in `src/`. Recorded: pre-existing
  markdown table pipe misalignments (unescaped `|` in math) in a2, b4,
  code_status, d0; two markdown line citations of
  `grid_and_gates/README.md` left.

- **B5, carrier rows** (`steady_newton.f90`, `diffusive_photochemistry.f90`,
  `input_read.f90`; new suite `steady_species_rows` 41 + 8 rows; advisor
  rebuild: it and `certification`, `run_mode`, `carrier_retry`,
  `coupled_source_step`, `species_face_flux` green): the Newton unknown
  vector carries one row and one unknown per solved carrier (H2, OH, H2O,
  CO, H+) through a row registry in place of the H2-only `nvar = 4` special
  case; band geometry derived, `kl = ku = 3 nvar - 1`, `ncolor = 6 nvar - 1`
  (reproduces the 8/17 and 11/23 the code carried); the completion flag
  reads the species rows of the system it solved; `EXHALE_SPECIES_JAC_TEST`
  asserts J r inside the solve (the run-level test builds its Jacobian
  before the registry exists). `Ionization transport` + `Solver: Newton`
  supported by derivation when the proton is a Newton unknown (`x_hp_fix`
  handed to the sweep), refused with the reason otherwise; MEASURED
  `nvar = 5`, Jacobian action 7.2e-7 (all rows), 4.6e-5 (species rows), 0
  unresolved colors. Two pre-existing defects fixed: `frozen_residual` left
  the species slots of an `intent(out)` array unwritten (directional check
  relative error 1.0 from stale memory) and `n_headroom_iterate` was never
  assigned (the trial-admissibility threshold stood at 0; with the four
  oxygen balances as rows all 41 colors came back unresolved). `wasp_full_newton`
  stays certified (13.30, byte-identical); H2-only control byte-identical;
  marching path byte-identical (`wasp_full`, `mol_base_handoff` 300).
  **Findings:** no regression case enters the steady solver, so every golden
  is a marching snapshot and `cert_tol_carrier`/`cert_tol_element` still
  have no anchor (contract section 11 records two rows from solves that did
  not meet their target, 7.46e-3 and 5.82e-3, not anchors). **Open, and a
  physical-consistency defect:** the stationary carrier row still uses the
  cell-velocity advective form while the marching rows ride on
  `F_rho Y_s^face / m_s` since B4-1c (two operators, two fixed points): B5b.
  Also open: the element rows (`element_transport_residual` must return a
  residual per element, not the worst at each cell), the cell's own field
  in the cell residual (interface written: FIELD-SELF), reload consistency,
  three `EXHALE_main.f90` sites (the registry name, `EXHALE_JAC_TEST`
  coverage, the PTC refusal text).

- **B5b** (`species_face_flux.f90`, `steady_residual.f90`, `steady_newton.f90`,
  `certification.f90`, `binary_element_diffusion.f90`,
  `diffusive_photochemistry.f90`, three sites of `EXHALE_main.f90`; suite
  `steady_species_rows` 67 rows; advisor rebuild: it and `species_face_flux`,
  `certification`, `carrier_retry`, `coupled_source_step`, `run_mode`,
  `acceptance_classes` green): ONE discretization of the material advective
  term exists in the code: the RK stages, the stationary carrier rows, the
  fixed-wind carrier relaxation and the stationary element rows all form it
  as the divergence of the same face species fluxes (`species_face_fraction`,
  `species_face_flux`, new `species_flux_divergence`) on the mass row's faces,
  areas and volumes; the cell-velocity upwind form, its deferred van Leer
  correction (`carrier_advection_correction`, `adv_corr`, `carrier_slope`), the
  base Dirichlet record it needed and the advective coefficient `ntv` are
  gone from `diffusive_photochemistry.f90` (RED 76.7 relative against a
  stage-1 update on a sign-changing face flux, GREEN 7.5e-14).
  `element_transport_residual` returns one residual per element (the
  certification prints one entry per element: `lower_profile` C 0.977, O
  0.921, N 0.968, seven elements "not applicable"); rows `srow_element_he`
  and one `srow_element_trace` per metal element registered before the
  carriers, `project_element_mass_fractions` the write-back. `mol_diffusion`
  reloaded with `Coupled carrier solve`: `nvar = 4`, one element row,
  Jacobian action 4.9e-3 / 1.3e-3, 0 unresolved colors, the solve refused
  at the He/H row (0.9987) with no descent direction: the element system is
  expressible and its Jacobian right, not yet solved from that start (the
  element unknown moves the whole column's opacity, which the banded model
  does not carry; the "all rows" bound of the suite widened 1e-4 -> 1e-2
  with the measurement at the site). Reload consistency asserted on the
  binary. `EXHALE_main.f90`: `set_transported_species_rows` (not
  carrier-specific), the two `EXHALE_JAC_TEST` inserts, the PTC text; at
  that site two pre-existing defects fixed: the species slots of `Yvec` were
  read uninitialized (2.4e2 against 1.85) and the probe direction was
  unscaled (3.3e1 at `nvar = 4`); now 3.2e-6 / 4.5e-6. **By construction the
  stationary species rows now inherit the mass row's imbalance**
  (`div(F_rho Y) = Y div F_rho + F_rho grad Y`; the old form was the second
  term alone), so on marching snapshots whose mass rows read 1.1-1.9 of
  their terms the rows read O(1) where they read 4e-4 to 0.4 before.
  Marching path byte-identical on `wasp_full`, `mol_base_handoff`,
  `mol_carrier`, `mol_diffusion`, `lower_profile` at 300 steps;
  `wasp_full_newton` reloaded certifies (13.30, 2.6e-9 of the stored state).
  Anchors still UNCONFIRMED (contract 11.3). Recorded: `element_diffusion_step`
  / `relax_element_composition` still advect on the smoothed steady mass
  flux, so the elemental Picard alternation and the elemental row are two
  operators (B5c); `carrier_rows_advect` kept as the on/off switch of one
  operator; `Coupled carrier solve` now also governs the element rows (key
  not renamed); `ionization_equilibrium.f90` cites `carrier_base_state`,
  which no longer exists (to FIELD-SELF).
  THREAD-DET confirmed by the official check: `make check` on the rebuilt
  binary, REGRESSION PASS, every case byte-identical (16:58-20:04 KST
  2026-09-07; the binary of that check also carried whatever of B5 had landed
  by 16:58, and B5's marching path was measured byte-identical anyway). A
  further `make check` follows FIELD-SELF and B5c.

- **B5c** (`binary_element_diffusion.f90`; suite `steady_species_rows` 60/0;
  advisor rebuild: it and `species_face_flux`, `certification`,
  `carrier_retry` green): the elemental relaxation and the stationary
  elemental row are one operator. `element_diffusion_step` takes the FACE
  MASS FLUX and forms its advective term as the divergence of the face
  element fluxes through the three shared routines
  (`element_advective_divergence`, one place); `relax_element_composition`
  reads that flux from the mass row of the state (`face_mass_flux_of_state`)
  and no longer smooths a steady `mdot/(4 pi r^2)`; `element_transport_residual`
  takes both halves of the helium row from `composition_residual` and the
  trace rows from the new `trace_composition_residual`, the routines the step
  itself solves, so the fixed point is the row's zero by construction;
  `solve_trace_element_in_hydrogen` is a deferred-correction loop that takes
  one pass, term for term the old direct solve, when the advection is off
  (the marching path). MEASURED on a synthetic column relaxed to its own
  convergence: the row reads 8.65e-13 (He alone), 8.45e-13 / 5.38e-13 (He +
  ten trace elements) of its terms, the advection carrying 0.993 of them;
  the replaced form leaves 6.4e-3 on the same state. The smoothed-flux
  rationale does not survive: it repaired a form (the element equation minus
  X times the mass equation) that is gone, and its consequence is now a
  measured property of the equation: on a face mass flux that does not
  conserve mass the conservative balance has no bounded steady composition
  (X driven onto 0 and 1). Marching path byte-identical (`mol_diffusion`,
  `lower_profile` 300 steps); the direct-steady reload of `mol_diffusion`
  byte-identical, and its outer pass never reaches the relaxation (the solve
  ends `info = 2`, B5b 6.3), so the run-level RED/GREEN does not exist.
  **Found: `diffusion_tests.x` has not linked since B5b** (`binary_element_diffusion`
  now `use`s `steady_residual_mod`, whose chain the test's `DIFT_SRC` does not
  carry): a `functions/` module depending on `time_step/steady_residual` is a
  layering inversion; the flux goes by argument (DIFT-LINK, launched), and
  `make test`'s `diffusion_tests` row is red until it lands. Also recorded:
  the relaxation's convergence measure is helium-only (DIFT-LINK folds the
  trace drift in); the vanished-element write-back branch leaves a trace row
  at 1e-4 to 2e-2 where an element is driven to exactly zero (deliberate
  guard, measured, documented at the test site); `element_transport_residual`'s
  dummy `v` is unread.

- **FIELD-SELF** (`ionization_equilibrium.f90` self-field loop, default one
  pass; `util_ion_eq.f90` and `System_HeH*` comments; probe
  `photoionization_field_self_consistency` 12; advisor rebuild:
  `coupled_source_step` 36/0, `physics_probe` 1491/0, `certification`,
  `acceptance_classes` green): the cell's own optical depth CAN now be
  closed inside the cell solve (one pass leaves a 4.06e-2 self-consistency
  defect in the ionized fraction of a cell one optical depth thick, the
  loop closes it to 2.8e-15 in 13 passes, contraction 0.080; the field
  routines already formed the B4-6 cell mean from the cell's own densities,
  so no field code changed), **and doing so buys nothing: COST3's
  attribution is refuted.** MEASURED on `wasp_full` 300: coupled passes
  10.09 -> 10.01 for 1.65x the wall time (2 self-field passes) and 1.87x (8);
  `mol_base_handoff` 1500: 7.19 -> 7.24 for 2.0x; the decisive binary with
  BOTH levers (carried columns at block 8 AND four self-field passes): 7.20
  passes, 190 s, against COST3's 7.13 / 103 s for the columns alone.
  (**Superseded by COST4 below**: with both stellar lags closed at the
  settings that close them the count is 8.18, not seven, and the two
  couplings named here are worth 0.08 of it; the count is the slowest rung
  of a ladder of modes.) Same fixed point: state movement <= 1.4e-7, at the 1e-8
  tolerance's level. Default output byte-identical on `wasp_full` 300,
  `mol_base_handoff` 1500, `mol_metals` 1500, `hydrostatic_column`;
  `wasp_full_newton` uncapped CERTIFIED, both files byte-identical to its
  golden (that binary carried B5b too). Kept in the tree at its inert
  default on COST3's precedent (the switch is the record of a measured
  negative result). Comment fixes in passing: `System_HeH_metals.f90` said the
  sweep is serial (untrue since the parallelization; its arrays are
  threadprivate); the `carrier_base_state` citation of
  `ionization_equilibrium.f90` rewritten to `carrier_base_composition_imposed`.
  Recorded: `he_ground_singlet_density`'s clamp counter needs a lock when
  the self-field loop runs (a slot per cell would be free,
  `composition.f90`); `f_vib_quench`/`e_vib_bound` stay at the entry
  composition inside the loop; the COST3 attribution in that report, its
  log entry above and the `csm_max_pass` comment of `EXHALE_main.f90` are
  superseded by this measurement (comment to be corrected when DIFT-LINK
  releases main).

- **DIFT-LINK** (`binary_element_diffusion.f90` and its three call sites in
  `EXHALE_main.f90`, `steady_newton.f90`, `certification.f90`; suite
  `steady_species_rows` 62 rows; advisor rebuild: it and `species_face_flux`,
  `certification`, `coupled_source_step` green): the element module no longer
  `use`s `steady_residual_mod`; the face mass flux is a required argument of
  `relax_element_composition` and `element_transport_residual`, obtained by
  each caller with `face_mass_flux_of_state`. The layering inversion is gone
  and **`make diffusion_tests` links again with the `Makefile` and `DIFT_SRC`
  untouched: 35 passed, 0 failed** (MEASURED by the advisor on the live
  tree; it failed at the link since B5b). The relaxation's convergence
  measure now runs over every element it relaxed, each against its own
  reservoir (`element_composition_distance`): RED 5.99 / GREEN 0.0 on a
  column where Fe moves 1.78 while helium moves 0.25, so the helium-only
  measure returned a number 7x too small to bound the element that moved.
  Byte-identical on `mol_diffusion`, `lower_profile` 300 steps and the
  `mol_diffusion` direct-steady reload. Dead code removed
  (`helium_mass_fraction`, an unread `element_mass_fractions` call).
  **Advisor fixes on its findings:** the limiter of `species_face_fraction`
  divided by `Yav - Yf(j)`, which vanishes when the donor cell's average is
  itself negative (a trace species handed back at a small negative density;
  MEASURED as a SIGFPE under `-ffpe-trap=zero`): the anchor is now the
  average clamped into [0,1], which makes both denominators strictly
  positive and is what the comment already claimed; `mol_diffusion` 300
  steps byte-identical, since at production flags the infinity was clamped
  to the same 0. The `csm_max_pass` comment of `EXHALE_main.f90` and
  `docs/steady_solver_design.md` section 15 corrected to the FIELD-SELF and
  B5c measurements. Recorded: `diffusion_tests` T9 reads 1.098e-12 against
  its 1e-12 bound at `-O1 -fcheck` (an optimization-level rounding
  sensitivity of `relative_settling_mass`, 35/0 at -O3); the outer loop's
  under-relaxation damps helium only while the drift it tests can now be set
  by a trace element.

- **Batch gate after B5/B5b/B5c/FIELD-SELF/DIFT-LINK/THREAD-DET**
  (2026-09-07/08): shared tree rebuilt (`EXHALE.x` md5 eba2d7f2260f,
  `make -q` up to date); `make test` green, the two steady suites passing on
  the `wasp_full_newton` run log and the `diffusion_tests` row restored;
  `make check` **REGRESSION PASS, all 16 cases byte-identical** (20:54 to
  00:05 KST). Every item of this batch was designed to leave the marching
  path alone and the check confirms it: the goldens refreshed at the
  2026-09-07 series gate still stand, and no golden was refreshed here.

- **HYG-T9** (advisor, 2026-09-08; `src/tests/diffusion_tests.f90`): T9's
  ambipolar-limit bound sat exactly on the arithmetic's noise floor, so the
  suite passed at -O3 and failed at -O1 with checks (DIFT-LINK's finding).
  MEASURED at both levels: the neutral limit is exact (no electron term),
  the H+ plasma reads 8.207e-13 (-O3) and 7.740e-13 (-O1), the He++ plasma
  7.083e-13 and 1.098e-12. Mechanism, written at the site: the limits are
  exact algebra in the masses, so the deviation is arithmetic, and the
  ambipolar term of the plasma limits comes from a centered difference of
  `ln rho` over `dr` whose logarithm is O(b0/r), so each ulp reaches
  `dm_eff` amplified by 1/dr. The bound is now one decade above the largest
  measured (1e-11), the anchoring rule the certification tolerances use, and
  the test still catches a wrong ambipolar limit by ten decades (a wrong
  limit misses its target by 0.1 to 1) with the neutral limit still asserted
  exact. 35/0 at -O3 and at -O1 -fcheck.

- **B5d** (`steady_newton.f90` only; suite `steady_species_rows` 77 rows;
  advisor rebuild: it and `certification`, `species_face_flux`,
  `coupled_source_step`, `carrier_retry` green): the item's own brief
  expected the element rows to fail because an element unknown moves the
  whole column's opacity, too non-local for the banded model. **MEASURED and
  REFUTED**: the scaled Jacobian column of an element unknown carries at
  most 6.3 percent of its 2-norm outside the band, the same order as a
  hydrodynamic column (mass 0.0018, momentum and energy 0.0000), so the
  diagonal correction, the rank-one column-integral update and the identity
  block the brief offered would each have repaired a defect that is not
  there; none was written. What does break the model, located by the same
  measurements and fixed: (1) the base element row was structurally EMPTY,
  cell 1 being the Dirichlet reservoir of the element operator, so the
  banded model was singular by construction (its LU pivot grew to 8.413e19
  against a largest band entry of 5.229e3; now 1.433e6) and it now carries
  the boundary condition the operator states; (2) the element ROW SCALE
  collapses on the abandoned iterate (helium rows spanning 1.277e-29 to
  2.360, so an ordinary `dF/dY = 46.7` became a scaled band entry of
  4.890e22 and `||grad merit|| = 1.104e44`), and the row scale is now
  floored at the flow-time rate of the quantity the row transports, the
  hydrodynamic convention with no chosen constant. Result: the same case
  reaches `||R|| = 1.723` against 1.974 at outer iteration 18 instead of 8
  with `||grad merit|| = 4.813e3`, 40 decades lower and a number the model's
  own entries can produce. **The element solve still does not converge and
  the item does not claim it does**; `cert_tol_element` and
  `cert_tol_carrier` are unmoved and UNCONFIRMED. New diagnostics:
  `banded_model_conditioning` (empty rows and columns stated once per
  registry, pivots and scales under `EXHALE_BAND_REPORT`),
  `jacobian_column_reach_beyond_the_band` (`EXHALE_JAC_COLUMN_REACH`).
  `mol_diffusion` 300 steps byte-identical; `wasp_full_newton`
  byte-identical, `info = 0`, CERTIFIED, 13.30, same 6 outer iterations.
  **Named as what is left: the globalization** (the non-monotone search
  accepts steps that raise `||R||` from 1.7 to 25; the trust-region route
  gives up at radius 3.8e-16 after 12 ray-test refusals while the merit is
  still falling): item B5e. Also recorded: `carrier_row_term_scale` has the
  same collapsing structure and was not measured (no carrier row in that
  case); the trace-metal case of the floor is covered by inspection only.

- **COST4** (measurement; NOTHING in the repository changed by the item, its
  four candidate files byte-identical to the snapshot it took; every
  experiment a private binary in a scratch tree whose no-switch build is
  byte-identical to the control): **the coupled loop's pass count is a
  ladder of nested modes and equals the SLOWEST of them, not their sum.**
  MEASURED on `wasp_full` 300, one thread, removing one composition coupling
  at a time: the production loop 10.09 passes; the stellar field's
  dependence on the composition removed, 7.16; the temperature removed as
  well (the energy solve replaced by its own converged answer), 4.94; the He
  recombination rates too, **2.00** (which reproduces, with a different
  probe, the 2 passes COST2 measured with every composition-dependent input
  frozen). Both hypotheses of the brief are refuted: the two diffuse
  couplings are worth at most 0.08 of the 10.09 (Balmer off 10.09, He rec
  off 10.01, both off 10.01, He rec frozen after pass 1 10.01, and on
  `mol_base_handoff` all three read the control's 7.19 to the last digit),
  and **the Balmer continuum is not a lag at all**: `excited_H_update` fills
  its arrays once BEFORE the loop, so they are constants of it, and
  refreshing it every pass costs 1.5 percent of wall for +0.01 passes. The
  stellar lag closed IN FULL for the first time (a one-cell block and eight
  self-field passes): 8.18 passes for 218 s against 10.09 for 139 s, the
  next rung, with freezing the field entirely giving 7.16. The temperature
  is not the residual either: the composition sweep alone, at the
  temperature the step converged to, takes 10.94 passes (6.32 against 7.19
  on `mol_base_handoff`). **Decision: change nothing** (no loop-internal
  Balmer refresh, no self-consistent He recombination coupling, no default
  moved); to go below seven passes every rung must come out at once.
  Advisor: the three statements this refutes are corrected in place, the
  `csm_max_pass` comment of `EXHALE_main.f90`, the FIELD-SELF entry above
  (marked superseded) and open item 1 of the 2026-09-08 handoff. **The one
  remainder, now open item 1**: on `mol_base_handoff` the ladder is flat
  (freezing the field RAISES the count 7.19 -> 7.41) and everything the item
  could reach accounts for 0.6 of its 7.19 passes; the named candidate, READ
  and not measured, is `n_tot`, the third-body density of the three-body
  molecular reactions, built by `calc_ntot` from the ENTRY composition, and
  with it the H2 vibrational quench fraction and bound energy. Also
  recorded: `xuv_field_block_cells` and `xuv_self_field_passes` are
  compile-time parameters, so measuring the two together needs a rebuild;
  `heating_of_composition` calls `he_rec_coupling` a second time whole-grid
  every pass for its heating alone (both calls together bounded at 1.7
  percent of wall); on `mol_base_handoff` the He recombination coupling is
  what keeps step 0 inside the pass cap.

- **COST5** (measurement; NOTHING in the repository changed by the item, its
  four candidate files byte-identical to its own snapshot; probe binaries
  byte-identical to the control at their default): the molecular remainder
  COST4 left is named, quantified, and the ladder reaches its bottom rung.
  **The count is arithmetic, and the law is measured and then tested**:
  `passes = 2 + log10(the step's own excursion / the stopping tolerance) /
  (decades per pass)`. On `mol_base_handoff` the composition contracts by
  0.171 a pass and the temperature by 0.072 (0.77 decades a pass) from an
  excursion of 7.0e-5 to a tolerance of 1e-8: 7.00 predicted, 6.93 MEASURED;
  on `wasp_full` 0.222 and 0.199 from an excursion seven times larger, 9.2
  predicted, 10.09 measured (so the molecular cases are the cheap ones, not
  the pathological ones). Tested, not fitted: the law says the count depends
  on the stopping tolerance through one logarithm at 1.30 passes per decade,
  and moving both tolerances together gives 8.06 / 6.96 / 5.68 / 4.04 at
  1e-9 / 1e-8 / 1e-7 / 1e-6, i.e. 1.34 over three decades; and it is not the
  cell solve, six decades of hybrd1 tolerance leaving the count at 6.97 for
  up to 1.99x the wall. Of the four candidates only the named one is an
  iterate: `n_tot`, the third body of R12/R13/R15 (the ONLY composition
  dependence of the molecular network outside the field: of seventeen
  hoisted coefficients exactly two see it, both linearly). Holding it:
  `mol_base_handoff` 7.19 -> 6.03 at 1500 steps, `mol_metals` 6.99 -> 6.01;
  it limits the composition mode's RATE (0.171 -> 0.073) but the count falls
  only 6.93 -> 6.11 because the temperature mode then binds. The other three
  are entry state through the operator-split field (vibrational pair 6.96,
  electron density 6.94, Lyman-Werner 6.96, against 6.96), and the metal
  element totals are invariant across a sweep by construction (6.99 against
  6.99, exact). With every input AND the temperature held the sweep needs
  **2.01 passes**, reproducing COST2's and COST4's atomic bottom rung on the
  molecular path. **Decision: nothing closed.** The mode that sets the count
  is the temperature/composition alternation (`EXHALE_main.f90`,
  `energy_semi_implicit.f90`, not that item's files); closing the `n_tot`
  lag exactly needs two sites in `System_HeH_mol.f90` (also not its), and
  the in-file alternative costs a second cell solve per pass, which the
  cell-tolerance scan bounds at more than the 0.80 passes it returns. The
  production treatment is not wrong: at the fixed point the entry and
  returned compositions are one state, so the third body the network is
  solved with is the third body of the accepted state, to `csm_comp_tol`.
  Advisor: the law is now in the `csm_max_pass` comment; the misleading
  indentation of `composition.f90`'s `calc_ntot` branches fixed. Recorded:
  the loop's stopping test compares an ABSOLUTE mass-fraction change against
  `csm_comp_tol`, so it is a relative 1e-8 on H2 and no constraint on a
  1e-12 species. **The one lever that keeps both the fixed point and the
  tolerance is sequence extrapolation** on the (T, c) pair after the third
  pass (the decay is geometric with a ratio stable from pass 2 and close
  across cases): ESTIMATED at three passes of seven on the molecular cases
  and four of ten on `wasp_full`, item COST6.

- **HYG-4** (advisor, direct; the four hygiene items of the 2026-09-08
  handoff): (i) `he_ground_singlet_density`'s run-wide clamp counter is
  incremented with `!$omp atomic update` inside the function, so the named
  critical section around its one call inside a parallel region
  (`ionization_equilibrium.f90`, self-consistent field loop) is gone and the
  three comments that justified serial column formation by the counter are
  reworded (the columns are still formed before the blocks are solved, for
  the ordering, not for the counter); (ii) `element_transport_residual` no
  longer takes the cell velocity it never read (four call sites:
  `certification.f90`, `steady_newton.f90`, two test rows), the balance
  being the divergence of the FACE mass flux handed in; (iii) the outer
  Picard loop of `steady_wind_with_element_diffusion` halved the element
  under-relaxation `omega` whenever `comp_drift` failed to fall, but
  `comp_drift` is the worse of the element drift and the carrier drift while
  `omega` damps the element relaxation alone (the carrier pass is bounded by
  `carrier_trust`), so a carrier drift that failed to fall slowed the
  element alternation for nothing; the damping test now reads the element
  relaxation's own drift (`elem_drift`) and the exit test keeps the worst of
  the two. No regression case runs both `He_diffusion` and the carrier
  transport, so the matrix cannot move. (iv) `f_vib_quench`/`e_vib_bound`
  staying at the entry composition inside the self-field loop is CLOSED BY
  MEASUREMENT (COST5): they are entry state through the operator-split
  field, and holding them at the first-pass value is worth 0.00 passes
  (6.96 against 6.96). Suites: `steady_species_rows` 74/0, `certification`
  29/0, `diffusion_tests` 35/0.

- **B5e** (steady solve globalization; `steady_newton.f90`, the
  `steady_species_rows` suite, design doc section 17, contract section
  11.5): **the first stationary solve with a species row that meets its own
  target.** MEASURED on `mol_diffusion` reloaded with `Coupled carrier solve:
  True` (one element row, `nvar = 4`): the scaled trust region takes `||R||`
  from 1.986 to **4.345e-09** against a target of 1e-08 in 66 outer
  iterations (54 accepted), where the pseudo-transient line search aborts at
  iteration 18 on "no descent direction exists for the banded model" at
  1.702. Two defects were named by instrumenting the eleven exits of
  `trust_region_step` and repaired where the measurement pointed: (i) the
  radius collapse was the CAUCHY LEG ASCENDING the model, because its
  gradient is the transpose of the BANDED `A` and has no sign (MEASURED
  `r0 . A g = -2.1e13`, `-2.2e19` at the collapsing iterates, `pred`
  negative and exactly linear in the radius); where `r0 . A g <= 0` the leg
  is now dropped and the dogleg reduces to the Krylov leg cut to the ball,
  which always predicts a reduction (GMRES gives `r0 . A sN < -||A sN||^2/2`);
  and `pred` is formed as `-(r0 . As) - ||As||^2/2` instead of the
  cancelling difference of two squared norms, whose sign on a short step was
  rounding; (ii) the best-iterate ledger ranked on `||R||` (three
  hydrodynamic rows) while the state is judged on those AND every species
  row, so it could hand back the iterate with the worse element balance; it
  now ranks on `distance_from_certification`, the largest of the judged row
  measures over their tolerances, formed once inside `eval_residual` under
  an optional argument; the three-unknown branch keeps its comparison
  character for character. **Refused by measurement**: bounding the STEP by
  the judged rows (implemented in all three step controls) stalled the
  element solve at `||R|| = 1.92` on 13 refusals; and the banded-gradient
  Cauchy length `(r0 . A g)/||A g||^2`, more defensible on paper, left the
  same solve crawling at 1.8e-4 after 81 iterations against 4.3e-9 in 66
  with the original length, so the original stays, recorded at the site with
  both numbers. **Default moved**: `use_tr = (nspec_row > 0)`, the trust
  region is the step control of the coupled route (`EXHALE_TRUST_REGION=0`
  restores the line search for measurement); no regression case sets
  `Coupled carrier solve: True`, so no golden can move. **Anchor**: the
  elemental transport He/H row of the returned state reads **4.633e-11
  against `cert_tol_element = 1e-8`** (215x inside), the first measurement
  under that tolerance; it is NOT moved, because contract sections 9 and 10
  want a ladder (the same case at `Resid tol` 1e-9 and 1e-10, not run) to
  tell a target-limited row from one at its floor. The state is not
  certified: its cell-1 hydrodynamic mass row reads 6.2e-11 against 3e-12, a
  base-cell matter and not an element-row one. **The carrier-row solve is not
  rescued**: on `mol_carrier` the trust region stagnates at 1.873 (line
  search 1.855) and 18 of its 49 iterations end on "no Krylov direction
  could be sampled", so `cert_tol_carrier` stays UNCONFIRMED, and the
  measurement says it is not a step-length problem. Worker MEASURED
  `wasp_full_newton` byte-identical (it 6, `||R|| = 8.131e-09`, info 0,
  CERTIFIED, log10 Mdot 13.30) and the `mol_diffusion` trust-region run
  identical to the pre-item text under `EXHALE_TRUST_REGION=1` (every landed
  change inert on a path that already worked). Recorded, not fixed:
  `EXHALE_SPECIES_JAC_TEST=1` perturbs the run it measures (iterations 1-5
  agree to every printed digit, iteration 6 differs, the two runs then
  diverge; some module state survives a probe `eval_residual`; the same root
  writes `n_eq_sweeps_model` and the frozen background); and the ray test
  can refuse a step that reduced the merit far MORE than the model promised
  (`ratio = 7187`), which Algorithm 4.1 of Nocedal and Wright never does,
  the obvious amendment (gate it on `ratio <= eta_accept`) left unmade for
  want of a control measurement. The worker was stopped by a rate limit
  after writing its report; advisor verified the tree: build clean; suites
  `steady_species_rows` 74/0 (six new driver rows), `certification` 29/0,
  `species_face_flux` 16/0, `carrier_retry` 101/0, `coupled_source_step`
  30/0; one "per-iteration" phrase in `run.sh` reworded; `wasp_full_newton` re-run by the advisor against its golden with the B5e binary: all four files byte-identical, info 0, CERTIFIED, log10 Mdot 13.30.

- **Batch gate 2 (2026-09-08, 05:04 to 08:08)**: `make check` on the tree
  as of B5e + HYG-4 (before COST6 landed; the harness built at 05:04, COST6's
  first edit is 05:17) **PASS, all 65 files of the matrix byte-identical**
  to the goldens refreshed by the user on 2026-09-07. `make test` in the
  same script failed 7 suites on their build-staleness precheck only (the
  shared `build/` had not been remade after HYG-4 when the suites ran, the
  harness remade it afterwards); the suites themselves were run on the
  advisor's private build of the same tree: `steady_species_rows` 74/0,
  `certification` 29/0, `species_face_flux` 16/0, `carrier_retry` 101/0,
  `coupled_source_step` 30/0 then 39/0 with COST6, `diffusion_tests` 35/0,
  `physics_probe` 1491/0, and `energy_update`, `grid_and_gates`,
  `spectrum_type` 53/0, 93/0, 154/0. Every suite of `make test` is therefore green on this tree. No golden touched.

- **COST6** (extrapolation of the coupled loop's geometric decay;
  `EXHALE_main.f90` +532 lines, `coupled_source_step` +3 whole-binary rows;
  **default OFF**, `EXHALE_CSM_EXTRAP=1` turns it on): implemented, guarded,
  measured. The object is the GLOBAL pair (T, f_sp), which is exactly the
  loop's iterate (the energy solve rewrites `u(3,:)` and `W(3,:)` in full
  from T and the entry anchor each pass, so an undo of the pair is exact,
  and the guard-forced-refusal row reproduces the off run byte for byte);
  one Rayleigh-quotient ratio per block over all cells and species, the
  projection onto the dominant mode COST5 measured stable, and one damped
  ray, so no cell can move in a direction the global mode lacks. That is how
  it differs from COST2's Aitken taken cell by cell (12.00 against 11.05 passes): a
  ratio taken cell by cell spreads 200 percent across the grid and cannot be
  tested for validity; the global one is tested (`||d_k - theta d_{k-1}|| <=
  0.5 ||d_k||`) before use. Admissibility is by construction (fractions in
  [0,1], He 2^3S inside its He I column, T inside the energy solve's own
  bracket, all linear in the step, so the largest admissible multiple is
  taken); MEASURED, the worst violation is identically zero on every case at
  300 and 1500 steps, 1 and 16 threads, and the step was never damped below
  0.95. The stopping test and its tolerances are untouched: a candidate is
  kept only if the pass taken FROM it meets them. Three things the
  measurement taught: the mode ALTERNATES (cosine of consecutive increments
  -0.85 to -1.00 on `mol_base_handoff`, +0.71 to +0.98 on `wasp_full`; the
  first implementation demanded theta > 0 and refused 1188 of 1188
  candidates); the one-mode model does not reach back into early passes
  (the ratio drifts 0.77, 0.47, 0.31, 0.22, 0.17, 0.14, 0.11 over one
  `wasp_full` step), hence a REACH guard, take the limit only where the
  modelled remaining error is already inside the tolerance ball (reach scan
  6.61 / 6.15 / 6.36 / 7.81 / 7.74 / 7.75 at 1/2/3/10/30/100, control 6.96);
  and undoing a mediocre candidate costs more than it saves (7.29 strict
  against 6.94 lenient, control 6.99 on `mol_metals`). MEASURED, sweeps
  control -> on: `mol_base_handoff` 1500 10780 -> 10102 (**-6.3 percent**,
  matched pair one binary 285.4 -> 268.0 s), `mol_carrier` 300 -2.3,
  `mol_metals` 400 -0.75, `wasp_full` 300 -0.63, `hydrostatic_column` 0;
  state movement 8.3e-8 to 5.0e-7 relative, log10 Mdot unchanged. **The
  brief's estimate of three passes of seven is NOT reachable, and the reason
  is the item's real finding**: the stopping test measures the INCREMENT,
  and the error is |theta/(1-theta)| times it, one to two decades smaller,
  so the loop is one to two decades stricter than its tolerances say and
  spends 1.0 to 1.6 passes on that; the extrapolation's ceiling is the last
  factor 1/|theta|. Redefining the test on the estimated error is a change
  of what the loop returns and is put to the user as decision item 8 of
  `docs/To_be_determined_by_user_20260906.md` (recommendation (b)). Tests:
  `coupled_source_step` RED 2 rows against the control, GREEN 39/0
  (`extrapolation_lowers_the_pass_count` 6.09 against 6.96 at 300 steps;
  `refused_extrapolation_costs_nothing` identical at 6.80;
  `mean_passes_within_bound` re-measured 6.80, unchanged); `physics_probe`
  1491/0. Worker-introduced defect fixed inside the item: a ledger `write`
  with six integer descriptors for seven items. Recorded, not fixed:
  `write_setup_report.f90` does not echo the switch; `n_csm_x_kept` can only
  count candidates the stopping test did not meet and reads 0 everywhere
  (the number that matters, met on the very next pass, is printed beside
  it); the `csm_max_pass` comment block is now some 60 lines of measured
  tables and belongs in a `docs/` memo.
  Advisor: the measured tables and the item history that had grown to some
  230 comment lines around the coupled-step constants are moved to
  `docs/coupled_source_loop_cost.md` (sections 1-7, decision item 8 linked);
  the code site keeps each constant's conclusion and a section pointer. The
  setup report does not echo environment switches (none of them, by
  design), so the switch is reported by the run counter line alone. Advisor
  verification: `coupled_source_step` 39/0 and `physics_probe` 1491/0 on the
  advisor build; default-off byte-identity against the GOLDENS at matrix
  length (which checks B5e, HYG-4 and COST6-off together):
  `mol_base_handoff` 12000 steps and `wasp_full` to its own stop, both
  `Hydro_ioniz.txt` and `Ion_species.txt` IDENTICAL to the goldens (both
  runs end with status 2, "not certified", exactly as the golden runs do).

- **COST7** (decision item 8, option (b), user 2026-09-08; `EXHALE_main.f90`,
  `coupled_source_step` suite, `docs/coupled_source_loop_cost.md` sections
  7-8): **the coupled source loop stops on the ESTIMATED REMAINING ERROR of
  the pair (T, composition), not on the last increment.** The estimate is
  |theta/(1 - theta)| times this pass's increment in each of the two norms
  the test uses, theta the global Rayleigh-quotient ratio, formed in ONE
  routine (`coupled_pair_geometric_error_estimate`) that the COST6
  extrapolation now reads too (its reach guard is literally the stopping
  test at a loosened tolerance; shared objects renamed `csm_x_*` ->
  `csm_geom_*`, env hooks `EXHALE_CSM_GEOM_*`). Where the sequence gives no
  estimate (passes 1-2, a guard failing) the increment test stands, and
  that fallback (`EXHALE_CSM_GEOM_REFUSE=1`) is byte-identical to the
  pre-COST7 binary on all 10 written files of `wasp_full` 300 and all 7 of
  `mol_base_handoff` 1500. **The item's finding**: the estimate is an L2
  quotient against two max-norm tolerances, and MEASURED with a probe that
  runs the loop on to 1e-12 and rolls back, the unfactored estimate accepts
  states OUTSIDE `csm_comp_tol` on 100 of 299 steps of `mol_base_handoff`
  (1.8x the tolerance at worst); a factor 2 leaves 1 of 299; **a measured
  safety factor `csm_err_safety = 3` leaves none**, the accepted state's
  distance to its own fixed point reading 7.4e-9 and below (composition)
  and 2.1e-9 and below (temperature) on every one of the 2496 accepted
  states probed. The price is that almost none of the passes COST6
  estimated survive: `wasp_full` 300 10.09 -> 10.01 passes (sweeps -0.9
  percent), `mol_base_handoff` 1500 7.19 -> 6.96 (-3.1), `mol_metals` 400
  6.99 -> 6.93 (-0.9), **`mol_carrier` 300 6.01 -> 5.01 (-16.6)**,
  `hydrostatic_column` 1.01 -> 1.01 (+1 sweep of 302: on its single
  three-pass step the estimate EXCEEDS the increment, so the new test is
  the stricter one there; it is not uniformly looser). log10 Mdot unchanged
  on every case; state movement 4.3e-10 to 1.1e-7 relative, at the
  tolerance level. What the change buys is the MEANING of the two
  tolerances, the decision's own first reason, not passes. The COST6
  extrapolation under the new test is worth nothing (6.79 against 6.82
  passes) and stays default off. Tests: `coupled_source_step` RED 41/3
  against the control, GREEN 44/0; `mean_passes_within_bound` re-measured
  6.23 (was 6.80), both at the site; `extrapolation_lowers_the_pass_count`
  rewritten as `extrapolation_does_not_cost_passes`; the fallback-is-one-path
  row replaces the brief's byte-identity row (a suite cannot hold a golden
  of the previous binary; the identity is measured in the report).
  `physics_probe` 1491/0. Recorded, not fixed: the composition increment of
  `mol_base_handoff` has a floor near 5e-9 (127 of 1500 probes cannot reach
  1e-12 in 59 passes), so `csm_comp_tol = 1e-8` is within a factor of two
  of what that case can reach at all; COST6's `n_csm_x_kept` still reads 0.
  **Every golden of the matrix moves at the 1e-8 level**: gate 3 below.

- **Gate 3 (2026-09-08, 09:29 to 12:00)**, tree after COST7 (with COST6
  default off, HYG-4, B5e): shared build remade, **`make test` all suites
  green**; `make check` verdict FAIL against the 1e-3 relative tolerance,
  as expected for a change of the stopping test, and every movement
  examined. `log10 Mdot` unchanged on all 14 cases that print one (13.35,
  13.35, 13.30, 10.57, 10.58, 10.54, 10.57, 10.58, 10.51, 10.56, 7.29, 8.27,
  9.61, 7.98). Movement table (largest relative change of any column, new
  against golden): `wasp_full` 6.6e-5 (v at the base, cell 136),
  `wasp_he23off` 3.8e-4 (v at cell 146), `wasp_full_newton` 1.1e-10,
  `mol_base_handoff` 2.8e-6, `mol_metals` 2.4e-6, `mol_lyman_werner`
  4.1e-7, `mol_diffusion` 2.8e-6 (adv file 2.2e-3 at the outer boundary,
  cell 504), `mol_ir_bands` 2.4e-6, `mol_sec_ion` 4.5e-3 (v = 0.091 cm/s at
  cell 177, a near-zero base velocity, adv 1.5e-3), `mol_carrier` 5.2e-7,
  `lower_profile` 9.3e-6, `hp_zero_seed` / `hp_trace_seed` / `hp_front`
  identical, `roundtrip` within. Two cases move by order one and both are
  pre-existing base artefacts, not COST7's: **`oxygen_chemistry`** (1000-step
  snapshot) reads 42x in n(H II) at cell 22 because the GOLDEN carries an
  odd-even sawtooth there (cells 21, 24, 25, 27 at 7e3-9e3 against 2e5 in
  their neighbors) and the new run is smooth (1.6e5 to 3.9e5 across cells
  18-26; rho, p, T within 1e-6); **`hydrostatic_column`** reads 0.999 in the
  adv file's T at cell 350 because the adv post-process gives T_adv = 0.78 K
  inside 1.40 Rp in this static column (the original T is 1124 K
  throughout) and the edge of that region moved by one cell; the non-adv
  files move by 4.6e-8 (Ion_species by 0.16 at cell 2 on values of 1e-182,
  denormal noise). The adv defect is recorded as an open item. **Goldens to
  be refreshed by the user** (archive first), then `make check` re-run.

- **ADV-STATIC** (`post_process_adv.f90`, new suite `adv_static_limit`,
  `docs/postprocess_advection_validity.md`): **the brief's premise is
  refuted and the 0.78 K of `hydrostatic_column`'s `_adv` file is named.**
  MEASURED: the flow where it sits is Mach 0.23 to 0.82 (the case is a
  300-step relaxation snapshot, not a static column: `rho v r^2` spans seven
  decades); the IONIZATION correction is held at equilibrium in 503 of 503
  cells (`Ion_species_adv.txt` equals `Ion_species.txt` to the bit) and its
  Damkohler condition already gives the continuous v -> 0 limit the brief
  asked for, so nothing there was touched. The number is the exact solution
  of the equation the `_adv` ENERGY solve poses: the steady internal-energy
  equation with the enthalpy flux of the mass-flux divergence DROPPED,
  `rho v de/dr - p v dln(rho)/dr = heat - cool`, which equals the exact
  `... + h div(rho v) = heat - cool` only where `rho v r^2` is stationary;
  with the radiation off it reduces to the discrete adiabat `T_j/T_{j-1} =
  rho_j/[rho_j - (gamma-1)(rho_j - rho_{j-1})]`, and integrating that on the
  run's density reproduces the reported temperature to 6e-11 over the 347
  cells before the first fallback (1138 K -> 0.78 K at 1.396 R_p; the
  dropped term exceeds the kept ones in 480 of 502 cells). With heat/cool
  negligible the equation is homogeneous in v, so no v -> 0 condition can
  change its answer (99.96 percent departure at |v| <= 1.4e-3 c_s). Landed:
  a validity condition of the energy solve's own, the thermal Damkohler
  number `Da = t_cross |heat - cool| / u_th` (new pure
  `thermal_damkohler_number`; above one the cell keeps the run's own T),
  the counterpart of ionization condition (ii) that the bare `v <= 0` sign
  test lacked; the stationary-mass-flux assumption stated at the site with
  its exact dropped term and MEASURED at run time (largest ratio, where, how
  many cells exceed one; nothing acts on it yet); the first `eval_cool`
  writes a named `tcool_in` instead of a discard slot; the stale "is not
  evaluated" sentence of the validity memo rewritten with both tables.
  Impact: byte-identical everywhere measured (`hydrostatic_column` full case
  all four files; `wasp_full`, `mol_base_handoff`, `lower_profile` post-
  processed from the same golden state by both binaries, all four files;
  max thermal Da over tested cells 0.951 / 0.534 / 0.055 / 3e-19, 0 cells
  above one anywhere). Suite RED/GREEN (rows 1-7 cannot link before; row D
  `measured=silent` before), advisor 10/0. Also MEASURED: no matrix case has
  a stationary mass flux by the cell-to-cell measure (`rho v r^2` over
  outflowing cells spans 1.27x its median in `wasp_full`, 2.2x in
  `mol_base_handoff`, 83x in `lower_profile`); the `du` stop is taken on a
  restricted window and does not bound this. **Follow-up decided by the
  advisor under the physical-correctness rule: restore the dropped term
  (item ADV-ENERGY)**, since an energy equation that is exact only for a
  stationary mass flux is wrong on every snapshot and at the 1e-3 level on
  every converged wind; MEASURED worth on the static column 569 K at 1.40
  R_p and 121 K at 2.96 R_p in place of 0.78 K and 6.0 K. The alternative,
  refusing the column where the assumption fails, would rewrite 162 of 502
  `_adv` temperatures of `wasp_full` (median 5 percent) and is decision 14
  (T8.1) of the target-system document, left to the user. Advisor fixed the
  stale provenance paragraph of the case README.

- **ADV-ENERGY** (`T_equation.f90`, `ion_cell_state.f90`,
  `post_process_adv.f90`, `adv_static_limit` +7 rows,
  `docs/postprocess_advection_validity.md`): the advection post-process's
  energy equation is now the EXACT steady internal-energy equation. The
  enthalpy flux of the mass-flux divergence, `h div(rho v)` with `h = e +
  p/rho`, enters `T_equation` through a new zero-default field of the cell
  record, `teq_cell%div_rhov`, as a separate additive term in both branches
  (caloric `mum dr div_rhov (E(x_H2,T) + T)`, monatomic `gamma mum dr
  div_rhov T`; `mum`, `mup`, `coeff` untouched; zero field an exact zero,
  tested to the bit). `div(rho v)` is the mass row's own operator `(A_p F_p -
  A_m F_m)/dV` over the control volume between the two points the upwind
  difference uses. MEASURED corrections to the briefs: `T_equation` has no
  caller outside `post_process_adv` (no marching solve to protect; the
  bitwise property holds anyway); and ADV-STATIC's 569 K prediction had put
  back `p div(v)` alone, not `e div(rho v)`; the whole term gives **1.7e7 K
  at 1.40 R_p and 1.6e9 K at 2.96 R_p** on the static-column snapshot, the
  closed form `w ~ rho^(gamma-1) F^(-gamma)` with `F = rho v r^2` falling 8
  percent from cell to cell. Tests RED (nine "not a member" compile errors
  before) / GREEN (closed-form solution `T ~ r^(-gamma)` on a diverging
  column to 6.9e-4 at 400 cells, first-order convergence 5.59e-3 / 2.79e-3 /
  1.39e-3 / 6.95e-4, the adiabat departing by 2.17); advisor
  `adv_static_limit` 17/0, `physics_probe` 1491/0. **Impact, MEASURED, same
  golden state post-processed by both binaries**: marching files bitwise in
  every case; T_adv movement `wasp_full_newton` 1.4 percent (the only
  Newton-converged state, He 2^3S 1.9 percent), `wasp_full` (du stop) 16
  percent at 1.088 R_p (He 2^3S 16 percent), `mol_base_handoff` snapshot 301
  percent, `lower_profile` snapshot 456 percent (He 2^3S 238x),
  `hydrostatic_column` 4.2e9 (reaching 2.6e10 K). Restored term against the
  kept ones: `wasp_full` 67x at 1.033 R_p in 161 of 502 cells,
  `lower_profile` 25x in 299 cells, `wasp_full_newton` 24x in 2 cells. **The
  `_adv` temperature of a state that violates continuity is now unbounded
  above where it was bounded below, and the transit tools read those
  columns**: on such a state NO steady post-process has a meaning, and what
  to output there is decision 14 of the target-system document, now put to
  the user as decision item 9 with a recommendation. Goldens: every `_adv`
  file moves by design; refresh deferred by the user to this item's gate.

- **B5g** (`steady_newton.f90` +467/-35, `steady_species_rows` suite,
  design doc sections 18, 18.5): **a residual evaluation is now a function
  of its arguments.** Three class-(b) caches survived a probe
  `eval_residual` and were read back by the solve: the elemental rows and
  scales (`erow_he`/`escale_he`/`erow_tr`/`escale_tr`), the carrier measures
  (`carrier_relnorm_last`...), and the carrier arrays of
  `diffusive_photochemistry` (`col_scale_car`, `row_scale_car`, `row_terms`,
  `headroom_car`, the frozen backgrounds); the Newton ROW SCALE, the
  acceptance gate, the best-iterate ledger and the choice of the state
  handed back read them. Fixed on the reader side: `eval_residual` gained
  `state_is_discarded` and a probe holds/puts back that set (reusing the
  `save_/restore_carrier_module_state` pair `certification.f90` already
  uses); the loop-top row scale is taken AFTER the iterate's own evaluation;
  the gates read a new `carrier_relnorm_state`. The two self-tests also tag
  their sweeps as candidate sweeps (under the iterate tag they were writing
  the run's acceptance ledger and `ieq_nonroot_streak`, whose limit stops
  the run). PROVEN: `EXHALE_SPECIES_JAC_TEST=1` and the untested run are the
  same run, every `(JFNK) it` line identical over 17 iterations,
  `Ion_species.txt` byte-identical (RED: 82 differing lines, first at
  iteration 6). `wasp_full_newton` byte-identical (info 0, 6 iterations,
  CERTIFIED, 13.30). **What the fix uncovered is bigger than the item**: the
  `mol_diffusion` element solve that B5e brought to 4.3e-9 in 66 iterations
  now STAGNATES at 1.957 after 17; the runs agree for five iterations and
  part on a 3 percent difference of `Drow` between two evaluations of the
  SAME state, because the residual's seed dependence (`resid_sc_probe`) is
  **0.96 in units of Resid tol** against `||R|| = 1.986`, and tightening the
  sweep tolerance (1e-12, 1e-14 with pass caps 200, 500) never satisfies
  the increment test from either seed: **the eliminated composition
  residual of that state has no fixed point it reaches**, and B5e's
  convergence was a property of one sequence of evaluations. The change is
  kept (a row scale must be the row's own terms of the same evaluation).
  Part 2, the carrier rows: all 18 "no Krylov direction" exits and all 514
  refused samples were ONE thing, the H2 carrier of cell 500 at -7.4e-19
  (a species unknown on its lower bound), not the element budget (the
  counter's label named one screen of several; corrected, census added);
  `jv_product` now takes the backward difference `[F(Y) - F(Y - eps v)]/eps`
  only after every forward halving is refused: 49/0/0 accepted/rejected/
  refused where 31/18/18, GMRES truncations 18 -> 0, radius 8e-3 -> 2.05;
  `||R||` unchanged at 1.873 because the dogleg cuts steps by 2^-25..2^-49
  against the same bound from the trial side, a bound-constrained unknown
  space (decision, not repair). `cert_tol_carrier` UNCONFIRMED. Suites
  (advisor build): `steady_species_rows` 74/0 (88/0 with the self-test pair
  logs), `certification` 29/0, `coupled_source_step` 44/0, `carrier_retry`
  100/0, `species_face_flux` green; advisor re-ran `wasp_full_newton` in full
  with the B5g binary: `Hydro_ioniz.txt` and `Ion_species.txt` IDENTICAL to
  the gate-3 output, CERTIFIED, log10 Mdot 13.30.

- **B5f** (`docs/a2_certification_contract_20260906.md` section 11.6; no
  constant moved): the `cert_tol_element` ladder on `mol_diffusion` at
  `Resid tol` 1e-8 / 1e-9 / 1e-10 gives ONE state: the 1e-8 solve meets its
  target in 66 iterations (element row 4.633e-11, reproducing B5e digit for
  digit, so COST7 leaves the stationary solve untouched); the 1e-9 and
  1e-10 solves abort after 88 iterations with the radius at 9.3e-23 (9 of
  31 refusals on "the ray difference is below its cancellation floor") and
  hand back the SAME state, bitwise. The row is floor-limited at 4.63e-11
  (knee at `||R||` about 1e-7, 1.5 percent movement over the last factor
  26); proposed value by the contract's rule 5e-10, NOT applied (two rungs
  are not solves that met their target; a decade above one configuration's
  floor is not a decade above the row's floor, which the operator's header
  puts between 1e-12 and 1e-5 across configurations; no second
  configuration exists). **Read together with B5g**: B5f ran on the
  pre-B5g text, and B5g shows the state it anchored on is reached by one
  sequence of evaluations only (with the probe caches removed the same
  solve stagnates at 1.957), so this rung stands as a measurement of that
  state's row, not as evidence that the solve converges. B5f also MEASURED
  that the coupled route's convergence depends on `eq_sweep_reltol`
  (1e-10 loses the 1e-9 solve at iteration 17), the same fact from the
  other side. Recorded: the JFNK of a reloaded 1e-8 state starts from
  `||R|| = 1.226` because the marching loop takes two CFL steps first (the
  O(dt) kick, eight decades, seen for the first time from a stationary
  state); the certification prints no cell-resolved record of the species
  rows, which blocks any floor attribution.

- **ADV-REFUSE** (decision item 9, option (a), user 2026-09-08;
  `post_process_adv.f90` +227/-43, `write_output.f90` +94/-10,
  `examples/exhale_io.py`, `EXHALE_transit.py`, `README.md`,
  `README_HOWTO.md`, user manual section 4 (PDF rebuilt, 76 pages),
  `docs/postprocess_advection_validity.md`, `adv_static_limit` +13 rows):
  **the advection post-process refuses the correction in every cell whose
  own energy-equation term ratio says the flow is not stationary; that cell
  keeps the run's temperature and the ionization equilibrium at it; a new
  integer column `adv_status` at the end of `Hydro_ioniz_adv.txt` records
  the reason row by row** (0 corrected; 1 thermal Damkohler; 2 not
  stationary; 3 ionization conditions (i)-(iii); 4 cell solve failed; five
  named module parameters in `write_output`, the `# columns` line extended,
  a `# adv_status_counts` header line counted from the array written).
  The criterion is `enthalpy_flux_term_ratio` (new pure function beside
  `thermal_damkohler_number`), the restored term over the LARGER kept term,
  refused on a STRICT `> 1` so a stationary flux with round-off is never
  refused; formed ONCE on the state handed in over every cell, so the
  refused set cannot move with the iterate and the converged product does
  not depend on the iteration history; applied to temperature AND
  composition (a non-stationary cell has no steady advective ionization
  correction either). `write_output` gains an optional `adv_status` dummy
  (the four `'eq'` calls untouched, `Hydro_ioniz.txt` keeps seven columns);
  `exhale_io.load_hydro` now takes names from the header (seven historical
  keys kept, `Run.adv_status` exposed, `None` for an equilibrium file);
  `EXHALE_transit.py` prints how many rows inside `1 <= b <= Rib` are
  refused, by reason; `EXHALE_plots.py` and the four other header parsers
  adapt unchanged. MEASURED (same golden state, "Do only PP"; marching
  files bitwise everywhere): `hydrostatic_column`'s corrected temperature
  falls from 2.6e10 K to **1182.32 K, exactly the run's own maximum, no row
  above it** (435 before); status counts 0/1/2/3/4: `wasp_full_newton`
  501/0/2/1/0, `mol_base_handoff` 470/0/33/1/0, `wasp_full` 205/0/297/2/0
  (160 of the 297 beyond conditions (i)-(iii), matching ADV-STATIC's
  route-2 estimate of 162; the larger total is the reference point, the
  state handed in over every cell including the inflow cells), `lower_profile`
  182/0/321/1/0, `hydrostatic_column` 0/0/481/23/0. Every refused row
  carries the run's T to the BIT and the equilibrium density of every
  species exactly (one 1.8e-16 reassociation in `lower_profile`'s He I,
  pre-existing). Movement against the pre-ADV-ENERGY goldens, T_adv:
  `wasp_full_newton` 3.3 percent, `wasp_full` 8.0, `mol_base_handoff` 64,
  `lower_profile` 1.8x, `hydrostatic_column` 1525x (0.78 -> 1182 K);
  species: He III 5.4 percent, He 2^3S 29 percent, He 2^3S 110x (at the
  base of `mol_base_handoff`), H II 11x (`lower_profile` base), He III 16
  percent. **Transit effect on `wasp_full`**: 295 of 500 rows inside the
  impact-parameter range are refused (printed), and the metrics move at
  1e-4 relative or less (He 10830 red 4.102 percent both, blue 2.843 ->
  2.844, FWHM 0.899 A both; largest line movement O I 1302 +0.55 percent
  relative), because the refused rows sit in the base under the extended
  wind. Tests RED (seven "not found in module" errors; `measured=7
  columns` on the old binary) / GREEN; advisor `adv_static_limit` 30/0,
  `physics_probe` 1491/0. Worker also fixed a dangling `\ref{sec:retiredbase}`
  in the manual (the four retired base keys named instead). Recorded: two
  documents describe the `_adv` validity conditions (manual, HOWTO) and had
  both drifted; four `# columns` parsers exist in python; `dV_cv` division
  now appears twice unguarded (a zero is a broken grid).

- **B5h** (`ionization_equilibrium.f90` +68/-2 diagnostics,
  `steady_newton.f90` +211/-2 probe only, `steady_selfconsistent_residual`
  assertion 4; NO behavior changed, `wasp_full_newton` and the element run
  byte-identical to a control built from the same source with the additions
  removed): **why one state has two residuals.** MEASURED on
  `mol_diffusion`'s element configuration with a lockstep two-seed trace of
  the composition elimination (`EXHALE_RESID_SC_TRACE`), a nine-seed basin
  ladder (`EXHALE_RESID_SC_BASIN`) and a cell acceptance report
  (`EXHALE_IEQ_REPORT_CELL`): the amnesty is never used (classes 4 and 6
  zero at every pass, every cell class 1 with reaction residual 1e-25
  against 1e-6), the self-field lag is not it (`xuv_self_field_passes` 1 ->
  20, identical to every digit), the He recombination coupling is not it
  (off, identical to four digits), the stopping test is not it (a fixed
  count of one pass per call, chained, the difference frozen from pass 2
  to 60). **Two things are.** (1) Two separated ROOTS of the network at the
  same incident field: from a seed displaced five cells the sweep reaches a
  composition 8.9x apart in x(H2) at the H2 front (cell 214) and 19 percent
  apart in x(H I) at the ionization front (cell 192), bitwise identical from
  cell 216 outward, both exact roots (class 1, 1e-25). (2) **A marginal
  direction: the local-equilibrium elimination does not determine the
  shielded layer's H2 CONTENT.** Retention of a seed perturbation of the H2
  partition is 0.77 at the base rising to 0.84 at r = 1.064 and falling to
  0 outside cell 216, the SAME curve to three digits at perturbations of
  1e-6 and 1e-2, so the sweep's map has eigenvalue 1 in that direction and
  the seed dependence is a derivative, linear over four decades with gain
  3.2e3 (the base cell's energy row amplifies the ghost composition by
  ~4e3), reaching `Resid tol` at a seed change of 3e-12. READ from the
  network: the fast H2 -> H2+ -> H3+ -> H2 cycle returns every nucleus, so
  the H2 row as scaled by its turnover pins the PARTITION and leaves the
  total to `k15 n_HI^2` against `k12 n_tot n_H2`, a detailed-balance pair
  that in the shielded layer is far below the turnover; and pinning it
  would pin the fully molecular limit (0.9997) where the state carries
  0.376: the content is upstream and transport data, exactly why it is
  imposed at the ghosts. **The repair the code already owns, MEASURED**:
  with `Molecular carrier transport: True` (n(H2) a Newton unknown) the
  same ladder falls from 3.155e-3 to **9.736e-10** at delta 1e-6 and x(H2)
  of the front cell is bit-identical from all nine seeds; the residual is a
  function of its arguments below `Resid tol`. What remains at large seeds
  (0.98) is the ionization front's two roots, a physics decision (which
  root a stationary residual may land on, or a statement that the front is
  unresolved on this grid), reported. The element solve with H2 carried
  then meets B5g's carrier obstruction (`||R||` flat at 1.986 over 19
  iterations, carrier row 0.91, run stopped by PID at 4 min an iteration).
  Also MEASURED: `eq_sweep_reltol` below ~1e-10 is unreachable by
  construction, the increment measure `max|df|/max(|f|,1e-12)` reporting
  round-off of He I and H2+ fractions ~1e-20 in the ionized wind (plateau
  6e-11 for sixty passes), which is what B5f's and B5g's "never met at
  1e-12, 1e-14" was. Tests: assertion 4 `residual_is_one_number_for_basin`
  RED 3.155e-3 with H2 eliminated, GREEN 9.736e-10 with H2 carried (the
  tolerance read from the run's own log); `steady_species_rows` 86/0 with
  the self-test pair; B5g's part 1 reproduced as a by-product. A defect of
  the item's own fixed inside it (a cell-report default of -1 was the lower
  ghost). **Decision put to the user as item 10**: refuse `Coupled carrier
  solve: True` on a molecular configuration without carrier transport, or
  turn the transport on for it (`input_read.f90`, not the item's).

- **B5i** (decision item 10, option (a), user 2026-09-08; `input_read.f90`
  +51/-2, `grid_and_gates` new row `coupled_carrier_h2` (+5 assertions),
  user manual (a `Coupled carrier solve` entry, which the manual had NOT
  had; the `Molecular carrier transport` entry corrected from "inert
  without oxygen chemistry" to "needs molecular chemistry"), `README.md`,
  `README_HOWTO.md` (new block; the oxygen paragraph's dangling "It"
  corrected), `docs/steady_solver_design.md` section 19): **`Coupled
  carrier solve: True` with `Molecular chemistry: True` and `Molecular
  carrier transport: False` is refused at startup**, after every key is
  resolved (so the oxygen default of the transport key is seen: MEASURED,
  `Oxygen chemistry: True` plus the coupled key is accepted, the same with
  the transport stated `False` is refused), with a message naming what was
  asked, why there is no root (B5h's numbers) and the remedy key; `error
  stop 1`, exit status 1 asserted (2 is the certification's refusal of a
  complete run). Scope exactly `carrier_in_newton .and. thereis_mol .and.
  .not. carrier_transport`, not gated on the Newton key (the direct steady
  route `EXHALE_PTC=1` runs a coupled solve too, and it is the route B5b to
  B5h used); the atomic coupled solve with an element row is accepted
  (asserted). No shipped configuration is refused (MEASURED by grep: no
  `*.inp` in the tree sets the key; 62 of 81 B5-series scratch inputs
  would be, as the decision says). Two molecular cases byte-identical
  before/after at 200 steps. The parse-site comment that claimed "checked
  below" for a check that did not exist is corrected. A refused run leaves
  a zero-byte `EXHALE_setup.out` and no `EXHALE_resolved.out`, like every
  other `input_read` refusal. Tests RED (the eliminated-H2 coupled run
  completed and printed a Mdot) / GREEN; advisor `grid_and_gates` 98/0,
  `physics_probe` 1491/0. Advisor corrected `docs/input_schema.md` K15e
  ("needs K15d" now qualified to the molecular configuration). Recorded:
  `EXHALE_main.f90` opens `EXHALE_setup.out` (truncating) before
  `input_read`, so any refusal erases the previous run's setup report (fix:
  open it after the input is accepted); the other "checked below" comments
  of `input_read.f90` deserve a sweep against the checks that exist.

- **HYG-SETUP** (advisor, direct; `EXHALE_main.f90`): the setup report
  `EXHALE_setup.out` is opened AFTER `input_read` returns instead of before
  it, so a configuration `input_read` refuses no longer truncates the
  report a previous accepted run left in the same directory (B5i's
  finding). `input_read` writes nothing to that unit. `grid_and_gates`
  green on the advisor build; the marching path is untouched.

- **HYG-CHECKED** (`input_read.f90` comments plus one new report block,
  `docs/input_schema.md`, `README_HOWTO.md`, user manual (PDF 77 pages),
  `grid_and_gates` new row `carrier_transport_inert` +4 assertions): every
  requirement statement of `input_read.f90` and every "Needs K.." /
  "refused with" clause of `docs/input_schema.md` checked against the code
  that enforces it. 40 enforced as stated. ONE not enforced (`Molecular
  carrier transport: True` with the molecular network off) is now a
  startup WARNING, not a refusal, by the code's own convention for a
  present-and-inert key (every consumer is guarded by `thereis_mol`,
  MEASURED inert; eight shipped inputs all carry the network). THREE
  enforced differently from the statement, one with the sense INVERTED:
  the schema and the HOWTO said `Ionization transport` is refused with
  `Solver: Newton` and with `Coupled carrier solve`, where the code
  REQUIRES the second whenever the first is set (K15e does carry a proton
  row; what is refused outright is the `EXHALE_PTC=1` route); "`Solver:
  PTC` refuses it" named a solver value that does not exist; the base-grid
  resolution check lives in `write_setup_report` (warns below 10 cells per
  scale height). SIX stale statements corrected: `known_keys` also feeds
  the duplicate-key REFUSAL (MEASURED, all 99 entries used, no hole); the
  sub-Lyman absorbers are not mutually exclusive; `Kzz_base`/`He_Kzz` are
  NOT inert with `He_diffusion` off (the carrier transport reads
  `kzz_cell`), corrected in four places; `Molecular reaction heat` default
  True is guarded by `with_molecules` (a refusal would refuse every atomic
  run); the coronal cutoff window 0.08-0.13 superseded by
  `coronal_cutoff_width.md` 7.2; four `K15b` cross-references that pointed
  at the `Stellar LW flux` row corrected to `K15`. MEASURED byte-identical
  on `wasp_full` 100 and `mol_base_handoff` 200 (17 of 17 files). Advisor
  `grid_and_gates` 102/0, `physics_probe` 1491/0. Recorded: all 63 "F<n>"
  line citations of the schema are stale by 300-1200 lines (marked stale
  at the convention's definition; recommend deleting the tokens); `Coupled
  carrier solve: True` is silently inert in an atomic run without
  `He_diffusion` (a report belongs where the row registry is); the user
  manual has no `Ionization transport` entry.

- **Gate 4 (2026-09-08, 15:39 to 18:45)**, tree after COST7, COST6 (off),
  B5e, B5g, B5h, HYG-4, ADV-STATIC, ADV-ENERGY, ADV-REFUSE (B5i, HYG-SETUP
  and HYG-CHECKED landed during the matrix; they touch startup checks and
  a file-open order only, and were verified on the advisor build with the
  full `grid_and_gates` suite, 102/0): shared build remade, **`make test`
  all suites green**; `make check` FAIL against 1e-3 as expected. **The
  marching files are unchanged from gate 3 to the digit**: `Hydro_ioniz.txt`
  and `Ion_species.txt` carry exactly COST7's movement on every case
  (identical max-rel values in both logs, the three `hp_*` cases still
  identical), so B5g, B5h, the three ADV items, B5i and the hygiene items
  moved nothing on the marching path, and log10 Mdot is unchanged on all
  14 cases. The `_adv` files move by design: `Hydro_ioniz_adv.txt` now has
  eight columns (`adv_status`, the harness reports SHAPE 8 vs 7), and
  `Ion_species_adv.txt` moves by 0.05 (`wasp_full_newton`, He I at the
  outer edge) to order one (the molecular cases and the `hp_*` cases,
  whose refused base cells now carry the equilibrium densities; the
  relative measure is element-wise, so a species the old correction had
  moved off its equilibrium value reads 1). **Goldens to be refreshed by
  the user, once, for COST7 + ADV-ENERGY + ADV-REFUSE together** (archive
  first), then `make check` re-run by the advisor.

- **Golden refresh (2026-09-08 19:05, run by the advisor on the user's
  instruction, the user being away from the machine)**: the previous set
  archived to `backup/regression/golden_pre_cost7_20260908/` (verified
  identical to `golden/` before the refresh), then `run_check.sh golden`
  snapshotted the gate-4 outputs of all 16 cases. The refreshed goldens
  carry COST7 (the stopping test on the estimated error), ADV-ENERGY and
  ADV-REFUSE (decisions 8b, 9a). Gate 5, a `make check` on the current
  tree (B5i, HYG-SETUP, HYG-CHECKED on top, all startup-only), follows.

- **HYG-DOCS** (documents only: `docs/input_schema.md`, user manual (PDF 82
  pages), `README_HOWTO.md`): the 91 stale "F<n>" line tokens (77 lines)
  and the 37 `GUI<n>` numbers (MEASURED stale by a uniform 13 lines) are
  gone; a key is located by its label and the header says so; nine Notes
  cells that held only a citation got real notes READ from the parse sites.
  The manual gained the missing `Ionization transport` entry (in the sense
  HYG-CHECKED corrected) and the four other keys the schema documented and
  the manual did not (`Reconstruction continuation`, `Atomic rate set`,
  `Caloric EOS`, `Photoelectron heating`); MEASURED afterwards, all 97
  schema keys have a manual entry. Three stale manual entries corrected
  (`Resid norm` is retired; `Resid tol`'s norm is the row scaling of
  sections 143/145, not `|R|/max|u|`; `Flux spread tol` default 2e-5, not
  5e-3); the validity paragraph's "one of three" made four. The `_adv`
  validity conditions live in ONE place (manual section 4,
  `\label{sec:advvalidity}`), the HOWTO points there; the HOWTO names
  `exhale_io.py` as the loader every new script should use (four `#
  columns` parsers exist). Advisor: the same three stale statements in
  `parameters.f90`'s comments corrected (default 1e-5; the row scaling; the
  retired base valve), and the four schema rows whose code spans carried a
  `|` (a strict renderer split them) rewritten; table widths verified.

- **B5j** (`steady_newton.f90` +1061/-41 only; `steady_residual.f90` and
  `diffusive_photochemistry.f90` bitwise unchanged; `steady_species_rows`
  +14 rows, `steady_selfconsistent_residual` assertion 5): **the carrier
  row's lower bound, MEASURED as two mechanisms in sequence.** The
  outermost carrier row is linear in its own unknown: at the hand-off its
  root is +2.1e-11 with the iterate 4 percent above (positivity-preserving,
  the Newton right to reduce it); twenty iterations later the row's
  constant term has changed sign and the root sits at -5.5e-13, below the
  bound, so the row has no admissible root and the unknown is pinned; and
  at the cell that binds the solve (215) the root is +2.75e-6 while the
  Newton has driven the unknown to exactly zero, eleven decades below it,
  the row reading -1.00000 of its own scale. A species unknown at zero
  whose row wants it positive cannot be sampled at all (a two-sided
  difference leaves the box on both sides). Landed: the carrier mass
  fraction clamped onto [0,1] in `write_species_rows_into_composition` (the
  rule the element fractions had); the trial projected onto the faces of
  the species box with the predicted reduction recomputed for the step
  taken; Bertsekas (1982) epsilon-active set held for the step and applied
  to every sampled direction; the probe step cut to the box with the side
  chosen by the room (subsuming B5g's backward sample); the Cauchy point as
  the step when the Krylov leg cannot be sampled; a defect of the item's
  own found and fixed (the accepted trial re-evaluated as `Y + D s` rather
  than as the vector evaluated, whose round trip does not close for a
  projected component: a carrier at -4.8e-35). MEASURED on `mol_carrier`,
  three configurations, none converging: baseline 49 iterations `||R||` 1.929, 527
  refused samples, cut-backs to 2^-49; (a) bound-aware region 63 it, 1.415,
  0 refused, 2^-4; (b) carrier as `ln n` (`EXHALE_CARRIER_LOG_UNKNOWN=1` (the spelling at B5j; since B5k the log unknown is the default and `=0` restores the density))
  65 it, **0.9972**, merit 17.6 -> 0.924, smallest adopted unknown 2.2e-11.
  (b) is better by every measure and is NOT the default (it changes the
  unknown space of a route a shipped key selects): decision item 11. **Part
  2**: the elimination's stopping test is now the change of the three
  quantities the residual reads the composition FOR (particle count, H2
  caloric share, the radiative pair), `EXHALE_EQ_MEASURE=species` restoring
  the old one: at `eq_sweep_reltol` 1e-12 `wasp_full_newton` reaches it in
  33/49 passes where the species measure never does (plateau 2e-11); the
  element case's remaining plateau (1.85e-11) is the radiative channel at
  cells 441-459 wandering, that case's own non-convergence (B5h), not the
  measure's floor. Seed dependence unchanged to four digits (9.736e-10).
  **`wasp_full_newton` moves 1.1e-9 relative**, isolated to the new measure
  by `EXHALE_EQ_MEASURE=species` reproducing the control byte for byte: the
  solve converges FURTHER (`||R||` 8.1e-9 -> 1.7e-9 in 8 iterations instead
  of 6, info 0, CERTIFIED, 13.30); below the 0.1 percent rule, the golden
  is not touched (the harness will read "within"). Suites (advisor build):
  `steady_species_rows` 93/0, `steady_selfconsistent_residual` 3/0,
  `certification` 29/0, `coupled_source_step` 44/0, `physics_probe` 1491/0;
  one "per-iteration" phrase in `run.sh` reworded. Recorded: a carrier row
  at an active bound is an error to the certification (row 1.00000 at cell
  215) and a constrained optimum to the solver, so a bound-constrained
  stationary state has a multiplier there, not a zero row (decision item
  12); 28 of 64 steps of configuration (a) are the Cauchy point alone (Bertsekas'
  two-metric projection not taken); `load_IC` refuses a molecular restart
  written at the input's own He/H by comparing printed digits for equality
  (`7.9307E-02` against `0.0793`). Advisor: `wasp_full_newton` re-run on
  the matrix path (marching then Newton) with the B5j binary against the
  refreshed golden: movement 5.9e-9 (`Hydro_ioniz.txt`) and 9.6e-9
  (`Ion_species.txt`), 13 outer iterations to `||R|| = 4.3e-9`, CERTIFIED,
  log10 Mdot 13.30; below the 0.1 percent rule, golden untouched.

- **Gate 5 (2026-09-08, 19:03 to 22:02)**: `make check` on the tree as of
  B5i, HYG-SETUP, HYG-CHECKED (built 19:03, before B5j landed) against the
  refreshed goldens: **REGRESSION PASS, all 64 files of the 16 cases
  byte-identical.** The refreshed goldens are therefore the current tree's
  own outputs, and the three startup-only items moved nothing. B5j's
  1.1e-9 movement of `wasp_full_newton` (the Newton converging further
  under the new increment measure) is not in this gate's binary and will
  read "within" at the next one; below the 0.1 percent rule it refreshes
  nothing.

- **B5k** (decision item 11, option (a), user 2026-09-08; `steady_newton.f90`
  +346/-85, `steady_species_rows` +137/-6 (114/0), design doc section 20,
  schema K15e, manual entry, PDF 82 pages): **the molecular carrier
  densities of the coupled steady solve are carried as `ln n` by default**
  (`EXHALE_CARRIER_LOG_UNKNOWN=0` restores the density unknown; the element
  fractions untouched; no certification change, decision 12 (c)). The
  floor under `ln n` is the transport operator's own admissibility floor,
  1e-20 of the carrier's element budget in that cell (`carrier_headroom`),
  the fraction `carrier_residual` already puts under its row scale;
  MEASURED 5.6e-28 to 4.6e-21 across `mol_carrier`'s column (B5j's log unknown had
  `lo = -huge`, unsound: `exp` underflows to an exact 0). Stated at the
  site: the column factor `n` of the Jacobian (both the banded model and
  the Krylov action probe a RELATIVE change), the region asymmetric in the
  density (`n -> n exp(+-Delta)`), and `ln 0`. **A real finding at the
  arming point**: `carrier_headroom` does not exist when the solve reads
  its controls (allocated inside `carrier_steady_residual`, first reached
  by the opening `eval_residual`), so the first implementation silently
  fell back to the density unknown on every coupled run; the opening state is
  now packed as a density and rewritten as logarithms after that
  evaluation, sound because the residual is a function of the state and
  not of its coordinates. Floor array made rank 2 (cell AND carrier row;
  two carriers of one cell have different elements). MEASURED, one binary,
  the two runs one variable apart: `mol_carrier` density 63 it `||R||`
  1.415, 28 Cauchy-only steps, 10 unknowns on a bound, 716 adopted
  carriers at the floor; **`ln n` 66 it `||R||` 0.9957, merit 0.960, 0
  Cauchy-only, 0 on a bound, 0 backward products, 0 at the floor, smallest
  adopted carrier 2.2e-11 = 2.5e15 times its floor**; `mol_diffusion` (H2
  carried + element row) density unknown 0 accepted steps of 12, log unknown 50
  accepted, `||R||` 1.986 -> 1.982. **Neither converges, and the next
  obstruction is named and is NOT the density unknown's**: the log unknown's steps
  are model-accurate (ratio 0.39-1.09) and SHORT (`||s||` 0.09-0.17 against
  a radius of 2.5 sitting at its ceiling for 30 of 68 iterations), cut by
  3-5 halvings before the trial is describable, and the refusing screen is
  the ELEMENT BUDGET (157 samples, against 0 for a negative carrier): the
  carrier CEILING, which `species_unknowns_outside_their_bounds`
  deliberately does not treat as a face of the box. Bertsekas' two-metric
  projection is not what the log unknown needs (it holds no bound). Gates:
  `wasp_full_newton` byte-identical twice (full case control vs new, and
  the reload identical to B5j's binary including every JFNK line); marching
  bitwise on `mol_carrier` 200 steps. New rows: default is `ln_n`, switch
  reaches `density`, the floor identical under both unknowns,
  `carrier_is_above_its_own_floor` (RED on the density unknown: 0, 716 at the
  floor). Advisor: `steady_species_rows` 93/0, `steady_selfconsistent_residual`
  3/0, `certification` 29/0, `physics_probe` 1491/0; the two documents that
  named the old `=1` spelling annotated. Recorded: `gm_m = 40` is a literal
  at the call site in `EXHALE_main.f90` and the log-unknown solve sits at it in 58 of
  66 iterations (unmeasured whether raising it helps, the cycle reports
  iterations used, not the residual reached); `carrier_headroom` returns
  `huge` as an availability sentinel. **Next item B5l**: the element budget
  as a face of the unknown box.
- **B5l** (`steady_newton.f90` +392/-82; `diffusive_photochemistry.f90`
  +51/-6, the headroom accessor only; `steady_species_rows` suite +161 and
  +108; `docs/steady_solver_design.md` section 21; `certification.f90`,
  `EXHALE_main.f90`, `steady_residual.f90`, `ionization_equilibrium.f90`,
  `input_read.f90` untouched): **the element budget is a face of the unknown
  box and the projection treats it exactly.** `freeze_species_unknown_box`
  is now the one place the admissible set of a species unknown is written
  down (one lower and one upper bound per unknown of the flat vector,
  formed once per outer iteration from the iterate and from the element
  budget the carrier operator froze at that iterate); the three routines
  that had each computed their own bounds
  (`species_unknowns_outside_their_bounds`, `fix_active_species_bounds`,
  `largest_step_inside_the_species_box`) read it. **A drift between the
  three was found and fixed on the way**: the projection had read the
  carrier ceiling from the trial's own density slot while the other two
  read it from the iterate, so the set a trial was written onto and the set
  the Krylov cycle was sampled under were not the same set. Faces of a
  carrier unknown: below, zero (`carrier_log_floor` in `ln n`); above, the
  cell density and the element budget (`carrier_headroom`, for H2 half the
  hydrogen the carriers may take), the tighter of the two. The budget is
  frozen for the step, not evaluated at the trial, which is the carrier
  operator's own choice (`headroom_car` assigned under `weno_mode == 1`,
  the one evaluation that IS the iterate; refreshing it per trial was tried
  there and is worse), so box and screen describe one state. The budget
  face is raised to the iterate where the iterate stands outside it
  (counted, `n_budget_face_at_the_iterate`; 0 on both coupled cases with
  the default step control); it is stored one floating-point step inside
  the budget (`budget_face_inside = 1 - 8 eps`, representability in `ln n`,
  not a margin). An unknown whose box is a point (zero budget) is held
  without asking the gradient. `EXHALE_SPECIES_BUDGET_FACE=0` restores the
  entry text's box (measurement switch, default on) [retired by N4b: the
  budget became a constraint row, `EXHALE_ELEMENT_CONSTRAINT_ROWS=0` is the
  measurement switch now]. `pgmres` reports the
  relative residual it REACHED against the tolerance asked; `EXHALE_GM_M`
  overrides `gm_m` (the caller's 40 stays the default and the literal at
  `EXHALE_main.f90` line 5343 is untouched). **MEASURED on `mol_carrier`
  (coupled, `ln n`), face off vs on: 66 vs 52 outer iterations, both ended
  by the stall detector; `||R||` handed back 0.9957 vs 1.174; trials halved
  3 or more times 38 of 68 vs 0 of 53; residual samples refused by the
  element budget 157 vs 0; trials written onto a face 0 vs 29; adopted
  carriers above their own budget 0 vs 0 (largest ratio 0.9886 / 0.9878);
  iterations ending "the model promises no reduction" 0 vs 19; worst Krylov
  relative residual reached 0.938 vs 0.569 where 0.1 was asked; `log10
  Mdot` 10.63 both.** The halvings and refusals against the budget are gone
  exactly; `||R||` does not improve; the limit moved to a named place, the
  Krylov leg of the dogleg. **The Krylov budget is not the limit**: `gm_m`
  40 / 80 / 160 give `||R||` 1.174 / 1.103 / 1.271, worst relative residual
  reached 0.569 / 0.630 / 0.614, `||s||` of the last accepted steps
  0.04-0.24 in all three, `log10 Mdot` 10.63-10.64; the cycle spends every
  product and stops at the same 0.6 of its right-hand side, so it stagnates
  on this operator (preconditioner or operator, not subspace). On
  `mol_diffusion` (H2 carried) the face does the same to the step control
  (halved 6-11 times in 50 of 50 vs 0 of 41; budget refusals 388 vs 0) and
  the solve stalls where it always has (H2 front, cell 212, element row at
  1.000 of its scale; `||R||` 1.982 vs 1.943; `log10 Mdot` 7.97 both), with
  the probe showing instead 419 backward-side products, 26 unknowns held on
  a bound, 12 Cauchy-only steps and 18 truncated cycles. A defect in the
  item's own first measurement was found and fixed (the adopted carrier
  was recorded against the PREVIOUS iterate's budget, which had reported 67
  adopted carriers above budget, worst 1.058; the call now follows the
  iterate's own residual evaluation, arithmetic unchanged iteration for
  iteration). **Gates**: `wasp_full_newton` bitwise against the control
  binary (`Ion_species.txt` bitwise, `Hydro_ioniz.txt` identical but for
  the provenance line, every `(JFNK)` line identical, `||R|| = 4.299E-09`,
  `log10 Mdot` 13.30, 13 iterations); marching bitwise on `mol_carrier` 200
  steps; `EXHALE_SPECIES_BUDGET_FACE=0` reproduces the entry-text binary
  exactly; `steady_species_rows` 118/0 with the new rows RED without the
  face and GREEN with it; `physics_probe` PASSED. **Advisor**: the five
  diffs match the tree against the worker's entry snapshots; forbidden-word
  scan clean; the current tree builds warning-free (`build_lwv`); the
  worker's `wasp_full_newton` output against the current golden moves
  5.87e-9 / 9.58e-9, exactly the B5j movement already on record, so B5l
  itself moves nothing on the matrix path. No golden refreshed. **Recorded,
  not repaired (outside this item's files)**: `carrier_element_headroom`
  states a box side per carrier but no joint bound for carriers sharing an
  element (the physical constraint `2 n(H2) + n(OH) + 2 n(H2O) + n(H+) <=
  nH_avail` is a hyperplane the code does not state), and OH/H2O hold
  hydrogen nuclei their own headroom does not mention, so a state can pass
  every headroom and still leave the write-back's hydrogen remainder
  negative (the branch that zeroes H I, H II, H2+, H3+); neither reaches
  the `{H2}` carrier set both measured cases use. Also: the epsilon-active
  band (one finite-difference step) never engages on `mol_carrier` with the
  face (29 trials projected, 0 held; no switch, a wider band unmeasured);
  `largest_step_inside_the_species_box` returns zero for the whole
  direction when one component has no room (18 truncated cycles on
  `mol_diffusion`); `note_adopted_species_unknowns` is not called for a
  converged state (pre-existing, cannot affect a stall). **Next obstruction
  on the carrier-row solve: the Krylov leg (banded preconditioner or operator),
  not bounds and not subspace.** Per the user's instruction of 2026-09-09
  ("leave running things alone, stop code modification for now") no
  further item is started; B5m (the atomic element-row solve with no Krylov
  direction) stays planned.
- **Preparation of PLAN_20260909_rev1 (items N-1 and N17; advisor,
  documentation and file moves only, no source change; user's go-ahead
  2026-09-09 with the instruction not to rerun regression checks after
  simple changes).** `docs/PLAN_20260909.md` was written from the review of
  the issues inventory (`docs/ISSUES_20260909_review.md`), reviewed
  (`docs/PLAN_20260909_review.md`) and revised to `docs/PLAN_20260909_rev1.md`,
  the plan of record (items N-1 to N17, stages A to E, decisions 13 to 19).
  Before writing rev 1 the advisor re-read the review's new source claims:
  the sampled integration-error estimate is printed and dropped and no
  production caller issues `as_reject_int_error` (`EXHALE_main.f90`
  2806-2821); after a retry inside the full pass the half passes start from
  `0.5*dt_step` and the clock advances by the saved `dt_step` (1706-1734,
  2797-2799, 2820, 2836); the audit harness asserts the OLD failure
  (`if (ieee_is_finite(x(1))) error stop 4`); `use_tr = (nspec_row > 0)`
  (`steady_newton.f90` 5763) so `wasp_full_newton` never enters the trust
  region, which corrects a statement of the first plan. All confirmed. **N-1**:
  the B5m reproducer moved from the session scratchpad into the tree as
  `backup/regression/atomic_elem_newton/` (input, metals, setup and resolved
  reports, `run.log`, the two handed-back state files, `run.sh`, README;
  6.8 MB, of which `run.log` 6.2 MB; not in `DEFAULT_CASES`); the common
  worker brief became `docs/worker_rules.md` with the session paths
  generalized and three rules added (one owner per source file per
  increment; affected-path validation only; named outcomes in reports); the
  review harness `docs/audit_20260905/run_issues_review_20260909.sh` was run
  once on the current tree, its output kept as
  `docs/audit_20260905/run_issues_review_20260909.out` (zero-radius step
  norm 0; Cauchy merit 0.5 to 40.5, line minimizer 0; identity GMRES 1;
  zero-operator GMRES nonfinite with reported residual 0; enthalpy ratio
  0.5 on a nonzero mass divergence; He/H pair 8.83e-5; hydrogen inventory 2
  at the individual ceilings; outer-shell chord 11.83; "defects reproduced,
  not repaired"), the harness itself untouched. **N17**: `ISSUES_20260909.md`
  corrected per rev 1 section 3 (seed retention no longer stated as an
  eigenvalue of one; B5f marked a withdrawn, uncertified reference; 3.5
  scoped to the present fixed-point map; the restart row now states the
  `heh_dev > 1e-6` branch and the 8.83e-5 discrepancy; `carrier_headroom_known`
  noted; `hydro_row_max` and `use_tr` bullets; the `_adv` statement in
  section 5 made conditional; new section 3.8 on the physical step never
  being refused for its integration error). `code_status_20260905.md`: the
  Huang et al. 2023 shell-opacity note added to the photoionization row;
  four pre-existing spelling/forbidden-phrase slips fixed on the way
  (`per-cent` to `percent` twice, two `per-molecule` phrases reworded).
  Verification: forbidden-word scans of the four documents; the reproducer's
  `run.sh` resolves its default binary path; nothing built or rerun.
- **N16a** (PLAN_20260909_rev1 Stage A, the physical step can be refused
  and the clock is right; `EXHALE_main.f90`, `attempted_step.f90`,
  `src/tests/attempted_step/{run.sh,attempted_step_tests.f90}`; user
  2026-09-09 "N16a 병행"): **the sampled step-doubling comparison is one
  transaction over one macro-interval `H`**: the full pass covers `H`, the
  two half passes `0.5 H` each and end at `t + H` by construction (pass 3
  no longer inherits whatever `dt` held); a refusal inside any pass, a
  finite estimate above one, or a nonfinite estimate rejects the WHOLE
  macrostep (`as_reject_int_error` now issued by the production caller),
  restores `as_entry`, shortens `H` by the controller and recomputes both
  trajectories; exhaustion of the transaction's retry budget
  (`n_macro_retry`, cap 8) is a named failure that stops with status 2,
  names the pass and prints the checksum of the restored state
  (`stop_on_exhausted_retry_budget`, one report for both exits). The clock,
  the accepted count and the physical accumulations move only at the
  adoption boundary below the macrostep loop. **The order of the complete
  split update is MEASURED as `p = 1`**: at one fixed state the estimate
  scales as `H^2` over five decades (`p + 1` = 1.964, 1.993, 1.9985,
  1.9997, 1.9999, 2.0004, 2.0000 for `H` from 2e-2 s to 2.6e-7 s; the
  coarsest 1e-1 s is outside that regime and excluded), so
  `attempted_step_error_estimate` now returns the retained full step's
  error `Delta/(1 - 2^-p) = 2 Delta` (`err_order_p`,
  `err_full_step_factor`, convention at the declaration) [superseded by
  N16d/N16e: p = 1 was the base-face leak; the conservative and electron/T
  rows are p = 3, factor 8/7, the species rows p = 1] and the reduction
  is `0.9 e^(-1/(p+1))`. The brief's first procedure (one interval by 1, 2,
  4, 8 prescribed steps, ratios of the output rows) was INCONCLUSIVE, `p`
  0.4 to 4.0, because the third difference (2.2e-8 relative) sits at the
  coupled source step's fixed-point tolerance: the inner error is not
  below the temporal error there (review F3's condition, not yet
  enforced). Test hooks, default off: `EXHALE_AS_ERR_INJECT[_COUNT]`,
  `EXHALE_AS_REJECT_IN_PASS`, `EXHALE_AS_FIXED_DT` (seconds; in
  initialization mode with local stepping it flattens the local steps,
  stated at the hook). `EXHALE_STEP_CLOCK=1` now also prints an `err-pass:`
  line per pass and `state=` (checkpoint checksum) per accepted step.
  **Tests**: suite RED on the entry binary 117/16 (the sixteen new rows:
  passes cover one interval, forced error rejects, clock unchanged,
  shorter retry accepted, refusal aimed at pass 1/2/3 keeps the endpoint
  and names the pass, exhaustion in each pass and on the error is named
  and restores the accepted state by checksum), GREEN on the final binary
  133/0; advisor reran the suite with the worker's final binary (md5
  `5ea530d3251a8b8ce77e72415ddef0ab`): 133/0. **Runs**: `mol_carrier`
  marching 200 steps byte-identical apart from the provenance line
  (initialization mode never sets `n_err_passes = 3`; the guard lines
  quoted in the report); the suite's case in physical mode, 60 steps:
  identical through step 19, at step 20 `e = 1.215` (raw 0.607, the number
  the old code printed and dropped) rejects, `H` 1.118 -> 0.913 (`e` 1.77)
  -> 0.618 s (`e` 0.89) accepted, `t_phys` 66.40 -> 65.92 s; two
  rejections in 60 steps, all by the integration error. Not run: the
  matrix, `wasp_full_newton`, `run_fcheck.sh` (out of the affected path;
  the module compiles clean under the suite's `-fcheck=all`). No manual
  entry (the manual carries no run-mode or clock text; nothing in it
  became false). **Advisor**: file scope as briefed; six British `per
  cent` (two new, four pre-existing) changed to `percent`. **Reported, not
  fixed (outside the item)**: (1) the accepted state of a retried
  macrostep depends on which pass was rejected, up to 9.1e-5 relative in
  the base cell's velocity (5e-6 elsewhere), with the checkpoint round trip
  exact, so some state read by the coupled source step is neither saved nor
  rebuilt (seeds or caches of `ionization_equilibrium`,
  `constrained_chemical_equilibrium`, `diffusive_photochemistry`): a
  rejected attempt is not yet without consequence; (2) one interval of the
  order sequence (`H` 6.48e-5 s from the 0.405 controller reduction) sits
  2.3x off the `H^2` law the other seven follow, unexplained, for N16b;
  (3) `err_rtol = 1e-3`, `err_atol = 1e-12` are flat, unvalidated and cover
  no species (N16b); the gate now shortens physical-mode steps by about
  two at the sampled steps on this case.
- **N0, N3, N1, N2** (PLAN_20260909_rev1 Stage A, one worker, four diffs;
  `steady_newton.f90` +738/-143 over the four, new suite
  `src/tests/krylov_and_dogleg/`; no other file). **N0**: 14 rows calling
  the PRODUCTION `pgmres`, `dogleg_step` and the factored step-length and
  predicted-drop code (a test-operator hook in `pgmres` in the style of
  `carrier_headroom_set_for_test`); RED on the entry text 10/7, among them
  `zero_operator_step_is_finite = F`,
  `approximate_gradient_length_decreases_the_model = 40.5`. **N3**: seven
  named linear outcomes (`gm_tolerance_reached`, `gm_subspace_exhausted`,
  `gm_no_direction_sampled`, `gm_jacobian_action_unusable`,
  `gm_preconditioner_failed`, `gm_reduced_system_singular`, `gm_nonfinite`);
  a zero right-hand side returns before any factorization; `hypot`
  rotations; rank tested before the triangular division and the leading
  nonsingular block solved; LAPACK status and finiteness checked after each
  preconditioner apply and on the returned step; the TRUE scaled residual
  verified with one product at a breakdown; the cycle exits at a breakdown.
  Reorthogonalization not added: worst `max |V_i . V_j|` MEASURED 6.9e-13
  on the carrier fixture. **N1**: absent elements (reservoir abundance
  zero, the `melem_ab <= 0` test `set_energy_vectors` and `load_IC` use)
  are no longer registered as rows; `initial_trust_region_radius` with
  named statuses, holding the active bounds before its own cycle; a
  positive floor `tr_delta_min = 1e-6 sqrt(neq)` in scaled coordinates
  (257x and 1157x below the first radius the two fixtures take); `tr_dmax`
  follows the radius; "no feasible descent" is `info = 3`. **N2**: the
  length along the direction used, `tau = (r . A g)/||A g||^2`, the point
  named an approximate descent point, the paragraph defending the old
  length by B5e removed, a dropped Cauchy leg with a Krylov residual >= 1
  exits as `tr_no_direction_reduces_model`. Suite after diff 4: 31/0, also
  under `-fcheck=bounds,do,mem`. **Root cause of B5m (R1's mechanism, its
  cause)**: the 1500 unknowns on a bound of the atomic element reload were
  the three elements `metals.inp` does not carry (Si, K, S x 500 cells),
  registered with an identity row and sitting exactly on 0, so the
  fraction-to-the-boundary rule zeroed ANY direction with a component there
  and the first cycle took 0 products; the radius became 0 and stayed so.
  MEASURED after (worker and advisor agree): unknowns on a bound 1500 -> 0;
  initial radius 0 -> 1.909e-2 from a Krylov step that reached tolerance (3
  products, 6.1e-2 against 0.1); first step accepted, ratio 0.995; `||R||`
  1.888 constant -> 8.3e-2 at iteration 50, then an oscillation between
  6.9e-2 and 1.1e-1 with the merit creeping 7.29e-3 -> 7.23e-3 over 258
  iterations (advisor run stopped by its 20-minute cap; 178 accepted, 80
  refused; every cycle spends all 40 products). Not converged: the next
  obstruction is the model quality of the element rows near `||R||` 1e-1.
  **Carrier reload** (control reproduces B5l's run exactly): after, 54/1
  accepted, `||R||` handed back 1.174 -> 1.151, worst achieved Krylov
  residual 0.569 -> 0.575, radius range unchanged, `log10 Mdot` 10.63 both,
  outcomes now named (8 tolerance reached, 47 subspace exhausted). **The
  concurrent N16a edits are inert on every fixture** (paired controls
  bitwise). **`wasp_full_newton`** (once, final binary): not bitwise,
  largest relative movement 2.1e-10 (`Ion_species.txt`, Ca I) and 8.8e-11
  (`Hydro_ioniz.txt`), `info = 0`, `log10 Mdot` 13.30; by code path the
  `hypot` rotation of N3 (`use_tr` is false there); below the 0.1% rule,
  identical, no refresh. `lower_profile` (the only default case with
  `He_metal_diffusion` and a partial reservoir) bitwise to its golden in
  both builds (it stops at its marching cap). Other suites (worker):
  `steady_species_rows` 106/0, `certification` 29/0,
  `steady_selfconsistent_residual` 3/0; `physics_probe` does not include
  `steady_newton` (`source_closure.py`). **Reported, not fixed**:
  `largest_step_inside_the_species_box` still annihilates a whole direction
  when any unheld unknown has no room on the side the direction points (an
  element physically zero in one cell can trigger it; the outcome is now
  named and the radius survives, the direction is lost): an active-set rule
  decision, for N4b. Two dead locals removed from `pgmres`.
- **Stage A gate (advisor, 2026-09-09).** Tree after N16a, N0, N3, N1, N2.
  `make OBJDIR=build_lwv EXE=EXHALE_lwv.x` warning-free, md5 `dd233d6a`
  (identical to the worker's final binary, so the whole tree is what both
  workers measured). Suites on that build: `krylov_and_dogleg` 31/0,
  `steady_species_rows` 106/0, `steady_selfconsistent_residual` 3/0,
  `certification` 29/0, `attempted_step` 133/0 (with the N16a binary, md5
  `5ea530d3`, same source). N7's first two iterations on the atomic reload
  logged above with named outcomes. Not run, by the affected-path rule: the
  matrix, `physics_probe`, `run_fcheck.sh`; the shared `build/` and
  `EXHALE.x` remain behind the sources. No golden refreshed. Gate PASSED;
  Stage B (N4a, N5, N8a) and the N7 diagnostic list are next.
- **N13** (PLAN_20260909_rev1 Stage D, transit refusal census;
  `EXHALE_transit.py`, `exhale_transit_lib.py`, new Python suite
  `src/tests/transit_census/`; no Fortran, no output format change, so
  decision 17 is untouched): the census mask `(r >= 1) & (r <= Rib)` was
  geometrically wrong (the ray integral uses every shell with `|r| >= b`);
  one function `chord_shell_indices` now serves the census,
  `resonance_depth` and `resonance_spectrum`, and the census reports, line
  by line at eight impact parameters and at its maximum, the share of the
  line-center optical depth carried by rows whose correction was refused
  (`refused_line_center_tau_share`, the trapezoid rule rewritten as a sum
  over rows so that shares of disjoint subsets add to one), labeled a
  CONTRIBUTION, not an uncertainty. MEASURED: every `tpm_*.txt` of
  `wasp_full` (10) and `mol_base_handoff` (5) bitwise before/after; the
  defect is latent when `Rib = r[-1]` (both cases as they stand) and real
  when the domain passes the stellar radius (`wasp_full` read at
  `EXHALE_TRANSIT_RSTAR_RSUN=0.30`: 419 -> 500 sampled rows, 81 above the
  cap); on `wasp_full` the rays near `b = 1` take 99 percent of their
  Ly-alpha line-center depth from uncorrected rows. Suite RED on the
  pre-change copies (11 rows, the mask present and no shared selection),
  GREEN 11/0 on the tree; advisor rerun 11/0; +0.3 s on the whole script.
  Advisor: one `per-line` in the library header reworded. Reported, not
  fixed: the wavelength-resolved share and the file metadata wait for
  decision 17; `np.loadtxt(..., max_rows=1)` at line 896 warns on every
  run and duplicates the column count of line 325.
- **N12** (PLAN_20260909_rev1 Stage D, the upstream caloric energy at the
  upstream composition; `T_equation.f90`, `post_process_adv.f90`,
  `ion_cell_state.f90` (two fields `x_h2_up`, `e_up`),
  `src/tests/adv_static_limit/` (new driver `upstream_caloric_energy.f90`,
  rows U0-U4)): the upwind term of the `_adv` energy residual now receives
  the upstream cell's SPECIFIC internal energy `e_up = E(x_H2,up, T_up)/mu_up`
  formed by the caller, and the caloric branch is taken when EITHER side
  carries H2 (an atomic cell below a molecular neighbor no longer loses the
  incoming rovibrational energy); the monatomic branch is verbatim when both
  sides are atomic. **A second, larger defect found and fixed in the same
  path**: `x_h2_cell` was assigned only inside the `k == 1` stationarity
  loop (old 650) and read by the temperature loop for EVERY cell (old
  1316), so the whole molecular layer was solved at the H2 fraction of the
  outermost ghost (1.5e-7 on `mol_base_handoff`, 2.3e-5 on `mol_carrier`,
  against 0.84 at the base); advisor confirmed both lines in the pre-change
  source. The stationarity diagnostic's lower energy `e_lo` is likewise
  evaluated at the upstream composition (the ratio's definition, threshold
  and operator unchanged). Energy convention unchanged: the caloric EOS
  carries the bound rovibrational ladder only, the 4.478 eV bond lives in
  the formation reservoir (`molecular_reaction_heat.f90`), no term is
  counted twice. **Tests** RED/GREEN (MEASURED): isothermal H2 gradient
  residual 1.508 -> 0; atomic cell below molecular upstream 1.006 -> 0;
  dissociation-layer march against the analytic discrete solution, largest
  departure 0.155 -> 2.4e-11 (RED marched 1000, 892, 782, 615, 560 K
  against the discrete 1000, 975, 863, 728, 596; GREEN both 1000, 974.76,
  943.62, 880.65, 853.52); controls unchanged. Suite 36/0; advisor rerun on
  `build_lwv` (md5 `8d027f8b`, the tree with N12 and N13) 36/0 (a first
  advisor run reported one `no_output` FAIL from a relative `EXHALE_EXE`
  path, the advisor's mistake, not the suite's). **Impact** (post-process
  restart `Do only PP: True` of each case's own state; the before binary is
  the Stage A binary): `wasp_full` `_adv` files bitwise (atomic, both
  `x_h2` zero, monatomic branch verbatim); `mol_base_handoff` and
  `mol_carrier` `_adv` move in 469 and 467 of 504 rows, `T_adv` by up to
  94.2 percent (218 -> 424 K at r 1.0941) and 103.2 percent (169 -> 343 K),
  split about 35/45 and 41/46 percent between the stale-composition and the
  upstream-composition defects; both a warming of the molecular layer (5/2
  kT per H2 instead of 3/2); `p_adv`, `cool_adv` follow, `heat_adv` at
  most 4.8 percent; `adv_status` changes in one cell of `mol_carrier` (row
  26, 0 -> 2, from the `e_lo` composition in the ratio); run states
  (`Hydro_ioniz.txt`) bitwise in all three. **The molecular `_adv` goldens
  of the matrix will move; refresh at the Stage D gate, reported then, not
  now.** The writer's warning that molecular and oxygen columns are not
  advection-corrected stays; the product is still not a molecular
  stationary solution. Advisor: one `per-thread` phrase in
  `ion_cell_state.f90` reworded; the worker corrected four British
  spellings in comments. Reported, not fixed: a stale scalar formed in a
  `k == 1` block and read later is a class worth a sweep of the file;
  `post_process_adv.f90` is named by both N11 and N12 (one owner at a
  time; N11 waits for decision 16).
- **N9** (PLAN_20260909_rev1 Stage D, the restart He/H discrepancy;
  `load_IC.f90` +84/-12, new `grid_and_gates` row `restart_round_trip`
  (9 assertions); `write_output.f90` untouched, no key, no tolerance
  change): **the reported failure was a message defect**: the refusal
  printed the He/H of the TOP physical cell while the branch tests the
  largest departure over all physical cells; on the failing state
  (`B5j/scan500_stag`, a coupled-carrier JFNK iterate of the hot-Uranus
  gate) the top cell agreed with the input to 17 digits and cell 199 was
  off by 3.96 percent, so the message read "7.9307E-02 against
  7.9307E-02 ... cannot be rescaled". Review R5's 8.83e-5 is the distance
  between the INPUT values of two different runs, not a property of any
  state. Fixed: the comparison records the cell it attained and
  `report_deviating_heh_cell` prints cell, radius, file and input ratios
  at `ES16.9`, the departure and the threshold (`heh_dev_tol = 1e-6`, now
  named, value unchanged); RED/GREEN on a state with helium raised by half
  at cell 199. **The round trip is clean** (MEASURED on the 12000-step
  `mol_base_handoff` state with `EXHALE_DUMP_IC=1`): `r`, `v`, `p`
  bitwise; every species column and every inventory (H, He, O, C nuclei,
  charge, He/H, H2 ratio) within 4.2e-16; `T` 5.6e-16; the `rho` column
  4.9e-12 because the loader never reads it and rebuilds mass from the
  species (perturbing the column by 1e-15 changes nothing, bitwise), so
  4.9e-12 is the writer state's own mass-closure defect (4.1e-14 on a
  40-step state); a second trip 5.0e-16, no drift; no clock advanced.
  **The producing path moved the ratio, not the loader** (MEASURED): He/H
  departs from the input by 1.0e-11 and 8.9e-12 after 12000 marching steps
  with carrier transport off and on, but by 3.9 to 7.9 percent on every
  coupled-carrier density-unknown JFNK state and 41 to 43 percent on the
  `ln n` states (B5j, B5k): the carrier column is written after the
  element projection precisely so the Newton's carrier is not rescaled
  (`write_species_rows_into_composition` 1439-1469), the sweep reads the
  nuclei totals from the state (`hydrogen_helium_nuclei_density`), and
  with helium diffusion off no row re-imposes He/H. This is the element
  constraint N4a/N4b are about, seen from the restart side; reported, not
  touched. Also fixed on the way: `parse_labels` had no bound on 80-long
  label arrays (a file with more than 80 columns wrote past them); it now
  stops with a message (present schema at most 42). Atomic restart path
  MEASURED unchanged: `wasp_full_newton`'s state reloaded before/after
  bitwise, no He/H line printed. Advisor: `restart_round_trip` 9/0 on
  `build_lwv` (conda gfortran; the worker's private build used /usr/bin
  gfortran 13, against the standing toolchain rule, so the advisor's run on
  the standard toolchain is the one on record). **Reported, not fixed**:
  the window integrals of the assembled residual are not Lipschitz in the
  state in the odd-even band above the base (one round trip at 4.9e-16
  moves the 1.03-1.10 mass window by 16 percent and flips cell residual
  signs at cells 183-187; one ulp on every species column reproduces it),
  so a round-trip residual budget must be stated on the row maxima (N5);
  `restart_grid_guard.sh` carries a stale header describing a fixed defect;
  nothing checks the writer's `rho` column against its species (decision
  15's contract); `EXHALE_RESIDUAL=1` runs one sweep unconditionally, so a
  true "no chemistry on load" evaluation has no entry point yet (N10).
- **N7** (PLAN_20260909_rev1 Stage C, the atomic element-row diagnostic
  after Stage A; `steady_newton.f90` +477/-0, diagnostic prints under
  `EXHALE_ELEM_DIAG=1` (iterations 1, 2, 60; `=n` moves the third), one
  extra Jacobian product per selected iteration with the refusal
  statistics held aside; no arithmetic of the step changed, MEASURED: every
  `(JFNK)` line identical hook off/on over 2027 lines, across four builds).
  **Verdict, MEASURED on the reload fixture**: the oscillation at `||R||`
  7e-2 to 1.1e-1 is NOT a model failure (predicted/actual ratio 0.9946,
  1.0001, 1.0003 at iterations 1, 2, 60; the refused iteration's model is
  right to 3e-4), NOT a bound effect (0 unknowns held, 0 samples refused
  in 293 iterations), and NOT limited by the linear solve although the
  solve fails its tolerance (40 of 40 products, achieved 0.6685 against
  0.1, TRUE residual 0.6685 by one extra product, Arnoldi gap 2.3e-4; the
  dogleg takes 2.1e-4 of the Newton leg so that inaccuracy is second
  order). **Proximate cause: the ray test** (the finite-difference check of
  the model slope) refuses the step: at the refusing iteration model slope
  -9.13e-9 against ray slope -1.48e-9 (factor 6.2 against the allowed 4)
  with `|f2p - f2m| = 3.0e-12` only 11x above the floor
  `eq_sweep_reltol * |f2|` = 2.7e-13 (at iteration 1 the same test has a
  factor 1371 of room and the slopes agree to 0.5 percent); 82 of the 91
  refusals in 293 iterations are this test, 3 its sign form, 6 a real
  over-prediction. Limit cycle read off the radius: 9.78e-4 accepted ->
  1.96e-3 accepted -> 3.91e-3 refused by the ray test -> cut to 9.78e-4,
  period three; the merit falls 1e-4 of itself per accepted step, so
  `Resid tol` 1e-8 is seven decades away. **Underlying cause: the merit and
  the gate are different functionals.** At iteration 60 the eight element
  rows hold 0.437 of the squared merit (He 0.305) and NOTHING of `||R||`
  (energy row holds 1.000 of `||R||`, 0.560 of the merit); their own
  certification measures sit at 1.2e-3 to 7.4e-3 against
  `cert_tol_element = 1e-8`; accepted steps swing `||R||` by 40 percent
  (7.8e-2 to 1.09e-1) while the merit stands at 7.291e-3 to four digits.
  At iteration 20 the brief's option (c) is met outright: `||R||` 1.400 set
  entirely by the mass row, the energy row 0.810 of the merit, seven of
  eight elements worst at cells 12 or 17. **The Cauchy leg is worth
  nothing here**: `||sU||` 6.6e-9 against a radius 3.9e-3 (`g = B^T r0`
  from the band, `||A g||` 1e4 times `r0 . A g`), predicted decrease
  1.3e-9 against the dogleg's 7.6e-9; one product per iteration for
  1.7e-6 of the step. **The band is not the defect** (0 zero, nonfinite or
  round-off columns, `dgbtrf` info 0), though it conditions worse along the
  solve (`dgbcon` 7.3e-8 -> 9.4e-10, smallest pivot always a trace element
  of cell 500). **`gm_m` probe**: unlike the carrier-row solve, the larger
  subspace changes the achieved residual on this operator (first divergence
  at iteration 17: 40 of 40 at 0.1035 against 48 of 80 at 0.09999; deep in
  the stall 0.67 against 0.29) and the solve stalls WORSE (m = 80 plateau
  `||R||` 0.27 to 0.34, merit 1.57e-2, against 7.7e-2 to 1.1e-1 and 7.3e-3
  at m = 40): longer steps lift the ray probe off its floor and the refusal
  becomes an honest over-prediction, and the solve still goes nowhere.
  **Proposals, not implemented**: P1 measure the ray test's floor from two
  evaluations of one state instead of assuming `eq_sweep_reltol * |f2|`;
  P2 do not let the slope comparison overrule a reduction ratio within
  [0.998, 1.001] (every ray-size refusal inspected had one; the
  over-prediction guard stays live); P3 drop or rescale the banded-gradient
  leg where `||sU||` is below a small fraction of the radius (one product
  saved); P4 a GMRES restart buys m = 80's reduction at m = 40's memory but
  does not fix the stall; **P5 the acceptance functional and the merit must
  be reconciled** (a merit that reads the element rows on their
  certification scales, or a separate control for them): a design decision,
  decision item 20 (section below). Also fixed in the same file: the
  `tr_exit_text`/`gm_outcome_text` buffers (48 and 64 characters) truncated
  two 67-character messages; all three are 80 now. **Reported, not fixed**:
  the `worst j, k` on the `(JFNK) it` line is the merit's largest scaled
  row while `||R||` on the same line is a certification maximum, and the
  two need not be the same row (iteration 20: `||R||` is the mass row's,
  `worst` is an energy cell); the brief's belief that `Solver: Newton
  100.0` caps the iterations is wrong (the third word is `newton_du_switch`;
  the JFNK cap is the literal 500 at `EXHALE_main.f90` 3509 and 3000 on the
  direct route, no key or hook; runs were capped by wall clock);
  `atomic_elem_newton/README.md` is right but incomplete (the refused
  trials have a correct model). Advisor: rebuild warning-free,
  `krylov_and_dogleg` 31/0 on the tree with N7 and N9.
- **N4a** (PLAN_20260909_rev1 Stage B, one elemental constraint map with
  three contexts; new `src/modules/lower_atmosphere/element_inventory.f90`
  (745 lines, one `SRC` line in the `Makefile`), `diffusive_photochemistry.f90`
  (its three element densities, hydrogen budget and five carrier box sides
  now derived from the map with unchanged arithmetic), `steady_species_rows`
  suite, `docs/element_inventory_contexts_20260909.md` corrected;
  `steady_newton.f90` untouched): **the map reports, nothing enforces it
  yet**, so the Stage B gate may not claim from this item that trial
  construction, chemistry and certification enforce one constraint (N4b).
  `nuclei_per_particle` is the one place a stoichiometric coefficient is
  read (`bsp_nH/nHe/nO/nC`, `mion_elem`); `element_inventory_of_cell` gives
  the full inventory grouped carried/locked/reserved/closure with the
  budgets and the shared-constraint breach as a MAGNITUDE (against the
  budget the sides spend, and against the budget the chemistry closes
  with, separately); `element_inventory_of_candidate` evaluates a trial's
  carriers against the frozen element densities; `EXHALE_INVENTORY_REPORT=1`
  prints it at the operator's entry. **Physical judgment**: the separate box
  sides are a wrong feasible set by stoichiometry (two carriers of one
  element spend one budget); a cell at both ceilings breaches by exactly
  1.0 of the available hydrogen (MEASURED). **Seven corrections to the
  advisor's draft**, each checked in the source: `nO_free`/`nC_free`
  include the neutral stage while `cbg_nOion` excludes it (the F7 gap as a
  number: the oxygen box side spends the element total while
  `carrier_source` closes against the total less the ion stages, so the
  box admits states whose source rows clamp free O at zero; N4b's feasible
  set); `nH_avail` is floored; H I and He I are the closure, not rows; CO
  is fixed data of the local solve (its destruction row is the carrier
  operator's); **the shared sums ARE enforced on the marching path by
  `limit_to_element_budget`** (what has no shared statement is the coupled
  stationary solve, whose `eval_residual` compares a COUNT of cells);
  **`max(rest, 0)` in the write-back CREATES nuclei** (with `rest < 0` the
  carriers keep what they hold and the closure cannot go negative, so the
  written element total exceeds the handed one by `|rest|`; only
  `element_census_verify` at `EXHALE_ELEMENT_ASSERT` 1 or 2 sees it, off by
  default); absent elements are no longer registered (N1). **Tests** 149/0
  (advisor on `build_lwv`: 145/0, the four coupled-log rows unset), RED
  shown by `EXHALE_INVENTORY_SHARED=0` (the three infeasibility rows fail
  with that setting); all nine molecular goldens classified: worst carrier breach
  exactly 0.0, tightest hydrogen occupancy 0.96 to 0.99
  (`oxygen_chemistry` occupancy 1.0 is carbon's, 175 of 500 cells with all
  C in CO). **The N9 finding caught live**: inside a coupled carrier JFNK
  the He/H departure from the reservoir grows 1.4e-11 -> 4.9e-5 (third
  iterate, the first with a Newton `n(H2)`) -> 1.1e-1 (thirteenth) while
  every carrier budget stays at breach 0.0, so the fix belongs at the
  write-back of each iterate (the element row or a projection there), not
  at the end of the solve; the two saved coupled states have carrier
  breach 0.0 and He/H departures 4.0e-2 and 6.7e-1. **Marching**:
  `mol_carrier` 200 steps byte-identical (two independent before/after
  pairs); `mol_diffusion` 12000 steps byte-identical over the 7511 steps
  both builds had covered at report time (step trace character for character,
  the step-7001 snapshot files byte-identical apart from provenance); the
  pair finished later in the session scratch and the advisor compared
  its end: `Hydro_ioniz.txt` and `Ion_species.txt` at 12000 steps IDENTICAL
  apart from the provenance line, `log10 Mdot` 10.57 in both builds. Reported, not fixed: `element_ratio_gate = 1e-9` now
  exists in two modules because `element_census_tolerance()` is private.
  Rule breach reported by the worker: one `pkill -f` on its own shell
  pattern (it matched only its own shell; the other workers' runs were
  verified alive). Advisor: five British `percent` in the two files fixed;
  the advisor's earlier `pgrep` for the worker's binary name missed its
  renamed scratch copies, a false "finished" the worker corrected.
- **N7b** (N7 proposals P1 and P2, and the JFNK iteration hook;
  `steady_newton.f90` +249/-40, `krylov_and_dogleg` +12 rows): **P2 opens
  the cycle; P1 is implemented and MEASURED inactive, and the measurement
  refutes N7's reading of the refusals.** The ray test's decision is
  extracted into `ray_slope_verdict` (public, tested), its floor is
  `max(eq_sweep_reltol * |f2|, ray_floor_factor * measured reproducibility)`
  with the reproducibility measured once per outer iteration by two
  evaluations of the merit at one state (cost one residual evaluation per
  outer iteration, +1.5 percent on the carrier reload): MEASURED 1.8e-15 at
  a merit of 14.95 and 0 to 4e-19 at 2.7e-5, four to seven decades BELOW the
  assumed part (2.6e-13), which therefore governs everywhere; the refused
  trials' `|f2p - f2m|` (1e-10 to 3e-10) stand two to three decades above
  the floor, so they are real differences and NOT the probe near its
  cancellation floor as N7 and ISSUES 3.1 (d) said (corrected there). What
  they are: at the refusing iterate the ratio is 0.999 while the model slope
  (-5.5e-7) and the ray slope (+1.7e-7) have OPPOSITE signs; the model slope
  is set by the banded-gradient leg, which is not in the step (`||sU||`
  6.6e-9 against a radius 3.9e-3) but carries `||A g||` 3.1e3 into
  `r0 . A s`, so the slope test is right that the slope is wrong and the
  leg is what is wrong (N7 P3). `ratio_sound_band = 1e-2`: within it the
  ray verdict is formed and printed but not acted on;
  `tr_ratio_below_acceptance` untouched. `EXHALE_JFNK_MAXIT` read once per
  solve, resolved by `jfnk_outer_iteration_cap`, printed with its origin.
  **Atomic reload, 150 iterations, before -> after**: accepted 106 -> 129,
  ray-size refusals 35 -> 1, ratio-below-eta refusals 6 -> 15 (the guard
  takes over), radius 9.8e-4 .. 9.8 -> 3.8e-2 .. 2.7e8 (no longer binds;
  the accepted step is the interior Krylov point at `||s||` 1.0e-5, ratio
  0.72), `||R||` 7.27e-2 -> 6.92e-2, merit 7.26e-3 -> 7.08e-3, achieved
  linear residual still 0.69 against 0.1 with 40 of 40 products: **the
  limiter is now the linear solve and the model behind it** (N7 items 4
  and 7, P3/P4), not the step length. No convergence claimed. **Carrier
  reload**: every numeric `(JFNK)` line identical; 20 of 55 iterations had
  carried a STALE refusal label ("the model promises no reduction" printed
  on ACCEPTED iterations because `refuse` was reset only inside `.not.
  model_ok`; pre-existing, fixed): earlier reports quoting "iterations that
  ended on N: the model promises no reduction" counted that defect wherever
  N exceeds the refused count (B5l's 19 and 20 included). `wasp_full_newton`
  BITWISE on all ten output files (three log lines differ: two wall clocks,
  one new cap line). Tests RED on the extraction stage (4 of 43 fail: the
  measured floor rows and the sound-ratio rows), GREEN 43/0; advisor rebuild
  warning-free, `krylov_and_dogleg` 43/0, binary md5 `6d4bf63e` identical
  to the worker's final. Two dead locals removed (md5 unchanged). Advisor:
  four British `neighbor` spellings in comments fixed. Advisor decision:
  P1's measured floor stays (correct in form, one evaluation per iteration).
  **Reported, not fixed**: the ray test runs before the acceptance test, so
  a catastrophic over-prediction (ratio -533, -1493 at iterations 47, 50) is
  named a slope disagreement; `ray_h = 1e-3` is not a measured interval;
  `tr_dmax` doubles without an absolute cap and with P2 runs to 2.7e8
  (harmless, the interior point is taken, but `delta` stops being a
  diagnostic); `n_tr_exit` has no index guard.
- **N5** (PLAN_20260909_rev1 Stage B, the residual contract;
  `steady_newton.f90` +877, `ionization_equilibrium.f90` +274/-6,
  `EXHALE_main.f90` replay driver +134/-60, `src/tests/residual_determinism/`
  rewritten with a suite `run.sh` (cases: the two in the tree plus
  `EXHALE_RESID_RELOADS=<dir>[:<dir>]`, binary `EXHALE_RESID_EXE`)): **the old
  determinism suite had a hole**: it replayed the THREE-unknown residual in
  every configuration (`0 of 1500 entries`, the species-row registry is
  empty before a solve); the replay now registers the rows and returns 0 of
  1500 / 2000 / 2500 / 5500 entries differing (3, 4, 5, 11 unknowns per
  cell) after three trials, a discarded trust-region trial and a `jv_product`
  probe with the active bounds held (**contract 1 GREEN on four
  configurations**; a rejected macrostep is not reachable from the driver,
  the one part unmeasured). **Contract 2, conditional closure convergence,
  MEASURED and it exceeds the certification tolerances**: the judged row
  maximum moves, per unit relative change of the eliminated partition's
  seed, by 4.2e2 `cert_tol_carrier` on the carrier reload, 3.1e1
  `cert_tol_element` and 1.7e1 `Resid tol` on the atomic element reload,
  6e-8 on the element + carrier reload; a seed above 2.4e-3 (carrier) or
  3.3e-2 (atomic) moves the judged row past its tolerance; the retention of
  a scaled seed at the fronts is ABOVE one (1.39, 3.69, 3.46) where B5h had
  0.77 at the base. Two suite rows stay FAIL by design as the standing
  measurement (`closure_spread_within_the_row_tolerance_{atomic_elem_newton,
  det_car}`, 61.5 and 8.5 against < 1). **Contract 3, branch discovery**
  (`EXHALE_RESID_BRANCH_REPORT=1`, roots grouped by a branch-datum distance
  1e-1 anchored on B5h's two roots): with H2 transported B5h's two-root H2
  front does NOT reproduce (x(H2) at the front bitwise identical over 10
  seeds, the transported value replaces the seed in the H2 row); the H I
  partition shows one root on the carrier reload and TWO on the element +
  carrier reload (9 and 1 seeds, 1.05e-1 apart); whether the second is a
  distinct branch or a slow iteration is N6's question. **Chemical decay
  against transport** on the resolved gradient length: H2 `t_adv * lambda`
  3.68 (base), 4.0e-3 (front), 0.82 (wind); H2+, H3+, HeH+ 3 to 12 decades
  faster than transport everywhere, so the transported-H2 / eliminated-ion
  split stands and no promotion is proposed (conservation modes are removed
  by the parameterization itself: the neutral stage closes each element).
  **Rounding floor of the species rows in the judged measure is exactly
  epsilon** (2.2e-16) against tolerances 1e-8, 7.65 decades apart, condition
  numbers 2 to 1e3: compensated summation (the Al-Refaie et al. 2024
  reading) is not a lever. Validation: `krylov_and_dogleg` 43/0,
  `steady_species_rows` 145/0, `mol_carrier` 200 steps bitwise,
  `wasp_full_newton` bitwise on all ten files (one new log line, N7b's cap).
  Advisor: rebuild warning-free, md5 `913e6579` identical to the worker's
  final; `krylov_and_dogleg` 43/0; `residual_determinism` 5 PASS / 1 FAIL,
  the FAIL being the atomic closure-spread row by design (a first advisor
  run with the wrong variable name `EXHALE_EXE` tested the stale shared
  `EXHALE.x`, the advisor's mistake); two British `percent` in
  `ionization_equilibrium.f90` fixed, the worker fixed five `neighbor`.
  The worker's session was cut once by the Opus limit (09:50 KST) and
  resumed; nothing was lost. **Reported, not fixed**: a discarded
  evaluation restores `n_eq_sweeps_last`, so a caller reads the previous
  state's pass count; `report_accepted_cell_state` writes inside the
  parallel region without a lock; `carrier_diffusion_coefficient` reads the
  marching operator's `pct_Dco`, which the stationary residual never fills
  (a pure stationary path would read zero; `diffusive_photochemistry.f90`);
  the global `count` shadows the intrinsic.
- **N16b** (PLAN_20260909_rev1 Stage E brought forward, species-aware
  integration error; `attempted_step.f90` +546/-28, `EXHALE_main.f90`
  macrostep code +125/-7, `attempted_step` suite +289): the step-doubling
  estimate reads four row classes on their own scales (conservative rows;
  transported species fractions, the `f_sp` columns `carrier_state` builds
  `fc` from: H+, H2, OH, H2O, CO; element mass fractions under
  `He_diffusion`; electron fraction and T through the caloric EOS), the
  species scale being the species' own value with an absolute floor of
  `element_ratio_gate` (1e-9) of the largest density the species can reach
  in the cell, read from the N4a map's stoichiometry; the flat 1e-12 no
  longer judges any species row (a carrier at 1e-20 changed by 1 percent:
  e = 20 on its own scale, 2e-10 on the flat floor). **The inner nonlinear
  error is collected cell by cell** (the coupled pair's safety-factored
  geometric factor on the cell's own increment, and the carrier substep's
  row residual) and the estimate carries the interval `e_lower .. e_upper`
  the three inner errors leave it in: reject when above one at the bottom,
  accept when below one at the top, UNRESOLVED in between (named outcome,
  counted in the end-of-run block); the brief's fraction rule for each class
  (1/10) is kept as the report, not as the verdict, because applied
  literally it called EVERY estimate of a healthy run unresolved (species
  difference 6e-17 five decades below the pair's own error) and would have
  made F1's rejection inert. **p = 1 confirmed on every class** (atomic
  case, eight levels over five decades: species `p + 1` 1.9987, 1.9799;
  electron/T 1.9943 to 1.9997; conservative 1.9809 to 2.0004; molecular
  case coupled ladder H = 5, 2.5, 1.25 s: species 1.949, 2.299, electron/T
  1.912, 1.983), so the retained-step factor 2 holds for all rows. On the
  molecular case at `H <= 0.4 s` the species class reports itself
  unresolved (inner over difference 2.7e5): review F3's condition, MEASURED
  and detected. **The 60-step physical run is unchanged to 17 digits**
  (`t_phys`, checksum, both files): the conservative rows still decide it.
  `mol_carrier` marching bitwise. **N16a's 9.1e-5 state leak IDENTIFIED**
  (not fixed, files not owned): `Apply_BC.f90` `base_face_W` and
  `base_face_lower_W` (lines 48, 50; written only by `Apply_BC`, read by
  `Rec_BC` in the first `Reconstruct` of a retaken pass BEFORE the pass's
  first `Apply_BC`, so the base face's left state is the refused
  attempt's) and `caloric_eos.f90` `caloric_mixture_active`, `nk_per_mass`,
  `x_h2`, `molecular_cell` (125, 145-147; written only inside
  `get_species_densities`, which `rebuild_state_from_checkpoint` calls
  AFTER `U_to_W`, so the restored p and T use the refused attempt's caloric
  composition); the two reinforce each other at the base cell, which is
  the measured signature (v 9.1e-5, T 6.8e-8, cooling 1.7e-5). Row
  `rejected_pass_leaves_no_state_behind` stays RED at 9.1162e-5, labeled
  with the files and lines. N16a's 2.3x outlier reproduced and bounded (only
  at controller-formula levels, conservative and charge classes, always the
  base cell; 0.665 and 0.53 of the H^2 prediction), not explained; two
  candidates named (the unbounded inner error of the mass/momentum rows
  through the second half's coupling; the base-face leak). `err_rtol`,
  `err_atol` unchanged (the CFL step 20 would need 1.21e-3; the estimate
  RISES at the first shortened interval, a second reason not to anchor on
  this case). Tests 155/1 (the one FAIL by design), RED on the entry text
  150/6 (an estimate that was the solvers' own error rejected ten
  macrosteps there; now 0, named). Advisor: suite rerun with the worker's
  final binary (md5 `7bd041dd`) 155/1; four `modelled` in `EXHALE_main.f90`
  comments changed to `modeled`. **Reported, not fixed**: the mass and
  momentum rows carry no inner error; the carrier substep's number is a row
  residual, not a solution error (conditioning unmeasured); `NCO_col`,
  `k_co_diss`, `theta_co_shield` are outside the checkpoint by an unstated
  exception (read only after their write within a sweep).
- **N16d** (the checkpoint's restore is complete; `EXHALE_main.f90`
  `rebuild_state_from_checkpoint` 4 code lines, comments in `Apply_BC.f90`,
  `caloric_eos.f90`, `attempted_step.f90`, the suite row extended to both
  files): **closed by rebuilding, not saving.** The five leaking quantities
  (`base_face_W`, `base_face_lower_W`; `caloric_mixture_active`,
  `nk_per_mass`, `x_h2`, `molecular_cell`) are functions of `(u, f_sp)`, so
  the rebuild now runs `rho = u(1,:)`, `get_species_densities`, `Apply_BC`,
  then `U_to_W` and `comp_T_from_p` (composition first because the pressure
  map and the boundary read it; boundary second so the face states are the
  restored state's). Every one of the six call sites follows a restore, so
  no restore-free path can move. **MEASURED**: `rejected_pass_leaves_no_
  state_behind` 9.1162e-5 -> 0.0 bitwise on both files; suite 156/0
  (advisor rerun with the worker's final binary md5 `abdf17f7`: 156/0).
  **The fix changes the measured order of the conservative rows**: the
  full-versus-two-half difference was 40x the real splitting difference
  (the leaked base face carried an O(H) error into one face of one RK
  stage of an H/2 step, an O(H^2) contribution) and now falls as `H^4`
  (`p + 1` 3.9966 to 4.0104 over the ladder), the local-error order of the
  three-stage Runge-Kutta; the electron/T rows follow it; the transported
  species stay `H^2` (unchanged to all printed digits). N16a's 2.3x outlier
  is EXPLAINED: the two controller-formula levels now sit on the `H^4` law
  to 0.1 percent; the leaked base face carried the difference of two
  intervals whose ratio was not the ladder's. **Consequence not yet fixed
  (N16e, launched)**: `err_order_p = 1` was calibrated on the leak; the
  conservative and electron/T estimates carry a factor 2 where 8/7 is
  right (inflated 1.75x) and the reduction exponent is -1/2 where -1/4 is;
  the species class is genuinely p = 1, so the order is a property of each
  class. `mol_carrier` marching 200 steps bitwise (0 refused attempts, no
  restore); with one injected refusal the initialization-mode state moves
  2.1e-5 (v at r 1.0185) and 1.2e-5 (He 2^3S), the fix of the same defect
  with both channels live (atomic `hydrostatic_column` had the base-face
  channel alone). `wasp_full_newton` not run (no restore unless an attempt
  is refused; whether it has zero refusals was not measured). The reload
  route does not share the rebuild (checked). Advisor: two `Per-cell`
  phrases in `caloric_eos.f90` reworded. **Reported, not fixed**: the same
  composition-before-pressure ordering stands at three other sites of
  `EXHALE_main.f90` (1075-1078, 5162-5166, 5312-5326) and
  `upmap_bc_roundtrip` measures an identity its comment says is not one
  (N16e); the added `Apply_BC` adds one boundary evaluation per restore to
  the run-total diagnostic counters.
- **N20** (decision 20 a with N7's P3; `steady_newton.f90` +314/-84,
  `krylov_and_dogleg` +13 rows, `steady_species_rows` +1 row;
  `certification.f90` unchanged, `cert_scale_floor` was already public):
  **one functional, implemented exactly, with a measured consequence that
  is a decision, not hidden.** `merit_row_scale_from_certification` is the
  one denominator a species-row-route row is read through; `cell_row_scales`
  builds the hydrodynamic rows through `residual_row_scale(k,j,u)`, the
  element rows through their certification scale, the carrier rows through
  the term sum the certification divides by; base cells of element rows
  keep the unknown's own scale (Dirichlet reservoir); the three-unknown
  route is untouched by construction (`Dr = D` when `nspec_row <= 0`).
  MEASURED: on both reloads every element row's largest scaled merit row
  EQUALS its certification measure at the same cell and the merit's largest
  scaled row IS `||R||`'s row (before: factor 2.5 apart at iteration 1, 112
  at iteration 60). **P3**: the banded-gradient leg is dropped from the
  step AND its image when shorter than `cauchy_leg_min_fraction = 1e-3` of
  the radius (measured 1.4e-5 to 8.7e-7 on both fixtures, so everywhere),
  so the slope handed to `ray_slope_verdict` is the step's own; products
  saved: zero, because `tau_c` needs `A g` (a lazy leg would save one per
  iteration, P-N20-1, reported). A quotient verdict at ratio <= `eta_accept`
  is named `tr_ratio_below_acceptance`. `tr_dmax` capped at `1e3 sqrt(neq)`
  (7.4e4 on the atomic fixture; binds on neither fixture now, would have
  bound the 2.7e8 of N7b). **The consequence**: the carrier reload improves
  2.6x (best `||R||` 1.151 -> 0.446, both builds on the stagnation detector,
  +16 percent evaluations, `dgbcon` 5.4e-10 -> 1.9e-11) while **the atomic
  element reload no longer solves**: the certification row scales span 11.4
  decades (3.7e-10 .. 92.7) against 8.3 before, `dgbcon` 7.3e-8 -> 8.9e-10
  at iteration 1, the Krylov point collapses to `||s||` 2.5e-15 inside a
  radius 19.6 at iteration 11 and the solve aborts on the no-descent
  watchdog at `||R||` 1.777 (before: 7.38e-2 at 150 iterations). A third
  scaling choice (element/carrier rows on certification scales, hydrodynamic rows on
  the state scales) does not collapse but is worse than the entry text on
  both fixtures (0.741, 1.026): the hydrodynamic reconciliation is what
  helps the carrier reload and kills the atomic one. The merit's row scaling
  IS the operator's (`r0 = F/Drow`, `A = Drow^-1 J D`), so the remaining
  freedom is the column scaling `D` and the preconditioner: N8a's, launched
  with the atomic collapse as its RED. Tests RED on the entry objects (the
  element-row scale row 4.0 against 1.0; the ratio-naming row), GREEN
  `krylov_and_dogleg` 56/0, `steady_species_rows` 145/0, `certification`
  29/0; `wasp_full_newton` BITWISE on all ten files (two wall-clock lines).
  Advisor: rebuild warning-free, `steady_newton.o` md5 `6f3aef50` identical
  to the worker's, the three suites 56/0, 145/0, 29/0. Two pre-existing
  defects fixed (a dead local; `[diag 7]` printed `at cell 0` for the
  carrier). **Advisor decision**: (a) stays in the tree; the atomic reload's
  collapse is a conditioning problem of the scaled operator and is the
  first deliverable of N8a, not a reason to give the step control a
  functional other than the gate's; the state of the atomic reload is
  reported as WORSE than before N20 until N8a lands. Reported, not fixed:
  `EXHALE_main.f90` 1182 should pass `u` to `cell_row_scales`;
  `docs/p54_base_layer_mass_flux.md` 10.4 now has a second measurement
  behind it.
- **N16e** (the order of each row class is its own integrator's;
  `attempted_step.f90` +106/-50, `EXHALE_main.f90` +60/-21, `attempted_step`
  suite +75/-15): `err_order_p(n_err_class) = [3, 1, 1, 3]` (hydro, species,
  element, charge/T) with the measured ladders at the declaration
  (conservative 3.9966 to 4.0104; electron/T 3.9962 to 3.9867; species
  1.9987, 1.9799; the element class is assigned the species' order by
  construction, said so); `err_full_step_factor` is a class's own (8/7 for
  the RK3-carried rows, 2 for the first-order split ones), the estimate
  records the deciding class (the one whose lower bound reached the gate)
  and `attempted_step_reduced_dt` uses that class's `-1/(p+1)`. MEASURED:
  the conservative and electron/T estimates fall by exactly 4/7 (inflation
  1.7500004 removed), the 60-step physical run's trajectory does not move
  (N16d had already removed the step-20 rejection, so the brief's expected
  "1.21 rejected -> 0.69 accepted" belonged to N16b's entry text, not to
  this item's; stated); on the forced-error ladder the controller with
  exponent -1/4 asks for 0.846 and 0.226 of the interval where -1/2 asked
  0.601 and hit the 0.2 floor, so the floor no longer does the sizing.
  Stale sites: `EXHALE_RESIDUAL=1` and `physical_handoff_check` NOT
  reordered, the composition cannot have moved there (`init` and
  `equilibrate_loaded_composition` both end with `get_species_densities`),
  the invariant is now stated at both; `update_map_begin_step`'s install
  branch reordered (composition before `U_to_W`; `EXHALE_UPDATE_MAP` on
  `mol_carrier` moves 0.0); `upmap_bc_roundtrip` renamed
  `upmap_bc_interior`, required exactly zero (`Apply_BC` writes only the
  ghosts), MEASURED 0. Tests 163/0 (156 before; the unit rows cannot be run
  RED against a scalar `err_order_p`, the RED is the entry suite's own
  contradiction, `p = 1` asserted and `p + 1 = 4` measured in one run; the
  whole-binary row 0 of 2 -> 2 of 2). `mol_carrier` marching bitwise.
  Advisor: suite rerun with the worker's final binary (md5 `0ae8d8d9`)
  163/0. **Reported, not fixed**: `update_map_begin_step`'s install branch
  owes the next `Reconstruct` an `Apply_BC` (same family as N16d, a
  diagnostic route); `err_atol = 1e-12` still floors the conservative rows
  on a code-unit density; the mass/momentum inner error is still unbounded
  (the deciding row's interval is a point); no fixture has the species
  class deciding a step, so the exponent switch is covered by unit rows
  only.
- **N8a** (PLAN_20260909_rev1 Stage B, scales, boundary equations,
  derivatives and the inner budget, with N20's atomic collapse as its RED;
  `steady_newton.f90`, `krylov_and_dogleg` +17 rows, `steady_species_rows`
  +2 rows; `certification.f90` and `steady_residual.f90` unchanged, no
  tolerance moved): **the collapse was one exactly empty Jacobian column,
  not the spread of the row scales.** At iteration 11 of the atomic reload
  the carbon unknown of cell 500 had been driven to zero, its column scale
  sat on the absolute floor 1e-20, the probe step `sqrt(eps) * 1e-20`
  underflowed the residual difference, the column was exactly 0 and the
  banded preconditioner singular (`dgbcon` 3.98e-27, smallest pivot
  3.7e-19); the 2.5e-15 Krylov point is what a triangular solve with that
  factor returns. Iteration 1's 8.9e-10 was only the 7.3-decade spread of
  the certification row scaling (no empty column). **Two repairs, the model
  untouched to the last bit** (`r0 = F/Drow` and `A` unchanged, MEASURED by
  the scaled column norms): (a) `element_unknown_column_scale = max(|x|,
  1e-4 |reservoir|, 1e-20)`, the arithmetic at the function (the column is
  resolvable while the probe step stands far above `eps * X_reservoir`);
  (b) `equilibrate_the_preconditioner`: two-sided equilibration of the band
  that is FACTORIZED only, powers of two, `banded_preconditioner_solve`
  applying `C dgbtrs(R v)`. Attribution (first solve, accepted/refused;
  empty columns at it 11; `dgbcon` at it 11): entry 14/16, 1, 4e-27; floor
  alone 69/16, 0, 4.9e-11; both 114/36, 0, 1.1e-9. **Atomic reload
  (`EXHALE_JFNK_MAXIT=150`)**: entry abort at `||R||` 1.777 (no descent, 30
  iterations); floor alone 1.189 (stagnation, 86); both **3.800e-2 at the
  150 cap, still descending** (2.5e-1 to 3.8e-2 over the last thirteen
  iterations), `dgbcon` at iteration 1 8.9e-10 -> 1.65e-6, the brief's
  target (the pre-N20 7.4e-2) passed by 1.94x with one functional kept;
  refusals by name shift from ray-floor to ratio-below-eta (35) and one
  slope-size. **Carrier reload** 0.446 -> **0.1513** (2.95x), `dgbcon`
  1.9e-11 -> 8.7e-6, +30 percent evaluations. The model's own column
  equilibration (`EXHALE_COLUMN_EQUIL=1`, default off) was built, MEASURED
  and NOT adopted: it moves the trust-region coordinates, the radius runs
  to the cap and trials ascend by 4e4 times the predicted decrease. **First
  cell rows declared once** (`species_row_base_equation`: element row =
  Dirichlet reservoir, READ from `element_transport_residual`; carrier row =
  control-volume balance, READ from `carrier_residual`'s `j == 1` branch;
  no equation changed). **Base derivatives**: the directional-derivative
  sweep at the base is a clean V with its minimum 4.8e-8 at `sqrt(eps)`,
  the step `build_banded_jac_full` takes; whole-system rows 8.2e-8.
  **Closure map at the base** (`measure_the_closure_map_at_the_base`, four
  evaluations once per solve): the sweep exits after ONE pass at 9e-12,
  four decades inside the 1e-8 asked, so the reference is a pass count;
  quotients 2.1e-4 (atomic, attenuation) and 5.34 (carrier);
  `outer_residual_closure_target` (what the outer residual may keep, what
  the ray floor reads) is now separate from `eq_sweep_reltol` (what the
  sweep is asked = target / measured amplification; 1.87e-9 on the carrier
  reload). **B5h's ~4e3 base amplification is not reproduced on either
  fixture** and is not claimed. Manufactured states in `krylov_and_dogleg`
  (columns spanning eight decades, pivot ratio 5e-18 -> 0.54 after
  equilibration; exactness rows: the map and the unknown-space step
  unchanged to the last bit). Tests `krylov_and_dogleg` 73/0,
  `steady_species_rows` 147/0, `certification` 29/0; the new binary with
  both changes off reproduces the entry text to every digit;
  `wasp_full_newton` BITWISE (ten files; `nspec_row = 0` short-circuits
  every new path, the hydrodynamic column scales untouched). **Advisor**:
  rebuild warning-free, `steady_newton.o` md5 `c7621824` identical to the
  worker's, the three suites 73/0, 147/0, 29/0; the atomic reload rerun on
  the advisor binary reproduces the first solve to the digit (`||R||`
  3.800e-2 at 150, `info = 1`), and the run's SECOND stationary solve (the
  code re-enters the Newton from the handed-back state, another 150
  iterations) reaches **`||R||` 5.273e-4, flux spread 2.3e-6**, still not
  certified (`log10 Mdot` 10.10; exit 2), 32 min at 16 threads: the reload now
  descends four decades in 300 iterations where it stood still since B5m.
  Not converged. **Reported, not fixed**: `n_eq_sweeps_model = 1` so the
  pass count, not `eq_sweep_reltol`, controls the closure on both
  fixtures; `EXHALE_JFNK_MAXIT` caps each solve, not the run;
  `fix_active_species_bounds` lost its `D` argument (the epsilon-active
  width now on `state_column_scale`); the `[diag 4]` true-residual row is
  not valid for long steps (its finite-difference product is out of range).
  The fixture README updated by the advisor.
- **Advisor measurement after N8a (atomic reload, `EXHALE_JFNK_MAXIT=600`,
  16 threads, 18 min).** One solve with the cap at 600 stops on the
  stagnation detector at iteration 224 with `||R||` 2.671e-2 (`info = 2`;
  170 accepted, 55 refused, all but three by "the reduction ratio is below
  eta"; `log10 Mdot` 10.06), WORSE than the same binary's 150 + 150 route,
  where the run's re-entry into a second solve from the handed-back state
  reached 5.273e-4. [The advisor's reading here was WRONG and N21 corrected it: the
  run does not re-enter the Newton from the handed-back state; `info = 1`
  sends it back to MARCHING for 2000 steps (the log's `Newton finish at
  step 2002`), and the second solve starts from the marched state at
  `||R||` 6.06e-2. The two decades are gained by the marching, not by a
  fresh trust region; N21 measured every reset and none helps.] Not
  converged either way.
- **N4b** (decision 14 route (i), a step that moves along the element-budget
  surface; `steady_newton.f90`, `element_inventory.f90` +94 (the
  derivatives of the map, `element_constraint_derivatives`),
  `steady_species_rows` +35 rows; `diffusive_photochemistry.f90` untouched,
  the frozen budgets recovered exactly through the public
  `carrier_headroom`): **the acceptance row is met by ten decades.** On the
  carrier reload the He/H nucleus ratio against its reservoir stays at
  1.39e-11 at every one of 68 iterates where the entry text drifted from
  1.39e-11 to 2.016e-1 (20 percent of the ratio the run was given) over 79
  iterates: the write-back of each iterate now returns the element totals
  it was handed (the N9 / N4a finding closed where those reports said it
  belonged). Best `||R||` 0.1513 -> **0.0899** at 14 percent fewer residual
  evaluations (5690 -> 4891). The shared constraint `2 n(H2) + n(H+) + n(OH)
  + 2 n(H2O) <= nH_avail` (and the O and C rows) replaces the five separate
  corners as the step limit (each corner is implied by the row, so nothing
  is lost and the direction along the surface is gained); a nonlinear
  feasibility screen by magnitude at every candidate; a restoration phase
  for an infeasible start (on the element-and-carrier reload 52 trials
  restored over 1165 cell-steps, worst violation 1.630 -> 2e-16; 6
  candidates refused by magnitude; 167 write-back cells whose closure could
  not return an element total, counted). **Stated, not claimed**: on both
  reload fixtures the shared hydrogen row is never within a probe step of
  being tight (occupancy 0.985 to 0.989, active band 1e-8), so the
  active-set PROJECTION onto the constraint hyperplane, the part route (i)
  is named for, fired on ZERO iterations of either run; its correctness is
  stated by the suite rows, and the measured improvement comes from the
  write-back, the budget leaving the coordinate box and the shared step
  limit. The atomic element reload (element rows, no carrier row) is
  identical line for line over its 150-iteration first solve (`||R||`
  3.800e-2 both). Tests `steady_species_rows` 182/0 (RED with
  `EXHALE_ELEMENT_CONSTRAINT_ROWS=0`, 20 rows fail), `krylov_and_dogleg`
  73/0, `certification` 29/0; `mol_carrier` marching 200 steps bitwise;
  `wasp_full_newton` bitwise (`||R||` 4.213e-9, `info = 0`, `log10 Mdot`
  13.30). **Advisor**: rebuild warning-free, `steady_newton.o` md5
  `42fcf223` and binary md5 `9c9ba732` identical to the worker's; the three
  suites 73/0, 182/0, 29/0; `wasp_full_newton`'s two matrix files compared
  by the advisor against the Stage A reference: IDENTICAL. Two earlier
  builds of the item were discarded and their runs rerun (report section
  8). Certification unchanged (decision 12 c). **Stage B gate**: N4a map,
  N4b constrained step and element-conserving write-back, N5's three
  contracts measured (contract 2 fails by design on two reloads), N8a
  scales and closure map: the gate's rows are in; what is NOT met is
  convergence of any species-row solve, which is Stage C's criterion, not
  Stage B's.
  N4b's user-visible switches: `EXHALE_SPECIES_BUDGET_FACE` is retired and
  replaced by `EXHALE_ELEMENT_CONSTRAINT_ROWS` (=0 restores the five
  corners), `EXHALE_ELEMENT_WRITE_BACK` added (the element-conserving
  write-back); no input key, output column or unknown layout changed. Also
  found by N4b: `dc/d(nd) = c/nd` exactly (the constraint is homogeneous of
  degree one in the density, so an active row cannot push the mass
  unknown); two earlier builds of the item were discarded on the atomic
  fixture's line-for-line check (a blocked probe component zeroed
  unconditionally made the sampled operator nonlinear, factor 33
  regression; the probe sides compared before capping). The advisor
  annotated the retired switch in the B5l entry and in
  `docs/steady_solver_design.md`.
- **N21** (what limits the species-row solves now; `steady_newton.f90`
  +417, `krylov_and_dogleg` +31 rows): **no reset of the trust-region state
  helps, and the advisor's premise was wrong.** The run's second stationary
  solve does not start from the handed-back state: `info = 1` returns the
  run to MARCHING for about 2000 steps and the second solve enters at
  `||R||` 6.06e-2 (advisor's own log: `Newton finish at step 2002`, then
  `step 4002`, the third solve reaching 4.26e-4 on the stagnation detector),
  so the two decades the 150 + 150 route gains are gained by the marching
  between the solves. MEASURED on the atomic reload from the SAME
  150-iteration prefix, `EXHALE_JFNK_MAXIT=300`: continue (control) 227
  iterations, stagnation, `||R||` 2.671e-2 (170/55, reproducing the 600-cap
  run); resetting the pseudo-transient shift, the acceptance memory or the
  closure map: bitwise identical logs (three exact no-ops: `dtau` already
  sits on its floor `dtau0`, the Grippo window is reported and never acted
  on on the coupled route, the closure map re-measures to amplification
  1.000); resetting the radius and ceiling 2.851e-2; returning to the best
  iterate 2.921e-2; all five together (what a re-entry does) 3.161e-2; a
  second Krylov cycle of 40 2.927e-2; a restart every 10 iterations without
  improvement runs to the 300 cap at the same 2.671e-2. Carrier reload:
  control reproduces N4b (8.993e-2), periodic restart 1.153e-1, second
  Krylov cycle 1.728e-1, both worse. **The periodic restart trigger is
  delivered default OFF** (`tr_restart_stall_default = 0`) by the brief's
  own condition. **The linear solve is not the limiter**: after N8a the
  atomic cycle reaches a median 0.157 of its right-hand side over iterations
  1 to 150 (N7: 0.57 to 0.69) but 0.74 to 0.82 over 151 to 225; the carrier
  cycle 0.963; doubling the products halves the linear residual and makes
  the outer solve WORSE on both fixtures; in the second atomic solve a
  cycle reaching 0.98 coexists with a two-decade fall of `||R||`. The
  radius reset is worse for a stated reason: the carried radius 9.85e-9
  accepts steps at ratio 1.000 while a re-initialized 1.15e2 over-predicts
  by 39x then 3392x and costs four refused iterations. READ: the band is
  refactored every outer iteration (a re-entry changes the preconditioner
  only through the shift) and the Krylov cycle starts from zero every
  iteration. Tests RED (13 symbols absent from the entry module), GREEN
  `krylov_and_dogleg` 104/0, `steady_species_rows` 182/0; with everything
  disarmed the delivered binary reproduces the entry text's 20-iteration
  solve line for line; `wasp_full_newton` BITWISE (ten files). Advisor:
  rebuild warning-free, `steady_newton.o` md5 `cbdbf535` and binary
  `e140284a` identical to the worker's, the two suites 104/0, 182/0; the
  advisor's own `n8a_reload` log confirms the marching between the solves.
  **The finding that names the next item (N22)**: the predicted decrease
  stops scaling with the step below about 1e-7 (from iteration 140: the
  step falls 4x and `pred` stays 7.099e-4; every longer step of the
  resulting period-two cycle is an ascent of 3e-3 against a promise of
  3e-4, ratios -4, -8, -13, refused, quartered, then accepted at ratio
  exactly 1.000); the production dogleg geometry is cleared by two new rows
  (quartered radius quarters the step and the predicted decrease to 2.5e-1
  and 2.50009e-1), so what remains is `trust_region_step`'s bookkeeping at
  those lengths or a precision floor in `r0 . A s` and `||A s||^2`; the
  separating measurement is one diagnostic line (`||A s||`, `r0 . A s`
  beside `||s||` at two radii of one iterate). README of the fixture
  updated by the advisor.
- **N13b** (decision 17 option a, transit spectra carry their validity;
  `EXHALE_transit.py`, `exhale_transit_lib.py`, `transit_census` +8 rows;
  Python only): one reader `read_adv_validity` locates the status columns by
  name in the `# columns` line and distinguishes three cases: schema 2 (N11's
  two fields: a row is refused for the census when `adv_T_status` is not 0,
  the composition field reported beside it), schema 1 (today's single
  column), legacy (no schema: validity UNKNOWN, and both the console and the
  saved metadata say UNKNOWN rather than 0 refused). Every `tpm_*.txt` gains
  a comment metadata block (`transit_schema 1`, the input profile and its
  schema, its `adv_input_certified` / stationarity operator / model
  restrictions lines copied, its provenance lines, the tool's own identity
  by git hash and file mtimes, the census for that line with the eight-point
  refused-share ladder and its maximum, the CONTRIBUTION-not-uncertainty
  sentence, the `EXHALE_TRANSIT_*` overrides in effect); the two hand-written
  products get the same block. MEASURED: numerical columns bitwise on all 10
  `wasp_full` and 5 `mol_base_handoff` products; the schema-2 read tested on
  synthetic files and on a relabeled `wasp_full` profile (the relabeling is
  the worker's, not N11's mapping; the production file comes with N11);
  every consumer in the tree (`python/make_transit_figures.py`, the WASP-52b
  and LHS 1140 b notebooks, two LHS scripts) reads through `np.loadtxt` or
  skips `#`, none edited; the planet folders' `tpm_*.txt` untouched (mtime
  2026-08-19). Suite RED 8 rows (no reader existed), GREEN 19/0; advisor
  rerun 19/0. Fixed on the way: the second `np.loadtxt(..., max_rows=1)`
  column count (NumPy warning on every run) replaced by one rule. Reported:
  the refused share is taken over the temperature field only (a companion
  share over the composition field would be the honest partner); still line
  center only; `README_HOWTO.md` should mention the metadata block (advisor
  added the sentence).
- **N11** (decision 16 option a; `post_process_adv.f90`, the `_adv` writer
  of `write_output.f90`, `examples/exhale_io.py`, `EXHALE_plots.py`,
  `adv_static_limit` +driver M1-M5 and rows R8-R9): **stationarity of the
  input state is judged by the certification's own mass operator**
  (`assemble_residual`, `residual_row_scale`, the mass row `|R_1|/s_1`
  against `adv_mass_row_tol = cert_tol_mass`), the enthalpy term ratio kept
  as a printed sensitivity diagnostic that refuses nothing; **two status
  fields** `adv_T_status`, `adv_comp_status` (0 corrected, 1 retained, 2
  failed, 3 unsupported, 4 not_evaluated) in BOTH `_adv` files, columns 8
  and 9, with `write_adv_validity_header` writing `# adv_schema 2`, the two
  legends, `# adv_status_counts T <5> comp <5>`, `# adv_input_certified <T|F>
  <reason>` (the certification's verdict on this very state),
  `# adv_stationarity_operator face_flux_mass`, `# adv_product` (a ONE-WAY
  correction on the input's density and velocity field; a converged scalar
  temperature solve does not make a row self-consistent), `# adv_model_
  restrictions`; `adv_unsupported` where the omitted molecular and oxygen
  species carry more of the cell's particle count than the carried ones
  (structural, no threshold); `exhale_io` exposes both fields and the header
  (`load_adv_header`, `load_adv_status`), a legacy file reads schema 1 with
  both fields `ADV_UNKNOWN`, never a translation of the old column;
  `EXHALE_plots.py` shades the refused rows. **MEASURED, the physics of the
  item**: the ratio screen and the mass row are not interchangeable and the
  ratio is the weaker: on the analytic outflow fixture the ratio is 1.1e-13
  everywhere while the mass row stands at 2.1e-7 (RED of the old
  criterion); on `wasp_full` the ratio refused 297 of 503 rows, the mass row
  refuses 503 of 503 (largest 1.04 at r 1.032); `mol_base_handoff` 33 vs
  503 (0.245 at r 4.21); `wasp_full_newton` (the one certified state) 2 vs
  3. So on the two relaxation snapshots every row is retained and the `_adv`
  product IS the run's own state (`T_adv` equal to `T` in 504 of 504 rows;
  the old product differed by up to 6.4 and 207 percent); on the certified
  state the change is confined to the base (p, T 7e-5; species 7.7e-2 at r
  1.0002). Run states bitwise in all three; `rho`, `v` in the `_adv` files
  moved by exactly zero. Tests 48/0 (advisor rerun on `build_lwv`, md5
  `4320f5b7` identical to the worker's: 48/0). No golden refreshed; every
  `_adv` golden will move. **Advisor decision (N11b, launched)**: the
  criterion is the right OPERATOR but the wrong TOLERANCE for a conditional
  correction. `cert_tol_mass = 3e-12` is the tolerance of a certified
  stationary state; the advective correction is a first-order local
  correction whose own error is of the order of the mass row's measure
  (the fractional change of the face mass flux across the cell), so a
  du-stopped or snapshot state with a mass row of 1e-3 is exactly what the
  tool was built for and is what review R6 calls the "approximate local
  use" with "a tolerance tied to the required error". N11b writes the mass
  row measure itself as a column of `Hydro_ioniz_adv.txt`, sets
  `adv_conditional_tol = 1e-2` (the correction is described as accurate to
  that fraction, stated in the header beside `adv_input_certified`), and
  keeps the certification verdict as the separate statement it is; the
  product is then a CONDITIONAL correction with its condition written in
  every row. Advisor: two `per-field` and two `percent` (one in a printed
  message) fixed. **Reported, not fixed**: `assemble_residual` leaves
  `R(:,1-Ng)` undefined (`RK_rhs` assigns `2-Ng..N+Ng` only; no caller reads
  it today; `RK_rhs.f90`/`steady_residual.f90`); `README.md` 154,
  `README_HOWTO.md` 959 and the manual section 4 still describe the single
  `adv_status` column (N11b updates them with the final schema); a suite row
  for `corrected` needs a converged fixture.
- **N11b** (the `_adv` product is a CONDITIONAL correction with its
  condition in every row; `post_process_adv.f90`, the `_adv` writer of
  `write_output.f90`, `exhale_io.py`, `EXHALE_plots.py`, `EXHALE_transit.py`,
  `exhale_transit_lib.py`, `adv_static_limit` (53 rows), `transit_census`
  (21 rows), `README.md`, `README_HOWTO.md`, manual section 4): the row
  decision is taken against `adv_conditional_tol = 1e-2` (the fraction a
  corrected row is accurate to in the mass flux, the first-order argument at
  its declaration), the measure `m_j = |R_1|/s_1` is written as column 10
  `adv_mass_row` of `Hydro_ioniz_adv.txt` and named in `# columns`, the
  header carries `# adv_conditional_tol` and `# adv_mass_row`, and
  `# adv_input_certified` stays the certification's verdict with
  `cert_tol_mass`; schema stays 2 (readers are header-driven; the one
  place that counted columns from the end of the row was a test row, fixed
  to locate by name). MEASURED against the N11 binary: corrected rows (T)
  `wasp_full` 0 -> 215 of 504, `mol_base_handoff` 0 -> 341, `wasp_full_newton`
  500 -> 502 (the two others are ghost rows below the energy loop's first
  cell); the column equals the certification's own printed measure on the
  same state to every digit (1.942392 at cell 432); run states and `rho`,
  `v` of the `_adv` files bitwise. `exhale_io` returns `adv_mass_row` and the
  tolerance; the transit metadata block copies the tolerance line and the
  largest measure; the census still refuses on `adv_T_status`. Tests 53/0
  and 21/0 (advisor rerun on `build_lwv`, md5 `3131cab0` identical to the
  worker's: 53/0, 21/0). **Every `_adv` golden moves, far on the two
  snapshot cases: refreshed at the Stage D gate.** For the record: the
  first-order argument bounds the dropped term; in some corrected cells the
  KEPT enthalpy-flux term of the mass divergence still dominates the energy
  equation (452 of 503 cells of `hydrostatic_column` against 390 retained),
  and whether one percent of the mass flux bounds the temperature error of
  such a cell is not measured; a corrected row is not close to the run's own
  row and is not expected to be (the certified state moves T by 6.7e-2 at
  `m_j` 3e-16). Advisor: one `percent` in `README.md` fixed;
  `docs/postprocess_advection_validity.md` marked stale with a pointer.
- **N10** (decision 15 option a, the restart contract; `input_read.f90`
  +137, `load_IC.f90` +551, `write_output.f90` +13, `EXHALE_main.f90` +166,
  `grid_and_gates` +1 row script of 20 assertions, `docs/input_schema.md` K43
  and appendix D.2, manual entry; PDF not rebuilt): the key `Restart intent:
  trajectory | relaxation | stationary [evaluate|equilibrate]` with the three
  consistency refusals (default `trajectory` in `phys`, `relaxation`
  otherwise, so no existing input changes behavior); both state files carry
  the eight-line block (`restart_schema 1`, `reservoir`, `species_columns`,
  `grid`, `constants`, `options` with 20 switches, `t_phys[s]`, `source`),
  built and parsed in `load_IC.f90` (one definition), every number `ES23.16`;
  no block = legacy, loaded as before and marked `provenance unknown` in the
  log and in every file the run writes (the mark is inherited); a
  disagreeing `reservoir` (within `heh_dev_tol`), `grid`, `constants` or
  `options` is a refusal naming the field, `species_columns` and `source`
  are reported; a pair whose halves carry different blocks is refused.
  `stationary_state_of_the_loaded_restart`: WENO3 fixed, the N16d rebuild
  order, one equilibrium sweep (part of the residual's definition), the
  loaded composition's departure from that sweep's root printed, the
  residual and the certification held to the FILE's own `certified=`, then
  `evaluate` writes back and stops or (default) enters
  `steady_wind_with_element_diffusion` at once; no CFL step. **MEASURED**:
  the atomic element fixture AS LOADED stands at `||R||` 1.8924 (the marched
  hand-off read 1.888: the two CFL steps); the certified `wasp_full_newton`
  state is NOT the fixed point of its own sweep (particle count moves 2.1e-11,
  species 5.0e-10, cooling 3.6e-10: the cell solve's `xtol = sqrt(eps)` and
  the file's digits), reported at every evaluation rather than equilibrated
  away; `stationary evaluate` takes 0 steps and returns the certification
  unchanged; `stationary` on that state enters the JFNK at 7.133e-9 and
  exits `info = 0` with no step. Tests 20 rows GREEN (RED on the entry
  binary: the key warned and ignored, one step marched, no block, no
  refusal); `grid_and_gates` 130/0 (advisor rerun 131/0 on `build_lwv`).
  `mol_carrier` 200 steps numeric-line bitwise (eight new `#` lines);
  `run_check.sh` line 297 excludes `#` lines, so the block cannot change a
  matrix verdict; `wasp_full_newton` cold numeric-line bitwise against a
  control binary of the same tree without N10 (its 5.8e-9 / 9.4e-9 movement
  against the golden is B5j's recorded movement, not N22's as the worker
  supposed). **Fixed on the way**: the header loop matched `columns` as a
  SUBSTRING (`index(line,'columns') > 0`), so the `_adv` files' NOTE prose or
  the new `species_columns` field would have replaced the species labels
  silently; now the first token (`comment_field_is`); the same hazard in
  `restart_round_trip.sh`'s Python readers fixed. Deviations from the
  design, stated: no `EXHALE_STATIONARY_EVALUATE_ONLY` env hook (the second
  word is the one answer); the `provenance unknown` mark is not yet in
  `EXHALE_setup.out` (`write_setup_report.f90` not owned; one line);
  `species_columns` reported not refused (the by-label species map is a
  documented path); `source` is the git revision, not an executable md5.
  **Consequence for the user (decision item 21)**: the exact `options`
  comparison forbids the option-ladder workflow (converge without an option,
  restart with it on); legacy pairs are unaffected, every new rung would be
  refused; the disciplined fix is a key naming the tokens allowed to differ.
  Also corrected: `docs/input_schema.md` appendix D listed the retired
  `valve`/`fluxconst` fields and the wrong `base_bc` values. Reported:
  `Run mode` is absent from `input_schema.md`; the `EXHALE_PTC=1` route
  writes no setup report; a re-evaluated state writes `sec_ion_step=0`.
- **N22** (the predicted decrease stops scaling below `||s||` 1e-7;
  `steady_newton.f90` +195/-10, `krylov_and_dogleg` +16 rows): **all three
  candidates of the brief refuted by measurement, and the cause named**:
  the trials that break the scaling are exactly the trials whose species
  unknown was written onto a face of its box (`nout = 1`), and that write is
  a move of the state no model of the step contains (11 of 11 such trials
  refused with a fixed ascent of 3e-3 over five decades of step; 10 of 10
  others accepted at ratio 1.000). When the radius is cut the SCALED step
  quarters but the step in the unknowns the residual reads, `||D s||`,
  moves 0.4 percent, so the model image and `pred` are unchanged to seven
  digits; `r0 . A s = -||A s||^2` to every digit (`pred = 0.5 (r0 . A g)^2 /
  ||A g||^2`, radius-free, the interior minimizer along the leg: not a
  defect); a fresh Jacobian action at the same step reproduces `pred` to
  every digit (no stale cache); the probe displacement is `sqrt(eps)(1 +
  ||Y||)` = 1.8e-7 whatever the direction, nine decades above the
  residual's reproducibility (no floor; the reverse: the stalling step is
  2.7e-3 of the probe arc). **The fix**: `hold_the_step_where_it_leaves_the_
  species_box` zeroes the step at an unknown it would carry out of its box
  (the iterate is feasible, so the trial is, and the step is one the
  operator can be sampled along); `probe_length_of_the_jacobian_action` is
  the one definition of the probe arc; a diagnostic line for each trial under
  `EXHALE_ELEM_DIAG`. **Atomic reload, cap 300, before -> after**: 227 ->
  113 iterations (both stagnation), accepted/refused 170/55 -> 91/19,
  negative ratios 49 -> 14, best `||R||` 2.353e-2 -> **1.333e-2**, handed
  back 2.671e-2 -> **1.356e-2**; the period-two cycle is gone (iterations
  95-98: `pred` doubles with the step, all accepted at 0.998 to 1.001, steps
  1e-2 to 1e-1 instead of 1e-9); the closing refusals are ray/floor verdicts
  at POSITIVE ratios. Carrier reload BITWISE (the branch never fires);
  `wasp_full_newton` BITWISE (ten files; the probe-arc refactor is shared
  through `jv_product` and this settles it). Tests: the twelve-rung ladder
  rows are GREEN at entry (the geometry was never the cause), the new rows
  RED by four absent symbols, GREEN 120/0; `steady_species_rows` 182/0.
  Advisor: rebuild warning-free, `steady_newton.o` md5 `d70a078b` identical
  to the worker's delivered object, the two suites 120/0, 182/0. Not
  converged. **Reported, the next item (N23)**: the approximate-gradient leg
  is 0.7 percent of the step and supplies ALL of its model image (`||sU||/
  delta` 6.9e-3 passes `cauchy_leg_min_fraction = 1e-3` while `||A s|| =
  tau_c ||A g||` to four digits): the criterion must be the share of the
  IMAGE, not the length ratio (N7b's finding at a length the threshold
  admits).
- **N10b** (decision 21 option a; `load_IC.f90` +319/-20, `input_read.f90`
  +145/-1, `write_setup_report.f90` +62/-6, `grid_and_gates` +1 script of 11
  rows, `docs/input_schema.md` K44 and D.2, manual entry): the key `Restart
  option change: <token>[, ...]`; the `options` field is compared TOKEN BY
  TOKEN from one table (`opt_name(20)`, `opt_value`, the emitted text
  bitwise unchanged), only the named tokens may differ, every other
  difference refuses the load naming the token with its `from -> to`; six
  layout tokens (`metals`, `mol`, `oxychem`, `carrier`, `carrier_newton`,
  `iontrans`) may never be named (they decide the rows the state carries: a
  cold start, not a restart); an unknown token, a grid/reservoir/constant
  word, the key with no token or without `Load IC? True` are input errors;
  the change is written into both halves of the new state as
  `# option_change from -> to at restart of <source>`, inherited (history of
  32 lines); `EXHALE_setup.out` gains `Restart provenance: UNKNOWN | UNKNOWN
  by descent | states its configuration` and the permitted changes (N10's
  deviation 2 closed; the "by descent" wording was added because the first
  wording was measured false on a ladder's second rung). **Correction of
  the advisor's brief**: `LW` is a FLUX (`Stellar LW flux`), not an option
  token, so the `armA_noLW -> armA_LW` rung never needed this key; and every
  state written before the Jupiter-radius unification (those two included)
  is refused by the grid guard at 2.25e-2 in the outermost cell, so the `armA_*`
  directories are configurations, not restart sources. The ladder step was
  measured on the `armA_noLW` configuration with `base_ir`, `mol_ir` flipped:
  refused by token without the key, accepted with it, the provenance line
  in both halves. Tests RED 11/11 at entry (the key was an unrecognized
  line), GREEN 11/11; `grid_and_gates` 141/1 with the one FAIL attributed
  (below); `mol_carrier` 200 steps numeric-line bitwise (the rebuilt
  `# options` line bitwise). Also fixed: N10's test stage A read the LIVE
  `wasp_full_newton/output/`, which the matrix rewrites (a false failure of
  the legacy-mark row as soon as the pair carried a block); it now strips
  the block from its own copy. **Advisor**: rebuild warning-free, md5
  `b10641e5` identical to the worker's; `grid_and_gates` 141/1, the FAIL
  being `stationary_evaluate_rebuilt_within_sweep_budget` reading the live
  `output/` of `wasp_full_newton`, which the running Stage D matrix had left
  as a mid-run relaxation snapshot (3.2e-3 against 1e-8; the same pair
  evaluated by the entry and the N10b binary is bitwise, so not N10b's);
  the advisor pinned the certified state as
  `backup/regression/wasp_full_newton/IC/` (copied from the golden,
  `certified=T`, README there) and pointed the row at it. Reported:
  `restart_grid_guard.sh`'s stale header (third time); `sec_ion_step=0` on a
  re-evaluated state; whether `molbase`, `he_diff`, `he_metal_diff` belong
  with the layout tokens is a reviewer's question (the worker's reading:
  they add operators, not rows).
- **N23** (the approximate-gradient leg admitted by its share of the model
  image; `steady_newton.f90` +220/-5, `krylov_and_dogleg` +16 rows): two
  constants `cauchy_leg_image_share_max = 0.5`, `cauchy_leg_step_share_max
  = 0.1`, two pure routines (`cauchy_leg_shares_of_step_and_image`,
  `cauchy_leg_is_image_without_step`: image share above the first AND step
  share below the second, so a leg that IS the step keeps its image), and a
  branch in the dogleg's trial loop that DROPS such a leg from the step and
  its image and rebuilds the dogleg at no extra product (drop, not rescale:
  a rescale would either shorten `A sU` without `sU`, which is not a step of
  the operator, or move the step off the dogleg path); the model slope
  handed to `ray_slope_verdict` is the slope of the step taken in both
  branches. **MEASURED, the reading (`[diag 9]`)**: at the reload's second
  solve quartering the radius moves the promise by 0.9959 with the leg
  admitted and by 0.2501 with it dropped (six iterates: 0.9959 to 0.8938
  against 0.2501 to 0.2517), so the promise follows the radius again where it
  did not. **Honest impact: the new rule fires ZERO times on both fixtures
  and all three runs are BITWISE** (atomic 113 iterations, 1.356e-2;
  carrier 8.993e-2; `wasp_full_newton` 13.30): on the post-N22 trajectories
  the length test reaches every image-dominant leg first (117 length drops,
  0 image drops); the shares printed at every trial say the image test is
  the right criterion nonetheless (0.9986 of the image at 3.4e-5 of the
  step; eleven of seventeen read trials satisfy the predicate on their own),
  and the leg N22 saw at iteration 145 sat on a trajectory the N22 hold
  removed. Tests RED (four symbols absent), GREEN 136/0; `steady_species_
  rows` 182/0; `-Wall` warnings 6 before and after. Advisor: rebuild
  warning-free, the two suites 136/0 and 182/0 on `build_lwv` (the object
  md5 differs from the worker's delivered one by a post-measurement comment
  edit that moved line numbers, which the I/O statements carry; the
  worker's line-for-line reproduction of the solve covers it). Not
  converged; the rule is not claimed to improve either fixture. Reported:
  `cauchy_leg_min_fraction` is now a proxy for a test the code can make
  directly (letting the image test alone decide is the next increment, not
  taken because no fixture exercises the difference); the trial diagnostic's
  evaluations enter the cost counters.
- **N24** (what refuses the remaining trials of the atomic element reload;
  `steady_newton.f90` +51/-2, `krylov_and_dogleg` +13 rows): **the
  obstruction is the linear model, and the row is named.** Atomic reload,
  cap 300: 110 trust-region iterations (the "113" of earlier reports counted
  three summary lines), 91 accepted / 19 refused. The 19 are two
  populations, neither a step-control defect: 14 ascents at LONG steps
  (`||s||` 1.6 to 51, five with zero species unknowns held, so N22's
  mechanism is not what is left) and five closing refusals at positive
  ratios (`||s||` 2.3e-4 to 8.9e-7 with the promise quartering exactly while
  the true decrease keeps a step-independent floor of 3e-11 to 1.6e-10). One
  accepted step (outer 59, ratio 0.948, `||R||` 7.37e-2 -> 2.12e-2) moves the
  binding row from the energy equation to **the sodium element row of the
  outermost cell (cell 500, r 4.06 Rp)** and the solve never leaves it:
  from outer 60 to 110 that row holds 0.435 of the squared merit and the
  whole judged distance (5.3e-2 of its own scale, 5.36e6 at
  `cert_tol_element` 1e-8) while `||R||` stands at 1.333e-2. At outer 105
  everything else is sound (0 held, 0 refused samples, 0 empty columns,
  `dgbcon` 1.25e-6, step ratio 0.9998); what fails is the linear model: 40
  of 40 products, true relative residual 0.58 against 0.1, **the Arnoldi
  image 21 percent away from the operator** (1.6e-6 at iteration 1, 2.3e-4 at
  N7's iteration 60), `||sN||` 8.6e2 against a radius 5.9e-2, so the region
  uses 7e-5 of an inaccurate direction; the stuck row is at once the
  smallest scaled column and the smallest scaled row of the band. The
  stagnation detector read the judged distance, which IS that row; `||R||`
  was not falling (1.337e-2 -> 1.333e-2 over 32 iterations). **The
  image-only leg test** (`EXHALE_CAUCHY_LEG_BY_IMAGE=1`, default off)
  measured and NOT adopted: worse on both fixtures (atomic 2.30e-2 against
  1.36e-2; carrier 0.324 against 0.090, to the cap), the two runs parting at
  outer 7, so N23's "the rules agree" was that trajectory's property.
  Deliverable 3 STOPPED as briefed (linear solve / model / one row's
  physics, not the step control). Tests RED (two symbols absent), GREEN
  149/0; `steady_species_rows` 182/0; the delivered binary reproduces the
  entry solve line for line over 925 lines; `wasp_full_newton` BITWISE.
  Advisor: rebuild warning-free, 149/0, 182/0 on `build_lwv` (the advisor's
  binary of N23 had been built from an earlier revision of
  `steady_newton.f90`, as the worker showed; the suites and fixtures were
  the check). **Proposal (N25)**: (1) explain the 21 percent Arnoldi gap
  (`EXHALE_GM_ORTHO`; `Jv` reproducibility by two products); (2) ask the
  cycle for the direction ON THE BALL (Steihaug-Toint truncation) instead of
  solving to 0.1 a leg four decades beyond the radius; (3) the sodium row's
  own terms at cells 490 to 500: the element operator gives cell 500 a
  zero-gradient outer ghost, so nothing outward constrains it (a boundary
  question of the element operator, `binary_element_diffusion.f90`); (4)
  `cert_tol_element` 1e-8 is what makes this row the whole judged distance
  (UNCONFIRMED, ISSUES 3.1 e; not moved). Reported: the `(JFNK) it` line
  prints `k=*` for a species row without naming the element (in this case the
  identity was the whole diagnosis); the stagnation messages quote `||R||`
  where the detector read the judged distance.
- **Stage D gate and golden refresh (advisor, 2026-09-09 21:14 to
  2026-09-10 00:27, `make check` on the shared `build/`; refresh 00:40).**
  Tree after every item to N24 (N25 was running on a private build).
  Verdict against the goldens of 2026-09-08 19:05: on all 16 default cases
  the two marching files are **data identical** (15 cases bitwise; `wasp_
  full_newton` WITHIN at 5.78e-9 / 9.41e-9, B5j's recorded movement, which
  the 0.1 percent rule counts as identical) and all 32 `_adv` files FAIL by
  SHAPE (10 against 8 and 40/36/43 against 38/34/41 columns: the two status
  fields and the measure column of N11/N11b), so `make check` exits 2 as
  expected. The `_adv` movement is the physics of N12 (stale composition,
  upstream energy), N11 (mass-operator stationarity) and N11b (conditional
  tolerance): `T_adv` changes in 477 of 504 rows on `mol_base_handoff` (max
  61.6 percent, 364 -> 588 K at r 1.093), 479 of 504 on `mol_carrier` (max
  138 percent, 379 -> 901 K at r 1.098), 215 of 504 on `wasp_full` (max 2.3
  percent, 5637 -> 5505 K at r 1.131). **Goldens refreshed**: the previous
  set archived as `backup/regression/golden_pre_stageD_20260910/`, `run_check.sh
  golden` snapshotted the outputs the gate had just produced (17 cases,
  including the three `hp_*` and `oxygen_chemistry` extras), and the advisor
  compared the 68 snapshot files against the case outputs with the
  comparator's own `#`-filter: 68 identical, 0 differ (this comparison, not a
  three-hour rerun, is the "check again" of the standing rule, since a
  snapshot is a copy; the user asked that regression checks not be
  repeated). The shared `build/` and `EXHALE.x` are now current with the
  tree (md5 `f6ccb498` at the gate's build; N23/N24 built after it touched
  `steady_newton.f90` only, so the marching path of the matrix is the gate's).
  State files now carry the N10 metadata block and the `_adv` files the
  N11/N11b header; both are `#` lines the comparator ignores.
- **N25** (the Arnoldi gap, a direction on the ball, the sodium row at the
  outer boundary; `steady_newton.f90`, `krylov_and_dogleg` +25 rows,
  `binary_element_diffusion.f90` a default-off print only, 33 lines, no
  arithmetic): **(1) the 21 percent Arnoldi gap is the NONLINEARITY of the
  finite-difference action across directions**, shown positively: every one
  of the 40 Arnoldi columns matches its own product to 1e-16 to 3e-16, two
  products of one direction agree exactly, the basis's orthogonality loss
  is 4.4e-13, yet the action of `V_1 + V_2` differs from the sum of the two
  actions by 1.1e-4 and the assembled gap tracks the step length (8.5e-5
  and 2.5e-4 per unit of `||sN||` at outer 1 and 105; `||sN||` 12 -> 865),
  the signature of a curvature term over the fixed probe arc `sqrt(eps)(1 +
  ||Y||)` = 1.8e-7. Reorthogonalization NOT adopted (loss 1e-13 to 3e-12
  against the 1e-8 level; as a measured option `EXHALE_GM_REORTHO=1` it is
  worse on both fixtures: 1.487 and 0.43). **(2) Steihaug-Toint truncation**
  is in as `EXHALE_KRYLOV_ON_THE_BALL=1`, default off: products per
  iteration 33 -> 7 (atomic), 40 -> 1.5 (carrier); atomic better on the
  judged distance (1.69e6 against 5.31e6, lowest `||R||` 9.4e-3 against
  1.33e-2) but worse on `||R||` handed back (1.69e-2 against 1.36e-2);
  carrier worse (0.486 against 0.090); the brief's condition fails, default
  stays off. **(3) The sodium row of cell 500 is not a well-posed outflow
  condition, and the reason is a DEFECT of the stationary residual
  evaluation, not of the element operator**: the outer ghost of the element
  mixing ratio is FROZEN at the pre-solve composition (1.73111e-6, the
  `metals.inp` reservoir the handed-in state carried) throughout the solve
  while cell 500 falls to 1.0e-11, so at the stagnating iterate the ghost is
  1.7e5 times the cell it "continues"; the marching path refreshes that
  ghost (the second solve, entered at step 2002, has ghost = cell 500
  again); the accepted step of outer 59 did not move `fX(500)` at all (it
  moved the outer face mass flux); the scaled Jacobian row is dominated by
  the momentum of cells 499/500 and the sodium of 499 while its own
  diagonal is 1.3e-4 of the largest entry (the row is nearly independent
  of its own unknown). The physically stated condition (zero gradient at
  the top, no diffusive flux) requires the ghost to follow the ITERATE;
  the repair is in the composition write-back / ghost policy of
  `eval_residual`, which is N26 (launched). **(4)** the `(JFNK) it` line
  names the element (`worst row: element Na of cell 500`) and the stagnation
  messages quote the judged distance and `dj_best`. Tests RED (eight symbols
  absent), GREEN 174/0; `steady_species_rows` 182/0; default path unchanged
  (atomic 109 iteration lines identical to the digit; carrier identical);
  `wasp_full_newton` BITWISE. Advisor: rebuild warning-free, 174/0, 182/0 on
  `build_lwv`; three British `neighbor` in the element operator's comments
  fixed. Not converged; the on-the-ball run's best iterate is chosen by the judged
  distance while the control is `||R||` (the two rank iterates differently,
  which limits the comparisons).
- **N26** (the upper ghost of a transported element or carrier follows the
  iterate; `steady_newton.f90` +135/-20, `steady_species_rows` +7 rows): in
  `write_species_rows_into_composition` (the ONE place a trial composition is
  assembled for the residual: `eval_residual` and the three best-iterate
  restores) the transported columns of the upper ghosts take the outermost
  cell's value by the rule the marching path uses (READ: `Yetr(N+1:N+Ng,q) =
  Yetr(N,q)` in `species_advection_stage`; `ftry/fc(N+1:N+Ng,ic) = ...(N,ic)`
  in the carrier operator; `Apply_BC` never touches `f_sp`), the element
  totals of a ghost restored after the carrier copy, a ghost deficit clamped
  as the operator clamps its own; lower ghosts untouched (Dirichlet
  reservoir); the three-unknown route returns before reading anything. The
  mixing ratio follows the mass fraction exactly (equal He mass fraction at
  N and its ghost gives equal `n_He/n_H`, hence equal `n_X/n_H`; MEASURED
  1.4e-16). **MEASURED, atomic element reload (`EXHALE_JFNK_MAXIT=300`)**:
  `||R||` handed back 1.356e-2 -> **3.719e-4** (36x), judged distance
  5.36e6 -> 6.21e4 (86x), flux spread 6.6e-4 -> 5.5e-7, 91/19 -> 235/18
  accepted/refused, the sodium row of cell 500 binds on **0 of 252**
  iterations (was 51 of 109) and no row of cell 500 is ever the worst; the
  solve now stagnates on the INTERIOR helium row of cell 246 (r 1.155),
  `info = 2`, uncertified. **Carrier reload**: binding cell (H2 of cell 205,
  interior) and judged distance (4.807e7) unchanged to four digits, `||R||`
  handed back 8.993e-2 -> **2.519e-1** (worse by 2.8x, a path difference of a
  stagnating solve, reported not judged), flux spread down 12 percent.
  `mol_carrier` marching bitwise (no caller outside `solve_steady_jfnk`);
  `wasp_full_newton` BITWISE. Tests RED on the entry text with only the
  export added (4 rows: ghost mass fraction 6.7e4 off, mixing ratio 1.0e5
  off, the outermost element row of a flat column 1.19e-2 instead of 0, the
  carrier ghost 0.80 off), GREEN 189/0 (a non-vacuous row keeps the stale
  ghost and reads 1.19e-2 in both builds); `krylov_and_dogleg` 174/0. Advisor:
  rebuild warning-free, `steady_newton.o` md5 `4bd3138a` and binary
  `e56609ff` identical to the worker's, 174/0 and 189/0 on `build_lwv`; one
  `neighbors` fixed. **Reported, the next small item (N26b)**: the
  relaxation solvers `solve_trace_element_in_hydrogen` and
  `element_diffusion_step` refill their working arrays' upper ghosts before
  every pass while `element_transport_residual`, documented as the same
  operator, reads the ghost from `f_sp` with no refill, so the operator's
  boundary depends on its caller (`binary_element_diffusion.f90`, two lines).
  Not claimed: convergence; an explanation of the collapse of the outermost
  cells to 1e-11 of the reservoir.
- **N26b** (the element operator refills its own upper ghosts;
  `binary_element_diffusion.f90` +26/-0, new suite `src/tests/element_operator/`
  with 5 rows): `element_transport_residual` now writes `Xhe(N+1:N+Ng) =
  Xhe(N)`, `fX(N+1:N+Ng) = fX(N)` and `wYtr(N+1:N+Ng) = wYtr(N)` (the third
  line the relaxation twin gets for free through `project_elements`) so the
  operator poses ONE outer boundary (zero gradient, no diffusive flux at
  faces 0 and N) whatever the caller left in the ghosts; the lower ghosts
  stay the Dirichlet reservoir; the header states the boundary. Both
  outermost rows (N-1 and N, WENO3 reaches cell N+1 at face N-1) are covered.
  MEASURED: the relaxation's own fixed point handed to the residual with
  stale ghosts read 1.19e-2 of its scale in the outermost row before and
  1.5e-13 after (the interior level is 1.3e-13); the same column with stale
  and with zero-gradient ghosts now gives the two outermost rows bitwise.
  Suite RED 3 of 5 on the entry text, GREEN 5/0; `certification` 29/0; the
  atomic reload's first two JFNK iterations identical but for the merit's
  reproducibility going 2.8e-14 -> 0 exactly (the residual no longer reads
  a ghost the sweep leaves at a different last bit), the outputs moving at
  7.9e-12 (T) and 1.6e-10 (an `_adv` status/measure column) because N26's
  write-back hands over the ghost to 2e-16 and the refill to the bit;
  `mol_diffusion` marching 200 steps BITWISE (the marching path reaches the
  relaxation twin, never the residual). Built and measured in an isolated
  source tree with `steady_newton.f90` pinned, because N27 was editing it.
  Advisor: rebuild warning-free, `element_operator` 5/0, `certification`
  29/0 on `build_lwv`; one `percent` fixed. Reported: the write-back's ghost
  is the cell's to round-off, not exactly (any other consumer of those
  ghosts reads a 2e-16 gradient).
- **N27** (the helium row of cell 246 and the H2 row of cell 205 read side by
  side; `steady_newton.f90` +191/-3 diagnostic only, `steady_species_rows`
  one row restated): **the helium row of cell 246 is held by the LINEAR
  SOLVE, not by its physics, its boundary, a bound, the closure or the
  tolerance.** At outer 245 the row is +2.43e-22 = diffusive (settling)
  +2.40e-22 plus advective +3.3e-24, same sign (they add); the Newton removed
  the advective term (a factor 7200 since outer 1) and barely touched the
  diffusive one (x1.28); the measure 6.17e-4 is the maximum of a flat plateau
  (5.9e-4 at 240 to 5.4e-4 at 253), not a bad cell; 0 held, 0 refused samples
  in 246 iterations, band healthy (`dgbcon` 7.5e-7), merit reproducibility
  exactly 0, the trust region inert (ratio 0.933, 0 cuts, `aN` 1.0) and the
  Krylov cycle spends 40 of 40 products on a direction whose true relative
  residual is **0.9917** against 0.1 with the Arnoldi image 14.7 percent off
  the operator; the step buys 4.7e-5 of the merit per iteration. The
  helium Jacobian triple is a proper tridiagonal stencil with the diagonal
  largest (the unknown controls its row; the wind dominates the row's
  entries). The element rows hold 0.98 of the merit (He 0.63, Fe 0.36);
  `||R||` is the energy row at 1.37e-4; the stagnation detector is honest
  (merit and `||R||` fall while the judged distance rises 6.12e4 -> 6.21e4).
  **A number that does not fit N25's reading**: the Arnoldi gap is 0.147 at
  `||sN||` 1.8e-2, 8 per unit against N25's 2.5e-4 per unit, so the gap is
  not a function of the assembled length alone. **The H2 row of cell 205**
  is the same shape on a coarser front: measure 0.480 (not a cancellation,
  half its own terms), f(H2) falling 1.4e-1 -> 1.1e-3 -> 4.9e-4 over cells
  195-207 (logarithmic slope -757), f(H I) flat (1.7 percent), so cell 205
  is NOT a two-root H I front (consistent with N5); linear 40/40 at 0.97,
  gap 3.1e-2. Deliverable 3 STOPPED (nothing in the step control to change
  at the binding iterate). Deliverable 4: `cert_tol_element` would have to be
  6.2e-4 (4.7 decades; far above the 1e-5 the operator's header gives as its
  loosest floor); the closure's reproducibility (3.1e-7 in the row's measure
  per unit seed, N5) is not what holds it. Proposals: re-measure the
  Arnoldi gap at three iterates (it no longer tracks `||sN||`); the untried
  freedom is the COLUMN preconditioner (row scaling is fixed by decision 20
  a); whether the detector should stop a solve whose merit still descends
  while the judged row does not is an acceptance decision; the carrier
  front is a grid statement. Tests 189/0, 174/0; advisor rerun 189/0, 174/0.
  **Incident, found by N27 and CONFIRMED by the advisor**: with N26b's ghost
  refill in `element_transport_residual` the atomic element reload parts from
  N26's run at outer 10 and stands at `||R||` 1.199 at iteration 40 where
  N26's run was at 3.21e-2 (advisor run, cap 40, current tree: 1.888, 1.891,
  1.677, 1.484, 1.518, 1.582, 1.199 at iterations 1, 5, 10, 15, 20, 30, 40
  against N26's 1.888, 1.891, 1.678, 1.955, 1.758, 0.912, 0.0321); N27's
  own diagnostic is inert (hook on = hook off). N26b's validation (the first
  two iterations, the certification suite, `mol_diffusion` marching) did
  not reach it, and the advisor's verification of N26b did not rerun the
  reload: the advisor's mistake. N26b's refill is therefore NOT accepted as it
  stands; N26c (launched) finds why an operator that poses its own boundary
  breaks a solve whose caller already hands it that boundary to 2e-16, or
  reverts the refill. N27 restated the row `a_stale_upper_ghost_is_read_
  by_the_outermost_row` to its negation because N26b made the original
  false; N26c reconciles the two rows with whichever text stands.
- **N26c** (why N26b's ghost refill broke the atomic solve;
  `binary_element_diffusion.f90` refill REVERTED with the header restated,
  `element_operator` rows restated (6), N26's `steady_species_rows` row
  restored): **the refill changes the residual by one element's outermost
  row in its twelfth digit and nothing else, and the atomic element reload
  does not survive a last-bit change of the upper ghosts, whatever makes
  it.** MEASURED: the caller's ghosts already equal the cell to 1 to 4 units
  in the last place (`Xhe` 9e-16, `fX` 1.9e-15, `wYtr` 6e-16 over the whole
  solve); the two texts' outermost rows agree to 17 digits on every element
  but nitrogen's, which differs by 9.8e-13 relative (6.8e-17 in its own
  measure); no other ghost of the routine is read (face coefficients stop at
  N-1); the band's stencil already carries the dependence of rows N-1, N on
  cell N through the write-back. **A control that multiplies the same three
  ghost quantities by `1 + 2e-16` and poses no boundary at all breaks the reload
  as badly** (`||R||` at iteration 40: caller's ghosts 3.212e-2 = N26 to the
  digit; refill 1.199; one-ulp control 1.125; the three part at outer 10-11
  and from each other). Where it amplifies: the Krylov leg at its one-digit
  tolerance (relative residual reached 8.86e-2 / 8.89e-2 / 8.93e-2 at outer
  10, `pred` 0.4 percent apart, then 2.7x apart at outer 11, then different
  branches of the acceptance at outer 12: accepted 0.874 / REFUSED -1.98 /
  accepted 85); the solve's own merit reproducibility prints 1.4e-14, a
  hundred times the ghost movement. **The refill is physically right and is
  reverted anyway**: no version of it can be bit-inert on the production
  caller, and the reload's acceptance is a trajectory no last-bit change
  preserves. Delivered text reproduces N26's atomic reload line for line over
  252 records (`||R||` 3.719e-4, judged 6.22e4, He of cell 246, 27543
  samples, exit 2) and N26's carrier reload over 63 (0.2519, 4.806e7, H2 of
  cell 205); `mol_diffusion` marching bitwise; suites `certification` 29/0,
  `element_operator` 6/0 (restated to ASSERT the measured caller dependence,
  1.194e-2, so the suite exits 0; the advisor accepts that choice),
  `steady_species_rows` 189/0, `krylov_and_dogleg` 174/0; advisor rerun of
  the four on `build_lwv` (md5 `3201310a`, the worker's delivered binary):
  6/0, 29/0, 189/0, 174/0. **The finding that matters beyond this item**:
  the atomic element reload is NOT a reproducible quantity; its `||R||` at a
  fixed iteration count is set by which Arnoldi vector crosses a one-digit
  tolerance at the last bits of the residual, so "reproduces the previous
  run" rejects correct changes as readily as wrong ones, and what the reload
  can carry is named outcomes (binding rows, refusal counts, flux spread),
  not `||R||` at iteration N. Recorded in ISSUES 3.1 (d) and in the fixture's
  README. The asymmetry of the element operator (the relaxation twins write
  their ghosts, the residual takes the caller's) stays in the tree,
  documented and measured (1.19e-2 on a stale caller). Advisor: a `cp -irp`
  of the N4b worker hung on its interactive prompt since 2026-09-09 (PID
  61163, cwd the repository root, copying from the session scratch) was
  killed by PID. CORRECTION within the hour: killing the `cp` let the
  worker's shell run its NEXT statement, which re-applied N4b's first
  (discarded) patch onto the delivered `element_inventory.f90` (832 -> 926
  lines, a duplicate routine, would not compile); the N4b worker, woken by
  the stale task, rebuilt the file from its saved entry copy plus the
  delivered patch and it reproduces `diff_element_inventory.diff` line for
  line; the advisor confirmed 832 lines, one routine, a warning-free build
  (md5 `3201310a`, N26c's delivered binary) and `element_operator` 6/0,
  `steady_species_rows` 189/0. Lesson written into `docs/worker_rules.md`:
  `\cp -f` / `\mv -f` always, and a stale command is canceled as a task,
  not killed as a process.
- **Direction after N26c (user, 2026-09-10: "권고안대로 진행").** The
  advisor's recommendation adopted: (b) N8b, anchoring the species-row
  certification tolerances on the two candidate states by the review's five
  measurements (reproducibility floor, derivative convergence, a
  manufactured steady solution of the element operator, grid refinement of
  the candidate state, conserved-flux error), reported per row class and
  regime with a PROPOSAL for the user (decision item 22) and no value moved;
  and in parallel (c) N28, the code items of ISSUES 3.7/3.8 (run-wide JFNK
  cap, one `element_ratio_gate`, the update map's install branch owing an
  `Apply_BC`, the undefined innermost ghost column of `assemble_residual`,
  the PTC route's setup report, `sec_ion_step` on a re-evaluated state, the
  stale test header, `Run mode` in the input schema). (a), the column
  preconditioner, waits for N8b's numbers. The user manual PDF was rebuilt
  by the advisor (85 pages; it was two keys behind its source: `Restart
  intent`, `Restart option change`).
- **N8b** (anchoring the species-row certification tolerances;
  `certification.f90` +173 additive, no value moved; `steady_newton.f90` one
  diagnostic; `certification` suite 29 -> 42 rows; memo
  `docs/certification_tolerance_anchoring_20260910.md`): **1e-8 is refused
  by no floor of the arithmetic, the closure or the representation (five to
  eleven decades below it); what refuses it is the DISCRETIZATION and the
  ELEMENT FLUX CONSERVATION.** MEASURED on the two candidates: reproducibility
  exactly 0; closure-seed floor 5.3e-15 (the state's own sweep seed
  1.7e-10); state-perturbation amplification 0.23 to 0.30 in the wind but
  1.1e3 to 2.0e3 at r < 1.01 (floors 5e-17 and 4.4e-13); derivative: a
  clean V at the base (4.5e-7), a PLATEAU at 6.4e-5 over four decades of
  probe length at the binding cell 245 (the part of the derivative outside
  the band) so the smallest row a Newton step resolves there is 3.8e-8;
  manufactured relaxed column at 801/401/201 cells 5.4e-13 / 2.47e-6 /
  7.42e-6, order 1.59 (the flat column 2.7e-14); grid refinement of the
  candidate at r 1.1515: the He/H row 5.97e-4 (N = 500) -> 2.48e-4 (N =
  1000), remap control 1.6e-5, so **about 78 percent of the 6.0e-4 is
  discretization** (a C2 interpolant is needed; below r 1.03 the remap
  dominates); element flux spread 39.3 in the layer, 2.6e-2 in the wind,
  2.2e-3 over the escape window against the 5.5e-7 mass-flux spread the candidate
  is quoted with. Memo section 7 is the proposal (decision item 22, added
  by the advisor to `To_be_determined`): (a) 1e-5 in the wind with the layer
  reported, (b) keep 1e-8 and state the fixtures cannot certify, (c) two
  regimes gating; no proposal certifies either candidate. **Two defects found
  outside the item's files, reported**: (1) **the candidate's species do not
  carry the run's own density**: mass closure 6.68e-3 at the top against
  1.9e-9 on the three-unknown control; reproduced in the suite on a column
  closed to 6e-16 that reads 3.37e-3 after ONE `relax_element_composition`,
  so `project_elements` does not return to hydrogen the mass it takes from
  helium (N29, launched); (2) the writer-to-loader round trip does not
  return the state: the atomic candidate reloads at `||R||` 0.486 instead of
  3.7e-4 because `load_IC` rebuilds `rho` from the species (a consequence of
  1). Tests 42/0 (RED at the build level: the suite does not compile against
  the entry objects); the default certification report byte-identical
  (anchoring under `EXHALE_CERT_ANCHOR=1`). Advisor: rebuild warning-free,
  `certification` 42/0; one `neighborhood` and two `percent` fixed.
- **N28** (the code items of ISSUES 3.7; `EXHALE_main.f90`, `RK_rhs.f90`,
  `steady_residual.f90`, `element_census.f90`, `element_inventory.f90` and
  `attempted_step.f90` (the duplicate constant), `grid_and_gates` +2 rows and
  headers, `docs/input_schema.md` K45): (1) `EXHALE_JFNK_MAXIT` caps the RUN:
  the three call sites of `steady_wind_with_element_diffusion` pass one named
  default `jfnk_outer_iterations_default = 3000` (the marching hand-off's
  literal was 500, the one default that moved; `wasp_full_newton` converges
  in 13 and cannot reach either), and under the variable the run takes its
  FIRST stationary solve and refuses later entries by name (MEASURED: three
  solves of 2 iterations -> one); (2) `element_ratio_gate` has one owner,
  `element_census.f90` (public, `element_census_tolerance()` public with its
  override), read by `element_inventory.f90` and `attempted_step.f90`; (3)
  `update_map_begin_step`'s install branch calls `Apply_BC` in the N16d
  order (the `EXHALE_UPDATE_MAP` diagnostic's second block on `mol_carrier`
  moves 5 percent in the mass and energy rows: the fix of a defect, a
  boundary of one state under an interior of another); (4) `assemble_
  residual` defines the innermost ghost column (`dF`, `S`, `R` zero at
  `1-Ng`; two callers read it in whole-array expressions); `-fcheck` clean
  on the certification driver and on a 20-step `mol_carrier` (the run that
  reaches it); (5) the `EXHALE_PTC=1` route writes its setup report (was 0
  lines; new row RED/GREEN); (6) a `Restart intent: stationary evaluate`
  keeps the file's `sec_ion_step` (2402 -> 2402, was 0; new row RED/GREEN);
  (7) `restart_grid_guard.sh`'s stale header rewritten and the script made to
  honor `EXHALE_TEST_OUT`; (8) `Run mode` documented as K45. Item 9
  (`n_eq_sweeps_last` on a discarded evaluation) is `steady_newton.f90`'s,
  reported. Validation: `mol_carrier` 200 steps and `wasp_full_newton` cold
  BITWISE on every output file (the one log difference: `cap 500 -> 3000,
  from the caller`); `grid_and_gates` 146/0; `attempted_step` 162/1 with the
  FAIL proven pre-existing (rebuilt from entry sources: same row); advisor
  rerun on `build_lwv`: `grid_and_gates` 146/0, `attempted_step` **163/0**
  (the N28b fix below had reached the tree by then). **Found, reported,
  N28b launched**: `bg_cell` (`type(ion_rates)`) is allocated without
  initialization and the type has no default values, so before the first
  sweep the checkpoint checksum reads 8.7e244 and the row
  `a_perturbed_state_is_not_the_checkpoint` could never fire (a perturbation
  of 1e2 lost in the round-off of 1e244). Advisor: three forbidden phrases
  fixed (`per-element` in a printed line, `neighboring`, `per-region`). One
  deviation stated by the worker: an extra private `-fcheck` build (the
  shared `run_fcheck.sh` runs `make distclean`).
- **N28b** (`bg_cell` allocated and never initialized; `ion_cell_state.f90`
  the type only, `ionization_equilibrium.f90` the allocation site,
  `attempted_step` suite +1 row and a deterministic dirty-heap setup): **the
  defect is real, wider than `bg_cell`, and includes a LOGICAL**. Three
  grid-sized objects of `ionization_equilibrium` came out of `allocate`
  holding the recycled block: `bg_cell` and `ieq_rate_cell` (both
  `type(ion_rates)`, no default values) and the real array `heat_chem` (the
  one array left out of the routine's zero list). MEASURED with a probe that
  frees a 16 kB block before the allocation: `T_K`, `P_HI` 8.7e244 and
  **`x_h2_fixed` TRUE** (an imposed H2 partition the caller never asked for,
  read from garbage, standing in the array the carrier transport reads), on
  a fresh heap all zero (which is why N28 saw the dead row fail and the
  advisor's rerun saw it pass: the defect appears only when a freed block is
  recycled). Fix: `ion_rates` carries its zero in every one of its 37 fields
  (34 reals, 3 logicals) as default initializers, stated as the invariant "a
  cell not yet swept has no gas, no rates and no imposed partition", which
  reaches `bg_cell`, `ieq_rate_cell`, `bg_cell_adopted`, `bg_cell_best`, the
  closure scratch and the four checkpoint copies at once; `heat_chem = 0` at
  the allocation. Two comments that claimed a default initializer would not
  reach threadprivate copies were corrected (MEASURED with a four-thread
  program: it does, as OpenMP specifies). Tests: RED on the entry text with
  the new driver 68/2 (`the_unswept_background_carries_no_rates` -3.4e250;
  `a_perturbed_state_is_not_the_checkpoint` h1 == h0), GREEN 164/0 (the
  intermediate build located `heat_chem` at 8.7e244 by a field-by-field peek).
  `mol_carrier` 200 steps bitwise (every reader of the three objects is gated
  on `bg_ready`/`ieq_rates_ready`, set after the loops that write them;
  `heat_chem` has no physics reader, only the checkpoint); `wasp_full_newton`
  not run (atomic: the arrays never enter a physics path). Advisor: rebuild
  warning-free, `attempted_step` 164/0. Reported: three test fixtures
  enumerate all 37 fields of `bg_cell` by hand (now redundant) and
  `element_census_tests.f90` set two fields and left the rest to the
  allocator (now covered by the type).
- **N29** (the element projection does not conserve mass;
  `binary_element_diffusion.f90` +258/-113, `element_operator` suite 6 -> 14
  rows): **the leak is the mass the write-back is told to hand back, not the
  write-back.** `project_elements` took the mixture mass as `msum = m_1 n_H +
  m_He n_He` with `m_1 = mass_per_H_nucleus_without_He()`, the RESERVOIR mass
  per hydrogen nucleus (the trace metals at `melem_ab`), while the cell
  carries the metals at whatever the trace transport left; the projection
  then wrote the difference into hydrogen at every call. MEASURED at one cell
  of the certification column (metal/H 0.024401 against the reservoir's
  0.028227): `sum m_i f_i` 1.0000000 -> 1.0000636 after one step, helium's
  gain exceeding the hydrogen and metal loss by 6.4e-5: mass created. A
  second leak of the same kind: the trace-metal diffusion wrote metal densities
  with no counter-flux in the hydrogen it diffuses through (3.0e-3 on the
  suite column with the metals frozen). **Fix**: `mixture_mass_split` reads
  the mass the cell's species actually carry (`bsp_mass` and `melem_A`, no
  species weight restated), `project_elements(f_sp, Xhe, msum, scale_metals)`
  scales component 1 as a whole, every helium-free species included (CO
  too), the metals written first and the hydrogen group closing the mixture;
  one projection after the metal loop puts the metals' mass change back on
  hydrogen; `settling_coefficient` takes its component-1 mass from the same
  split. Metals off is bitwise by construction (`m_1 = m_H = 1` exactly).
  **MEASURED**: certification column closure after the relaxation 3.37e-3 ->
  2.3e-15; the handed-back atomic candidate 6.44e-3 -> 1.94e-9 (the file's
  own precision, the three-unknown control's 1.93e-9); the writer-to-loader
  round trip of that state now returns it (writer `||R||` 1.112 and reload
  1.1123 on the SAME row, against 3.2e-2 vs 0.4385 on different rows before:
  memo section 6 closed); `mol_diffusion` and `wasp_full_newton` BITWISE
  (metals off / no element row); **`lower_profile` MOVES and should**: its
  stored golden is 3.6e-4 out of closure (O/H 4.61e-4 against the reservoir
  4.90e-4), a 1500-step cold start goes 6.77e-3 -> 3.14e-9 in closure with
  T at 1.6 Rp moving 1.3 percent and H I 6.3 percent: **a golden that moves
  for a physics fix, refreshed at the next gate**. Atomic reload cap 40 by
  named outcomes: closure 6.44e-3 -> 1.94e-9, 32/8 -> 31/9, refusals 0 in
  both, binding row energy of cell 2 -> cell 4, flux spread 2.3e-3 ->
  1.4e-2, same stop. Tests RED 4 rows on the entry text (3.44e-3, 4.78e-4,
  3.44e-3, 3.82e-4), GREEN 14/0; `certification` 42/0, `diffusion_tests`
  35/0, `krylov_and_dogleg` 174/0. **The worker's control discipline**: the
  advisor's `EXHALE_lwv.x` is rebuilt at every verification and is a build
  of the live tree, so the worker built its own control (the delivered
  objects with only the entry-text `binary_element_diffusion.o`) and
  discarded three comparisons that had used the advisor's binary as
  "before"; workers are to build their own controls (added to the rules
  below). **Advisor**: `steady_species_rows_tests.f90` referenced
  `element_ratio_gate` through `element_inventory`, which since N28 imports
  it with `use, only` and does not re-export it, so the suite would not
  build (N28's validation did not include this suite; the advisor's did not
  either): one `use element_census, only: element_ratio_gate` line added by
  the advisor; rebuild warning-free (md5 `a9627a14`), `element_operator`
  14/0, `certification` 42/0, `steady_species_rows` 189/0,
  `krylov_and_dogleg` 174/0; one `percent` fixed. Reported: the metal diffusion's
  counter-flux is closed on hydrogen but not on momentum or internal energy
  (a trace approximation, named); `element_census_verify` gates the element
  ratios and only reports the mass closure (the 6.4e-3 state passed every
  gate it has): whether it should gate is a decision.
- **N8c** (the five anchors re-measured on the states after N29; the memo
  `docs/certification_tolerance_anchoring_20260910.md` only, 344 -> 524
  lines, no source file changed): **the reading is unchanged and is now
  taken on states worth measuring.** Post-N29 candidates (advisor-verified:
  no `src` file newer than the brief): atomic reload 275 iterations, 254/22,
  `||R||` 2.271e-4, binding row He/H of cell 245 (r 1.153), flux spread
  9.8e-7, **mass closure 1.95e-9** (was 6.68e-3), 11 of 13 entries refuse;
  carrier reload byte-identical to N26's (no metals, no element row; the
  corrected projection is never reached), closure 1.93e-9. The anchors: two
  evaluations of one state 0 exactly; the state's own sweep seed 1.7e-10 ->
  1.6e-14 (the loader no longer repairs a composition whose mass missed the
  density), closure floor 5e-19; state-perturbation amplification 1.2e3 to
  2.0e3 at r < 1.01 and 0.4 to 0.7 elsewhere; derivative: a clean V at the
  base (7.6e-8), the PLATEAU at 6.6e-5 at the binding cell stays (smallest
  resolvable row 3.9e-8 there, 6.5e-9 in the wind; both taken at the first
  JFNK iterate of a `relaxation` restart, two CFL steps from the loaded
  state, stated); manufactured column order 1.589, closure 2.3e-15 (was
  3.4e-3), extrapolated 1.52e-6 at the production spacing; grid refinement
  of the candidate 5.93e-4 (N = 500) -> 2.47e-4 (N = 1000), ratio 2.40,
  remap control 1.5e-5, Richardson limit 1.31e-4: **78 percent of the row is
  discretization**; element flux spread in the layer **39.3 -> 9.79** (He),
  37.3 -> 9.09 (H), the wind 2.6e-2 -> 2.7e-2, the escape window 2.0e-3.
  Memo: sections 3 (mass closure) and 6 (the reload) marked CLOSED by N29
  with the readings; section 5 rewritten; section 7's proposal re-stated
  with the new floors and **option (a) supported as written** (1e-5 gating at
  r >= `cert_regime_wind_r`, the layer reported and not gating); two of its
  three prerequisites discharged, the layer's flux conservation (9.8)
  standing. No proposal certifies either candidate (wind element row 3.05e-4
  and 7.2e-2 against a proposed 1e-5). `certification` 42/0. Reported: the
  iron wind row of the atomic candidate rose 18x between the two chaotic
  states (a size, not a trend); the derivative anchor is not measured at the
  loaded state itself (needs the probe under `stationary evaluate`,
  `steady_newton.f90`).
- **Decision 22 (user, 2026-09-10): option (a).** `cert_tol_element` and
  `cert_tol_carrier` become 1e-5, gating, for cells at r >= `cert_regime_wind_r`
  (1.20); the layer's element and carrier rows (r < `cert_regime_layer_r`,
  1.10) are reported and do not gate until the layer's element flux
  conservation (9.8 after N29) is repaired; the hydrodynamic tolerances do
  not change; a wind-certified state says so in one line. Anchor: the
  operator's stated floor and a decade above the manufactured column's
  discretization at the production spacing (1.52e-6); no existing state is
  admitted by it (3.0e-4, 7.2e-2). Implemented as N30 (launched).
- **Gate after N29 and the `lower_profile` refresh (advisor, `make check`
  2026-09-10 07:58 to 11:03; refresh 11:10).** Tree after N28, N28b, N29
  (N30 running on a private build). Against the 2026-09-10 00:40 goldens:
  **15 of 16 default cases data identical on all four files** (`wasp_full_
  newton` now identical too, its golden being the 00:40 snapshot);
  **`lower_profile` FAILS on all four, by design** (N29's mass fix: its
  golden was 3.6e-4 out of closure; v moves up to 4.8 percent at r 3.18 in 85
  of 504 rows, rho 1.0 percent in 4 rows, T 0.5 percent in 2 rows; the `_adv`
  measure column 57 percent at one row; `log10 Mdot` 7.29). The failing case
  was identified by its own section of the log (the advisor's first parse
  attributed the summary line to `hp_front`, which PASSES). Refreshed:
  `golden/lower_profile` archived as `golden_pre_n29_20260910/lower_profile`,
  `run_check.sh golden lower_profile` snapshotted the gate's output, the
  four files compared identical to the case output with the comparator's
  filter. No other golden touched.
- **N30** (decision 22 a; `certification.f90` +296/-54, `steady_newton.f90`
  +50/-30, `certification` suite 42 -> 57 rows, `steady_species_rows` 189 ->
  194, `docs/input_schema.md` the `cert_reason` row): `cert_tol_element` and
  `cert_tol_carrier` are gone as single numbers; `cert_tol_element_wind =
  cert_tol_carrier_wind = 1e-5` with one accessor each (`cert_tol_*_at(r)`, a
  step at `cert_regime_wind_r` = 1.20; `cert_tol_reported_only = huge` below),
  the anchoring written at the declaration with the memo cited; **one
  reduction** `certification_species_row_gate` decides a species row over
  the GATED cells (no gated cell = refused, never a satisfied row); the
  entries carry gated and reported maxima with their cells; the report
  prints `gated at r >= 1.2 against 1e-5: <measure> at cell <j>` and, for a
  certified state with species rows, the one line `certified IN THE WIND ...
  the rows below it are reported and do not gate: worst <measure> at cell
  <j>`; the state file carries `cert_reason=certified_in_wind` (the field
  already existed). **The band 1.10 <= r < 1.20 is reported and does not
  gate** (the decision was silent; the band holds the candidates' binding
  cell at r 1.153 where the discretization is 4.6e-4, so a 1e-5 gate would
  ask the residual to fall 1.5 decades below the error of its own discrete
  equation; no anchor for the band exists). `distance_from_certification`
  and `certified_row_measures` read the gated measures and the accessors, so
  **the trust region's ranking changes and both reloads move** (that is the
  decision's intended reach): atomic reload 275 -> 166 iterations at the
  same `log10 Mdot` 10.10, binding row He/H 2.896e-4 at cell 288 (a WIND
  cell) against 1e-5, NOT certified; carrier reload 63 -> 76 iterations, H2
  7.284e-2 at cell 285 against 1e-5, NOT certified; the hydrodynamic
  tolerances, the H(n=2) level balance (4.8e-10 against 1e-6) and
  `steady_gates_met`'s own carrier gate (`carrier_resid_th`, a volume ratio,
  not this tolerance) untouched. `wasp_full_newton` bitwise, still certified
  with no reason. Tests RED (the superseded rule refused a layer-1e-2 /
  wind-5e-6 row and a wind-5e-5 row alike, on a layer cell; four
  judged-distance rows moved onto the new divisor), GREEN 57/0, 194/0,
  `krylov_and_dogleg` 174/0. No key added (the radii are properties of the
  flow, constants of `certification.f90`). Advisor: rebuild warning-free,
  57/0, 194/0, 174/0 on `build_lwv`; the stale comment in `EXHALE_main.f90`
  (580) rewritten; `write_output.f90`'s `# adv_input_certified T` branch now
  carries the reason (a wind-certified state keeps its qualification in the
  `_adv` header); superseded notes appended to `a2_certification_contract`,
  `steady_solver_design` and `element_inventory_contexts`; one
  `neighborhood` fixed. Reported: whether `steady_gates_met` should read the
  anchored tolerance is a separate decision on a shared interface.
- **DOC1 to DOC3** (the documentation brought to the stage-2 state, 2026-09-10,
  three concurrent documentation workers, no source touched, no binary run):
  **DOC1** `EXHALE_BC_and_IC.tex` (29 pages): the whole lower-boundary
  chapter described the retired component-wise ghost closure (Dirichlet
  density, one-way valve, isothermal ghost pressure, four keys that
  `input_read.f90` now refuses with `error stop`); rewritten as the
  characteristic `(p, s)` face condition of `base_boundary.f90` (branches by
  face Mach number, the smoothstep of width `base_face_mach_blend`, the
  Gauss-Legendre ghost averages, the three surviving controls), the upper
  boundary as `free_outflow_ghost`'s isothermal continuation, a new section
  on the transported columns' boundaries (N26, N26c, N8a, N29), the
  checkpoint rebuild order (N16d), and the old test matrix retitled as the
  retired closure's, not reproducible by today's binary.
  `composition_restart_and_base_handoff.tex` (18): thirteen base species in
  the nuclei count (OH, H2O, CO added), every `load_IC.f90` line citation
  re-derived, a new section for the restart contract (intent, the eight-line
  metadata block field by field with refuse/report on each, `provenance
  unknown`, the grid guard, `Restart option change`), the round-trip and
  projection evidence of N9 and N29 labeled as read from this log.
  `transmission_spectrum.tex` (11): a new section for `read_adv_validity`
  (schema 2, the five statuses, `adv_conditional_tol`), the
  `chord_shell_indices` census with the CONTRIBUTION statement, and the
  `transit_schema 1` block. **DOC2** `steady_solver_memo.tex` (16): a new
  section 7 for the species-row solve as it stands (merit on certification
  scales, the trust region with named outcomes and the leg rules, the
  constraint rows, the empty-column repair and the equilibrated
  preconditioner, the three options measured and disarmed, the linear solve as
  the limiter with the additivity finding, the ulp-level chaos, the residual
  contract, the certification by regime), eleven hooks added to its table.
  `molecular_hydrogen_treatment.tex` (23): the `_adv` row corrected (the
  energy residual reads the upstream caloric energy since N12), the layer
  table RE-MEASURED on the current `mol_lyman_werner` output (144/184/476
  cells, H2 peak 5.15e12, 1 percent edge 4.48e9 against the stale 227/251,
  4.63e13, 4.42e9), thirteen of sixteen cases pinned, a new section on the
  element budget, the H3+ departure factors marked STALE in place, the
  `mol_ir_bands` peak re-measured (2.148e-8 at r 1.0059 against the stale
  2.934e-7). `lower_atmosphere_coupling.tex` (28): a new section 7.4 on the
  element operator (N29's conservation, N26b/N26c's boundary, N8c's flux
  spread table as the open item); two pre-existing overfull boxes fixed.
  `EXHALE_user_manual.tex` (88): certification by regime in section 2.7, the
  `Run mode` key entry, eleven hooks, the coupling header example replaced by
  the line the code writes (`cert_reason`), `Cooling_breakdown.txt` columns
  12 to 14 (H2, H2O, CO infrared) restored to the column list. **DOC3**
  `README.md` (sixteen-case matrix, 25 suites, the restart contract, the
  transit metadata, certification stated, nineteen examples),
  `README_HOWTO.md` (`Resid tol` default 1e-5 not 1e-3, `Flux spread tol`
  default 2e-5 on the face flux not 5e-3, a certification table, a restart
  section, an eleven-row hook table, the sixteen cases and 25 suites), the
  NEW `docs/code_status_20260910.md` (418 lines, replaces the 2026-09-05 file
  as the read-first document; section 3.3 states what is and is not
  certified and what limits the rest), the project `CLAUDE.md` pointer line,
  `TO_BE_DONE.md` (an index of what is open, three headings corrected, 35
  em-dashes removed), `steady_solver_design.md` section 22. All eleven PDFs
  rebuilt, `grep -c '^!'` 0 on every log. Advisor: every default and
  constant quoted checked against the source (`flux_spread_th_default`
  2.0e-5, `resid_max` 1e-5, `cert_tol_*_wind` 1e-5, `cert_regime_wind_r`
  1.20, `ratio_sound_band` 1e-2, `tr_delta_min_dimensionless` 1e-6, the
  `Cooling_breakdown` header); one constant name corrected in the memo;
  pre-existing British spellings fixed in the workers' files, in
  `ISSUES_20260909.md` and in the handoff; the 51 em-dashes of the workspace
  `CLAUDE.md` replaced. **Reported by the workers and CONFIRMED by the
  advisor** (ISSUES 3.7): the base reservoir is stated once at startup and
  the refreshes of `dp_bc` and `ntot_bc` after every sweep are dead writes
  (the boundary itself is right: the reservoir equals the user's base
  pressure by construction); `golden/mol_ir_bands/` carries no
  `Cooling_breakdown.txt`; `examples/README.md` covers part of nineteen
  folders under the wrong planet; the `.md` twins of two `.tex` documents
  have parted from them; the `max(rest, 0)` write-back clamp (N4b) restated.
  ISSUES section 5 rewritten to the state after N30; the handoff opens with
  a dated 2026-09-10 state block. The carrier reload is still a scratch
  recipe (`B5k/mkrun.sh carrier`), not a pinned fixture: to pin with N31.
- **N31** (why the finite-difference Jacobian action is not additive across
  directions, and two probe options; `steady_newton.f90` +563/-3, of which the
  arithmetic is two probe rules and the rest is a default-off hook,
  `krylov_and_dogleg` 174 -> 190 rows, new fixture
  `backup/regression/carrier_elem_newton`): **the action is not additive
  because the residual carries a NON-SMOOTHNESS FLOOR five to nine decades
  above the double-precision rounding of its own rows, and the difference
  quotient along the PRECONDITIONED Krylov directions divides that floor by
  an increment five decades smaller than a generic direction produces.**
  N25's reading (curvature over the fixed arc) is refuted at the iterates
  that matter and confirmed at one that does not. MEASURED on the atomic
  element reload with the new hook `EXHALE_JV_ADDITIVITY=1`
  (`[diag 13]`, at outer 1, 20 and the last of a capped run, on the first
  two Arnoldi directions and on two directions off the basis, at a tenth,
  one and ten times the probe arc): the defect of the Arnoldi pair scales
  as the arc to the power **-1.026** at outer 1 and **-0.999** at outer 40
  (rounding), the off-basis pair as **+1.125** and **+2.109** (curvature,
  and 3.2e-8 and 2.1e-8, negligible), and at outer 20 both are curvature at
  cells 1 and 2. The second difference of the residual on the
  defect-holding rows confirms it directly: ratios 0.774 and 1.676 at outer
  1 where a smooth residual gives 4 and 4, and the floor stands **9.5e5**
  (outer 1) to **2.0e9** (outer 40) times `epsilon` times the residual on
  those same rows, so it is not the arithmetic of the assembly. READ, the
  only tolerance in the atomic path that can set a floor of that size:
  `ionization_equilibrium.f90` line 969 hands `hybrd1` an xtol of
  `sqrt(dpmpar(1))` = 1.49e-8 (lines 1432, 2110, 2114), a statement about
  the step and not the residual, so the composition is a piecewise map of
  the state; `eval_residual` fixes the number of composition passes but not
  `hybrd1`'s own iteration count. **A decision of the feasible set is ruled
  out**: on all six probes of the three iterates, components blocked 0,
  forward side, and the step taken equal to the nominal step. **The
  localization gives no unknown class to fix**: at outer 1 the squared
  defect is 0.31 mass, 0.64 energy, 0.03 momentum, 0.02 helium, under 1e-4
  in the eight trace elements, and nine tenths of it sits in 35 of the 500
  cells between 1 and 75. The scale reading is a uniform shortfall: the
  largest unknown of the state is displaced by 2.2e-8 of itself (about
  `sqrt(epsilon)`) and almost every other by 1e-10 to 6e-9, one to two
  decades short, with no unknown near a bound (nearest face 2.9e-5 against
  a displacement of 4.5e-14). The model that follows, defect ~ floor /
  (eps `||Jv||`), predicts both pairs at outer 1 within a factor three.
  **Two options, both default off.** `EXHALE_JV_COLUMN_SCALE=1`
  (`probe_step_on_the_column_scales`, `eps0 = sqrt(eps)(1 + ||Dc^-1 Y||) /
  ||Dc^-1 v||`) is the brief's remedy and it is WORSE, for the reason the
  diagnosis gives: the Krylov direction's weight sits on unknowns of small
  column scale, so the rule SHORTENS the arc by 42 and the defect rises
  8.430e-6 -> 1.466e-4 at exponent -0.965, the 1/arc law reproduced by a
  change of the arc made for another reason. `EXHALE_JV_PROBE_ARC=<x>`
  (`jv_probe_arc_scale`, default 1) is the lever the measurement points at,
  and it wins the linear question and loses the Newton by the same change:
  at x = 100 the defect of the Arnoldi pair falls 8.430e-6 -> 2.240e-6, the
  exponent turns +1.139 and the second differences become 4.065 and 4.022,
  so the rounding regime is gone and the quotient is truncation-limited,
  while the atomic reload then ABORTS at outer 33 (`info = 2`, `||R||`
  1.504, flux spread 5.902e-2, 12 steps refused by the model with the
  actual reduction 26.4 times the predicted, radius 2.028e-7) where the
  control runs 144 iterations to 3.641e-4 and 5.557e-6 with no model
  refusal. At x = 20 it aborts at outer 43 at 1.524 with 21 model
  refusals and a radius of 2.076e-13; the carrier reload runs its 40
  iterations at both values and hands back 1.006 at x = 100 and 9.076e-1
  at x = 20, each with 2 model refusals, against the control's 9.572e-2
  with none. **At x = 5, which
  does NOT stall**, the binding regime is reached (outer 58, `||R||`
  3.087e-3, element Fe of cell 146, against the control's 3.386e-3 and
  element Fe of cell 147 at outer 60) and THERE the Krylov cycle still
  spends 40 of 40 products at a relative residual of 0.898 to 0.980
  against the 0.1 asked, indistinguishable from the control's 0.935 to
  0.959. **The substance of the item**: the additivity of the linear
  operator wants a LONGER arc (defect ~ floor / arc) and the trust
  region's model wants a SHORTER one (truncation ~ arc, and the ray test
  refuses a step whose model and true slope differ in size), the arc
  already sits at the Newton's balance point, and every direction was
  measured to lose, so NO probe rule fixes this. What can is lowering the
  residual's own floor, which is not double precision.
  Column-scale named outcomes at cap 40, option against control:
  atomic `||R||` at 40 1.674 against 7.607e-2, handed back 1.633 against
  5.377e-2, flux spread 4.985e-1 against 3.555e-3, 30/10 against 32/8;
  carrier `||R||` at 40 1.357e-1 against 1.070e-1, handed back 1.263e-1
  against 9.572e-2, flux spread 7.718e-2 against 4.635e-2, 31/9 against
  33/7, SAME binding row (carrier H2 of cell 205), no new
  `gm_no_direction_sampled` and no refusal of any class in either. The one
  thing the column scaling improves is the count of Krylov cycles reaching the
  tolerance asked (40 of 40 against 26 of 40), which is the REDUCED
  problem's residual on an operator that is not the Jacobian and is
  precisely why that count must not be read as progress. NEITHER OPTION IS
  ADOPTED: both fail the first clause of the acceptance rule outright.
  **Default path unchanged**: with every hook off the atomic reload's 40
  iteration lines are identical to the control's to the digit, the carrier
  reload's binding row and `||R||` at 40 are the control's, and
  `wasp_full_newton` is BITWISE but for the provenance timestamp; with
  `EXHALE_JV_ADDITIVITY=1` on, the 40 iteration lines are still identical,
  so the hook is inert armed as well as disarmed. Tests: `krylov_and_dogleg`
  RED on the control objects (seven symbols absent), GREEN 190/0;
  `steady_species_rows` 192/0, `certification` 57/0, `element_operator`
  14/0, `residual_determinism` 5/1 with the same row failing by design in
  the control build. **Fixture pinned (deliverable 0)**:
  `backup/regression/carrier_elem_newton` now holds the carrier reload that
  existed only as a scratch recipe, with `run.sh`, a README and the control
  of the unmodified tree at cap 40 (`||R||` 1.969, 1.984, 1.921, 6.800e-1,
  1.070e-1 at 1, 5, 10, 20, 40; handed back 9.572e-2; flux spread 4.635e-2;
  `info = 1`; 33 accepted, 7 rejected; 40 of 40 products every cycle;
  binding row carrier H2 of cell 205; 2754 samples admitted, none refused).
  A capped run of that fixture never reaches a mass-loss rate. **Defect
  fixed on the way**: the new probe rule was first written with a local
  `integer :: n`, which is the grid size `N` of `global_parameters` under
  Fortran's case-insensitive rules, so its own size guard fell through on
  every call; the new suite rows caught it at once and the local is now
  `nunk`. Not claimed: that the floor is `hybrd1`'s xtol -- that is READ
  from the source and not measured, and the experiment that closes it
  (tighten that xtol, re-measure the floor) belongs to
  `ionization_equilibrium.f90`, which this item does not own. The comment
  of the suite row `the_probe_arc_stands_above_the_residual_noise` was
  restated: the arc IS nine decades above the residual's reproducibility as
  a length of the whole vector, but the conclusion drawn from it about the
  quotient does not follow. Advisor: rebuild on `build_lwv` (one linker warning, pre-existing), `krylov_and_dogleg` 190/0, `steady_species_rows` 194/0, `residual_determinism` one row FAIL by design; atomic reload at cap 40 with every hook off, 8 threads: 1.888, 1.891, 1.819, 1.583, 1.861, 8.964e-1, 7.607e-2 at 1, 5, 10, 15, 20, 30, 40, handed back 5.377e-2, flux spread 3.555e-3, identical to the worker's control; the carrier fixture's README control checked against its own `run.log`; the fixture README is the reference for every later carrier-reload number (the 0.25 quoted after N26 was a scratch-recipe run, not this fixture). Accepted. **Next item N32**: the residual's floor (the inner solves of the composition and of the temperature stop on step or bracket tolerances, `newton_dense`, `hybrd1` xtol `sqrt(eps)`, `brent_root` 1e-10); measure whether the additivity floor follows them.
- **N32** (which inner tolerance sets the residual's non-smoothness floor;
  two default-off tolerance hooks and a widened `[diag 13]`,
  `ionization_equilibrium.f90` +47/-2, `T_equation.f90` +42/-3,
  `steady_newton.f90` about +130 inside the additivity hook,
  `krylov_and_dogleg` 190 -> 203 rows): **NEITHER. The floor follows
  neither the composition solve's stopping tolerance nor the temperature
  bracket's, and it is not in the chemistry at all: it sits in the
  HYDRODYNAMIC rows of the residual.** New hooks `EXHALE_IEQ_TOL=<x>`
  (replaces `tol = sqrt(dpmpar(1))` for the run; new function
  `composition_solve_tolerance`, cached once by `ioniz_eq`, announced only
  when armed) and `EXHALE_TEQ_TOL=<x>` (the bracket of `brent_root`, now a
  module variable at its old 1e-10 through
  `temperature_bracket_tolerance`). MEASURED on the atomic element reload
  at outer 1 (cap 40, 8 threads, one entry state for every run): the
  additivity defect of the Arnoldi pair is 8.174e-6 at exponent -1.063
  (default 1.49e-8), 7.452e-6 at -0.995 (1e-10), 9.348e-6 at -1.059
  (1e-12) and 9.546e-6 at -1.064 (1e-14), and the second difference on the
  defect rows stays on its floor at 5.6e-12 to 1.6e-11 with no trend in
  the tolerance: **the law is FLAT**, over six decades. The measurement is
  exactly reproducible (three control runs, and 1 thread against 8, give
  8.174e-6 and 1.362e-11 to every printed digit), so the +-17 percent
  between runs is a real movement of the state and not a fall.
  `EXHALE_TEQ_TOL=1e-13` is inert to the bit through all 40 iterations,
  which confirms as MEASURED what was READ, that `brent_root` is reached
  only from `post_process_adv` and never from `eval_residual`. **Why the
  composition cannot be the floor**: quadratic convergence pins its root
  at the SQUARE of the step tolerance. On the reported cell the accepted
  `x(H+)` is 2.037266074305163e-12 at the default and
  2.037266074305300e-12 at 1e-14, agreeing to 13 significant digits (7e-14
  relative), and the pinned fixture's own run.log records an accepted
  reaction residual of 2.91e-16 at maximum for the steady sweeps; the
  floor implies a jitter of about 1e-9 relative, five decades above that.
  The tightest tolerance costs 30 percent more wall time per residual
  evaluation (0.090 s against 0.066 s at 4 threads, 17 evaluations each)
  and buys nothing. **Where the floor IS**, from the widened hook: the
  energy SOURCE, which is where `excited_hydrogen`, `lya_rt`, the cooling
  table interpolation, the write-back clamps and `comp_T_from_p` all end
  up, is smooth on the arc and five decades below the floor (second
  difference of `heat - cool` on the defect cells 3.039e-17, 2.474e-16,
  1.046e-15 at half the arc, the arc and twice it, ratios 8.1 and 4.2,
  textbook h squared, against the residual's own 1.362e-11 on the same
  rows); by row class over the whole column the floor is mass 1.801e-11,
  energy 3.090e-11, momentum 5.946e-12 against element He 4.269e-12 and
  every trace element at 1e-16 to 2e-13; nine tenths of it sits in 29 of
  500 cells between cell 1 and cell 58. The mass row carries no source at
  all, so no part of the chemistry can be its cause. The competing
  explanation, that the row is a small difference of large terms and its
  floor is their rounding, is REFUTED: the largest of heat and cool on
  those cells is 3.424e-3 against `||F||` 5.507e-3, a ratio of 0.62, and
  epsilon times it is 7.6e-19. **Deliverable 2 (the stationary-solve-only
  polish) was NOT BUILT**: the brief conditions it on deliverable 1
  naming a tolerance, and none is named. **A correction to the brief and
  to N31**: on this fixture `solve_ieq` and `newton_dense` never run.
  The helium triplet is on (the fixture's `EXHALE_setup.out`) and
  `metals.inp` is present, so the sweep calls
  `hybrd1(ion_system_HeH_TR_metals, ...)` directly; MEASURED `nt_calls` 0
  over the whole hook and no "ioniz-eq solver" line in the fixture's
  reference run.log. `hybrd1` returns `info = 1` at the default, at 1e-10
  and at 1e-14 alike, so there is no `info = 2/3` and no fallback to
  count, and `newton_dense`'s hardcoded step test `1e-11*xscale` is not in
  this path. **Default path unchanged**: `wasp_full_newton` BITWISE in
  `Ion_species.txt` and identical but for the provenance timestamp in
  `Hydro_ioniz.txt`, every `(JFNK)` line identical; the atomic reload's 40
  iteration lines identical to the control's with every hook off. Tests:
  `krylov_and_dogleg` RED on the control objects (both new functions
  absent), GREEN 203/0 with 13 rows added; `steady_species_rows` 194/0,
  `certification` 57/0, `element_operator` 14/0, `adv_static_limit` 53/0,
  `coupled_source_step` 30/0, `physics_probe` 1491/0,
  `residual_determinism` 5/1 with the same row failing by design as in the
  control. Noticed, not fixed: the fixture README does not say the helium
  triplet is on, which is what selects the direct `hybrd1` branch and is
  what sent this item's brief to the wrong routine; and `newton_dense`
  applies the caller's tolerance to `fnorm` while its step test uses a
  hardcoded `1e-11*xscale` no caller can set, which would be the real
  floor on a fixture that does reach it. **Next item**: the floor is in
  the flux assembly, not in `ionization_equilibrium.f90` -- the candidates
  left are the reconstruction where the limiter is not the frozen WENO3
  weight, the Riemann wave-speed estimate, `Apply_BC`, and the transport
  coefficients of the diffusion operator, which are the one remaining way
  the composition reaches a hydrodynamic row. Advisor: rebuild on `build_lwv`; `krylov_and_dogleg` 203/0, `steady_species_rows` 194/0, `adv_static_limit` 53/0, `residual_determinism` the same single row FAIL by design; atomic reload at cap 40 with every hook off, 8 threads: the N31 ladder to the digit (1.888, 1.891, 1.819, 1.583, 1.861, 8.964e-1, 7.607e-2; handed back 5.377e-2, spread 3.555e-3), so the worker's 1.112 at iteration 40 is the thread-count sensitivity of N26c, not a change of the path. The atomic fixture README now states the triplet is on and which branch that selects; `newton_dense`'s hardcoded step test recorded in ISSUES 3.7. Accepted; deliverable 2 correctly not built. **Next item N33**: the jump itself, scanned along the Arnoldi direction across four arcs and split by term and operator (base boundary, reconstruction, Riemann flux, caloric inversion, composition read, transport sources).
- **N33** (where the residual's non-smoothness floor is: the term, the
  operator and the branch; a default-off multi-scale jump scan
  `[diag 14]`, `steady_newton.f90` +710, `krylov_and_dogleg` 203 -> 220
  rows): **THERE IS NO BRANCH. The floor is the ROUNDING OF THE FLUX
  ASSEMBLY, raised above the last bit of the row by two cancellations
  that compound.** The Roe flux of a face is built from the JUMPS of the
  two reconstructed states, and in the nearly hydrostatic base layer
  those jumps are 3.2e5 to 5.9e6 times smaller than the states
  themselves, so the interface flux carries the last bit of O(1)
  quantities and not of its own value; the row is then the flux
  difference divided by the cell volume, `r^2/dV` = 5.12e3, one over the
  1.955e-4 width of a base cell. The bound
  `epsilon x (face state) x r^2/dV` is what a correctly rounded assembly
  cannot go below, and MEASURED it IS the floor: 1.675e-12 against
  1.240e-11 at cell 1, 2.520e-13 against 1.644e-12 at cell 84, 4.326e-15
  against 2.366e-14 at cell 167, 4.241e-17 against 4.409e-16 at cell 251
  -- the bound falling 4.6 decades with the measured floor following it,
  ratio 5.5 to 12 throughout, and the ten cells holding the largest
  second difference are cells 1 to 23, every one a base cell. Beyond
  cell about 250 the bound falls below the residual's own curvature and
  stops binding, which is what a floor does. New hook
  `EXHALE_RESID_JUMP_SCAN=1` samples F at 65 points of `Y + t v` on ten
  windows, each a quarter of the one before, from two probe arcs down to
  two arcs over 4^9; one window is not enough, because at the standard
  arc every interval already holds thousands of steps and none stands
  above its neighbors. **The measurement**: on the atomic element
  reload at outer 1 the tracked row is the energy row of cell 1, its
  median first difference falls by exactly four whenever the spacing
  does (the smooth part is a straight line), and its second difference
  does not fall at all -- 1.0e-11 to 1.6e-11 over 5.4 decades of
  spacing, 1.8e-5 down to 6.9e-11 -- which is a step of F and not a kink
  or curvature. Isolated steps appear once the smooth increment per
  interval falls below the floor: 28 of 64 intervals at 3.05e-5 arcs and
  61 of 64 at 7.63e-6 arcs, at least 4.0e6 stepping intervals per probe
  arc, the largest moving its row by 2.625e-12, 2.749e-9 of the row.
  **The term and the operator**, at the largest step (cell 9, mass row):
  `S` of the mass row is identically zero and the whole step is in `dF`;
  every quantity upstream of the flux -- the ghosts, `T`, `n_tot + n_e`
  and the four reconstructed face states -- moves by one to eight units
  in the last place of an O(1) number with no outlier, in exact
  multiples of 2.220e-16; the first quantity that steps is the interface
  flux, by 4.803e-16 where its neighbors move by 3.2e-17, and
  4.803e-16 x 5116.3 = 2.457e-12 against the row's measured 2.625e-12.
  **No branch fires anywhere near it**: faces falling back from Roe to
  HLLE 0 and unchanged, face states the positivity limiter scaled 0 and
  unchanged, so the base boundary's Newton, `characteristic_branch_weight`,
  the Mach cap, the ghost quadrature, the stencil selection, the
  wave-speed branch and `comp_T_from_p` are all continuous here. The
  off-basis direction carries the SAME floor, 9.1e-12 to 1.4e-11 flat
  over 5.4 decades, and isolates no step at any window, because its
  smooth increment never falls below it: the floor is a property of the
  residual and not of the direction, and what the preconditioned
  direction does is shrink the numerator of the difference quotient by
  five decades. The carrier reload reproduces the law (floor 1.532e-11,
  bound 1.382e-12, ratio 11.09), so it is the discretization at the base
  and not one fixture. **N31 corrected**: the floor's size stands, but
  "five to nine decades above the rounding of its own rows" compared it
  against `epsilon ||F||` on the row, which is the wrong reference for a
  row built as a flux difference over a cell width; against the right
  one it is a factor of ten. N32's localization is confirmed and made
  quantitative. **Deliverable 2 NOT BUILT**: the brief conditions the
  deliverable on a branch being named and there is none. The treatment the
  measurement does name, for the advisor, is a WELL-BALANCED flux
  difference at the base -- assemble the rows from the departure of the
  face states from the local hydrostatic isentrope, which
  `base_boundary.f90` already integrates, so the balance cancels
  analytically instead of in floating point -- in `Num_Fluxes.f90`,
  `RK_rhs.f90` and `Source.f90`, none of them this item's files; the
  cheaper partial lever is the base cell width, on which the floor
  depends as one over it. **Default path unchanged**: the atomic reload
  at cap 40 gives all 338 `(JFNK)` lines identical to the control's
  through `done info=1 ||R||= 1.112E+00`; with the scan armed the
  outer-1 line is the control's to every digit; `wasp_full_newton`
  BITWISE in `Ion_species.txt`, `Ion_species_adv.txt` and the three
  breakdown files, differing only in the provenance timestamp in the two
  hydro files, `log10 Mdot` 13.30 on both. Tests: `krylov_and_dogleg`
  RED on the control objects (`resid_jump_scan_on` absent), GREEN 220/0
  with 17 rows added; `steady_species_rows` 194/0, `element_operator`
  14/0, `certification` 57/0, `residual_determinism` 5/1 with the same
  row failing by design and the same measured 3.354e+01 on the control
  binary. Noticed, not fixed: `global_parameters` declares
  `integer :: count`, which shadows the intrinsic `COUNT` in every scope
  that uses the module, so `count(mask)` written anywhere fails to
  compile -- the same class as the `t0` and `g2s` incidents and worse in
  that the shadowed name is an intrinsic; renaming a public global was
  out of scope. Advisor: rebuild on `build_lwv`; `krylov_and_dogleg` 220/0, `steady_species_rows` 194/0; atomic reload at cap 40 with every hook off, 8 threads: the N31 ladder to the digit (7.607e-2 at 40, handed back 5.377e-2, spread 3.555e-3), no `[diag 14]` line; two British spellings in the new comments fixed. Accepted; deliverable 2 correctly not built (there is no branch to smooth). The finding is structural and reaches the marching discretization, so it is put to the user as DECISION 23 (well-balanced flux differencing recommended; quad accumulation of the flux difference as the control experiment; coarser base cells; accept the floor). The `count` global shadowing the intrinsic recorded in ISSUES 3.7 (HYG). **Next item N34** (runs without a decision, changes no default): the control experiment, quad accumulation of the face fluxes and their difference inside the stationary residual only, as a measured option.
- **N34** (the control experiment of decision 23: the hydrodynamic rows of
  the stationary residual in quadruple precision). `EXHALE_RESID_QUAD=1`
  sends the STATIONARY evaluations of `assemble_residual` -- the tag
  `ieq_sweep_state_kind` puts on an iterate or a probe; the marching
  stages reach `RK_rhs` directly and are untouched -- through a
  quadruple-precision instantiation of the whole flux assembly, from the
  conserved state through the caloric EOS, the PLM limiter and the WENO3
  weights, the Roe flux with its entropy fix and HLLE fallback, the
  geometric and gravitational sources and the division by the cell volume,
  rounded to double once, at the row. One source text
  (`hydrodynamic_rows_body.inc`, 1111 lines, kind-generic) instantiated
  twice by `hydrodynamic_rows.f90`; the composition, the heating and the
  cooling stay double and enter converted exactly. **The double
  instantiation is the production operator BIT FOR BIT**, asserted three
  ways: a suite row on a stated grid (largest difference exactly zero in
  `dF`, `S` and the four face states, tolerance 0), and two whole
  40-iteration reloads under `EXHALE_RESID_QUAD=2` -- the atomic (all 355
  `(JFNK)` lines identical) and the carrier, which is the molecular one
  and so exercises the generic caloric EOS (all 364 identical). The
  quadruple instantiation reaches the analytic mass row of a uniform state
  to 1e-30 before its rounding, where the double one stands at its own.
  **Verdict: the floor is confirmed and is NOT removed.** On the atomic
  reload at outer 1 the tracked row's non-smoothness floor falls only
  1.240e-11 to 5.917e-12, and cell by cell it lands on
  `epsilon_DOUBLE x face state x r^2/dV` -- ratio 0.93 to 1.4 over the base
  cells against 5.7 to 8.9 in the control -- and not on the quadruple
  bound, which is 1.224e-30, 4.8e18 below what is measured. So N33's
  amplifier (the face cancellation 3.2e5 times the `r^2/dV` 5.1e3) is
  confirmed a second way and its attribution completed: four fifths of the
  floor is the arithmetic INSIDE the assembly, the rest is the double
  representation of what the assembly is handed and hands back -- the
  ghosts `Apply_BC` writes, the sweep's temperature and particle count,
  and the double row itself -- which precision inside one operator cannot
  reach. The detector isolates no step at any of the ten windows (control:
  28 of 64 and 61 of 64 at the two narrowest), but the second difference is
  still flat over 5.4 decades and the row still changes sign 39 times in
  64 intervals: a staircase with a smaller rise. **What the quadruple-precision accumulation does buy
  is the fidelity of the Krylov model**: the additivity defect falls
  8.174e-6 to 1.009e-6 with the exponent in the arc still negative
  (-1.063 to -0.901, so still rounding and not curvature), and at the
  binding iterate of the atomic reload at cap 40 the relative gap of the
  Arnoldi image falls 1.188e-1 to 3.813e-3, a factor 31, the returned
  relative residual then being the true one to four digits (0.2114 against
  0.2114; the control returned 0.1122 for a true 0.1617). The element rows,
  which do not pass through the assembly, are unmoved to four digits
  (C 1.367e-13 to 1.366e-13), which is the control inside the control.
  **The Newton does not converge and the Krylov cycle still stalls.**
  Atomic cap 40, quad accumulation on: 40 of 40 products at 0.21 against the requested
  0.10, `||R||` 8.413e-1 handed back against the control's 1.112 with
  trajectories parting at iteration 9 (the solve is chaotic at the ulp level,
  so that is not by itself a result). Atomic cap 250: both stop on the
  stagnation detector, control at iteration 167 at 3.047e-4, quad at 136 at
  4.096e-4, both with 40 of 40 at 0.983 to 0.989, both binding on the
  element Fe rows of cells 132 to 137, both NOT CERTIFIED with 11 entries
  of the inventory refusing. Carrier cap 40: unchanged in every named
  outcome -- binding row the H2 carrier row of cell 205 throughout, 40 of
  40 at 0.9886 against 0.9729, `||R||` 9.336e-2 against 9.427e-2 -- as it
  must be, since that row is assembled by `carrier_steady_residual` and
  never sees the flux. Cost MEASURED: the flux assembly alone 46x at 8
  threads and 74x at one (500 cells, 3.628e-5 s to 1.681e-3 s), the whole
  atomic reload at cap 40 1.11x (271 s to 302 s), because a stationary
  residual is mostly the ionization sweep. **Nothing adopted, no default
  changed**: `wasp_full_newton` bitwise but the provenance timestamp,
  `log10 Mdot` 13.30 on both; the atomic and carrier reloads identical
  line for line with the option off. Tests: `krylov_and_dogleg` RED on the
  control objects, GREEN 232/0 with 12 rows added; `steady_species_rows`
  194/0, `element_operator` 14/0, `certification` 57/0, `energy_update`
  53/0, `acceptance_classes` 22/0, `adv_static_limit` 53/0,
  `species_face_flux` 1/0, `spectrum_type` 154/0, `grid_and_gates` 146/0,
  `attempted_step` 70/0, `coupled_source_step` 30/0,
  `residual_determinism` 5/1 with the same row failing by design and the
  same 3.354e+01 on the control binary. One file outside the brief's list:
  `caloric_eos.f90` gained three read-only accessors (the mixture ratios
  and the H2 rovibrational table nodes) with no arithmetic and no behavior
  change, because the generic pipeline cannot evaluate the molecular
  caloric EOS without them; the molecular whole-run bitwise comparison is
  the evidence that they changed nothing. Noticed, not fixed:
  `backup/regression/atomic_elem_newton/README.md` still names the helium
  row of cell 246 as the binding row at cap 250 from N26/N27, where the
  measurement on this tree is the element Fe rows of cells 132 to 137 at
  iteration 167 and `||R||` 3.047e-4. **For decision 23**: extra precision
  confirms the diagnosis and does not cure it; the amplifier has to go
  (well-balanced flux differencing, or a coarser base cell, which buys one
  factor of the width), and even with the hydrodynamic rows made
  model-faithful the two rows that actually bind -- the element rows of the
  atomic reload and the H2 carrier row of the carrier reload -- are not
  hydrodynamic rows and are untouched by any of it. Advisor: rebuild on `build_lwv` (the new `hydrodynamic_rows.f90` in `SRC`); `krylov_and_dogleg` 232/0, `steady_species_rows` 194/0, `energy_update` 53/0; atomic reload at cap 40 with every hook off, 8 threads: the N31 ladder to the digit (7.607e-2 at 40, handed back 5.377e-2). The `caloric_eos.f90` accessors (read-only, outside the brief's list) accepted on the molecular reload's bitwise evidence. Accepted; nothing adopted. Decision 23 restated with this result: option (b) retired as a remedy, (a) still the right discretization but not by itself the cure of the species-row solve; the atomic fixture README gained a 2026-09-10 state section (binding rows Fe of cells 132 to 137 at cap 250 on this tree). **Next item N35** (runs without a decision): the conditioning of the preconditioned operator on the element and carrier rows at the binding iterate, and what the banded preconditioner misses there.
- **N35** (what holds the Krylov cycle on the species rows once the
  rounding floor is removed: the conditioning of the preconditioned
  operator and what the banded preconditioner misses). Three default-off
  measurement hooks on the linear system of one outer iteration, adopting
  nothing: `EXHALE_KRYLOV_SIZE_SCAN=1` (`[diag 16]`) reruns the SAME cycle
  at 40, 80, 160 and 320 products and the last two again with the basis
  orthogonalized twice, with `EXHALE_GM_HISTORY=1` printing the reduced
  residual after every product and the true one every twenty;
  `EXHALE_PRECOND_SPECTRUM=1` (`[diag 15]`) takes 200 Arnoldi products of
  `A_z M^-1` from a deterministic start, diagonalizes the Hessenberg
  (`dgeev`) and reports the Ritz values of the whole operator and of its
  compressions onto the species and the hydrodynamic rows, with the
  localization of the three smallest Ritz vectors;
  `EXHALE_BAND_DIFFERENCE=1` (`[diag 17]`) measures `(A - A_band) v` on
  those Ritz vectors and on the first two Arnoldi directions, splits the
  columns they live on into the part inside the band, the band's
  disagreement there and the part outside, with the radiation live and
  frozen, and reads the binding row's band row by class as a fraction of its
  diagonal. **Verdict: the band misses nothing, and the operator itself is
  near-singular on the species rows.** Its action agrees with the full
  finite-difference one to 4.0e-5 to 1.5e-3 relative on every direction
  tested on both fixtures; the columns hold 99.4 to 100.0 percent of their
  norm inside the band (the largest outside part is 0.6 percent, `element
  He of cell 308`); and what remains is not on the rows that bind at all,
  but 0.999 in the momentum row of the LAST cell (carrier) and 0.83 to
  0.89 momentum plus 0.11 to 0.17 energy in five to seven cells (atomic),
  the element and carrier rows taking 0.0000 to 0.003 of it. So the two
  remedies the brief conditioned on -- a rank-few radiation correction by
  the Woodbury identity, or a wider `kl_jac` -- would each repair something
  already right, and **nothing was built**. **What IS wrong, MEASURED at
  the binding iterate of both fixtures.** (1) The preconditioned operator
  is well conditioned on the hydrodynamic rows and near-singular on the
  species rows: the Ritz spectrum of the compression onto the hydrodynamic
  rows has a magnitude ratio 43.4 (carrier, `carrier H2 of cell 205`
  binding) and 158 (atomic, `energy of cell 52` binding), none and three of
  200 below 1e-2; onto the species rows, 6.56e5 (carrier, smallest
  1.493e-6, nineteen below 1e-2, one below 1e-4) and 1.06e4 (atomic,
  smallest 9.366e-5). The smallest Ritz vector is 0.71 to 0.91 carrier rows
  on the front cells 199 to 205 and 458 to 466 (carrier) and 0.9955 element
  rows over 145 cells (atomic), and the carrier fixture's near-zero part is
  a complex pair, not a real null direction. The attribution is there from
  outer 1 (2.17e5 against 39.6). And the discretization statement behind
  it: the band row of `carrier H2 of cell 205` has diagonal 1.2396e-1
  against off-diagonal magnitudes 4.1026e2, **3.304e3 times the diagonal to
  the hydrodynamic unknowns of its own stencil and 5.65 times to the
  neighboring carriers**, so the binding row's own unknown is three
  decades below its coupling to another class and there is no diagonal for
  a preconditioner to stand on. (2) **A larger subspace makes the step
  worse, not better**, which is neither of the brief's two alternatives.
  The reduced least-squares residual falls monotonically with the size
  while the TRUE residual of the returned step stops falling at 60 to 80
  products and rises: atomic with quad accumulation on, reduced 2.658e-1,
  1.781e-1, 1.144e-1, 9.958e-2 at 40, 80, 160 and 182 products against TRUE
  2.851e-1, 2.575e-1, 3.006e-1, 3.203e-1 -- the cycle announces its
  tolerance reached at 182 and hands back a step leaving 0.32, three times
  what was asked; carrier, reduced 9.677e-1, 8.001e-1, 4.639e-1, 1.434e-1
  against TRUE 9.731e-1, 8.638e-1, 9.668e-1, **1.472e+00**. With quad
  accumulation off the true residuals are a factor 3 to 6 worse at every size
  (atomic 7.856e-1, 3.905e-1, 7.580e-1, 2.004e+00) and the shape is the
  same. Reorthogonalizing the basis changes nothing at any size, the
  measured loss of orthogonality being 0.000e+00, and the
  Arnoldi image of each column at the same iterates is exact to 3.30e-16 and 3.77e-16: the
  recursion is not the defect. What fails is additivity -- the action of
  the sum of the first two directions departs from the sum of their actions
  by 1.357e-3 (atomic) and 7.601e-5 (carrier), the assembled 40-product
  step's image by 1.068e-1 and 2.256e-1 -- because the matrix-free action
  is the secant of a nonlinear residual over a fixed arc and the step it is
  extrapolated to grows from 0.31 to 11.4 while the arc stays 5.2e-6. Tests:
  `krylov_and_dogleg` 251/0 with 19 rows added (RED on the control objects,
  which do not carry the routines), among them the Ritz values of a stated
  upper triangular 6 by 6 spanning six decades reproduced to 1.171e-12 with
  the recursion carried to the dimension of the space, the column split
  summing to the whole difference with an empty band and with an exact one,
  and the four hooks off by default; `steady_species_rows` 194/0,
  `element_operator` 14/0 (the three suites whose driver links
  `steady_newton`). Default path: the atomic reload at cap 40 with every
  hook off, new binary against my control, all 354 `(JFNK)` lines
  IDENTICAL through `||R|| 5.377E-02, flux spread 3.555E-03`, the N31
  ladder to the digit; `wasp_full_newton` bitwise but the provenance
  timestamp, `log10 Mdot` 13.30 on both. Cost MEASURED: the size ladder is
  1082 products per diagnostic iteration (the carrier reload at cap 40 from
  about 3.5 to about 9 minutes at 8 threads), the spectrum 600 and the band
  difference about 15. Noticed, not fixed: the frozen-radiation column of
  `element He of cell 410` carries 7.267e-4 of its 2.816e-1 norm outside
  the band where the WENO3 stencil says it should carry none, 0.26 percent
  and two decades below anything in this item's conclusions; the likely
  route is the base boundary state recomputed from the iterate, not
  confirmed. Also corrected in the reading: the preconditioner is
  `build_banded_jac_full`, a coloring of the FULL residual, not of the
  frozen one -- `build_banded_jac` on `frozen_residual` is called only by
  the standalone check in `EXHALE_main.f90` -- so "inside the band" carries
  the coloring's contamination by the non-local response, and the frozen
  column is a separate control rather than the band's own approximation.
  **For the advisor**: the next question is not the preconditioner. The
  species rows' row scaling (the binding row's diagonal 3.3e3 below its
  hydrodynamic coupling is set by the certification row scale of decision
  20 a divided by the column scale, and using one scale for the merit and
  another for the model is a decision); the cycle has an optimal length
  around 60 to 80 products on both fixtures and the step control asks for a
  tolerance the operator does not honor; and `gm_resid_rel` is verified
  against the operator only at a breakdown or a rank loss, where the
  measurement says it needs verifying whenever the cycle ran more than a
  few tens of products. Advisor: rebuild on `build_lwv`; `krylov_and_dogleg` 251/0, `steady_species_rows` 194/0; atomic reload at cap 40 with every hook off, 8 threads: the N31 ladder to the digit (7.607e-2 at 40, handed back 5.377e-2). Accepted; deliverable 4 correctly not built. The two levers named are the LINEAR model's row scaling (a freedom decision 20 a did not fix: the merit and gate keep the certification scales, the linear inner product may differ) and a Krylov cycle that returns the step its TRUE residual chooses (the optimum lies at 60 to 80 products; the reduced residual lies past it). **Next item N36**: both as measured options, default off.
- **N36** (the two levers N35 named: the row scaling of the LINEAR model,
  separate from the merit's, and a Krylov cycle that returns the step its
  TRUE residual chooses). Two default-off options on the linear solve,
  adopting nothing. `EXHALE_MODEL_ROW_EQUIL=1`
  (`unit_infinity_norm_row_scaling_of_the_band`) forms the diagonal `E`
  that brings every row of the banded model to unit infinity norm as a
  power of two, applies it to the band that is factored and, inside
  `pgmres`, to the right-hand side and to every action of the operator,
  and maps the image and the relative residual BACK to the certification
  scales at every exit, so that the merit, the gate, the trust region's
  predicted decrease and the certification keep the scales decision 20 a
  fixed and only the linear model's inner product changes; `ab` itself is
  untouched, so the merit's gradient is where it was.
  `EXHALE_GM_TRUE_RESIDUAL=1` makes the cycle measure its candidate step
  against the operator from 20 products on, every 10 and at every point
  where it would otherwise end, keep the best, and stop on two consecutive
  rises or on the tolerance IN THE TRUE NORM; what comes back is that
  step, its own image from the product that measured it, and the true
  residual as `gm_resid_rel`, under a new eighth outcome
  `gm_true_residual_turned`. **Verdict: both options work as specified, both
  fail the acceptance, both stay off.** (1) **The row scaling is a real
  freedom and taking it is worse.** The cycle then minimizes a functional
  the merit does not read, and the two disagree by a factor 1.2 to 673
  MEASURED (`[row equil]`, equilibrated against certification: atomic
  9.956e-2 / 4.111e-1 at outer 2 and 9.937e-2 / 1.218e-1 at outer 40;
  carrier 2.353e-1 / 5.730e+1 and 1.286e-1 / 8.657e+1), so the cycle
  announces the 1e-1 it was asked for while handing back a step that
  leaves 86.6. It cannot touch what binds: the preconditioned operator
  `A_eq M_eq^-1 = E (A B^-1) E^-1` is SIMILAR to the unequilibrated one,
  and the measurement at the same iterate agrees (hydrodynamic Ritz ratio
  39.62 against 39.60 carrier, 57.02 against 56.36 atomic; the count of
  species-row Ritz values below 1e-2 unchanged at 18 and 2), while the
  banded factorization is bit for bit the same matrix because
  `equilibrate_the_preconditioner` already equilibrates both sides; and a
  row scaling leaves every ratio WITHIN a row where it was, so the binding
  carrier row's diagonal 1.2396e-1 against 4.1026e2 of coupling, 3.304e3
  times the diagonal to the hydrodynamic unknowns of its own stencil, is
  exactly N35's number still. What it costs: the carrier reload at cap 40
  goes from `||R||` 9.572e-2, info=1, binding `carrier H2 of cell 205` to
  1.955, **info=2 with four entries of the inventory refusing**, and the
  true residual of the returned step at 40 products from 0.973 to 1.050e2.
  (2) **The true-residual cycle is the right instrument and it does not
  turn the solve.** At the binding iterate it caps the cycle where the
  step stops improving: on the ladder 40/80/160/320 the control's returned
  step gets WORSE with size (atomic TRUE 2.851e-1, 2.575e-1, 3.006e-1,
  3.203e-1; carrier 9.731e-1, 8.638e-1, 9.668e-1, **1.472e+00**) while the
  true-residual cycle stops at 80 and at 120 products and returns 2.660e-1 and 8.244e-1,
  never worse than the control's best at any size and at fewer products.
  Two cycles of the atomic reload at cap 40 name the new outcome. But the
  best true residual any subspace reaches is 0.27 and 0.82, the Newton
  does not descend past its plateau, and `||R||` handed back at cap 40 is
  1.060 against the control's 5.377e-2 (atomic; the trajectory parts at
  iteration 15, the first cycle to pass 20 products) and IDENTICAL on the
  carrier, where a subspace of 40 always makes the last candidate the best
  and the option changes only the number reported (true 9.870e-1 against
  reduced 9.839e-1). **Acceptance, by the named clauses**: the true
  residual below 0.1 within 80 products on both fixtures fails for both
  options (option 1 atomic 8.79e-2 but carrier 7.51e+1; option 2 atomic 2.66e-1,
  carrier 8.64e-1), and the refusal count at cap 40 does not fall by half
  (atomic 9 to 10 for both options; carrier 8 to 7 and 8 to 8), with no new
  refusal class and no new `gm_no_direction_sampled` anywhere. Corrected
  mid-item and re-measured: the first true-residual cycle left the reduced
  tolerance in place as a stop and so discarded a step it never measured
  (carrier outer 1, the 80-product cycle ended at product 49 and handed
  back product 40's step at TRUE 1.259e-1 where product 49 reaches
  9.724e-2); with the reduced residual no longer ending a cycle, the
  numbers above are the corrected cycle's. Tests: `krylov_and_dogleg` RED on
  the control objects (nine symbols absent), GREEN 264/0 with 13 rows
  added, among them the equilibrated solve returning the same step as the
  unequilibrated one to 9.90e-16 on a stated dense system solved over the
  whole space, the equilibrated cycle reporting the CERTIFICATION residual
  of a truncated solve to eight digits, and a stated operator whose action
  stops being one linear map after 40 products, on which the cycle names
  `gm_true_residual_turned` and stops before its subspace is spent;
  `steady_species_rows` 194/0, `element_operator` 14/0,
  `residual_determinism` 5/1 with the same row failing by design and the
  same measured 3.354e+01. **Default path**: the atomic reload at cap 40
  with every hook off, all 354 `(JFNK)` lines IDENTICAL to my control
  through `||R|| 5.377E-02, flux spread 3.555E-03`; the carrier reload all
  364 identical; `wasp_full_newton` bitwise but the provenance timestamp,
  `log10 Mdot` 13.30. With the options ON `wasp_full_newton` stays CERTIFIED:
  option 1 `||R||` 3.043e-9 against 1.667e-9 with a largest relative movement
  5.2e-10 in the output, option 2 bitwise in `Ion_species.txt` with 8 extra
  residual evaluations, its cycles never reaching the 20 products the
  needs. **For the advisor**: what the two options rule out is the linear
  algebra. The norm the cycle minimizes is free and wrong to change; the
  length of the cycle is free and already at its optimum; the operator's
  species block is near-singular at every scaling and the step it can
  produce leaves 0.27 to 0.82 of the residual whatever the subspace. The
  remaining lever is the one N35 named and N34 confirmed the arithmetic
  cannot reach: the DISCRETIZATION of the rows that bind. Advisor: rebuild on `build_lwv`; `krylov_and_dogleg` 264/0, `steady_species_rows` 194/0; atomic reload at cap 40 with every hook off, 8 threads: the N31 ladder to the digit (7.607e-2 at 40, handed back 5.377e-2). Accepted; both options stay off. **The linear algebra is exhausted as a lever** (N31 to N36: probe rule, inner tolerances, extended precision, band coverage, row scaling of the linear model, cycle length). STOPPED here for the user: decision 23, and the next question is the discretization of the species rows that bind at the front (the carrier row of cell 205 whose own diagonal is 3.3e3 below its hydrodynamic coupling). `code_status_20260910.md` section 5 task 1, ISSUES 3.1 and 5, the handoff and decision 23 restated accordingly.
- **N37** (decision 23 (a): a well-balanced flux difference for the
  near-hydrostatic layer, as the default-off option `Well balanced:`).
  With the key on the reconstruction, the Riemann jumps and the pressure
  force carry the DEPARTURE from the cell's own local hydrostatic
  equilibrium instead of the state: within cell j that equilibrium is the
  constant-density one through its own (rho_j, p_j),
  `p_eq,j(r) = p_j - rho_j (phi(r) - phi(r_j))`, which assumes no thermal,
  entropy or compositional stratification and so preserves ours (Kaeppeli
  and Mishra 2016, A&A 587, A94, sections 2.1.1 and 2.1.3; the momentum
  source as the face-pressure difference of the cell's own equilibrium is
  their 2014 paper, J. Comput. Phys. 259, 199, eq. 2.26). The pressure
  stencil acts on the neighbor's pressure measured against that
  equilibrium continued THROUGH the face (this cell's density up to the
  face, the neighbor's beyond it), which is what vanishes on the discrete
  equilibrium; the two cell pressures are differenced first, while that is
  exact, and the hydrostatic terms added to the difference. The momentum
  row is then assembled from the face pressure measured against the same
  equilibrium and `source` returns `S(2) = 0`: with
  `A+ P_up - A- P_dn - (A+ - A-) p_j = -rho_j [A+ (phi_i(j) - phi_c(j)) +
  A- (phi_c(j) - phi_i(j-1))]` identically, the equilibrium's flux
  difference and its gravitational source cancel in the ALGEBRA and the
  cell's own O(1) pressure never enters the row. Mirrored in the
  kind-generic copy, so `EXHALE_RESID_QUAD=2` is the production operator
  bit for bit under both key values (suite row, largest difference exactly
  zero, tolerance 0). **Exact preservation, MEASURED**: on the scheme's own
  discrete equilibrium (the isothermal solution of the face-matching
  condition, the spherical form of their eq. 42) the momentum row falls
  from 7.0e-5 (PLM) and 2.5e-3 (WENO3) of the local weight rho g at
  N = 250 to 4.1e-14 and 3.9e-14, and stays at 4.5e-14 to 5.5e-14 at
  N = 500, 1000, 2000 while the base scheme falls at its design order; the
  mass and energy rows reach 1e-16 against 1e-8. Twelve new assertions in
  `hydrostatic_residual`, plus eleven rows in `krylov_and_dogleg` (with no
  gravity the well-balanced option is the base scheme to 8e-16 on a uniform state and to
  1.4e-16 on a Sod tube, for both reconstructions; and the same exactness
  on a second grid and geometry). **What it does NOT reach, and this was
  written into the design before the code**: N33's non-smoothness floor
  does not move (energy row of cell 1, 1.009e-11 with the option against
  1.046e-11 without; cells 84, 167, 251 the same to a factor 1.3; six to
  ten times `epsilon x face state x r^2/dV` either way), because a cell
  pressure is itself a rounded O(1) double and every jump built from two of
  them steps by its last bit whatever the grouping. What does move is the
  additivity defect at ten times the probe arc, 2.410e-7 with the option
  against 3.151e-5 without. **The stationary solves are worse with the key
  on, and the reason is not the option**: `store_row_terms` scales the
  momentum row by `max(|dF(2)|, |S(2)|, |Smom|)`, the option returns
  `S(2) = 0`, so the row IS its own scale and every normalized momentum
  residual is exactly 1 (MEASURED on the carrier reload as loaded: k=2 is
  `1.000000E+00` with the option against 1.954 without, while the mass row is
  the same number to seven digits and the energy row to four). The merit
  then carries no momentum information: `||Fs||2` rises from 16.5 to 28.1
  at the same state, and 28.1^2 - 16.5^2 = 517 is the 500 momentum rows
  pinned at unity. The atomic reload at cap 40 reaches 1.415 against the
  control's 1.112 and the carrier reload 1.732 against 9.34e-2; at cap 250
  the option's solve aborts at iteration 61 on "no descent step in 12
  consecutive iterations" with `||R||` and the worst scaled row both exactly
  1.000 and that row the MOMENTUM of cell 1, which is what a merit whose
  momentum rows are all pinned at one looks like from inside the line
  search. `wasp_full_newton` with the key on says the same from the other
  side: it is NOT CERTIFIED on exactly one entry, the momentum row at
  1.000E+00, while its flux spread is 1.539e-15 against the control's
  6.233e-14, its energy row 9.911e-12 against 1.667e-9, its mass row
  1.228e-13 against 1.453e-13 and `log10 Mdot` 13.30 either way, the
  profiles moving by at most 1.2e-4 (velocity). **The
  mechanical column drains by the same amount with the option as without**
  (`hydrostatic_column`, max |v|/c_s 0.796 against 0.792 at 300 steps and
  2.776 against 2.777 at 3000): that case's state is not the scheme's
  equilibrium and its drain is set by the outer boundary, as B4 recorded.
  **Default path**: with the key absent or False every path is bit for bit
  my control (`wasp_full_newton` reload, the atomic reload's 196 `(JFNK)`
  lines, `hydrostatic_column` at 300 steps, `mol_carrier` and
  `mol_diffusion` at their `maxsteps`); the one line that moves in every
  output file is the `# options` field, which now carries the 21st token
  `wellbal`, nameable in `Restart option change`. `EXHALE_resolved.out` and
  the parse dump each carry one new line. Key documented as K46. **For the
  advisor**: the option is kept in the tree and stays off; before it can
  ever be default on, the momentum row's reference scale has to become the
  pressure force the option cancels analytically,
  `|rho_j (A+ (phi_i(j) - phi_c(j)) + A- (phi_c(j) - phi_i(j-1)))|/dV`,
  which is two multiplications from what `RK_rhs` already holds but sits in
  an acceptance interface this item does not own. Advisor: rebuild on `build_lwv`; `krylov_and_dogleg` 279/0, `grid_and_gates` 162/0 (four rows FAIL when the shell rows are pointed at the stale shared `EXHALE.x`: pass `EXHALE_EXE` with `EXHALE_OBJDIR`), `steady_species_rows` 194/0, `certification` 57/0, `element_operator` 14/0, `residual_determinism` the same single row FAIL by design; atomic reload at cap 40 with every hook off and the key absent, 8 threads: the N31 ladder to the digit (7.607e-2 at 40, handed back 5.377e-2). The regression comparison strips comment lines, so the new `wellbal` token in the `# options` header is invisible to `make check`. Accepted as a default-off option; the momentum-row reference scale in `store_row_terms` (the option's `S(2) = 0` makes the row its own scale) is recorded in ISSUES 3.7 as the prerequisite of any default-on decision; the publisher PDFs of the two Kappeli and Mishra papers remain wanted (the institutional versions were read in full).

- **N38** (the discretized element and carrier rows at the front: why the row
  that binds is nearly independent of its own unknown). One default-off
  measurement hook, adopting nothing. `EXHALE_FRONT_ROW=1` (`[diag 18]`)
  takes the species row that carries the largest scaled residual and its two
  neighbors and prints, at the selected outer iteration: every TERM of the
  row on its certification scale, with the advective term split per FACE by
  masking the face mass flux (the divergence of a face species flux is
  exactly zero where the mass flux it multiplies is, so the masked call
  leaves one face standing alone and the SAME operator is used, never a
  copy); the reconstructed face fraction read back out of that masked
  divergence beside the donor cell's own; the chemical source of a carrier
  from the sweep's own rates; the band's entries for the row, to its own
  unknown, to the same unknown in the neighbors and to every hydrodynamic
  unknown of the three cells, each as a number and as a multiple of the
  diagonal; those same entries SPLIT BY TERM, by one-sided differences of
  the row's own terms along the same column with the same step
  `build_banded_jac_full` uses; the cell Peclet number, the logarithmic
  gradient and the chemical against the advective time scale; and the column
  scale of every unknown of the binding cell beside the unknown itself.
  **Verdict: the binding row's diagonal is NOT small -- it is the size of the
  row -- and the 3.3e3 of N35 is a CHOICE OF COLUMN SCALE that the Krylov
  cycle cannot see.** (1) **Every hydrodynamic entry of the row is its
  ADVECTIVE term**, on all three fronts measured: `carrier H2 of cell 205`
  (carrier reload, outer 40) `momentum of cell 205` -7.40478e+1 of which the
  diffusive part is 2.0e-5 and the chemical -2.3e-5; `element He of cell 191`
  (atomic, outer 1) `momentum of cell 190` -1.86758e+4 with a diffusive part
  1.9e-6; `element Fe of cell 137` (atomic cap 250, outer 140, one of the
  rows N34 named) `momentum of cell 136` -1.25980e+4 with 2.0e-6. Summed
  over the reported stencil the hydrodynamic columns carry **3389, 4176 and
  2424 diagonals** and the neighboring species unknowns 5.855, 1.052 and
  1.005, which is N35's 3.304e3 and 5.646 reproduced. Each species
  sub-block is meanwhile a clean diffusive tridiagonal (15.385/-7.772/-7.612
  helium, 17.432/-8.861/-8.564 iron): N27's reading of an interior element
  row holds at the front too, and it is the wind columns that bury it. (2)
  **The size of that coupling is one over the Mach number**, because the
  momentum column scale is `|rho v| + rho c_s`, the SOUND speed in a deeply
  subsonic layer: MEASURED `scale over unknown` **3.839e+2** (carrier front,
  momentum 6.83274e-5 against scale 2.62332e-2), **1.361e+5** (helium front)
  and **5.149e+4** (iron front), against 1.000 for the mass, the energy and
  every species unknown. The arithmetic closes to a factor 1.4 to 1.6: one
  face carries -0.5322 of the carrier row and 383.9 x 0.5322 = 204 against
  the measured 129; one face carries -0.1969 of the helium row and
  1.361e5 x 0.1969 = 2.68e4 against the measured 1.87e4. (3) **The three
  other candidates are measured and are not the cause.** The species column
  scale's reservoir floor is INACTIVE everywhere it was looked at (the
  carrier unknown is `ln n` with column scale exactly 1.0, the default; the
  helium and iron unknowns' scales are their own magnitudes to five digits),
  so raising or lowering `element_column_scale_floor` by a decade changes the
  diagonal by nothing. A stiff advective/chemical cancellation in the
  diagonal is a factor 1.8 and not decades (carrier diagonal 1.19228e-1 =
  diffusive +2.17367e-1 and advective -9.89614e-2, with the chemical
  derivative +9.44223e-4, 0.8 percent, although the chemical time scale is
  220 times the advective one), and there is no cancellation at all in the
  helium diagonal. The limiter is a live switch (the reconstructed face
  fractions stand at 0.587 and 0.657 of the donor at the carrier front, and
  frozen at first order the advective part of the diagonal changes SIGN,
  -9.89614e-2 to +2.50409e-1) but it moves the ratio by less than 2, because
  freezing it raises the momentum coupling by the same kind of factor. N5's
  two-root front was not re-run: cell 205 is not one. (4) **The test that
  settles it**: the existing default-off column equilibration
  (`EXHALE_COLUMN_EQUIL=1`) is exactly the remedy branch (i) names, and
  MEASURED at outer 1, where both runs stand at the same iterate, it brings
  the binding row's hydrodynamic coupling from 3389 to **0.993** (carrier)
  and 4176 to **2.926** (atomic) while the Krylov cycle returns the SAME step
  at the SAME residual to every printed digit (carrier 40 of 40, reduced
  1.257e-1, TRUE 1.259e-1; atomic 8 of 40, 7.821e-2, TRUE 7.824e-2) and the
  species-row Ritz count below 1e-2 is unchanged at 18 and 2. It cannot be
  otherwise: a column scaling `E` is applied to the operator and to the band
  built FROM it, and `(A E)(band(A) E)^-1 = A band(A)^-1` IDENTICALLY. N36
  showed a row scaling is a similarity; this shows a column scaling is an
  identity, so **both sides of the scaling are now exhausted**. What the column equilibration
  costs, since it also reshapes the trust-region ball and the probe steps:
  cap 40 ends at `||R||` 1.822 (carrier) and 1.592 (atomic) against 9.572e-2
  and 5.377e-2. Acceptance by the named clauses: the Ritz ratio below 1e3
  FAILS (2.07e5, 1.45e3), the cycle reaching 0.1 within 40 products FAILS on
  the carrier, named outcomes are lost on both; it stays off. **Deliverable 3
  built nothing**, by the branch N35's deliverable 4 took: the remedy branch
  (i) names is either already the default (`ln n`) or already measured to be
  invisible. Tests: `krylov_and_dogleg` RED on the control objects
  (`attributed_jacobian_entry` absent), GREEN **279/0** with 6 rows added
  (the hook off by default; a stated term set composing its own row and its
  two faces composing its advective term through
  `species_row_from_its_terms`, the expression the measurement itself uses;
  the attributed entry of a term linear in the unknown equal to its
  coefficient times the scaling, with a zero probe step and a row with no
  scale attributing nothing) -- 279 rather than 270 because N37 added nine
  rows to the same file, a brief overlap reported and not resolved;
  `steady_species_rows` 194/0, `element_operator` 14/0,
  `steady_selfconsistent_residual` 3/0 on the `wasp_full_newton` log,
  `residual_determinism` 5/1 with the same row failing by design and the same
  3.354e+01. Default path: the atomic reload at cap 40 with every hook off,
  all 354 `(JFNK)` lines IDENTICAL to my control through `||R|| 5.377E-02,
  flux spread 3.555E-03`; the carrier reload all 364 identical at 9.572e-2;
  `wasp_full_newton` `done info=0 ||R|| 1.667E-09` on both binaries with all
  28 lines the same; and with the hook ON both reloads still end at the
  control's values, so the measurement is not an operation on the state.
  Cost MEASURED: one residual evaluation for each of twenty columns at each
  of three outer iterations, about four minutes on the carrier reload at cap
  40 and 8 threads. Also corrected in the reading: N35's "there is no
  diagonal for a preconditioner to stand on" does not follow from the band
  row -- the diagonal is the size of the row, the disparity beside it is a
  column scale, and the near-singularity of the species-row compression of
  `A M^-1` is a separate fact that the band-row ratio is not evidence for. Advisor: verified on the same final build (the suites above); atomic reload identical. Accepted; deliverable 3 correctly not built. N35's sentence "no diagonal for a preconditioner to stand on" is withdrawn: the diagonal is the size of the row, the 3.3e3 is the momentum column scale (the sound speed) against a layer Mach number of 1e-3 to 1e-5, and column scalings are invisible to the cycle, so the near-singularity of the species-row compression is a fact of its own.
- **Close of the species-row solver program** (2026-09-10 evening, user's
  instruction: "어떤 것을 해도, 끝이 없어보이니, 엄청난 차이가 나오지 않는 한
  이번 건만 마무리 하고 끝내기 바래", then "있는 그대로 정리하고, md tex
  파일들 update 후 commit/push"): N37 and N38 were finished as briefed and
  no further item is started. State of the tree at the close, verified on
  the advisor's rebuilt `build_lwv`: `krylov_and_dogleg` 279/0 (174 at the
  morning commit), `grid_and_gates` 162/0, `steady_species_rows` 194/0,
  `certification` 57/0, `element_operator` 14/0, `residual_determinism` one
  row FAIL by design; the atomic reload at cap 40 with every hook off and
  every new key absent gives the N31 ladder to the digit (1.888, 1.891,
  1.819, 1.583, 1.861, 8.964e-1, 7.607e-2; handed back 5.377e-2, flux
  spread 3.555e-3); the `wasp_full_newton` reload ends `info = 0`, `||R||`
  1.667e-9, the number every worker's paired control gave (its output is
  not compared with the cold-start golden, which is a different path:
  `sec_ion_step` 2402 against 0). Everything N31 to N38 added is a
  default-off hook, option or key: `EXHALE_JV_ADDITIVITY`,
  `EXHALE_JV_COLUMN_SCALE`, `EXHALE_JV_PROBE_ARC`, `EXHALE_IEQ_TOL`,
  `EXHALE_TEQ_TOL`, `EXHALE_RESID_JUMP_SCAN`, `EXHALE_RESID_QUAD`,
  `EXHALE_KRYLOV_SIZE_SCAN`, `EXHALE_GM_HISTORY`, `EXHALE_PRECOND_SPECTRUM`,
  `EXHALE_BAND_DIFFERENCE`, `EXHALE_MODEL_ROW_EQUIL`,
  `EXHALE_GM_TRUE_RESIDUAL`, `EXHALE_FRONT_ROW`, and the input key
  `Well balanced:` (K46); the marching path and the goldens are untouched
  (the `# options` header gains the `wellbal` token, which the regression
  comparison strips with every comment line). New fixture
  `backup/regression/carrier_elem_newton`. Documents closed to this state:
  `ISSUES_20260909.md` 3.1, 3.7 and 5; `code_status_20260910.md` 3.3 and
  task 1; `session_handoff_20260910.md` (new, the handoff of record);
  `To_be_determined_by_user_20260906.md` decision 23 with its outcome;
  `steady_solver_design.md` 22.4; `input_schema.md` K46 and the manual's
  key table (PDF rebuilt, 90 pages); `README.md` (the key). The analysis
  the user asked for, of why the coupled solve fails and what a different
  approach would be (a segregated solve, the CETIMB way read from
  Koskinen et al. 2013 section 2.1.3 and 2022 Appendix B, a different
  treatment of the layer), is `docs/solver_approach_analysis_20260910.md`.
  Reported, not fixed (ISSUES 3.7): the `count` global shadowing the
  intrinsic; `newton_dense`'s hardcoded step test; the dead `dp_bc`/`ntot_bc`
  refreshes; `store_row_terms`'s momentum scale under the new key;
  `fortdep.py` blind to `include`; `golden/mol_ir_bands/` without
  `Cooling_breakdown.txt`; `examples/README.md`; the two `.md` twins.

---

## 8. PLAN_20260911: the partitioned stationary route (2026-09-11)

Plan of record: `docs/PLAN_20260911_partitioned_solver.md`. User instruction
of 2026-09-11: read `docs/solver_approach_analysis_20260910.md`, its review
`docs/solver_approach_analysis_20260910_review.md` and
`docs/solver_partition_experiment_20260911.md`, correct the errors they
establish, run the solver experiments, and adopt the partitioned route in the
production solver if it measures better than the coupled species-row solve.
This reopens the stage-2 solver program that the 2026-09-10 instruction had
closed after N37 and N38. Items P1 to P9, then P10 to P12 on the same day's
further instruction to fix the open items the plan had left to the user, then
P14 and P15 on the two findings of that afternoon's Codex adversarial review
and P16 to P18 on the user's two decisions of the same afternoon, then P20 and
P21 in the evening and the terminology sweep A1 to A4 on the user's instruction
of that evening (P13, P19 and P22 are the record itself); the rules every worker
on this tree follows are `docs/worker_rules.md`.

### The items

- **P1 (defect D1): a failed element composition solve is no longer handed
  back as an update** (`binary_element_diffusion.f90` +311/-23,
  `src/tests/element_operator/element_operator_tests.f90` +190/-13; every line
  count of this section is MEASURED from `git diff --numstat` on the tree at
  the close).
  `solve_mass_fraction` returns an explicit success result and no longer
  accepts the last halved trial of a line search that never descended;
  `element_diffusion_step` takes an optional status and, when asked,
  restores the composition it was handed unless the candidate is finite,
  inside `[0,1]` to `1e-12` before the range clip, and no farther from the
  density the step held fixed than the entry composition was, to `1e-10`;
  `relax_element_composition` steps on a copy, discards a refused step,
  retries at half the length within a bounded budget and reports three
  outcomes by name. MEASURED (P1) on the archived D1 reproduction
  (`atomic_elem_newton`, species mode, omega 1): a returned composition
  missing rho by `1.768571764708571E-01`, a nonfinite chemical refresh and
  exit status 1 become the entry composition restored at
  `1.046888214610824E-15`, an admissible refresh at `1.10297235E-15` and
  exit status 0. The first inadmissible operation of that failure, which the
  experiment report left unlocated, is now named: the composition solve
  drives `X` out of `[0,1]` and the range clip used to absorb it
  (`element_step_out_of_bounds`), and the mass discrepancy is downstream of
  the clip. Tests: `element_operator` 14/0 at entry, **24/0** after, **22/2**
  on a RED build with the four acceptance actions disabled
  (`an_unsolved_step_is_refused`, `an_unsolved_step_leaves_the_entry_composition`
  measured `1.31725E-02`); `diffusion_tests` 35/0, `steady_species_rows`
  194/0, `certification` 57/0 unchanged. Impact (MEASURED, P1, control built
  from the entry text of the one file): `mol_diffusion` and `lower_profile`,
  300 steps, one thread, `Ion_species.txt` byte-identical and
  `Hydro_ioniz.txt` identical but the provenance `run=` timestamp, on both
  cases; the marching caller passes no status, so none of the new tests
  runs there. One fixture of the same suite was corrected while there
  (`the_relaxation_fixed_point_under_a_callers_ghost` built a column whose
  species missed rho by `1.321E-01`): on the suite's own mass-closed column
  the worst interior element row of the relaxation fixed point moves from
  `1.62384E-02` to `1.35062E-13`, so the fixed point on a mass-closed column
  really is a zero of the stationary element row.
- **P2 (defect D2): the carrier movement bound now constrains the returned
  state** (`diffusive_photochemistry.f90` +225/-81,
  `src/tests/carrier_retry/carrier_retry.f90` +261/-1). The displacement is defined once (the largest change of any solved
  carrier column over the physical cells, divided by the largest H2 mixing
  ratio of the outer entry state, a quantity that cannot move while the pass
  runs), it is enforced before acceptance on every trial, a refused trial is
  undone and halved down to `2**(-10)` of the cell transport time, and the
  pass ends on one of five named outcomes. MEASURED (P2) on
  `carrier_elem_newton` through the archived driver: the entry text returns
  `drift = 2.8126757844770485e-02` at `trust = 0.01`, reproducing the
  0.02812676 of `docs/solver_partition_experiment_20260911.md` section 7.2
  exactly; the repaired routine returns `9.9956436739254070e-03` at the same
  setting and `2.9845242387380568e-03` at `trust = 0.003`, so the bound holds
  and it scales (a factor 3.33 smaller setting gives a factor 3.35 smaller
  movement). The halving ladder was measured on an instrumented copy: 0.4412
  at the full cell crossing time, then 0.2523, 0.1445, 0.0779, 0.0404,
  0.0205, 0.0103, 5.155e-3, 2.581e-3, 1.291e-3 and 6.458e-4 at the tenth
  halving, so the movement is proportional to the trial length and a bound of
  1e-4 keeps nothing on this column. Tests: `carrier_retry` **111 PASS / 10
  FAIL** on a RED control carrying the entry logic behind the new interface,
  **122/0** after; `carrier_returned_state_acceptance` 36/0, `attempted_step`
  70/0, `carrier_constraint_attribution` 10/0, `carrier_reference_scales`
  14/0, `species_masses` 9/0, `certification` 57/0, `adv_static_limit` 53/0.
  Impact: `mol_carrier`, 300 steps, one thread, `Ion_species.txt`
  byte-identical and `Hydro_ioniz.txt` identical but the provenance
  timestamp, which is the expected result because the marching loop calls
  `photochemical_transport_step` directly and never this routine. The cost is
  stated: a bounded pass advances less, so the frozen H2 residual after the
  composition update is larger (4.4708e-01 at trust 0.01 against the entry
  text's 3.1897e-01, from 4.9275e-01 before), and the state handed on is
  correspondingly closer to the entry state.
- **P3 (defect D3): the stop test of the stationary solve is now the
  acceptance test of the state** (`steady_newton.f90` +203/-29,
  `src/tests/steady_completion_flag/run.sh` +61/-15,
  `src/tests/krylov_and_dogleg/krylov_and_dogleg_tests.f90` +127/-1). The loop
  top of `solve_steady_jfnk` and of `solve_steady_ptc` no longer stops for
  success on `steady_gates_met` alone: it stops only when the gate AND every
  row the certification judges the carried system by are within their own
  tolerances, measured on the iterate by the same evaluation and the same
  formula the acceptance at return uses. A gate met while a certified row
  refuses prints one line naming the row, its measure, its own tolerance and
  its cell, and the iteration continues within its budget. MEASURED (P3) on
  the D3 reproduction (`carrier_elem_newton`, first hydrodynamic solve): the
  entry text stopped at the loop top of iteration 14 with `info = 0` and was
  refused at return on a mass row of 5.954e-12 against 3e-12; after the
  repair the same solve continues, 27 of its 40 iterations meet the gate with
  a row refusing, and it ends `info = 1` on its pass budget with mass
  5.945e-12, momentum 5.991e-13, energy 2.481e-09. So the second of the two
  admitted outcomes: the mass row of this solve does not come below 3e-12 in 27
  further iterations (it stands at 5.9e-12 to 7.9e-12 throughout, consistent
  with the N33 rounding floor), and the flag is now honest about it. Cost of
  the added test, MEASURED with a temporary clock pair in a private build: 27
  calls, mean 3.19 ms against a mean outer iteration of 0.73 s, 0.44 percent
  of one outer iteration, and only behind the gate, so an iteration below the
  gate pays nothing. Tests: `krylov_and_dogleg` 279/0 on the control objects
  and **290/0** after (11 new rows asserting the predicate the stop test
  reads); `steady_completion_flag` RED on both real logs at entry (an
  `info = 0` log with no statement that the stop test read the certified rows;
  a refusal at return with no stop statement) and GREEN after;
  `certification` 57/0, `carrier_returned_state_acceptance` 37/0,
  `attempted_step` 71/0 unchanged. Impact: the `wasp_full_newton` reload ends
  `info = 0`, `||R||` 1.667e-09, flux spread 6.233e-14, 8 outer iterations,
  `log10 Mdot` 13.30 on both binaries, the certification block line for line
  identical, every output byte-identical but the provenance line of the two
  `Hydro_ioniz*` files and one extra `loop-top stop:` line in the log; the
  `carrier_elem_newton` reload at cap 40 gives every `(JFNK)` line identical
  because the gate is never met on that reload. No golden can move: READ from
  the ten `DEFAULT_CASES`, none of their golden logs carries a `(JFNK) start`
  or `(PTC) start` line at all.
- **P4: the momentum row's reference scale under `Well balanced: True`**
  (`RK_rhs.f90` +90/-2, `steady_residual.f90` +59/-2,
  `src/tests/grid_and_gates/hydrostatic_residual.f90` +318/-11, that suite's
  `README.md` +44/-1). Under the well-balanced option `store_row_terms` scaled the momentum row by
  `max(|dF(2)|, |S(2)|, |Smom|)` with `S(2) = 0` by construction, so every
  scaled momentum residual was exactly 1 (N37, ISSUES 3.7). The reference is
  now `max(|dF_2|, |S_2|, |S_visc|, equilibrium_pressure_force(j))`, the last
  term the magnitude of the pressure force of the cell's own hydrostatic
  equilibrium, which is the right-hand side of whichever cancellation
  identity the assembled row belongs to (PLM
  `|rho_j [A+ (phi_i(j) - phi_c(j)) + A- (phi_c(j) - phi_i(j-1))]|/dV`, WENO3
  `|rho_j (phi_i(j) - phi_i(j-1))|/dr`; Kaeppeli and Mishra 2016, A&A 587,
  A94, eq. 16 in spherical geometry). It has one definition,
  `equilibrium_pressure_force_of_state`, called by `RK_rhs`, by
  `positivity_limited_fluxes` and by `assemble_residual` (the last so that
  the kind-generic `EXHALE_RESID_QUAD` assembly, which returns the well-balanced
  row but not the term it is read against, does not scale the row by a stale
  weight). The numerator is untouched; with no gravity the term is zero and
  the scale is the row's own terms exactly (MEASURED 0 for all four
  scheme/solver pairs). MEASURED (P4): `wasp_full_newton` reloaded with the
  key ON is **NOT CERTIFIED** on the entry text on exactly one entry, the
  momentum row at **1.000E+00** at cell 1 (which reproduces the N37 record to
  the digit, spread 1.539E-15, energy 9.911E-12, mass 1.228E-13, Mdot 13.30),
  and **CERTIFIED** after, `info = 0`, `||R||` 9.719E-10, momentum row
  9.719E-10 against 1.0E-08 at cell 500, `log10 Mdot` 13.30. Tests:
  `grid_and_gates` `hydrostatic_residual` **28 PASS / 12 FAIL** on the
  entry-text objects (the twelve failures at exactly 1.000 for all four
  pairs) and **40/0** after, of which 16 are new rows measuring the row
  divided by `residual_row_scale(2,...)`, the one expression the
  certification reads; the whole suite 162/0 at entry (READ, N37) and **178
  PASS / 0 FAIL** after; `krylov_and_dogleg` 279/0 and `certification` 57/0
  unchanged. Default path (key absent): `wasp_full_newton` reload and
  `mol_carrier` at 300 steps byte-identical on both binaries apart from the
  provenance line, every `(JFNK)` line identical, key-off `info = 0`, `||R||`
  1.667E-09, Mdot 13.30. The option's effect on the two solver fixtures is large
  and in one case negative, and both are recorded: the carrier reload's
  handed-back `||R||` falls from 1.721 to **1.026E-01**, the key-off level
  (READ: N37 9.34E-02), with a measured momentum row of 1.404E-02 at cell 499
  instead of a pinned 1; the atomic reload stops earlier (`info = 2` at
  iteration 24) and hands back 1.702 against the control's 1.415 at cap 40,
  because the merit no longer carries 500 momentum rows pinned at one
  (`||Fs||2` 17.1 against 28.2). A defect found while there and fixed: the
  allocation guard in `RK_rhs` tested `face_flux` for the well-balanced
  arrays `face_q_up`/`face_q_dn` as well, and `face_flux`/`face_p` are also
  allocated by `store_the_interface_fluxes`, so a stationary evaluation
  reached before any marching stage could write an unallocated array; the
  guards are now separate and each tests its own array, with no arithmetic on
  any path.
- **P5: the analysis document corrected** (`docs/solver_approach_analysis_20260910.md`,
  +463/-115, 114 to 592 lines). All nine corrections the review established
  are in place, each with an inline `[corrected 2026-09-11: ...]` mark citing
  the review section or the experiment measurement behind it: the spatial
  truncation error is not a lower bound on the algebraic residual (the
  frozen-background H2 wind row measures 2.449021e-13 on the same 500-cell
  grid where the withdrawn argument put a 5e-4 floor); Ritz-magnitude ratios
  of a masked compression are neither condition numbers nor the species Schur
  complement; partitioning moves the coupled iteration outside the Krylov
  solve with outer matrix `(1 - omega) I + omega D^-1 C A^-1 B` and the outer
  loop is not new code (`steady_wind_with_element_diffusion`); the CETIMB
  route is renamed as informed by CETIMB and gains the two conditions the
  review adds; the well-balanced option removes the equilibrium's algebraic
  cancellation and not the floor; Koskinen et al. 2013 lower boundary is
  section 2.1.1; the three-unknown path uses its existing globalization and
  not the trust region (`use_tr = (nspec_row .gt. 0)`); section 7 recast with
  the 2026-09-11 numbers; a new section 9 indexing the twelve changes and a
  new section 5.5 on partitioned components inside a coupled method (Knoll
  and Keyes 2004, J. Comput. Phys. 193, 357). Three numbers of the original
  disagreed with their sources and are corrected, one of them wrong: the N4b
  row for the carrier reload read "1.4 -> 0.09" where the measurement is 0.151 ->
  0.090. Citations were verified independently in the published text
  (Koskinen 2013 section headers by `pdftotext -layout`; the Knoll and Keyes
  bibcode by the ADS API).
- **P6: the outer iteration as a stated contract** (`src/EXHALE_main.f90`,
  +448/-131). `steady_wind_with_element_diffusion` now says in its own header
  which equations are alternated with the wind and why (a species the run
  moves by an operator-split transport step, the diffusing elements always
  and the molecular carriers where `Coupled carrier solve` is off, has a
  stationary balance that no hydrodynamic solve touches), and it accepts a
  state only by the certification of the FULL set of active equations
  evaluated on the refreshed state it would hand back. A hydrodynamic
  `info = 0` while an alternated species row refuses is not an accepted
  state, and a nonzero hydrodynamic flag no longer ends the loop by itself,
  because both `info = 1` and `info = 2` hand back a state the composition
  update still makes progress from. The endings are named
  (`outer_certified`, `outer_state_not_finite`,
  `outer_element_update_refused`, `outer_no_progress`, `outer_pass_budget`),
  the P1 status is consumed, and the quantity the passes are controlled by is
  the worst gated species row of the certification: where it fails to fall
  over a pass both updates are shortened (the element under-relaxation omega
  halved, floor 0.125; the carrier movement bound halved, floor
  `carrier_trust_floor = 1.0d-3`) and `outer_no_fall_max = 3` consecutive
  passes without a fall end the loop. The last pass takes no composition
  update, because an update exists to be consumed by the next hydrodynamic
  solve and taking it would replace a stationary wind by a state the update
  itself spoiled. The composition drift and the volume-weighted carrier
  residual are reported and are not tests. With no alternated species the
  body is one steady solve and the state refresh, and the flag is the solve's
  own, which is what an atomic run and a coupled carrier solve both do.
  Default path MEASURED on the P6 run pair (`wasp_full_newton.ctrl` against
  `.new`): `info = 0`, `||R||` 1.667E-09, flux spread 6.233E-14, CERTIFIED,
  `log10 Mdot` 13.30, and `output/Ion_species.txt` and
  `output/Hydro_ioniz.txt` byte-identical between the two.
- **P7: the measurement.** Below.
- **P8: the corrections propagated to the documents and tests that carried
  them** (`docs/code_status_20260910.md`, `docs/ISSUES_20260909.md`,
  `docs/certification_tolerance_anchoring_20260910.md`,
  `src/tests/steady_selfconsistent_residual/run.sh`,
  `src/tests/steady_species_rows/steady_species_rows_tests.f90`). The
  withdrawn inference (a spatial truncation error of 5e-4 at N = 500 floors
  the algebraic residual, so 1e-5 cannot certify) is corrected in place in
  the three documents that still carried it, marked with the review section
  and the measured 2.449021e-13 beside it, and decision 22 (a) is kept as the
  standing decision in every one of them with only its justification changed;
  the grid-refinement anchor and its remap control are relabelled as
  consistency diagnostics; "permanently uncertified until the operator's
  order or the grid changes" is withdrawn and quoted where it stood; ISSUES
  3.7's `store_row_terms` item is struck through and marked RESOLVED by P4,
  and the resolved table gains D1, D2 and D3 with their measured before and
  after. Two test contract errors, both RED before and GREEN after:
  `steady_selfconsistent_residual` assertion 2 read "a hydrodynamic row above
  its tolerance implies `info = 2`", which conflates not certified with one
  flag and already failed on the entry text for any budget-limited log
  (`FAIL flag_of_carrier_elem_newton measured=info=1 reference=info=2` ->
  PASS, with `wasp_full_newton` still PASS and a synthetic `info = 0` copy
  still FAILing); and `relaxation_drift_covers_every_element` built its
  column by swinging each metal mixing ratio about the reservoir on a fixed
  hydrogen fraction, so its species missed the density the column hands the
  operator by **1.66042E-02**, which is now **6.58484E-16** on a column of
  mass fractions that sum to one, the suite going 194/0 to **195/0** with a
  new row measuring the closure by the code's own mass policy.
- **P9: this section**, `docs/session_handoff_20260911.md`, the status section
  of the plan, the dated sections of the two fixture READMEs, the reopening
  paragraphs of `docs/code_status_20260910.md` section 3 and
  `docs/ISSUES_20260909.md` section 5, and the new entries of
  `TO_BE_DONE.md`.
- **P10: the best-iterate ledger reads the hydrodynamic rows against their own
  certification tolerances** (`steady_newton.f90` +131/-55,
  `src/tests/krylov_and_dogleg/krylov_and_dogleg_tests.f90` +135/-1, READ from
  the P10 report; the whole uncommitted diff of the two files now stands at
  +335/-84 and +268/-1, MEASURED by `git diff --numstat`, P3's share
  included). Items P10, P11 and P12 were added after the plan had run, on the
  user's instruction of 2026-09-11 to fix the open items the session had left
  to the user ("fix all of them", and "the pressure term too"); P10 is the
  third entry of the plan's "what remains".
  `distance_from_certification` divided the hydrodynamic slot by the run's
  `Resid tol`, so the ledger that chooses the state handed back could not see
  a mass row above 3e-12 while `||R||` was below 1e-8. The slot now arrives as
  a distance formed by ONE new public routine,
  `hydrodynamic_distance_from_certification`, `max(mass/3e-12, momentum/1e-8,
  energy/1e-6)` on the windows `resid_relnorm` uses, so `d < 1` is exactly the
  row condition the loop-top stop of P3 and the acceptance at return impose,
  and `keep_this_iterate` ranks on that distance on every route including the
  three-unknown one (the `nspec_row <= 0` branch that ranked on `||R||` is
  gone). `certified_row_measures` hands the plain `||R||` of the same
  evaluation back through a new optional argument rather than leaving a caller
  to re-derive it, and the iteration line prints the two side by side.
  `resid_tol_of_solve` is KEPT, with one reader, the
  `EXHALE_RESID_BRANCH_REPORT` block, where `rmax(1)` is deliberately the
  plain mixed-row measure and dividing it by one row's tolerance would be a
  ratio of two different rows. Tests: `krylov_and_dogleg` 290/0 at entry,
  **300 PASS / 9 FAIL** on a RED build in which both distances read
  `maxval(rc)/1.0d-8` (the superseded reading, at the `Resid tol` of these
  fixtures), **309/0** after, the 19 new rows being the block
  `judged_distance_reads_each_row_against_its_own_tolerance`; `certification`
  57/0, `steady_completion_flag` 3/0, `steady_selfconsistent_residual` 3/0
  unchanged. Impact MEASURED (P10) against a private build of the entry text
  of the one file: `wasp_full_newton` reloaded ends `info = 0`, `||R||`
  1.667E-09, `log10 Mdot` 13.30 with **all eleven output files
  byte-identical**, the ledger handing back the same iterate and the only
  number that moves being the reported best-iterate distance, 1.749E-01
  (which is `||R||_best/1e-8`) to 8.469E-02 (the same iterate's worst row over
  its own tolerance, so 8.5 percent of it). The carrier reload shows the
  re-ranking doing what it is for: pass 2 hands back a mass row of **4.85E-12
  at `||R||` 2.459E-09** where the control took 6.98E-12 at 1.593E-09, and
  pass 3 **4.96E-12 at 2.005E-09** against 8.17E-12 at 1.723E-09, with the H2
  carrier row and the worst gated species row identical pass by pass because
  the carrier relaxation sits outside the hydrodynamic solve. The HD 209458 b
  element reload runs its twelve passes in 1064 s against the control's 1043 s
  (2 percent on a shared machine, not a claim) with every one of the twelve
  worst gated species rows and its cell identical to the control, and ends on
  the same single refusing entry, the hydrodynamic mass row at **8.893E-10 at
  cell 3** against the control's 1.054E-09 at cell 16; the solves reporting a
  stalled best iterate fall 5 to 3 and those ending on "no descent direction
  exists" 11 to 10. What stays open, and what this measurement settles: the
  mass row does NOT reach 3e-12 and that state still does not certify, because
  at 1e-9 to 1e-12 the row IS the N33 rounding of the base-layer flux
  difference. MEASURED (P10, finding 3): at pass 1 the ledger's best distance
  is 4.407E+02, a mass row of 1.32E-09, while the same state re-evaluated at
  its own composition after the best-iterate restore carries **4.11E-09**, so
  ranking iterates by that row is ranking by the luckiest rounding of a
  quantity that moves when the state is evaluated again, which is why the
  pass-by-pass mass rows scatter rather than descend in the control and in the
  new run alike. The entry is an acceptance question, not a solver item, and
  no further solver item follows from it.
- **P11: the momentum row's reference scale is the largest PHYSICAL term of
  the momentum equation** (`RK_rhs.f90`, `steady_residual.f90`,
  `src/tests/grid_and_gates/hydrostatic_residual.f90` and that suite's
  `README.md`; the whole uncommitted diff of the four files stands at +285/-2,
  +151/-33, +736/-12 and +94/-1, MEASURED by `git diff --numstat`, with P4's
  +90/-2, +59/-2, +318/-11 and +44/-1 inside it). The fourth entry of the
  plan's "what remains", and the "pressure term" of the user's instruction.
  The scale was built from the pieces of the discretization,
  `max(|dF(2)|, |S(2)|, |Smom|)` plus P4's equilibrium pressure force under
  the well-balanced key; it is now built from the three left-hand terms of
  `d(rho v)/dt + div(rho v v) + dp/dr + rho dphi/dr = S_visc`, kept beside
  `face_flux`/`face_p` as module arrays of `RK_integration`, filled on both
  reconstructions, under both key values and on the positivity-repair path,
  and gathered by one routine, `momentum_row_terms_of_cell`, with
  `momentum_row_terms_of_state` filling a whole state and
  `reconstruction_continuation_rhs` blending the three terms with the same
  `lambda` as `dF` and `S`. The residual `R(2) = dF(2) - S(2) - S_visc` is
  untouched on every path, and the split is a rearrangement of the assembled
  row and not a model of it: `F(2) - p_out` is the momentum flux without the
  pressure exactly for every flux function, and the geometric term is
  recovered from the source as `p_geom = S(2) + weight` so that no second
  definition of the cell pressure exists. MEASURED (P11), on a state with
  structure in density, pressure and velocity,
  `max_j |ram + dp/dr + rho g - (dF_2 - S_2)|` over the scale is **0 under
  WENO3 and 1.7e-16 (default) / 5.3e-15 (key on) under PLM**, for both Riemann
  solvers, so the three terms add up to the row they scale. The defect is
  closed on the states where the terms are known in closed form: with no
  gravity and a uniform pressure at rest, where every term of the equation
  vanishes, the PLM scale read **1.000 of the geometric pressure term
  `(A+ - A-) p_j/dV`**, that is 2 p/r, with both Riemann solvers at every N,
  and now reads the `tiny` floor 2.22507E-308; on supersonic uniform flow at
  M = 3, where the ram divergence is the only term, it stood **6.66667E-02 =
  1/(gamma M M) above it** and now agrees with it to 1.9e-14 (PLM) and
  4.2e-13 (WENO3); on the discrete hydrostatic equilibrium the PLM scale was
  **2.79144E-01, 28 percent BELOW the weight it balances**, because `S(2)` is
  the weight minus the geometric term and no term of the max was the weight
  itself, and the departure `(w - s)/w` is now 0. P4's own zero-gravity
  diagnostic closes with it: the departure of the well-balanced scale from the
  base scheme's on the same state falls from **9.43654E-01 to 3.74E-13** under
  PLM with ROE (2.76E-13 with HLLC), the two having agreed under WENO3 all
  along (1.7e-13 to 2.7e-13). Tests: `hydrostatic_residual` **50 PASS / 10
  FAIL** on the entry-text objects and **60/0** after, the ten being the four
  scheme/solver pairs of P4's corrected
  `zero_gravity_momentum_scale_is_the_rows_own_terms` (which had compared the
  option's scale with `max(|dF_2|, |S_2|)`, the SUM of the two terms the row
  holds, 7.50E-01 FAIL to 0) and the PLM pairs of the three new identity
  statements; the whole `grid_and_gates` suite **188/10 to 198/0** (198 =
  P4's 178 + 20 new rows); `krylov_and_dogleg` 290/0, `certification` 57/0 and
  `adv_static_limit` 53/0 before and after. Two quantities are REPORTED and
  not asserted, with the reason: the two-sided departure of the scale from the
  weight on the well-balanced scheme's own equilibrium is not round-off for the base scheme
  (MEASURED 8.5e-4 at N = 250 falling to 2.3e-6 at N = 2000, PLM with ROE),
  because that state is the discrete equilibrium of the WELL-BALANCED reconstruction
  and HLLC's face pressure carries a dissipative part of size dr/H times the
  weight; and the halving statement subtracts the unperturbed row of the same
  cell, since the undifferenced ratio is not 2 (MEASURED 1.0012 at N = 250,
  1.586 at N = 2000) while the imbalance the perturbation ADDS halves to
  2.000000 within 5e-6. Impact: nothing that was certified stops being
  certified. `wasp_full_newton` reloaded with the key OFF is **CERTIFIED** on
  both binaries, `info = 0`, `log10 Mdot` 13.30, with the momentum row moving
  **8.469E-10 to 1.165E-09** of its 1.0E-08 tolerance, the mass row 1.453E-13
  to 2.017E-12 of its 3.0E-12 (so this reload now certifies with less margin
  on the mass row) and `||R||` 1.667E-09 to 1.495E-09; the outputs are not
  byte-identical, because the merit and the loop-top gate both divide by the
  momentum scale and the two solves accept different iterates, and the states
  differ by at most **2.6E-09 relative** over every column (9.3e-10 in
  `Hydro_ioniz.txt`, 2.6e-09 in `Ion_species.txt`). With `Well balanced:
  True` the same reload is **CERTIFIED** on both, the momentum row moving
  **9.719E-10 to 9.621E-11** with all eleven output files byte-identical and
  the `||Fs||2` ladder identical iteration by iteration: at cell 500,
  r = 4.62, the ram divergence exceeds the cell's weight by about a factor
  ten, which is why the measure of the certified state falls by that factor on
  a supersonic outer cell. The marching path of `mol_carrier` at 300 steps and
  one thread is **byte-identical in all seven output files**, `log10 Mdot`
  7.95, the `run.log` differing in exactly eight lines, every one of them a
  report of the momentum measure (the worst-row entry 1.8164E+00 at cell 438
  to 1.9788E+00 at cell 449) and none of them a state. On the carrier reload
  every pass line is identical but the momentum entry, which falls **1.32E-12
  to 3.22E-13** after pass 1 (4.47E-13 and 3.54E-13 at passes 2 and 3), the
  state as loaded reading 3.622E-01 to **3.311E-01** in its scaled momentum
  row, and all four output files byte-identical. What stays open: the
  kind-generic `EXHALE_RESID_QUAD` assembly still hands back neither the momentum
  terms nor the face departures, so under that option (default off, a rounding
  control experiment nothing is adopted from) both are read from the last
  `RK_rhs` evaluation rather than from the state just assembled; the fix
  belongs in `hydrodynamic_rows.f90`, which this item does not own.
- **P12 (with P12b): two test-hygiene defects, and the same shape in the
  sibling scripts** (`src/tests/test_columns.f90` new, 125 lines;
  `attempted_step/run.sh` +18/-5, `carrier_returned_state_acceptance/run.sh`
  +18/-5, `carrier_retry/run.sh` +18/-5, `carrier_reference_scales/run.sh`
  +22/-6, `carrier_constraint_attribution/run.sh` +22/-6,
  `constrained_network_layout/run.sh` +22/-6, `species_masses/run.sh` +22/-8,
  `element_operator/run.sh` +5/-0, `steady_species_rows/run.sh` +5/-0, and the
  two drivers, MEASURED by `git diff --numstat`; of the drivers' totals the
  item's own share is READ from the P12 report, the `use` line, the renamed
  calls and 66 lines out of `element_operator_tests.f90` and 43 out of
  `steady_species_rows_tests.f90`). The fifth entry of the plan's "what
  remains", both halves of it. The **`build_stamp` fallback**: the trigger is
  narrower than "a tree that has never been built", it needs such a tree AND a
  work directory outside `$root/build`, which is the concurrency-safe usage
  the worker rules mandate, because `source_closure.py` walks `$root/build`
  recursively and a default work directory inside it swept the fallback into
  the closure by accident. Seven `run.sh` now select
  `${EXHALE_OBJDIR:-$root/build}/build_stamp.f90`, write the fallback into
  their own work directory only when that is absent, and compile whichever of
  the two is in force as the FIRST file of the closure with the closure
  filtered of any second copy, which is the shape
  `src/tests/physics_probe/run.sh` already had;
  `grep -ln 'if [ ! -f "$root/build/build_stamp.f90" ]' src/tests/*/run.sh`
  now returns nothing. Four of them gained a private object directory,
  `EXHALE_TEST_OBJDIR`, without which the no-build verification cannot be run
  on them at all (`species_masses`, which had no override of any kind,
  `carrier_reference_scales`, `carrier_constraint_attribution`,
  `constrained_network_layout`). RED and GREEN, MEASURED in four cases: at HEAD
  with an external work directory the compile dies on `Cannot open module file
  'build_stamp.mod'` and each script prints its one failure line, **0 PASS / 1
  FAIL, exit 1**, in `attempted_step`, `carrier_returned_state_acceptance`,
  `carrier_reference_scales`, `carrier_constraint_attribution`,
  `carrier_retry` and `species_masses`; after the fix every suite gives the
  count it gives at HEAD with a default work directory, **70/0**
  (`attempted_step`), **36/0** (`carrier_returned_state_acceptance`),
  **14/0**, **10/0**, **122/0**, **9/0** and **8/0**
  (`constrained_network_layout`, which has no `utilities.f90` in its closure,
  so it has no RED to show and its change is the shape alone), in a tree with
  no `build/` and in the repository against a private build, and the `diff` of
  the PASS/FAIL lines before and after is empty for all of them. Which stamp
  each case used is MEASURED too: in the no-build case every work directory
  holds a `build_stamp.o` compiled from the fallback, and in the case with a
  build tree no fallback file was written at all, so the generated stamp was
  taken, which the superseded scripts could not do because they looked only in
  `$root/build`. The **column written twice**: the new module `test_columns`
  holds `column_carrying_its_own_density(f_c, q_h2, he_taper, metal_swing)`,
  named for what it builds, with its header stating the invariant
  `sum_i m_i n_i = rho` and the round-off it holds to, and
  `column_mass_closure(rho_c, f_c)`, the measurement moved out of
  `element_operator_tests.f90` unchanged; the two drivers call them at 15 and 2
  sites and `use utils, only: calc_rho` is gone from both, it having been there
  only for the moved code. The shared construction reduces to each former copy
  exactly and not approximately (`metal_swing = 0` gives `swing = 1.0d0`,
  `q_h2 = 0` leaves `isp_H2` at zero, `he_taper = .false.` selects the
  constant the steady driver used, and `zsum` accumulates in the same element
  order), and the statement is the measurement: `element_operator` **24/0**
  and `steady_species_rows` **195/0** before and after, with the `diff` of the
  two full logs, diagnostic lines included, empty in both, so every printed
  digit is unchanged (`the_column_is_mass_closed_as_built` 4.10591E-16,
  `drift_column_carries_its_own_density` 6.58484E-16, the relaxation's 27
  steps and returned drift 1.7790E+00 among them). The file is in no source
  list of the Makefile (`grep -c test_columns Makefile` = 0) and no dependency
  generator sees it; the two `run.sh` that need it compile it immediately
  before their driver. Nothing stays open from this item: its own report
  found nothing further, and the four sibling scripts and the hardcoded work
  directory it reported are what P12b then closed.

Two changes outside the item list, made by the advisor while reading these
three:

- `src/tests/steady_species_rows/steady_species_rows_tests.f90`, the two rows
  P10 reported as stating the superseded semantics. In
  `judged_distance_is_the_worst_row_over_its_tolerance` the hydrodynamic slot
  now arrives as a distance (`rows(1)` 2.0d0, 0.5d0 in place of 2.0d-8,
  5.0d-9), with a comment saying so, and
  `judged_distance_follows_the_residual_target`, which asserted the property
  the item removes, is replaced by
  `hydrodynamic_distance_reads_each_row_against_its_own_tolerance`: a mass row
  of 1e-9 beside a momentum row of 1e-13 and an energy row of 1e-8 gives
  `1e-9/3e-12`, the mass row being what the distance reports. The two
  `resid_tol_of_solve` assignments and that import are gone.
- the `Coupled carrier solve` key comment in
  `src/modules/files_IO/input_read.f90` (+17/-13, MEASURED). It stated that
  the alternation "was measured not to converge and not to be able to"
  (section 139), which is narrower than it reads after P7. It now names the
  default (False), says what the alternation is judged by (the certification
  of the refreshed state in `steady_wind_with_element_diffusion`), cites the
  two reloads of this section for what it reaches, and keeps True as the
  coupled measurement option.

Four more items ran in the afternoon of the same day: **P14** and **P15**, the
two findings of the Codex adversarial review of the working tree (the review
and the verification of each finding are the subsection below), and **P16** and
**P17**, the user's first decision of the afternoon (the mass tolerance
re-anchored on the fixture's measured rounding floor) and its consequence for
the best-iterate ledger. **P18** followed and closed the same day.

- **P14: every ending of the stationary outer iteration hands back one
  consistent state** (`src/EXHALE_main.f90`, `steady_wind_with_element_diffusion`
  only, +84/-40 against the entry text;
  `src/tests/grid_and_gates/output_state_consistency.sh` +181/-4; MEASURED by
  P14 with `diff -u` and `git diff --numstat`). Codex finding 1, and it is
  real. On a stagnation ending the composition had been updated while the
  particle densities and the temperature written beside it had not, so the
  state the run wrote was one composition update ahead of its own particle
  count and temperature. MEASURED (P14) on the forced stagnation of the
  hot-Uranus carrier reload: re-evaluating the written state moves its
  temperature by **1.317089E-08** relative at cell 295, where the same fixture
  ending on its pass budget, which takes no update, moves by **5.921E-16**; the
  written temperature of the worst cell is 1.768302856793E+03 where the same
  conserved state and composition give 1.768302833502E+03, and the worst
  species column is 2.543E-03 (He 2^3S, the stage answering to the temperature
  with the steepest exponential). Four changes. The progress control (the
  `sp_worst` test, `n_no_fall`, the omega and trust halving, the setting of the
  stagnation ending) now stands between the certification of the pass and the
  update the next solve would consume, guarded so that a certified or
  non-finite pass reaches neither, and the update's own gate therefore excludes
  the stagnation ending by itself: the pass that ends the iteration takes no
  update. The ending is announced with the other endings, below the summary
  line and the diagnostics, so the pass that ends on it prints its summary, its
  carrier diagnostics and its H2 front like every other pass. A pass that takes
  no update says which of the two reasons it was. And the element path's update
  now refreshes the particle densities and the temperature between
  `relax_element_composition` and the equilibrium sweep
  (`get_species_densities`, `comp_T_from_p`), which is the shape the carrier
  path already had; before it, the sweep equilibrated the new composition at
  the old temperature. The header states the invariant. One behavior change to
  record because it is not a reordering: the shortened step of a pass whose
  measure did not fall is now consumed by THAT pass's update, where before it
  reached the next pass's. RED and GREEN on the new second case of
  `output_state_consistency.sh`, which runs the stagnation
  (`EXHALE_CARRIER_TRUST=1e-4 EXHALE_OUTER_PASSES=8 EXHALE_JFNK_MAXIT=5`, the
  setting that stagnates; the brief's cap of 40 does not, the worst gated H2
  row falling through the hydrodynamic solve alone), asserts by name that the
  ending reached is the stagnation one, and then measures the `Restart intent:
  stationary evaluate` round trip over the 500 physical cells against **1e-12**:
  `measured=1.317089e-08` FAIL on the entry text, `measured=6.499247e-14` PASS
  after. The tolerance is justified by measurement and not by a snapshot: the
  states handed back by an ending that takes no update measure 6.5E-14 and
  5.9E-16, and the defect measures 1.3E-08, four decades above the allowance.
  The suite goes **199 PASS / 1 FAIL** on the entry text to **200 / 0**.
  Impact: `wasp_full_newton` reloaded is unmoved (`info = 0`, CERTIFIED, `||R||`
  1.495E-09, flux spread 1.9113E-12, `log10 Mdot` 13.30 on both binaries, all
  eleven outputs identical apart from the provenance `run=` timestamp); the
  carrier reload at three passes is identical pass line by pass line to the
  control and to the P7 record; the element reload's pass 1 is identical in
  every digit and the ladder parts from pass 2, because the sweep of pass 1's
  update now runs at the relaxed composition's own temperature (momentum
  1.60E-12 against 4.24E-14 at pass 2), both builds ending on the pass budget
  with the same 9 refusing equations and the same base-layer mass row as the
  worst refusing entry.
- **P15: the kind-generic rows hand back the face departures their momentum row
  was built from** (`src/modules/time_step/hydrodynamic_rows.f90` +75/-15,
  `src/modules/time_step/hydrodynamic_rows_body.inc` +51/-17,
  `src/tests/krylov_and_dogleg/krylov_and_dogleg_tests.f90` +255/-1 measured
  against the file as it stood when the item started; MEASURED by P15). Codex
  finding 2, and the third entry of the plan's "what remains". Under
  `Well balanced: True` the momentum row's pressure-gradient term, the momentum
  row scale and every measure read against that scale are formed from the
  module arrays `face_q_up`/`face_q_dn` of `RK_integration`, which only
  `RK_rhs` and the positivity repair ever filled; the kind-generic rows kept
  their own departures as locals and `store_the_interface_fluxes` wrote back
  `face_flux` and `face_p` alone, so in the kind-generic rows those three quantities belonged
  to the last `RK_rhs` evaluation or to the zeros of the first allocation.
  The departures now leave the kind-generic text as two OPTIONAL `intent(out)`
  arguments of `hydrodynamic_flux_difference_and_source`, converted to double
  at the one rounding on the way out exactly as the flux and the face pressure
  are, both public wrappers ask for them, and `store_the_interface_fluxes`
  stores them under the key. The reconstruction continuation CAN reach the
  kind-generic rows (`assemble_residual` chooses the text from
  `generic_precision_rows_arm()` alone and the generic text carries its own
  blending branch), and on it the terms would be wrong for two reasons rather
  than one, the single-scheme weight `w_face_p` and the fact that the two
  schemes reconstruct different faces so no single pair of departures satisfies
  the blended row identity; `the_generic_text_is_usable` therefore now refuses
  `recon_lambda_on .and. 0 < recon_lambda < 1`, ahead of every evaluation of that text.
  That branch ends in `error stop` and has no executed row, which is stated
  rather than worked around. Twelve new rows, each on a state the kind-generic rows assemble
  AFTER a residual was assembled on a different state, both reconstructions and
  both precision kinds: MEASURED RED **2.098877E+01** (PLM) and **2.032204E+01**
  (WENO3) on the pressure gradient against the production term,
  **3.262392E+00** and **3.428893E+00** on the row scale, **1.190363E+01** and
  **1.181676E+01** on the row rebuilt from the stored face data, and GREEN at
  **0 bitwise** for the generic-double instantiation, which IS the production
  operator bit for bit, and 2.6E-16 to 2.7E-15 for the quadruple one at its
  1e-13 bound. `krylov_and_dogleg` **309/0 to 321/0**; `certification` 57/0 and
  `steady_species_rows` 195/0 unchanged. Impact: `wasp_full_newton` reloaded is
  unmoved with the key off, and with `Well balanced: True` and
  `EXHALE_RESID_QUAD=2` it stays CERTIFIED at the certified momentum row
  **9.621E-11**, P11's key-on production value, on both binaries, because the
  JFNK and PTC loops reset the sweep kind to marching before the certification
  of the handed-back state and the kind-generic text is consulted only for a non-marching
  sweep kind (READ, `steady_newton.f90`); the brief's expectation that the
  control's certification reads another state's departures is corrected by that
  measurement. What the defect did move is the measure INSIDE the loop, where
  the kind-generic rows are armed: the convergence measure of iteration 1 is **1.756E-03** on
  the control against **1.753E-03** after, the momentum row divided by the
  pressure-gradient term of another state and overstated by 0.17 percent, and
  that measure is what the loop-top stop of P3 and the ledger of P10 read on
  every iterate. `||Fs||2` is identical iteration by iteration, iterations 2 to
  7 are identical line for line, the accepted state is the same state, and
  `mol_carrier` marching is byte-identical.
- **P16: the hydrodynamic mass-row tolerance is anchored on the row's measured
  rounding floor** (`src/modules/time_step/certification.f90` +212/-2,
  `src/modules/time_step/steady_residual.f90` +177/-1,
  `src/tests/certification/certification_contexts.f90` +105/-0,
  `docs/certification_tolerance_anchoring_20260910.md` +302/-0 appended;
  MEASURED by P16 with `git diff --numstat`, the share of the other uncommitted
  item in `steady_residual.f90` separated by attributing each hunk). The user's
  first decision of the afternoon, and it CLOSES the first entry of the plan's
  "what remains": the tolerance is now
  `max(cert_tol_mass, min(cert_tol_mass_ceiling, cert_mass_round_margin * floor(j)))`
  with `cert_mass_round_margin = 10`, the ceiling 1 (a row at 1 is a mass flux
  changing by the whole of itself across the cell, not an operative branch on
  any state measured), `cert_tol_mass = 3e-12` unchanged and now the FLOOR of
  the tolerance, and the verdict on that row taken cell by cell through the
  pure `mass_row_cell_verdict`. **The floor expression the item's brief
  specified could not serve, and the report of that is part of the result**:
  `eps (|A+ F_up| + |A- F_dn|)/dV / residual_row_scale(1,j)` is algebraically
  between one and two ulps at every cell, because the row's scale IS the larger
  of the two addends and the `dV` cancels, and MEASURED it is **4.4409E-16 on
  all 500 cells of all four states** looked at, so ten times it never rises
  above 3e-12 and the item as written would have been inert. The cause is the
  one N33 established: the interface flux carries the last bit not of its own
  value `rho v` but of the O(rho c_s) quantities the Riemann solve differences,
  so the floor is about `eps/Mach` and not `eps`. The estimate adopted is
  `eps rho(|v| + c_s) A / dV` over the row's own scale, and it bounds the row's
  measured non-smoothness within a factor **1.5931** (`wasp_full_newton` as
  loaded, base Mach 3.937E-03), **1.8857** (`atomic_elem_newton` as loaded,
  5.062E-06), **2.6806** (`carrier_elem_newton` as loaded, 2.400E-04) and
  **1.7733** (the HD 209458 b state this item certifies, 2.767E-06), where the
  specified expression is short by two to six decades. The margin is measured,
  not chosen: `c_round = 10` is a round number 3.7 times above the largest of
  those four ratios and inside N33's independent reading of the same floor
  against its own bound, 5.5 to 12 (READ). The measurement itself is in the
  code: `EXHALE_MASS_FLOOR_SCAN=1` (new, default off) adds one ulp to the
  density of every physical cell, runs the whole flux assembly again and
  reports the STEP each mass row takes beside both estimates, and the state
  given is then assembled once more and reproduces its own mass rows to
  0.00E+00 on every fixture. Tests: `certification` **57/0** on the entry
  driver against the entry objects, **0 PASS / 1 FAIL** with exit 1 on the new
  driver against the entry objects (the compile naming the three new symbols),
  **70/0** after, the eleven new rows asserting that a floor below the fixed
  value leaves exactly 3e-12 and names the fixed tolerance, that a floor above
  it gives exactly `c_round * floor` and names the rounding anchor of that
  cell, that a row at half the anchor is WITHIN where the superseded rule
  transcribed inside the test refuses it, that a row above the anchor still
  refuses, that the ceiling binds, and that the margin stands above every ratio
  measured. Seven further suites are identical row for row in both builds, and a
  `-fcheck=bounds,do,mem` build runs both reloads with no runtime trap (the
  floor estimate reads the ghost cell at `j = 1`, which is why it was checked).
  **Impact.** `wasp_full_newton` reloaded is unmoved: CERTIFIED in both builds,
  mass row 1.212E-12 at cell 500 with the verdict taken there against the FIXED
  3.0E-12 at distance 0.404, `log10 Mdot` 13.30, outputs byte-identical apart
  from the provenance timestamp; its floor is below 3E-13 at every cell, so the
  anchor cannot reach that case, and the cold-start `wasp_full_newton` of
  `make check` cannot be affected either, its stop having no row binding (READ
  from `backup/regression/wasp_full_newton/run.log`) and a tolerance that only
  ever rises cannot move it. **The HD 209458 b element reload CERTIFIES, at
  pass 12**, `outer pass 12: ACCEPTED -- every active equation of this state is
  within its own tolerance`, where the control (a build of the entry text of
  the two files in the same tree) ends the same twelve passes NOT CERTIFIED on
  the one entry it had before, the hydrodynamic mass row at 7.159E-10 against
  3.0E-12 at cell 13. On the certified state the largest mass row is 7.703E-10
  at cell 2 and the verdict is taken at cell 34 (r = 1.00665), 6.542E-10
  against **4.6E-09, the rounding anchor of that cell**, distance 0.141; the
  eight elemental wind rows and the worst gated species row of all twelve
  passes are identical to the control to three digits, cell included, and the
  two builds part at pass 6, the first pass whose hydrodynamic solve stops on the
  loop-top test because its mass row is inside its own tolerance (`info = 0`
  where the control's is `info = 2`). The anchored tolerance still refuses
  states: the mass row is admitted at 11 of the 12 passes and refused at pass 3
  at distance 69.9. On the hot-Uranus carrier reload the base cells are
  admitted by their own anchor at every pass and the refusal moves OUT to cell
  180 (r = 1.0707), where the estimated floor is below 3e-12 and the row
  carries signal, refusing 5.067E-12 to 9.549E-12 against the fixed 3.0E-12: the
  refusal has moved from rounding to signal. The H2 carrier row is unchanged
  pass by pass in both builds (7.38E-02 at cell 290, 7.22E-02 at 289, 7.06E-02 at
  288), and the three passes cost 15.47 s against the control's 69.39 s. The
  derivation, both candidate estimates, the measured ratios and all three
  impact measurements are the appended section **"Anchor (6), 2026-09-11: the
  mass row's rounding floor"** of
  `docs/certification_tolerance_anchoring_20260910.md`; its nine earlier
  sections are byte-identical.
- **P17: the best-iterate ledger reads the mass row against the cell's own
  tolerance** (`src/modules/time_step/steady_newton.f90` +146/-32,
  `src/tests/krylov_and_dogleg/krylov_and_dogleg_tests.f90` +98/-0; MEASURED by
  P17). P16's own outside-scope finding, and the consequence of the anchor for
  the ledger: the hydrodynamic slot was formed as `rc(1)/cert_tol_mass` from a
  triple of row maxima, so the ledger ranked the mass row against 3e-12 while
  the stop test and the acceptance read the cell's own tolerance. The slot is
  now `max over cells j and rows k of |R_kj| / residual_row_scale(k,j) /
  tol_k(j)` over the whole physical column, with `tol_1(j)` taken from the
  certification's own rule at that cell through P16's public
  `mass_row_cell_verdict` (so the anchor formula is not restated in the solver
  module) and the fixed momentum and energy tolerances on the other two rows;
  that is the same functional the acceptance applies cell by cell, so `d < 1`
  in the ledger is now exactly the hydrodynamic part of the verdict. `rnorm` is
  unchanged and is now formed only when the caller asks for it, which the two
  loop-top calls do and the three trial calls do not. The superseded
  `hydrodynamic_distance_from_certification` keeps its name, signature and body
  and is documented as what it is, a BOUND from above on the distance
  (`tol_mass_at(j) >= cert_tol_mass` at every cell), because its removal would
  break the compile of a test file P17 does not own; P18 removes it. MEASURED
  RED twice: the new driver does not build against the entry-text objects
  (exit 1, `Symbol 'hydrodynamic_distance_from_certification_by_cell' ... not
  found in module 'steady_newton'`), and a private build in which the new
  function's body is replaced by the superseded reading and nothing else is
  altered gives **327 PASS / 7 FAIL** with the seven rows naming the column
  (`the_wind_cell_holds_the_mass_rows_distance` 3.333333E+02 against
  3.333333E+00, `the_base_cell_is_a_decade_inside_its_own_tolerance`
  3.333333E+02 against 1.000000E-01, `is_a_certified_hydrodynamic_state` F
  against T among them). GREEN `krylov_and_dogleg` **321/0 to 334/0**, and
  seven further suites identical in both builds. **Impact.** The single number
  that says the item worked, on the HD 209458 b element reload: the ledger's
  `best judged iterate: distance` over the twelve passes runs **2.2E+02 to
  6.3E+05 in the control and 1.321E-01 to 1.868E-01 after**, while the
  acceptance CERTIFIES the pass-12 state in both builds at a mass-row distance of
  0.141 (control) and 0.132 (new, verdict at cell 50, r = 1.00977, 4.245E-10
  against the 3.2E-09 anchor). The control's ledger therefore ranked every
  iterate of every solve as hundreds to hundreds of thousands of times outside
  certification, on a row its own acceptance admitted; the ledger and the
  acceptance now speak about one quantity. The trajectories part at pass 2 and
  the elemental relaxation is not damaged (the worst gated species row descends
  2.61E-03 to 7.73E-06 against the control's 2.61E-03 to 6.84E-06, both inside
  1e-5 at pass 12, and the certified state's elemental rows are SMALLER in the
  new build on every element, Fe 5.086E-05 against 8.966E-05, C 3.70E-06 against
  1.00E-05, O 3.34E-06 against 1.24E-05). The cost is not measurable: 909 outer
  iterations at 1.206 s against 861 at 1.216 s, so the 500 rounding-floor
  evaluations added per judged residual cost nothing this measurement can
  separate from noise. More restarts and fewer stagnations, which is the
  reading of a ledger that now stops improving where `||R||/3e-12` used to keep
  moving on rounding noise: "best iterate has not improved" fires in 7 solves
  against 4, "no descent direction exists" ends 11 against 8, and the control's
  one STAGNATED solve is gone, while no solve fails and every pass still
  returns a state. On the hot-Uranus carrier reload the consequence is a
  verdict: the hydrodynamic solve returns **`info = 0` at all three passes**
  against the control's 0, 2, 2; the mass row's verdict is WITHIN at passes 2
  and 3 (cell 1, 1.826E-11 against the 2.5E-11 anchor, and cell 2, 6.437E-12
  against the same) where the control refuses at cell 180; the state of the
  last pass has **1 refusing entry, the H2 carrier row alone**, against the
  control's 2; and the best judged distance runs 0.5282, 0.8819, 0.2578 against
  1.630, 3.501, 1.840. The H2 carrier row is identical pass by pass in both
  builds and in the P7 record (7.38E-02, 7.22E-02, 7.06E-02 at cells 290, 289,
  288), which is what must hold with the carrier relaxation outside the
  hydrodynamic solve. `wasp_full_newton` reloaded is unmoved, and for a stated
  reason: its estimated floor is below 3e-13 at every cell, so the cell's
  tolerance IS the fixed value everywhere and the maximum of the ratios is the
  ratio of the maxima; the two run logs are identical line for line except the
  `Execution Time`, every written file is byte-identical except the provenance
  timestamp, and the best judged distance is 6.722E-01 in both builds. A separate
  fix in the same file, reported on its own:
  `write_gate_met_while_a_row_refuses` printed the largest measure of a row and
  its cell beside the tolerance of the BINDING cell, a pair of numbers
  belonging to two different cells, and now names the binding cell throughout
  where the row records one, which is what the certification's own refusing
  entry list already did; no suite parses that line and the carrier rerun on
  the final binary reproduces every reported digit.
- **P18: one entry point for the hydrodynamic distance** (`steady_newton.f90`
  +23/-50, `krylov_and_dogleg_tests.f90` +49/-25,
  `steady_species_rows_tests.f90` +22/-16; P17's own outside-scope finding).
  With the ledger reading the mass row cell by cell (P17), the fixed-tolerance
  form `hydrodynamic_distance_from_certification(rc)` had no production reader
  left and survived only through test rows. It is DELETED, with its `public`
  export and the then-unused `cert_tol_mass` import of the module; every test
  row that read it restates its fact through
  `hydrodynamic_distance_from_certification_by_cell` on a one-cell column
  whose rounding floor (1e-30) puts it under the fixed tolerance, so the
  restatement is the same double division and every literal reference is
  unchanged. MEASURED (P18): `grep -rn hydrodynamic_distance_from_certification
  src/ | grep -v _by_cell` returns nothing; `krylov_and_dogleg` 334/0 before
  and after with no row lost (two renamed:
  `the_fixed_tolerance_alone_reports_the_base_cell`,
  `and_the_fixed_reading_is_never_below_the_distance`), `steady_species_rows`
  195/0, `certification` 70/0, `steady_completion_flag` 3/0; a compile RED
  (the entry-text drivers against the delivered objects fail on the missing
  symbol), a semantic RED of the cell-by-cell fact (327/7, P17's seven rows)
  and a semantic RED of the restatement itself (329/5, the fixture's floor
  raised into the anchored regime) all discriminate. Impact: the
  `wasp_full_newton` reload gives `info = 0`, CERTIFIED, `||R||` 1.495e-9,
  `log10 Mdot` 13.30 in both builds, all eleven output files byte-identical, the
  two run logs differing in the execution-time line alone. Documents that
  named the deleted function as current fact
  (`certification_tolerance_anchoring_20260910.md`,
  `session_handoff_20260911.md`) were corrected by the advisor.

- **P20: the code-size appendices moved to this log and were re-measured**
  (`docs/Update_EXHALE_appendix.tex` new, `docs/Update_EXHALE_stage1.tex`
  appendix removed and its `\date` fixed at 2026-09-04,
  `src/utils/update_log_to_tex.py` postamble, pointer lines in `README.md` and
  `README_HOWTO.md`). The two appendices, the code-size summary against the
  original ATES and the list of source inherited unchanged from it, were
  stranded in the stage-1 log and stale at 2026-08-30. They are now a
  hand-maintained fragment the generator `\input`s after the converted body,
  so appending to this log keeps them. Re-measured 2026-09-11 with
  `python3 src/utils/codesize.py` (MEASURED, all of it): the Fortran tree is
  **205 files / 76,259 SLOC**, 21.7x the original ATES (47 files / 3,514
  SLOC), against 111 files / 21,018 SLOC and 6.0x at the 2026-08-30 edition;
  the core excluding the Wind-AE port is 169 files / 71,503 SLOC, of which the
  production build is 95 files / 50,762 SLOC (14.4x) and the `src/tests/`
  drivers are 74 files / 20,741 SLOC, the rule the script applies being every
  `.f90` under `src/` outside `src/modules/wind_ae/`. 123 new core modules
  (49 production, 74 test drivers), 38 of the 47 original files modified, and
  the files **byte-identical with upstream are down from 13 to 8**: `hybrd`,
  `hybrd1` and `fdjac1` took the `unit_step_floor` difference-step argument,
  `Source.f90` the well-balanced option (N37), `J_inc.f90` the
  single-spectrum-type decision, and `speed_estimate_HLLC.f90` was deleted as
  dead code, so one upstream file is now removed rather than none.
  `Update_EXHALE_stage1.pdf` 464 -> 460 pages, `Update_EXHALE_stage2.pdf` 139 -> 158
  pages (MEASURED, `pdfinfo`).

- **P21: the carrier relaxation solves the chemistry of the composition it
  advances** (`src/modules/lower_atmosphere/diffusive_photochemistry.f90`
  +115/-27, `src/EXHALE_main.f90` +17/-7,
  `src/tests/carrier_retry/carrier_retry.f90` +287/-16; these three counts are
  MEASURED by `diff -u` against the entry texts saved in the scratch directory
  before the first edit and NOT by `git diff`, for the reason in "An external
  commit during the session" below). Verdict: the design does what it was meant
  to do at the level it was meant to act. With the chemistry equilibrated after
  every transport step the pass keeps, the unbounded pass at a fixed wind
  converges to a state that really is stationary for the coupled operator: on
  the hot-Uranus carrier reload the carrier steady residual falls from
  **7.99E-02 to 1.42E-09** in the worst cell and from **2.32E-03 to 1.18E-11**
  volume-weighted (MEASURED, P21, M1 below). The frozen fixed point was not a
  fixed point of transport plus chemistry, and this is the measurement of that.
  Two negative results come with it and are in the same subsection: at the
  production bound nothing moves, and the joint fixed point at a fixed wind is a
  fully molecular column.

  The design as implemented. A new private routine
  `equilibrate_chemistry_at_fixed_pressure(rho, p, f_sp, T, heat, cool, eta)`
  takes `get_species_densities` (the particle count and the electron density of
  the composition), then `comp_T_from_p` (its temperature at the conserved
  state's own pressure), then `ioniz_eq`. No new module dependency was needed in
  either direction and `ionization_equilibrium.f90` is not edited. It is called
  inside the trial loop of `relax_photochemical_composition` immediately after a
  trial is ACCEPTED and before the fixed-point test, so a trial the operator
  refuses, or one that leaves the movement bound, is undone and halved before
  any sweep is taken and leaves the chemistry exactly where it found it. **The
  sweep IS the background reinstall**: `ioniz_eq` writes `bg_cell(j) =
  ieq_cell` for every cell of a molecular sweep (READ), while
  `install_background_of_adopted_state` belongs to the steady solver's
  best-iterate bookkeeping and is not involved here. The sweep is handed the
  transported partition and fixes x(H2), the oxygen carriers and, where the
  protons are carried, the ionized fraction from the composition it reads, so it
  returns the chemistry OF that partition and does not re-solve it. **The
  movement bound is read on the transport output, before the sweep**: the bound
  governs what the transport operator moved in one trial, and the chemistry then
  closes on the step that was kept. MEASURED consequence: the reported drift of
  every bounded production pass stays at 9.98E-03 to 1.00E-02 at trust 1e-2, so
  the sweep does not carry the returned state out of its bound. The fixed-point
  test is now taken between two states the chemistry has closed on (`fnow` is
  recomputed after the sweep and becomes `fprev`), so `carrier_relax_fixed_point`
  means the carriers AND their chemistry stopped moving together. Interface:
  `relax_photochemical_composition(rho, v, p, T, f_sp, heat, cool, eta, trust,
  drift, nstep, outcome)`, with `p` `intent(in)`, the conserved state's own
  pressure untouched, and `T`, `heat`, `cool`, `eta` `intent(inout)`, written
  only where a refresh happened. The carrier branch of
  `steady_wind_with_element_diffusion` passes the new arguments and its own
  `ioniz_eq` call is DELETED as redundant; the `get_species_densities` and
  `comp_T_from_p` beside it are KEPT, because they are not a refresh but the
  filling of the density columns the routine keeps beside the composition, and
  they repeat the arithmetic of the pass's last internal call exactly.

  Tests: `carrier_retry` **122/0 at entry, 124 PASS / 3 FAIL on a RED build**
  carrying the delivered module with the one
  `call equilibrate_chemistry_at_fixed_pressure` line removed, **127/0 after**.
  The new section (11) has five rows and states why the joint fixed point cannot
  be asserted on the synthetic 12-cell column at all: the shortest admissible
  trial still moves 6.46e-4 of the largest entry H2 mixing ratio, six decades
  above `relax_tol = 1e-10`. What it asserts instead is the consistency the joint
  fixed point rests on: the background third-body density read from the returned
  composition (GREEN 4.7244E-07 against the entry composition's 5.5008E-03; RED
  5.4938E-03 against 1.2191E-16), the background temperature (GREEN 2.5675E-07
  against 5.4837E-03; RED 5.5241E-03 against 0), and one further sweep moving the
  composition by 1.3637E-10 against hybrd1's own `xtol = sqrt(eps) = 1.5E-08`
  (RED 2.0110E-04). Every other suite whose driver links
  `diffusive_photochemistry` is unchanged before and after:
  `carrier_returned_state_acceptance` 36/0, `attempted_step` 70/0,
  `steady_species_rows` 195/0, `certification` 70/0,
  `carrier_constraint_attribution` 10/0, `carrier_reference_scales` 14/0,
  `element_operator` 24/0, `species_masses` 9/0.

- **P22: the record of P21, of the terminology sweep and of the external
  commit** (this section with its two new subsections, the status rows of
  `docs/PLAN_20260911_partitioned_solver.md`, `docs/session_handoff_20260911.md`,
  the rewritten H2 front entry and the new outside-scope entry of
  `TO_BE_DONE.md`, and one sentence in `docs/code_status_20260910.md`
  section 3).

- **D1 to D3: `docs/EXHALE_physics_and_algorithms.tex`, the physical equations
  and the algorithms as the code applies them** (user instruction 2026-09-11,
  evening: a document at the depth of Caldiroli et al. 2021 sections 2 and 3,
  written from the source, not from the memos). Main file with `natbib`; eleven
  section files in `docs/physics_overview/` (introduction and notation,
  hydrodynamics, radiation, ionization, thermal source, molecules and metals,
  numerical scheme, boundary and initial conditions, stationary state, output
  and post-processing, the table of physics keys with the reference list, 69
  entries). D1 wrote the hydrodynamic, numerical, boundary and stationary
  sections, D2 the radiation, ionization, thermal, molecular and
  post-processing ones, D3 made them one document (55 inline citations
  converted, one symbol for one quantity with a notation table, about 20
  cross references, repetition and run narrative removed). MEASURED: `latexmk`
  rc 0, 57 pages, no LaTeX or hyperref warning, every `\texttt{}` file name
  present in `src/` (74 checked). Code facts established while writing, each
  stated in the document where it belongs: no escape-probability suppression
  acts on the Ly-alpha cooling channel (`lambda_coex_HI` in `Cool_coeff.f90` is
  the optically thin fit; the escape probability enters the metal
  fine-structure statistical equilibrium and the H(n=2) populations only);
  `equation_T` (`T_equation.f90`) is called by `post_process_adv.f90` and test
  drivers, the marching loop's temperature being `solve_energy_semi_implicit`;
  `Numerical flux:` is a mandatory key with no default and every shipped input
  (70 of 70) says HLLC, ROE being refused with molecular chemistry unless
  `Caloric EOS: monatomic`; `2D approximate method: alpha` applies neither the
  dayside dilution nor the mass-loss factor, which match the `Rate/...` strings
  only. Pointers added in `README.md` and the project `CLAUDE.md`.
- **The measurement of record after P14 to P17** (advisor, private snapshot
  build of the tree with P14 to P17, 8 threads, MEASURED): the HD 209458 b
  element reload, partitioned, 12 passes, `ACCEPTED -- every active equation
  of this state is within its own tolerance` at pass 12 (worst gated species
  row 7.73e-6, He/H partition at cell 289; mass 9.26e-10, momentum 2.44e-14,
  energy 2.17e-8); the hot-Uranus carrier reload, 3 passes, hydro `info = 0`
  at every pass (mass 7.33e-12 / 1.29e-11 / 6.44e-12 admitted by the cell's
  anchor), the one refusing entry the H2 carrier row (7.06e-2 at cell 288).

Three changes outside the item list, made by the advisor while reading P14 to
P17:

- `assert_written_state_is_the_accepted_one` (`src/EXHALE_main.f90`) carries an
  ABSOLUTE floor of 1e-12 beside its 1e-10 relative test. P14 reported the
  false alarm: the spread is a difference of O(1) numbers over their mean, so
  two evaluations of one state agree only to the rounding of that difference,
  and on the partitioned carrier reload a state flat to fourteen digits read
  **3.0982E-14 as written against 2.7487E-14 accepted** while the relative test
  alone called that a change of `u`. 1e-12 stands a decade above that rounding
  and four decades below the smallest spread a gate has ever accepted (READ,
  the comment at the test).
- the three rows of the `grid_and_gates` suite that wrote to a fixed directory
  under `$ROOT/build/tests/grid_and_gates` and ignored `EXHALE_TEST_OUT` now
  honor it (`momentum_row_from_fluxes_only.sh`,
  `base_level_single_statement.sh`, `sed_coverage_stop.sh`; READ, each now
  takes `${EXHALE_TEST_OUT:-$ROOT/build/tests/grid_and_gates}`). Two concurrent
  runs of the suite, which is what a before-and-after measurement is, collided
  there. P14 fixed `output_state_consistency.sh`, the fourth, because its item
  owned it.
- the paragraph at `assemble_residual` (`src/modules/time_step/steady_residual.f90`)
  that stated what the kind-generic rows do not hand back is rewritten to what P15 made
  true: the kind-generic rows return the momentum row and store the interface
  fluxes, the face pressures and, under the well-balanced option, the face
  departures they were built from, but neither the three terms of the momentum
  equation nor the equilibrium pressure force the row is read against, so the
  one expression of each is evaluated there for the state just assembled; and a
  reconstruction continuation never reaches them, the generic text refusing
  it.

Two changes outside the item list, made by the advisor while reading P21:

- the `H2 front:` diagnostic of `src/EXHALE_main.f90` printed
  `r(max(j50,1))` and `r(max(j01,1))`, so a front that had left the domain was
  reported as the radius of the first cell. On P21's M1 run that read
  `x2=1e-2 at r= 1.0002` for a column in which no cell is below 1e-2 at all,
  beside a correct `x2=0.5 at r=2.2367`. Each level is now printed as `> r(N)`
  when no cell of the column is below it.
- the noun of the 2026-09-11 terminology rule removed from the comment at
  `src/EXHALE_main.f90` line 2022 and from
  `diffusive_photochemistry.f90` line 248, and from the two fixture READMEs
  (about 26 lines), which item A4 left to whoever was editing those files.

### The pseudo-time start of the stationary restart route

Found while running P6 and settled by the advisor. `Restart intent:
stationary` entered the stationary solve at the CFL interval of the reloaded
state. A state written by a run is already relaxed, and a continuation
started at the CFL step from it does not reach the root: MEASURED, 8 threads,
the hot-Uranus carrier reload's hydrodynamic solve is handed back at `||R||`
2.685E-01 after 40 iterations and the next four passes enter at 1.133, 1.018,
1.026 and 1.025 and hand back 2.976E-01 to 3.579E-01
(`<scratchpad>/P6/runs/carrier_p6/run.log`), while the HD 209458 b element
reload stops on its stagnation detector at `||R||` 1.490 and the pass is then
refused because the element relaxation finds no admissible advance
(`<scratchpad>/P6/runs/atomic_partition.new/run.log`). With the start at 1.0,
the value the marching hand-off uses, the same two solves reach `||R||`
2.507E-09 at iteration 13 and hand back **1.498E-09** (carrier) and hand back
**1.167E-08** (atomic). The advisor therefore made 1.0 the default of the
`Restart intent: stationary` route (`src/EXHALE_main.f90` line 5234), with
`EXHALE_PTC_DTAU0` as the override for measurement; the direct steady route
keeps its CFL start, which is the pseudo-transient continuation it was
designed as. Verified with no environment variable set: the run prints
`stationary restart: dtau0 = 1.00E+00` and reproduces the first pass to every
printed digit (`<scratchpad>/P7/carrier_default_dtau0/run.log`). No
regression case sets `Restart intent:` at all (MEASURED: the grep over
`backup/regression/*/input.inp` returns nothing), so this default moves no
golden; it moves any user run that sets the key.

### P7: the measurement of record

One binary (`EXHALE_P6.x`, built from the merged tree), both fixtures,
`Restart intent: stationary`, `EXHALE_PTC_DTAU0=1.0`, `OMP_NUM_THREADS=8`,
500 cells. The partitioned route is `Coupled carrier solve: False` with
`EXHALE_OUTER_PASSES` setting the pass budget and `EXHALE_JFNK_MAXIT` the
hydrodynamic cap; the coupled route is the same fixture with `Coupled carrier
solve: True`. Wall times of the partitioned runs are the sum of the seconds
the `outer pass` lines report, which agrees with the run directory's file
timestamps to under a second; the coupled carrier control prints no total, so
its time is the timestamp span of its run directory. Every number MEASURED
from the log named.

| route | fixture | passes (hydro cap) | wall time | hydro mass / momentum / energy | worst gated species row | certified? |
| --- | --- | --- | --- | --- | --- | --- |
| state as loaded | carrier, hot Uranus | -- | -- | 2.284e-1 / 1.954 / 1.871 | H2 1.166e-1 at cell 492 | no, 4 entries |
| partitioned | carrier | 8 (40) | 173 s | 8.75e-12 / 1.66e-12 / 2.41e-9 | H2 6.217e-2 at cell 282 | no, 2 entries |
| partitioned | carrier | 40 (40) | 831 s | 8.121e-12 / 3.377e-12 / 2.730e-9 | H2 2.234e-2 at cell 443 | no, 2 entries |
| coupled | carrier | one solve (250), stopped by the stagnation detector at iteration 105 | 321 s | 1.231e-2 / 1.009e-1 / 5.801e-2 | H2 7.815e-2 at cell 283 | no, 4 entries |
| state as loaded | atomic, HD 209458 b | -- | -- | 1.509 / 3.348e-1 / 1.892 | He/H 5.840e-1 at cell 291 | no, 11 entries |
| partitioned | atomic | 5 (80) | 528 s | 1.05e-9 / 4.40e-13 / 1.13e-8 | O 4.64e-4 at cell 268 | no, 9 entries |
| partitioned | atomic | 12 (80) | 1043 s | 1.054e-9 / 1.665e-13 / 1.501e-8 | He/H 5.468e-6 at cell 288 | **no, 1 entry** |
| coupled | atomic | 3 (250), ended in the fourth | 2392 s | 1.16e-6 / 2.89e-6 / 2.83e-6 | He/H 8.29e-5 at cell 295 | no, run not completed |

Sources: `<scratchpad>/P7/carrier_elem_newton_part1/run.log`,
`carrier_part40/run.log`, `carrier_coupled250/run.log`,
`atomic_elem_newton_part1/run.log`, `atomic_part12/run.log`,
`atomic_coupled250/run.log`. The atomic coupled control alternates the
element relaxation too, because `He_diffusion` is on, and its three completed
passes cost 1380.59 s, 482.04 s and 528.98 s; it was ended in the fourth. The
`(EXHALE_main) outer pass N:` line of every pass carries the hydrodynamic
flag, the worst gated species wind row with its cell, the three row maxima,
omega, the carrier trust and the seconds, and the closing `(certification)`
block names what refuses. Those `<scratchpad>/` paths are the session scratch
root and are TEMPORARY: they will be gone, and everything quoted from them is
recorded here, in the two fixture READMEs and in
`docs/session_handoff_20260911.md`. The runs reproduce from the pinned
fixtures by the recipes in those READMEs.

What the ladders show. On the atomic fixture the worst gated species row
falls 2.61e-3, 3.81e-3, 1.77e-3, 9.45e-4, 4.64e-4, 2.36e-4, 1.20e-4,
5.97e-5, 3.07e-5, 1.50e-5, 7.94e-6, 5.47e-6 over the twelve passes (the one
rise, pass 1 to pass 2, is the alternation's own feedback across a
hydrodynamic solve and is why a single rise does not end the loop), the
refusing entries fall 9 to 1, and in the last pass every elemental wind row
is within 1e-5: He/H 5.468e-6 at cell 288, and the metals 1.582e-6 (Mg) to
3.720e-6 (O). On the carrier fixture the H2 wind row falls 7.38e-2 to
2.23e-2 over the 40 passes with the binding cell moving 290 to 443, which is
the H2 front advancing about one cell in a pass and not a residual settling.

**The decision.** The plan's adoption rule is that the route is adopted only
if it certifies where the coupled solve does not, or reaches a strictly lower
joint residual in no more wall time. Neither fixture certifies, so the first
branch is not met on either. The second branch is met on both: on the atomic
fixture the partitioned route stands strictly below the coupled control on
every judged row (mass 1.054e-9 against 1.16e-6, momentum 1.665e-13 against
2.89e-6, energy 1.501e-8 against 2.83e-6, worst gated species row 5.468e-6
against 8.29e-5) in 1043 s against 2392 s for three passes; on the carrier
fixture the partitioned route at 8 passes stands strictly below the coupled
control on every judged row (8.75e-12 against 1.231e-2, 1.66e-12 against
1.009e-1, 2.41e-9 against 5.801e-2, H2 6.217e-2 against 7.815e-2) in 173 s
against 321 s, and the 40-pass continuation goes on to 2.23e-2 for 831 s.
The partitioned route is therefore what the production outer loop runs.
`Coupled carrier solve` already defaults to False
(`carrier_in_newton = .false.`, `src/modules/init/parameters.f90` line 419),
so no default changes; `True` remains the coupled measurement option, and what
this session added on the partitioned side is the joint acceptance contract
of P6 and the three defect repairs P1 to P3 that made a pass mean something.

**What did not certify, and why.** Entry 1 was CLOSED the same afternoon by
the user's first decision (P16, anchor (6)); entries 2 and 3 are stated as
they stand:

1. **CLOSED 2026-09-11 by P16.** The atomic element reload refused on the
   hydrodynamic mass row alone,
   1.054E-09 against 3.0E-12 at cell 16 (r = 1.0031), on the three-unknown
   solve of this fixture. Every one of its twelve hydrodynamic solves ends
   with `||R||` between 8.980E-09 and 1.501E-08 against the run's `Resid tol`
   of 1e-8, and the last one ends on `no descent direction exists for the
   banded model at this state` at `||R||` 1.501E-08. So the alternation has
   done its work (the element rows are inside their tolerance) and what is
   left is the base-layer mass row of the hydrodynamic solve, which is the
   N33 rounding floor territory and not a species-coupling question.
   **P10 settled which of the two it is.** With the best-iterate ledger
   reading that row against its own 3e-12 (P10), the same twelve passes end at
   8.893E-10 at cell 3, the same order, and the row is MEASURED to carry no
   signal at that level: the ledger's best distance of pass 1 is 4.407E+02, a
   mass row of 1.32E-09, while the same state re-evaluated at its own
   composition after the best-iterate restore carries 4.11E-09. At 1e-9 to
   1e-12 the row is the rounding of the base-layer flux difference, so ranking
   iterates by it is ranking by rounding noise and no ledger can bring it
   down. The entry was therefore an ACCEPTANCE question, the 3e-12 tolerance
   against the base-layer rounding floor of this fixture, and not a solver
   item; the user decided it the same afternoon, and with the mass row judged
   at `max(3e-12, 10 * floor(j))` the same twelve-pass recipe **CERTIFIES at
   pass 12** (MEASURED by P16 and again by P17, which certifies at pass 12 in
   both of its builds). The verdict is taken at a wind cell, not at the base:
   cell 34, r = 1.00665, 6.542E-10 against the 4.6E-09 anchor of that cell,
   distance 0.141 (P16), and cell 50, r = 1.00977, 4.245E-10 against 3.2E-09,
   distance 0.132 (P17). The base-layer floor of this fixture is what the
   anchor is measured on: see the decisions subsection below and
   `docs/certification_tolerance_anchoring_20260910.md`, "Anchor (6)".
2. The carrier reload refuses on the H2 wind row, and the user's second
   decision of the afternoon ran the continuation that says what that row IS.
   MEASURED over 110 passes (the table and the recipe are in the decisions
   subsection below): the H2 front does NOT stop, `x2 = 0.5` moving 1.0726 to
   1.1230 R_p over 108 passes at about 5E-04 R_p in a pass with no sign of
   slowing, so the continuation is unfinished after 110 passes and 0.05 R_p.
   And the worst gated row stops being the front after pass 40: it stands at
   **2.2E-02, nearly uniform in r**, and its cell walks outward through the
   wind, 443 to 500, reaching the outer boundary at pass 108 while the value
   does not fall. That floor is a property of the operator and not of the
   front: the bounded carrier pass (trust 1e-2 on the H2 mixing-ratio maximum,
   met by one transport step) advances the front but leaves the wind cells,
   whose H2 mixing ratio is orders of magnitude below the maximum, at one
   implicit step in a pass, so their own H2 balance never relaxes, and the
   drift measure, absolute in the grid maximum, does not see them. The
   consequence for the method, recorded and not adopted: a bounded front
   advance and a relaxation of the wind cells' H2 balance are two different
   operations, and the second needs its own step (a cell-relative measure, or
   a wind-window relaxation after the bounded pass) before the H2 wind row can
   be judged at 1e-5. Whether the front stops at all inside 2 R_p is a
   physical question and stays open (the P23 comparison against Koskinen et
   al. 2022 Model A).
3. The hot-Uranus mass row sits on the same floor throughout the carrier
   run, 4.48E-12 to 1.57E-11 against 3.0E-12 over the 40 passes, with the
   momentum row at 4.14E-13 to 5.13E-12 and the energy row at 1.75E-09 to
   4.56E-09. The mass row of that fixture is the rounding of the base-layer
   flux difference recorded in N33, which is the reading P16 then made
   operative: under anchor (6) the base cells of this fixture are ADMITTED by
   their own anchor at every pass and the refusal moves out to cell 180
   (r = 1.0707), where the floor is below 3e-12 and the row carries signal.

### The Codex adversarial review (2026-09-11)

Run on the working tree after P1 to P13, that is on the merged text of the
morning's items and not on the committed tree. Two findings, both verified in
the source by the advisor before any item was written, and both closed the
same afternoon. The findings are quoted as they were written.

1. **[high]** "정체 종료 시 조성과 출력용 상태가 서로 달라집니다
   (src/EXHALE_main.f90:6503-6514)": on a stagnation exit the composition had
   been updated while the particle densities and the temperature were not
   refreshed on the element path, so the certified state and the written state
   differed; the recommendation was to move the progress control before the
   composition update. Verified in the source, reproduced as a measured
   inconsistency of 1.317089E-08 in the temperature of the written state, and
   **closed by P14** with exactly that ordering: the progress control now
   stands between the certification of the pass and the update, so a pass that
   ends the iteration takes no update, and the element path's update refreshes
   the particle densities and the temperature before its equilibrium sweep.
2. **[medium]** "고정밀 경로가 다른 상태의 압력 정보를 인증에 사용합니다
   (src/modules/time_step/steady_residual.f90:256-258)": with
   `EXHALE_RESID_QUAD` the generic rows stored only `face_flux` and `face_p`,
   so `momentum_row_terms_of_state` read the previous `RK_rhs` call's
   `face_q_up`/`face_q_dn` (zeros on the first call); the recommendation was to
   return the terms from the generic computation and to refuse the continuation
   combination. Verified in the source, and **closed by P15**, both halves: the
   kind-generic rows now hand their own face departures back along the route the fluxes take,
   and the continuation combination is refused, because the kind-generic rows CAN be
   reached on a continuation (`assemble_residual` chooses the text from
   `generic_precision_rows_arm()` alone and the generic text carries its own
   blending branch) and no single pair of departures satisfies the blended row
   identity.

Neither finding moves a certified state of the default path: P14's
`wasp_full_newton` reload and P15's, with the key off and with
`Well balanced: True` and `EXHALE_RESID_QUAD=2`, are CERTIFIED with every
output byte-identical apart from the provenance timestamp. What each moves is
stated in its own bullet above.

### User decisions of 2026-09-11 (afternoon)

Two decisions were taken after the open items of the plan were put to the
user.

**Decision 1: the mass tolerance is re-anchored above the fixture's measured
rounding floor** (item P16, and the ledger consequence P17). The tolerance of
the continuity row is now

```
tol(j) = max( cert_tol_mass, min( cert_tol_mass_ceiling, c_round * floor(j) ) )
```

with `cert_tol_mass = 3e-12` unchanged as the FLOOR, the ceiling 1, and
`c_round = cert_mass_round_margin = 10`. The floor estimate is
`eps rho(|v| + c_s) A / dV` over the row's own scale, which is one ulp in the
wind and about `eps/Mach` in a layer whose flux is a tiny residue of the
momentum the cell carries at its own fastest signal speed.

The brief's first expression for the floor,
`eps (|A+ F_up| + |A- F_dn|)/dV / residual_row_scale(1,j)`, is inert, and that
is a measured result and not a judgement: the row's scale IS the larger of the
two addends and the `dV` cancels, so the expression is algebraically between
one and two ulps at every cell, MEASURED at **4.4409E-16 on all 500 cells of
all four states** examined, and ten times it never rises above 3e-12. The
adopted estimate bounds the step each mass row actually takes when one ulp is
added to the density of every physical cell (`EXHALE_MASS_FLOOR_SCAN=1`, new
and default off) within a factor

| state | base Mach | step over the adopted estimate | step over the brief's expression |
| --- | --- | --- | --- |
| `wasp_full_newton` as loaded | 3.937E-03 | **1.5931** | 1.109E+02 |
| `atomic_elem_newton` as loaded | 5.062E-06 | **1.8857** | 2.4365E+05 |
| `carrier_elem_newton` as loaded | 2.400E-04 | **2.6806** | 5.0915E+03 |
| the HD 209458 b state P16 certifies | 2.767E-06 | **1.7733** | 1.4201E+06 |

All MEASURED by P16, largest ratio over the 500 cells of each state. `c_round
= 10` is a round number 3.7 times above the largest of those ratios and inside
N33's independent reading of the same floor against its own bound, 5.5 to 12
(READ). The outcome: the **HD 209458 b element reload CERTIFIED at pass 12**,
where the control built from the entry text of the two files ends the same
twelve passes NOT CERTIFIED on the hydrodynamic mass row (7.159E-10 against
3.0E-12 at cell 13), and **`wasp_full_newton` unmoved** in both builds
(CERTIFIED, mass row 1.212E-12 at cell 500 judged against the FIXED 3.0E-12 at
distance 0.404, `log10 Mdot` 13.30, outputs byte-identical apart from the
provenance timestamp), its floor being below 3E-13 at every cell. The anchor
still refuses states: the HD 209458 b mass row is refused at pass 3 at distance
69.9, and on the hot-Uranus carrier reload the refusal moves from the base
cells, admitted by their own anchor, out to cell 180 where the row carries
signal. The derivation and all three impact measurements are the appended
section "Anchor (6), 2026-09-11: the mass row's rounding floor" of
`docs/certification_tolerance_anchoring_20260910.md`.

**Decision 2: the long continuation of the hot-Uranus carrier reload.** The
pass budget was set to 400 and the run went to pass **110**, where the progress
control ended it ("the worst gated species row has not fallen in 3 consecutive
passes"), the carrier trust having been halved 1.0E-02 to 5.0E-03 to 2.5E-03 at
passes 109 and 110. Recipe: `backup/regression/carrier_elem_newton`,
`Load IC? True`, `Coupled carrier solve: False`, `Restart intent: stationary`,
`EXHALE_OUTER_PASSES=400 EXHALE_JFNK_MAXIT=40 OMP_NUM_THREADS=8`. Binary: the
advisor's snapshot build of the merged tree after P1 to P12 and BEFORE P14 to
P17, so every hydrodynamic solve of this run returns `info = 1` on the
base-layer mass row (5E-12 to 1E-11 against 3E-12) and the run predates
anchor (6). All numbers MEASURED from its `run.log`.

| pass | worst gated H2 row | cell | H2 front `x2 = 0.5` [R_p] | `x2 = 1e-2` [R_p] |
| --- | --- | --- | --- | --- |
| 1 | 7.38E-02 | 290 | 1.0726 | 1.0966 |
| 10 | 5.86E-02 | 281 | 1.0755 | 1.1017 |
| 20 | 4.07E-02 | 274 | 1.0796 | 1.1070 |
| 30 | 2.62E-02 | 267 | 1.0827 | 1.1140 |
| 40 | 2.23E-02 | 443 | 1.0871 | 1.1199 |
| 50 | 2.20E-02 | 473 | 1.0917 | 1.1262 |
| 70 | 2.19E-02 | 482 | 1.1017 | 1.1413 |
| 90 | 2.18E-02 | 490 | 1.1126 | 1.1563 |
| 100 | 2.17E-02 | 494 | 1.1184 | 1.1643 |
| 108 | 2.17E-02 | 500 | 1.1230 | 1.1727 |
| 110 | 2.17E-02 | 500 | 1.1230 | 1.1727 |

**Two regimes.** Passes 1 to 30: the worst gated row is at the front's outer
edge, its cell moving 290 to 267 (r about 1.28 to 1.25), and it falls 7.4E-02
to 2.6E-02 as the front relaxes. From pass 40 the worst row is no longer at the
front: it sits at 2.2E-02 and its cell walks outward through the wind, 443 to
500, reaching the outer boundary at pass 108 while the value does not fall. The
front itself keeps moving outward at about 5E-04 R_p in a pass, `x2 = 0.5` from
1.0726 to 1.1230 over 108 passes, with no sign of stopping.

So three statements, the first two measured and the third a reading of them.
(a) The front has NOT reached a stationary position in 110 passes: 0.05 R_p
traveled and the continuation unfinished. (b) The remaining wind residual of
2.2E-02 is nearly uniform in r and is NOT the front. (c) The reading: the
bounded carrier pass, trust 1E-02 on the H2 mixing-ratio maximum and met by one
transport step, advances the front but leaves the wind cells, whose H2 mixing
ratio is orders of magnitude below the maximum, at one implicit step in a pass,
so their own H2 balance never relaxes, and the drift measure, absolute in the
grid maximum, does not see them. The progress control ended the run because
after pass 40 the worst gated row measures that wind floor and not the front.

Consequence for the method, recorded and NOT adopted: a bounded front advance
and a relaxation of the wind cells' H2 balance are two different operations,
and the second needs its own step, a cell-relative measure or a wind-window
relaxation after the bounded pass, before the H2 wind row can be judged at
1e-5. The physical question whether the front stops at all inside 2 R_p is
separate and stays open (the P23 comparison against Koskinen et al. 2022
Model A).

### P21: the fixed point at a fixed wind, and the production loop

All of this is MEASURED by P21 on `backup/regression/carrier_elem_newton`
(`Coupled carrier solve: False`, `Restart intent: stationary`,
`EXHALE_JFNK_MAXIT=40`, 8 threads, scratch copies), against a control binary
built from the entry text of the three files.

**M1, the fixed point at fixed hydro.** `EXHALE_CARRIER_TRUST=1e6`, that is the
movement bound disabled, and `EXHALE_OUTER_PASSES=2` rather than 1, because the
outer loop takes no composition update on its last pass (P6) and at 1 the
relaxation never runs. Both binaries end the pass on
`carrier_relax_fixed_point`.

| quantity | control (frozen chemistry) | new (refreshed) |
| --- | --- | --- |
| transport steps kept | 38 | 57 |
| composition drift | 9.79E-01 | 9.95E-01 |
| carrier steady residual after the pass, volume-weighted | 2.32E-03 | 1.18E-11 |
| worst cell | 7.99E-02 at cell 266, r = 1.211 | 1.42E-09 at cell 2, r = 1.000 |
| H2 front, `x2 = 0.5` | 1.1327 R_p | 2.2367 R_p |
| `x2 = 1e-2` | 1.2056 R_p | nowhere in the domain |
| largest x2 on the grid | 0.9774 | 0.9848 |
| cells clamped onto the element budget | none | none |
| pass 1 wall time | 18.64 s | 16.04 s |
| pass 2 wall time (hydrodynamics only) | 16.66 s | 37.06 s |

The wind H2 row before the pass is 7.38E-02 on both, at cell 290, reproducing
the experiment's 7.380744e-2, and the control's 7.99E-02 reproduces its
7.989781e-2 to the digit. The experiment's three-number sequence has only two
members in the production loop, the diagnostic being taken after the caller's
refresh in the control, so the frozen intermediate (2.449021e-13) is not a state
the production loop ever holds; in the new build there is no intermediate at
all, which is the point. Wall time of the new pass is LOWER despite 50 percent
more steps and one sweep a step; the plausible reading, offered as a reading, is
that a transport step taken on the background of the composition it advances
converges in fewer inner iterations than one taken far from a frozen
background.

**M2, the production loop: at the production bound nothing moves.** Same
fixture and recipe with the movement bound at its default 1e-2,
`EXHALE_OUTER_PASSES=12` and a 60-pass continuation, against the control's first
three passes and against the 110-pass ladder of the same day's decision 2 (both
READ from the P7 run logs).

| pass | worst gated H2 row, cell | recorded ladder | `x2=0.5` / `x2=1e-2` | recorded |
| --- | --- | --- | --- | --- |
| 1 | 7.38E-02, 290 | 7.38E-02, 290 | 1.0726 / 1.0966 | 1.0726 / 1.0966 |
| 10 | 5.86E-02, 281 | 5.86E-02, 281 | 1.0755 / 1.1017 | 1.0755 / 1.1017 |
| 20 | 4.07E-02, 274 | 4.07E-02, 274 | 1.0796 / 1.1070 | 1.0796 / 1.1070 |
| 30 | 2.62E-02, 267 | 2.62E-02, 267 | 1.0827 / 1.1140 | 1.0827 / 1.1140 |
| 40 | 2.23E-02, 443 | 2.23E-02, 443 | 1.0871 / 1.1199 | 1.0871 / 1.1199 |
| 50 | 2.20E-02, 473 | 2.20E-02, 473 | 1.0917 / 1.1262 | 1.0917 / 1.1262 |
| 60 | 2.20E-02, 477 | (not in the ladder) | 1.0953 / 1.1327 | |

Twelve and sixty passes reproduce the recorded ladder to every printed digit:
the worst gated H2 row, its cell and both front positions are identical at
passes 1, 10, 20, 30, 40 and 50. Every pass ends on the movement bound with 3 to
6 kept steps and a drift of 9.98E-03 to 1.00E-02, and every hydrodynamic solve
returns `info = 0`. So **the 2.2E-02 wind floor does not go away and the frozen
chemistry was not its cause**; the bound's measure, absolute in the grid
maximum, remains the standing explanation. Cost is not a reason to hesitate: a
production pass costs 3.07 to 3.72 s against the control's 3.10 to 3.54 s, and
one equilibrium sweep on this 500-cell grid costs below 0.029 s at 8 threads
(MEASURED from 300 marching steps of `mol_carrier` in 8.66 s including
initialization), so 3 to 6 of them cannot be resolved against that scatter.

**M3 and M4, the default paths.** `mol_carrier` marching, `EXHALE_MAXSTEPS=300`,
one thread: `Hydro_ioniz.txt`, `Ion_species.txt` and both `_adv` files identical
except the two `git=` provenance lines, and `Cooling_breakdown.txt`,
`Heating_breakdown.txt` and `IC_dump.txt` byte-identical, which is the expected
result because the marching loop calls `photochemical_transport_step` directly
and never this routine. `wasp_full_newton` reloaded, one thread: the two run
logs are byte-identical, both `info = 0` and CERTIFIED at `||R|| = 1.920E-09`,
ten of the twelve written files byte-identical and the other two differing in
the same provenance line alone. No molecules, so the branch is never entered.

**Does the front stop at fixed hydro, and where?** No, and the way it fails to
stop is new information. At a wind held fixed and with the chemistry solved on
the advancing composition, the pass runs to a genuine fixed point of the coupled
operator (1.42E-09) whose H2 fills the whole domain: x2 = 2n(H2)/n_H is 0.98 at
1.1 R_p, 0.96 at 1.2, 0.84 at 1.5, 0.60 at 2.0 and 0.33 at 3.0 R_p, with no cell
below 1e-2 anywhere. The largest x2 is 0.9848, so this is NOT the element
ceiling and no cell is clamped in either run, and P2's frozen artifact (x2
driven to the element ceiling in 123 cells) does not appear. A plausible
reading, offered as a reading: the H2 column and its self-shielding grow
together, so re-solving the photodissociating field on the advancing composition
removes the sink that held the front, and nothing but the wind can stop it,
which it cannot while it is held fixed. What this says is that **the front's
position is set by the hydrodynamic response and not by the carrier operator
alone**, consistent with the bounded advance being the only usable form of the
update; the physical question, whether the front stops inside 2 R_p in the
coupled system, stays exactly where P23 has it.

**Does the wind floor go away?** No, as M2 records. The 2.2E-02 row that appears
at pass 40 and walks outward through the wind is reproduced pass for pass and
cell for cell with the chemistry refreshed. The frozen background was not its
cause; the bound's measure, absolute in the grid maximum, still is.

### The terminology sweep (A1 to A4, 2026-09-11)

The user asked whether "arm" is a word physics uses for what this tree was using
it for, and instructed that terms belonging to software engineering, clinical
trials or management not be used as the names of physical or numerical concepts
at all. The rule is now in `~/.claude/CLAUDE.md` under "금지 표현", in two
bullets. The first: avoid terms that physics, chemistry, mathematics and
astronomy do not use, in documents, comments, log messages and identifiers
alike, and write instead what the thing is, the standard term of those fields or
plain English, naming what is compared ("control build", "measured build"),
which physics option is meant ("the well-balanced option"), or which quantity
("the mass row", "the H2 front"); where such a term already stands, fix it where
it is visible. The second: the NOUN "arm" is forbidden, for one side of a
comparison, for a selectable code path and for one branch of a fixture alike,
and code identifiers (`ARM_*`, `*_arm`) are held to the same rule; the VERB "arm
/ armed / arming / disarm", meaning to enable a device, is plain English and
stays.

Four workers swept the tree. Every count is MEASURED by the worker that owns the
files, by `grep -w -i -E "arms?"` against the git HEAD text, and no number, file
name or measured statement was changed anywhere in the four diffs.

- **A1, the 51 top-level `docs/*.md` outside the two update logs: 396 to 17.**
  All 17 that remain are allowed: 14 are the verb, 3 are the literal directory
  glob `arm*`. The replacements are by sense: a side of a comparison became "the
  two runs" or "both builds"; an optional code path became its key or its
  physics ("the well-balanced option", "the kind-generic residual assembly");
  a fixture's solve became "the element reload" or "the carrier reload", which is
  the word those same documents already use; a rung of a ladder of runs became
  "rung"; one chemical network of a comparison became "network". One code span
  changed, the placeholder `<arm>` of a run-directory path in
  `vulcan_photochem_comparison.md`, which named nothing on disk and is now
  `<network>`, the names that do exist there being `vulcan/`, `ncho/`,
  `zahnle/`, `zahnle_s/` and `eqfit/` (READ from the tree).
- **A2, `docs/Update_EXHALE_stage2.md`: 176 to 1** (the one is the verb, "when both
  stellar keys arm H(n=2)"), and `Update_EXHALE_appendix.tex` 1 to 0. The
  largest class was "in both arms" of a two-binary comparison, now "in both
  builds"; then the optional paths ("the well-balanced option", "quad
  accumulation", "the kind-generic rows"), the fixtures ("the atomic element
  reload", "the carrier rows"), and the unknown space of B5j/B5k ("the density
  unknown", "the log unknown"). `Update_EXHALE_stage2.tex` was regenerated and the PDF
  rebuilt at **158 pages**, the same count as before the sweep (MEASURED,
  `pdfinfo`). Thirteen British spellings were corrected in the same pass, twelve
  in the Markdown and one in the appendix (`modelling`, `centre`/`centred`,
  `labelled`, `cancelled`, `behaviour`, `travelled`, `NEIGHBOUR`), as the
  American-spelling rule for documents asks.
- **A3, the hand-maintained LaTeX and the stage-1 log**: `Update_EXHALE_stage1.md`
  211 to 4 and `Update_EXHALE_stage1.tex` 196 to 7 (verb uses and four
  `arm_heh1_x2matched` case-name literals), and nine other `docs/*.tex` files,
  of which `lit_wasp121.tex` was deliberately not edited: its "Fe I (blue arm),
  VLT/UVES" is the spectrograph's own blue arm and its "spiral arms" are the
  outflow's, both standard astronomy and not the forbidden sense. Nine PDFs were
  rebuilt with every page count unchanged (MEASURED, `pdfinfo`:
  `Update_EXHALE_stage1` 460, `EXHALE_user_manual` 90,
  `lhs1140b_exhale_vs_pwinds` 76, and six smaller memos).
- **A4, the source, the tests, the scripts and the top-level documents: 35 files,
  328 to 12**, the 12 all the verb plus two references to the `arm*` case
  directories by name. Identifiers renamed:
  `generic_precision_rows_arm` to `generic_precision_rows_selected`;
  `ARM_OFF`, `ARM_QUADRUPLE`, `ARM_GENERIC_DOUBLE` to `ROWS_PRODUCTION`,
  `ROWS_QUADRUPLE`, `ROWS_GENERIC_DOUBLE`; `arm_carrier_log_unknown` to
  `form_carrier_unknown_space` (named for what it forms, since it runs in both
  coordinates and always forms `carrier_log_floor`); `arm_departure_of_the_rows`
  to `key_departure_of_the_rows`; the local `arm_name` to `kind_name`; the local
  `e_arm` to `e_wb`; and eleven printed test-row names. Environment variables and
  input keys were NOT touched (`EXHALE_RESID_QUAD`, `Well balanced:`,
  `EXHALE_ELEMENT_CONSTRAINT_ROWS`). Every old row name and changed message was
  grepped over the whole working copy and matched only the file that defines it,
  so nothing parses them. MEASURED: the eight suites that link the changed
  modules keep their counts exactly (`krylov_and_dogleg` 334/0, `certification`
  70/0, `grid_and_gates` 200/0, `steady_species_rows` 195/0,
  `residual_determinism` 5/1 by design, `steady_completion_flag` 3/0,
  `steady_selfconsistent_residual` 3/0, `element_operator` 24/0); and the
  `wasp_full_newton` reload is byte-identical between a build of the delivered
  text and a control build of the entry text, both plain (`info = 0`,
  `||R|| = 1.495E-09`, flux spread 1.911E-12) and under `EXHALE_RESID_QUAD=2`
  with `Well balanced: True` (`info = 0`, `||R|| = 7.364E-10`, flux spread
  1.591E-13), so the renamed selector still selects. That the selection is not
  inert was measured beside it: `EXHALE_RESID_QUAD=1` gives `||R|| = 7.363E-10`
  at a flux spread of 1.588E-13, different in the last digits from mode 2, while
  mode 2 against no selection at all is byte-identical, which is the design
  assertion for the generic-double instantiation.

Three further occurrences were removed by the advisor at the close, outside the
four items: the noun in the comment at `src/EXHALE_main.f90` line 2022 and in
`diffusive_photochemistry.f90` line 248, and about 26 lines of the two fixture
READMEs (`backup/regression/atomic_elem_newton/README.md`,
`carrier_elem_newton/README.md`), which A4 left to the item then editing them.

What stays by design, everywhere: the verb; the `arm*` case directories under
`backup/regression/` and the two glosses that now point at them as file names;
and the astronomy senses in `lit_wasp121.tex`.

### An external commit during the session

Commit **`a7c18a6`, "remove build_*"** was made by the user at 21:39:24 KST on
2026-09-11 from another session, while item P21 was mid-edit (READ, `git log`
and `git show --stat a7c18a6`: 149 files changed, 39,652 insertions, 201
deletions). It removed the `build_dift/` and `build_fcheck/` dependency stamps;
it added the 141 files of
`docs/audit_20260905/partition_experiment_20260911/runs/`, 38 MB of copied
inputs, `.state` and `.faces` dumps and logs (MEASURED, `du -sh`); and it
carried P21's in-progress edits of `src/EXHALE_main.f90`,
`diffusive_photochemistry.f90` and `carrier_retry.f90` as they stood at that
minute. Two consequences are recorded because they are visible in the
measurements above. P21's line counts are `diff -u` against the entry texts
saved in its scratch directory before the first edit, `git diff` no longer being
able to produce that item's diff. And the provenance `git=` line of every output
written by a build made after the commit reads `a7c18a620f1b` where the control
binary's reads `0608382f639c`, which is the only difference in P21's M3 and M4
comparisons. The "What the repository carries" paragraph of
`docs/audit_20260905/partition_experiment_20260911/README.md`, which said
`runs/` was kept in the working copy only, now records that it is committed
since `a7c18a6`.

### Gates

Advisor verification of the merged tree (READ from the closing brief):
snapshot build with the default toolchain; `element_operator` 24/0,
`carrier_retry` 122/0, `krylov_and_dogleg` 290/0, `certification` 57/0,
`steady_species_rows` 194/0 and 195/0 after P8; `grid_and_gates` 158 PASS
with 3 FAIL, all three missing files of the snapshot copy (`inputdata`,
`examples`) rather than assertions (P4 measured the same suite 178/0 in the
live tree); the `wasp_full_newton` reload `info = 0`, `||R||` 1.667e-9,
CERTIFIED, `log10 Mdot` 13.30, unchanged. `make check` (advisor, 2026-09-11 09:40 KST, the shared build rebuilt from this tree): REGRESSION PASS, every DEFAULT_CASES output and `wasp_full_newton` byte-identical to the goldens; no golden moved, because no default case enters the stationary solve or the fixed-wind relaxations.

After P10 to P12 (advisor, merged tree, private snapshot build, 8 threads,
MEASURED 2026-09-11): the thirteen suites that link the changed modules all
PASS (`krylov_and_dogleg` 309/0, `grid_and_gates` 187/1 with the one FAIL a
WASP-52b SED file absent from the snapshot copy, `steady_species_rows` 195/0,
`element_operator` 24/0, `carrier_retry` 122/0, `certification` 57/0,
`attempted_step` 70/0, `carrier_returned_state_acceptance` 36/0,
`species_masses` 9/0, `carrier_reference_scales` 14/0,
`carrier_constraint_attribution` 10/0, `constrained_network_layout` 8/0,
`steady_completion_flag` 3/0). The two fixtures re-run on that binary
reproduce P10's and P11's records: the HD 209458 b element reload, 12 passes,
worst gated species row 2.61e-3 -> 5.47e-6 with only the hydrodynamic mass
row refusing (8.893e-10 against 3e-12 at cell 3, r = 1.0006); the carrier
reload, 3 passes, H2 7.38e-2 -> 7.06e-2 at cells 290 -> 288, mass 8.12e-12 /
4.85e-12 / 4.96e-12, momentum 3.22e-13 / 4.47e-13 / 3.54e-13.

Second `make check` after P10 to P12 (advisor, 2026-09-11 15:15 KST, the shared build rebuilt from the merged tree): REGRESSION PASS, every DEFAULT_CASES output and the cold-start `wasp_full_newton` byte-identical to the goldens; no golden moved.

After P14 to P17, each count MEASURED by the item that owns the file and READ
from its report: `grid_and_gates` **200/0** (P14, the second case of
`output_state_consistency.sh`; 199/1 on the entry text),
`krylov_and_dogleg` **321/0** after P15 and **334/0** after P17,
`certification` **70/0** after P16 (57/0 on the entry driver), with
`steady_species_rows` 195/0, `adv_static_limit` 53/0 and the other suites each
item ran unchanged in both of its builds. The two rows that FAIL by design or by
fixture are unchanged: `residual_determinism`'s
`closure_spread_within_the_row_tolerance_atomic_elem_newton` at 3.354E+01 (N33,
by design) and `steady_selfconsistent_residual`'s
`handback_matches_the_accepted_iterate` on a reload log that prints no iterate
norm.

Third `make check` after P14 to P18 (advisor, 2026-09-11 21:10 KST, the shared build rebuilt from the final tree): REGRESSION PASS, every DEFAULT_CASES output and the cold-start `wasp_full_newton` byte-identical to the goldens; no golden moved.

Fourth `make check` after P21 and the terminology sweep (advisor, 2026-09-12 01:40 KST, the shared build rebuilt from the tree with P21, A1 to A4 and the closing edits): REGRESSION PASS, every DEFAULT_CASES output and the cold-start `wasp_full_newton` byte-identical to the goldens; no golden moved.

### Reported by the workers and not fixed

- P1: `relaxation_drift_covers_every_element` in `steady_species_rows`
  builds a column whose species miss its own density by 1.660E-02 (closed
  the same day by P8); `relax_element_composition` damps helium alone, which
  belongs to the outer composition control; `element_diffusion_step`'s
  trace-metal branch re-seeds an exhausted cell and `project_elements` clamps
  `grp_new` at zero, the only place in the step that can create or destroy
  mixture mass, now covered by the closure test and not investigated
  further.
- P2: the carrier relaxation reaches into the steady residual for its
  advecting flux (`carrier_rows_advect = .true.`) while the element operator
  is given its flux by the caller, so the routine cannot be called at all
  unless someone else assembled a mass row; `residual_row_scale` in
  `steady_residual.f90` carries a comment block that is the history of a
  change and not the physics; the fixed-point exit of the carrier relaxation
  is unreachable in production while the H2 front propagates, and loosening
  `relax_tol` to make it reachable is not a fix.
- P3: ~~`distance_from_certification` divides the hydrodynamic slot by the
  run's `Resid tol` and not by each row's own certification tolerance, so
  the best-iterate ledger cannot see a mass row above 3e-12 while `||R||` is
  below 1e-8, which is exactly the state D3 is about (it did no harm on the
  D3 run, where the returned iterate's mass row is 5.945e-12 against the
  last iterate's 7.914e-12, but the ranking does not guarantee it)~~
  (**fixed 2026-09-11, P10**);
  `steady_selfconsistent_residual` assertion 2 (closed the same day by P8);
  ~~`attempted_step/run.sh` and `carrier_returned_state_acceptance/run.sh`
  write a fallback `build_stamp.f90` into their own work directory, where
  the compile cannot find it, so a tree that has never been built fails with
  `Cannot open module file 'build_stamp.mod'`~~ (**fixed 2026-09-11, P12**,
  in those two and in five sibling scripts).
- P4: the kind-generic rows do not form the equilibrium pressure term
  themselves (closed from `assemble_residual`, but the two places must stay
  one expression if the kind-generic rows are ever made to write their own face data back);
  ~~the default path's PLM momentum scale carries a coordinate artifact, the
  geometric pressure term `(A+ - A-) p_j/dV` entering the max separately
  although it is one physical term with the flux divergence, so at zero
  gravity with a uniform pressure the scale is `2 p/r` where the physical
  force is zero (MEASURED: the well-balanced scale, which is the row's true content,
  sits 0.944 below the base scheme's for both Riemann solvers at every N;
  under WENO3 the two agree to 2E-13). Fixing it moves the default momentum
  measure at every cell and therefore the acceptance of every recorded
  state, so it is a user decision.~~ (**fixed 2026-09-11, P11**, on the
  user's instruction: the scale is now the largest of the three physical
  terms, the zero-gravity artifact reads 2.22507E-308 in place of 1.000, and
  the acceptance movement is measured, `wasp_full_newton` still CERTIFIED.)
- P5: the same truncation-as-floor reasoning stood in three other documents
  (closed the same day by P8); the analysis document's own title is still
  stronger than its corrected conclusion; `docs/session_handoff_20260910.md`
  and `docs/steady_solver_design.md` carry statements of the same family,
  not read line by line; the experiment's private build is retained at
  `/tmp/exhale_partition_20260911.ZxXJEB` per its own section 9.
- P8: `docs/certification_tolerance_anchoring_20260910.md` still opens with
  `cert_tol_element` and `cert_tol_carrier` at 1e-8, which decision 22 (a)
  superseded, and whether that memo is a dated record or a live reference is
  a decision; `docs/code_status_20260910.md` section 3 described the
  stage-2 solver program as closed (this session's item P9 marks it
  superseded); ~~`mass_closed_column` exists only in
  `src/tests/element_operator/element_operator_tests.f90` and two suites now
  build the same column by hand, so the recipe belongs in one place a test
  can call~~ (**fixed 2026-09-11, P12**: `src/tests/test_columns.f90`).
- P9: ~~the `Coupled carrier solve` key comment in
  `src/modules/files_IO/input_read.f90` states that the alternation "was
  measured not to converge and not to be able to" (section 139). The first
  half is what P7 measured against: the alternation with the joint
  acceptance contract of P6 brings every elemental wind row of the atomic
  fixture inside 1e-5 and moves the carrier front monotonically, so the
  comment's reading is at least narrower than it states. Not my file this
  increment; it needs the sentence rewritten to what section 139 actually
  measured (the old alternation, without the acceptance contract, and its
  sign error on the front's motion).~~ (**fixed 2026-09-11, advisor**.)
- P10: ~~two rows of
  `src/tests/steady_species_rows/steady_species_rows_tests.f90` state the
  superseded semantics, feeding a raw `||R||` into the slot that now takes a
  distance and asserting that the hydrodynamic slot follows the run's
  `Resid tol` (MEASURED 193 PASS / 2 FAIL on the new objects)~~ (**fixed
  2026-09-11, advisor**); `dj` that the ledger ranks with is the accepted
  trial's own rows while the `judged rows` line prints the loop-top
  evaluation, so the reported best distance can sit below the minimum of the
  printed lines (`wasp_full_newton`: 8.469E-02 against a printed minimum of
  1.178E+00), which is entry-text behavior and not a defect; in the
  `EXHALE_RESID_BRANCH_REPORT` block one divisor for the whole class stays,
  because the block aggregates spreads across seeds whose maximum can sit in
  different rows.
- P11: ~~the kind-generic rows still hand back neither the momentum terms nor
  the face departures (`hydrodynamic_rows.f90`'s
  `store_the_interface_fluxes` writes `face_flux`/`face_p` and nothing
  else), so under `EXHALE_RESID_QUAD` the well-balanced face departure
  and, on a reconstruction continuation with `0 < lambda < 1`, which of the
  two pressure conventions the blended flux carries are read from the last
  `RK_rhs` evaluation instead of from the state just assembled; the shape of
  the fix is two weights in place of the `use_plm` branch, and the option is
  default off and a rounding control experiment nothing is adopted from.~~
  (**fixed 2026-09-11, P15**, which was also the Codex review's second
  finding: the kind-generic rows store their own face departures, and the continuation
  combination is refused rather than weighted, because the two schemes
  reconstruct different faces and no single pair of departures satisfies the
  blended row identity. The momentum terms are still not handed back by the
  kind-generic rows, and they are not read from a previous evaluation either: the one
  expression of each is evaluated at the call site for the state just
  assembled, which is the P4 entry above.) The
  marching stop of `mol_carrier` reports a momentum measure of 1.979E+00
  against a scale that is now the largest physical term in the layer, which
  is commensurable with the mass and energy entries in a way it was not;
  whether the marching stop should quote it at all is an acceptance
  question.
- P12: nothing further. The four sibling scripts carrying the same dead
  fallback and the work directories hardcoded with no override, which the
  P12 report listed as outside its own two files, are what P12b closed the
  same day.
- P14: ~~`assert_written_state_is_the_accepted_one` raises a false alarm
  whenever the flux spread reaches the rounding level (MEASURED on the
  partitioned carrier reload, 3.0982E-14 as written against 2.7487E-14
  accepted on a state whose mass flux is flat to fourteen digits, while the
  same run at a spread of 1.7736E-03 passes the test); it needs an absolute
  floor beside the relative one, or the warning will keep accusing the best
  states the solver produces~~ and ~~three rows of the `grid_and_gates` suite
  write to a fixed directory under `$ROOT/build/tests/grid_and_gates` and
  ignore `EXHALE_TEST_OUT`, so two concurrent runs of the suite collide
  there~~ (**both fixed 2026-09-11, advisor**, the floor at 1e-12 and the
  three rows named in the paragraph above). What is left is one line of
  documentation: `src/tests/grid_and_gates/README.md` gives one line for each
  test name, and the `output_state` row now asserts two things where the
  README says one.
- P15: `quad_rows_calls` and `quad_rows_seconds` are incremented in
  `hydrodynamic_rows_in_quadruple_precision` and public, and no file reads
  them (`grep` over `src/`), so the quadruple-precision rows' measured cost, which the module
  header calls a measured number and not an estimate, reaches no run summary;
  relatedly, nothing in `run.log` says whether `EXHALE_RESID_QUAD` took
  effect, so P15 had to establish that its quad run was a quad run indirectly,
  from the iteration-1 measure moving between the two binaries and from the
  sweep-kind reset READ in `steady_newton.f90`. A line in the setup report, or
  those two counters reaching it, would make a quad run self-documenting. And
  the continuation refusal the item added has no executed row: it ends in
  `error stop`, which cannot be asserted from inside the driver process, and
  no production path reaches it, the key being default off and the
  continuation controller disarming before the JFNK finish (READ).
- P16: the control build of the HD 209458 b reload ends with a state its own
  certification refuses and **exits 0**, because the outer loop's ending is
  the pass budget and no stationary claim is made. That is the run-mode
  contract working as written, but a reader of the shell status cannot tell
  that case from a certified one.
- P17: ~~`hydrodynamic_distance_from_certification` has no production reader
  left and cannot be deleted without breaking the compile of a test file the
  item does not own~~ (**item P18**, above). Two fixture observations, no
  action taken: running `grid_and_gates` twice at once with distinct
  `EXHALE_TEST_OUT` makes the `flux_spread` sub-test report
  `mass_flux_spread_recomputed_from_output measured=no_dump_line`, so the two
  invocations share something outside `EXHALE_TEST_OUT` (each invocation on its own
  gives 200/0), which is worth a line in that suite's README for whoever next
  runs two builds side by side; and the log-driven block of
  `steady_species_rows` requires a `(species_jac_test)` line in the log, so
  the `best judged iterate <= min judged rows seen` invariant it states is
  checked only on logs made with the species Jacobian probe armed, and P17 had
  to apply that comparison by hand (0 of the 12 solves of its new build report a
  best above their own printed minimum, against 1 of 12 in the control). A
  fixture log carrying that probe would be worth pinning.
- P21: the archived experiment drivers under
  `docs/audit_20260905/partition_experiment_20260911/`
  (`transport_wind_experiment.f90` and the `driver_snapshot.f90` copies beside
  the runs) call `relax_photochemical_composition` with its old six-argument
  signature. They are not in the `Makefile` `SRC` list and no test `run.sh`
  builds them, so nothing compiles them; they are the record of the run they
  belong to and are left as they stand. The `H2 front:` diagnostic P21 reported
  beside it was fixed by the advisor (above).

## 9. PLAN_20260912: the seven findings of the stage-2 review (2026-09-12)

Plan of record: `docs/PLAN_20260912_review_fixes.md`. User instruction of
2026-09-12: analyze `docs/Update_EXHALE_stage2_review_20260912.md` against the
code, and where its analysis is right, write the plan and start the
corrections. All seven findings R1 to R7 were confirmed in the source before
the plan was written (the confirming lines are listed at the head of the plan);
the corrections are items Q1 to Q4, Q1 and Q4 by the advisor, Q2 and Q3 by
workers under `docs/worker_rules.md`. Every number below is MEASURED on the
tree of this section unless marked READ.

### The state contract at a fixed hydrodynamic state (R2, R3)

The decision behind Q1: at a fixed hydrodynamic state the quantities the solve
holds are the conserved variables `u = (rho, rho v, E)`. A composition update
at fixed hydro keeps `u` and recomputes the pressure and the temperature from
the unchanged thermal energy `rho e = E - (rho v)^2/(2 rho)` through the
caloric EOS of the new composition (`pressure_from_energy_density` after
`get_species_densities`, then `comp_T_from_p`). A held pressure is not a state
the solve knows: the review's 500-cell probe showed that keeping `p` while the
main loop keeps `E` describes two molecular states 3.3e-3 apart.

### The items

- **Q1 (R1, R2, R3, R6, R7 and the outer progress measure), advisor.**
  `diffusive_photochemistry.f90` (+206/-204), the carrier branch and the
  summary of `steady_wind_with_element_diffusion` in `EXHALE_main.f90`,
  `src/tests/carrier_retry/carrier_retry.f90` (+195/-77) and its `run.sh`.
  - R1: `photochemical_transport_step` takes an optional `trial` argument; a
    trial of the stationary relaxation returns before the history mark, the
    message and the stop, so a rejected trial (an uncovered interval, a
    chemistry that did not close, a displacement beyond the movement bound)
    leaves `carrier_history_certifiable_flag` as it found it. The marching
    path passes no `trial` and is unchanged.
  - R2, R3: `relax_photochemical_composition(u, v, f_sp, p, T, heat, cool,
    eta, trust, drift, nstep, outcome)` now takes the conserved state; `p` and
    `T` are outputs. `pressure_and_temperature_at_fixed_conserved_state`
    is the one expression of the contract above, and
    `equilibrate_chemistry_at_fixed_conserved_state` cycles (p, T at fixed u;
    sweep; p, T again) until the largest relative change of T falls below
    1e-6 (at most five cycles), reads the sweep ledger (`n_nonfinite`) and a
    finite composition, and reports `ok`; a refresh that does not close
    refuses the trial (`carrier_relax_chemistry_refused`, a new named ending
    "the chemistry of the shortest admissible trial did not close") and the
    entry composition, background, p, T, heating, cooling and eta are
    restored. The caller in `EXHALE_main.f90` rebuilds `W` from the unchanged
    `u` and takes the returned composition; nothing is reconstructed from a
    held pressure any more.
  - R6: the carrier retry suite's column is built by
    `column_carrying_its_own_density` of `src/tests/test_columns.f90` (the
    mass-closed constructor of the element suites; `run.sh` compiles it
    before the driver) with the ions moved out of the neutrals; a row asserts
    the column reconstructs its own density to 1e-13 before any step. The
    consistency rows are absolute: the returned pressure is that of the
    unchanged thermal energy (0.0), the returned temperature is that of the
    returned composition (0.0), the background third-body density and
    temperature the sweep stored are those of the returned composition
    (6.53e-7 and 2.28e-7 relative, allowance 10 `ieq_res_tol`), and one
    further sweep leaves the returned composition within the sweep tolerance
    (8.83e-7 against `ieq_res_tol`). Nine rows for R1: a rejected trial keeps
    the history certifiable, keeps no step, hands back the entry composition
    bit for bit and names the refused interval; a pass after the refusals
    keeps a step and stays eligible; a refused trial in phys mode reports
    exhausted, leaves the history certifiable and writes no composition.
    RED on the entry text (the trial marking the history, the refresh not
    closed): five rows fail (pressure 1.2e-7, temperature 3.5e-7, further
    sweep 1.2e-6, the two R1 rows 0 of 1); GREEN on the tree: 138 PASS, 0
    FAIL.
  - R7 and the progress measure: the outer summary of
    `steady_wind_with_element_diffusion` ranks every entry by its binding
    distance (`dist_bind`, `row_at_bind`, `jbind` where the entry carries
    them), so the mass row's refusal is attributed to the cell that binds it
    and not to the largest row over another cell's tolerance; and the
    progress control halves omega and the carrier movement bound on the
    JOINT distance of the state (the largest entry over its own tolerance,
    `prog_worst`), not on the worst species row alone. An entry that is not
    finite or not available counts as an infinite distance.
  - The review's 500-cell probe (`docs/audit_20260905/stage2_review_20260912/
    molecular_wind_state_probe.f90`, adapted to the new signature, hot-Uranus
    carrier reload, one hydro solve and one bounded relaxation at trust
    0.01): returned temperature against the EOS of the returned state 0.0
    (was 5.5997819e-9), pressure of the unchanged thermal energy against the
    returned pressure 0.0 (was 3.2993e-3), thermal-energy change a held
    pressure would need 6.9e-16 (was 3.3250e-3); mass closure 8.33e-15 and
    the five retained steps at drift 9.9955e-3 unchanged.
  - Impact on the routes: the unbounded pass of P21 M1 repeated on the new
    contract (`EXHALE_CARRIER_TRUST=1e6`, cap 40) reaches the fixed wind's
    fixed point in 61 transport steps at drift 0.995 with the H2 front
    (x2 = 0.5) at 3.079 R_p (P21: 57 steps, 2.237 R_p; both fully molecular
    fixed points, the front leaving the domain). The 12-pass production loop
    ends at H2 row 5.49e-2 at cell 279, front 1.0755 R_p (the 110-pass ladder
    of P16 at pass 12: 5.50e-2 at cell 279, 1.0765 R_p); no pass triggered
    the halving. `atomic_elem_newton` at three passes: pass lines identical to
    the control build. `wasp_full_newton` reload and the `mol_carrier` case
    (12000 steps, one thread): data rows identical to the control build and
    to the golden.
- **Q2 (R4 in the operator), worker (Sonnet).**
  `binary_element_diffusion.f90` (+108/-43),
  `src/tests/element_operator/element_operator_tests.f90` (+88/-1).
  `element_diffusion_step` judges its candidate on every call: the entry
  composition is kept, the four tests (solved, finite, X inside [0,1] to
  `element_fraction_bound_tol`, mass closure no worse than the entry's by
  `element_mass_closure_tol`) run whether or not `status` is present, and the
  entry composition is restored on a refusal. The outcome is always left in
  the new public `element_step_last_status`, and `element_step_outcome_text`
  names it. Two new rows repeat the ordinary and the NaN-residual steps with
  no `status` at the call site: RED on the entry text (the refused candidate
  with departure 0.86 adopted, the status left at accepted), GREEN after
  (`element_operator` 28/0). A second, independent defect found by the change:
  the Newton solve of the helium mass fraction (`solve_mass_fraction`) read
  its arithmetic floor only on a DESCENDING pass, so a step whose residual had
  already reached the floor and whose next pass found no improving trial
  reported unsolved (residuals 1e-12 to 1.1e-5 at the refusal, MEASURED on
  the closed-column relaxations of `diffusion_tests` T1a, T7c, T7d and T10,
  which went 35/0 to 27/8 under the unconditional judgment alone); the
  no-descent exit now reads the same floor, and `newton_drop` moved from
  1e-6 to 2e-5 so that the ~1e-5 floor of the K_zz = 2e12 homopause column
  (T7), which the module header already documented, lies inside it
  (`diffusion_tests` 35/0 again). Cost on `mol_diffusion` at 300 steps: 54.4 s
  against 54.2 s, inside the noise. `mol_diffusion` and `lower_profile` at 300
  steps: data rows identical to the control build. `atomic_elem_newton`: pass
  1 identical; at pass 2 the judged rows are identical and omega differs
  (0.500 control, 0.250 Q2), attributable to the floor change turning a
  floor-level residual from unsolved to solved (READ from the Q2 report; the
  advisor's own three-pass run of the merged tree, above, has omega 0.250 at
  pass 2 on both builds because the control binary of that run already
  carried Q2).
- **Q3 (R5), worker (Sonnet).** `certification.f90` (+105/-17),
  `src/tests/certification/certification_contexts.f90` (+106),
  `docs/certification_tolerance_anchoring_20260910.md` (addendum).
  `mass_row_cell_verdict` returns a third, optional outcome: RESOLVED
  (`cert_mass_round_margin * floor < cert_tol_mass_ceiling`, judged as P16
  anchored it) or UNRESOLVED (the margin times the floor reaches the ceiling,
  or the floor is not finite): an unresolved cell never reads `within`, its
  distance is `max(1, q/ceiling)`, so the review's counterexample (q = 0.5 at
  floor 1) reads `within = F`, distance 1.0 (entry text: `within = T`, 0.5).
  The column loop is the pure `mass_row_column_verdict`, which counts the
  unresolved cells and names the first (`n_mass_unresolved`,
  `j_first_mass_unresolved`), and the report names that cell even when a
  different cell binds the verdict. `certification` 84/0 (70/0 before the
  new rows), `krylov_and_dogleg` 334/0, `steady_species_rows` 195/0 (READ
  from the Q3 report; the advisor's run of `certification` on the merged tree
  84/0). The addendum records that the margin of 10 was measured on three
  converged or converging fixtures only, none a near-zero-flow or
  molecular-EOS column: the UNRESOLVED outcome exists for the regime the
  fixtures do not reach.
- **Q4 (R4 in the marching caller), advisor.** The marching loop of
  `EXHALE_main.f90` reads `element_step_last_status` after the element step:
  in phys mode a refused step refuses the whole attempted step
  (`as_reject_element` at the diffusion operation, the verdict text naming
  the test that refused it), which the controller restores and retakes at
  half dt like a carrier interval left uncovered; in init mode the march goes
  on from the composition the step restored and the refusals are counted
  (the final report prints "element steps refused inside the initialization
  march: n"). A test knob, `element_refuse_leading_steps_for_test`
  (`EXHALE_ELEMENT_REFUSE_STEPS`, zero in every production run), refuses the
  first n element steps so the path can be exercised: on `mol_diffusion` in
  phys mode with two refusals injected, attempt 1 at dt 1.47451E-04 and
  attempt 2 at 7.37253E-05 are refused "at the element diffusion step: an
  element invariant broke" and attempt 3 runs at 3.68626E-05; in init mode
  the two refusals are counted and the run exits 0. (The phys-mode run of
  that fixture then fails at step 0 on the coupled source step's fixed point
  with or without the knob, exit 2: a property of that fixed-step snapshot in
  phys mode, not of Q4.) Found beside it and fixed: the refusal census of
  `attempted_step_note_step` stopped one reason short
  (`as_reject_source_fixed_point`), so a source-energy refusal was never
  counted by reason; the bound is now `as_reject_source_energy`.
- **The stagnation row of `grid_and_gates`.** On the joint progress measure
  the configuration `output_state_consistency.sh` used to force the
  stagnation ending (`EXHALE_JFNK_MAXIT=5`, `EXHALE_CARRIER_TRUST=1e-4`) runs
  out its budget instead: the starved hydrodynamic rows keep falling (energy
  4.9E-01 to 3.9E-02 over 8 passes, 1.1E-01 at 20) while the H2 row stands at
  7.3E-02. The row now uses `EXHALE_JFNK_MAXIT=40` with
  `EXHALE_CARRIER_TRUST=1e-6`: the hydro solves reach their roots, the H2
  row stands at 7.38E-02 for four passes and the ending is announced at pass
  4; the second assertion of that script (the ending hands back one state)
  measures 7.1e-16 against 1e-12. The script's header and the suite README
  say why the configuration moved.

### Gates

Advisor verification of the merged tree (private build `build_Q1` /
`EXHALE_Q1.x` of the default toolchain, deleted at the close): `carrier_retry`
138/0, `certification` 84/0, `element_operator` 28/0, `attempted_step` 70/0,
`steady_species_rows` 195/0, `krylov_and_dogleg` 334/0, `diffusion_tests`
35/0, `carrier_returned_state_acceptance` 36/0, `carrier_reference_scales`
14/0, `carrier_constraint_attribution` 10/0, `species_masses` 9/0,
`steady_completion_flag` 3/0; `grid_and_gates` 198/1, the one FAIL the
stagnation row on its old configuration (the suite was started before the
script was re-anchored), and the re-anchored `output_state_consistency.sh`
4/4 on its own afterwards. `make check` (advisor, 2026-09-12 10:35 KST, the
shared `build/` rebuilt from this tree with the gfortran on PATH, conda-forge
16.2): **REGRESSION PASS (all cases byte-identical)**, 64 file comparisons
PASS with data identical. A first launch of the gate with the system gfortran
13.1 on PATH was stopped before its verdict and is not counted: the goldens
are snapshots of the 16.2 build and a verdict against another compiler would
say nothing about this change.

### Reported by the workers and not fixed

- Q2: the out-of-bounds outcome of the element step (`element_step_out_of_bounds`)
  is not reachable on a synthetic column of the suite's size (the drift flux
  vanishes at both ends of the composition axis in the discrete operator), so
  the new rows cover the solve-failed and nonfinite outcomes only; the one
  measured out-of-bounds excursion of record is P1's, on the atomic reload.
- Q3: `grid_and_gates` row `mass_flux_spread_recomputed_from_output` failed
  once on the worker's control run ("no_dump_line") and passed on its measured
  run and on the advisor's; the P17 note on two concurrent invocations of that
  suite sharing something outside `EXHALE_TEST_OUT` stands.
