# EXHALE: state of the code, what physics is in and verified, what is certified, where the problems begin, and what comes next (2026-09-10)

Written to restart from in a fresh session. Read this first; then
`docs/ISSUES_20260909.md` (the problem-centered account of stage 2: what was
done, every problem resolved, every problem open with its evidence),
`docs/PLAN_20260909_rev1.md` (the plan those items belong to),
`docs/session_handoff_20260910.md` (the last handoff), and
`docs/Update_EXHALE_stage2.md` section 7 (the dated record of items B5e to N30).
This document replaces `docs/code_status_20260905.md`, which is kept as the
state at the opening of stage 2; every row below is the state after item N30.

Every statement was verified against the source or against a measurement of
the sessions of 2026-09-05 to 09-10. A number is labeled MEASURED (this
document ran it), READ (from the source or from a file) or LOGGED (copied
from `docs/Update_EXHALE_stage2.md` or `docs/ISSUES_20260909.md`, which state their
own measurement). Statements carried over from before 2026-09-05 without
re-measurement are marked "per the earlier record".

## 0. The one-paragraph verdict

**The atomic line, H, He (with the 2^3S metastable) and the ten trace metals,
with Roche or spherical gravity, the two-stage march and the Newton finish,
and the transit post-processing, is implemented and validated** (Huang et al.
2023 WASP-121 b reproduced in mass-loss rate, metal cooling, n=2 physics and
transit depths, per the earlier record), and its three-unknown stationary
solve is the one route in the code that converges and CERTIFIES. **The
molecular line, H2, H2+, H3+, HeH+, the oxygen cycle, the molecular carriers
and the transported ionization state, is implemented, regression-guarded and
only partly validated.** **No stationary solve that carries a species row,
an element row from `He_diffusion` or a carrier row from the molecular
carriers, converges on any route**, and both species-row certification
tolerances were anchored by measurement only on 2026-09-10 (decision 22 a,
item N30): 1e-5 gating at r >= 1.20, the rows below reported. The two
candidate states stand at 2.9e-4 (atomic element reload, wind) and 7.3e-2
(carrier reload) against that gate, so neither is certified. The problems begin
exactly there (section 4).

## 1. The tree

* Branch `v1.00`, HEAD `35d9dd5` ("Heating and cooling assembled from the
  post-sweep composition"). **367 paths are uncommitted** (274 modified, 88
  untracked, 5 deleted; MEASURED `git status --porcelain`); nothing since HEAD
  has been committed or pushed. `backup/` (regression harness, goldens) is not
  in the git remote.
* Toolchain: plain `make` builds with the conda-forge gfortran 16.2.0 on PATH
  (MEASURED, `/opt/miniconda3/bin/gfortran`) and the OpenBLAS of the same
  prefix (rpath'd); the compiler's path and version are part of the rebuild
  stamp. `make FC=ifx OBJDIR=build_ifx EXE=EXHALE_ifx.x` builds with ifx
  2026.1 + MKL. Do not prepend `/usr/bin` to PATH.
* Regression: `backup/regression/run_check.sh check`, **sixteen default cases**
  (`wasp_full`, `wasp_he23off`, `wasp_full_newton`, `mol_base_handoff`,
  `mol_metals`, `mol_lyman_werner`, `mol_diffusion`, `mol_ir_bands`,
  `mol_sec_ion`, `mol_carrier`, `lower_profile`, `hydrostatic_column`,
  `oxygen_chemistry`, `hp_zero_seed`, `hp_trace_seed`, `hp_front`; READ,
  `DEFAULT_CASES`), plus the `roundtrip` golden. Bitwise first, else within
  `REGRESSION_REL_TOL` (1e-3). **The goldens are current**: fifteen cases
  re-snapshotted at the Stage D gate and `roundtrip` with them (MEASURED file
  times 2026-09-10 00:30 and 00:15; previous set
  `golden_pre_stageD_20260910/`), and `lower_profile` re-snapshotted after
  N29 (MEASURED 2026-09-10 11:06; previous set `golden_pre_n29_20260910/`).
* Two pinned reload fixtures live outside the matrix and are not goldens:
  `backup/regression/wasp_full_newton/IC/` (the certified atomic state, read
  by the `restart_intent` rows of `src/tests/grid_and_gates/`) and
  `backup/regression/atomic_elem_newton/` (HD 209458 b with eight element
  rows, the diagnostic entry point of the atomic element reload, with its own
  `IC/` reload pair).
* Assertion suites: `make test` builds and runs `element_census_tests`,
  `diffusion_tests`, `residual_determinism` and **every executable
  `src/tests/*/run.sh`, of which there are 25** (MEASURED): `acceptance_classes`,
  `adv_static_limit`, `attempted_step`, `carrier_constraint_attribution`,
  `carrier_reference_scales`, `carrier_retry`,
  `carrier_returned_state_acceptance`, `certification`,
  `constrained_network_layout`, `coupled_source_step`, `element_operator`,
  `energy_update`, `fuv_band_ledger`, `grid_and_gates`, `krylov_and_dogleg`,
  `physics_probe`, `residual_determinism`, `run_mode`, `species_face_flux`,
  `species_masses`, `spectrum_type`, `steady_completion_flag`,
  `steady_selfconsistent_residual`, `steady_species_rows`, `transit_census`.
  Eight of them carry a `README.md`; the rest state what they compare in the
  driver's header. `run_fcheck.sh` (a `-fcheck` build on one bounded
  HD 209458 b case) was last CLEAN at the 2026-09-07 gate and has not been
  re-run since (LOGGED, ISSUES 3.8).
* Documentation: `docs/Update_EXHALE_stage2.md` (stage 2, sections restart at 1;
  section 7 holds every item from the PLAN rev2 series to N30) and
  `docs/Update_EXHALE_stage1.md` (26,086 lines, sections 1-171); 146 memos in
  `docs/`; `docs/input_schema.md` (every input key, keys K43 to K45 added by
  N10, N10b and N28); `README.md`, `README_HOWTO.md`,
  `docs/EXHALE_user_manual.tex`.

## 2. Physics inventory

Legend: **V** verified against a published solution or an analytic limit;
**I** implemented and regression-guarded, not validated against an
external result; **L** implemented with a known limitation stated in the
row; **X** not implemented.

### 2.1 Hydrodynamics and the steady state

| piece | status | evidence / limitation |
|---|---|---|
| Godunov (HLLC) + RK3, PLM and WENO3 reconstruction, Mixed grid, runtime cell count | V | ATES heritage; the regression matrix. The WENO3 stencil coefficients were swapped until item WENO3-ORDER (`C2(j)` on `dWp` in WL, `C1(j-1)` on `dWm` in WR); the measured convergence rate moved 2.991 to 3.000. The positivity limiter clamps to `epsilon*W_avg` and can no longer return a zero density or pressure |
| gravity: Roche potential truncated at L1 (default) or spherical to a stated radius (`Domain mode: Spherical`, `Outer radius [R_p]`) | V | Huang 2023 Case B (Roche, per the earlier record); Koskinen 2022 Model A velocity and neutral densities (spherical) |
| lower boundary: fixed rho or p at the base, isothermal ghost, base velocity from the mass flux; characteristic face condition | L | the molecular base "breathes": a marching limit cycle at the default CFL (absent below CFL 0.15-0.3). It keeps `du` from converging in molecular runs and distorts fixed-step snapshots. Separately, the base cell amplifies any composition inaccuracy by about 4e3 through `Apply_BC`'s continuous-temperature ghost (ISSUES 3.3) |
| outer boundary: outflow. The upper ghost of a transported element or carrier column follows the iterate (zero gradient), the marching path's own rule | I | N26. Before it the ghost was frozen at the pre-solve composition and the outermost element row was unconstrained from outward; the change moved the atomic element reload from `\|\|R\|\|` 1.36e-2 to 3.7e-4 (LOGGED) |
| viscosity, conduction (Watson et al. 1981 kappa) | I | keys `Viscosity`, `Conduction`; conduction negligible on Model A, as they state |
| caloric EOS: H2 rovibrational ladder, atoms/ions/electrons 3/2 k; `Caloric EOS: monatomic` (comparison option) | I | the Koskinen gate's outer temperature is insensitive to it (20 K) |
| the material advective term: ONE discretization shared by the RK stages, the stationary rows and both relaxations | I | B4-1, B4-1c |
| convergence: two-stage `du` (flux spread) hand-off, JFNK Newton finish with the scaled trust region, PTC | L | see section 3.3. The atomic three-unknown solve converges and certifies (`wasp_full_newton`: info 0, `\|\|R\|\|` 4.3e-9). A `du` stop alone leaves Mdot path-dependent at the percent level: **always Newton-finish an atomic run before quoting Mdot** |
| the physical step's integration error | L | one transaction over one macro-interval, the step-doubling estimate rejects, the clock moves only at adoption (N16a); the order of the complete split update is MEASURED as 1. Still open: `err_rtol` 1e-3 and `err_atol` 1e-12 are flat and unvalidated and cover no transported species, and a retried macrostep's accepted state depends on which pass was rejected (up to 9.1e-5 relative in the base cell's velocity). Until that is closed, no physical-mode trajectory is quoted as error-controlled for its species (ISSUES 3.8) |
| restart contract: `Restart intent: trajectory \| relaxation \| stationary [evaluate\|equilibrate]`, an eight-line `restart_schema 1` metadata block in both state files, refusal on a disagreeing grid, reservoir, constants or options, `Restart option change: <token>[, ...]` to permit named option differences | I | N10, N10b. A file with no block is legacy, loaded as before and marked `provenance unknown`, and the mark is inherited by everything the run writes. Six layout tokens (`metals`, `mol`, `oxychem`, `carrier`, `carrier_newton`, `iontrans`) may never be named: they decide the rows the state carries, so changing one is a cold start, not a restart |
| OpenMP over cells (ionization sweep, cooling, hydro stages, carrier Jacobian) | V | one binary gives the same bits at every thread count (THREAD-DET; the 6.6e-14 difference between 1 and 16 threads is gone, fixed block length 32) |

### 2.2 Radiation and heating

| piece | status | evidence / limitation |
|---|---|---|
| spectrum: one type builds every band (power law, Planck, or a loaded table), the XUV and the part below 13.6 eV alike; the Balmer continuum follows it | V | decision of 2026-09-05. Before it the sub-Lyman field was 500x too faint for He 2^3S. A loaded table that stops above the lowest threshold of an active absorber stops the run |
| photon grid: every active threshold is a bin edge | V | 2c-QUAD, 2c-METALEDGE. Threshold bins had been integrated across an edge, +9.5 percent in the H I photoionization rate |
| dayside average `2D approximate method: Rate/4 + Mdot` | V | the same as Model A's "flux divided by 4" |
| photoionization H, He, He+, He 2^3S, metals (Verner et al. 1996 and Verner & Yakovlev 1995; opacity tables optional) | V | the seventeen metal stages take each fit only where it is valid and carry the SHELL SUM, outer plus subshells, reproducing Verner's own `phfit2` to 1e-10 at 17 ions x 12 energies. Before item REF-METALS the fits were evaluated above their `E_max` (Fe I 49x the truth at 1240 eV) and the inner shells were absent; the correction raised base heating 21 to 26 percent on the metal cases |
| metal inner-shell absorption: the ion advances one stage, one electron of h nu - E_th(outer) to the cascade | L | decision 7 b. The true Auger event leaves the atom two or more stages up, which the three-stage metal ladder cannot hold; the one electron delivers 0.95-1.09x the heat and 1.10-1.37x the secondary H I ionizations of the two the real event has. Fluorescence neglected and quantified (Krause 1979 Tables 3 and 5) |
| H2 photoabsorption: four channels (H2+; H+ + H; 2H+; neutral dissociation window) sharing one cross section | I | no external check; `f_di` from Chung 1993 Table II |
| photoelectron energy: h nu - I heats (default); secondary ionization (Shull & van Steenberg, staged after first convergence, or `Immediate`, or off) | V/L | standard for the atomic line; the Koskinen comparison shows Model A deposits 2.5-3x h nu - I at the same ionization rate, a property of their code. `Photoelectron heating: full` or `<fraction>` are comparison-only options |
| He recombination photons ionizing H (`He_rec_coupling`, default on) | I | off for the Model A runs |
| H(n=2): Balmer photoionization, Ly-alpha pumping, in-situ Ly-alpha emissivity from LaRT | V | Huang 2023 phases 3a/3b, per the earlier record |
| Ly-alpha escape probability: `beta = 1/(1 + 0.909316 B)`, Neufeld 1990 eq. 3.27, face partition eq. 2.25, `A_2p1s` = 6.2649e8 | V | item LYA-BETA corrected the coefficient; `n_2p` moved by up to 84x at the base, so **every H-alpha, H-beta and Mg II transit number in the repository from a `Jlya escape-prob: True` run is superseded** |
| Lyman-Werner H2 photodissociation, band 912-1201 A, Cloudy-derived self-shielding table, fluorescence heating | I | decision 6 B. The table and the `Stellar LW flux` key had been normalized on different intervals; the LW-only regression case moved 343.0 to 480.9. Default off unless a flux is stated or a spectrum is loaded; a stated zero switches it off |

### 2.3 Ionization, chemistry, cooling

| piece | status | evidence / limitation |
|---|---|---|
| coupled H/He (+ He 2^3S, + metals, + molecules) local-equilibrium solve per cell (`System_HeH_*`, analytic-Jacobian Newton, hybrd1 second attempt) | V | Huang 2023 phases 2-3 |
| rates: Badnell case B recombination, Voronov collisional ionization; `Atomic rate set: Koskinen2022` | V | the switch moved the Koskinen gate by less than 5 percent. The Ca I radiative recombination row was 6x off until REF-METALS |
| metals C, N, O, Mg, Si, Ca, Na, K, S, Fe; closed-form CHIANTI cooling; `eos_metals` in the mass/electron/particle budget | V | Huang 2023 WASP-121 b Fe II / Mg II cooling and Mg II, Na I, Ca II transit depths, per the earlier record. An element `metals.inp` does not carry is no longer registered as a stationary unknown (N1) |
| He 2^3S: Penning ionization (Taylor 2025 fit), triplet cooling in the energy solver, charge exchange, He 10830 transit | V | Falorca 2026 review items closed, per the earlier record |
| molecular network H2, H2+, H3+, HeH+ (Koskinen 2022 Table 1 R1-R23); collisional reaction heat (default on); H2 caloric EOS | I | the Koskinen gate; the H2 loss rate of one molecule at 1.5 r_base agrees with Model A (1.0e-5 against 1.2e-5 s^-1) |
| oxygen cycle OH, H2O, CO, four FUV photolysis bands, CO reservoir, `Oxygen_chemistry.txt` | I | no external check. CO now has destruction channels (`He+ + CO` 1.6e-9, UMIST RATE22 4068, with its He+ sink; shielded photodissociation on the LW beam, Visser et al. 2009 Table 6) and the thermodynamic ceiling that stood above the layer is deleted |
| cooling: recombination, collisional excitation (Ly-alpha and He), bremsstrahlung, H3+ (Miller 2013), H2/H2O/CO infrared bands, the infrared field of the atmosphere below | I/L | H3+ per molecule is about 3x what Model A's Figure 9 implies at 3 r_base. The fine-structure line escape argument was off by sqrt(pi) until REF-LODI (de Jong 1980 eq. B-7, frequency-integrated depth), and the infrared tabulation grid ended below the wind temperatures until it was extended to about 30100 K |
| the sweep returns the heating and cooling OF THE COMPOSITION IT RETURNS | V | fixed 2026-09-05 |
| acceptance of a cell: five classes, one definition, with a cap on an amnestied cell (`ieq_nonroot_res_cap` = 1e2, class 6 keeps the entry composition), and a class-4 cell written but reported "NOT certified" | I | ACCEPT-2, decision 2 a; A1 |
| the coupled source loop stops on the ESTIMATED error of the accepted state, with a measured safety factor 3 | V | COST7, decision 8. Stopping on the last increment left a third of the accepted states outside `csm_comp_tol` |
| the unswept background carries no gas, no rates and no imposed partition | V | N28b: `ion_rates` carries a default zero in all 37 fields. Three grid-sized objects had come out of `allocate` holding a recycled heap block, including a LOGICAL `x_h2_fixed` read TRUE from garbage |

### 2.4 Transport of the composition

| piece | status | evidence / limitation |
|---|---|---|
| advection correction of the ionization state and the temperature in post-processing (`_adv` profiles) | V/L | see 2.6 and section 4. The energy equation is now the exact steady one with the enthalpy flux term restored (ADV-STATIC, ADV-ENERGY), and each row says whether its temperature and its composition were corrected or retained |
| carrier transport operator in the hydro: implicit advection-diffusion of H2 (and OH, H2O, CO), 2nd-order limited, base inflow, Blanc's-law molecular diffusion, eddy `Kzz` | I | Picard alternation with the wind |
| `Ionization transport: True`: H+ as a fifth carrier, the sweep handed the transported fraction | V | reproduces the advected estimate where P r/v < 1; the Koskinen gate's n_e excess fell from 5-10x to 1.5x. **Molecular branch only** (requires `Molecular chemistry` and `Molecular carrier transport`, refused otherwise) |
| He/H binary element diffusion (`He_diffusion`); trace metals diffusing independently | I | the element projection took the mixture mass from the RESERVOIR metal/H rather than the cell's and created mass at every call (6.4e-3 on the atomic candidate) until N29; with the cell's own mass the closure is round-off (1.9e-9, the file's precision) and the writer-to-loader round trip returns the state. `lower_profile` moved (T 1.3 percent, H I 6.3 percent at 1.6 R_p) and its golden was refreshed |
| **He+, H2+, H3+, HeH+ (and every metal stage) are local equilibrium in the hydro** | X | where the Koskinen comparison's remaining ionization excess sits: He+ 3x theirs above 1.5 r_base. A local root overstates an ion where P r/v < 1, as it did for H+ |
| **a stationary solve that keeps a transported H+, or any other transported balance** | L | the row EXISTS: the stationary Newton carries a row and an unknown for every transported balance (B5, B5b, B5c), and `Ionization transport` with `Solver: Newton` is supported when `Coupled carrier solve: True` is set and refused without it. What does not exist is a CONVERGED one: no species-row solve converges on any route (section 4.1) |

### 2.5 Lower atmosphere and hand-off

| piece | status |
|---|---|
| Koskinen 2022 analytic column (Tier 1); `base.inp` scalar handoff (`p_base`, `T_base`, `HeH_base`, `Kzz_base`, `q_H2_base`); lower-atmosphere PROFILE handoff from VULCAN/photochem (Tier 3, `Lower atmosphere profile:`, the `lower_profile` regression case) | I; the Koskinen runs use the scalar handoff at 1 microbar |
| a coupled steady solve on a molecular configuration WITHOUT `Molecular carrier transport` | refused at startup (decision 10 a): an eliminated H2 content is one the local equilibrium cannot determine, so the stationary residual would not be a function of its unknowns |

### 2.6 Post-processing

| piece | status |
|---|---|
| transit spectra (He 10830, Ly-alpha, H-alpha/beta, Mg II, Ca II, Na I, O I 1302/1304/1306), Roche geometry, LaRT coupling (`EXHALE_transit.py`) | V for the atomic line (Huang 2023 phase 5, per the earlier record). Every `tpm_*.txt` now carries a `transit_schema 1` metadata block: the input profile and its `_adv` schema, the certification and provenance lines copied from it, the tool's own identity, the refused-row census of that line with its eight-point ladder, and the `EXHALE_TRANSIT_*` overrides in effect (N13, N13b, decision 17 a). Still line center only, and the refused share is taken over the temperature field alone |
| `_adv` profiles | L, a CONDITIONAL correction with its condition written in every row: the stationarity of the input is judged by the certification's own mass operator, the measure `adv_mass_row` is a column, the decision is taken against `adv_conditional_tol` = 1e-2, and `adv_T_status` / `adv_comp_status` say per row whether the temperature and the composition were corrected (0), retained (1), failed (2), unsupported (3) or not evaluated (4). `# adv_input_certified` is the separate certification verdict, and the product is declared a ONE-WAY correction on a fixed density and velocity field. Molecular and oxygen columns are the equilibrium sweep's values |

## 3. Validation and certification status

### 3.1 The atomic + metals line: implemented and validated, with four qualifications

Reproduced (per the earlier record): Huang et al. (2023) WASP-121 b, Mdot
(0.056 / 0.337 for Cases A / B against 0.052 / 0.32), the Fe II / Mg II
cooling structure, the n=2 proton source and Balmer heating, the Mg II,
Na I, Ca II transit depths; Newton-converged winds for HD 209458 b,
HD 189733 b, WASP-52 b, WASP-121 b.

1. **Those runs predate the physics of stage 2.** The Ly-alpha escape
   probability, the metal photoionization shell sum, the photon-grid
   quadrature, the cell width `dr_j`, the Jupiter radius and the spectrum
   assembly all moved. The re-convergence of the planet folders and
   `benchmarks/` is user-gated (section 5.9).
2. **Advected ionization is post-processed, not in the hydro**, for the
   atomic line. Where P r/v >> 1 in the wind, the hot Jupiters validated,
   local equilibrium is right; the in-hydro `Ionization transport` exists
   only with molecular chemistry on.
3. **The paper and poster tables** stand on a superseded binary
   (`docs/paper_materials_vintage_audit.md`); re-convergence is
   instruction-gated.
4. The photoelectron accounting (h nu - I) is the standard one and is what
   the validation used.

### 3.2 The molecular line: implemented, partly validated

Target: Koskinen et al. (2022) Model A. Within the plotted points: velocity,
H2/H/He densities, the layer temperature (1.05-1.30 r_base, within 100 K
after 4e5 steps). Decided by the heat deposited per ionization: the outer
temperature, velocity and ionization. Their Figure 9 has 2.5-3x the h nu - I
rate at the same ionization rate; with a stated 0.55 of the photon energy
(comparison only) the heating and temperature match, and n_e stays 1.4-1.6x,
H3+ 0.35-0.5x, He+ 3x, f(H2) 0.06-0.08 low. Mdot 0.66-0.78 of theirs, from
the unsettled layer. `docs/koskinen2022_model_a_comparison.tex` section 4.1
lists all ten differences with their sizes.

### 3.3 What is certified, and what limits the rest

The certification (`src/modules/time_step/certification.f90`) states one
condition per active balance, evaluated by one evaluator in four contexts
(probe, stationary trial, physical step, certification). The tolerances, READ
from the declarations:

| row | tolerance | how it was fixed |
|---|---|---|
| mass | 3e-12 | measured, contract section 9 |
| momentum | 1e-8 | measured, contract section 9 |
| energy | 1e-6 | A2tol2, a decade above the largest energy row of any state whose solve met its own target |
| closure, level balances | `ieq_res_tol` | the value the cell sweep already judges a root by |
| element, carrier | **1e-5 for r >= 1.20 R_p; below that reported and not gating** | decision 22 a (N30), anchored in `docs/certification_tolerance_anchoring_20260910.md` |
| time-discrete hydrodynamic row | 1.0, diagnostic | an order-one admissibility statement, not a convergence statement; it does not gate |

The species-row tolerance is a function of radius through
`cert_tol_element_at(r)` and `cert_tol_carrier_at(r)`, because an element or
carrier balance is a cancellation between transport terms whose ratio the
flow sets: in the wind the advective divergence carries the row, below the
homopause the eddy and settling terms carry it and each is orders of
magnitude larger than their difference. The same equation on the same state
reads 2 to 11 times smaller in the wind than in the layer (MEASURED, memo
section 1). Two reporting radii bound the windows, `cert_regime_layer_r` =
1.10 and `cert_regime_wind_r` = 1.20; the band 1.10 to 1.20 is reported and
does not gate, because it holds the candidates' binding cell at r = 1.153
where the operator's own discretization error is 4.6e-4 and no state of this
implementation has come near 1e-5 there.

[corrected 2026-09-11: the reason given here read "and a 1e-5 gate would ask
the residual to fall 1.5 decades below the error of its own discrete
equation". That inference is withdrawn. The algebraic residual of the
discrete equations and the truncation error of those equations against the
continuum are different quantities, and the first is not bounded below by the
second (review section 3.1 of `docs/solver_approach_analysis_20260910_review.md`);
an algebraic error well below the discretization error is what makes a
grid-convergence study readable in the first place. MEASURED, 2026-09-11: on
this same 500-cell grid the existing H2 transport routine took the wind H2
row from 7.380744e-2 to **2.449021e-13** on a held background
(`docs/solver_partition_experiment_20260911.md` section 4), nine decades
below the 4.6e-4 the withdrawn sentence treated as a floor on a species row
of this grid. Decision 22 (a) stands unchanged
(1e-5 gating for r >= 1.20, the layer and the 1.10 to 1.20 band reported and
not gating); only its justification changes, from "the discretization forbids
it" to "no state this implementation produces reaches it there, and the rows
below 1.20 are a cancellation of terms whose conservation is itself measured
at 2e-3 in the wind and 9.8 in the layer".]

**What certifies today:**

* **The atomic three-unknown stationary solve.** `wasp_full_newton`: info 0,
  `||R||` 4.3e-9, CERTIFIED, log Mdot 13.30 (LOGGED). It is in the regression
  matrix and its state is pinned as a reload fixture.

**What does not certify:**

* **The atomic element reload** (HD 209458 b, `He_diffusion`,
  `He_metal_diffusion`, eight element rows). Handed back at `||R||` 2.3e-4
  after N29; its binding row in the gated window is the helium row at 2.9e-4
  of its scale against 1e-5 (LOGGED, N30). NOT certified.
* **The molecular carrier reload** (the hot-Uranus column with H2 transported).
  Handed back at `||R||` 0.25; the H2 row of the binding cell reads 7.3e-2
  against 1e-5 (LOGGED, N30). NOT certified.
* **No other configuration with a species row** has converged on any route.

**What limits them**, in one sentence each, all MEASURED by the items named:
the linear solve, which cannot reach 0.1 of its right-hand side in 40 Krylov
products on this operator and whose Arnoldi image is 15 to 21 percent off it
(N24, N27, and N21 agreeing); the atomic reload's ulp-level chaos, a one-ulp
change of the upper ghosts moving `||R||` at iteration 40 from 3.2e-2 to 1.1,
because the Krylov leg's one-digit tolerance decides which vector crosses it
at the last bits of the residual (N26c); and, for the layer, the element
flux, conserved to 2.6e-2 in the wind and not at all below (memo section 9.8).
The consequence for acceptance, stated by N26c: **the reload's `||R||` at a fixed
iteration is not an acceptance quantity**; named outcomes are (binding rows,
refusals, flux spread). UPDATE 2026-09-10 evening (N31 to N36): the linear
algebra is exhausted as a lever. The Arnoldi gap is a ROUNDING FLOOR of the
residual, not curvature: the flux difference over a base cell amplified by
the near-hydrostatic cancellation at the face (N33; bound `epsilon x face
state x r^2/dV`), following no inner tolerance (N32); a quadruple-precision
evaluation of the hydrodynamic rows brings the image within 3.8e-3 of the
operator and the cycle still stalls (N34). The banded preconditioner misses
no coupling; the preconditioned operator is near-singular ON THE SPECIES
ROWS (Ritz ratio 6.6e5 carrier, 1.1e4 atomic; the binding carrier row's
diagonal 3.3e3 below its coupling to the hydrodynamic unknowns of its own
stencil) and a longer Krylov cycle returns a WORSE true step (N35); the
linear model's row scaling and a true-residual cycle are measured and both
fail acceptance (N36). N37 put the well-balanced flux difference in as the default-off key
`Well balanced:` (exact on the discrete equilibrium; the rounding floor does
not move; a prerequisite in `store_row_terms` before any default-on), and
N38 showed the binding row's small diagonal is a column-scale reading, not
physics. ~~The program is CLOSED (2026-09-10 evening, user)~~ **SUPERSEDED
2026-09-11** (see the dated paragraph after the bracket below); the analysis
of a different approach is `docs/solver_approach_analysis_20260910.md`.

[bounded 2026-09-11, review sections 3.4 and 4.1 to 4.3: both statements of
this paragraph are measurements of the present implementation and neither is
a property of the equations. The rounding floor is the floor of THIS flux
assembly in double precision, and N34 measured it moving (down by 2.1, the
Arnoldi gap by 31, from 1.19e-1 to 3.81e-3) under a quadruple evaluation of
the same rows, so it bounds no other formulation of the same physics. "Near
singular on the species rows" is a Ritz-magnitude ratio of the masked,
right-preconditioned action over 200 products, not a measured condition
number of the coupled system; what it establishes is that the directions
associated with the transported species are where the present preconditioned
Krylov solve makes no progress. Neither number says that a discrete stationary root does not
exist or that no method can certify these rows.]

**[reopened 2026-09-11]** The stage-2 solver program was reopened the next
day by the user's instruction, as the PARTITIONED route rather than the
coupled species-row solve, and it delivered: the plan is
`docs/PLAN_20260911_partitioned_solver.md` and the dated record with every
number is `docs/Update_EXHALE_stage2.md` **section 8**. Three defects of the
partitioned route were found and repaired (a failed element composition
solve handed back as an update; a carrier movement bound tested after the
step was applied; a stop test weaker than the acceptance test at return),
the outer iteration became a stated contract accepting a state only by the
certification of the FULL set of active equations on the refreshed state,
the momentum reference scale under `Well balanced: True` was corrected (with
the key on `wasp_full_newton` now CERTIFIES, `info = 0`, `||R||` 9.719e-10,
where its momentum row stood at exactly 1.000), and the pseudo-time start of
`Restart intent: stationary` was set to 1.0 instead of the CFL interval.
MEASURED outcome on the two fixtures, 8 threads, 500 cells: the partitioned
route reaches a strictly lower joint residual than the coupled solve in no
more wall time on both, which is the plan's adoption rule, and it is what
the production outer loop runs (`Coupled carrier solve` already defaults to
False). Neither fixture certifies. The HD 209458 b element reload refuses on
the hydrodynamic mass row ALONE, 1.054e-9 against 3.0e-12 at cell 16, with
every elemental wind row inside 1e-5 after twelve passes; the hot-Uranus
carrier reload refuses on an H2 wind row still moving, 2.234e-2 against
1e-5, its binding cell walking 290 to 443 over 40 passes as the front
advances about one cell in a pass. So what this paragraph's "no
configuration with a species row converges" reports is now bounded further:
the alternated species rows of the atomic configuration DO come inside
their tolerances, and what refuses that state is the base-layer mass row of
the three-unknown hydrodynamic solve, the N33 rounding floor. **Closed the
same afternoon (user decision, item P16, `certification_tolerance_anchoring_20260910.md`
anchor (6)):** the mass row's tolerance is now the larger of 3e-12 and ten
times the cell's own estimated rounding floor, and the HD 209458 b element
reload CERTIFIES at pass 12 (MEASURED, P16 and P17). The hot-Uranus carrier
reload remains uncertified on its H2 wind row: a 110-pass continuation
(decision 2) shows the front still moving (x2 = 0.5 from 1.0726 to
1.1230 R_p) and a uniform wind floor of 2.2e-2 that the bounded carrier
pass does not relax (Update_EXHALE_stage2.md section 8, "User decisions"). Item P21,
which made the carrier relaxation solve the chemistry of the composition it
advances, then measured what that front is: at a fixed wind the joint fixed
point of transport plus chemistry is a FULLY MOLECULAR column (x2 = 0.98 at
1.1 R_p and 0.33 at 3.0 R_p, the front outside the domain, no cell clamped onto
the element budget), so the front's position is set by the hydrodynamic
response and not by the carrier operator alone, while the 2.2e-2 wind floor is
reproduced digit for digit with the chemistry refreshed and is therefore the
bound's measure, absolute in the grid maximum, and not the frozen chemistry
(MEASURED, P21; Update_EXHALE_stage2.md section 8, "P21: the fixed point at a fixed
wind, and the production loop").

## 4. Where the problems begin (ordered by size and depth)

The full account of each, with its evidence and its next item, is
`docs/ISSUES_20260909.md`. This is the ordering.

1. **No stationary solve with a species row converges** (ISSUES 3.1). Four
   distinct obstructions were named and three removed: the carrier's bounds
   (resolved by the `ln n` unknown, decision 11 a, and by making the element
   budget a face of the box, B5l); the shared element constraint, which the
   coupled-carrier Newton violated by 4 to 43 percent in He/H (resolved by
   N4b, route i of decision 14: the write-back returns the element totals it
   was handed, 1.4e-11 over 68 iterates); the trace-element rows of absent
   elements, which held the trust radius at zero (N1); and the outer ghost of
   a transported column, frozen at the pre-solve composition (N26). What is
   left is the linear solve and the model behind it (section 3.3). Blocks: a
   converged Model A comparison, a converged hot-Uranus Mdot, any molecular
   case that is not a fixed-step snapshot, and any element-transport state of
   HD 209458 b, which has never been obtained on any route.
2. **The advection post-process on a state that is not stationary** (ISSUES
   3.2). No matrix case has a stationary mass flux by the cell-to-cell
   measure: `rho v r^2` over outflowing cells spans 1.27x its median in
   `wasp_full`, 2.2x in `mol_base_handoff` and 83x in `lower_profile`. The
   product is now a conditional correction that says row by row which it is,
   so this is a fact about the runs rather than a defect of the tool, but
   **every `_adv` product and transit spectrum in the repository predates the
   N11/N11b/N12 definition** and the molecular ones also predate the N12
   physics fix (T_adv up to 94 and 103 percent higher in the molecular
   layer).
3. **Local-equilibrium ions in a transported wind** (He+, the molecular ions;
   for the atomic line every stage): the same error H+ had before the carrier
   route existed, fixable with the existing operator.
4. **The base cell amplifies any composition inaccuracy by about 4e3**
   (ISSUES 3.3), so any acceptance measure taken as a maximum over cells
   reports cell 1 on a molecular configuration. Whether the base cell's
   hydrodynamic rows should be judged on the interior's tolerance is an
   acceptance question not yet put to the user.
5. **The O(dt) kick on a reloaded stationary state** (ISSUES 3.4): a state
   reported stationary to `||R||` 8.4e-9 read 1.226 after the two CFL steps
   the marching loop used to take before the Newton. `Restart intent:
   stationary` now enters the JFNK with no CFL step (N10), which is the
   hand-off that item asked for; re-measuring the kick itself at `dt`,
   `dt/2`, `dt/4` from a seed-independent evaluator is still open.
6. **Heating per ionization, external** (Koskinen Model A). Their heating
   profile is not what their stated accounting and spectrum produce. EXHALE
   keeps h nu - I; the comparison is a bracket until the authors' spectrum
   and heating implementation are known. Do not tune the fraction silently.
7. **The base limit cycle / breathing** at the default CFL: numerics, not
   physics; it spoils fixed-step molecular snapshots and `du` convergence.
8. **The ionization front's second root** (ISSUES 3.1 c). With H2 transported
   the H2 front has one root, but the H I partition at the front has two on
   the element + carrier configuration (9 seeds and 1 seed, 1.05e-1 apart).
   Which root a stationary residual may land on, or a statement that the
   front is unresolved on this grid, is a physics rule the code does not
   have. Study item, no owner.
9. **The physical step is never refused for its integration error, for its
   species** (ISSUES 3.8): `err_rtol` and `err_atol` are flat and unvalidated
   and cover no transported species, and a retried macrostep's accepted state
   depends on which pass was rejected.
10. **H3+ non-LTE cooling** differs from the published detailed-balance
    treatment by about 3x per molecule at wind densities.
11. Smaller, all in ISSUES 3.6 and 3.7: `Load IC` cases written before the
    Jupiter-radius unification cannot be reloaded (the grid guard refuses
    them); the spectrum file has nothing below 5 A; four independent
    `# columns` parsers in python; `carrier_element_headroom` states no joint
    bound for carriers sharing an element; the epsilon-active band of the
    projection never engages on the present fixtures.

## 5. Tasks, in priority order, with acceptance criteria

1. ~~The column preconditioner of the species-row linear solve~~ DONE and
   REFUTED as a lever (N31 to N36, 2026-09-10): probe rule, inner
   tolerances, extended precision, band coverage, row scaling of the linear
   model and cycle length are all measured; none frees the species-row
   linear solve. Replaced by: **decision 23** (a well-balanced flux
   difference in the near-hydrostatic base layer, the physically right
   discretization; not by itself the cure) and **the discretization of the
   species rows that bind** (the carrier row of the front cell whose own
   diagonal is 3.3e3 below its hydrodynamic coupling; the element rows of
   the front cells): a statement of the discretized element and carrier
   equations at the front, their Peclet number and what makes the row
   nearly independent of its own unknown. DONE by N38 (the term is the
   advective coupling to the momentum column, scaled by the sound speed);
   ~~the program is CLOSED by the user~~ **SUPERSEDED 2026-09-11**: it was
   reopened as the partitioned route and the discriminating measurement was
   made (`docs/PLAN_20260911_partitioned_solver.md`,
   `docs/Update_EXHALE_stage2.md` section 8). The corrected form of that section 7
   is in the document itself; the analysis document's own errors were
   corrected the same day.
2. **The layer's element flux conservation** (memo section 9.8): conserved to
   2.6e-2 in the wind and not at all in the layer (39). Acceptance: a stated
   conservation measure of the element operator in the layer, and, if it does
   not improve, a written statement of what the layer's tolerance can be.
3. **The element operator's discretization** (order 1.59 on a manufactured
   column; about 5e-4 at the binding radius at N = 500, and a remapping
   estimate that puts about 78 percent of what the atomic candidate stands at
   in that direction). Acceptance: the measured order at N = 500 and
   N = 1000, on independently relaxed grids or against a manufactured
   solution with known forcing, and a statement of the ACCURACY that order
   buys at the gate radius. [corrected 2026-09-11, review sections 3.1 and
   3.2: this is an accuracy question, not a bound on the algebraic residual
   the solve can reach, and the 78 percent is a Richardson-style estimate on
   a remapped, unconverged candidate, so it is a consistency diagnostic and
   not a split of the residual into an unavoidable and an avoidable part.
   The acceptance no longer asks "is the gate above the discretization"; a
   frozen-background H2 wind row of 2.449021e-13 was MEASURED on the same
   500-cell grid, `docs/solver_partition_experiment_20260911.md` section 4.]
4. **Transport He+ (then H2+, H3+, HeH+) as carriers**; `ic_Hp` is the
   template. Acceptance: He+ on the Koskinen gate within the plotted points;
   element census closed to 1e-10; key-off bitwise.
5. **N16b's remainder**: validate `err_rtol` and `err_atol`, extend the
   estimate to the transported species, and close the retried-macrostep path
   dependence (some state the coupled source step reads is neither saved nor
   rebuilt by the checkpoint). Acceptance: a physical-mode trajectory that
   can be quoted as error-controlled for its species.
6. **Ask the authors** (user decision) for the Koskinen mean solar spectrum
   file and the heating deposited per ionization, or Figure 9 in numbers.
7. **H3+ non-LTE factor**: the Koskinen et al. (2009) Table 2 ratios as an
   alternative table. Acceptance: cooling of each molecule at 3 r_base within
   1.5x of what their Figures 8 and 9 imply.
8. **The base limit cycle**: local time-step control in the base cells (CFL
   below 0.15-0.3 there). Acceptance: `mol_sec_ion` base velocity free of the
   odd-even pattern; goldens refreshed once.
9. **Instruction-gated (never autonomous)**: the D-physical reruns of the
   planet folders and `benchmarks/` on the loaded SEDs with the new physics,
   Newton-finished; `Load IC` case regeneration on the current grid; a
   converged run reaching 3000 K for the CO domain; paper and poster
   re-convergence; the LHS 1140 b write-up; the GitHub public switch.

## 6. Conventions and pitfalls worth remembering

* A comparison run states every condition and is judged on plotted profiles,
  not on two numbers.
* Physics defects are fixed when found; goldens follow (refresh, report the
  movement); technical movement below 0.1 percent needs no refresh.
* Goldens are never refreshed before a fix series. The control during a
  series is a scratch baseline (`REGRESSION_GOLDEN_DIR`); the refresh happens
  once, at the end, and is reported.
* Comparison-only options (`Photoelectron heating: full|<f>`, `Caloric EOS:
  monatomic`, `Atomic rate set`) are labeled so in the setup report and in
  `EXHALE_resolved.out`; a run with one must say so in any report.
* A stated `Stellar LW flux: 0.0` turns the band off even with a spectrum
  file; an absent key with a spectrum file integrates the band.
* An exact `options` comparison guards every restart, so the old restart-ladder
  workflow (converge without an option, restart with it on) needs
  `Restart option change:` naming the tokens allowed to differ.
* Kill by PID after checking `/proc/<pid>/cwd`; never `pkill`. Use `\cp -f`
  and `\mv -f`, and `tail -n N`; `$PWD` inside `( cd dir && ... )` is `dir`.
* Workers build privately (`make OBJDIR=build_<item> EXE=EXHALE_<item>.x`)
  and never run the shared `make check` while another worker is in the tree;
  `run_check.sh` takes a lock and refuses a second run.
* The auto-mode classifier refuses multi-case golden refreshes and compound
  commands that include one.

## 7. Where things are

* Problems, resolved and open: `docs/ISSUES_20260909.md`; its review
  `docs/ISSUES_20260909_review.md`; the plan `docs/PLAN_20260909_rev1.md`;
  the common rules every worker follows, `docs/worker_rules.md`.
* Decisions 1 to 22 with their options and what was chosen:
  `docs/To_be_determined_by_user_20260906.md`.
* The species-row tolerance anchoring, with the five measured anchors:
  `docs/certification_tolerance_anchoring_20260910.md`.
* The certification and run-mode contracts:
  `docs/a2_certification_contract_20260906.md`,
  `docs/a0_run_mode_contract_20260906.md`.
* The stationary solver's design, section by section, with section 22 the
  stage-2 outcome: `docs/steady_solver_design.md`.
* Koskinen comparison runs: `benchmarks/koskinen2022_model_a/`; figures and
  table `docs/k22_model_a_figures.py`; digitized published profiles
  `docs/p23_published_profiles.md`.
* Reload fixtures: `backup/regression/wasp_full_newton/IC/` (certified
  atomic) and `backup/regression/atomic_elem_newton/` (eight element rows).
* Open-items list: `TO_BE_DONE.md`.
* Memory index:
  `~/.claude/projects/-nfs-mocafe-kiseon-RT-Codes-ExoAtmosphere/memory/MEMORY.md`.
