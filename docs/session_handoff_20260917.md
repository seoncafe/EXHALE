# Session handoff, 2026-09-17: the 2026-09-16/17 series on `docs/PLAN_20260916_rev3.md`

Written 2026-09-17 at the user's request, at the close of the series. The
previous handoff is `docs/session_handoff_20260915_rev1.md` (its sections 1 to
21 carry the campaign and the items that preceded this one); the ones before it
are `docs/session_handoff_20260912.md` and earlier.

**Every number in this document is READ from the memo, run log or plan row
named beside it. Nothing was re-measured for this document**, with three
exceptions, each labeled MEASURED where it appears: the regenerated tables of
`LHS1140b/MODELS.md` (section 3), the page count of
`docs/lhs1140b_exhale_vs_pwinds.pdf` (section 7) and the file lists of section 2
(read from the memos and checked against the modification times of the tree).

**The final binary.** `EXHALE.x` md5 `59bfdb3fc4d0104fc2e9c3734596d2f6`,
built 2026-09-17 10:14 with bare `make` (conda-forge gfortran 16.2.0, the
OpenBLAS of its prefix), reproducible from `make distclean`; manifest
`LHS1140b/models/BINARY_MANIFEST_59bfdb3fc4d0.txt`. The golden table is in
section 2.3 and in `docs/Update_EXHALE_stage2.md` section 11 ("The golden
refresh of 2026-09-17").

Until that block is filled, the binary of record in this document is
`c2e9c9990b9f14f1be8cd77abca68945` (the catalog binary of 2026-09-16 11:05),
which is what every catalog state on disk was made on and what the memos of
this series call "the tree binary".

---

## 1. What the series did and what it found

Eleven items of `docs/PLAN_20260916_rev3.md` were executed between 2026-09-16
17:37 and 2026-09-17 08:00. In one page:

**L26, the base boundary was discontinuous, and the repair changes no
physics.** `wind_window_mass_flux` returned `have_F = .false.` at exactly zero
mean window flux, and the blend then fell back to the cell-1 weight, so the
boundary density jumped at a zero window; the jump is directional (a uniform
window and a nonuniform one reach the same zero point with two different
limits) and measures 37.407 and 37.468 per cent of `rho_b` on the two LHS
1140 b wind states and 0.411 per cent on the cold start, equal to
`|w_i - 1/2| |rho_res - rho_rev|` on every state; the `max`/`min` spread had a
converged slope jump of 242 at the argmax switch (L26 memo sections 1, 4, 8).
The repair makes the window's say `s_wind = A(M_wind) C(d_window)` with
`d_window` the relative standard deviation (no argmax, so no slope jump;
threshold recalibrated to `2e-3`, the geometric middle of the measured gap
between 302 certified states at 4.06e-05 to 1.82e-04 and the tightest non-wind
state at 1.39e-02) and `A` a cubic smoothstep in the window Mach number that
vanishes quadratically at zero flux. Three purely local closures were measured
first and all three are refused by the states: every reading of the base mass
flux built from cells 1 to k is negative on the four LHS wind states whose
branch is inflow (L26 memo section R3). **Nothing the code meets moves**:
`rho_b`, both ghosts and `w_rev` are bit for bit the entry text's on five
probed states; `wasp_full`, `wasp_he23off`, `mol_base_handoff` and
`lower_profile` reproduce every data line over 12000+ steps, the two provenance
lines excepted; the fiducial re-solved returns `info = 0`, CERTIFIED, and
reproduces the catalog profile to 1.83e-16 in density and 2.05e-12 in velocity
(L26 memo section R7, R6.1).

**L25, the low-XUV stalls were the HLLC flux's resolution failure on the
500-cell grid.** The unmodified HLLC flux stands at `||R|| = 1.779` after 60
Newton iterations on `x0.10_kzz1e9/HeH9.7` at 500 cells; the same unmodified
HLLC flux CERTIFIES the same case on a 1000-cell grid (pass 11, `||R||`
3.552e-07), and that state is the 500-cell Roe state (Mdot +1.03 per cent, EW
-0.40 per cent, T and rho within 0.1 to 3 per cent above 1.005 R_p), so the
stall is a resolution failure of that flux and not a different wind (L25 memo
section 3.1, 3.5). With `Numerical flux: ROE` as the only change three of the
seven unsolved low-XUV cases certify directly, and the four remaining certify
by continuation in XUV, each rung seeded from the certified rung above it: all
eighteen rungs certified and not one ladder was stopped (L25 memo section 4.5).
**The continuation does not change the wind**: the ladder's target state and
step 3's uncertified direct solve of the same case agree to 4.4e-07 in T, rho
and v on all 500 cells, and what moved is the gated elemental He/H partition
row, 2.234e-04 refused against 7.20e-06 accepted, so what those cases needed
was decades of elemental relaxation inside the pass budget and not a flux, a
ramp or a source change (L25 memo section 4). **The flux-family systematic**
where both fluxes solve: Mdot -0.31 and -0.62 per cent, He I 10830 EW -0.31 and
-0.44 per cent, first-cell T -5.70 and -7.21 per cent and first-cell rho +5.86
and +7.52 per cent, the mass-flux spread unchanged to three digits (L25 memo
section 3.4). The published Rieper (2011) `min(Ma,1)` velocity-jump factor was
implemented behind `Low Mach velocity jump:` (default False) and MEASURED TO
PREVENT the solve: it is exactly zero where two neighboring velocities cancel
in the Roe average, so it removes the damping of the odd-even mode at exactly
the faces where that mode lives, and Rieper's premise (local and global Mach
comparable) is violated by eight decades on this column (L25 memo section 2).

**L7g, five physically wrong statements of the molecular layer, corrected from
published data.** (i) R5 (`H2+ + e`) and R16 (`HeH+ + e`) deposited their whole
enthalpy as gas kinetic energy; one H atom leaves in n = 2, so the gas gets
0.749 and 1.554 eV instead of 10.948 and 11.753 eV and 10.199 eV goes to the
H(n = 2) population (Takagi 2002; Giusti-Suzor et al. 1983; Guberman 1994).
(ii) R6 (`H3+ + e`) deposited its whole 9.250 eV promptly; 2.4795 eV is
internal energy of the H2 fragment and branches through the quench fraction
(Kokoouline et al. 2001; Strasser et al. 2001). (iii) R12 and R15, the
three-body association and its reverse, applied the M = H2 rate coefficient to
the TOTAL heavy-particle density; the collider sum is now formed with the
published per-collider coefficients, helium given argon's as an upper bound
(Cohen & Westberg 1983). (iv) the quench fraction's radiative rate used the
v = 1 total decay rate 8.3e-07 s^-1 while the code's own validity claim
asserted all levels; the all-level maximum over the 302-level ladder is
5.5939e-06 s^-1, a factor 6.56, cross-checked against Wolniewicz et al. (1998)
to 1.1e-04. (v) the collisional side was two 1979 fits and helium was absent;
three published quantum calculations replace them, one per collider (Lique
2015 for H, Jozwiak et al. 2024 for He, Le Bourlot et al. 1999 for H2), the old
H2-H2 fit being 21 to 500 times the measured rate (L7g model memo section 0, 1;
L7h addendum). **The base heat is down 22 per cent**: `heat_mol_chem` at cell 1
of the certified `molecular_scalar_gj1132_kzz1e9/HeH2.13` falls 21.6 per cent
(1.11205700e-06 to 8.71940708e-07 erg cm^-3 s^-1) and 9.4 per cent at cell 25;
total heating falls 21.3 per cent at cell 1; about 95 per cent of the movement
is (iii) (L7g model memo section 3.1). **The molecular catalog no longer
roots**: that certified state is not a solution of the corrected equations, its
energy residual sitting at 4.27e-01 through two complete outer passes with
every other row inside tolerance, and L7h adds that its H2 carrier balance row
moves from 5.645e-07 to 7.817e-01 once the collider-resolved rate reaches the
composition solver too. Re-solving it is a campaign and not a bounded pass (L7g
model memo section 3.5).

**L22, the movement bound is not the defect.** Relocating the bound-controlling
cells only moves the bound's control to the mask edge, and the refusing cell's
row measure rises in every waived run (`wellmixed/HeH0.083` cell 500: 2.462e-02
control to 3.345e-02, against 4.537e-02 narrow waive and 2.362e-01 far waive;
L22 memo step 2). The two named competing explanations were tested by
intervention and both excluded: stating an outer diffusive face makes the row
worse by factors 1.12 to 182, because H2 is heavier than the mean particle of
an ionized hydrogen wind and the settling drift such a ghost admits is INWARD,
thirty times the row's largest term; and the deferred reconstruction terms of
the carrier Jacobian reproduce the central difference to a mean relative error
1.542e-03 and the pass run with them in is the control pass to every digit (L22
memo step 2b). **The alternation fails where the composition update leaves the
wind's linear range**: the one run whose carrier block reached its own fixed
point (residual 7.32e-10 against 7.01e-03 in the control) moved the base
particle count by 6.893e-02, 6.9 times the bound, and the hydrodynamic solve
that had to follow it plateaued in a limit cycle, `||R||` between 1.100 and
1.162 from iteration 100 to 205 with all 500 mass rows outside tolerance (L22
memo step 2b (C)). **The coupled block as designed fails on the Krylov solve**:
on the certified `kzz1e9/HeH2.13` state the alternation CERTIFIED at outer pass
12 (worst gated carrier row 7.82e-01 to 8.18e-06) while the block over a longer
wall clock stood at trust-region iteration 14 with `||R||` moved from 4.267e-01
to 4.241e-01, 0.6 per cent, and every one of its thirteen linear solves ended
with the subspace exhausted (40 products of 40, relative residual 2.597e-01
against the 1.00e-01 asked for) (L22 memo step 3, I3). **The ghost rule was
unified**: both ghost rules of the carrier column now live in
`carrier_face_mass_fraction`, so the certification row and the relaxation row
are one function of the interior state; the rule is `X_ghost = X_N` in both
directions of the face mass flux, and no diffusive, eddy or settling flux
crosses either end of the column. Test row at tolerance zero: 2.1724e-02
before, exactly 0.0 after; the three frozen states are bitwise unchanged and
`mol_carrier`, `mol_base_handoff` and `mol_diffusion` are byte-identical over
12000 steps but for the provenance line (L22 memo step 2c).

**L9, the old post-processing pass wrote PLM snapshots over solved states.**
The runner's `Do only PP: True` with `CFL: 1.0e-12` was found to be doing three
unrequested things: reconstructing with PLM where the solution was reached with
WENO3; overwriting the solved state's own `Hydro_ioniz.txt` with an uncertified
relaxation snapshot (`certified=F cert_reason=no_stationary_claim`), so every
`LHS1140b/models/*/output/Hydro_ioniz.txt` on disk today is such a snapshot with
the certified state surviving only in the `_IC` copy; and writing a `# coupling:`
header, which implies a certification claim, onto the `_adv` files, which
describe a different composition. `Restart intent: stationary evaluate` now
names three states and writes every product from the work state with no step
taken. The measure the operator itself weights by, `adv_mass_row`, falls from
1.5731e+00 to 1.7338e-09 at cell 1 and from 3.0840e-05 to 3.9957e-14 at cell
500, nine decades at the base and four in the wind; the He I 10830 red-pair EW
moves by 5.1e-06 relative (L9 memo sections 5, 6).

**L23, the certification pair of a restart header.** `(certified,
cert_reason)` is now imported-state metadata carried unchanged by a reload that
writes without measuring and overwritten by any run that evaluates the state;
an over-length reason, a duplicate key and a pair whose halves disagree refuse
the load. A defect found on the way: the coupling line was selected by a
substring test that also matched `# mapped-from-coupling:`, so eight kept seeds
under `.L4h/`, `.L14/` and `.L18/probe_seedmap/` were loaded as certified when
their own line said `certified=F` (L23 memo sections 4, 6a, 6b).

**L15 closed with no thread dependence.** The whole 105-iteration trajectory of
the element reload, 152 Krylov products and 923 recorded quantities, is bitwise
identical at 1, 8 and 16 threads, one md5 `bd12e7c8ad917beda538d789080e10e3`
across all traces; the two molecular configurations (the certified fiducial and
the stalled well-mixed state the campaign spent hours on) likewise; the closing
measurement on the catalog binary reproduces both. The 2026-09-15 divergence is
not reproducible on the present tree and is attributed to nothing. This is
non-reproduction, not a proof that no race exists (L15 memo, stages 2, 5, 6).

**L12b, the derivation.** The L12a design's stage flux double-counted the eddy
term; the corrected face flux is
`F_k(f) = x_k(f) N_el(f) - n_el(f) K(f) [x_k(j+1) - x_k(j)] / dr(f)`, and
`sum_k F_k(f) = N_el(f)` holds to rounding under three stated conditions
(one shared element flux, face fractions closing the simplex exactly, and the
closing stage's eddy term equal to minus the sum of the others'). Fifteen rows
in the new suite `src/tests/ionization_stage_flux/`, 15 passed, with four
deliberately mismatched rows confirming the test discriminates (L12b memo).

**L24, the `atomic_elem_newton` fixture.** Neither entry pseudo-time certifies.
The default entry descends `||R||` 1.674e+00 to 4.892e-05 and stops on the
stagnation detector at iteration 90 with 4 of 11 entries still refusing; the
`1e8` entry does not descend at all (stagnation at iteration 64, residual
1.903, 11 entries refusing), so the L14 recipe is EXCLUDED as the repair. What
holds the default-entry run is the linear solve, 62 of 91 Krylov cycles
exhausting the 40-vector subspace at 0.91 to 0.999 of the right-hand side on
the helium element row of cell 335; admissibility and every state-dependent
residual branch are excluded (0 refusals, 0 activations). The fixture is not
re-pinned (L24 memo).

**The renaming.** Sixteen `backup/regression/` case directories carrying the
banned noun were renamed (twelve active, four quarantined). None appears in
`golden/` or in any dated `golden_*/` snapshot and none is in `DEFAULT_CASES`,
which is why the rename provably could not move a number; only directory names
and `NOTE.txt` prose were edited, with every citing document updated
(`docs/named_case_audit.md` section 6 carries the mapping).

---

## 2. The state of the tree

### 2.1 Source files changed in this series, and by which item

READ from the memos, checked against the modification times of `src/`. The
series began at 2026-09-16 17:37; anything with an earlier time belongs to the
previous session.

| file | item |
|---|---|
| `src/modules/states/base_boundary.f90` | L26 repair (-97 / +183) |
| `src/modules/flux/Num_Fluxes.f90` | L25 step 2 (the `Low Mach velocity jump:` branch) |
| `src/modules/files_IO/input_read.f90` | L25 step 2 (the new key), L22 step 3 I2 (`Coupled carrier solve:` three-valued) |
| `src/modules/files_IO/write_setup_report.f90` | L25 step 2, L22 step 3 I2 (the echoes) |
| `src/modules/files_IO/load_IC.f90` | L23 (+215 / -18, the certification claim), L22 step 3 I3 (`opt_is_route`) |
| `src/modules/files_IO/write_output.f90` | L9 (+40, `write_derived_provenance_header`), L7h (the R12/R15 diagnostic and the FUV band ledger) |
| `src/EXHALE_main.f90` | L9 (the evaluate route, about +260), L23 (+9 / -1), L22 step 1 (the four progress measures) and step 3 I2 (the handover) |
| `src/modules/time_step/steady_newton.f90` | L15 (+339, instrumentation, off by default), L22 step 3 I1 (the coupled Jacobian probe) |
| `src/modules/lower_atmosphere/diffusive_photochemistry.f90` | L22 steps 2, 2b, 2c (the mask, the interventions, the unified ghost rule) |
| `src/modules/functions/binary_element_diffusion.f90` | L22 step 1 (`element_transport_residual_norm`) |
| `src/modules/lower_atmosphere/molecular_reaction_heat.f90` | L7g model (+303 / -35) |
| `src/modules/lower_atmosphere/mol_rates.f90` | L7g model (+535 / -50, the collider-resolved association) |
| `src/modules/lower_atmosphere/h2_vibrational_relaxation.f90` | L7g model (+464 / -61), L7h (`gamma_10_H2` from the para-H2 file) |
| `src/modules/lower_atmosphere/lyman_werner.f90` | L7g model (comment only, +11 / -1) |
| `src/modules/radiation/excited_hydrogen.f90` | L7g model (+83 / -11, the chemical n = 2 source) |
| `src/modules/radiation/ionization_equilibrium.f90` | L7h (+271 / -6, He ground singlet into `f_vib_quench`) |
| `src/modules/radiation/util_ion_eq.f90` | L7h (+43 / -12, the `nheiS` argument) |
| `src/modules/nonlinear_system_solver/System_HeH_mol.f90` | L7h (+210 / -51, `mk15` from the collider sum in the composition rows) |

Test files added or changed: `src/tests/low_mach_stress_energy/` (new, L25 step
1), `src/tests/physics_probe/roe_low_mach_velocity_jump.f90` (new, 24
assertions, L25 step 2), `src/tests/physics_probe/molecular_energy_recipients.f90`
(new, 424 lines, 62 assertions, L7g), `src/tests/physics_probe/molecular_third_body_in_the_chemistry.f90`
(new, 252 lines, 9 assertions, L7h), `src/tests/physics_probe/run.sh`,
`src/tests/ionization_stage_flux/` (new, 15 rows, L12b),
`src/tests/carrier_outer_boundary/` (new in L22 step 2b, rewritten to 16 rows in
step 2c), `src/tests/coupled_block_jacobian/` (new, 8 rows, L22 step 3 I1),
`src/tests/grid_and_gates/base_boundary_continuity_probe.f90` (probe in L26,
grown to a 32-row suite in the repair),
`src/tests/grid_and_gates/base_branch_discriminant.f90` (L26, 6 rows),
`src/tests/grid_and_gates/stationary_evaluate_products.sh` (new, 20 assertions,
L9), `src/tests/grid_and_gates/restart_intent_and_metadata.sh` (L23, sixteen new
rows), `src/tests/grid_and_gates/restart_option_change.sh` (L22 I3, five new
rows), `src/tests/grid_and_gates/coupled_carrier_h2_row.sh` (L22 I2, five new
rows), `src/tests/grid_and_gates/run.sh`,
`src/tests/element_operator/element_operator_tests.f90` (L22 step 1, 6 rows),
`src/tests/carrier_retry/carrier_retry.f90` (L22 step 1, 9 rows),
`src/tests/steady_species_rows/steady_species_rows_tests.f90` (L22 step 2c).

Catalog scripts: `LHS1140b/models/run_case.sh` (L9, the post-processing pass on
the evaluate route), `make_models.py`, `write_reproduce.py`, `status.py` and
`LHS1140b/MODELS.md` (the L25 decision: the flux column and the frozen
seven-case ROE list), and four new scalings
`LHS1140b/sed/lhs1140_sed_gj1132_at_b_xuv0p{18,16,28,26}.txt` written by the
rule the existing grid was written by and checked by regenerating five existing
files to all 29992 rows each (L25 step 4).

### 2.2 New keys and their defaults

| key | default | item |
|---|---|---|
| `Low Mach velocity jump:` (input) | `False` | L25 step 2; byte-identical off, and measured to PREVENT the solve when on |
| `Coupled carrier solve:` (input) | `False` (unchanged); now three-valued `False` / `True` / `On stall`, and a word it has no meaning for STOPS the run | L22 step 3 I2 |
| `EXHALE_REACTION_HEAT_RECIPIENTS` | unset = corrected physics; `0` restores the old recipients and the total-density third body | L7g model |
| `EXHALE_BASE_WIND_AMPLITUDE` | `1e-12` (the code constant `base_wind_window_amplitude`) | L26 repair (new) |
| `EXHALE_BASE_WIND_SPREAD` | recalibrated `1e-2` to `2e-3` | L26 repair |
| `EXHALE_BASE_WIND_SAY_ON_MACH` | REMOVED | L26 repair |
| `EXHALE_EVAL_STATE_DUMP` | off | L9 |
| `EXHALE_ELEMENT_DRIFT_IS_MAP_DISTANCE` | off | L22 step 1 (restores the superseded element measure) |
| `EXHALE_L22_MASK`, `EXHALE_L22_MASK_MODE` | off | L22 step 2 (diagnostic) |
| `EXHALE_L22B_JAC_RECON`, `_JAC_ACTION`, `_JAC_CELLS`, `_DISPLACEMENT` | off | L22 step 2b (measurements, kept) |
| `EXHALE_L22B_OUTER_DIFF`, `EXHALE_L22B_OUTER_ADV` | REMOVED | L22 step 2c |
| `EXHALE_COUPLED_JAC_ACTION`, `_CELLS`, `_DIR` | off | L22 step 3 I1 (writes and stops, takes no step) |
| `EXHALE_L15_TRACE`, `_DUMP_AT`, `_REPEAT` | off | L15 stage 2 |

No existing default changed, and no golden was refreshed by any item of this
series.

### 2.3 The binary and the goldens

`EXHALE.x` md5 `59bfdb3fc4d0104fc2e9c3734596d2f6` (predecessor, the catalog
binary, `c2e9c9990b9f14f1be8cd77abca68945`); manifest
`LHS1140b/models/BINARY_MANIFEST_59bfdb3fc4d0.txt` (the md5 of every source
file, the items carried beyond the predecessor, the new keys and their
defaults). MEASURED by the closing worker on 2026-09-17, single-threaded
(`make check` on the 16 default cases plus `roundtrip`, tolerance 1e-3):

| case | against the 2026-09-16 goldens | what moved |
|---|---|---|
| `wasp_full`, `wasp_he23off`, `wasp_full_newton`, `lower_profile`, `hydrostatic_column` | byte-identical data lines | nothing: no molecular network runs in them |
| `mol_base_handoff`, `mol_metals`, `mol_lyman_werner`, `mol_diffusion`, `mol_ir_bands`, `mol_sec_ion`, `mol_carrier` | beyond tolerance | the H2 front moved outward: n(H2) at r = 1.08 to 1.10 R_p 5 to 190 times the golden, H3+ following; rho and T by 8 to 20 percent at r = 4.1 to 4.3 R_p (item L7g and its follow-up) |
| `oxygen_chemistry` | beyond tolerance, order unity | the same physics on a 1000-step snapshot with oxygen carriers; not apportionable |
| `hp_zero_seed`, `hp_trace_seed`, `hp_front` | beyond tolerance by 4e-3 to 2e-2 | the same physics, small because the proton transport cases sit mostly outside the front |
| `roundtrip` | beyond tolerance in the data (the same hot-Uranus gate) | L7g; and its second half, the coupling header round trip, now PASSES (`cert_reason` preserved, item L23), which failed on the predecessor |

Attribution: every case that runs the molecular network moved and no other
case moved, the gate item L7g asserted. Nothing else moved a case: L9 changes
an `_adv` header only, L22 keeps the default route, L25's and L15's keys
default off, L26 leaves the five non-molecular cases byte-identical.

The goldens were refreshed ONCE: the 2026-09-16 set is kept entire as
`backup/regression/golden_L7g_20260917/` with a `NOTE.txt`; after the refresh
`make check` and the roundtrip check read byte-identical, 17 of 17
directories. `backup/regression/run_fcheck.sh` (bounds-checked build, one HD
209458 b case): clean, no trap. No pinned fixture was refreshed.

Everything on the catalog disk today was made on `c2e9c9990b9f14f1be8cd77abca68945`,
which predates every source change of this series.

---

## 3. The catalog

### 3.1 What stands

95 cases in `LHS1140b/models/`, solved on `c2e9c9990b9f` with the OLD
post-processing products. MEASURED 2026-09-17 from the tables
`models/status.py --write` regenerated: **82 rows certified** (79 atomic
prescribed and closure rungs, plus the three molecular cases
`molecular_scalar_gj1132_kzz1e9/HeH0.55`, `HeH2.13` and `HeH9.7`) and 13 not.
Of the 13, seven are the low-XUV atomic cases, which carry
`Numerical flux: ROE` in their `input.inp` and a `seed_from_ladder.txt` but no
`output/` state yet, their certified states living in `.L25/` and
`.L25/ladder/`; the other six are molecular and are listed in section 3.4.
One caution when reading the tables: `molecular_scalar_gj1132_wellmixed/HeH2.13`
is printed as `running`, which is what the files on disk say; the case was in
fact never started, as its own `not_solved.md` states in the reason column
beside it.

**Every `LHS1140b/models/*/output/Hydro_ioniz.txt` on disk is an uncertified
relaxation snapshot** written by the old post-processing pass, with the
certified state surviving only in the `_IC` copy (L9 memo section 1). That is
true of all 95 cases and is one of the reasons for the re-run below.

### 3.2 What must be re-run, and why

1. **Every molecular case, on the new physics.** The L7g corrections move the
   base chemical heat by 21.6 per cent and put the certified fiducial's energy
   residual at 4.27e-01 and its H2 carrier row at 7.817e-01; the archived
   molecular states are not roots of the present equations (L7g model memo
   section 3.5, L7h). This is the largest item and it is a campaign.
2. **The seven low-XUV atomic cases, with `Numerical flux: ROE`.** Their
   catalog `input.inp` files already carry the key (the 2026-09-17 decision);
   their `output/` states do not exist yet. Each is seeded from the rung named
   in its own `seed_from_ladder.txt`, or from the certified `.L25/` state where
   one was solved directly (L25 memo sections 3.3, 4.5).
3. **Every case's post-processing pass, through the L9 evaluate route.** The
   old pass overwrote the solved `Hydro_ioniz.txt` with a PLM relaxation
   snapshot and wrote a certification-claiming header onto the `_adv` files;
   `run_case.sh` now takes the evaluate route, and only a re-run leaves a
   certified state in both the main file and the `_IC` copy (L9 memo section 1).
4. **The `Hydro_ioniz.txt` snapshots the old pass overwrote**, for every case
   that is not otherwise re-solved: the certified state is in the `_IC` copy and
   the evaluate route rewrites the pair consistently.

### 3.3 The recipe of the re-run, step by step

1. `cd LHS1140b/models && ./run_campaign.sh` with the case list of the group
   being re-run; the molecular branch takes `SEED=` (a certified state of the
   same group) and the ROE cases take their `seed_from_ladder.txt` rung. The
   runner's post-processing pass is the L9 evaluate route and no longer sets
   `CFL: 1.0e-12`.
2. `python3 models/status.py --write` regenerates sections 7 and 8 of
   `LHS1140b/MODELS.md` from the files as they stand (the reason column and the
   flux column).
3. `python3 models/write_reproduce.py` refreshes every case's `REPRODUCE.md`
   (it records the two answers of the evaluate route and the host).
4. `python3 LHS1140b/make_memo_figures.py` regenerates
   `docs/figures/lhs1140b_*`; then `latexmk -pdf` in `docs/` for
   `lhs1140b_exhale_vs_pwinds.tex`.
5. The long solves go on `lart4` (the user's instruction of 2026-09-17: a new
   long solve goes to the idle host). The I4 runs measured a wall clock per
   pass about 8.5 times shorter there than on the loaded `lart3` (1161 s to
   136 s and 1821 s to 213 s on the first passes; L22 memo step 3, I4).

### 3.4 The not-solved list after the series

MEASURED 2026-09-17 from the `not_solved.md` files and the regenerated
`MODELS.md` tables:

| case | reason as its `not_solved.md` states it |
|---|---|
| `molecular_photochem_gj1132_kzzprofile/HeH2.09` | 21 outer passes on `c2e9c9990b9f` without an accepted state, the carrier balance row of the H2 front holding it; stopped 2026-09-17 because that binary carries the molecular physics of before L7g |
| `molecular_photochem_gj1132_kzzprofile/HeH9` | the run ended after 2 outer passes with no stationary claim (a relaxation snapshot written `certified=F`) |
| `molecular_scalar_gj1132_kzz1e9/HeH0.083` | one scalar movement bound over a column holding a slow H2 front and a far wind; the refusing row is the H2 carrier balance at cell 500, r = 29.031 R_p, 2.972e-03 against 1.0e-05, a band and not a fall (2.07e-03 at pass 11, 3.29e-03 at 24, 2.972e-03 at 40) |
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | the same bound; refused at pass 24 by the stagnation rule, the row 2.462e-02 at cell 500, the chemistry 3 per cent of it and the two advective face fluxes cancelling |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | the same bound; forty passes out, the row 2.477e-02 at cell 306, r = 1.9517 R_p, the refusing cell migrating outward at the pace of the front |
| `molecular_scalar_gj1132_wellmixed/HeH2.13` | never started: a molecular well-mixed case is seeded from another certified molecular well-mixed state and neither of the two that would supply one is solved |

The four low-XUV cases that had a `not_solved.md` before this series no longer
have one: the `.L25/ladder/` continuation certified them and each catalog
directory carries a one-line `seed_from_ladder.txt` naming its rung.

---

## 4. What is not going well, in detail

### 4.1 The two well-mixed molecular cases and `kzz1e9/HeH0.083` (the I4 result)

Three runs under `Coupled carrier solve: On stall` were restarted on `lart4` at
04:30 KST and stopped by the advisor on 2026-09-17 at 10:10 KST on the user's
instruction. READ from `LHS1140b/models/.L22/i4_*/run.log`:

| run | last pass | worst gated species row | where | hydro |
|---|---|---|---|---|
| `i4_wm0083` (`wellmixed/HeH0.083`) | 37 | 4.79e-02 of 1.0e-05, carrier balance H2 | cell 500 | `info = 2` at pass 37, mass 3.40e-03; `info = 0` and mass 1e-08 at passes 26 to 36 |
| `i4_wm055` (`wellmixed/HeH0.55`) | 16 | 2.84e-02 of 1.0e-05, carrier balance H2 | cell 445, the cell walking outward 350, 359, 369, 379, 389, 400, 411, 423, 435, 441, 444, 445 over passes 5 to 16 | `info = 2` at passes 13 to 16, mass 5.72e-03 to 9.89e-02 |
| `i4_kz0083` (`kzz1e9/HeH0.083`) | 40 (the cap) | 8.19e-04 of 1.0e-05, carrier balance H2 | cell 499 | `info = 0`, mass 8.8e-09; `state of the last outer pass, NOT accepted` |

The carrier H2 row falls about 10 per cent a pass in the two well-mixed runs
(8.16e-02 at pass 26 to 4.79e-02 at pass 37; 3.74e-02 at pass 13 to 2.84e-02 at
pass 16) and about 1 per cent a pass in the K_zz case. **No handover fired in
any of the three**, and in the one that ran furthest the reason is that it is
not stalling: neither of the handover's two counts reached three (L22 memo step
3, I4). So the criterion was not exercised, which is a statement about the runs
and not about the criterion. At 10 per cent a pass from 4.79e-02 a run would
need about 80 further passes to reach 1.0e-05; that is an extrapolation of a
rate that has not been constant and not a forecast. None of the three is a
state whose carrier row alone refuses it: the entry certification puts the
hydrodynamic energy row five to six decades outside its tolerance on all four
states the block and the alternation were started from (L22 memo step 3, I3-I5).

### 4.2 The coupled block's Krylov exhaustion, and what a follow-up needs

Section 1 has the numbers: 13 of 13 linear solves exhausting the 40-vector
subspace at relative residual 2.597e-01 against 1.00e-01, while the alternation
certified the same state at outer pass 12. The species box is not what holds it
(0 unknowns on a bound; I1 measured the fraction-to-the-boundary length at 13.8
times the direction). The block's own Jacobian action is sound: it reproduces a
central difference of the full residual to a mean relative error 3.692e-04 over
the rows a carrier direction moves and 3.792e-03 over the rows an energy
direction moves, the band holds all thirteen two-cell reconstruction entries
that the carrier matrix drops, and the outer ghost row the design expected to
FAIL passes bitwise, which is step 2c holding (L22 memo step 3, I1). Cost as
designed: about 64 full residual evaluations per block iteration (24 for the
banded preconditioner at 23 colors plus 40 products) against 1 to 2 products per
iteration of the alternation's three-unknown solve.

**What a follow-up would need is a preconditioner that holds the species rows.**
The design memo's risk table (`docs/lhs1140b_stationary_L22_step3_design_20260916.md`
section 8.3) names this risk first and names where it was met before: L24, 62 of
91 Krylov cycles at 0.91 to 0.999 of the right-hand side on the helium element
row of cell 335. The older record is N26c/N27 (`docs/ISSUES_20260909.md`): N27
found the helium row of cell 246 held by the LINEAR SOLVE, 40 of 40 products at
0.99 of the right-hand side with the Arnoldi image 15 per cent off and no longer
tracking the step length, and N26c found the atomic reload chaotic at the ulp
level, a one-ulp change of the upper ghosts moving `||R||` at iteration 40 from
3.2e-2 to 1.1, because the Krylov leg's one-digit tolerance decides which vector
crosses it at the last bits of the residual. I3 met the same failure mode on the
case the design expected to be the easy one. The design's fallback, a base-only
block, is kept as a named fallback to be opened only on a measured cost, and was
not designed further.

### 4.3 The base face flux is not the wind's on the LHS states

READ from the L26 memo section R9: the Riemann-solve base face mass flux, in
units of the window's mean flux, reads +32.38 on the fiducial, +336.55 on the
certified 0.10, -404.32 on the certified 0.03, -584.45 on the 0.02 transient and
-4735.68 on the cold start; the faces just above the base read -18.6, -203.7,
-153.9 and -232.5 before settling near +1 at faces 4 to 8. This refutes item
L21's statement of 2026-09-15 that the face mass flux is +1.00000 F_wind at
every face including the base, on the states as they now stand; L21's memo was
deliberately not edited, and the correction lives in L26. Mass is conserved, so
a stationary state should carry one flux through every face. The mismatch is
attributed to the base odd-even mode of `docs/p44_base_sawtooth.md` and L25
section 3, explicitly not to L26's closure, and it is reported and left
unfixed. **It passes certification**: the mass-row certification below r = 1.2
is gated only against a rounding anchor for cells of density 1e14, so an
imbalance hundreds of times the wind's own flux is not caught.

### 4.4 The 0.02 transient's energy row under Roe

The 0.02 transient was re-solved under the repaired boundary (which is bitwise
inert on that state), 40 passes each with `Numerical flux: HLLC` and `ROE`, in
`LHS1140b/models/.L26/x002_hllc` and `.L26/x002_roe`. **Neither certifies,
`info = 1` in both** (L26 memo section R6.4):

| row | HLLC | ROE |
|---|---|---|
| mass | 2.906E-01 at cell 1, tolerance 1.2E-06 | 5.947E-06 at cell 382, tolerance 3.0E-12 |
| momentum | 1.441E-08 at cell 2, tolerance 1.0E-08 | 6.076E-07 at cell 500, tolerance 1.0E-08 |
| energy | 6.299E-01 at cell 1, tolerance 1.0E-06 | **1.773E+00 at cell 6, tolerance 1.0E-06** |

Under ROE the base mass row DOES close, the worst mass row moving to cell 382
and falling from 2.9E-01 to 5.9E-06 and the momentum row moving to the outer
boundary; what remains is the energy row at cell 6, with the base at a different
thermal structure (the ghost at 226.02 K, `base.inp`'s `T_base` to the digit,
cells 1 to 5 at 220.4 to 203.5 K against 470 to 526 K under HLLC). That contrast
is the flux family of L25 step 3 seen on this transient, and it is that item's
question.

### 4.5 The two `molecular_photochem` cases stopped

`HeH2.09` was stopped at outer pass 21 and `HeH9` ended after 2 passes with no
stationary claim; both were running on the catalog binary, which carries the
molecular physics of before L7g, so a solution of it would have been superseded
by the re-run in any case. Their last written states are the seeds for the
re-run (their `not_solved.md` files).

---

## 5. Deferred items, with their reasons

- **The L7g cascade (the Bossion tables).** The level-resolved (v, J) cascade
  that would replace the bounded two-parameter form is not buildable now: the
  H2-H collisional rates above v = 4, together with the state-resolved
  three-body collisional dissociation for M = H whose detailed-balance inverse
  gives the nascent distribution, exist as a calculation (Bossion, Scribano,
  Lique & Parlant 2018, MNRAS 480, 3718, extended to all bound states by
  Nesterenok et al. 2019) but the rate tables are not distributed with the
  paper. The code keeps `q(R15) f_quench` with its validity stated where it is
  used. Because `f_quench` is 1 to a part in 1e6 to 1e7 in this layer, the
  cascade's absence does not matter here; it would matter for a shallower base
  or a thinner, cooler layer (L7g inventory memo sections 6, 8; the entry is
  already in `docs/TO_BE_DONE.md`).
- **L12 stages B to F.** Stage B is blocked: `binary_element_diffusion.f90`
  exposes neither `drift_and_gradient_face_coefficients`, nor `element_face_flux`,
  nor the element nucleus face flux itself, **which exists nowhere, because the
  operator forms two divergences and never a single face flux**; and the two
  divergences use two different discrete geometries (the diffusive half
  `Kj = 1/(rp^2 (rep_+ - rep_-))`, the advective half a volume
  `(rep_+^3 - rep_-^3)/3`) which agree only to O(dr^2). Stage B cannot be
  written without duplicating the element flux, and reconciling the two
  geometries moves every element row of every state. Stages C to F are
  sequenced after B and not started (L12b memo section 7).
- **L4e's two proposals** for the linear solve and the line search
  (`steady_newton.f90`), unchanged since 2026-09-14: the judged-row scaling
  exists as an option, is off, and was measured not to be the repair on its own;
  the ledger-ranked line search is not written.
- **`grid_and_gates/output_state_consistency.sh` row
  `outer_iteration_ending_is_the_stagnation_one`** expects `no_progress` and
  measures `pass_budget`. It runs with `Coupled carrier solve: False` and FAILS
  IDENTICALLY on every control binary of this series (L23, L9, L22 I1/I2 and
  I3-I5 all report it), so it is the outer loop's ending having moved under the
  row and not any item's doing. The test says what to do: pick a setting that
  reaches the ending again rather than loosen the assertion.
- **`steady_selfconsistent_residual`'s row `handback_matches_the_accepted_iterate`**
  reads 2.83103e-08 against a reference of 1; whether that is a defect or a
  property of a partially converged solve is not established (L22 step 3 I3-I5).
- **`global_parameters` declares `integer :: count`**, shadowing the intrinsic
  wherever a module uses that block without an `only` list. Reported by L15 as a
  latent compile hazard and not fixed, because the rename touches many files.
- **The h2 module notes.** `carrier_mass_fractions` in
  `binary_element_diffusion.f90`, the stage path's spelling of the same routine,
  states the inner ghost rule and no outer one (L22 step 2c, reported in a file
  that item could not edit). `low_mach_dissipation.f90`'s header still carries a
  nonnegative-dissipation claim that L25 step 1 proved wrong as written, and the
  correctly signed replacement `-M^-1 B^T K B` was derived and measured but not
  inserted (L25 step 1). `docs/input_schema.md` counts four channels under the
  `EXHALE_REACTION_HEAT_RECIPIENTS` key where the L7g memo counts three grouped
  corrections: the same physics, two counting conventions.

---

## 6. Decisions

### 6.1 Decisions the user made this series

| date | decision |
|---|---|
| 2026-09-16 | The renaming of the sixteen `backup/regression/` case directories carrying the banned noun, per `PLAN_20260916_rev3` section 10 |
| 2026-09-17 | **L25, option 2**: the Roe flux for the seven low-XUV atomic cases on the catalog grid, applied to `make_models.py`, `MODELS.md`, `write_reproduce.py`, `status.py` and the memo tex |
| 2026-09-17 | **L25 step 4**: seed each still-unsolved low-XUV case from a solved case at a higher XUV and lower the XUV in small steps, each rung seeded from the previous certified rung |
| 2026-09-17 | **L22 step 3**: the coupled-block design approved with all nine decisions D1 to D9 as recommended |
| 2026-09-17 04:30 | A new long solve goes to the idle host (`lart4`); the three I4 runs were restarted there |
| 2026-09-17 10:10 | The two I4 well-mixed runs stopped at outer passes 37 and 16 |

### 6.2 Decisions still open, each with a recommendation

| open decision | recommendation |
|---|---|
| Whether to re-solve the whole molecular catalog on the corrected L7g physics, and in what order | Recommend yes and soon: the archived molecular states are not roots of the present equations, so every molecular number in the memo and the update log stands on superseded physics. Order: the fiducial `kzz1e9/HeH2.13` first (it is the seed of the group), then `HeH0.55` and `HeH9.7`, then the two `photochem` cases from their stopped states, and the three well-mixed ones last, since L22 says they need something the alternation does not have |
| Whether to adopt the Roe flux beyond the seven low-XUV cases | Recommend NOT yet. R37 requires re-solving and re-judging every certified case on its own residual, conservation and grid check, and Roe's own grid convergence is not measured (only one case, only a factor 2, and the Roe leg at 1000 cells had not converged). The seven-case list is the right scope until that measurement exists |
| Whether the 1000-cell grid, rather than the flux, should be the catalog's answer to the low-XUV stall | Recommend measuring before deciding: the pair (ROE-500, ROE-1000) is missing, and the doubled grid moves both ends (`r_max` by 1.9 per cent), so it is not a pure refinement |
| Whether to open a preconditioner item for the species rows | Recommend yes, as the successor to L22: three separate measurements (N27, L24, L22 I3) now put the same failure at the same place, and nothing else in the solver has been shown to bind |
| Whether the base face flux mismatch (section 4.3) is a certification gap to close | Recommend opening it as its own item: a mass row hundreds of times the wind's own flux passes today, which is a statement about the gate as much as about the base |
| Whether `_adv` products should be regenerated for every archived case or only for re-solved ones | Recommend regenerating for every case in the same pass as the re-run, since the old pass also overwrote the solved `Hydro_ioniz.txt` |

---

## 7. Paths

**Memos of this series** (all under `docs/`):
`lhs1140b_stationary_L7g_inventory_20260916.md`,
`lhs1140b_stationary_L7g_model_20260916.md`,
`lhs1140b_stationary_L9_20260916.md`,
`lhs1140b_stationary_L12b_derivation_20260916.md`,
`lhs1140b_stationary_L15_20260916.md`,
`lhs1140b_stationary_L22_20260916.md`,
`lhs1140b_stationary_L22_step3_design_20260916.md`,
`lhs1140b_stationary_L23_20260916.md`,
`lhs1140b_stationary_L24_20260916.md`,
`lhs1140b_stationary_L25_20260916.md`,
`lhs1140b_stationary_L26_20260916.md`,
`cloudy_h2_model_reference_20260916.md`,
`named_case_audit.md` section 6 (the renaming).

**Plan and its reviews**: `docs/PLAN_20260916_rev3.md` (the plan of record;
`PLAN_20260916.md`, `_rev1.md`, `_rev2.md` are frozen and referenced by SHA-256
from rev3, and `PLAN_20260916_review.md`, `_review1.md`, `_review2.md` are the
three Codex reviews). The item rows are in
`docs/PLAN_20260913_lhs_stationary.md`.

**Audit directories**: `docs/audit_20260905/L15_stage2_20260916/` (the traces,
checkpoints and comparisons), `L22_step2_20260916/`, `L22_step2b_20260916/`,
`L22_step3_20260916/`, and the standalone probe
`docs/audit_20260905/plan_20260916_rev1_boundary_limit.py`.

**Run directories** (all under `LHS1140b/models/`):
- `.L22/`: `EXHALE_L22e.x`, `i3_alt/` (the alternation run of I3, certified and
  post-processed, with five `tpm_*.txt`), `i3_block/` (the coupled block, no
  certified state), `i4_wm0083/`, `i4_wm055/`, `i4_kz0083/` (each with
  `run.log` and the `run_lart3.log` of the host before the restart),
  `launch.sh`, `launch_i4_lart4.sh`, `pids.txt`, `pids_lart4.txt`,
  `pp_i3_alt.sh`.
- `.L25/`: step 3's seven direct ROE solves (`roe_s010_HeH9.7`,
  `roe_p020_HeH9.7`, `roe_s001_HeH2.13` certified; `roe_p001_HeH9.7`,
  `roe_p015_HeH9.7`, `roe_p025_HeH9.7`, `roe_s001_HeH9.7`), the grid study
  `grid1000_probe/`, and the drivers `launch_*.sh`, `post.sh`, `queue.sh`.
- `.L25/ladder/L{1,2,3,4}/x<tag>/`: the eighteen rungs of step 4, every one
  certified, with `input.inp`, `run.log`, `output/`, `REPRODUCE.md`,
  `seed.log`; `x<tag>_pp/` beside each target rung carries the post-processed
  copy and the He I 10830 measurement; ladder logs `L{1,2,3,4}.ladder.log`.
- `.L26/`: `fid_resolve/` (the fiducial re-solved under the repaired boundary,
  `info = 0`, certified), `x002_hllc/` and `x002_roe/` (the 0.02 transient,
  both `info = 1`).
- `.L21/L15/t{1,8}_r{1..5}`: the ten stage-1 repeats of L15.
- `.L14/x003_HeH2.13`, `.L14/x002_HeH2.13`: the frozen states the L26 probes
  read (read-only, never written to).

**Reference PDFs added this series**, all in
`/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/references/`:
`Rieper_2009JCP_228_2918.pdf`, `Rieper_2010JCP_229_221.pdf`,
`Rieper_2011JCP_230_5263.pdf` (L25); and for L7g and L7h:
`Bossion_2018_MNRAS_480_3718.pdf`, `Lique_2015_MNRAS_453_810.pdf` with its
supplementary directory `Lique_2015MNRAS_453_810_data/Rates_H_H2.dat`,
`Jozwiak_2024A&A_685_A113.pdf` with the CDS tables
`Jozwiak_2024J_A+A_685_A113/`, `LeBourlot_1999_MNRAS_305_802.pdf`,
`Wolniewicz_1998_ApJS_115_293.pdf`, `Wrathmall_2007_MNRAS_382_133.pdf`,
`Flower_2007_MNRAS_377_705.pdf`, `Glover_2008_MNRAS_388_1627.pdf`,
`Lepp_1983_ApJ_270_578.pdf`, `Palla_1983_ApJ_271_632.pdf`,
`Balakrishnan_1999_ApJ_514_520.pdf`, `Nesterenok_2019MNRAS_489_4520.pdf`,
`Orel_1987JCP_87_314.pdf`, `Schwenke_1988JCP_89_2076.pdf`,
`Esposito_2009_JPCA_113_15307.pdf`, `Paolini_2011PRA_83_042713.pdf`,
`Kokoouline_2001Nature_412_891.pdf`, `Strasser_2001PhysRevLett_86_779.pdf`,
`Guberman_1994PRA_49_R4277.pdf`, `Giusti-Suzor_1983PRA_28_682.pdf`,
`Takagi_2002_PhST_96_52.pdf`, `Mack_2006PA_74_052718.pdf`,
`Ohlinger_2007PRA_76_042712.pdf`, `Lee_2008_ApJ_689_1105.pdf`,
`Jozwiak_2024A&A_685_A113.pdf`.

**Data directories**: `LHS1140b/sed/` (the spectra, four new scalings this
series), `LHS1140b/models/` (the catalog), `LHS1140b/models_20260915_db87/` and
`models_20260914_preL21/` (the preserved earlier trees),
`observational_data/`, `cooling_data/`.

**The memo document**: `docs/lhs1140b_exhale_vs_pwinds.{tex,pdf}`. MEASURED
2026-09-17: the PDF was regenerated with `latexmk -pdf` so that the added L25
paragraph is in it, and it is **77 pages** (73 at the previous handoff). The
figures were NOT regenerated: the catalog states are unchanged.

---

## 8. Terminology

Four words this series used constantly, each with one meaning, so that the next
session does not invent a second one. `LHS1140b/MODELS.md` section 2 is the
naming contract.

- **rung**: one step of a continuation ladder, and nothing else. A composition
  rung is a step in reservoir He/H (the atomic ladders were built with steps of
  at most 0.3 dex, an empirical continuation limit); an XUV rung is a step in
  the spectrum scaling (L25 step 4's four ladders have three or six rungs each).
  A rung is not a case and not a pass; a rung's own directory carries a full
  `input.inp` and `output/`.
- **photochem**: the lower-boundary field of a group name, meaning the Photochem
  column handed over as `Lower atmosphere profile:` (T(p), K_zz(p), q_H2 and the
  element reservoirs at p_match = 1 microbar). It is a boundary treatment, NOT
  the code's molecular photochemistry module and NOT the `diffusive_photochemistry.f90`
  source file. A group can be `atomic_photochem_*` (atomic chemistry on a
  photochemical lower boundary), which reads as a contradiction only if
  "photochem" is taken for the chemistry.
- **the XUV grid**: the set of cases whose spectrum is the GJ 1132 proxy scaled
  by a factor, `gj1132x0.01` through `gj1132x0.33`. It is a grid in the spectrum
  scaling at fixed composition, not a grid in space and not the computational
  grid. When the number of cells is meant, say "the 500-cell grid".
- **pass**: one outer pass of the stationary alternation (the hydrodynamic solve
  followed by the composition relaxation). It is never a Newton iteration, never
  a Krylov product and never a marching step. `EXHALE_OUTER_PASSES` counts
  passes; `EXHALE_JFNK_MAXIT` counts iterations. "Pass 14" of a certification
  means fourteen alternations, each of which may hold hundreds of iterations.

---

## 9. Request list for external data

One item is blocking and the rest are wanted:

1. **BLOCKING for the L7g cascade.** The H2-H state-to-state rate tables of
   Bossion, Scribano, Lique & Parlant (2018), MNRAS 480, 3718, as extended to
   all bound states by Nesterenok et al. (2019), MNRAS 489, 4520, together with
   the state-resolved three-body collisional dissociation for M = H. The paper
   carries no supporting information, no repository and no data-availability
   statement, and the follow-on cooling function is deferred by its authors to a
   paper in preparation. The named authors to ask are Nesterenok, Bossion,
   Scribano and Lique. Without them the code keeps the bounded
   `q(R15) f_quench` form, which is measured valid in this layer and would not
   be in a shallower or cooler one.
2. The origin of Cloudy's `coll_rates_He_ORNL.dat`, cited in Cloudy only as
   "Lee et al. 2007, ApJ, in preparation", which ADS does not carry. Three
   plausible published leads were identified and none confirmed
   (`docs/cloudy_h2_model_reference_20260916.md`).
3. Moses & Bass (2000), JGR 105, 7013, still outstanding from the previous
   handoff: the source of the R17 Arrhenius branch above about 700 K, where
   nothing is measured.
4. A published He collider coefficient for the three-body H2 association. None
   exists, and the code currently uses argon's as a stated upper bound (Cohen &
   Westberg 1983); if a measurement or calculation can be obtained the bound
   becomes a value.
