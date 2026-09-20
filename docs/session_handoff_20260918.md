# Session handoff, 2026-09-18: the series on `docs/PLAN_20260917.md`, items L27 to L37

Written 2026-09-18 at the close of the series. The previous handoff is
`docs/session_handoff_20260917_rev1.md` (the revision written after
`docs/session_handoff_20260917_review.md`, the Codex review of the first
2026-09-17 text); its format is the format of this one. The plan this series
executed is `docs/PLAN_20260917.md`. The plan that follows from this handoff
is `docs/PLAN_20260918_rev1.md` (itself the revision of `docs/PLAN_20260918.md`
after `docs/PLAN_20260918_review.md`).

**Every number in this document is READ from the item report, the memo, the
run log or the source line named beside it. Nothing was measured for this
document, and no binary was run for it.** Where a worker labeled a number
MEASURED in its own report, it is quoted here as READ from that report and the
report is named; the measurement is that worker's.

**The binaries.** Three identifiers appear below and they are not one thing.
The plan stood on `EXHALE.x` md5 `59bfdb3fc4d0104fc2e9c3734596d2f6` (manifest
`LHS1140b/models/BINARY_MANIFEST_59bfdb3fc4d0.txt`), which is the binary item
L29 reclassified the archived claims on. The catalog re-run of item L34 was
made on `LHS1140b/models/EXHALE_3146d11b.x`, md5
`3146d11b4090306dcea75bb9718edd22`, manifest
`LHS1140b/models/BINARY_MANIFEST_3146d11b4090.txt`; every catalog state written
on 2026-09-18 carries that md5 in its `REPRODUCE.md`. The tree's `EXHALE.x`
was rebuilt by the regression harness at 05:30 to
`267e0b2e85be71e7d66c221ecc58cf34` (READ, the L30c and L35 reports) while items
were still landing, and the sources have moved past it again since.

**FINAL BINARY: `74b96cdcf887dcfee301ec547f012d88`** (MEASURED by the advisor 2026-09-18: `make distclean` then bare `make`, `make -q` returning 0; the same md5 the rename increment's measured build reported; manifest `LHS1140b/models/BINARY_MANIFEST_74b96cdcf887.txt`)

**GOLDEN REFRESH: none in this series (the user's instruction of 2026-09-18 14:10: the corrections of `docs/PLAN_20260918_rev2.md` continue on the same tree, so a refresh now would be repeated within days; the goldens stay the 2026-09-17 set, every item's control is its own entry-text build, and one refresh follows when the rev2 items have landed).**

No worker refreshed a golden, and no worker wrote under `backup/regression/`.
The update-log section that records the refresh is written PENDING and carries
a marked placeholder block for the refresh date and the final md5 (READ, the
update-log report).

---

## 0. What the series did, in one page

**Where it started.** The revised handoff of 2026-09-17 folded the review of
the first handoff text and carried its verification table: which of the
review's findings were read in the source or in its retained logs and which
were accepted as its judgment. Three of them set this series going. The L26
base face flux readings (+32.38 and +336.55 of the wind's flux) were PLM
evaluations of states solved under WENO3, so the certification-loophole claim
was withdrawn and a diagnostic contract was needed instead (item L27). The
stationary log reader joined records of different solves (item L28). And the
element and carrier transport operators divided advection and diffusion by two
different discrete geometries, which is a conservation inconsistency of the
implemented equations (item L30).

**What the plan asked for, and what landed.** `docs/PLAN_20260917.md` holds
items L27 to L36 and a bookkeeping section. All of them landed, and four of
them grew follow-ups that landed in the same series: L30b and L30c after L30,
L34a, L34b and L34c for the catalog, L36c, L36d and L36e after L36, and L37,
opened on the finding L34b reported.

- **L27**: one routine installs the stationary operator and one entry point
  evaluates a stored state; the certification report carries a face-flux
  budget as a report line. The probe now reads 0.99998 of the wind-window mean
  where it read 32.376.
- **L28**: every stationary record names its solve; the reader reads one
  completed solve. The row that read 2.83103e-08 reads 0.842987.
- **L29**: 119 archived certification claims labeled under the binary of
  record, 101 reproduced and 4 refused.
- **L30 / L30b**: one spherical face area and one shell volume for every
  transport term and for the hydrodynamic rows, the element face flux exposed
  as one public object, and the Runge-Kutta species conservation row repaired
  in its fixture rather than in the operator.
- **L30c**: the order-unity movement of `oxygen_chemistry` diagnosed as a
  discrete branch selection amplifying a rounding seed on a snapshot that is
  not a solution.
- **L31**: the molecular energy cycles audited with no double count, the four
  approximations relabeled with their validity, and a reduced level ladder
  built that brackets the scalar quench fraction.
- **L32**: the coupled block's directional Jacobian action measured and found
  not to be a linear map of its direction, so no preconditioner is the item.
- **L33**: the three stalled outer relaxations classified from the existing
  logs; the handover route was never exercised.
- **L34a / L34b / L34c**: 79 atomic states reproduced, seven low-XUV Roe cases
  certified, the molecular reference solution and `HeH9.7` certified, and six
  molecular cases not certified, two of whose seeding routes the restart
  contract forbids.
- **L35**: the matched-domain flux and grid study, and the local base closure
  compared against the retained remote one.
- **L36 / L36c / L36d / L36e**: ionization stage transport built, given a
  caller, extended to both elements, and carried in an atomic run, where the
  state certifies.
- **L37**: the molecular base mass row's round trip measured, its band and its
  mechanism named.
- **Bookkeeping**: the `low_mach_dissipation.f90` header corrected, the
  stagnation-ending row made reachable again, the `TO_BE_DONE` L7g entry
  extended, and `global_parameters`' `count` renamed `marching_step`.

**What is left as decisions.** `docs/PLAN_20260918_rev1.md` holds nine items (D8, the runner's continuation and state-pair contract, and D9, the pass budget of a transported-ionization solve and the four key-on goldens, were added from the L34c and L36e reports),
five of them the first text's corrected (D1 the default base cell width and
its input key; D2 the base boundary at zero and reversed flow; D3 the
transported He II row and its acceptance norm; D4 the thermochemical closure's
returned state; D5 one boundary-state operation and the restart that reproduces
it) and two added by the review and by L36e (D6, the helium held in HeH+ is not
reserved from the helium ion-stage simplex, a conservation defect in the
expressions; D7, metal charge exchange enters the sweep's rows but not the
transported stage sources). D6 and D7 are items to execute rather than
decisions to take; the decisions are in section 6.2.

---

## 1. The items, one by one

### L27: the diagnostics evaluate the production stationary operator

A new module `src/modules/states/stationary_operator.f90` (200 lines) holds
`select_stationary_reconstruction`, which installs `WENO3` and returns the
previous selection, and `stationary_face_mass_flux`, which installs the
composition and caloric state, applies the boundary to a copy of the conserved
array, reconstructs and returns the face mass flux on every face together with
the operator identity. The nine literal WENO3 assignments of
`src/EXHALE_main.f90` (entry-text lines 1231, 1357, 1497, 1531, 1604, 1938,
3657, 3728, 5464) call it; the marching `PLM` at 1873 and the residual probe's
save and restore at 4617 to 4627 keep their own meaning. The certification
report gained an optional face-flux block that gates nothing and changes no
tolerance. READ from the L27 report: the rewritten probe prints a base face
flux of 0.9999843726323676 of the window mean on `.L26/fid_resolve` and
0.9999748281788378 on the certified `x0.10_kzz1e9/HeH2.13` state, against
32.376280 and 336.551906 on the one explicitly labeled non-stationary row, and
the production evaluate route prints the same bits as the probe on both states,
face indices included. The refactor is byte-identical on the data lines of
`wasp_full` (12251 steps), `mol_base_handoff` and `lower_profile`. Suites on
the measured build: `grid_and_gates` 253 PASS with the one pre-existing
stagnation-row FAIL, `certification` 84, `krylov_and_dogleg` 334,
`steady_species_rows` 195, `adv_static_limit` 53, `element_operator` 40.

### L28: a stationary record names its solve

A module counter in `src/modules/time_step/steady_newton.f90` (+26 / -8, print
statements only) appends ` solve=<n>` to the JFNK and PTC iteration lines, the
best-iterate restore line and both completion lines. The reader of
`src/tests/steady_selfconsistent_residual/run.sh` (+204 / -25) names the solve
by the last completion line, requires a record to lie between the previous
completion line and this one and, where the token exists, to carry it, matches
the certification block positionally, and prints `incomplete_evidence` rather
than a ratio built across solves. READ from the L28 report: the row
`handback_matches_the_accepted_iterate` on `LHS1140b/models/.L22/i3_alt/run.log`
moves from 2.83103e-08 (FAIL) to 0.842987 (PASS), which is 1.208e-08 over
1.433e-08, the two numbers standing at lines 2284 and 2232 of that log. Four
synthetic fixtures cover multiple solves, a restore only in an earlier solve,
no restore in the final solve and a truncated final solve. All 437 archived
`run.log` files under `LHS1140b/models/` were re-read with both readers: 156
carry no stationary solve, 142 give the same ratio, 7 change, 41 are runs
stopped inside a solve and now answer `incomplete_evidence`, 16 have no iterate
record in the chosen solve, and 75 have none anywhere. No interval or tolerance
was changed.

### L29: the archived certification claims, labeled

119 archived states were evaluated on `59bfdb3fc4d0104fc2e9c3734596d2f6`
through the zero-step evaluate route, on copies, with nothing written into a
case directory. READ from the L29 report: **101 reproduced, 4 refused, 6 state
no claim, 7 have no state, 1 is not evaluable.** The four refused are
`.L14/x003_HeH2.13` on the hydrodynamic mass row, 9.984E-01 at cell 2 against
3.2E-10, and the three certified `molecular_scalar_gj1132_kzz1e9` states on the
hydrodynamic energy row, 4.634E-01 at cell 257 (`HeH0.55`), 4.267E-01 at cell
203 (`HeH2.13`) and 3.597E-01 at cell 1 (`HeH9.7`), all against 1.0E-06. The
x003 evaluation is line-identical to the review's retained log. New:
`LHS1140b/models/reclassify_claims.sh`, which refuses any other binary, and
`LHS1140b/models/CLAIMS_59bfdb3fc4d0.md` with the 119-row table; `status.py`
gained one column read from that table and nothing else. Of the 119, 90 carry a
`Hydro_ioniz.txt` whose own coupling line says `certified=F
cert_reason=no_stationary_claim`, which is the L9 snapshot; the claim was read
from the `_IC` copy wherever `Hydro_ioniz.txt` says `recon=PLM` (110 rows). The
one state the route could not evaluate is
`molecular_scalar_gj1132_wellmixed/HeH2.13`, whose mapped seed carries
`he_diff: T` while the case does not.

### L30: one spherical geometry, and the element face flux exposed

`grid_construction` now owns `spherical_face_area_and_cell_volume`, with
`face_area(f) = r_edg(f)^2` and `cell_volume(j)` in the factored form
`(r_+ - r_-)(r_+^2 + r_+ r_- + r_-^2)/3`; the four element diffusion sites and
the two carrier diffusion sites divide by that volume instead of
`r_j^2 max(dr_j, 1 cm)`, the 1 cm floor on the width is gone, and a nonpositive
volume is refused by the grid constructor. `element_nucleus_face_flux` is the
one public object returning the advective element mass face flux, the diffusive
face flux, the face element nucleus density and the face mass per nucleus.
READ from the L30 report, with only the six denominators reverted so that the
numbers isolate the geometry: `element_row_is_the_divergence_of_the_exposed_flux`
2.40865e-11 to 0.0 bitwise (tolerance 1e-14) and
`element_column_sum_is_the_boundary_flux_difference` 1.38185e-09 to 2.03150e-19
(tolerance 1e-13), both relative to the operator's own row scale; the second is
the conservation defect itself. Impact on three states through the evaluate
route: **no certification verdict flips**; the largest movement of a gated row
is the molecular fiducial's elemental transport row, 5.725e-09 at cell 217 to
7.245e-08 at cell 315, which sits 138 times below its tolerance on both builds.
Regression movement: `mol_diffusion` `H3p` 1.026e-04, `mol_carrier` `H3p`
8.039e-06, `lower_profile` `OIII` 4.421e-07, `mol_base_handoff` and `wasp_full`
byte-identical (neither runs a changed operator). Three conservation rows of
`diffusion_tests` had to follow the operator, because the fixture summed
helium mass on the old diffusive weight; 35 of 35 pass after.

### L30b: one spelling of the volume, the failing budget row, the flux writer

Eleven further sites spelled the shell volume out as the difference of cubes,
all of them the same exact volume but in the cancellation-prone form. A public
`spherical_cell_volume(j)` now holds the one expression and the array form
fills its `cell_volume` from it; seven of the eleven assemble rows over the
padded ghost range and two are scalar functions called once for each cell, so
the array form could not serve them. READ from the L30b report: the two spellings differ
by up to 6.4e-13 at the base of the refined `hydrostatic_residual` fixture,
measured. The failing row
`h2_nucleus_total_changes_only_by_boundary_flux` was diagnosed as the fixture
and not the operator: it built its own face composition and continued its own
cell N+1 value past the outflow face where `carrier_mass_fractions` states
`X_ghost = X_N`. With the fixture reading the operator's two routines and the
production objects unchanged, it moves from 9.9166E-09 (FAIL) to 1.3739E-15
(PASS) against 1.0E-12. `element_flux_profile.txt` is now written from
`element_nucleus_face_flux` for both halves, after `project_elements`, so the
file is the flux of a state that exists, and `F_H + F_He = Mdot_face` holds by
construction. Movement of task A on the regression cases: `hydrostatic_column`
8.46e-15 of the column scale, `wasp_full` 4.57e-13 over all 12251 steps (both
builds stopping at the same step with the same final line),
`mol_base_handoff` 1.19e-05 and `mol_diffusion` 2.06e-05 of the column scale.
Cost measured on two cases: 0.998 and 0.999 of the control wall clock. Reported
on the way: `src/utils/fortdep.py` did not scan `include`, so
`hydrodynamic_rows.o` did not depend on the `.inc` file it is compiled from;
that has since been fixed (`fortdep.py` now emits the included file as a
prerequisite and scans it for its own `use` statements).

### L30c: what moves `oxygen_chemistry` by order unity

Verdict: a discrete branch selection amplifying a rounding-level seed on a
snapshot that is not a solution. READ from the L30c report: the seed is the
restatement of the shell volume that L30b carried into the hydrodynamic rows,
whose two spellings differ by at most 5.7707e-13 relative on that case's own
grid (median 1.1853e-14). Two bitwise experiments name it: a build with only
the carrier diffusion denominators reverted is bit for bit the measured build,
and a build with only the shell-volume restatement reverted is bit for bit the
control, in all four output files and at a cap of one step as well. The first
difference above 1e-12 is inside step 0, in the first ionization sweep, in the
accepted composition; `du` at step 1 is already apart by 6.4e-04 and at step
1000 by 2.6e-02, a factor 41 over 1000 steps, so it is a jump and not a growth.
The integral the case does pin is stable: `log10 Mdot` 9.61 on both builds.
`mol_sec_ion` is the other class: the first step is bitwise identical, the
second is apart by 9e-12, the separation rises smoothly to 2.3e-08 at 1000
steps, and at 12000 steps both builds stop at the same step with the same line
and the matrix movement is 3.2102e-04 in `v`, which is 6.31e-07 of that
column's own scale. The root ladder does not fire at all in that case.

### L31: the molecular validity contract

READ from the L31 report and its memo
`docs/lhs1140b_stationary_L31_energy_cycles_20260917.md`. Step 1, the
reaction-cycle audit: five complete cycles close on the production assembly to
2.2e-13 eV or better, the nascent H2+ vibrational excitation of R23 deposited
once and not again at R5; an injected 0.55 eV double count, exactly the v = 2
energy the R23 comment names, is caught by three cycle rows at exactly that
size. No double count was found in a file the item does not own. Steps 2, 3, 5:
the R5/R16 n = 2 recipient is relabeled a LOW-STATE RECIPIENT APPROXIMATION
with a startup validity line, the R6 internal energy a MODEL ESTIMATE from the
position of a published peak, and the helium third body an ARGON-BASED ESTIMATE
over 77 to 5000 K with Cohen and Westberg's +/-0.3 in log; the claims to be
bounds are removed. Step 4: the new suite `src/tests/h2_level_ladder/` solves
54 rovibrational levels with one explicitly uncertain group above the data, 32
of 32 assertions passing (Boltzmann departure 2.97e-14 / 2.92e-15 / 6.12e-15 at
808 / 2128 / 5083 K, the steady partition closing to 3.7e-13 relative). Its
bracket on the radiated share `1 - f` against the scalar form: 6.793198e-07
against at most 5.140654e-08 at the base, 1.006850e-05 against 2.004984e-06 at
the front, 5.963722e-04 against 2.146044e-04 in the dilute column, ratios
13.21, 5.02 and 2.78, so the scalar radiates more than any member of the
bracket and replacing it could only add heat; the scalar form stays the default
and no default changed. Step 5's replacement of the argument "above the front
the H2 has gone": the R15 source at the 76 cells above 5000 K is 2.3404e-08 to
3.5145e-06 of the other three H2 formation channels on the certified fiducial,
and at most 1.6745e-13 of the hydrogen nuclei on the hottest atomic wind of the
catalog (`atomic_scalar_gj699_wellmixed/HeH1000`, maximum T 9517.55 K). Step 6,
the uncertainty record on three cells, as fractions of the nominal molecular
chemical heat: the **helium third body +41 / -21 per cent at the molecular
depth**, the recipient +4.6 per cent at the front and +130 per cent in the
dilute column where the heat is 3e-12 erg cm^-3 s^-1, and the R6 shares and the
quench fraction at 1e-07 or below everywhere. Step 7: the collider sum agrees
bitwise in the chemistry and in the carrier row and to 1e-14 relative in the
heat ledger. Impact: every output file of `mol_base_handoff` byte-identical
between the control and the measured build apart from the provenance
timestamp; the only change to a run is three validity lines of standard output.

### L32: the frozen coupled Newton system

READ from the L32 report and `docs/lhs1140b_stationary_L32_20260917.md`.
Verdict (i): **the directional action is not a linear map of its direction**, so
the Krylov solve has no theory to stand on and neither a wider subspace nor a
preconditioner can repair it. The residual and the action are bitwise
repeatable across an interleaved unrelated state; `A(3v) = 3A(v)` to 1.400e-14
and each Arnoldi column is a correct product on its own direction to 3.523e-16;
but `||A(v1+v2) - A(v1) - A(v2)||` is 1.867e+02 of `||A(v1)+A(v2)||` on the
first two Arnoldi directions, with the defect scaling as the arc to the power
2.291, which is curvature and not a floor, nine tenths of it in cells 1 and 2
and 87 per cent in the momentum rows. The reduced least-squares residual falls
monotonically (0.985, 0.372, 0.259 at products 1, 20, 40) while the true
residual of the same iterates does not (0.373 at 20, 0.511 at 40), the Arnoldi
image being 42 per cent off the operator's image of the step; the factor
between what the block reports and what its step has is 1.98. Subspaces of 80,
160 and 320 vectors all stop at 70 products at the same 2.81563e-01. The
chemistry stopping tolerance from 1e-6 to 1e-12 leaves the floor at 3.189e-01
to four digits. The banded model omits no long-range coupling: the part of six
sampled columns outside the band is exactly zero and the band's own entries are
0.09 to 3.5 per cent off. Two corrections to the record: `2.597e-01` is the
FIRST linear solve's reduced residual, the largest of the seventeen legs, and
trust-region iteration 14 was an ACCEPTED step at `1.337e-01`, the rejected ones
being 13 and 16. The item also delivered the loop-top-stop iterate record that
L28's reader needs, print only, with byte-identical data lines on
`carrier_model_a_newton`.

### L33: slow outer relaxation or a stalled mode

Read-only, from the existing logs; no solve was run. READ from the L33 report
and `docs/lhs1140b_stationary_L33_20260917.md`. `i4_wm0083` and `i4_wm055` are
**front relaxation** (the refusing cell of `i4_wm055` walks outward 315 to 445,
that is 2.1136 to 11.7471 R_p, and the row decays from pass 5); `i4_kz0083` is
**neither of the plan's two**, a relaxation throttled by the composition
movement bound with a fixed cell and unbounded row. The control that makes the
reading possible is `i3_alt`, the only one of the four in which the carrier
relaxation reached its own fixed point (2.87e-09 to 4.70e-09 at every pass,
ending "the carriers stopped moving at the fixed wind" 11 of 11 times); in the
three refusing runs every relaxation of every pass was cut off by the movement
bound. The three progress factors of the revised handoff are confirmed
(0.95273, 0.91232, and 1.04 per cent a pass over passes 23 to 40), with two
readings they leave out: `i4_wm0083` ROSE by a factor 1.0611 a pass over passes
1 to 22 before it fell, and the rate follows the movement bound (2.49 per cent
a pass at a bound of 2.5e-03 against 1.04 per cent at 1.25e-03, the decrement
falling to 0.42 of itself where the bound halved). **The handover was never
exercised**, and the revised handoff's wording is corrected: the bound-ending
count reached three in all three runs and stood at 37 of 37, 16 of 16 and 39 of
39; the count that never reached three is the consecutive-no-progress one,
whose longest run was 1, 2 and 1. The item also designed a five-rung
continuation in `q_H2_base` (0.190110 to 0.476179 in four steps of 0.0715) from
the certified fiducial, and did not run it.

### L34a: the atomic catalog re-run

READ from the L34a report and `docs/lhs1140b_stationary_L34a_20260918.md`.
**79 archived atomic states reproduced their own stationary claim on
`3146d11b4090306dcea75bb9718edd22`, 0 refused, 0 re-solved**, so the re-solve
branch never ran. The states did not move: mass density bit-identical in every
cell of every state, velocity moving in 32 of 79 by at most 2.219e-16 relative,
temperature by a median 5.836e-14 and at most 3.807e-12, and `4 pi rho v r^2` of
the outermost physical cell unchanged to the last bit in all 79. The He I 10830
equivalent width moved by a median 5.4e-06 relative and at most 3.6e-05, which
is the size L9 measured for the change of post-processing route. The gated
elemental He/H partition row moved by a median 0.9 per cent and at most 16.2
per cent against the value each case's own solve accepted on the predecessor
binary, with the same refusing cell in 56 of the 64 states that run element
diffusion; every state still certifies, and the movement is NOT attributed to
L30 alone, because the two numbers come from two binaries and two operations.
The L27 face-flux budget over the 79: base face within 3.5e-05 of the
wind-window mean, face-to-face spread at most 5.05e-08 and median 1.96e-09,
operator identity `WENO3, HLLC, well balanced: T` on all 79. **The seven
low-XUV cases with `Numerical flux: ROE` all certified at outer pass 1**,
`info = 0`, `||R||` 6.245e-08 to 8.546e-07, each seeded from the rung its
`seed_from_ladder.txt` names or from its certified `.L25/` state; against the
L25 numbers the mass-loss rate agrees to at most 5.6e-04 relative and the
equivalent width to at most 8.3e-05. Nine flux-closure rung `input.inp` files
were repaired (`Solver: Newton` and `Restart intent: stationary` restored,
`Do only PP` set False), each recorded in its own `REPRODUCE.md`.

### L34b: the molecular reference solution

READ from the L34b report and `docs/lhs1140b_stationary_L34b_20260918.md`.
**`molecular_scalar_gj1132_kzz1e9/HeH2.13` CERTIFIED at outer pass 12**,
`info = 0`, `||R|| = 1.417e-08`, in 40 m 35 s at 8 threads on `lart4`, seeded
from its own archived certified state (which the current binary refuses on the
energy row at 4.267e-01, exactly as L29 records). Every gated row inside its
tolerance: momentum 5.588e-12 of 1e-08, energy 1.417e-08 of 1e-06, carrier
balance H2 8.172e-06 of 1e-05, elemental transport He/H 6.379e-08 of 1e-05, He
2^3S level 1.824e-19 of 1e-06, eliminated-species closure 1.665e-16 of 1e-06.
The L27 budget on the solved state: the face mass flux is one number across the
whole column, base face included, to 1.39e-08 of the wind's. The energy balance
by the code's own operator with no step taken closes at a cell-wise maximum of
1.4464e-08 at cell 31 against 1.0e-06. `HeH9.7` certified at pass 10. What the
corrected L7g physics did: `log10 Mdot` 7.90907 to 7.90239 (-1.53 per cent in
Mdot), He I 10830 equivalent width 1.57583 to 1.56032 %A (-0.99 per cent), base
T 808.319 to 778.792 K (-3.65 per cent), base heating rate -5.60 per cent,
`2 n(H2)/n_H` at cell 1 -5.61 per cent, the H2 front one cell inward; on
`HeH9.7` the same corrections give -0.84 per cent in Mdot, -0.62 per cent in
the equivalent width, -7.99 per cent in base temperature and -12.80 per cent in
base heating, so **the colder the base the larger the movement, and the line
forms far above the layer that moved.** The L31 bracket on the solved state
puts the helium third body at +42.0 / -21.1 per cent of the base chemical heat,
**eight times the change the re-solve itself made in it.** One acceptance item
does not come back clean and is the origin of item L37: the no-step evaluate
route refuses the state on one entry of seven, the cell-1 hydrodynamic mass row
at 1.248e-08 against that cell's own anchor 7.8e-09, a factor 2.39 across the
re-entry against the 0.78 to 1.71 band L18 measured on 120 certified atomic
states.

### L34c: the remaining molecular cases, and the seeding rules

READ from the L34c report and `docs/lhs1140b_stationary_L34c_20260918.md`.
**None of the six certified.** One new certified solution was produced as a
seed, `.L34/atomic_wm2.13`, `info = 0`, `||R|| = 2.476e-08`, `log10 Mdot` 7.86,
equivalent width 1.9533 %A. The substantive finding is that two of the three
seeding routes the brief offered are inadmissible by the restart contract.
First, a metal-free scalar molecular state cannot seed a `photochem` case:
`load_IC` refuses it because that group carries C/H, N/H and O/H reservoirs at
the matching level and a scalar case states one ratio. Second, a diffused
molecular state cannot seed a well-mixed run even with `Restart option change:
he_diff` in place, because the state carries HeH+, whose nucleus of each element
cannot be rescaled by one factor, and the reference state's He/H runs from 2.13
at the base to 6.042513123E-01 at cell 500; the same clause forbids a
composition step inside the well-mixed family, so the only entry to a molecular
well-mixed case is the atomic-to-molecular conversion of an atomic well-mixed
wind at that composition. In the three well-mixed cases the stationary wind is
found and held (`hydro info = 0` from pass 3 to 5 onward) and the H2 carrier
balance alone refuses, with every relaxation of every pass ending on the
composition movement bound: that isolates the carrier transport as the single
obstruction there. `wellmixed/HeH0.55` reached 3.76e-02 on the gated carrier
row, the nearest any molecular well-mixed case has come to 1.0e-05. The L33
ladder's rung 1, a 0.179 dex composition step from the certified reference,
ended its first solve at outer pass 2 of a budget of 25 because the element
composition relaxation found no admissible advance and restored its entry
composition, which ends the ladder by its own acceptance.

### L35: matched-domain flux and grid, and the local base closure

READ from the L35 report and `docs/lhs1140b_stationary_L35_20260918.md`. The
three grids were matched to the bit: `r_edg(0)` and `r_edg(N)` agree to 0.0 at
500, 1000 and 2000 cells, solved for the two keys on a Python replica that
reproduces the binary's radius column exactly. Holding the physical faces fixed
also holds the first cell width fixed (1.933188e-04 in all three), so a
matched-domain refinement cannot refine the base. **The grid does relieve the
stall**: on the matched domain the unmodified HLLC flux has no stationary root
of `atomic_scalar_gj1132x0.10_kzz1e9/HeH9.7` at 500 cells (mass row 1.801e-01
and energy row 2.406e-01 at cell 1, base face carrying 0.700465914 of the
wind's own mass flux with a 36 per cent face spread), and at 1000 and 2000
cells the same flux brings the base face to within 4e-07 of unity at the first
outer pass. **It does not remove the flux-family difference, which is
converged**: T(1) -10.32 per cent, rho(1) +11.08 per cent, Mdot -0.996 per
cent, He I 10830 equivalent width -1.79 per cent (Roe against HLLC), identical
at 1000 and 2000 to three digits, and the profile difference is 10 to 24 per
cent below 1.01 R_p and 0.5 to 0.8 per cent out to the top of the domain.
**Refining at the fixed domain costs the certification**: none of the four
refined runs certifies, all four carrying a hydrodynamic energy row of 2.7e-04
to 8.2e-04 against 1.0e-06 at r = 4.3 to 5.0 R_p, at both fluxes, and the only
certified state of the six is Roe at 500 cells. What the item adds to the open
R37 question is the Roe flux's own grid convergence: it is converged at the
catalog's 500 cells, moving by 2.5e-06 in T(1) and 1.05e-04 in Mdot at a
doubling. Step 2: the local characteristic closure keyed on the production base
face flux and the retained remote-window closure install a BIT-IDENTICAL ghost,
all three primitive components, on both certified states, on a base whose
velocity is cut by a hundred, on a reversed base, on a perturbation confined to
the remote window and on the hydrostatic column, and the local one is a fixed
point of its own map in one iteration. They differ on two kinds of state only,
and neither is one a converged run holds: at rest, where the window carries
exactly zero flux and the smoothstep returns 1/2 so the base density is the
average of the two isentropes (3.79 per cent below the reservoir on the
fiducial, 3.07e-07 on the hydrostatic column), and on a 1.5x extended domain
(-7.30e-02). The remote closure's nonlocality is not a reason for an item: over
a factor ten in the window location the window mean moves by 1.3e-05 and the
face state does not move at all.

### L36: ionization stage transport, stage B

READ from the L36 report and `docs/lhs1140b_stationary_L36_20260918.md`. One
new module `src/modules/functions/ionization_stage_transport.f90` (518 lines)
holds the stage face flux on the element operator's own
`element_nucleus_face_flux`, its divergence on the L30 geometry, its
tridiagonal Jacobian rows and the simplex projection. The suite grew from 15
rows to 25, 25 passed, with the production sum identity at 2.0245e-16, the
identity at the boundary faces 1.0621e-16, after the simplex projection
2.1653e-16, and the stage divergence column sum 1.2421e-17; three rows are RED
by construction and measure the three mistakes the construction avoids (a
cell-centered volume 8.0694e-08, a full eddy term in every stage 1.4761e-04,
bounding the closing face value 3.6208e-04). At that point the module had no
caller in a run. The one-velocity approximation was measured and written at the
site: `|w/v| < 0.1` only inside 1.220 R_p (H) and 1.243 R_p (He) without the
resonant charge-exchange channel and 1.231 and 1.273 with it, reaching 0.625
and 0.430 over r >= 1.2 R_p, while `t_coll/t_flow < 0.1` out to 14.51 and 8.25
R_p, so the stages are collisionally locked in time scale and still drift,
because the drift is set by the ambipolar field against the friction. Stage F
landed on the transit side: `EXHALE_TRANSIT_STATE` selects the input PAIR and
every product header carries the selection, the files and the measured distance
to the other state; the He I 10830 equivalent width on the certified atomic
fiducial is 1.5544 %A from the `_adv` profile against 20.3950 from the
local-equilibrium solution, a ratio of 13.12.

### L36c: the stage row solved in the run

READ from the L36c report and `docs/lhs1140b_stationary_L36c_20260918.md`. The
H II row of the transport-chemistry operator is now the stage row, written on
the element operator's own hydrogen nucleus flux, and the bulk-mass `ic_Hp`
spelling is retired, so the module has a caller. The element's diffusive half
enters only where the run transports the elements against each other; the
stage's own eddy term is carried in every configuration. Key-off byte identity
on `wasp_full` (400 steps), `mol_diffusion` and `mol_carrier` (12000 each);
key-on movement 3.02e-05 (`hp_front`), 4.22e-03 (`hp_zero_seed`) and 3.03e-03
(`hp_trace_seed`) in `Hydro_ioniz.txt`, where `x(H II)` itself moves by 1.3e-10
to 1.4e-07 and what moves is the helium partition and the neutral hydrogen at
the front. Two of those goldens were already stale against the entry text by
3.3e-03 and 3.9e-03. On `carrier_model_a_newton` the stage sum identity holds
on real solved states at 2.268e-16 and 2.196e-16; the state is not certified on
either build, the `carrier balance H+` row standing at 6.432e-03 measured
against 6.392e-03 control, two and a half decades above its gate in both
spellings. A proposed anchor of 1e-13 for the non-gating stage-sum entry is
recorded as a proposal, not set.

### L36d: both elements' stages, and the anchoring measurement

READ from the L36d report and `docs/lhs1140b_stationary_L36d_20260918.md`. The
carried set became x(H II) per hydrogen nucleus and x(He II), x(He III) per
helium nucleus, each on its own element's nucleus flux, the two helium stages
sharing one simplex and one closing stage, with one stage-sum certification
entry per element. **The anchoring measurement that L12a and L36c owed was
made**: a manufactured ionization column whose faces are the same two radii at
every resolution gives a row measure of 7.2833e-06, 1.0135e-06 and 1.3327e-07
at 250, 500 and 1000 cells, observed order 2.845 / 2.927 / 2.886; extrapolated
to the spacing at which the carrier gate was anchored it reads 7.80e-07 and a
decade above it is 7.8e-06, so the rule of the anchoring memo applied to an
IONIZATION row returns the 1e-5 the stage rows already carry. **No tolerance
was set or changed**, and the reading at the gating radius (5.04e-06, which
would give a looser 5e-5) is reported and not adopted. The atomic gas was still
refused at this point, for three reasons read in the source: no caller (fourteen
gate sites in `src/EXHALE_main.f90` on `thereis_mol .and. carrier_transport`),
no frozen background, and no source outside the molecular network. The key-on
solve on `carrier_model_a_newton` is not certified and the refusing row is
`carrier balance He+`, 9.99e-01 at pass 1 to 7.74e-01 at pass 25 at cell 248,
while the control's H+ row falls 2.17e-01 to 5.71e-04. Read term by term the
He+ row is carried entirely by its chemistry (net source 68 times the stage
flux divergence at the gated cell and 5.0e4 times it in the shielded layer)
while the two builds' helium compositions agree to 0.20 per cent, and the same
operator on the control's own state reads 4.151e-02 instead of 7.736e-01: the
row is not broken, and even the locally solved composition is three and a half
decades above the gate. `hp_front` moved by 1.561e-01 in He III with the key
on, the elemental He/H ratio holding to 2.945e-12 across the column, so the
write-back moves charge and not nuclei.

### L36e: the atomic run carries the stages, and certifies

READ from the L36e report and `docs/lhs1140b_stationary_L36e_20260918.md`.
**On the certified atomic fiducial with `Ionization transport: True` the outer
iteration ACCEPTED at pass 56 and the state is CERTIFIED**, `info = 0`, the
three stage rows at 8.620e-07 (H+), 7.147e-07 (He+) and 2.360e-06 (He++)
against 1e-05, the two stage-sum entries at 2.161e-16 and 3.078e-16, the
elemental He/H partition at 8.679e-06 and every hydrodynamic row inside.
**Twenty-five passes are not enough, and the reason is a relaxation front and
not a refusing row**: over passes 1 to 25 the worst gated species row falls only
3.51e-02 to 1.76e-02 while the cell that binds it marches OUTWARD 217 to 491,
about 3.3 cells a pass, because the restart's local-equilibrium ionization is
being replaced by the transported partition and the front has to cross the
column; it reaches the outer boundary at pass 26 and the measure then falls
monotonically (1.27e-02 at 40, 2.03e-03 at 50, 8.68e-06 at 56). **Stage D is
closed on the negative side**: on the first state whose three ionization
fractions meet their own gate, the post-process moves the He 2^3S radial column
by 0.127 per cent, twenty times below the 5 per cent promotion rule, so the
level is not promoted to the carried set; the same column of the
local-equilibrium closure differs from its `_adv` by a factor 6.15. **Stage F
has its result**: the He I 10830 equivalent width from the certified
transported composition is 1.3678 %A (its own `_adv` 1.4053) against 20.3950
from the local-equilibrium solution and 1.5544 from that solution's `_adv`, so
the factor 13.12 is a property of the local-equilibrium closure and carrying
the stages removes it. Behind it, x(H II) reaches 0.0583 on the certified
transported state against 0.6657 on the local-equilibrium one, and
x(He II) + x(He III) 0.1888 against 0.9117. Six defects were found and fixed as
a separate item, among them an out-of-bounds write in the update-map diagnostic,
the gas particle density left at zero in every atomic run (`calc_ntot` called
only under `thereis_mol`, with `wasp_full` byte-identical across the fix), the
carrier row record printing a scale the acceptance does not divide by (a factor
101 at one binding cell), `iontrans` wrongly classified as a layout-changing
restart token, and **the He to H charge exchange missing from the stage row
sources**, which is the 1.2479e-04 key-on movement of `carrier_model_a_newton`.

### L37: the round trip of the molecular base mass row

READ from the L37 report and `docs/lhs1140b_stationary_L37_20260918.md`. **The
band**: on the six certified molecular states in hand the cell-1 hydrodynamic
mass row comes back from one no-step re-entry at 2.39 to 10.09 times the value
the run printed, and at 16.0 to 31.2 times that cell's own rounding floor, where
the tolerance is ten floors; on three certified atomic states of the same planet
and grid the ratio is 0.68 to 1.35 and the row stands at 0.21 to 1.52 floors,
which is L18's band met again. **It is a transient of the re-entry, not a
property of the state**: four of the six molecular states reach a reading that
repeats exactly by the third or fourth re-entry, and at that reading all six
stand at 0.08 to 7.5 floors, inside the anchor the code carries. **The mechanism
is the base ghost cells and the caloric equation of state**: the physical cells
round-trip exactly (`rho`, `v`, `p` of cells 1 and 2 to all 17 digits, the
composition of cell 1 to 1.7e-16), while the two inner ghost rows move by parts
in 1e-8 (ghost-0 H I, 4.835e-08) to 9e-6 (the He 2^3S column of ghost -1); the
molecular ghost is 99.998 per cent H2, so scaling its H2 column by `1 + eps`
moves the cell-1 row linearly with a slope of 0.30 to 0.33 of the row's own
scale, while the same perturbation up to 1e-06 moves an atomic base by nothing
at the printed precision, because there the pressure map carries no
composition. Exchanging the ghost composition rows between the archived and the
settled state reproduces the archived reading exactly, which is what identifies
the carrier. The composition refresh of the physical cells is innocent: it
moves the particle count by 1.295e-13, which at that slope is 4e-14 of the row,
five decades below the 6.9e-09 the re-entry moves it. The item proposes no
change of `c_round`: applying the anchoring rule to the round-trip
perturbation would give 100, but that perturbation is not one the arithmetic
forces and two further re-entries remove it, so anchoring on it would widen a
gate to admit a restart that does not reproduce its own base ghost. The three
ways out are in section 6.2 under D5. The item also made
`src/tests/physics_probe/molecular_energy_recipients.f90` read its three
representative cells from a state (`EXHALE_PROBE_STATE`): pointed at the
archived fiducial the distance from the compiled values is 0.00000E+00 and
pointed at the new solution it is 3.41525E-01, and pointed at an atomic state
or at a directory with no pair it stops rather than reading a zero.

### The bookkeeping items of the plan's section 11

READ from the bookkeeping report. (i) `src/modules/flux/low_mach_dissipation.f90`:
the nonnegative-dissipation claim is replaced by the indefinite-sign statement
with the measured eigenvalues of the symmetric part (+7.2e-17 on a uniform grid
at constant coefficient, -6.9e-03 with the gate closing at one face, -1.1e-10
on the LHS 1140 b 0.02-XUV state), the correctly signed form
`-M^-1 B^T K B` and its own measurements (-4.5e-17 and -5.9e-18), and the
statement that this module does not apply it; comment lines only, no executable
line moved. (ii) `src/tests/grid_and_gates/output_state_consistency.sh`: the row
`outer_iteration_ending_is_the_stagnation_one` read `pass_budget` against
`no_progress` on every control build of the previous series; the cause is that
the composition elimination sweep nudges the state every pass so the joint
distance keeps falling at full precision, and the ending is reached again at
`EXHALE_OUTER_PASSES=50`, where the first two consecutive passes without a fall
are 42 and 43 and pass 44 announces the refusal. The assertion was not touched;
the state handed back at that ending re-evaluates to itself at 6.4e-16 against
its 1e-12 allowance. Two other settings were measured and rejected and are
recorded in the row's comment. (iii) `docs/TO_BE_DONE.md`: the L7g cascade entry
keeps its blocking data request and gains the L31 reduced ladder as the bounded
substitute with the three cells' bracket.

### The rename

READ from the count_rename report. `global_parameters`' `integer :: count` is
now `integer :: marching_step` and `count_max` is `marching_step_max`. Of 541
occurrences of the token in `src/`, 95 are the variable (four of them inside
live `!$omp if(...)` clauses that a comment-skipping scan would have missed),
3 are the `system_clock` keyword, 1 is a genuine intrinsic call and 442 are
plain English in comments and printed strings; no local anywhere in `src/` is
named `count`. The printed labels are deliberately unchanged (`count=` in the
log, `count_max` in `EXHALE_setup.out`), because `src/tests/attempted_step/run.sh`
parses the first and `docs/input_schema.md` documents the second. RED and GREEN
on the hazard itself: a program that writes `use global_parameters` with no
`only` list and calls `count(mask)` fails to parse against the entry text and
compiles and prints `true elements = 3` against the delivered source. Eleven
source files changed; `wasp_full` (12251 steps), `mol_base_handoff` and
`hydrostatic_column` give data lines byte-identical to the control build, and
all 34 test suites give verdict lines byte-identical to the control, 3459 of
them, with the same six suites exiting nonzero on both for the entry text's own
reasons. Four comments that explained a loop by the shadowing were removed, the
loops left as they are.

### The update log

READ from the update-log report. One new section was appended to
`docs/Update_EXHALE_stage2.md`, `## 12. PLAN_20260917: the corrections that
follow the review of the 2026-09-17 handoff (2026-09-17 to 2026-09-18)`,
Markdown lines 11561 to 12516, with twelve subsections in plan order and 14
tables. The typeset twin was regenerated with
`python3 src/utils/update_log_to_tex.py` and `latexmk -pdf`, both exit 0:
`docs/Update_EXHALE_stage2.pdf` went from 218 to 235 pages. A second append
the same day (READ from the second update-log report) added nine subsections
for L36d, L36e, L37, L34b, L34c, L35, L30c, the `count` rename and the
further corrections, the file going from 12516 to 13614 lines and the PDF
to 257 pages, so every item of the series is in section 12. The closing
subsection stays PENDING until the advisor fills in the golden refresh: the
matrix at `REGRESSION_REL_TOL=0` on `3146d11b4090306dcea75bb9718edd22` moved
every one of the sixteen cases (its table is there), and the second matrix
on the final binary `74b96cdcf887dcfee301ec547f012d88` is the one the refresh
stands on.

---

## 2. The state of the tree

### 2.1 Source files changed in this series, and by which item

READ from the item reports.

| file | item |
|---|---|
| `src/modules/states/stationary_operator.f90` | L27 (NEW, 200 lines) |
| `src/EXHALE_main.f90` | L27 (+48 / -33, the operator selection and the evaluate route), L32 (the loop-top iterate record is in `steady_newton.f90`; here the outer loop's prints), L36e (+24 / -18, thirteen gate sites on the new predicate, one out-of-bounds allocation), L36c (a stale refusal reason, comment), the rename |
| `src/modules/time_step/certification.f90` | L27 (+59 / -3, the face-flux block), L36c (the stage-sum entry), L36d (one entry per element), L36e (+6 / -5, four gate sites) |
| `src/modules/time_step/steady_newton.f90` | L28 (+26 / -8, the `solve=` token), L32 (the loop-top-stop iterate record, print only), the rename (comment) |
| `src/modules/init/define_grid.f90` | L30 (+79 / -6, the shared geometry), L30b (+51 / -6, `spherical_cell_volume`) |
| `src/modules/functions/binary_element_diffusion.f90` | L30 (the shared geometry, `element_nucleus_face_flux`, the ghost rules), L30b (+116 / -58, the flux writer), L36 (the face-counter guard, `n_one`), L36c (+7, `element_nucleus_counts` exported) |
| `src/modules/lower_atmosphere/diffusive_photochemistry.f90` | L30 (+85, the shared geometry), L36c (+492 / -31, the stage row), L36d (the helium stages), L36e (+135 / -31, the predicate and the charge exchange), the rename |
| `src/modules/flux/species_face_flux.f90` | L30b (+9 / -2) |
| `src/modules/time_step/RK_rhs.f90` | L30b (+5 / -4, four sites) |
| `src/modules/states/Reconstruction.f90` | L30b (+10 / -5) |
| `src/modules/time_step/steady_residual.f90` | L30b (+4 / -5, three sites) |
| `src/modules/time_step/hydrodynamic_rows_body.inc` | L30b (+18 / -3, the factored form in the generic kind) |
| `src/modules/lower_atmosphere/molecular_reaction_heat.f90` | L31 (+151 / -23, the validity text and the startup report) |
| `src/modules/lower_atmosphere/h2_vibrational_relaxation.f90` | L31 (+58 / -22, comments) |
| `src/modules/lower_atmosphere/mol_rates.f90` | L31 (+57 / -25, comments) |
| `src/modules/functions/ionization_stage_transport.f90` | L36 (NEW, 518 lines), L36c (the halves a run carries) |
| `src/modules/nonlinear_system_solver/ion_cell_state.f90` | L36d (+22) |
| `src/modules/nonlinear_system_solver/ion_residual_core.f90` | L36d (+6 / -3, rows 2 and 3) |
| `src/modules/nonlinear_system_solver/constrained_chemical_equilibrium.f90` | L36d (+32 / -1), the rename |
| `src/modules/radiation/ionization_equilibrium.f90` | L36d (+43 / -18), L36e (+124 / -88, the predicate and `calc_ntot`), the rename |
| `src/modules/files_IO/input_read.f90` | L36c (+40 / -19), L36d, L36e (+42 / -33, the atomic refusal lifted), the rename |
| `src/modules/files_IO/write_setup_report.f90` | L36c (+11 / -3), L36d, L36e (+34 / -20), the rename |
| `src/modules/files_IO/load_IC.f90` | L36e (+21 / -13, `iontrans` reclassified) |
| `src/modules/flux/low_mach_dissipation.f90` | bookkeeping (+63 / -11, header only) |
| `src/modules/init/parameters.f90` | the rename (the declaration) |
| `src/modules/init/init.f90`, `src/modules/time_step/attempted_step.f90`, `src/modules/radiation/sed_read.f90` | the rename |
| `Makefile` | L27 (+1 in `SRC`), L30 (+1 in `DIFT_SRC`), L36 (+1 in `SRC`) |
| `src/utils/fortdep.py` | bookkeeping (scans `include` and emits it as a prerequisite) |
| `src/utils/map_state_to_grid.py` | bookkeeping (`cert_reason` dropped with `certified=` in a mapped seed) |
| `EXHALE_transit.py`, `exhale_transit_lib.py` | L36 (stage F, `EXHALE_TRANSIT_STATE`) |
| `LHS1140b/models/status.py` | L29 (one column), L34c (the states and the `superseded` reading) |
| `LHS1140b/models/write_reproduce.py` | bookkeeping (the seed and cold-start paragraphs) |

Test files added or changed: `src/tests/grid_and_gates/base_boundary_continuity_probe.f90`
(L27, rewritten on the new entry point), `stationary_evaluate_products.sh`
(L27, stage C), `run.sh` (L27, L36e), `output_state_consistency.sh`
(bookkeeping), `hydrostatic_residual.f90` (L30b),
`ionization_transport_atomic.sh` (L36e, NEW, 5 rows),
`restart_option_change.sh` (L36e, the route block re-pointed),
`base_closure_candidates_probe.f90` (L35, NEW, deliberately not registered);
`src/tests/steady_selfconsistent_residual/run.sh` and four new synthetic log
fixtures (L28); `src/tests/element_operator/element_operator_tests.f90` (+267,
L30); `src/tests/carrier_boundary_jacobian/` (+47, L30);
`src/tests/diffusion_tests.f90` (L30); `src/tests/species_face_flux/` (L30b);
`src/tests/physics_probe/molecular_energy_recipients.f90` (L31 blocks F, G, H,
62 to 81 assertions; L37 the state reader); `src/tests/h2_level_ladder/` (L31,
NEW, 32 assertions); `src/tests/coupled_block_jacobian/coupled_block_linear_system.f90`
and `run_linear_system.sh` (L32, NEW); `src/tests/ionization_stage_flux/`
(L36 15 to 25 rows, L36c to 27, L36d to 38);
`src/tests/ionization_imposed_fractions/` (L36d, 31 to 48);
`src/tests/carrier_reference_scales/`, `src/tests/steady_species_rows/`,
`src/tests/carrier_retry/`, `src/tests/carrier_returned_state_acceptance/`
(L36c, L36d, L36e).

### 2.2 New keys and their defaults

**No input key was added by this series**, and no default of an existing input
key moved. Two environment variables were added, both for a diagnostic or a
tool and both inert when unset:

| key | default | item |
|---|---|---|
| `EXHALE_PROBE_STATE` | unset; the compiled representative cells stand and no file is read | L37 (`src/tests/physics_probe/molecular_energy_recipients.f90`) |
| `EXHALE_TRANSIT_STATE` | unset, which is `adv`, the historical input pair | L36 (stage F, `EXHALE_transit.py`) |

What did change in the input contract, without a new key: `Ionization
transport: True` is now ACCEPTED in an atomic gas (L36e) and REFUSED with
`Coupled carrier solve: True` (L36c) and in a molecular gas whose carriers are
not transported; `iontrans` became a nameable `Restart option change:` token
rather than a layout-changing one (L36e); and the carried set behind the key is
now x(H II), x(He II) and x(He III), which `write_setup_report` echoes.

### 2.3 The binaries and the goldens

The catalog re-run was made on `LHS1140b/models/EXHALE_3146d11b.x`, md5
`3146d11b4090306dcea75bb9718edd22`, manifest
`LHS1140b/models/BINARY_MANIFEST_3146d11b4090.txt`; it carries L27, L28, L30,
L30b and L31, and it is the md5 every `REPRODUCE.md` written on 2026-09-18
names. The L29 reclassification was made on
`59bfdb3fc4d0104fc2e9c3734596d2f6`, the binary `docs/PLAN_20260917.md` stands
on. Items landing after L34 were measured each against a control build of its
own entry text, so the md5 chain of the series is
`59bfdb3fc4d0` to `3146d11b4090` (L34) to `d15a4b54fd5f` (L36) to
`91a0dca3daa3` (L36c) to `cf3d2941d1eb` (L36d) to `6a36529386102` (L36e) to
`74b96cdcf887` (the rename), each READ from the report that built it.

**FINAL BINARY: `74b96cdcf887dcfee301ec547f012d88`** (with the manifest
`LHS1140b/models/BINARY_MANIFEST_74b96cdcf887.txt`, which
`docs/PLAN_20260918_rev1.md` section 0.2 makes a precondition of every D-item
measurement).

**GOLDEN REFRESH: none in this series**, by the user's instruction of 2026-09-18 14:10 (the rev2 corrections continue on the same tree; the goldens stay the 2026-09-17 set; every item's control is its own entry-text build; one refresh when the rev2 items have landed). The strict matrix on `74b96cdcf887` was started and stopped on the same instruction; no golden was refreshed by any worker and
nothing was written under `backup/regression/`. What is known about the
movement, READ from the reports: the matrix at `REGRESSION_REL_TOL=0` on
`3146d11b` moved all sixteen cases (update-log report); `oxygen_chemistry`
moves by order unity for the reason L30c gives, and that movement is not
apportionable between the items of a series; `hp_zero_seed` and `hp_trace_seed`
were already stale against the entry text by 3.9e-03 and 3.3e-03 before L36c
touched them, `hp_front` is now stale by 1.6e-01 in the helium columns (L36d)
and `carrier_model_a_newton` by 1.2e-04 (L36e), and all four are key-on cases
whose refresh is one refresh. `run_fcheck.sh` was not re-run by any worker; the
bounds-checked build of L36c ran `hp_front` and a `He_diffusion` variant for 40
steps each with no runtime trap.

---

## 3. The catalog

### 3.1 What stands

READ from the L34a, L34b and L34c reports. 95 cases in `LHS1140b/models/`.

- **79 atomic states reproduced** their own claim in place on `3146d11b`
  through the L9 evaluate route: the 70 prescribed-composition cases that carry
  a state and the 9 last rungs of the flux-closure ladders. None refused, so
  none was re-solved.
- **The seven low-XUV cases are solved and certified**, each at outer pass 1
  with `Numerical flux: ROE`, seeded from the `.L25/` rung or direct solve its
  own `seed_from_ladder.txt` names. The catalog directories that had no
  `output/` state now have one.
- **Two molecular cases are certified**: the reference solution
  `molecular_scalar_gj1132_kzz1e9/HeH2.13` at outer pass 12 and
  `molecular_scalar_gj1132_kzz1e9/HeH9.7` at pass 10, both on `3146d11b`, both
  with the full acceptance list of section 1.
- **One certified atomic well-mixed seed was solved for the purpose**,
  `.L34/atomic_wm2.13`, because the catalog's atomic ladder has no He/H = 2.13
  rung and a molecular well-mixed case can be entered no other way.

### 3.2 The six molecular cases that did not certify, and where each stopped

READ from the L34b and L34c reports and from each case's own `not_solved.md`.

| case | seed | passes | where it stopped |
|---|---|---|---|
| `molecular_scalar_gj1132_kzz1e9/HeH0.55` | the certified reference at `--reservoir He/H 0.55` | 6 | carrier row 8.40e-01 to 7.09e-02 at pass 5 and back to 1.19e-01; hydrodynamic solve `info = 2` at every pass; stopped at the six-hour ceiling. (An earlier attempt from its own archived state gave a band 8.92e-03 to 3.11e-03 over six passes with the wind losing convergence from pass 5.) |
| `molecular_photochem_gj1132_kzzprofile/HeH9` | the certified wind of its atomic pair `.../HeH9/k04` | 1 at each of two pseudo-time starts | REFUSED: the element composition relaxation found no admissible advance and restored its entry composition, so every later pass is a copy of that one; mass row 2.848e-01 at cell 2, energy 3.167e-01 at cell 43 |
| `molecular_photochem_gj1132_kzzprofile/HeH2.09` | the certified wind of `.../HeH2.09/k05` | 12 | carrier row a band, 1.85e-01 to 4.42e-01, ending 3.07e-01; stopped at the ceiling |
| `molecular_scalar_gj1132_wellmixed/HeH2.13` | the atomic-to-molecular conversion of `.L34/atomic_wm2.13` | 7 | wind converged from pass 3, carrier row flat at 4.29e-01 from pass 5; stopped at the ceiling. **Started for the first time**, where `MODELS.md` recorded it as never started |
| `molecular_scalar_gj1132_wellmixed/HeH0.55` | the conversion of its certified atomic pair | 31 | wind converged from pass 5, carrier row to 3.76e-02 at pass 11, up to 1.56e-01 at 21, down to 6.48e-02 at 31; stopped at the ceiling |
| `molecular_scalar_gj1132_wellmixed/HeH0.083` | the conversion of its certified atomic pair | 40 + 40 | REFUSED on its own terms, both budgets spent: every hydrodynamic row inside tolerance (mass 2.005e-08 of 2.8e-08, momentum 5.509e-13, energy 3.110e-08 of 1.0e-06) and carrier balance H2 at 9.264e-02 of 1.0e-05 |

**The seeding rules the restart contract imposes**, which is the finding of
L34c and which the next campaign must plan around: a metal-free scalar
molecular state cannot seed a `photochem` case, because that group states C/H,
N/H and O/H at the matching level and no factor makes a carbon reservoir out of
a state that has none; and a state carrying HeH+ cannot be rescaled to a
different He/H by one factor, so a diffused molecular state cannot seed a
well-mixed run and a composition step inside the well-mixed family is forbidden
as well. The only entry to a molecular well-mixed case is therefore the
atomic-to-molecular conversion of an atomic well-mixed wind at that
composition. The measured bound on carrying a well-mixed molecular state at
all: the mapper accepts the 2026-09-16 `wellmixed/HeH0.55` state carried to
He/H = 1.0, and the mapped column then holds He/H between 0.999997254 and
1.000000000, a worst departure of 2.746e-06 against the loader's 1.0e-06
threshold.

### 3.3 `status.py`, `MODELS.md` and the memo

READ from the L29, L34a, L34b and L34c reports. `status.py` gained two things.
(i) A claim column: L29 added `claim on 59bfdb3f`, read from
`LHS1140b/models/CLAIMS_59bfdb3fc4d0.md` and from nothing else; L34c renamed it
`archived claim on 59bfdb3f` and made a row read `superseded` where the state
in the case directory has been written since that table was made, comparing the
`run=` stamp in the provenance header of the file the table read against the
table's own modification time (2026-09-17 23:27:50 KST). 14 of the 95 rows read
`superseded`: the seven low-XUV cases the L34 re-run solved and seven molecular
cases. (ii) The run states: `running` is no longer inferred from a missing
state file, so a case with no solver verdict prints `not solved` when a
`not_solved.md` is beside it, `not started` when neither is there, and
`marching-stop` unchanged for a state with no verdict; that removes the standing
false `running` on `molecular_scalar_gj1132_wellmixed/HeH2.13` and on
`molecular_photochem_gj1132_kzzprofile/HeH2.09`. `MODELS.md` sections 7 and 8
were regenerated with `python3 models/status.py --write`, and section 3 was
corrected by L34c to state the seeding rule the restart contract imposes and
what the three well-mixed cases did on it. The figures were regenerated with
`python3 LHS1140b/make_memo_figures.py` and `latexmk -pdf` built
`docs/lhs1140b_exhale_vs_pwinds.pdf` cleanly at **80 pages**, the same count
before and after L34c (77 at the previous handoff).

---

## 4. What is not going well, in detail

### 4.1 The molecular well-mixed carrier row

Three well-mixed cases were run to their ceilings on the current binary and
none certified, with the same shape in all three: **the stationary wind is
found and held and the H2 carrier balance alone refuses**, every relaxation of
every pass ending on the composition movement bound and none on its own
residual (READ, L34c). That is a sharper statement than the series began with,
because it removes the seed and the hydrodynamics as candidates. The best any
of them reached is 3.76e-02 on the gated row (`HeH0.55` at pass 11), against a
gate of 1.0e-05, and it rose again afterwards. L33's classification of the
earlier runs stands beside it: two of the three are front relaxation with the
refusing cell walking outward, the third is throttled by a movement bound that
only ever ratchets down within a run, and the handover route the design
provides has never fired because the consecutive-no-progress count never
reached three.

### 4.2 The two `molecular_photochem` cases

`HeH9` is REFUSED rather than slow: at both pseudo-time starts the element
composition relaxation found no admissible advance and restored its entry
composition, so the outer loop ends by construction and more passes would be
copies of one (READ, L34c). `HeH2.09` spends its passes in a band, 1.85e-01 to
4.42e-01 over twelve passes, and was stopped at the ceiling; its earlier run on
the predecessor binary needed 21 passes at 1614 to 11891 s each with the gated
row RISING over the last four (READ, L34b). Neither can be seeded from a scalar
molecular state, for the reservoir reason of section 3.2, so their seed is the
certified wind of their own atomic closure rung.

### 4.3 The molecular base re-entry

The cell-1 hydrodynamic mass row of a certified molecular state is not
reproduced by a no-step re-entry of the state the run wrote: the band is 2.39 to
10.09 of the run's own value and 16.0 to 31.2 of the cell's rounding floor,
against 0.68 to 1.35 and 0.21 to 1.52 on atomic states of the same planet and
grid (READ, L37). It settles: four of the six molecular states repeat exactly by
the third or fourth re-entry and all six are then inside the anchor. The carrier
is the pair of inner ghost rows, which the loader re-derives and projects and
which the file does not carry to the precision this row needs; the molecular
ghost is 99.998 per cent H2, so the caloric equation of state makes the base
face flux a function of its H2 count with a slope of 0.30 to 0.33, and a ghost
H2 count off by 1e-08 is four rounding floors of the row that gates the state.
Two further readings belong with it: two certifications inside ONE run, of one
conserved state, already differ by a factor 4.3 at the base cells while agreeing
to four digits at the verdict cell; and `EXHALE_RESIDUAL=1` is not a second
reading of the same operation, because its `reload_equilibrate` defaults on, so
L34b's "third reading" is one settling step further along the same ladder. This
is item D5 of the next plan.

### 4.4 The He+ row scale

On the hot-Uranus fixture `carrier_model_a_newton` the transported
`carrier balance He+` row stands at 9.99e-01 to 7.74e-01 over 25 passes while
the two builds' helium compositions agree to 0.20 per cent, and the same
operator evaluating the same row on the control's own locally equilibrated state
reads 4.151e-02, still three and a half decades above the 1e-05 gate (READ,
L36d). The reason is the denominator: the certification's carrier scale is
`row_terms_phys`, which contains the absolute NET chemical source, and for a
stage whose production and loss cancel to three decades the net is the same size
as the residual. On the atomic fiducial the scale does not bite at all, the same
row certifying four decades under its gate (READ, L36e), so this is a property
of a regime and not of the row. Measuring a fast stage against its production
PLUS its loss needs the split at the one place the terms are written
(`mol_heh_rows` rows 2 and 3, `heh_tr_rows` all four) and is a change to a
shared acceptance interface with its own anchoring. This is item D3.

### 4.5 The 25-pass front of transported ionization

A campaign that judges a transported-ionization configuration at 25 outer passes
reads a front and not a solution: the restart's local-equilibrium partition is
replaced from the inside outward at about 3.3 cells a pass, so on a 500-cell
domain no gated row can begin to fall before roughly pass 30, and the fiducial
accepted at pass 56 (READ, L36e). Over passes 1 to 25 the worst gated row fell
only from 3.51e-02 to 1.76e-02 while its binding cell marched 217 to 491. Any
later plan that budgets passes for this configuration has to budget for the
crossing first.

### 4.6 The continuation can overwrite a better state

`run_case.sh` continues at a second pseudo-time start whenever the first solve
does not return `info = 0`, including when the hydrodynamic rows are already
inside their tolerances and only the composition refuses. MEASURED by L34c on
`wellmixed/HeH0.083`: the first solve's forty passes reached a gated carrier row
of 4.13e-02 and the continuation's forty took it back to 7.37e-02, **and the
state on disk is the continuation's.** The composition movement bound is part of
it: it never recovers within a run (it only halves), and the continuation starts
it at 1.0e-02 again, where the row rises, so simply restoring the bound is not
by itself the remedy.

### 4.7 Two case directories hold a state pair from two runs

A solve writes `Hydro_ioniz.txt` only when it ends, so a run stopped inside one
leaves the previous run's state beside its own `*_IC.txt` seed. READ from L34c:
`molecular_scalar_gj1132_kzz1e9/HeH0.55` and
`molecular_scalar_gj1132_wellmixed/HeH0.55` are in that state, and both also
carry a `run_dtau0_first.log` dated 2026-09-16 beside a `run.log` of 2026-09-18,
because those runs never reached the branch that renames the log. Anything read
out of either directory has to name which file it came from.

---

## 5. Deferred items, with their reasons

- **The L7g cascade (the Bossion tables).** Unchanged and still blocking: the
  H2-H state-to-state rates above v = 4 and the state-resolved three-body
  collisional dissociation for M = H are not distributed with their papers. L31
  built the bounded substitute (`src/tests/h2_level_ladder/`) and measured the
  scalar form to keep more energy out of the gas than any member of its
  bracket, so the default stands and the absence does not matter in this layer;
  the entry in `docs/TO_BE_DONE.md` now says so.
- **The element operator's own derivative probe.** L30 delivered the
  conservation, divergence and order rows but not a direct probe of the
  element tridiagonal Jacobian, which is private and not reachable from
  outside; what stands in its place is indirect (the Newton of
  `solve_mass_fraction` still converging on every column of `diffusion_tests`,
  largest relative residual 2.496e-09, and the fixed-wind relaxation still
  reaching its fixed point). The carrier DIFFUSIVE weight likewise has no unit
  row, because the public route to `carrier_residual` needs a full molecular
  state no light driver builds; it is established by reading the two sites and
  by the regression movement.
- **Metal charge exchange in the stage sources (D7).** `cx_add_to_fvec` writes
  group A into row 1 and groups C and D into rows 2 and 3 of every metals
  ionization system, and `carrier_source` does not call it, so a run with
  metals and `Ionization transport: True` would carry stage sources that differ
  from the sweep's rows by those terms. Adding them needs the frozen metal
  densities of each stage and the `cx_metal_base` layout, which the transport
  operator does not hold. The fiducial L36e certified is metal-free. The plan
  makes the interim refusal and the addition item D7.
- **The production and loss split of a fast ionization stage.** Section 4.4;
  deferred because it changes a shared acceptance interface and needs its own
  anchoring. Item D3.
- **`src/tests/carrier_retry/`'s knife-edge row.**
  **CLOSED 2026-09-18 by D4a and D4b** (`docs/lhs1140b_stationary_D4a_20260918.md`, `_D4b_`): the row compared an abundance movement of 1.1e-6 against a residual tolerance that bounds it only at 1.8e+6 (HeH+ row sensitivity `|R|/eps` 5.5e-13); the returned composition IS a root of its network at the returned temperature (residual 2.6e-14). The row is replaced by `every_species_row_of_the_returned_state_is_a_root_at_the_returned_temperature` (5.3e-4 on the entry text, 4.4e-16 after), the abundance movement stays as a non-gating line with its conditioning, and the refusal branch of `relax_photochemical_composition` now restores the rate state, the caloric maps, `dp_bc`, `n_part_cell1` and the carrier checkpoint (nine differences, each exactly zero after). No tolerance moved.
  `one_further_sweep_leaves_the_returned_composition_within_the_sweep_tolerance`
  reads 1.100521341448987e-06 against `ieq_res_tol` = 1e-6 on the measured and
  the control build alike, bound by HeH+ at 1.3e-12 of the gas in the base cell;
  L36d showed the chain H2+ to H3+ to HeH+ amplifies by a factor 7 to 10 a link,
  that H3+ carries 53 per cent of that cell's electrons, that a third sweep moves
  the composition by 1.0085e-10 (so it settles and is not noise), and that
  without the imposition the same sweep moves it by 1.5067e+08 (so the
  imposition holds the state rather than breaking it). No tolerance and no
  presence threshold was touched. Item D4.
- **The helium held in HeH+ is not reserved from the helium ion-stage simplex
  (D6).** The projection imposes `x_HeII + x_HeIII <= 1` while the fractions are
  normalized by the TOTAL helium nucleus count, which includes HeH+, and the
  write-back forms neutral helium as the remainder minus HeH+; the review's
  executed counterexample (HeH+ 0.10, x(He II) 0.60, x(He III) 0.35) passes the
  simplex at 0.95 and leaves an inventory of 1.05. Reachability in a production
  run was not tested; the defect is in the expressions, and it gates D3.
- **`coupled_block_jacobian`'s two rows.** `the_species_rows_scale_is_the_
  certifications` reads `2_of_13` and then `3_of_13`, and
  `the_assembled_action_reproduces_the_central_difference_carrier` 8.95e-01
  against 1e-2, both on the CONTROL build as well; they are about the coupled
  route, which the ionization key now refuses, and L32's verdict says what that
  route's action is. Not fixed by any item.
- **`residual_determinism`'s row
  `closure_spread_within_the_row_tolerance_atomic_elem_newton`** reads
  1.715e+04 on the measured and the control build alike (L36d). Pre-existing.
- **Small documentation and tooling debts.** `docs/input_schema.md` line 171
  gives `count_max` in both its default and its variable column, and the
  variable is now `marching_step_max`; the printed labels `count=` and
  `count_max` no longer match the variable names, and bringing them into line
  means editing `src/tests/attempted_step/run.sh` and the schema in one change;
  three loops in `EXHALE_main.f90`, `steady_newton.f90` and `sed_read.f90` exist
  only because the name was shadowed and can now be written with the intrinsic
  (their false explanations were removed, the loops left);
  `src/tests/physics_probe/README.md` documents neither
  `molecular_energy_recipients.f90` nor `molecular_third_body_in_the_chemistry.f90`
  and does not list `EXHALE_PROBE_STATE`;
  `src/tests/grid_and_gates/README.md` section 4 is stale on three counts, one of
  them a temperature round trip quoted at 6.5e-14 where the bookkeeping worker
  measures 6.4e-16.
- **Two fixtures that are not what their names suggest.**
  `backup/regression/hydrostatic_column` is a 300-step relaxation snapshot whose
  first cell carries an inflow at Mach 0.2257, so a future well-balancedness
  test must build its own state at rest; `backup/regression/oxygen_chemistry` has
  no `NOTE.txt` recording that its comparison against a golden is a
  reproducibility statement and nothing else, which two bitwise experiments now
  establish. Changing either is a golden refresh, so neither was changed.
- **`src/tests/grid_and_gates/base_closure_candidates_probe.f90` is not
  registered** in that suite's `run.sh`; the L35 memo states the one block that
  would register it behind an environment variable.

---

## 6. Decisions

### 6.1 Decisions the user made this series

| date | decision |
|---|---|
| 2026-09-18 | Proceed with `docs/PLAN_20260917.md` as written, item by item, in its stated order and dependencies |
| 2026-09-18 | The review of the first 2026-09-17 handoff is folded into `docs/session_handoff_20260917_rev1.md` and its corrections are executed as plan items rather than by editing the memos of the previous series; the memos of that series carry head notes only (L26 sections R3 and R9, L22 step 3) |
| 2026-09-18 | The goldens are refreshed ONCE, at the close of the series, by the advisor, after L30c's diagnosis of the `oxygen_chemistry` movement; no worker refreshes one |
| 2026-09-18 | D2a: **A** at the stationary contact (the reservoir owns the level's thermodynamic state, including at rest), the rest of the Codex D2a text kept; D8: the Codex generation structure as reviewed in `docs/DECISION_D2a_D8_review.md`, with legacy import and D8a first; the golden refresh deferred and the strict matrix stopped; the rev2 items to proceed in the order of `docs/PLAN_20260918_rev2.md` section 0.1 |

### 6.2 Decisions still open, each with a recommendation

From `docs/PLAN_20260918_rev1.md`'s decision table, with the evidence each
needs before it is taken.

| open decision | recommendation | what must exist first |
|---|---|---|
| **D1b**, the historical grids and the new default base cell width | Pin the resolved old grid inputs in every reproduction record and adopt `2.0d-4` for NEW configurations only; no archive is mapped for this reason | D1a: one named default, a provenance flag set by the reader, `write_resolved_config` recording the width at round-trip precision, and coordinate equality and restart compatibility measured on the affected `Mixed` grids |
| **D2a**, the boundary model at zero and reversed flow | State the contact and reservoir model FIRST; no blending formula is chosen before it, and the first text's premise (that no characteristic leaves at rest) is withdrawn, since at rest the speeds are `-c_s, 0, +c_s` | the characteristic count, a manufactured DISCRETE hydrostatic state, a stationary pressure-balanced contact with differing entropy, and weak-flow directional derivatives around zero |
| **D3b**, the He II acceptance norm | Keep the gate while D3a measures; add complementary criteria only if D3a shows the need, and compare DIMENSIONAL residuals (`tolerance x scale`) and never normalized numbers | D3a: gross production and loss by channel from the production chemistry, the balances at the worst cells, the cancellation floor of `P - L`, the conditioning, stiff manufactured cases, and the structural tests. D3a depends on D6 |
| **D4b**, the closure contract | A directly evaluated returned-state contract (a chemical residual or increment criterion at the returned temperature), with settling sweeps kept as a SEPARATELY named diagnostic that never overwrites the returned-state assertion | D4a: the chemical and thermal residuals at the returned state, one complete fixed-conserved-energy closure increment, the rollback check, and the conditioning that relates `ieq_res_tol` to an abundance movement |
| **D5b**, the boundary-state operation | Derive the ghosts from the declared physical state and inputs in ONE operation with an explicit contract, and serialize only genuinely independent reservoir variables; do not widen `c_round` to cover a restart that does not reproduce its own base ghost | D5a: the boundary at each stage of the evaluate route, serialization measured separately from closure, boundary independence and idempotence, and same-state operator agreement as a signed dimensional difference and not as a ratio |
| **catalog replacement after D5b** | Validate on copies and re-solve where the corrected boundary moves the root; keep the three cases apart (reproduction under the same operator, measurement under the corrected one, a new solve where the root moves) | same-operator against changed-operator provenance and a reviewed result, never an old header. A corrected operator may refuse a formerly certified state, and that is reported, not absorbed |

**D6 and D7 are items to execute, not decisions.** D6 is a conservation repair:
the constraint becomes `x_HeII + x_HeIII <= 1 - n_HeH+/n_He,nuc` applied
identically in trial admissibility, finite-difference perturbations, source
evaluation, projection and write-back, with no clipping of a negative neutral
remainder while the positive species stand; its tests measure the actual
returned species (total H and He nuclei, mass, charge) with nonzero HeH+ and
nearly exhausted neutral helium, RED on the counterexample state and GREEN
after, key-off byte identity on the four cases of L36e. D7 adds the metal
charge-exchange terms to the stage sources through the same routine that writes
them into the sweep's rows (one spelling), and until it lands `input_read.f90`
must REFUSE `Ionization transport: True` with metals on, the message naming the
item.

---

## 7. Paths

**Memos of this series** (all under `docs/`):
`lhs1140b_stationary_L30c_20260918.md`,
`lhs1140b_stationary_L31_energy_cycles_20260917.md`,
`lhs1140b_stationary_L32_20260917.md`,
`lhs1140b_stationary_L33_20260917.md`,
`lhs1140b_stationary_L34a_20260918.md`,
`lhs1140b_stationary_L34b_20260918.md`,
`lhs1140b_stationary_L34c_20260918.md`,
`lhs1140b_stationary_L35_20260918.md`,
`lhs1140b_stationary_L36_20260918.md`,
`lhs1140b_stationary_L36c_20260918.md`,
`lhs1140b_stationary_L36d_20260918.md`,
`lhs1140b_stationary_L36e_20260918.md`,
`lhs1140b_stationary_L37_20260918.md`.
**L27, L28, L29, L30 and L30b wrote no memo**: their record is the artifact
that carries the statement, namely the head notes at
`lhs1140b_stationary_L26_20260916.md` sections R3 and R9 and at
`lhs1140b_stationary_L22_20260916.md` step 3, the table
`LHS1140b/models/CLAIMS_59bfdb3fc4d0.md`, and the corrected source headers and
new suite rows themselves.

**Plans and reviews**: `docs/PLAN_20260917.md` (the plan this series executed);
`docs/PLAN_20260918.md` (frozen, SHA-256
`be4184564f966bd5c38cad42c5617571d953768c7a1d2a86d1057c1444e16e69`),
`docs/PLAN_20260918_review.md` (the Codex review of it, with its retained
arithmetic checks `docs/audit_20260905/plan_20260918_review_checks.{py,log}`)
and `docs/PLAN_20260918_rev1.md` (the plan that follows this handoff).
The previous plan and its three reviews are `docs/PLAN_20260916_rev3.md` and
`docs/PLAN_20260916_review{,1,2}.md`; the item rows are in
`docs/PLAN_20260913_lhs_stationary.md`. The handoff this series started from is
`docs/session_handoff_20260917_rev1.md`, with `docs/session_handoff_20260917.md`
kept as the text its review read and `docs/session_handoff_20260917_review.md`
the review.

**Audit directories**: `docs/audit_20260905/L32_20260917/` (the frozen
checkpoint of item L32, its `identities.txt` and every run log; the checkpoint
files' md5s are in the L32 report and in the memo),
`docs/audit_20260905/L15_stage2_20260916/`, `L22_step2_20260916/`,
`L22_step2b_20260916/`, `L22_step3_20260916/`, and the review's retained
evidence `docs/audit_20260905/handoff_20260917_*`.

**Run directories** (under `LHS1140b/models/`):
- `.L34/`: `atomic_wm2.13` (the certified atomic well-mixed seed solved for the
  purpose), `rung1` (rung 1 of the L33 ladder, which did not certify and so
  ended it), the created and unrun rungs 2 to 4, and `rung_wm2.13` (the refused
  load that established the well-mixed seeding rule).
- `.L35/`: the six matched-domain solves `g{500,1000,2000}_{hllc,roe}`, three
  grid probes, `EXHALE_L35.x`, `launch_all.sh`, `solve_and_post.sh`,
  `grid_replica.py`, `solve_grid.py`, `collect.py`, `table.py`, `profiles.py`,
  `make_figure.py`; each stopped 40-pass attempt is kept as
  `<case>/run_passes40_stopped.log`.
- `.L22/`: `i3_alt/` (the alternation run that certified), `i3_block/` (the
  coupled block, the source of L32's checkpoint), `i4_wm0083/`, `i4_wm055/`,
  `i4_kz0083/` (the three runs L33 read), and the reader L33 wrote,
  `outer_pass_history.py`.
- `.L25/`: the seven direct Roe solves, the grid study `grid1000_probe/`, and
  `ladder/L{1,2,3,4}/x<tag>/`, the eighteen certified rungs the seven catalog
  cases were seeded from.
- `.L26/`: `fid_resolve/` (the certified fiducial, the state L27, L35 and L36e
  measured on), `x002_hllc/` and `x002_roe/`.
- `.L14/x003_HeH2.13`, `.L14/x002_HeH2.13`: the frozen states L29 and the
  review's probes read, never written to.
- Each catalog case touched by L34 carries its previous products in
  `output_pre_L34/` beside it, created once and never overwritten, and the three
  cases whose `output/` held L34b's stopped state also carry
  `output_L34b_stopped/`.

**The L32 checkpoint**: `docs/audit_20260905/L32_20260917/checkpoint/`
(`input.inp`, `base.inp`, `output/Hydro_ioniz_IC.txt`,
`output/Ion_species_IC.txt`), the entry state of the coupled block on the
certified `molecular_scalar_gj1132_kzz1e9/HeH2.13` case, with the exact command
of every diagnostic run recorded beside it.

**Figures added this series**: `docs/figures/lhs1140b_L34b_reference.png` (T,
rho, v, x(H2), x(H3+) and the He I 10830 line, new against seed, with an inset
on the base layer) and `docs/figures/lhs1140b_L35_flux_grid.png`.

**Reference PDFs used** (all already in
`/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/references/`; none was added this
series): for L31, `Takagi_2002_PhST_96_52.pdf`,
`Giusti-Suzor_1983PRA_28_682.pdf`, `Guberman_1994PRA_49_R4277.pdf`,
`Kokoouline_2001Nature_412_891.pdf`, `Strasser_2001PhysRevLett_86_779.pdf`,
`Cohen_1983JPCRD_12_531.pdf`, `Paolini_2011PRA_83_042713.pdf`,
`Lique_2015_MNRAS_453_810.pdf` with `Lique_2015MNRAS_453_810_data/`,
`Jozwiak_2024A&A_685_A113.pdf` with `Jozwiak_2024J_A+A_685_A113/`,
`LeBourlot_1999_MNRAS_305_802.pdf`, `Wolniewicz_1998_ApJS_115_293.pdf`, and
`Bossion_2018_MNRAS_480_3718.pdf` with `Nesterenok_2019MNRAS_489_4520.pdf` for
the deferred cascade; for L36e, `Huang_2023_ApJ_951_123.pdf` Table 4 (the He to
H charge exchange); for the L35 and D2 reading of the boundary,
`Thompson_1987JCP_68_1.pdf`; the Rieper papers of L25
(`Rieper_2009JCP_228_2918.pdf`, `Rieper_2010JCP_229_221.pdf`,
`Rieper_2011JCP_230_5263.pdf`) are cited by the low-Mach header the bookkeeping
item corrected.

---

## 8. Terminology

`LHS1140b/MODELS.md` section 2 is the naming contract. The four words of the
previous handoff keep their meanings, and four more were used constantly enough
this series to be fixed here.

- **rung**: one step of a continuation ladder, and nothing else. A composition
  rung is a step in reservoir He/H, an XUV rung a step in the spectrum scaling;
  a rung is not a case and not a pass, and a rung's own directory carries a full
  `input.inp` and `output/`.
- **photochem**: the lower-boundary field of a group name, meaning the Photochem
  column handed over as `Lower atmosphere profile:`. It is a boundary treatment,
  NOT the code's molecular photochemistry module and NOT
  `diffusive_photochemistry.f90`.
- **the XUV grid**: the set of cases whose spectrum is the GJ 1132 proxy scaled
  by a factor. It is a grid in the spectrum scaling at fixed composition; when
  the number of cells is meant, say "the 500-cell grid".
- **pass**: one outer pass of the stationary alternation (the hydrodynamic solve
  followed by the composition relaxation). Never a Newton iteration, never a
  Krylov product, never a marching step.
- **stage transport**: the transport of the ionization stages x(H II), x(He II)
  and x(He III) as fractions per element NUCLEUS, on the element operator's own
  nucleus face flux, in `ionization_stage_transport.f90`. It is not the
  molecular carrier transport (which carries species mass fractions on the bulk
  mass flux) and it is not the elemental diffusion (which carries the elements
  against each other). The key is `Ionization transport:`.
- **transported row**: a row of the transport-chemistry operator whose unknown
  is carried by that operator rather than eliminated by the local sweep. A
  transported row is not a certified row: `carrier balance He+` is a transported
  row that stands three decades above its gate on the hot-Uranus fixture and
  four decades under it on the atomic fiducial.
- **re-entry**: one no-step evaluation of a state through
  `Restart intent: stationary evaluate`, handed the state a previous evaluation
  wrote. E1 is the first, E4 the fourth; a reading that repeats under re-entry
  is SETTLED. A re-entry is not a re-solve, and a settled reading is not a
  certification.
- **the archived claim**: the `(certified, cert_reason)` pair in the header of a
  stored state, which is a claim made by the executable that wrote it and is
  never a current acceptance. `reproduced` and `refused` are the two verdicts
  the evaluate route gives it under a named binary and operator; `superseded`
  says the state in the directory has been written since the claim was
  tabulated.

---

## 9. Request list for external data

One item was raised to the head of the list by measurement this series; the
rest are unchanged.

1. **A published helium collider coefficient for the three-body association
   `H + H + He -> H2 + He`, over the 200 to 1600 K the molecular layer spans
   (200 to 5000 K would settle the extrapolation too). RAISED TO THE HEAD BY
   L31.** The code uses ARGON's coefficient (Cohen and Westberg 1983) as an
   argon-based estimate, with that paper's own +/-0.3 in log as the stated
   uncertainty and the distance from argon to helium explicitly outside it.
   READ from the L31 report: at the certified molecular base the helium
   efficiency moves the chemical heat by **+41 / -21 per cent for +/-0.3 dex**,
   four to six orders above every other approximation of the layer, and L34b
   measured the same bracket on the solved state as +42.0 / -21.1 per cent,
   **eight times the change the whole L7g correction made in that heat.** So the
   base heat the re-solved molecular catalog quotes carries this one uncertainty
   and nothing else of comparable size. Paolini, Ohlinger and Forrey (2011),
   Phys. Rev. A 83, 042713 plot a helium total over 0 to 350 K and tabulate
   nothing; a controlled extraction of that curve, or a helium measurement or
   calculation over the layer's range, would turn the estimate into a value with
   a stated range.
2. **BLOCKING for the L7g cascade.** The H2-H state-to-state rate tables of
   Bossion, Scribano, Lique and Parlant (2018), MNRAS 480, 3718, as extended to
   all bound states by Nesterenok et al. (2019), MNRAS 489, 4520, together with
   the state-resolved three-body collisional dissociation for M = H. The paper
   carries no supporting information, no repository and no data-availability
   statement. The named authors to ask are Nesterenok, Bossion, Scribano and
   Lique. Without them the code keeps the bounded `q(R15) f_quench` form, which
   L31 measured valid in this layer (the scalar radiates 2.8 to 13.2 times what
   the most radiative member of the reduced ladder does, so it never overstates
   the heat) and which would not be valid in a shallower or cooler one.
3. **Moses and Bass (2000), JGR 105, 7013**, still outstanding from the previous
   two handoffs: the source of the R17 Arrhenius branch above about 700 K, where
   nothing is measured.
4. **The origin of Cloudy's `coll_rates_He_ORNL.dat`**, cited in Cloudy only as
   "Lee et al. 2007, ApJ, in preparation", which ADS does not carry. Three
   plausible published leads were identified and none confirmed
   (`docs/cloudy_h2_model_reference_20260916.md`).
