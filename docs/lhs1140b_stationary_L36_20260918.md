# LHS 1140 b, item L36: the ionization-stage flux built on the shared element flux

Item L36 of `docs/PLAN_20260917.md` section 10, which carries stages B to F
of `docs/lhs1140b_stationary_L12a_design_20260913.md` as
`docs/lhs1140b_stationary_L12b_derivation_20260916.md` sequenced them.  It
follows L30 and L30b, which exposed `element_nucleus_face_flux` and put the
element and carrier transport on one spherical geometry.

Every number is MEASURED (run or computed here) unless marked READ.

---

## 0. Verdict, before the detail

1. **Stage B is built and tested.**  The stage face flux, its divergence on
   the L30 geometry, its tridiagonal Jacobian rows and its simplex
   projection are one new module,
   `src/modules/functions/ionization_stage_transport.f90`, written on the
   element operator's own `element_nucleus_face_flux`.  The identity
   `sum_k F_k(f) = N_el(f)` holds at every face to 2.0e-16 of the terms,
   at the boundary faces to 1.1e-16, after the simplex projection to
   2.2e-16; the stage divergences sum to the element nucleus divergence to
   4.5e-16 and the column sum of `V_j sum_k D_k(j)` is the difference of
   the two boundary nucleus fluxes to 1.2e-17 (all MEASURED, 500-cell
   production Mixed grid).  The analytic rows match a central difference of
   the operator they differentiate to 1.8e-11.

2. **Stages C and E cannot be built in this item and are not.**  The rows
   are solved inside the carrier operator
   (`src/modules/lower_atmosphere/diffusive_photochemistry.f90`) and the
   configuration is gated in `src/modules/files_IO/input_read.f90`; neither
   file belongs to this item.  The new module is therefore the ONE spelling
   of the stage flux, complete and tested, with no caller in a run yet.
   Section 6 states exactly what the owner of those two files has to do.

3. **The one-velocity approximation has a measured validity range, and it is
   narrow.**  MEASURED on the certified atomic fiducial
   `LHS1140b/models/.L26/fid_resolve`: the ambipolar drift of an ion against
   its own neutral stays below 0.1 of the bulk velocity only inside
   1.22 R_p (hydrogen) and 1.24 R_p (helium), and over `r >= 1.2 R_p` it
   reaches 0.62 and 0.43 without resonant charge exchange and 0.52 and 0.34
   with it.  The momentum transfer itself is fast (`t_coll/t_flow < 0.1`
   inside 14.5 R_p for hydrogen and 8.2 R_p for helium), so the stages are
   collisionally locked in time scale and still drift: the drift is set by
   the ambipolar field against the friction, not by the friction alone.
   Section 3.  This is written at the site, in the module header.

4. **The stage transport REPLACES the existing proton carrier; it is not a
   second transport.**  Section 1.

5. **Key-off byte identity: MEASURED on `wasp_full`, `mol_diffusion` and
   `hp_front`** (section 5).  No golden was refreshed and `EXHALE.x` was
   never touched.

6. **Stage F is built on the transit side and reports a factor 13.1.**  The
   tool's input pair is selectable (`EXHALE_TRANSIT_STATE`), the pair is
   named in every product header, and the measured distance between the two
   states travels with it.  MEASURED on the fiducial: the He I 10830
   equivalent width is 20.3950 percent A from the solved state against
   1.5544 percent A from the `_adv` profile, a factor 13.12 (L10 reported
   13.3 on 15.5116 against 1.1702 percent A for a different model, READ).
   Section 4.

7. **Stage D's promotion rule cannot be evaluated yet, and the half of it
   that can be measured was.**  Its precondition is a state whose three
   ionization fractions already agree to their own gate, which needs stage
   C.  The systematic it would be judged against is MEASURED: the He 2^3S
   radial column of the fiducial differs between the two closures by a
   factor 6.15 (4.9292e11 against 8.0114e10 cm^-2).  Section 4.3.

---

## 1. What this is, against the proton carrier that already exists

The key is unchanged: `Ionization transport: True`
(`input_read.f90` line 667, `ionization_transport` in `parameters.f90` line
464, READ).  Today it declares ONE extra carrier of the molecular carrier
operator, `ic_Hp` (`diffusive_photochemistry.f90` lines 1385 to 1388, READ),
and refuses any configuration without `Molecular chemistry: True` and
`Molecular carrier transport: True` (`input_read.f90` lines 2185 to 2199,
READ).

That proton carrier and the hydrogen row of this module are the SAME
physical object written in two variables:

| | the proton carrier today | the stage row |
|---|---|---|
| variable | `f_sp(:,isp_HII)`, a fraction per unit MASS | `x(H II)`, a fraction per hydrogen NUCLEUS |
| advected on | the BULK face mass flux | the element NUCLEUS face flux |
| eddy term | the whole mixing ratio flux `-n_tot K d(n_HII/n_tot)/dr` | `-n_el K dx/dr` alone |
| molecular diffusion | zero (`carrier_diffusivities`, the `ic_Hp` branch) | not carried |
| helium | not carried | `x(He II)`, `x(He III)` |

The present spelling is correct exactly while the bulk mass flux carries no
eddy term of its own, which is true today; it stops being correct the moment
the advective half becomes the element flux, which does carry one (L12b
section 1, last paragraph, READ).  **So this module is the replacement.**
The intended end state is one transport: `ic_Hp` and the two helium rows
formed from equation (1) of the module header, and no second copy of the
stage flux anywhere.  Until the carrier operator is re-pointed at it
(section 6), the tree holds the new spelling with no caller and the old
spelling behind the key, which is a state to leave as soon as possible and
not a design.

The He 2^3S level is NOT a stage of this partition: it is a sublevel inside
He I (`species_table.f90` `bsp_is_excited_level`, READ) and stays with the
sweep, for the Damkohler reason of L12a section 3.

---

## 2. What was built

### 2.1 The two L30b notes on `element_nucleus_face_flux`

Both were fixed in place (L30b report section 7, items 3 and 6, READ):

- **The face counters are restored inside the routine.**  The advective half
  calls `species_face_fraction`, which increments `n_species_faces_bounded`
  and raises `species_face_excursion`; those count the faces the RUN had to
  bound, so an evaluation of the flux of a state already in hand must not
  add to them.  The save and restore now sit in
  `element_nucleus_face_flux` itself, where every caller inherits them,
  instead of at each call site.
- **The units of the advective output are stated correctly.**  The header
  said CODE units; the routine multiplies the caller's `Frho` by a
  dimensionless face mass fraction and therefore returns whatever units
  `Frho` was in.  The header now says so, with both cases named (the code
  face mass flux of the mass row, and the g cm^-2 s^-1 flux the element
  flux profile writer passes).

One output was added: the optional `n_one(0:N)`, the HYDROGEN nucleus face
density, on the same arithmetic face rule as the helium one the routine
already returned.  A row written per hydrogen nucleus needs it, and
rebuilding it in the caller would be a second copy of a face rule.  It is
optional, so no existing caller changes.

### 2.2 `src/modules/functions/ionization_stage_transport.f90`

Public:

| routine | what it returns |
|---|---|
| `ionization_stage_face_fractions` | the face fractions, closing the simplex AT THE FACE: the carried stages through `species_face_fraction` with the donor side taken from the sign of the ELEMENT flux, the closing stage one minus their sum |
| `ionization_stage_face_flux` | `F_k(f) = x_k(f) N_el(f) - n_el(f) K(f) [x_k(j+1) - x_k(j)]/dr(f)` for the carried stages and the closing one, the closing stage's eddy term being minus the sum of the others', zero eddy term at f = 0 and f = N |
| `ionization_stage_divergence` | `[A_+ F_k(j) - A_- F_k(j-1)]/V_j` on `spherical_face_area_and_cell_volume`, the one geometry of L30 |
| `ionization_stage_face_jacobian` | `dF_k(f)/dx_k(j)` and `dF_k(f)/dx_k(j+1)` at the donor-cell face value |
| `ionization_stage_row_jacobian` | the tridiagonal `aa, bb, cc` of the divergence, built from those face derivatives and the SAME areas and volume the divergence divides by |
| `stage_simplex_projection` | the carried fractions returned to the simplex by scaling, with what it had to move measured into `stage_simplex_sum_over_one` and `stage_fraction_under_zero` |
| `hydrogen_and_helium_nucleus_face_flux` | `N_H(f)`, `N_He(f)`, `n_H(f)`, `n_He(f)`, `K(f)`, `dr(f)` from the element operator's one public flux, with no face coefficient rebuilt |

The three conditions of L12b section 3 are structural, not conventions the
caller can break: C1 because the element flux is an argument and is never
rebuilt; C2 because the closing stage is not reconstructed; C3 because
`n_el`, `K` and `dr` are single arrays and the closing eddy term is
`-sum_k E_k`.

The Jacobian omits the limiter and the second-order part of the
reconstruction, exactly as `element_advective_face_coefficients` omits them
from the element row, and the size of the omission is MEASURED and printed
by the suite: 9.49e-3 relative, largest over the column, between the
reconstructed and the donor-value stage divergence.

### 2.3 Boundaries

The element operator carries no diffusive flux through f = 0 or f = N (its
face coefficient loop runs `j = 1, N-1`, READ), so the stage eddy term is
zero at both and the boundary face carries `x_k(f) N_el(f)` alone.  Which
cell that fraction comes from is the sign of `N_el(f)` and nothing else,
which is the rule `species_face_fraction` already applies: at the base the
reservoir's partition enters when the flux is inward and the domain's own
leaves when the base breathes; at the top the zero-gradient ghost repeats
cell N, so nothing is created there in either direction.  No stage flux is
imposed at an outflow that the elemental flux does not carry.

---

## 3. The validity range of "ions and neutrals move together"

MEASURED on `LHS1140b/models/.L26/fid_resolve` (LHS 1140 b, He/H = 2.13,
`He_Kzz = 1e9`, `He_diffusion: True`, 500 cells to 29.03 R_p), with
`src/utils/collisional_validity.py`'s transcription of the code's own four
collision limits, its structure scale `L`, and the operator's own ambipolar
field term `G = dln(n_e T)/dr`.  `w = D_in G` is the drift of an ion against
its own neutral; `t_in` is the non-resonant momentum-transfer time the
operator's coefficients give and `t_cx` the resonant charge exchange they
exclude, at `sigma_cx = 2e-15 cm^2` (READ,
`docs/collisional_validity.md` section 1).  `t_coll` is their harmonic sum.

| r [R_p] | T [K] | v [cm/s] | x(H II) | t_coll/t_flow (H) | \|w/v\| (H) | \|w/v\| (H, with CX) | t_coll/t_flow (He) | \|w/v\| (He) | \|w/v\| (He, with CX) |
|---|---|---|---|---|---|---|---|---|---|
| 1.02 | 1135 | 1.17e+00 | 1.40e-04 | 5.5e-11 | 0.763 | 0.611 | 6.5e-11 | 0.501 | 0.378 |
| 1.20 | 5000 | 2.27e+02 | 1.18e-02 | 1.6e-07 | 0.030 | 0.023 | 2.7e-07 | 0.018 | 0.013 |
| 2.00 | 4520 | 3.19e+03 | 4.29e-02 | 1.2e-05 | 0.493 | 0.368 | 4.1e-05 | 0.294 | 0.209 |
| 4.01 | 2235 | 1.58e+04 | 1.28e-01 | 5.6e-04 | 0.500 | 0.387 | 2.8e-03 | 0.313 | 0.229 |
| 7.95 | 1133 | 4.73e+04 | 2.90e-01 | 1.3e-02 | 0.588 | 0.471 | 8.5e-02 | 0.386 | 0.291 |
| 12.93 | 736 | 8.01e+04 | 3.95e-01 | 7.0e-02 | 0.621 | 0.509 | 6.99e-01 | 0.421 | 0.324 |
| 20.10 | 534 | 1.09e+05 | 5.32e-01 | 2.5e-01 | 0.619 | 0.516 | 3.53e+00 | 0.429 | 0.336 |
| 29.03 | 458 | 6.55e-01 | 6.55e-01 | 4.97e-01 | 0.569 | 0.478 | 9.22e+00 | 0.399 | 0.315 |

Crossings, MEASURED:

- `t_coll/t_flow = 0.1` at 14.51 R_p (H) and 8.25 R_p (He); the helium
  ratio reaches 1 at 14.15 R_p and the hydrogen one never does inside the
  domain.
- `|w/v| = 0.1` at 1.220 R_p (H) and 1.243 R_p (He) without the resonant
  channel, 1.231 and 1.273 R_p with it, 1.238 and 1.297 R_p with it taken
  on the same Chapman-Enskog convention as the rest (a factor 1.6,
  `docs/collisional_validity.md` section 1, READ).
- Largest over `r >= 1.2 R_p`: 0.625 (H) and 0.430 (He) without the
  resonant channel, 0.518 and 0.336 with it, 0.470 and 0.298 on the
  Chapman-Enskog convention.

**The statement.** The one-velocity approximation is exact only inside about
1.2 R_p, which is also where the Damkohler number returns the local
ionization root whatever the transport does, so where it is exact it does
not matter.  Outward of it the neglected drift is a systematic of order
0.02 to 0.6 of the advective stage flux, and it is not smaller than the
stage transport the option exists to carry.  Above about 13 R_p the
continuum stage equation is itself unvalidated by the run's own Knudsen
measure (the tool's verdict on this state: HYDRODYNAMIC RESULT UNVALIDATED,
max Kn 0.88 in the heating and acceleration region, MEASURED).  This is the
same picture L12b section 6 measured on the 45 R_p state (0.72 and 0.53
there against 0.62 and 0.43 here), on a different state and with the
resonant channel now added rather than estimated.

An ambipolar stage drift is a separate physical term.  It is not folded into
this item, and the approximation is carried at the site with its measured
size, in the module header.

The script is `<scratchpad>/L36/stage_drift_validity.py`; it reads a run
directory and changes nothing.

---

## 4. Stage F, and what stage D can and cannot say

### 4.1 The selectable pair

`EXHALE_transit.py` read `output/Hydro_ioniz_adv.txt` and
`output/Ion_species_adv.txt` and nothing else (READ, its lines 79 and 80
before this item).  The pair is now resolved by
`exhale_transit_lib.transit_state_files(path, EXHALE_TRANSIT_STATE)`:

- `adv` (the DEFAULT, unchanged behavior) the advection-corrected profile;
- `solution` the state the wind solver converged and the certification
  judged, `output/Hydro_ioniz.txt` and `output/Ion_species.txt`.

It is a PAIR: the temperature of one state with the composition of the other
solves neither set of equations, so both files follow one selection.  Every
saved product now carries `state_selection`, `state_files` and
`state_difference` in its header, the last being the largest relative
distance between the two states over the rows the tool reads, in T, n(H I)
and n(He 2^3S), or the statement that the other state was not read and the
distance is UNKNOWN.

`post_process_adv` is NOT retired: it remains the independent
discretization (first-order upwind marching at the bulk velocity, no eddy
term, no element drift) the solution is measured against, and every quantity
it alone produces stays.

### 4.2 The He I 10830 equivalent width

MEASURED on the certified atomic fiducial (the `_adv` twin of the same state
lives in `LHS1140b/models/atomic_scalar_gj1132_kzz1e9/HeH2.13/output/`;
both give the same solution curve, so they are one state):

| state | EW [percent A] | max depth (rotation and instrument convolved) |
|---|---|---|
| `_adv` | 1.5544 | 5.20 percent |
| solution | 20.3950 | 66.12 percent |

Ratio 13.12 (MEASURED).  L10 reported 13.3 on 15.5116 against
1.1702 percent A for a different model (READ), so this is the same factor
reproduced on the fiducial and now reproducible from the tool itself.

**What that number is and is not.**  It is the distance between the
advection-corrected composition and the composition the solve carries.  On
this run the solve's composition is the LOCAL ionization equilibrium, because
the stage rows are not carried.  What stages C and D will change is what the
solved composition IS; the factor above is the size of the thing they move,
not yet the answer they give.

The measured distance between the two states over the 500 rows: T 7.37e-1,
n(H I) 6.31e-1, n(He 2^3S) 9.86e-1 (MEASURED, largest relative).

### 4.3 The He 2^3S promotion rule

The rule fixed in L12a section 3 (READ) is: if the post-process moves the
He 2^3S column by more than 5 percent **when the three ionization fractions
already agree to their own gate**, the level joins the carried set.  That
precondition needs stage C, so the rule cannot be evaluated and the level is
NOT carried, which is also the default the rule's own Damkohler measurement
supports.

The systematic it would be judged against is MEASURED on the fiducial:

| quantity | solution | `_adv` | `_adv`/solution |
|---|---|---|---|
| He 2^3S radial column over 1 to 29 R_p [cm^-2] | 4.9292e+11 | 8.0114e+10 | 0.1625 |
| n(He 2^3S) at 10.03 R_p | | | 3.03e-02 |
| n(He 2^3S) at 20.10 R_p | | | 1.71e-02 |
| n(He 2^3S) at 29.03 R_p | | | 1.38e-02 |
| n_e at 10.03 R_p | | | 2.02e-01 |
| n_e at 29.03 R_p | | | 1.44e-01 |

The column moves by a factor 6.15, and it follows the electron density the
metastable's recombination source is proportional to, not the level equation
(L12a section 3 measured 6.045 and 3.05e-2 / 0.195 at 10 R_p on the 45 R_p
state; the agreement says the two measurements are of the same state).  **The
line is a statement about the ionization closure.**

---

## 5. Tests and impact

### 5.1 `src/tests/ionization_stage_flux/`: 15 rows to 25, all GREEN

The 15 rows of L12b are unchanged and still pass at the values that memo
reported.  Ten rows are new, and they run the PRODUCTION module on the
production grid constructor (Mixed, 500 cells, 50 base cells of 2e-4 R_p
under a geometric stretch to 30 R_p), on two columns (a smooth one and one
with helium driven onto its own simplex faces) crossed with three winds
(outflowing, reversed, alternating face by face).

| row | measured | tolerance |
|---|---|---|
| `production_stage_flux_sums_to_the_element_flux` | 2.0245e-16 | 100 eps |
| `sum_identity_holds_at_the_boundary_faces` | 1.0621e-16 | 100 eps |
| `sum_identity_survives_the_simplex_projection` | 2.1653e-16 | 100 eps |
| `stage_divergence_sums_to_the_element_nucleus_divergence` | 4.4998e-16 | 100 eps |
| `stage_divergence_column_sum_is_the_boundary_flux_difference` | 1.2421e-17 | 1e-13 |
| `the_simplex_projection_returns_an_admissible_state` | 2.2204e-16 | 8 eps |
| `stage_row_jacobian_matches_a_central_difference` | 1.8490e-11 | 1e-7 |

Three rows are RED by construction, each a rule the derivation forbids:

| row | measured | must exceed |
|---|---|---|
| `a_cell_centred_volume_breaks_the_column_sum` | 8.0694e-08 | 1000 eps |
| `a_full_eddy_term_in_every_stage_breaks_the_identity` | 1.4761e-04 | 1000 eps |
| `bounding_the_closing_face_value_breaks_the_identity` | 3.6208e-04 | 1000 eps |

The second is the derivation's central point measured directly: replacing
the stage term by the whole mixing-ratio eddy flux, which is the `Phi_x` of
the L12a design, leaves `sum_k F_k` off the element flux by one copy of the
element's own eddy flux.  The first is the L30 geometry in the stage rows:
weighting the divergence by `r_j^2 dr_j` instead of the shell volume breaks
the column sum on the stretched grid.

The simplex projection is exercised rather than watched: the state handed to
it is deliberately inadmissible (every third cell pushed 1.3 past the
simplex, every seventh driven to -1e-3) and it returns an admissible one.

RED before: the seven GREEN rows could not be written before this item,
because the module they call did not exist; the three RED rows are what
shows they are not vacuous.

### 5.2 Key-off byte identity

Control build: `EXHALE_L36c.x`, md5 `16f5345d07ce70cdac09227839d5188b`, a
build of the ENTRY TEXT of every file of this item in
`<scratchpad>/L36/ctl_src/` (not the tree's `EXHALE.x`).  Measured build:
`EXHALE_L36.x`.

The key is off in all three, so each is a control on a path this item must
not move.  Bounded runs, `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`,
identical bounds on the two builds: `wasp_full` 400 steps, `mol_diffusion`
its pinned 12000, `hp_front` its pinned 100.

| case | `Hydro_ioniz.txt` | `Ion_species.txt` |
|---|---|---|
| `wasp_full` | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |
| `mol_diffusion` | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |
| `hp_front` | data rows BYTE-IDENTICAL | data rows BYTE-IDENTICAL |

The full files differ in two header lines only, and for a reason outside
this item: the control binary was built in a scratch directory that is not a
git work tree, so its build stamp reads `git=unknown tree=clean` where the
measured build reads `git=3c73905ca8a7 tree=dirty`.  Every data row is
identical byte for byte.

**The one run product that actually calls the routine this item changed.**
`element_nucleus_face_flux` is reached in a run only through
`write_element_flux_profile`, which fires on `EXHALE_DIFFUSION_CHECK=1` or
with a lower-atmosphere profile in use, so none of the three cases above
exercises it.  MEASURED separately: `mol_diffusion` bounded to 300 steps
with `EXHALE_DIFFUSION_CHECK=1` on both builds gives BYTE-IDENTICAL data
rows in `Hydro_ioniz.txt`, `Ion_species.txt` AND
`element_flux_profile.txt`.  The counter restore is therefore a no-op on
that path today, which is what was expected: the writer already restored the
counters around its own call, and the restore was moved into the routine so
that the next caller inherits it.  Nothing in a run reads either counter
(MEASURED by reading: the only references outside `species_face_flux.f90`
and the suites are the save and restore themselves).

### 5.3 The suites that link the changed objects

`element_operator` (40 rows, the one suite that calls
`element_nucleus_face_flux` directly) and `species_face_flux` (15 rows) were
run on BOTH builds: every row's verdict and every printed value is
identical, and every row passes.

---

## 6. What the next item has to do, and where

Stages C, D and E need two files this item does not own.

1. **`src/modules/lower_atmosphere/diffusive_photochemistry.f90`**: change
   `ic_Hp` to a fraction per hydrogen nucleus and add `ic_Hep`, `ic_Hepp`
   as fractions per helium nucleus; form their face flux, their divergence
   and their Jacobian rows by CALLING `ionization_stage_transport`, so the
   old face rule of `carrier_face_flux` is not used for them; take the
   element nucleus fluxes from `hydrogen_and_helium_nucleus_face_flux`;
   replace the carrier headroom of an ionization row by
   `stage_simplex_projection`; lift the `thereis_mol` and
   `carrier_transport` early returns for a set that is ionization only;
   extend `carrier_source` to the three charges with the trial electron sum.
2. **`src/modules/files_IO/input_read.f90`** lines 2185 to 2199: the two
   refusals are what make `Ionization transport: True` an error in an atomic
   gas.  They are CORRECT today, because the operator that would carry the
   row in an atomic gas does not exist, and they must be lifted in the same
   increment that builds it, not before.  `write_setup_report.f90` then
   echoes the carried set.
3. Then stage E: the certification entries, one per carried row, with the
   tolerance anchored by a measurement (L12a section 2.4).  **No tolerance
   was added here**, because no anchoring measurement exists without a
   solved state.
4. `src/modules/post_process/post_process_adv.f90` writes the measured
   difference between the two closures to `pp.log` and to the file header
   (L12b section 8).  The transit half of that statement is done (section
   4.1); the `pp.log` half is in that file, which this item does not own.

Until 1 and 2 land, the module has no caller in a run.  It is linked into
`EXHALE.x` (one line in `SRC`) and exercised only by the acceptance suite.

---

## 7. Noticed outside this item's scope

- **Three uses of "pipeline" in files this item owns were removed**
  (`exhale_transit_lib.py` in the two doubled-line docstrings and
  `EXHALE_transit.py` line 1040), each replaced by what the code does
  ("the same integration as the He and Lya lines", "the same spherical
  chord integration as resonance_depth").  The word survives in
  `src/tests/krylov_and_dogleg/krylov_and_dogleg_tests.f90` (L30b report
  section 7 item 5 already records it); that file is not owned here.
- **`docs/input_schema.md` K15f is still stale** on the two counts L12a
  section 7 and L12b section 9 record: it states a `Coupled carrier solve`
  requirement and a fatal `error stop` that was removed on 2026-09-12.
  Stage C rewrites K15f anyway, so it is left for that increment.
- **`element_nucleus_face_flux` returns the HELIUM flux and the mass per
  hydrogen nucleus, not the hydrogen flux.**  The hydrogen nucleus flux is
  formed here from the binary closure
  (`N_H = [F_rho - F_rho Y_He - J]/m_1`), which needs the caller to pass a
  face mass flux in the same units the diffusive half is in.  That is
  stated in both headers.  A cleaner shape would have the element operator
  return both elements' nucleus fluxes directly; it is a change to the
  public signature L30 has just landed and to the `element_operator` suite
  that reads it, so it is reported and not made.
