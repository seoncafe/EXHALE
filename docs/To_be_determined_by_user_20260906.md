# Decisions the user has to make (state at the close of Phase 1, 2026-09-06)

Five items. Each gives the background, the options, what each option moves,
and the advisor's recommendation. Nothing here has been implemented; the code
state each item refers to is that of `docs/Update_EXHALE.md` section 6
(binary `b2c15e53c888d09fa90ffdad80443c3a`, `make check` byte-identical on
sixteen cases, `run_fcheck.sh` clean).

## 1. Decision 4: what happens when the carrier Newton fails on a step

**Background.** The molecular carrier system (H2, OH, H2O, CO, and H+ with
`Ionization transport`) is advanced by a backward-Euler step solved with a
damped Newton (`diffusive_photochemistry.f90`, `solve_carriers`). Batch 2b
made a failed status stop the run. On 2026-09-06 `oxygen_chemistry` hit that
stop on the 89-bin photon grid; the diagnosis (report
`oxygen_carrier_newton_report.md`, summarized in `Update_EXHALE.md` section 6)
found the Newton correct and a root present, but the damped line search
needed 141 iterations to cross a stiff transient (an H2 destruction time of
1.4 s against a 2.85 s step) while the budget was 30. The budget is now 200
and the case runs; a step with no root would still stop the run.

**Options.**
- (a) Keep the present state: budget 200, stop on failure. Simple; a rootless
  step, which is physically possible, still stops the run.
- (b) Reject the step and retry on a shorter dt, the structure the hydro
  step already has (`retry_step` in `EXHALE_main.f90`). A rootless step is
  solvable at a small enough dt, so the stop disappears. Needs a retry cap
  and a dt floor; costs run time when it fires.
- (c) Accept the unconverged state and continue (the pre-2b behavior).
  Silent accumulation of unconverged rows; not recommended.

**Movement.** (b) is byte-identical on every regression case that never
fails a step (all seventeen today).

**Recommendation: (b).**

**Decided 2026-09-06: (b).** A failed energy update refuses the attempted step; the controller retakes it at half dt; exhaustion exits 2 with the diagnostics. Item ACCEPT-1/2 of `Update_EXHALE.md` section 7.

## 2. The relaxation amnesty of the ionization sweep

**Background.** The ionization and chemistry equilibrium of a cell is a
nonlinear system. When it does not converge, acceptance class 4 of
`ionization_equilibrium.f90` (the "relaxation amnesty", for early transient
steps) accepts the non-root and continues. In `oxygen_chemistry` on the
89-bin grid, cell 7 accepted a state with normalized residual 5.2e4 at step
12, which moved n_e from 1.1e3 to 3.5e10 cm^-3 and n(H+) from 0 to 2.6e8 in
one step and made the carrier system stiff (item 1). On the 20-bin grid the
same cell never took that path, so which grid stops is close to chance.

**Options.**
- (a) Keep the amnesty but cap the residual it may accept (a measured
  bound, of order 1); above it keep the previous step's composition.
- (b) Remove the amnesty: a non-root always keeps the previous state, or
  rejects the step together with item 1 (b).
- (c) Leave as is.

**Movement.** This is an acceptance policy, i.e. what the code counts as an
answer. The fixed-step relaxation snapshots (`mol_*` at 12000 steps) can move
where the amnesty fired in their early steps.

**Recommendation: (a) first**, with the bound measured on the matrix, and
"keep the previous state" rather than a stop. It belongs before the Phase 3
coupled update.

**Decided 2026-09-06: (a), cap = 1e2.** The amnesty keeps a non-root only up to a residual cap measured on the regression matrix; above it the cell keeps the composition it entered the sweep with (acceptance class 6, counted and reported like class 4). MEASURED on the sixteen default cases: every amnesty event of the matrix is in the step-0 sweep of `oxygen_chemistry`, 111 of them between 1e-6 and 8.9 and 4 between 8.6e4 and 1.5e5, with nothing in between; the four large ones are at cells 85 to 88, which are exactly the cells whose energy update then fails. The cap is the smallest power of ten at least 3x above the first population. Item ACCEPT-1/2 of `Update_EXHALE.md` section 7.

## 3. Acceptance of the D0 document (the Phase 2 gate)

**Background.** `docs/d0_governing_system_20260906.md` (1908 lines) writes
down the system the code solves before discretization: unknowns and state,
hydrodynamics, radiation, species balance, energy sources and sinks, the
operator split and time levels, the steady solver, a term table with a
consistency flag, and 15 open questions. Every later change is judged against
it, so it needs the user's reading and acceptance. Items C4, C5 and C16 were
closed on 2026-09-06 and are marked so in the file.

**Structural consistency items to read first.**
- C1 and C2: after the composition update the pressure is rebuilt at the same
  T from the new particle count and written back to `u(3)` with no physical
  source, and the chemistry stage overwrites `rho`, so the mass row is not
  conservative; the steady residual carries neither term.
- C3: the semi-implicit energy update runs two Newton iterations and
  re-evaluates cooling once, so the returned `cool` belongs to the previous
  iterate; `heat` is frozen; there is no residual test.
- C7: Lyman-alpha cooling charges the full 10.2 eV as escaping in cells where
  `lya_rt` computes an escape probability far below one.
- C8: the H/He/metal cooling sum is defined twice (`util_ion_eq.f90` and
  `T_equation.f90`), and the `_adv` temperature is solved with a different
  channel set.
- C17 and C19: the steady residual has no species rows, and `Ionization
  transport` refuses `Solver: Newton`, so that configuration has no steady
  solver.

**What is asked.** Read sections 6 (operator split) and 8 (term table), and
answer the open questions of section 9 where the code admits two readings
(the `_adv` channel set; the trace-metal diffusion counter-flux; the loaded
SED quadrature at the table's own resolution). Phase 3 code changes start
after acceptance, not before.

## 4. The WASP-121 b stellar spectrum

**Background.** Huang et al. 2023 (ApJ 951, 123, section 2.2) built their
spectrum from a solar SED below 1700 A rescaled to F_XUV = 1.6e6 erg cm^-2
s^-1 (Salz et al. 2019: the band is lambda < 912 A) and Lyman-alpha 1.0e5,
with LLmodels above 1700 A. LLmodels are not available here. Two candidate
files exist in `inputdata/sed/`; a third is proposed.

| Candidate | Construction | For | Against |
|---|---|---|---|
| A `wasp121b_solar_huang2023.txt` | WHI 2008 solar shape below 1700 A, 6459 K blackbody above | reproduces Huang's 700-1170 A band shape to 3 percent; continuous at 1700 A | 1700-2583 A, where He 2^3S is photoionized, is 1.8x brighter than the measured F6V photosphere |
| B `wasp121b_wasp17_huang2023.txt` | MUSCLES WASP-17 (F6V) with the same normalizations | measured photosphere in the near-ultraviolet | a factor 61 step at 1700 A; the FUV continuum below it is unphysical |
| C (not built) | A below 1700 A, B's photosphere above | closest to what Huang did (solar XUV joined to a stellar photosphere); removes the blackbody excess in the 1700-2583 A band | one more file to build; the WASP-17 photosphere is a proxy (Teff, gravity, blanketing not those of WASP-121), so its own bias is unquantified |

**Recommendation: C**, with A and B kept for sensitivity runs.

**Decided 2026-09-07: C.** `wasp121b_composite_huang2023.txt` is the spectrum of record for `WASP-121b/` and `benchmarks/wasp121/` (already the file both `input.inp` name since the 2d switch); A and B stay as sensitivity inputs.

## 5. Step 2d: switching the planet inputs to their loaded spectra

**Background.** Decision 16 puts every planet folder and benchmark on
`Spectrum type: Load` with its literature SED. The files are built and each
was read by the binary; `WASP-52b/input.inp` already loads its eps Eri table, every other `input.inp` of record still states a power law.

**What to decide.**
- (i) Switch the inputs now (a configuration change recorded in each folder
  README; the reruns stay a separate instruction).
- (ii) HD 189733 b: `hd189733b_epseri_salz2016.txt` (eps Eri normalized to
  Salz et al. 2016 Table 3, the modelling paper the comparison is against;
  a = 0.031 AU already set) or `hd189733b_bourrier2020.txt` (the star's
  measured XUV from the MOVES III tables, mean of Visits B-E: 1.4x brighter
  in the EUV, 30 percent fainter in the Lyman-Werner band, same photosphere).
- (iii) WASP-121 b: the choice of item 4.
- (iv) `VULCAN/atm/stellar_flux/sflux-HD189_B2020.txt` (a third-party rebin
  of the Bourrier spectrum, EUV 21 percent below the published average,
  epoch unstated): leave or remove.

**Recommendation.** Switch now; HD 189733 b on the Salz file as the run of
record and the Bourrier file as the sensitivity run; WASP-121 b on candidate
C; leave the VULCAN file untouched (it is superseded, and its removal is a
separate decision about that directory).

## 6. The wavelength band of the H2 Lyman-Werner pumping (LW-NORM, 2026-09-06)

**What was found.** The FUV band ledger (LEDGER-FUV) showed the Lyman-Werner
absorbers taking 45 percent more photons than the 912-1110 A beam loses.
LW-NORM located the cause: the self-shielding table
(`h2_self_shielding_table.f90`) was built from a line list that runs to
1200 A (the lines longward of 1110 A start from vibrationally excited levels,
which are populated at the 700-3200 K of the grid) but is normalized per
photon of the 912-1110 A band, while the beam loss of `lyman_werner.f90`
is the DB96 equation (39) equivalent width of 912-1110 A alone. MEASURED: the
column integral of the pump cross section reaches 1.515 at the top of the
column axis, against 1.519 for the photon content of 912-1200 A over
912-1110 A of a flat spectrum; the shipped `oxygen_chemistry` ledger ratio
1.454 is that integral at its base column. The 1110-1200 A lines are real
absorbers: at saturation they carry about a third of the dissociations the
table rates. The generator inputs are now in the tree
(`src/utils/h2_shielding_lbl/`) and regenerate the shipped table byte for
byte, so either repair is a script edit plus a regeneration.

**What to decide.**
- (A) Rebuild the table on 912-1110 A (`wl_max = 1110.0`). The beam and the
  table then describe one band and the ledger closes; the 1110-1200 A lines
  are dropped, so the dissociation rate of a hot thick layer falls (the P38
  memo section 4.4 measured `f_shield` 1.03x at 1e19 cm^-2 and 3.7x at
  4.4e21 cm^-2 between the two bands at 1300 K). No user-facing change: the
  `Stellar LW flux` key keeps its 912-1110 A meaning.
- (B) Widen the Lyman-Werner beam to 912-1200 A: the table normalized per
  photon of 912-1200 A (`WL_HI = 1200.0`), the band-B1 (1110-1201 A) flux
  joined to the beam in `water_photolysis.f90` and `util_ion_eq.f90`, and the
  H2O and OH continuum absorbers of B1 sharing that beam the way H I and H2O
  share band B2. This keeps every absorber the line list has. It redefines
  the `Stellar LW flux` key (the band the user's number covers) and the
  `Lyman_Werner.txt` and `FUV_bands.txt` columns, and `mol_lyman_werner`'s
  input value 343.0 (a 912-1110 A flux) would have to be restated.

**Recommendation.** (B), because it is the physical description: the
photons between 1110 and 1200 A exist in every SED the code loads, the hot H2
absorbs them, and (A) would remove an absorption the code already knows how
to compute to make a ledger close. (A) is the fallback if the key
redefinition is not wanted now; it is the memo's own earlier ruling and a
one-line script change.

**Decided 2026-09-06: (B).** The Lyman-Werner band becomes 912-1201 A (band B1 merged into it, the `Stellar FUV B1 flux` key retired), the table regenerated per photon of that band, and the beam loses what the table absorbs. Item LW-NORM-B of `Update_EXHALE.md` section 7.

## 7. Inner-shell photoionization of the metals (REF-METALS, 2026-09-07)

**What was found.** The metal photoionization fits (Verner et al. 1996 Table
1) are outer-shell fits with a stated `E_max` per ion (34 to 558 eV); above
it Verner's own routine `phfit2` hands the outer shell to the 1995 subshell
fit and ADDS the inner shells (K, L). The photon grid runs to 1240 eV. Item
METALS-OUTER adopts the handover for the outer shell (so no fit is evaluated
outside its validity; Fe I's row otherwise rises to 49x the truth at 1240 eV).
The inner shells remain absent: MEASURED after the handover (METALS-OUTER),
at 1240 eV the outer shell alone is 0.1 percent (C II) to 7.0 percent (Mg I)
of an ion's true total (Fe I 0.34, Ca II 0.32 percent); flux-weighted over
the column of `wasp_full` and `mol_metals` the metals carry 2.2 percent of
the 100-1240 eV photoabsorption (14 percent at the one radius where they
peak), and that whole band is 0.37 percent of the 13.6-1240 eV absorption,
so even the largest inner-shell correction moves a run's total absorption
by a few tenths of a percent, in a band the wind is thin in. Inner-shell
photoionization ejects an inner electron and is followed by an Auger cascade
that leaves the atom two or more stages up, with the photoelectron and Auger
electrons' energies as heat.

**What to decide.**
- (a) Leave the inner shells out, stated at the code site with the measured
  share (the metal contribution to the total XUV opacity above 100 eV is
  small against H and He in every case run so far; measure and record it).
- (b) Include the inner-shell absorption as opacity and heating (the photon
  removed from the beam, `h nu - E_th(inner)` plus the Auger electron
  energy to heat) with the ion advanced ONE stage (the ladder's limit),
  stating the under-count of the charge state.
- (c) Extend the metal ladder so that Auger products (two stages up) are
  representable, then include the full cascade (Verner's `phfit2` gives the
  cross sections; the Auger yields need a source, e.g. Kaastra and Mewe 1993).

**Published precedent (READ 2026-09-07 from the user-supplied PDFs).**
Cecchi-Pestellini et al. 2009 (A&A 496, 863) sum the photoionization
heating over "all ionization channels k" of "individual heavy elements and
their ions given in Verner et al. (1993)" (their Eq. 1), i.e. inner shells
included, with the primary photoelectron `E - I_k` degraded by the Dalgarno
et al. 1999 cascade; the chemistry there is H/He and the metals supply
opacity and heat only. Locci et al. 2022 (PSJ 3, 1) carry 128 species of
H/He/C/N/O whose neutrals have "singly charged, positive counterparts" only
(no doubly charged ions), and their secondary-electron source (their Eq. 3)
"includes also the contribution of Auger electrons" (Locci et al. 2018); the
K-shell edges of C (0.28 keV), N (0.40 keV), O (0.53 keV) appear as spikes
in their ionization profiles. So the published photochemistry models take
the inner-shell absorption as opacity, put the photoelectron and Auger
electrons into the secondary cascade, and advance the ion one stage: option
(b) exactly, with the charge under-count accepted by construction. No 1-D
escape hydrodynamics code (ATES, AIOLOS, Wind-AE, CETIMB) includes inner
shells at all (option (a)). Fox 2008 (SSRv 139, 3) reviews the solar-system
ionosphere practice, where Auger electrons are routine.

**Recommendation (revised).** (b), following the published exoplanet
photochemistry models: the Verner subshell cross sections above each
inner-shell threshold added to the metal opacity (`phfit2` gives them and
the code now carries the 1995 table), the photon's energy minus the
inner-shell potential handed to the secondary-electron degradation the code
already has (`electron_energy_degradation`, the Auger electron's energy with
it), the ion advanced one stage with the under-count stated at the code
site. The measured size (2.2 percent of the 100-1240 eV absorption, 0.37
percent of the total) says the effect is small either way; (b) is preferred
over (a) because it is what the field does and it removes a known omission
rather than recording one. (c) stays a ladder redesign for later.

**Decided 2026-09-07: (b), the revised recommendation.** Item METALS-INNER of `Update_EXHALE.md` section 7.

**Landed 2026-09-07.** `metal_photoion_sigma` returns the shell sum and
reproduces `phfit2`'s sum to better than 1e-10 at 17 ions x 12 energies from
20 to 1240 eV (MEASURED; 130 of those 204 rows disagreed before). The
subshell rows come from `phfit2.f`'s own `PH1` DATA statements, whose
thresholds differ from the printed 1995 table in six rows of five ions
(MEASURED). Every subshell turn-on is now a bin edge of the photon grid, so
the thresholds the grid is cut at go from 17 to 60 in `wasp_full`, 20 to 63
in `mol_metals` and 10 to 22 in `lower_profile`, while the bin COUNT is
unchanged (309, 289, 289). The handover at `E_max` is continuous where it
should be: the outer shell alone FELL by 84x for Fe I at 66 eV and 152x for
Fe II at 76 eV, because the 1996 row there stands for the whole 3d+4s
complex and the 1995 row is 4s only; with the subshells in the sum those
steps are 1.23 and 1.13, while the genuine K edges rise (C I 21x at 291 eV).
Impact, bounded snapshots on scratch copies with `OMP_NUM_THREADS=1`:
`mol_base_handoff` (metals off) byte-identical; `wasp_full` 300 steps
T +9.5 percent at most inside the wind and base heating +26 percent;
`mol_metals` 400 steps base heating +26 percent, T within 0.14 percent;
`lower_profile` 200 steps base heating +21 percent, T within 2e-5; log10
Mdot unchanged to the printed 0.01 dex in all four (7.21, 7.94, 8.38, 7.95).
The metal share of the 100-1240 eV absorption at fixed composition goes from
2.23 to 3.49 percent in `wasp_full` (2.20 to 3.63 in `mol_metals`, 2.30 to
2.67 in `lower_profile`), and that band is 0.38 percent of the whole. The
charge-state under-count, the one-electron cascade and the neglected
fluorescence are stated at the table in `cross_sec.f90` and in
`docs/photoion_cross_sections.tex`. The fluorescence neglect is now
quantified from the user-supplied Krause 1979 scan: with `e_top` = 1240 eV
only C, N, O and Na I can have a K hole and their omega_K are 0.28, 0.52,
0.83 and 2.3 percent (his Table 3, p. 315), while every reachable L-shell
event is radiationless to better than 0.7 percent (effective yields, his
Table 5, p. 320, and omega_3 of Table 3). Goldens NOT refreshed: `wasp_full`,
`mol_metals` and `lower_profile` will move when they are.

## 8. What the coupled source loop's stopping test measures (COST6, 2026-09-08)

**DECIDED 2026-09-08: (b).** Implemented as item COST7.

**The question.** The coupled source loop (`coupled_source:` in
`EXHALE_main.f90`) stops when one pass MOVES the pair (T, composition) by
at most `csm_T_tol` = 1e-8 (relative) and `csm_comp_tol` = 1e-8 (absolute
mass fraction). COST5 MEASURED that the sequence contracts geometrically at
a ratio |theta| of 0.07 to 0.22, so the ERROR of the accepted state, its
distance to the fixed point, is |theta/(1 - theta)| times the last
increment: one to two decades BELOW the 1e-8 the test names. The loop is
therefore one to two decades stricter than its tolerances say, and it pays
for that in passes: COST6 ESTIMATED 1.6 passes of 7 on `mol_base_handoff`
and 1.0 to 1.5 of 10 on `wasp_full` (the sweep is 84 percent of marching
time, so that is 15 to 23 percent of it), against the 6.3 percent COST6's
extrapolation, which keeps the present test, could take on the same case.

**The options.**
- (a) Keep the test as it is: the tolerances name the last increment and
  the accepted state is one to two decades inside them. Nothing moves.
- (b) Redefine the test on the ESTIMATED ERROR: stop when
  |theta/(1 - theta)| times the increment is at most the same 1e-8, with
  theta the global Rayleigh-quotient ratio COST6 already computes and
  guards (both increments non-zero, |theta| < 0.8, the model residual within
  tolerance; fall back to the present test where the guards fail). The
  accepted state then sits at the tolerance rather than below it; every
  output of every case moves at the 1e-8 level (a golden refresh, reported,
  user-run) and the `coupled_source_step` suite's pass-count bounds are
  re-measured.
- (c) Same as (b) but at 1e-9, which returns roughly the present accuracy
  for roughly the present cost, and only makes the statement honest.

**Recommendation: (b).** The tolerances were set (B3c) as statements about
the accepted state, and the state IS the object the certification and the
users read; a test that measures the increment instead is an implementation
accident that the geometric decay, now measured, lets us correct without a
model of the physics. (c) is the conservative alternative if the 1e-8
movement of every golden is unwelcome now. Under (a) the COST6 extrapolation
(default off, `EXHALE_CSM_EXTRAP=1`) is the only remaining lever and buys
6.3 percent on the molecular gate.

## 9. The advection post-process on a state that is not stationary (ADV-ENERGY, 2026-09-08)

**DECIDED 2026-09-08: (a).** Implemented as item ADV-REFUSE (decision 14 / T8.1 of the target-system document is thereby taken).

**What was found.** The `_adv` post-process re-solves the steady advection
equations along the recorded flow. Its energy equation had dropped the
enthalpy flux of the mass-flux divergence (ADV-STATIC), which is exact only
where `rho v r^2` is stationary; on a snapshot it gave the adiabat, 0.78 K
in a 1124 K column. ADV-ENERGY restored the term, so the equation is now
exact, and on the same snapshot it gives 1.7e7 K: with `F = rho v r^2`
falling 8 percent from cell to cell, the exact steady equation reads that
as a compression at the same rate (`w ~ rho^(gamma-1) F^(-gamma)`). MEASURED
movement of T_adv: the Newton-converged `wasp_full_newton` 1.4 percent; the
`du`-stopped `wasp_full` 16 percent at 1.088 R_p (He 2^3S column 16
percent) because its breathing base violates continuity in 161 of 502 cells
(restored term up to 67x the kept ones); the relaxation snapshots 300 to
450 percent, `lower_profile`'s He 2^3S by 238x. Neither the old nor the new
number is a temperature the gas has where continuity fails: a steady
post-process has no meaning on a non-stationary cell. The transit tools
(`EXHALE_transit.py`) read the `_adv` columns.

**The options** (decision 14 / T8.1 of `docs/b1_target_system_20260906.md`,
owner unassigned until now).
- (a) **Exact equation, and the correction is REFUSED where its assumption
  fails**: a cell where the restored term exceeds the kept terms (the state
  is not stationary there, the ratio the run already reports) keeps the
  run's own T and equilibrium composition, exactly as the thermal Damkohler
  condition already does where radiation dominates. The `_adv` product then
  means "the steady advective correction wherever a steady correction
  exists, the run's state elsewhere", and it is bounded. MEASURED cost
  (ADV-STATIC, route 2): rewrites 162 of 502 `_adv` temperatures of
  `wasp_full` by a median 5 percent (max 8.7), all in the breathing base;
  a Newton-converged state is essentially untouched (2 cells).
- (b) Exact equation with no refusal (the code as it stands): correct on a
  stationary state, unbounded on any other; `hydrostatic_column`'s `_adv`
  golden would pin 2.6e10 K.
- (c) Revert to the dropped-term form as the "stationary-flux projection",
  marked at the site: bounded, but it is not the equation of any state and
  it was 16 percent off on `wasp_full`'s base too.

**Recommendation: (a).** It is the only option under which every `_adv`
number is either the exact steady correction or the run's own state, and
the refusal criterion is the equation's own term ratio, not a tunable. The
transit tools should be told which cells were refused (a flag column, or
the count in the header), so a spectrum built on a breathing base is
labelled as such. Under (a) `wasp_full`'s He 10830 prediction will move at
the percent level; the Newton-converged states will not.

## 10. A coupled steady solve on a molecular configuration that eliminates H2 (B5h, 2026-09-08)

**DECIDED 2026-09-08: (a).** Implemented as item B5i.

**What was found.** With `Coupled carrier solve: True`, `Molecular
chemistry: True` and `Molecular carrier transport: False` the stationary
residual is not a function of its unknowns: the local-equilibrium
elimination does not determine the shielded layer's H2 content (a marginal
direction of the sweep, retention 0.77 of any seed perturbation, MEASURED
at 1e-6 and 1e-2 alike), and the base cell's energy row amplifies that by
~4e3, so two evaluations of ONE state differ by 3.2e5 times `Resid tol` per
unit relative seed change (0.96 `Resid tol` for the seeds a Newton actually
produces). Every convergence claim on `mol_diffusion`'s element solve (B5e,
B5f) rests on one sequence of evaluations. With `Molecular carrier
transport: True`, n(H2) is a Newton unknown and the same measurement reads
9.7e-10, below `Resid tol`. The content is upstream and transport data,
which is exactly why the code imposes it at the ghosts.

**The options** (`input_read.f90`, where the keys are resolved).
- (a) **Refuse the combination at startup** with a message naming the key
  to set (`Molecular carrier transport: True`), the way the SED-coverage
  and retired-base-key refusals already work. Nothing silently changes; a
  user who asked for a coupled steady solve on a molecular layer is told
  what it needs.
- (b) Turn the carrier transport on for that combination automatically and
  say so in the setup report.
- (c) Leave it: accept states at `Resid tol` on a residual whose seed
  dependence is 3.2e5 times that tolerance (the present behavior; what
  must not stay, in B5h's words).

**Recommendation: (a).** The rule elsewhere in this code is that a
configuration the method cannot solve is refused with the remedy named,
not repaired behind the user's back; (b) changes the unknown set of the
solve from a key the user did not set. Under (a) the `mol_diffusion`
element configuration used by B5b to B5h is refused as it stands and has to
carry the H2 row; the marching path is untouched (the refusal is on the
coupled steady route only).

**Also open from the same item, physics, not this decision**: the
ionization front of that state has two exact roots of the network at the
same incident field, 19 percent apart in x(H I) at r = 1.086, reached from
seeds far apart. Which root a stationary residual may land on (or a
statement that the front is unresolved on this grid) is a physics rule the
code does not have.

## 11. The carrier unknown of the coupled steady solve: density or its logarithm (B5j, 2026-09-08)

**DECIDED 2026-09-08: (a).** Implemented as item B5k.

**What was found.** The molecular carrier densities are non-negative
unknowns of the stationary system and the solver treated them as
unconstrained; on `mol_carrier` the outermost carrier is driven to and past
zero and the cell that binds the solve has its carrier at exactly zero
where its own row wants +2.75e-6, so no direction can be sampled there.
B5j MEASURED three arms on `mol_carrier`, none converging: the entry text
(`||R||` 1.929, 527 refused samples, step cut-backs to 2^-49); a
bound-aware trust region, now the default (1.415, 0 refused, cut-backs to
2^-4); and the carrier carried as `ln n` (`EXHALE_CARRIER_LOG_UNKNOWN=1` (the spelling at B5j; since B5k the log unknown is the default and `=0` restores the density):
**0.9972**, merit 17.6 -> 0.924, the smallest adopted unknown 2.2e-11
instead of an exact zero). The logarithm removes the bound by
construction, at the price of a floor under `ln n` and an asymmetric
region, and of changing the unknown space of the route that `Coupled
carrier solve: True` selects.

**The options.**
- (a) **Make the logarithmic carrier unknown the default** of the coupled
  route (one line), keep the bound-aware region for the element fractions
  (which are bounded on both sides and stay linear), and state the floor at
  the code site. Every measured number favours it; nothing converges yet
  either way, so no golden moves (no regression case runs the coupled
  route).
- (b) Keep the density unknown with the bound-aware region (the present
  default) and complete it with Bertsekas' two-metric projection (28 of 64
  steps are the Cauchy point alone, first-order descent).
- (c) Both: (a) now, (b) later if the log route stalls on a different face.

**Recommendation: (a).** A density is a positive quantity; a solver that
has to be told at every step not to make it negative is solving the wrong
variable. (b) is the published completion of the present default and is
not excluded by (a).

## 12. What a stationary state with a carrier at its bound means to the certification (B5j, 2026-09-08)

**DECIDED 2026-09-08: (c).** No certification change; revisit only if the logarithmic route stalls on a face.

**What was found.** With the bound-aware region the projected iterates pin
the carrier of cell 215 at zero, and that cell's carrier row reads 1.00000
of its own scale: to the certification an error (the row is not zero), to
the solver a constrained optimum (the row cannot be zero with the unknown
non-negative). A bound-constrained stationary state carries a multiplier on
the active face, not a zero row, and the certification has no notion of
one.

**The options.**
- (a) **The certification judges a carrier row at an active bound by its
  sign** (a KKT condition): a row that pushes the unknown further below the
  bound is stationary at the bound and is accepted, one that pushes it
  inward is an error and refused; the count of active-bound cells is
  printed and the state is marked as bound-constrained. Applies only to
  species rows (the hydrodynamic rows have no bounds).
- (b) Treat any active bound as a failure to certify (the present
  behavior): a stationary state may not have a carrier at zero. Under (a)
  of item 11 the bound is never reached, so (b) then costs nothing and
  says the physics: a wind carries a positive density of every carrier.
- (c) Reformulate so no bound exists (item 11 (a)) and leave the
  certification as it is.

**Recommendation: (c), and revisit only if the log route also stalls.**
Under item 11 (a) an active bound cannot occur; the certification then
needs no KKT branch, and adding one for a case that cannot arise would be
a rule nobody exercises.

## 13 to 19. Decisions of PLAN_20260909_rev1 section 9

Items 13 to 19 are stated in `docs/PLAN_20260909_rev1.md` section 9 (13 the
code pause, 14 the N4b route, 15 restart contract, 16 `_adv` status schema,
17 transit metadata, 18 two branches, 19 physical-time policy). Decided so
far: 13 (the user's "쭉 진행" of 2026-09-09, Stage A onward, with N16a in
parallel); 14 below; **15, 16, 17 DECIDED 2026-09-09, all option (a)**: 15 one
restart intent selector plus a versioned metadata block (design in
`docs/restart_contract_design_20260909.md`, items N10 and N9's metadata
half); 16 two `_adv` status fields with a versioned schema, legacy files read
as unknown validity (item N11); 17 comment metadata in `tpm_*.txt`,
numerical columns unchanged (item N13b).

**DECIDED 2026-09-09 (decision 14): route (i)**, the linearized constraint set
on the present unknowns (the advisor's recommendation in
`docs/constrained_step_design_20260909.md`). Implemented as item N4b.

## 20. The merit the step control descends and the gate the solve is judged by are different functionals (N7, 2026-09-09)

MEASURED on the atomic element reload (`backup/regression/atomic_elem_newton`,
N7 report): at outer iteration 60 the eight element rows hold 0.437 of the
squared merit `||F/Drow||^2` and nothing of `||R||` (a maximum of each
hydrodynamic row against its own certification scale), while the energy row
holds 1.000 of `||R||`; accepted steps swing `||R||` by 40 percent while the
merit moves by 1e-4 of itself; the element rows' own certification measures
sit at 1e-3 against `cert_tol_element = 1e-8`. The step control therefore
descends a functional that is not a proxy for the gate, and the proximate
refusals (the ray test near its floor, N7 P1/P2) are only what stops it
first. Options:

(a) **One functional.** The merit reads every row, hydrodynamic and element
and carrier, on its certification scale (`Drow` replaced by the
certification scaling for the element rows too), so that descending the
merit is descending the gate. Cleanest; changes the step lengths of every
species-row solve (the three-unknown route is untouched since its rows
already share one scaling); B5j to B5l measurements become historical.

(b) **Two controls.** Keep the merit for the step and add the element rows'
certification measures to the acceptance test of a trial (a trial that
worsens an element row's measure beyond a stated factor is refused even if
the merit falls). Smaller change; two functionals remain and can disagree.

(c) **Report only.** Keep both, print the gap (N7's item 7 does), and let N8a
(row scales, boundary derivatives) decide after its measurements.

Advisor recommendation: (a), preceded by P1 and P2 of the N7 report (the
ray-test floor measured, the slope test not overruling a sound ratio), which
are numerical repairs independent of the choice and can be done first.
Whichever option is chosen, `cert_tol_element` and `cert_tol_carrier` stay
UNCONFIRMED until N8b.

**DECIDED 2026-09-09: (a)**, one functional, with N7's P3 folded in. Implemented as item N20.

## 21. A restart that changes one physics option (N10, 2026-09-09)

The restart contract (decision 15 a) refuses a state whose `options` block
differs from the run's, exactly. That is right for a stationary state and it
forbids the arm-ladder workflow this project uses (converge without an
option, restart with it on, converge again: `armA_noLW -> armA_LW`, the He/H
rungs, the H2 channel rungs). Legacy pairs (no block) are unaffected; every
rung produced from now on would be refused. Options:

(a) **A key naming the tokens allowed to differ**, e.g. `Restart option
change: LW, sec_ion`, which permits only those tokens to differ and still
refuses everything else; the change is written into the new state's
provenance. Explicit, auditable. Advisor recommendation.

(b) Refuse only the tokens that change the unknown layout or the grid
(carriers, rows, N) and report the rest. Silent about which physics changed.

(c) Keep the exact refusal; ladders start from `Load IC? False` or strip the
block by hand. The workflow loses its restart.

**DECIDED 2026-09-09: (a)**, a key naming the tokens allowed to differ. Implemented as item N10b.

## 22. The species-row certification tolerances (N8b, 2026-09-10)

`cert_tol_element` and `cert_tol_carrier` are 1e-8, UNCONFIRMED. N8b MEASURED
the five anchors on the two candidate states (memo
`docs/certification_tolerance_anchoring_20260910.md`, section 7 has the full
text): the arithmetic, closure and representation floors are five to eleven
decades below 1e-8 (reproducibility exactly 0; closure seed 5e-15; state
perturbation 5e-17 in the wind, 4e-13 at the base); the derivative a Newton
step can resolve is 3.8e-8 at the binding cell; the DISCRETIZATION is about
5e-4 at the binding radius on the production grid (78 percent of the 6.0e-4
the atomic candidate stands at; N = 1000 halves it; the manufactured column
converges at order 1.59) and the ELEMENT FLUX conservation is 2.6e-2 in the
wind and 39 in the layer. No proposal certifies either candidate. Options:

(a) **Adopt the anchored values**: 1e-5 in the wind (the operator's own
stated floor and the discretization at N = 500), the layer reported and not
gating until its conservation is repaired (section 3 of the memo).
(b) **Keep 1e-8 and state that the fixtures cannot certify** on this grid.
(c) **A two-regime tolerance, both gating** (1e-5 wind, a layer value to be
anchored after the conservation repair).

**Re-measured after N29 (N8c, 2026-09-10)** on states whose mass closure is
round-off (1.9e-9): the floors are 0 / 9e-17 / 5e-19, the derivative resolves
3.9e-8 at the binding cell and 6.5e-9 in the wind, the discretization is
4.6e-4 at N = 500 (78 percent of the refused row; the manufactured column
gives 1.52e-6 at the production spacing, order 1.59), the layer's element
flux spread is 9.8 (was 39). The reading does not move.

**Advisor recommendation: (a)**, 1e-5 gating for cells at r >=
`cert_regime_wind_r`, the layer's element and carrier rows REPORTED and not
gating until its flux conservation (9.8) is repaired, the report saying in
one line that a wind-certified state is not layer-certified. No candidate
certifies under it (3.0e-4 and 7.2e-2 against 1e-5), so nothing that exists
is admitted by the change; what it does is make the gate say something true
about the wind. (b) remains honest if the user prefers a gate that claims
nothing.

**DECIDED 2026-09-10: (a)**. Implemented as item N30.

## 23. The rounding floor of the flux difference in the base layer (N33, 2026-09-10)

N31 to N33 traced why the species-row stationary solves are held by their
linear solve (40 of 40 Krylov products at a true relative residual 0.99
against 0.1 at the binding iterate, N27). The finite-difference Jacobian
action is not additive along the preconditioned Krylov directions because
the residual carries a non-smoothness floor (N31); the floor follows no
inner tolerance and is not in the chemistry (N32); and N33 MEASURED what
it is: the ROUNDING of the flux assembly, raised above the last bit of the
row by two cancellations. In the nearly hydrostatic base layer the Roe
flux is built from jumps of the reconstructed states 3e5 to 6e6 times
smaller than the states themselves, so the interface flux carries the last
bit of O(1) quantities; the row then divides the flux difference by the
cell volume, `r^2/dV` = 5.1e3 for a base cell of width 1.955e-4 (`Base
grid` key K33b, `dr_base = 2.0e-4`). The bound `epsilon x (face state) x
r^2/dV` is what a correctly rounded assembly cannot go below, and the
measured floor follows it over 4.6 decades along the column at a ratio of
5.5 to 12, on both fixtures (atomic reload and carrier reload), with no
branch, no fallback and no inner iteration flipping anywhere. The floor is
irreducible in this discretization: no probe rule (N31), no inner
tolerance (N32) and no arm inside the solver can lower it. Its consequence
is confined to the stationary solves with species rows (the three-unknown
route certifies; the marching path is unaffected, its steps never divide
by the cell width at this precision). Options:

(a) **Well-balanced flux differencing in the base layer** (RECOMMENDED):
assemble the hydrodynamic rows from the DEPARTURE of the face states from
the local hydrostatic isentrope, which `base_boundary.f90` already
integrates for the ghost construction, so the balance that now cancels in
floating point cancels analytically (hydrostatic reconstruction in the
sense of Kappeli and Mishra 2016, A&A 587, A94; LeVeque 1998, JCP 146,
346; the published references are to be read in full before the design is
fixed; `docs/b4_spatial_operator_design_20260906.md` sections on the
well-balanced pressure-gravity form are the in-house starting point). A change to `Num_Fluxes.f90`, `RK_rhs.f90` and `Source.f90`, so it
touches the MARCHING path too: built as a default-off arm, measured on both
reload fixtures by named outcomes (the floor must fall to the epsilon level,
the additivity defect's exponent must turn positive, the Krylov cycle must
reach 0.1 within 40 products at the binding iterate) and on the regression
matrix for its movement, then adopted as default only at a gate with the
goldens refreshed and the movement reported. This is the physically right
discretization of a layer in near-hydrostatic balance regardless of the
solver, which is why it is recommended over (b) and (c).

(b) **Extended precision for the hydrodynamic rows of the stationary
residual alone**: evaluate reconstruction, Riemann flux, geometric and
gravitational sources and the flux difference in `real(kind=16)` (quad)
from the double-precision state, inside the stationary solve only (the
face-state jumps already carry the rounding, so the whole pipeline from the
state to the row has to be in quad, not the difference alone). Lowers the
floor by about eighteen decades at a cost of roughly ten in wall time on
those lines, leaves the marching path and the goldens untouched, and treats
the symptom: the discretization stays as it is. As the CONTROL experiment
of (a) it is run first (N34) regardless of the choice.

(c) **Coarser base cells**: the floor is proportional to one over the cell
width, so a wider `dr_base` lowers it in proportion (and a finer one, the
opposite: the 2x base-grid key that fixed the HD 189733 b base wall makes
this floor twice as high). Trades resolution of the layer for a lower
floor; the base sound-wave physics of the layer (P44, P54) set the present
width, so this is a measurement, not a remedy.

(d) **Accept the floor** and judge the species-row solves by the
certification only, which the anchored tolerances already do (decision 22):
the Krylov cycle stays at 0.99, the arms stay chaotic at the last bit
(N26c), and no species-row configuration will certify below the floor's
reach in the layer. Honest, and it closes the program at "not converged".

Recommendation: (a), with (b) as the quick control experiment inside the
same item (if quad accumulation of the flux difference restores the Krylov
cycle's convergence, the diagnosis is confirmed before the discretization
is changed; if it does not, something else is still in the way and (a) is
not yet justified).

**N34 (the control experiment, 2026-09-10) changes the weighing.** With the
hydrodynamic rows of the stationary residual evaluated in quadruple
precision (default-off arm `EXHALE_RESID_QUAD=1`), the floor fell only by a
factor 2.1 and landed on `epsilon_double x face state x r^2/dV` exactly
(ratio 0.93 to 1.4): four fifths of the floor was the assembly's own
arithmetic, the rest is the double representation of what the assembly is
handed (the ghosts, the sweep's temperature and particle count) and of the
row itself, amplified by the same face cancellation and the same `r^2/dV`.
So N33's amplifier is confirmed a second way, and ONLY removing the
amplifier (a well-balanced form, option (a)) or widening the whole state
reaches the floor; option (b) cannot. But the Krylov stall is NOT held by
the floor alone: where the binding row is hydrodynamic (atomic reload, cap
40) the arm makes the Arnoldi image faithful (gap 1.19e-1 to 3.81e-3) and
the cycle still exhausts 40 products at 0.21 against 0.10; at cap 250 the
binding rows are the element Fe rows of the front cells (atomic) and the H2
carrier row of cell 205 (carrier), which do not pass through the flux
assembly, and their cycles stand at 0.98 to 0.99 with the arm on. What
holds them is the CONDITIONING of the preconditioned operator on the
species rows, not measured yet.

Revised recommendation: (a) remains the physically right discretization of
the near-hydrostatic layer and is still recommended, but it will not by
itself make the species-row solves converge, so it should not be sold as
that. The item that decides the species-row program is a measurement of the
preconditioned operator's conditioning on the element and carrier rows at
the binding iterate (the spectrum or a few singular values of the
preconditioned Jacobian restricted to those rows, and what in the banded
preconditioner misses their coupling: the diffusion coupling across cells,
the coupling of an element row to the hydrodynamic rows through
`mass_per_H`, the front). That measurement (N35) runs without a decision;
(a) waits for this one. Option (b) is retired as a remedy (kept as the
measured arm it is); (c) and (d) stand as written.

**N35/N36 (2026-09-10 evening) close the linear-algebra branch.** The band
misses no coupling; the preconditioned operator is near-singular on the
species rows themselves (Ritz ratio 6.6e5 carrier, 1.1e4 atomic), the
binding carrier row's diagonal is 3.3e3 below its coupling to the
hydrodynamic unknowns of its own stencil, a longer Krylov cycle returns a
worse true step, and neither a different row scaling of the linear model
nor a true-residual cycle passes acceptance. So (a) stands as the right
discretization of the layer, and the question that decides the species-row
program is now the discretized element and carrier equations at the front
(what makes a row nearly independent of its own unknown). Recommended
order: the user decides 23; the advisor then briefs the front-row
diagnosis (a measurement item, no default changed) whatever the choice.

**DECIDED 2026-09-10 (user: "(1)", read as the recommendation (a)): (a)**, a
well-balanced flux difference in the near-hydrostatic layer, built as a
default-off option and measured by named outcomes before any gate.
Implemented as item N37.

OUTCOME (N37, 2026-09-10): implemented as `Well balanced:` (K46), default
off, exact on the discrete equilibrium; the rounding floor of N33 does not
move under it (the cell pressure is a rounded double); the key stays off,
with the momentum-row scale of `store_row_terms` as the prerequisite of any
gate decision. The program is closed.
