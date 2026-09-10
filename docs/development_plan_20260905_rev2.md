# EXHALE development plan, revision 2 (2026-09-05)

Supersedes `docs/development_plan_20260905_rev1.md` after the second review,
`docs/development_plan_20260905_review2.md`. Inputs, in order of authority:
the source at HEAD `35d9dd5`, the published papers in `references/`, the two
reviews, the audit (`docs/physics_numerics_audit_20260905.md`), and
`docs/code_status_20260905.md`. Every source and paper claim of the second
review was re-checked before adoption (section 1), and its diagnostics were
re-run from `docs/audit_20260905/run_review2_checks.sh`. Rev 0 and rev 1 are
kept as the record of what each review addressed.

Nothing here is implemented. No atmosphere was run.

## 1. What the second review established, and what was verified today

| Item | Claim | Verified today | Disposition |
|---|---|---|---|
| Q1 | `h2_photo_channels.f90` reads Chung et al. (1993) Table II as a proton-count ratio `(N_S + 2 N_D)/N_M`; the published detector records one pulse for a double event, so the column is `(N_S + N_D)/N_M` | Module header lines 39-60 state `P = N_s + 2 N_d` and call the paper's normalization an overcount. `pdftotext` of `references/Chung_1993JCP_99_885.pdf`, p. 886: "detector responds with a single pulse (in double ionization both ions arrive almost simultaneously and cannot be resolved)". Re-run: at 80 eV the routine gives M, S, D = 0.7972, 0.1803, 0.0225 against the detector inversion 0.7796, 0.1763, 0.0441; the E1 channel-sum and stoichiometry checks pass in both. | **New confirmed source error.** Corrected in Phase 1 (A8) and tested in E1 (section 3.3). |
| Q2 | Rev 1's energy residual `Delta u_th + sum eps^0 Delta n - dt [H_rad - C]` charges the chemical energy twice when `H_rad` is the already partitioned photoelectron heat | Arithmetic re-run: H photoionization at 20 eV gives -7.197 eV instead of +6.402 eV; H2 double channel at 80 eV gives 16.65 eV instead of 48.33 eV | Adopted. D0 item 3 rewritten with the full external-exchange convention (section 4.2). |
| Q3 | 51.4 eV is Yan et al.'s vertical threshold while Chung et al. observe a double-ionization appearance near 47-48 eV; Chung et al. report excited neutral fragments H(2p) + H(2s); the neutral channel in `util_ion_eq.f90` gives the whole remainder to heat | Yan p. 1049: "has a vertical threshold of 51.4 eV". Chung p. 886: "above the double ionization threshold (~47-48 eV)"; p. 887: "which dissociates into H(2p) + H(2s)". `util_ion_eq.f90` neutral-channel comment: "the whole remainder of the photon goes into the kinetic energy of the two H atoms" | Adopted. A5 extended (section 4.3). |
| Q4 | The CO ceiling in `limit_to_element_budget` rescales OH, H2O, H2 and H+ as well; distance to the ceiling is not a constrained reaction-transport residual | `diffusive_photochemistry.f90` 2003-2133: CO clamped to `co_equilibrium_density`, then carbon clamp, then the water family scaled into the oxygen left, then H2, OH, H2O and H+ scaled together against `hydrogen_available_to_carriers` | Adopted. CO leaves Phase 1; its closure is a D0 derivation with a user decision (sections 4.2, 8). **Closed 2026-09-06:** the CO row carries published destruction rates (item B3b-CO) and the ceiling is deleted (item CEILING-DEL); `limit_to_element_budget` applies conservation only. |
| Q5 | The thermal enthalpy flux must exclude formation energy; `U(3)` is kinetic plus thermal; proton molecular diffusion is set to zero with the ambipolar term named as missing | `UW_conversions.f90` builds `U(3)` from `energy_density_from_pressure` plus kinetic energy. `diffusive_photochemistry.f90` 1196-1214: "the molecular diffusion of the proton is set to zero", "adding the ambipolar coefficient is a physics decision of its own" | Adopted. D0 items 3 and 5 rewritten (section 4.2). |
| Q6 | `speed_estimate_ROE.f90` enters its rarefaction/shock estimates only when `Q = p_max/p_min > 2`; an equal-pressure expansion stays in the PVRS guess and returns NaN | Lines 51-59: `Q_user = 2.0`. Re-run with `WL = (1,-1,1)`, `WR = (1,1,1)`, gamma 5/3: `u_star = 0`, `cL_star = cR_star = NaN`; exact `p_* = 0.2246`, `c_* = 0.9577` | **New reproduced source defect.** A2 extended (section 4.1). |
| Q7 | The non-LTE factor jumps from the table's 0.9985 to 1 at the last density column; the LTE fit follows the unscaled Table 4 emission | `h3p_cooling.f90` 136-139 returns `s = 1` for `ln >= sLogN(nNs)` while `sTab(11,5) = 0.9985`. Re-run at 5000 K: 0.99850 just below 1e14, 1.0 at 1e14; `E_LTE(5000) = 3.427e-18`, Table 4 unscaled 3.4119e-18, scaled 4.7688e-18 (both columns present in the PDF) | Adopted. B2 rewritten (section 4.4). |
| Q8 | Backward-Euler species and energy residuals are first order; a split fixed point is not the unsplit stationary equation; step rejection must restore every cached state | Design argument, consistent with the marching loop read for rev 0 (Lie split: hydro RK3, then carriers, then sweep, then energy) | Adopted. Sections 4.2 item 8, 4.3, 4.6. |
| Q9 | `assemble_residual` calls `Reconstruct` from U and, in the viscous/conduction branch, `U_to_W` again with `Tc = W(3,:)/n_part`; under the caloric EOS that pressure can differ from the seed-derived one | `steady_residual.f90` 109-136 (`Reconstruct(u_in, ...)`), 175-176 (`call U_to_W(u, W)`, `Tc = W(3,:)/n_part`) | Adopted. Section 3.1 restated. |
| Q10 | The transported-H+ constraint is applied only when `bg_ready` or a restart is loaded, so a zero-seed input can be replaced by an equilibrium state before the tested path runs | `ionization_equilibrium.f90` 1203-1211: `... .and. (bg_ready .or. do_load_IC)`; `bg_ready` set true at 1896 after the first sweep | Adopted. Test construction rules in section 3.3. |
| Q11 | `run_check.sh golden` overwrites `golden/<case>/` in place; `REGRESSION_REL_TOL` defaults to 1e-3; a positivity-repair counter already exists and is reported | `run_check.sh` 242-254 (`cp` into `golden/$c/`), 158 (`REL_TOL=${REGRESSION_REL_TOL:-1e-3}`); `RK_rhs.f90` 25, 319-320 (`n_faces_flux_positivity_limited`); `EXHALE_main.f90` 2166-2168 prints it when positive | Adopted. Section 3.2 and A1 corrected. |
| Q12 | The H2 partition function must agree with the populations assumed by chemistry, electron degradation and infrared transfer; the trace-species EOS assumption needs a runtime condition; the monochromatic cell photon balance belongs with the Phase 3 source interface | Design arguments | Adopted. Sections 4.3, 4.4, 6. |

The review also read Dalgarno et al. (1999), Glassgold et al. (2012), Koskinen
et al. (2022) Appendix B and Schulik and Booth (2023) for the points in Q3,
Q5 and Q10. Those readings were not repeated today; the plan cites them as
the review does and marks them for confirmation when the corresponding
derivation is written.

## 2. Principle of the order (unchanged from rev 1, sharpened)

Derive first; implement the bounded algebraic fixes that need no shared
derivation; then thermodynamics and the coupled source update; then the
spatial operator; then the steady solver on that operator; then validation
and extension.

```
Phase 0  evidence, provenance, gates, test construction
Phase 1  bounded corrections needing no shared derivation      (A1, A2, A3-scale, A8)
Phase 2  the common physical system, on paper                   (derivation D0)
Phase 3  thermodynamics and the coupled source update           (A4+B3, A5, A6, B1, B2, B4)
Phase 4  the spatial operator and its validity conditions       (D1, D2, D3, mixture and charged transport)
Phase 5  the steady solver on that operator                     (C1, C2, C3, C4)
Phase 6  coupled validation and extension                       (Koskinen, matching, CO kinetics, IR, observables)
```

What moved since rev 1: the CO constraint model leaves Phase 1 and becomes a
D0 derivation with a decision; the Chung detector inversion (A8) enters Phase 1
as a bounded correction; A2 gains branch admissibility and a vacuum limit;
the elementary photon-conservation test and the local energy-conservation
test move from Phase 6 to the increment that first claims them (Phase 3).

## 3. Phase 0: evidence, provenance, gates, test construction

### 3.1 What each code path holds fixed (precision after Q9)

- Marching step (`EXHALE_main.f90` 1168-1212): from the post-hydro U, T from
  p and the current composition (`comp_T_from_p`); carriers transported and
  the sweep solved at that T; particle count refreshed; p rebuilt at the
  unchanged T (`comp_p_from_T`) and U from it (`W_to_U`); then the energy
  step. The rebuild at fixed T is audit N1.
- Steady residual (`steady_newton.f90` 555-606): from the unknowns U, T from
  the seed composition; sweep at that T; particle count refreshed; ghosts
  rewritten with the new composition; `assemble_residual(u, n_tot+ne, heat,
  cool, R)` with the interior U unchanged. Inside it, `Reconstruct` derives
  primitive states from U with the current caloric composition, and the
  viscous/conduction branch calls `U_to_W` again and forms `Tc =
  W(3,:)/n_part`. Under constant gamma the pressure at fixed U is composition
  independent and the two pressures agree; under the caloric EOS the
  recomputed pressure can differ from the seed-derived one. The
  seed-dependent quantities are the optical depths and rate coefficients
  evaluated before the local nonlinear solve; the local chemistry evaluates
  trial populations, and the returned heating and cooling are of the
  returned composition at the sweep temperature (section 170).
- Final closure loop (2636-2655): post-sweep composition fed back as the seed
  until the residual norm moves by less than 1e-3.

Consequences stand as in rev 1: the steady-side work is a self-consistent
evaluation, not a mechanical substitution; preservation of the converged
atomic solutions under Phases 3-5 is a hypothesis.

### 3.2 Baseline provenance and vocabulary (Q11)

- `run_check.sh golden` overwrites the reference in place. The order is
  therefore: (1) copy `backup/regression/golden/` to
  `backup/regression/golden_pre170_<date>/` and verify the copy (`cmp`);
  (2) write its README from what the headers carry (`git=6d07d48afd41
  tree=dirty`, run timestamps, input checksums) and mark toolchain and
  stopping reason **unknown** where the header does not carry them; (3) only
  then refresh. Provenance is not reconstructed from memory.
- Promotion of the post-170 case outputs to references is not justified by
  the fact that all twelve differ; it is justified once each output's
  identity (input checksum, binary md5, `tree=clean` at the intended
  revision) and purpose are recorded. The user performs the refresh.
- Vocabulary: byte-identical; within the declared relative tolerance (1e-3
  by default, and that is a tolerance, not round-off); moved. "Round-off"
  is not used for a 0.1% agreement. A change of a physical constant is a
  model change even when its effect is below the tolerance. No small
  movement is an argument for accepting a physical defect.
- `make check` runs the whole matrix; a bounded run is
  `backup/regression/run_check.sh check <case...>`.

### 3.3 Test construction (Q1, Q6, Q10)

- **Chung detector inversion** (new): with `sigma_i` the ionizing absorption
  cross section, `r` the measured signal ratio and `q_D = N_D/(N_S + N_D)`,
  `sigma_M = sigma_i/(1+r)`, `sigma_S = (1-q_D) r sigma_i/(1+r)`, `sigma_D =
  q_D r sigma_i/(1+r)`; the test reconstructs `(S+D)/M = r` from the returned
  channels and checks the proton and electron multiplicities (S one, D two)
  from the stoichiometric table, not from the channel weights. Added to
  `src/tests/e1_h2/`, which today checks channel sum and stoichiometry only.
- **Equal-pressure ROE expansion** (new): `WL = (1,-1,1)`, `WR = (1,1,1)`,
  gamma 5/3; finite, positive star pressure and densities, and `c_*` within
  a stated tolerance of 0.9577 (an approximate estimator is compared at the
  accuracy of its approximation, not to machine precision).
- **Transported H+** (three cases, H/He only): `hp_zero_seed`,
  `hp_trace_seed` (1e-12 n_H), `hp_front`. Because the carrier constraint is
  applied only when `bg_ready` or a restart is loaded, an input abundance is
  not evidence of the state the solver receives. Each case is constructed
  through a verified restart file whose H+ column holds the intended seed,
  and the test asserts H+, n_e, the element budgets and the base state
  immediately before `carrier_source` and the Jacobian assembly. The
  tiny-abundance derivative test and the transport-front test are separate.
- **Hydrostatic column** (mechanical): fixed composition, cooling and
  chemistry off, conduction off, a thermodynamically stationary state with
  compatible boundaries; the quantity measured is the discrete
  pressure-gravity residual and the spurious velocity. Radiative-chemical
  equilibrium is a separate test. Zero incident radiation alone is not an
  equilibrium condition.
- The audit and review diagnostics (`docs/audit_20260905/*.f90`, `*.py`) move
  to `src/tests/physics_probe/`, are built through the Makefile's dependency
  generator, and each printed quantity acquires a tolerance and a verdict;
  nonzero exit on failure.
- `make test` aggregates `element_census_tests`, `diffusion_tests`, the E1
  checks, the A2 kinetics checks, `residual_determinism` and the probe.

## 4. Phases 1-5

### 4.1 Phase 1: bounded corrections (A1, A2, A3-scale, A8)

**A1, LLF signal speed** (`Num_Fluxes.f90` 287): `a1 = max(abs(vL)+aL,
abs(vR)+aR)`. Reached by `Numerical flux: LLF` and by
`positivity_limited_fluxes`. The existing counter
`n_faces_flux_positivity_limited` (already printed when positive) is the
record of whether the repair path was reached; what is added is a
zero-valued line in the summary and a count of direct LLF use, so absence of
the message is not read as absence of the path. The counter shows a path was
reached; it does not predict the size of the resulting change.

**A2, ROE star states** (`speed_estimate_ROE.f90`): the exponent `1/z` on the
two-rarefaction bracket and the upstream densities on the shock compression
ratios remain; in addition, the branch selection is made on admissibility,
not on `Q > 2` alone: a PVRS guess with `p_* <= 0` or a nonpositive star
density is rejected and the two-rarefaction estimate is used; an explicit
vacuum test (`(2/(gamma-1))(c_L + c_R) <= v_R - v_L`) returns the vacuum
speeds. Tests: the equal-pressure expansion of 3.3, symmetric and asymmetric
rarefactions and shocks, a stationary contact, normalization invariance, and
positivity of every returned star state (finiteness alone passes negative
pressure over negative density and is not the criterion). `input_read`
refuses `Numerical flux: ROE` with the caloric EOS; the header states the
constant-gamma validity.

**A3-scale, carrier reference scales and returned-state diagnostics**
(`diffusive_photochemistry.f90` 376-377, 1680-1686, 1830-1838, 1889-1893):
one element-reference table for every carrier used by the Jacobian step and
the residual floor; `pct_newton_resid` evaluated on the state after
`limit_to_element_budget`; `solve_carriers` returns a status; a failed status
in a carrier without an active constraint rejects the step. What the CO
ceiling means is **not** settled here (Q4): until D0 supplies the constraint
model, a run in which `pct_co_ceiling > 0` is labeled in the setup report
and in `EXHALE_resolved.out` as "CO ceiling active: constrained model not
derived", and its outputs are not used for quantitative C/O results. This is
disclosure, not closure. **Overtaken by events on 2026-09-06:** the ceiling
was deleted (item CEILING-DEL) after the CO row was given published
destruction rates (item B3b-CO), so no such label is written and none is
needed; `EXHALE_resolved.out` carries seven `co_domain_*` keys reporting where
the destruction-only row is inside its stated domain.

**A8, Chung detector inversion** (new, `h2_photo_channels.f90` 39-60,
`h2_channel_cross_sections`, `frac_H2_double_of_protons`, the `yan_rho`
model): rewrite the header and the inversion to the detector observable
`(N_S + N_D)/N_M`; redefine the double fraction as the event fraction `q_D
= N_D/(N_S + N_D)` (Chung's "about 20% at 80 eV" is read in that observable);
rederive `yan_rho` from Yan's double-to-single ratio `rho_D = N_D/(N_M +
N_S)` via `sigma_D/sigma_i = rho_D/(1 + rho_D)` combined with the detector
ratio, with a nonnegativity check on S; recheck the Table II interpolation.
The paper's normalization `sigma(H+) + sigma(H2+) = sigma_i` is consistent
with the detector response and the "overcount" note in the header is
removed. Energy per channel does not change here; it belongs to A5. Gate:
the detector test of 3.3; the four channel fractions still sum to one and
close nuclei and charge. Expected movement: every molecular case with the
default `chung80` model above 51.4 eV; measured, then reported.

Independent and delegable in parallel: A1, A2, A3-scale, A8, the test
construction of 3.3, the baseline archive of 3.2.

### 4.2 Phase 2: the common physical system on paper (D0)

One document, `docs/governing_system_<date>.md`, reviewed by the user before
Phase 3 starts. Its items, corrected after the second review:

1. **Species and metadata**: mass, elements, charge, formation energy
   `eps_s^0`, internal-state model and the populations it assumes, transport
   class. He(2^3S) is a level of He I.
2. **Reactions and photoevents as one ledger**: stoichiometric matrix,
   progress rates, `A nu = 0`, `z^T nu = 0` per reaction; for each event the
   recipients: translational energy of each product, internal excitation
   (level named), escaping or reabsorbed radiation, fast-electron energy and
   its degradation. Each recipient is assigned once. The present
   hand-maintained lists (`molecular_reaction_heat`, the collisional
   ionization potentials in `eval_cool`, the recombination-cooling fits, the
   carrier `carrier_source`) are mapped onto the ledger, and the mapping of
   each existing term is derived, not assumed: removing the ionization
   potential from `coio` does not turn a thermal recombination fit into the
   full escaping photon energy.
3. **Energy convention** (Q2, Q5). Evolved variable: `U(3)` = kinetic plus
   thermal internal energy, as today; formation energies are not stored. The
   accounting identity for a homogeneous cell is
   `Delta u_th + Delta E_chem = integral of Q_external,net dt`, with `E_chem
   = sum_s eps_s^0 n_s` and `Q_external,net` the energy entering the
   represented material from radiation and external agents minus the energy
   leaving it. For a photoevent the entering energy is the full photon
   energy, including the part that ends in chemical storage; for radiative
   recombination the leaving energy is the escaping photon energy, which
   includes the released chemical energy as well as the captured electron's.
   Thermal photoheating already partitioned (h nu minus I) is therefore
   **not** placed on the right-hand side of this identity; it is a derived
   quantity, equal to `Q_external - Delta E_chem` for that event. The
   alternative formulation (thermal sources with chemical terms only where
   the thermal definitions do not already contain them) is recorded and
   requires the reaction-by-reaction derivation of item 2 before use. If fast
   electrons or excited states carry their own reservoirs, their increments
   are on the left-hand side or their elimination is derived. Restart
   metadata state the convention; loading an old state is verified.
4. **Thermal enthalpy flux and the chemical flux** (Q5): with `J_s` the
   number fluxes in the barycentric frame, `F_chem = E_chem v + sum_s eps_s^0
   J_s` and `d_t E_chem + div F_chem = sum_s eps_s^0 S_s`; subtracting this
   from the formation-inclusive equation gives the thermal equation, whose
   diffusive enthalpy uses `h_s^th = h_s^full - eps_s^0` and whose chemical
   source is `- sum_s eps_s^0 S_s`. The finite-volume form uses the same
   chemical face flux as the species rows.
5. **Charged transport closure** (Q5): `sum_s m_s J_s = 0` does not close the
   charged species. D0 chooses between (a) a justified common-velocity
   approximation with negligible charged drift, with the condition stated
   and checked at runtime, and (b) an explicit closure: zero net current
   `sum_s z_s e J_s = 0` with the ambipolar field, electron transport, and
   the matching momentum and thermal exchange (Koskinen et al. 2022 Appendix
   B, eqs. B9-B10, to be confirmed in the PDF when written). Today's code
   sets proton molecular diffusion to zero and is neither.
6. **Species continuity in the face-flux form**: advective flux from the
   Riemann mass flux and a reconstructed mass fraction; diffusive flux `J_s`
   with the constraints of items 4 and 5; spherical volume `(r_R^3 -
   r_L^3)/3`; base condition as a species flux through the actual face.
7. **Boundary constraint count** (Q10): count for the full differential
   system: incoming hyperbolic characteristics of the multispecies Euler
   system already include composition modes (no separate condition per
   advected species on top of them); diffusion, viscosity and conduction
   add their own requirements; flow reversal and the outer boundary are
   covered.
8. **Time integration and the stationary residual** (Q8): one continuous
   source function `S(n, T, radiation)` and one spatial operator `L(U, n)`;
   the physical time integrator is derived from them with its ordering
   (transport, chemistry, radiation, energy) and its expected global order
   stated. Backward Euler on the sources is first order; the RK3 hydro stage
   does not raise a Lie-split scheme above first order; if a higher order is
   claimed it is by a stated scheme (IMEX or Strang) and measured. The
   stationary residual is the unsplit `L(U, n) + S(n, T) = 0`; a finite-step
   map `(Phi_dt(x) - x)/dt` may accelerate convergence but its fixed point
   is tested against the unsplit residual. Step rejection restores the
   entire attempted state: species, energy, cached columns and optical
   depths, and any mutable background the rates read. A common physical
   time step is distinguished from local pseudo-time.
9. **CO constraint model** (Q4): either the one-sided instantaneous
   destruction model `d_t n_CO + div F_CO = -lambda`, `lambda >= 0`, `n_CO <=
   n_eq`, `lambda (n_eq - n_CO) = 0`, with the implied sink evaluated and its
   sign checked at the active bound, the released C and O routed to named
   species, the dissociation energy taken from the ledger, and the active
   set and its derivatives following the same equations as the returned
   residual; or the exclusion of configurations that activate the ceiling
   until the kinetic network of Phase 6 exists. The choice is the user's
   (section 8); D0 writes both.
10. **Validity diagnostics**: Knudsen number, ion-neutral and electron-heavy
    equilibration times, chemical relaxation times from the constrained
    Jacobian, infrared critical densities, optical depths, and the H2
    internal-population assumption (item 1) against the density at which
    the assumed levels cease to be thermalized.

### 4.3 Phase 3, sources: A4+B3, A5, A6

**A4+B3, the coupled species-energy source step**, from D0 items 2, 3, 8.
Residual for a cell over a source step, in the convention of item 3:

- species: `n_s^{n+1} - n_s^* - dt S_s(n^{n+1}, T^{n+1}) = 0` for integrated
  species; algebraic equilibrium for eliminated species;
- energy: `u_th(T^{n+1}, n^{n+1}) - u_th^* + sum_s eps_s^0 (n_s^{n+1} - n_s^*)
  - Q_ext^{n -> n+1} = 0`, where `Q_ext` is the integrated external exchange
  of item 3 (absorbed photon energy in full, minus escaping radiation from
  recombination, line and continuum emission, minus energy carried by
  escaping fast particles).

The derivative of the energy residual includes the composition derivatives
through the chemical Jacobian and the temperature derivatives of `Q_ext`; the
present `|dC/dT|` damping may serve as a trial iteration matrix with a
residual test and a step fallback. Root selection: the root continuous with
`T^*` under step halving. A failed cell subdivides; the count is reported.
The two-iteration solve of `energy_semi_implicit.f90` is replaced by this
step; in the atomic, composition-eliminated limit it reduces to a scalar T
solve with rates re-evaluated at the returned T and the residual tested (what
rev 0 called A7).

Increments and what each may claim (Q8 corrected rev 1's "nothing before
(iii)"):

- (i) scalar T solve with residual test and returned-state cooling; claims
  the energy-solve gate, not conservation;
- (ii) equilibrium composition solved **together** with the energy
  constraint of item 3 (the sweep at `T^{n+1}`, iterated with the energy
  identity to tolerance); this increment may claim the closed-cell gate for
  its explicitly defined equilibrium model, because an equilibrium
  elimination solved with the energy constraint conserves material energy;
- (iii) finite-rate integration of the species D0 item 10 does not
  eliminate; claims the time-resolved trajectory gates.

Measured before replacing: the present projection energy of audit N1
integrated over a run for one atomic and one molecular case.

Gates (Q2, Q12):

- isolated H photoionization at threshold and at 20 eV: chemical storage
  13.598 eV, thermal 6.402 eV under stated local deposition;
- radiative recombination with escaping energy; collisional ionization and
  its reverse; a photon-driven cycle with unchanged final composition (zero
  net species change with nonzero throughput and radiation exchange);
- closed reacting cell: radiation switched off means photoionization off
  **and** radiative recombination and spontaneous emission either off or
  their escaping energy included in the measured total; total thermal plus
  chemical energy constant to the nonlinear tolerance; equilibrium reached;
- elementary photon balance for one cell and one frequency: absorbed equals
  entering minus leaving, and each absorbed event enters the material
  ledger once (moved here from Phase 6);
- source-step refinement: the declared order measured under a common
  physical time step;
- the four Newton-finished atomic planets: Mdot and profiles before and
  after, plotted; movement reported.

**A5, H2 photoevent ledger** (`util_ion_eq.f90` 937-968 and the neutral
channel), extended by Q3:

- Double channel: impose the exact event sum `h nu = Delta E_chem +
  E_electrons + E_fragments + E_internal + E_radiated`, with `Delta E_chem =
  2 I(H) + D0(H2) = 31.675 eV`. Adopt a stated partition (vertical picture:
  fragments carry the vertical-minus-asymptotic difference; electrons share
  `h nu - 51.4 eV`) as an approximation with its range; the degradation
  fraction at the mean electron energy is an approximation to an integral
  over the electron distribution and is labeled with that distribution.
- Threshold region: the hard zero below 51.4 eV is part of the adopted
  vertical approximation; the observed appearance near 47-48 eV is recorded
  and the threshold region investigated before 51.4 eV is described as a
  measured onset. Closing the energy sum does not validate the threshold.
- Neutral and dissociative channels: the experiment reports excited
  fragments (H(2p) + H(2s) in the neutral window; an excited H fragment
  permitted in dissociative ionization). The ledger gains an
  internal-excitation recipient per channel: a represented level (the
  existing H(n=2) population where the code carries it), a controlled
  eliminated level, radiation, or collisional quenching, assigned once;
  where the branching is unknown a bounded approximation and its evidence
  are recorded. "Whole remainder to heat" is replaced by that assignment.
- Deposition: the stopping length of energetic fragments **and** of primary
  and secondary electrons against the scale height, where local deposition
  is assumed; recorded where it fails.
- Degradation, later reaction heat and the thermal EOS are reconciled so no
  excitation energy is counted twice (Dalgarno et al. 1999 section 5.3;
  Glassgold et al. 2012; to be confirmed in the PDFs when written).

**A6, reaction-heat terms**: the He(2^3S)/HeH+ associative term enters the
ledger; the oxygen reactions and the O(1D) elimination enter together. If
the oxygen ledger is incomplete when the phase closes, the option is labeled
and excluded from oxygen-rich production runs; disclosure, not closure.

### 4.4 Phase 3, thermodynamics: B1, B2, B4

**B1, one H2 thermodynamic model** (Q12 added): the `molecular_infrared_data`
level list is the candidate single source, subject to a completeness check
over 100-10000 K; from it `Z_rv`, `u_rv = k T^2 d ln Z/dT`, `c_v`, and the
standard-state `H^0`, `S^0`, `G^0` with translation, standard pressure,
formation reference, statistical weights and the code's mass convention
explicit, before the oxygen module's Shomate H2 entry is replaced.
Ortho/para: conversion timescale in the target gas stated; equilibrium or
frozen ratio as the argument supports. **Population consistency**: the same
internal populations are assumed by the EOS, the chemistry rates that
presume LTE internal states, the electron-degradation model and the infrared
emission; where the gas is below the thermalization density of the assumed
levels, either a derived non-LTE internal energy or a clearly limited
approximation with a runtime flag (D0 item 10). Trace-species EOS: the
Phase 3 accepted configurations satisfy a runtime condition on the
internal-energy and heat-capacity contribution of H2O, CO, OH and the
molecular ions; outside it the run is labeled and not accepted as a
production result until Phase 6 supplies those heat capacities.

**B2, H3+ cooling**, complete specification (Q7):

- Emission definition: Miller et al. (2013) Table 4 has unscaled and scaled
  columns (3.4119e-18 and 4.7688e-18 W molecule^-1 sr^-1 at 5000 K); the
  current fit follows the unscaled column (3.427e-18). The choice is made
  explicitly, and the non-LTE factor of Table 6 is normalized compatibly
  (its high-density emission approaches the scaled value).
- Temperature: one smooth positive function fitted to the emission data
  (Table 4 or the underlying line list), not to Table 5's polynomial outputs
  at their joins; continuity by construction, derivative continuity where the
  solver needs it; fit accuracy reported against the data with the
  interpolation uncertainty between sparse samples.
- Density: low-density collisional asymptote `s proportional to n(H2)`
  (exponent 1 imposed, normalization from a stated rate model or matching
  approximation), smooth transition to the table; **high-density end**: the
  table's last column is used as the value at 1e14 (0.9985 at 5000 K) with a
  smooth approach to 1 beyond it, replacing the present jump to 1 at
  `sLogN(nNs)`; derivative behavior at both table boundaries tested.
- Scope: the linear gate applies to collisionally powered cooling under the
  chosen population assumptions; radiative pumping, formation pumping and
  other colliders are separate recipients in the ledger and are not forced
  to zero with n(H2).
- The Koskinen et al. (2009) table as a keyed alternative is Phase 6.

**B4, rate validity metadata**: range and collider in every rate routine;
out-of-range counters. No behavior change.

### 4.5 Phase 4: the spatial operator (D1, D2, D3, mixture and charged transport)

From D0 items 4-7 and 9.

**D1**: carrier advective flux from the face mass flux and a reconstructed
mass fraction; common spherical volume; diffusive fluxes with `sum_s m_s J_s
= 0`, the thermal enthalpy flux of item 4, the chemical face flux shared with
the species rows, and the charged closure of item 5; base species inflow
from the actual face flux. Tests: frozen passive scalar (speed, order,
integrated conservation); column element budgets with boundary fluxes; a
diffusion-only column that closes mass, elements, the specified current
constraint and thermal plus chemical energy with the chosen boundary fluxes.

**D3**: the constraint count of item 7 documented in `base_boundary.f90` and
adjusted where the imposed set exceeds it.

**D2**: the discrete hydrostatic residual is evaluated first on the
mechanical column of 3.3; a spatial pressure-gravity imbalance is corrected
before any time-integration remedy is tried (a local time step cannot repair
a spatial imbalance). Then the base cycle is measured on `mol_sec_ion`, and
local time-step control or the face condition is revisited only if the
residual test passes and the cycle persists. Gate: hydrostatic residual and
spurious velocity decrease with resolution; no odd-even pattern.

Validity diagnostics of item 10 written per cell; the run summary states
where they were exceeded and which conclusions remain supported there.

### 4.6 Phase 5: the steady solver on that operator (C1-C4)

**C1**: species rows as Newton unknowns, variable in number; the residual is
the unsplit stationary equation of D0 item 8 built from the Phase 4 spatial
operator and the Phase 3 continuous source function, not from the finite-step
map. The `input_read` refusals (1845-1881) lifted one at a time, each
remaining one naming the phase that lifts it. Gates: a Newton-converged
solution against a time-converged reference of the same equations on the
same grid; the `hp_*` cases; and an isolation check that the option-off
atomic path is byte-identical to the immediately preceding, physically
corrected build (not to pre-Phase-3 results).

**C2**: elimination by relaxation time from the constrained chemical
Jacobian against advection, diffusion, temperature-change and
irradiation-change times; the `P r/v` rule retired. He+ transported where
required. The digitized Koskinen He+ points are Phase 6 validation.

**C3**: composition convergence gate on species, T, optical depths and
diffuse rates; the reload test. Where the atomic branch cannot equilibrate
on the flow time, transport off is a stated approximation and the run
summary says so; the user may select the restricted experiment, and the
report states that the restriction is not an equivalent physical option.

**C4**: the Jacobian action is a finite difference with step, rounding and
inner-solve error; directional derivatives verified over a range of `eps`;
independent, frozen and eliminated variables named; inner
radiation/composition closure verified separately and deterministically;
what the banded preconditioner omits stated. The `wasp_full_newton` base
energy row (||R|| 5e-4 against 1e-5) is attacked here.

## 5. Phase 6: coupled validation and extension

- Koskinen 2022 Model A rebuilt term by term in the audit's order. The local
  paper suffices to begin the equation and process comparison now (Appendix
  B states the diffusion formulation and the photodissociation limitations);
  the missing spectrum file limits an exact numerical reproduction only. The
  request to the authors remains a user-authorized action.
- Lower/upper matching with a matching-pressure sweep and flux/energy
  continuity.
- `_adv` post-processing from the same non-equilibrium state, or restricted
  to a domain where the neglected molecules are shown irrelevant.
- Photon conservation across bands and overlaps (the LW/oxygen overlap,
  competing absorbers, refinement at thresholds, recombination photons
  counted once); the elementary one-cell balance is already Phase 3.
- CO kinetics with photodissociation and shielding (replacing the ceiling
  model of D0 item 9), heat capacities of H2O, CO, OH and the molecular ions,
  a net-exchange infrared form zero by construction at equilibrium, band
  optical depth.
- Multidimensional geometry, magnetic confinement and stellar-wind
  interaction stay outside this plan.

## 6. Ownership of every closure

| Closure | Phase | Gate |
|---|---|---|
| Reaction and photoevent invariants (elements, charge, one recipient per energy term) | 2 (D0), 3 | one check for each reaction and channel in `make test` |
| Chung detector inversion | 1 (A8) | `(S+D)/M = r`; multiplicities from stoichiometry |
| ROE admissibility and vacuum | 1 (A2) | equal-pressure expansion, positivity of star states |
| Thermodynamic identity, standard-state functions, K_eq across modules | 3 (B1) | identity to 1e-6; cycle test |
| Population consistency and trace-species EOS condition | 3 (B1), runtime | flag and label |
| Isolated photoionization, recombination, photon-driven cycle | 3 (A4+B3) | exact event sums |
| Closed reacting cell | 3 (increment ii onward) | energy constant to tolerance |
| Excited fragments assigned once | 3 (A5) | ledger check per channel |
| Elementary one-cell photon balance | 3 | absorbed = in - out |
| Source-step refinement at declared order | 3 | measured order |
| H3+ surface: low-density exponent, both joins, high-density end, emission normalization | 3 (B2) | four independent checks |
| Frozen passive scalar, species conservation | 4 (D1) | order and budget |
| Diffusion-only column: mass, elements, current, thermal plus chemical energy | 4 (D1) | budgets close |
| Charged transport closure | 2 (D0), 4 | chosen closure's constraint satisfied |
| Hydrostatic residual (mechanical) | 4 (D2) | decreases with resolution |
| Boundary constraint count | 2 (D0), 4 (D3) | documented and matched |
| CO active bound | 2 (D0), 4 or 6 | implied sink sign, element and energy closure, or exclusion |
| Unsplit stationary residual vs finite-step fixed point | 5 (C1) | both evaluated, difference reported |
| Final-state residual, deterministic inner closure | 5 (C3, C4) | reload reproduces the state |
| Momentum/thermal coupling validity | 4, runtime | exceedance and supported conclusions reported |
| Band photon conservation, IR cavity, spatial refinement, lower/upper interface, post-processing consistency, molecular benchmark | 6 | as listed in section 5 |

## 7. Parallelism and effort

Delegable now, in parallel: A1, A2, A3-scale, A8, the test construction of
3.3, the baseline archive of 3.2, the projection-energy diagnostic, B4. The
CO constraint model is not in this list. Everything from D0 onward is
sequential at the interface level. No calendar estimate is given; each
increment is sized after its derivation exists and ends at its gate. A
schedule is never an argument for keeping incomplete physics.

## 8. Decisions for the user

1. Golden archive and refresh in the order of 3.2; unknown provenance marked
   unknown.
2. Energy convention: `U(3)` stays kinetic plus thermal with the full
   external-exchange identity of D0 item 3 (this plan), or a
   formation-inclusive total energy. D0 is written for the former and records
   the latter.
3. Charged transport closure (D0 item 5): common-velocity approximation with
   a checked condition, or the explicit zero-current closure.
4. CO constraint (D0 item 9): derive and validate the instantaneous
   destruction model, or exclude ceiling-active configurations until Phase 6.
5. D0 review and acceptance before Phase 3.
6. Oxygen ledger incomplete at the end of Phase 3: label and exclude, or hold.
7. Atomic-line in-hydro ionization transport where C2's criterion requires
   it: the restricted experiment is selectable, and is reported as an
   approximation.
8. Request to Koskinen et al. for the spectrum file and the heating per
   ionization.
9. Instruction-gated items stay gated: paper and poster re-convergence, the
   LHS 1140 b write-up, the public switch of the repository.

## 9. What this plan does not claim

- No defect's effect on Mdot, on the base cycle or on any profile is
  asserted; each is measured at its gate. The channel-probability change of
  A8 and the H3+ high-density change of B2 in particular are unmeasured
  until run.
- Preservation of the converged atomic solutions under Phases 3-5 is a
  hypothesis.
- The Koskinen heating factor is attributed to neither code.
- Existence and stability of a molecular steady state with transported H+
  is what Phase 5 tests.
- The paper readings the second review reports for Dalgarno et al.,
  Glassgold et al., Koskinen et al. Appendix B and Schulik and Booth were not
  repeated today; the citations above are the review's, marked for
  confirmation in D0.
- Acceptance of a review finding in this plan does not repair the source;
  the two new source errors (Q1, Q6) remain in HEAD `35d9dd5` until Phase 1.
