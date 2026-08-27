# Phase E: the profile handoff and the elemental-flux closure

**Status: implemented, 2026-08-27 — milestones E1 through E5, section 8.
The six decisions of section 9 were confirmed as proposed and are in the code
as written, except that the helium-flux arm of decision 5 was answered the
other way: the closure iterates helium as well as hydrogen. Sections 2-6 now
describe code rather than a proposal, and every test of section 7 has been
run. Two of them are answered against their first form and the answer is the
result, not a defect: T-E5 holds in the steady-flux window and fails in the
overlap window (section 3.4), and the overlap window itself is measured and
rejected on physical grounds. The changelog record is
`docs/Update_EXHALE.md` sections 76-79; the application is
`docs/lhs1140b_lower_atmosphere_plan_new.md`, Phase E.**

This document is the implementation design that
`docs/composition_restart_and_base_handoff.tex` (Part 2, section "Out of
scope here") deliberately left unwritten: the schema of the profile handoff,
how EXHALE consumes it, the Photochem adapter that produces it, the climate
step that fixes the cold trap, and the iteration that makes the elemental
fluxes of the two models agree.

Inputs it is written against, all read for this design:

- `docs/lhs1140b_lower_atmosphere_plan_new.md`, Phase E (the two parts and
  their acceptance criteria) and baseline rows 4, 9;
- `docs/composition_restart_and_base_handoff.tex`, Part 2 — the decision of
  2026-08-26 in favor of profiles, the field list, and the reasons;
- `docs/vulcan_photochem_comparison.md`, section "2026-08-26 - Phase P1" —
  the Photochem environment, the gas-giant workflow, and the traps;
- `docs/Update_EXHALE.md` section 74 — the five `base.inp` key categories,
  the `<El>_H_base` elemental keys, `src/utils/element_budget.py`, and the
  `EXHALE_resolved.out` extension;
- `docs/oxygen_chemistry_new_plan.md` P1-P4 and section 4.1 (the physics
  ranking);
- `docs/binary_diffusion_design.md` sections 7.3, 7.4;
- `LHS1140b/kzz_decision.md` (the adopted `He_Kzz = 1.0e9` and the judgment
  that a constant K_zz is the wrong shape).

---

## 1. Goal and acceptance

The goal is the one stated in the plan: make the composition an output.

> **Acceptance (plan, Phase E).** Lower and upper models agree on elemental H
> and He fluxes at the match; the cold-trap water abundance is documented;
> the converged H:He is insensitive to the initial trial flux.

Why this is not already true. `HeH` is an input scalar
(`input_read.f90:1240-1243` sets it from `HeH_base`), it is the Dirichlet
reservoir the diffusion operator imposes at the base
(`binary_element_diffusion.f90:410-425`), and nothing in the handoff carries
a gradient, so no flux condition can be stated. A steady lower-atmosphere
chemistry solution cannot determine its own hydrogen supply independently of
its upper-boundary flux; without a flux-continuity condition the composition
is a prescribed input and "the wind returns H:He" is circular.

Two quantitative anchors this design uses as tolerances rather than
inventing new ones:

- **The reproducibility floor.** On HD 209458 b, JFNK solutions of the same
  configuration differing by 0.1% in their boundary data spread by
  **0.025 dex in log10 Mdot** (`docs/vulcan_photochem_comparison.md`, P1.6,
  "The same measurement at the deeper handoff"). No flux-closure residual
  tighter than that is meaningful on a wind-level quantity.
- **The flatness criterion.** Radial invariance of the elemental face flux
  `4 pi r^2 (rho X v + J)` measured against the code's own `du` spread is
  test T8, already implemented as the diagnostic
  `write_element_flux_profile` (`binary_element_diffusion.f90:1641-1700`,
  file `output/element_flux_profile.txt`, gated on `EXHALE_DIFFUSION_CHECK=1`
  or on a profile being in use; it was `./diffusion_faceflux.txt` before
  2026-08-27).

## 2. The profile file

### 2.1 Name and key

Proposed file name **`lower_atmosphere_profile.dat`** in the run directory,
and one new opt-in key in `input.inp`:

```
Lower atmosphere profile: lower_atmosphere_profile.dat
```

The key is named for what the file is (the lower atmosphere's solution over
an interval), not for its role in the code. Absent key: present behavior,
unchanged, including the scalar `base.inp` path. Parsed in the
keyword-extension block beside `Lower atmosphere` (`input_read.f90:442-451`)
and consumed at the same point as the scalar handoff, i.e. the
`read_base_inp` call site (`input_read.f90:856`) — after
`run_lower_atm_prestep` and before the derived constants, so overrides reach
everything computed from them. This placement is also what makes the
elemental reservoirs work: `thereis_metals`, `melem_ab` and
`thereis_lowIP_metal` are derived after that call
(`input_read.f90:889-910`), a reordering made in P2 precisely so a handoff
element could be heard (`docs/Update_EXHALE.md` section 74).

### 2.2 Format

Plain text, one header block of `#` lines followed by a fixed-column table,
deep to shallow (decreasing pressure). Same conventions as every other
EXHALE product: a `# columns:` schema line that the Python loaders read, SI-
free CGS units stated per column, full double precision.

Header (machine-readable `key value` lines, all required unless marked):

| key | meaning |
|---|---|
| `solution_id` | sha256 over the lower model's configuration, mechanism file, thermodynamic data, stellar flux file and the elemental abundance vector. The fingerprint that makes "same solution" checkable. |
| `source_code` | `photochem` / `vulcan` / `analytic` |
| `source_version` | e.g. `photochem 0.8.4`, `photochem_clima_data 0.3.1` |
| `mechanism` | mechanism file name and its own sha256 |
| `stellar_flux` | flux file name and sha256, plus the dilution applied |
| `p_match_bar` | the matching pressure: where the EXHALE base is placed |
| `p_top_bar` | the shallowest level the file carries; must be `< p_match_bar` so that an overlap exists |
| `p_deep_bar` | the deepest level carried (diagnostic; EXHALE never reads below the match) |
| `trial_flux_H` | the elemental H flux imposed at the lower model's upper boundary for this solution, `g/s` outward positive |
| `trial_flux_He` | the same for helium (usually 0 for the first arm) |
| `iteration` | closure iteration index k that produced this file |
| `reached_steady_state` | `T`/`F` from the chemistry solver |
| `notes` (optional) | free text |

Columns:

| # | name | unit | note |
|---|---|---|---|
| 1 | `p` | bar | monotonically decreasing |
| 2 | `r` | R_J | from the lower model's own hypsometric integration, same convention as `vulcan_to_base.py:139-154` |
| 3 | `T` | K | |
| 4 | `n_tot` | cm^-3 | total particle density, molecules counted once |
| 5 | `rho` | g cm^-3 | |
| 6 | `Kzz` | cm^2 s^-1 | the profile that replaces the constant |
| 7 | `q_H2` | - | H2 volume mixing ratio, the molecular/atomic H partition |
| 8 | `q_H` | - | atomic H volume mixing ratio |
| 9.. | `X_He`, `X_C`, `X_N`, `X_O`, `X_S`, ... | - | elemental El/H **nuclei** ratios, summed over every carrier, one column per element the file carries |
| .. | `F_He`, `F_C`, ... | g s^-1 | upward elemental fluxes at that level, outward positive |
| .. | `F_H` | g s^-1 | upward elemental hydrogen flux — the closure variable |
| .. | `q_H2O`, `q_CO`, ... | - | the molecular abundances the EXHALE network carries or that the cold-trap statement needs |

Element columns are summed over all carriers, as `vulcan_to_base.py` already
does for `C, N, O, S` and for the `HeH_base` correction
(`src/utils/vulcan_to_base.py:82-118, 163-164`):
photochemistry moves nuclei between molecules without creating or destroying
them, so the nuclei ratio is the conserved quantity and the only one
`melem_ab` can hold.

### 2.3 Overlap and interpolation

- **The matching level** is `p_match_bar`. EXHALE's base cell is placed
  there: `p_base_bar <- p_match_bar`, and `T0`, `R0` are the profile's `T`
  and `r` interpolated at that pressure.
- **The overlap interval** is `[p_match_bar, p_top_bar]` in pressure, i.e.
  from the EXHALE base outward to the radius where the lower model stops.
  Both models are defined there; it is where flux continuity is measured
  (section 3.4). A file with `p_top_bar >= p_match_bar` has no overlap and
  is refused.
- **Interpolation onto the EXHALE grid**: linear in `log p` for every
  intensive quantity (`T`, `Kzz`, mixing ratios, elemental ratios), which is
  the convention the existing adapters already use
  (`vulcan_to_base.py:136-137`). Extrapolation is never performed: EXHALE
  cells above `p_top_bar` take the shallowest profile value for `Kzz` (the
  eddy coefficient is a lower-atmosphere property and no statement about it
  exists above the file), and take nothing else. The choice of `log p` as
  the interpolation variable, rather than radius, is deliberate: the two
  models' radius scales are only equal if their hydrostatic integrations
  agree, and pressure is what both solve on.

### 2.4 Relation to the scalar `base.inp`

The scalar file stays, as the single-layer special case
(`composition_restart_and_base_handoff.tex`, section "Decision", item 2);
four regression cases pin it. What it may **not** be is a second source of
the same quantities.

**Proposed rule (single source).** When `Lower atmosphere profile:` is set
and the named file exists, every `base.inp` key of the EOS-boundary,
elemental-reservoir and boundary-constraint categories is **refused** with an
`error stop` naming the offending key. Provenance comments are allowed and
diagnostic keys are ignored as they are today. All of

| scalar key | replaced by |
|---|---|
| `T_base` | column `T` at `p_match_bar` |
| `r_base` | column `r` at `p_match_bar` |
| `p_base` | header `p_match_bar` |
| `q_H2_base` | column `q_H2` at `p_match_bar` |
| `HeH_base` | column `X_He` at `p_match_bar` |
| `<El>_H_base` | column `X_<El>` at `p_match_bar` |
| `Kzz_base` | column `Kzz`, as a profile (section 3.2) |

come from the profile. This removes the same-solution problem by
construction for these keys, rather than checking it.

**The enforced same-solution check** is then about the pair as a whole, and
it is a refusal, not a message. Today the two consistency statements are
printed and nothing else (`input_read.f90:1272-1290`). Proposed:

1. Every file produced by an adapter carries `solution_id`, including a
   `# solution_id <hash>` comment line in a `base.inp` written by the same
   adapter.
2. If both a profile and a `base.inp` are present and either lacks a
   `solution_id`, or the two differ, EXHALE stops with the two ids printed.
3. If the profile's `iteration` and `trial_flux_H` are absent, the run is
   permitted but the closure driver refuses to use it (a hand-written
   profile is a legitimate one-shot input; it just cannot be iterated).

This is the "enforce, not print" requirement of
`composition_restart_and_base_handoff.tex`, section "What the profile file
has to carry", closing paragraph.

## 3. What EXHALE does with it

### 3.1 Base state

`T0`, `R0`, `p_base_bar`, `q_h2_base`, `HeH` and `X_<El>` are set from the
profile at the matching level, through the same doors the scalar keys use —
`set_element_abundance` (`input_read.f90:1294-1313`) for the elements, direct
assignment for the rest — so that `comp_mass_per_H`, `comp_ntot_bc` and
`comp_rho_bc` (`src/modules/functions/composition.f90:126-173, 239-242`)
build the base EOS from exactly the values the lower model reported. Nothing
in those functions changes.

`EXHALE_resolved.out` gains the profile's provenance
(`solution_id`, `source_code`, `source_version`, `p_match_bar`,
`trial_flux_H`, `iteration`) alongside the abundances it already writes
(`write_setup_report.f90:417-456`), so the closure driver and
`element_budget.py` read one authority.

### 3.2 The eddy coefficient becomes a profile

Today `he_kzz` is one scalar (`parameters.f90:154`) used in three places
inside the diffusion operator:

- `binary_element_diffusion.f90:660` — the diffusive time scale that sets
  the relaxation step;
- `:1357` — `Agrad = rhof*(Df + he_kzz)`, the binary element operator's
  gradient-plus-eddy face coefficient;
- `:1572` — `DK = Df + he_kzz` in the trace-metal kernel.

`LHS1140b/kzz_decision.md` and `docs/eddy_diffusion_kzz.tex` already record
the judgment that a constant is the wrong shape: the molecular coefficient
rises by nearly five decades between the base and 1.2 R_p, so no single
constant can both mix the base and follow the molecular coefficient outward,
and the homopause — the one consequential thing K_zz does — is where
`K_zz = D`, a property of two profiles. The profile file is what resolves
this.

**Proposed change.** A module-level array `kzz_cell(1-Ng:N+Ng)` filled once
after the grid exists: from the interpolated profile column when one is
given, and from the scalar `he_kzz` in every cell otherwise. The three sites
above take `kzz_cell(j)` (and its face average `0.5*(kzz_cell(j) +
kzz_cell(j+1))` where they take a face value, matching how `Df` and `Gf` are
already face-averaged at `:1352-1354`).

Byte-identity of the constant case is exact, not approximate: for a uniform
array, `0.5d0*(a+a) = a` in IEEE double, and `Df + kzz_face` is then the same
sum of the same two numbers as `Df + he_kzz`. The three regression cases with
diffusion off return at the operator's early exit and never reach this code;
`mol_diffusion`, which has `He_Kzz: 1.0e9`, is the case that pins the
identity.

`he_kzz` keeps its name and its key `He_Kzz`: it is the constant a run states
when it has no profile, and section 9 lists the alternative of renaming it.

### 3.3 The Dirichlet reservoir

`X_base = m_He_amu*HeH/(m_1 + m_He_amu*HeH)`
(`binary_element_diffusion.f90:410-411, 638`) already reads the global `HeH`,
which the profile now sets. Nothing changes in the operator. What changes is
the meaning: `HeH` is no longer a user input but the converged output of the
closure, and the run record must say which iteration produced it.

The trace-metal reservoirs a profile states are Dirichlet in the same sense,
and a restart has to respect them too. `load_IC` restores the metal densities
the restart file carries, so a seed written at an earlier reservoir would hold
that older value everywhere above the base cell while the boundary condition
already uses the new one -- and the elemental budget then misses by the amount
the reservoir moved. Each element the handoff states is therefore renormalized
as it is loaded, by one factor common to every ionization stage,

```
r_El = (stated El/H) / (El/H of the loaded state at the base cell),
```

which changes the normalization and nothing else: the ionization split and the
shape of the loaded profile are preserved exactly. An element the handoff does
not state, and a restart with no handoff at all, are left untouched
(`docs/Update_EXHALE.md` section 80).

### 3.4 Measuring the elemental flux over the overlap

The measurement the closure needs already exists. `write_element_flux_profile`
(`binary_element_diffusion.f90:1641-1700`) writes, per face,

```
F_He(r_f) = 4 pi r_f^2 ( rho_f v_f X_upwind + J_f )   [g/s]
```

together with `4 pi r_f^2 rho_f v_f`, and the T8 criterion is that the two
have the same radial spread. Phase E needs three things from it that it does
not do today:

1. **Hydrogen as well as helium.** `F_H = 4 pi r^2 (rho (1-X) v - J)` on the
   same faces, from the same `X` and `J` the step used. It is the closure
   variable and it must not be reconstructed elsewhere.
2. **Written unconditionally when a profile is in use**, not only under
   `EXHALE_DIFFUSION_CHECK=1`, and into `output/` with the other products
   rather than into the run directory root.
3. **A single reported number per element**: the flux over the overlap
   window, reported as median and relative spread over the faces between the
   base and `p_top_bar`, so the closure compares one number and can say how
   well-defined it is. The spread is the honesty term — a mismatch inside
   the spread is not a mismatch.

Note the standing hazard the same module records: below ~1.02 R_p the base
carries a standing sound wave and the spread of `r^2 rho v` is 10^2 to 10^4
times its own median (`binary_diffusion_design.md` section 7.3). The
overlap window must therefore be taken outside the base cells, on the same
`[j_min:N]` escape window the solver uses to declare the wind steady
(`EXHALE_main.f90:687-710`), intersected with the profile's coverage. If that
intersection is empty — the lower model stops below `j_min` — the closure
cannot be measured and the run must say so rather than quote a base-cell
number.

**Measured 2026-08-27 (E1), and it is the normal case, not the exception.**
`j_min` is the first cell with `r >= r_esc` (`define_grid.f90:193-198`), and
`r_esc` is 2 R_p in every case that carries diffusion. A lower-atmosphere
profile reaching from the microbar match to 1e-8 bar covers about 0.01 R_p
above the base — two decades of pressure is a few scale heights, not a
radius doubling — so `[j_min:N]` and the profile's coverage never intersect,
and E1 reported `lower_profile_flux_state window_empty` on every
configuration tried.

### The two windows, as E4 implemented and measured them

E4 replaced the single `[j_min:N]` window with two, because measurement said
one interval cannot answer the question:

**(a) The overlap window**, the interval both models describe. Its lower edge
is now *measured*, not assumed: the first face at which the 5-face moving
spread of `r^2 rho v` falls below 10 per cent, which is where the standing
base sound wave stops dominating; when no face qualifies the edge falls back
to 1.02 R_p and `lower_profile_flux_r_lo_source` says which of the two it
was. Its upper edge is `r(p_top)` from the profile's own hydrostatic column.
Both edges, the face count and the reduction go into `EXHALE_resolved.out`.

**On LHS 1140 b this window is empty, and it is empty structurally.** The
measured edges are `r_lo = 1.009 999 8 R_p` (from the spread rule) and
`r_hi = 1.009 840 R_p`: the profile stops 1.6e-4 R_p *below* the radius at
which the wind first becomes flux-flat. Extending the profile is not the
fix, and the reason is a measurement, not a judgment. Between the match and
the profile top the column covers 1.96 decades of pressure in 0.0098 R_p,
i.e. about 0.005 R_p per decade; the wind is still 6x above its far-field
flux at 1.05 R_p and does not settle to within 10 per cent until ~1.5 R_p.
Reaching 1.5 R_p on that scale would ask the photochemical column for of
order 100 further decades of pressure. Mapping `p_top` through EXHALE's own
`r(p)` instead of the profile's — the two hydrostatic scales differ, the
wind being hot where the cold column is not — puts the upper edge at
1.0734 R_p rather than 1.0098 and makes the window non-empty, but does not
rescue it: the flux over `[1.010, 1.073] R_p` has a radial spread of 5.6
(F_H) to 7.3 (F_He), against a closure tolerance of 0.05. Both mappings say
the same thing, that the interval the two models share lies wholly inside
the region where the wind's mass flux is not yet its own.

**(b) The steady-flux window**, `r >= r_esc`. At a steady state the elemental
flux `4 pi r^2 (rho X v + J)` does not depend on radius, so a flux measured
there *is* the flux through the matching level; that identity, and not a
convenience, is what makes this window a legitimate statement about the
handoff. On the converged LHS 1140 b profile run it holds 191 faces and
gives `F_H = 1.8143e7 g/s` and `F_He = 1.3204e7 g/s` with radial spreads of
0.44 and 0.46 per cent — flat to the same tolerance as the mass flux itself
(0.45 per cent), which is test T-E5.

That last number depends on how far the wind was converged, and the
dependence was measured rather than assumed: at `Resid tol: 1.0e-3`, the
value `finish_case.sh` uses, the same solution's steady-window spread is
7.4 per cent (F_H 7.1, F_He 7.9) and would fail section 6.2's own validity
test; continuing the same solution to `Resid tol: 1.0e-4` brings it to
0.50 per cent and moves `log10 Mdot` from 7.47 to 7.50. The 7 per cent was
under-convergence, not a floor. The closure driver therefore runs the wind
at `1.0e-4`.

## 4. The Photochem adapter

**`src/utils/photochem_to_lower_profile.py`**, with
`src/utils/vulcan_to_lower_profile.py` the cross-check arm behind the same
command-line interface. The VULCAN arm is a cross-check, not a second
production path: P1 left the code choice to the user, and P1's evidence table
gives Photochem the climate model, which is the only path away from a
prescribed T(p) (`vulcan_photochem_comparison.md` P1.7).

The driving pattern is already written and validated in
`docs/p1_matched_comparison.py`; the adapter should reuse it rather than
re-derive it:

- environment `/home/kiseon/.conda/envs/photochem_cmp/bin/python`, photochem
  0.8.4, `photochem_clima_data` 0.3.1 (P1.0);
- `photochem.extensions.gasgiants.EvoAtmosphereGasGiant(mech, flux_file, mp,
  rp, solar_zenith_angle=..., thermo_file=..., data_dir=...)`
  (`p1_matched_comparison.py:160-162`);
- mechanism from `zahnle_rx_and_thermo_files(atoms_names=['H','He','N','O','C'
  (,'S')], remove_reaction_particles=True)` (`:111-118`);
- stellar flux converted by dilution only, `(R_star/a)^2`, mW/m^2/nm at the
  planet (`:85-96`);
- steady state through `initialize_robust_stepper` / `robust_step` in blocks,
  recording `reached_steady_state` and `gave_up` (`:174-192`).

### 4.1 Configuration generated

Planet mass and radius (from `input.inp` or the closure driver), the stellar
SED, the mechanism atom list, `T(p)` and `Kzz(p)` (section 5), the elemental
abundances mapped **by name** into `gdat.gas.molfracs_atoms_sun`
(`:168-171`), and the upper-boundary elemental H flux of the current closure
iteration.

### 4.2 Failure policy

| condition | action |
|---|---|
| `reached_steady_state = False` or `gave_up = True` | no profile is written; the driver stops the iteration and reports the block log. A non-steady chemistry solution has no elemental flux to hand over. |
| `set_temperature` raises the thermodynamic-range error (P1.8 trap 2: NH2 above ~1900 K) | the T(p) input is truncated at the highest level the data admit and the truncation is recorded in the header `notes`; if the truncation falls below `p_match_bar` the run is refused, because the match would then sit outside the chemistry solution. |
| `3*TOA_pressure_avg >= P_climate_top` (P1.8 trap 1) | the supplied climate grid is cut at `3.05*TOA` exactly as `read_climate` does (`p1_matched_comparison.py:99-109`); the adapter must set the TOA *and* cut the grid, never one alone. |
| `p_top_bar >= p_match_bar` after all cuts | refused: no overlap. |
| element columns that do not sum to the input abundances to within 1e-10 | refused; this is the adapter's own conservation check and it is cheap. |

### 4.3 Traps carried forward

- `gdat.gas.atoms_names` is ordered differently from run to run (P1.8
  trap 4): map by name, never by position. This applies to the profile's
  element columns as well — they are named in the schema header and read by
  name.
- An EXHALE run directory must contain `output/` before the run starts
  (P1.8 trap 3), which the closure driver creates.
- The tool shell's `cp` is interactive; the driver uses `\cp -f`, as
  `finish_case.sh` already does.

## 5. Climate and the cold trap

**Status: implemented and measured on LHS 1140 b, 2026-08-27
(`src/utils/radiative_convective_column.py`,
`photochem_to_lower_profile.py --climate`, `docs/Update_EXHALE.md`
section 78). Items 1, 2, 4 and 5 below are as proposed; item 3 is corrected
by what the codes actually do, and the correction is recorded with it.**

Photochem's chemistry API takes `T(p)` and `Kzz(p)` as **inputs**
(`initialize_to_climate_equilibrium_PT(P, T, Kzz, ...)`,
`p1_matched_comparison.py:174`), so the climate solution is a separate,
explicit step and not a byproduct. It is
`src/utils/radiative_convective_column.py`: Photochem's `clima`
(`AdiabatClimate`) with `solve_for_T_trop` on, so the stratospheric
temperature is the skin temperature of the solution rather than a number
chosen for it, and with the background gas picked by abundance rather than
by name -- on the helium-rich atmosphere the LHS 1140 b line implies, the
background gas is He and not H2.

The procedure, as implemented:

1. Solve the radiative-convective column on the planet's parameters and the
   Phase-A stellar spectrum, from a stated deep boundary up to the
   photochemical model top.
2. The tropopause is where that adiabat meets the isothermal stratosphere,
   i.e. the temperature minimum; `clima` states it as `P_trop`/`T_trop` and
   both go into the profile header `notes`.
3. **Corrected.** The proposal said to set the chemistry model's rainout to
   that tropopause and enable the mechanism's H2O particle. Measured, both
   halves are wrong for a gas giant, and in opposite directions:
   - the H2O particle needs no enabling. The Zahnle set restricted to
     H/He/N/O/C already carries `H2Oaer`, and asking for
     `water-condensation: true` in the gas-giant settings is **refused**:
     *"Either "fix-water-in-troposphere" or "water-condensation" is turned
     on in the settings file, but the reaction mechanism already implements
     H2O condensation via a particle."* The cold trap is therefore applied
     by the chemistry solver already, and the settings switch must stay
     `false`;
   - rainout is a surface process. `gas-rainout: true` is refused for want
     of a `rainfall-rate`, and a gas-giant model has no surface for rain to
     fall to. It is not the route by which a hydrogen or helium atmosphere
     traps its water, and it is not used.
4. `f_H2O` at the tropopause is reported as the cold-trap water abundance,
   carried in the `q_H2O` column and quoted in `notes`, together with the
   elemental sums that bracket it: the profile always computes the El/H
   ratios twice, once over the gas phase and once including the condensed
   carriers, and their difference is the cold trap.
5. `Kzz(p)`: the climate solution provides none, so the run states one and
   the header says which. On LHS 1140 b that is the adopted constant
   `1.0e9 cm^2/s` (`LHS1140b/kzz_decision.md` section 0), written into the
   `Kzz` column.

One approximation is stated in the code with its validity range: the climate
solve needs a composition before any chemistry has run, so the elemental
vector is partitioned into the carriers a cool hydrogen/helium atmosphere
holds in equilibrium -- O in H2O, C in CH4, N in N2, He atomic, the rest of
the hydrogen H2. That is the equilibrium partition only below the CO/CH4 and
N2/NH3 transitions, and the solve is **refused** above a stated
`t_partition_max` (1000 K) rather than returning a column whose opacity
carriers are the wrong molecules. Nothing downstream inherits it: the
photochemistry computes its own composition from the same elemental vector.

### Measured on LHS 1140 b (2026-08-27)

`He/H = 2.09` (the diffusion-limit crossing of `kzz_decision.md` section 6),
Lodders C/N/O, GJ 1132 SED at the b orbit, `K_zz = 1e9`, `p_match = 1e-6`
bar. The deep boundary of the climate solve is its one free parameter and it
was scanned:

| deep boundary [bar] | `T_deep` [K] | `P_trop` [bar] | `T_trop` [K] | `f_H2O` above the tropopause |
|---|---|---|---|---|
| 1 | 242.4 | 0.488 | 186.0 | 3.33e-7 |
| 3 | 308.5 | 0.765 | 185.8 | 2.06e-7 |
| 10 | 439.9 | 0.974 | 185.4 | 1.49e-7 |
| 20 | 555.9 | 1.031 | 185.0 | 1.32e-7 |
| 30, 50 | refused: the pseudoadiabat leaves the thermodynamic data ("Failed to compute heat capacity") | | | |

The stratospheric temperature is the robust part: 185-186 K over a factor of
20 in the deep boundary, against the skin temperature `T_eq/2^(1/4)` = 190 K
for `T_eq = 226 K`. The tropopause *pressure* is not: it moves by a factor
of two over the same range, and the water abundance moves with it, because
at a fixed `T_trop` the trapped mixing ratio is `p_sat(T_trop)/P_trop`.

20 bar is the configuration adopted for the profile, and not for a physical
reason: it is the shallowest boundary at which the photochemical model
bottom fits inside the climate column. On an atmosphere this cold the
chemistry never equilibrates within any column `clima` can reach, so
Photochem's `determine_quench_levels` returns the deepest level of the grid
as the quench level and its default `BOA_pressure_factor = 5` then asks for a
model bottom five times deeper than the column it was given. The adapter
takes that factor as an option and states 1 here, with the reason in the
refusal text.

**Against the published estimate.** Cherubim et al. (2026) quote a cold trap
at a 0.1 bar tropopause with `T_skin = 194 K` and `f_H2O ~ 7 ppm` for
`T_eq = 226 +/- 4 K`, from closed-form estimation
(`docs/lhs1140b_lower_atmosphere_plan_new.md` section 2). The climate
solution gives **0.13-0.33 ppm** at a tropopause of **0.5-1.0 bar**, i.e. a
factor 20-50 less water at a tropopause 5-10x deeper. The two statements are
consistent in their saturation physics -- `p_sat(185 K)/1 bar` is 1.3e-7 and
`p_sat(194 K)/0.1 bar` is about 5e-6, so each number follows from its own
tropopause -- and they differ in where the tropopause is put. Ours is a
solved level in a helium-dominated column (`mu = 3.6` at `He/H = 2.09`); 0.1
bar is an assumed one. This is reported as a result about the climate model,
not tuned toward the published number.

**What the trap does to the handoff.** The elemental oxygen the gas phase
carries falls from `O/H = 6.062e-4` at the deep boundary to **4.957e-7** at
the matching level, a factor of 1223, while `C/H`, `N/H` and `He/H` pass the
tropopause unchanged (CH4 and N2 do not condense at 185 K). At the match
itself the gas-phase and all-carrier sums are equal, because the condensate
is left behind at the tropopause and not carried to 1 microbar: the cold trap
shows in this handoff as a depletion with height, not as a condensed
reservoir sitting at the base.

## 6. The flux-continuity iteration

**Status: implemented 2026-08-27 as `src/utils/element_flux_closure.py`,
with the trial fluxes actually imposed on the chemistry (E2 only stated
them). Both hydrogen and helium are iterated, which is the one place this
section's original text was overruled: it proposed holding helium at zero.
What is written below is the code, with the measured departures marked.**

### 6.1 The loop

```
k = 0:  Phi_H^(0) = trial elemental H flux (section 6.3)
repeat:
  1. Photochem run with Phi_H^(k) as the upper-boundary elemental H flux
     -> lower_atmosphere_profile.dat, iteration k
  2. EXHALE wind from that profile, He_diffusion on, co-converged:
     marching -> JFNK (EXHALE_PTC) -> post-processing pass
  3. measure Phi_H^EX(k), Phi_He^EX(k) over the overlap window (section 3.4),
     with their radial spreads
  4. residual  eps_k = |Phi_H^EX(k) - Phi_H^(k)| / |Phi_H^EX(k)|
  5. if eps_k <= tol: done
  6. Phi_H^(k+1) = Phi_H^(k) + omega ( Phi_H^EX(k) - Phi_H^(k) )
until k = k_max
```

`omega = 0.5` to start, halved (floor 0.125) on any iteration whose residual
failed to fall — the same damping rule, and the same reason, as the diffusion
outer loop (`binary_diffusion_design.md` section 7.3): nothing in a Picard
iteration of two solves keeps them from chasing each other. As with that
loop, **the convergence test reads the undamped residual**, so a small
`omega` cannot buy a false convergence. A secant update may replace the
under-relaxation once three iterates exist; it is an accelerator, not a
different criterion, and it falls back to under-relaxation whenever it
proposes a non-positive flux.

**Both elements are iterated**, each with its own trial flux, its own
undamped residual and the shared `omega`; convergence is on
`max(eps_H, eps_He)`. Holding helium at zero, as this section first
proposed, would have made the helium column of the closure a diagnostic
rather than a condition, and on a planet whose whole question is the H:He
partition that is the wrong half to leave open.

**How the trial flux reaches the chemistry.** `set_upper_bc(species,
bc_type='flux', flux=...)` on the Photochem model, with the elemental flux
converted at the model's own top radius,

```
phi_El [nuclei/cm^2/s] = Phi_El [g/s] / m_El / (4 pi r_top^2)
```

and imposed on the dominant carrier of each element: helium on `He`, and
hydrogen on `H2` at `phi_H/2`, two nuclei to the molecule. The dominant-
carrier step is measured and not assumed — the atomic-H share of the
hydrogen nuclei at the model top is 1.0e-4 in the converged LHS 1140 b
solution, and the adapter refuses the solution if that share rises past
1 per cent, because the flux would then have been put on the wrong carrier.
The sign convention was read off Photochem's own right-hand side
(`photochem_evoatmosphere_rhs.f90`, `rhs(k) = ... - var%upper_flux/dz`):
positive is a loss at the model top, i.e. outward, the same convention this
document uses. A zero trial flux calls `set_upper_bc` not at all, so the
one-way profiles of E2 and E3 remain reproducible — verified by reproducing
`solution_id 51302c76...` bitwise.

That the boundary condition takes is checked every iteration rather than
assumed: the adapter measures the model's own top-of-atmosphere elemental
fluxes with `gas_fluxes()` and writes them into the profile header as
`measured_flux_H` / `measured_flux_He`. On the LHS 1140 b reference arm the
imposed 1.800000e7 g/s comes back as 1.800162e7 g/s, a departure of 9e-5.
The closed-top solution's own numerical floor is 3.0 g/s in hydrogen and
18.2 g/s in helium, seven decades below the fluxes at issue.

### 6.2 Tolerance and iteration count

`tol = 0.05` on the elemental H flux — chosen against the measured
reproducibility floor rather than picked: the same configuration reproduces
`log10 Mdot` to about 0.025 dex, i.e. ~6%, between JFNK solutions differing
negligibly in boundary data (P1.6). A flux residual below that floor is not
distinguishable from the solver's own noise, and the closure must not claim
it.

**The validity test, as implemented.** The flux is defined only to within
its own radial spread over the window it was measured on, so a window whose
spread exceeds `tol` cannot decide a residual at `tol` whatever that residual
comes out to be. The driver therefore tests the *window* every iteration and
independently of how far the iterate is: `spread > tol` stops the run as
UNRESOLVED with both spreads printed. This is a statement about the
measurement, not about the iterate, and it is the reason the wind is
converged to `Resid tol: 1.0e-4` rather than `1.0e-3` (section 3.4).

`k_max = 8`, with every iteration's `(Phi_H^(k), Phi_He^(k), Phi_H^EX(k),
Phi_He^EX(k), eps_H, eps_He, omega, window used, window spreads, HeH_match,
log10 Mdot, solution_id, EXHALE info)` written to
`closure_history.txt` in the driver's run record, and the narrative,
including every refusal, to `closure.log`. The `solution_id` column is the
fingerprint chain: each iteration's profile is identified by the hash of the
configuration that produced it, so the sequence is auditable after the fact.

**Failure policy**, each recorded and then stopped on, never retried
silently: a nonzero adapter exit or a non-steady chemistry solution (the
adapter log tail is printed); an EXHALE run that does not report
`done info=0` (the run log tail is printed); no usable window; a window
spread above `tol`; and `k` reaching `k_max` with the last residual printed.

### 6.3 Insensitivity to the initial trial flux

The plan's third acceptance criterion is a designed test, not an observation:
run the loop from three starts spanning an order of magnitude, and require
the converged `HeH` at the match to agree across them to within the same 5%
and the converged `log10 Mdot` to within 0.025 dex. If the three starts
converge to different fixed points, the result is that the closure is
multivalued on this planet — a reportable finding, not a failure to be tuned
away.

**Run on LHS 1140 b, 2026-08-27**, from `Phi_ref = 1.8e7 g/s` in hydrogen and
`1.3e7 g/s` in helium (the fluxes the neighboring solution carries) and from
`0.3 x` and `3 x` those. The three arms are
`LHS1140b/exhale/flux_closure/{ref,lo,hi}`.

| start | k at convergence | converged `F_H` [g/s] | `F_He` [g/s] | `He/H` at the match | `log10 Mdot` | He 10830 red depth [%] | EW [mA] |
|---|---|---|---|---|---|---|---|
| `1.0 x` | 0 | 1.81433e7 | 1.32040e7 | 2.0923516 | 7.500 | 1.699 | 4.743 |
| `0.3 x` | 5 | 1.82031e7 | 1.34266e7 | 2.0923486 | 7.500 | 1.721 | 4.804 |
| `3.0 x` | 6 | 1.82119e7 | 1.34710e7 | 2.0923602 | 7.500 | 1.725 | 4.818 |

The three converged `He/H` agree to **5.5e-6**, `log10 Mdot` to better than
the 0.005 the log is printed at, `F_H` to 3.8e-3 and `F_He` to 2.0e-2, all far
inside the 5 per cent criterion. **T-E7 passes**, and the closure is
single-valued on this planet rather than multivalued.

Two things the run says about the loop itself. `omega` stayed at 0.5 in every
arm: the residual fell on every iteration and the halving rule never fired,
so the Picard map is a plain contraction here and the damping was not tested
by this planet. And the iteration count is exactly what a contraction with
`omega = 0.5` predicts — the residual halves each step, so the arm starting
1.98 away needs seven iterations and the one starting 0.70 away needs six,
both inside `k_max = 8`, with no margin to spare for the 3x arm. A start
further than an order of magnitude out would need a larger `k_max` or a
secant acceleration.

### 6.4 The driver

**`src/utils/element_flux_closure.py`**, run from the case directory:

```
python3 src/utils/element_flux_closure.py <case_dir> --config <json>
        --phi0-H <g/s> [--phi0-He <g/s>] --seed <converged output/ dir>
        [--omega 0.5] [--tol 0.05] [--kmax 8] [--resume] [--dry-run]
```

It owns the loop and nothing else: the chemistry is
`photochem_to_lower_profile.py`, the wind is EXHALE driven the way
`LHS1140b/exhale/finish_case.sh` already drives it (seed -> JFNK
under `EXHALE_PTC=1 EXHALE_PTC_JFNK=1` -> post-processing pass), and the
budget check is `src/utils/element_budget.py`. Two constraints from
`binary_diffusion_design.md` sections 7.3-7.4 that the loop depends on and
that are recorded there as done: the diffusion loop must exist on the direct
steady route (the `EXHALE_PTC` path every LHS 1140 b case uses), and
`load_IC` must keep each cell's element split with `He_diffusion` on, or the
restart at step 2 destroys the diffused state.

Everything that is a property of the planet rather than of the iteration —
both command lines' fixed arguments, the binary, the `input.inp` template,
the thread count — lives in one JSON configuration
(`--print-config-template` emits one), so the driver carries no planet
constants. The LHS 1140 b configuration is
`LHS1140b/exhale/flux_closure/lhs1140b_closure.json`.

Two departures from `finish_case.sh`, both measured rather than
precautionary:

- `Resid tol: 1.0e-4` instead of `1.0e-3`, for the reason in section 3.4:
  at `1.0e-3` the far-field flux of this planet is flat only to 7.4 per
  cent, which section 6.2's validity test rejects.
- the post-processing pass runs **without** the `EXHALE_PTC*` variables.
  `finish_case.sh` sets them inline on the JFNK command alone; carried into
  the second pass they make the binary re-enter the steady solver and stop
  before writing the advection-corrected profiles, so `*_adv.txt` and the
  mass-loss rate never appear. This was found by the first LHS 1140 b
  closure run reporting `log10 Mdot = nan`.

## 7. Tests and gates

| id | test | criterion |
|---|---|---|
| T-E1 | schema round trip | a profile written by the adapter, read by the EXHALE reader and re-emitted by the Python loader reproduces every column to double precision; unknown columns are carried, not dropped |
| T-E2 | scalar fallback | with no `Lower atmosphere profile:` key, `make check` is **6/6 byte-identical**; goldens untouched |
| T-E3 | constant K_zz identity | a profile whose `Kzz` column is uniformly `1.0e9`, on the `mol_diffusion` case, reproduces that case's golden bitwise (section 3.2's IEEE argument, verified rather than asserted) |
| T-E4 | pair enforcement | a profile plus a `base.inp` carrying an EOS/reservoir/constraint key stops with the key named; a pair with mismatched `solution_id` stops with both ids printed; a legacy `base.inp` alone still runs |
| T-E5 | overlap flux | `F_H` and `F_He` over the overlap window are flat to the same tolerance as the mass flux (T8 criterion), on a converged LHS 1140 b run with diffusion on -- **measured 2026-08-27, and the answer is that they are not, in the overlap window.** Where that window is non-empty at all it holds 12-18 faces between 1.007 and 1.0098 R_p and the flux there has a radial spread of 13 (F_H) to 35 (F_He), with a median 4.5x the flux the wind actually carries. The criterion is met in the steady-flux window instead, and there it is met exactly as stated: `F_H` 0.44%, `F_He` 0.46%, mass flux 0.45% on the reference arm. Section 3.4 records why the overlap interval cannot do better |
| T-E6 | closure residual | the loop of section 6.1 reaches `eps_k <= 0.05` within `k_max`, and the residual is above neither the flux spread nor the reproducibility floor -- **passes 2026-08-27** on all three LHS 1140 b arms (`k` = 0, 5, 6 against `k_max` = 8); the converged residuals are 0.015, 0.031 and 0.029 against window spreads of 0.004-0.009 |
| T-E7 | initial-value insensitivity | section 6.3: three starts, converged `HeH` within 5%, `log10 Mdot` within 0.025 dex -- **passes 2026-08-27**: `He/H` agrees to 5.5e-6 and `log10 Mdot` to better than the printed 0.005 across starts spanning 10x |
| T-E8 | element budget | `src/utils/element_budget.py <case>` closes for H, He and every element the profile carries, at its default tolerance, against `EXHALE_resolved.out` (the pre-existing He 2^3S double count of `calc_rho`, `docs/Update_EXHALE.md` section 74, is the known exception and is reported, not hidden) -- **closed 2026-08-27** on the LHS 1140 b profile run: C 1.44e-14, N 1.47e-14, O 1.45e-14 against 1e-8. With `He_diffusion` on, He/H is a solved profile and not a column invariant (2.0921 at the base, 0.167 at 30 R_p), so the H and He rows are stated at the base cell, where the reservoir is the boundary condition, and the separation is reported beside them; `EXHALE_resolved.out` now carries `he_diffusion` so the checker knows which of the two it is testing |
| T-E9 | cross-code | the VULCAN adapter run at the same match on the same planet gives an `X_He` and an `F_H` whose difference is stated beside the P1 spread (network 3.95, code 1.70 on a 864 K base; 1.02/1.08 on a 2331 K base) |
| T-E10 | regression case | one new case with the profile route on, snapshotted at the end of the series, so the reader and the K_zz profile path are pinned from then on -- **added 2026-08-27** as `backup/regression/lower_profile`, the seventh case of the default matrix. It is `examples/17_lower_profile` (HD 209458 b, He 2^3S on, `He_diffusion` and `He_metal_diffusion` on, `Solver: Newton`) with a provenance-only `base.inp` beside the profile and a 12000-step cap, the relaxation-snapshot convention `mol_diffusion` uses. One run exercises the reader, the matching-level base state (`T0 = 1450.0 K`, `R0 = 1.401 R_J`, `p_base = 1e-6 bar`, `q_H2 = 0.11961`, `He/H = 0.083333`), the elemental reservoirs the profile carries (`C/H = 2.700e-4`, `O/H = 4.900e-4`, `N/H = 6.800e-5`, so the case is metals-on with no `metals.inp`), `K_zz(p)` interpolated onto the grid (`1.130e9` to `5.623e9 cm^2/s`, not a constant), the accepting branch of the `solution_id` pair enforcement, and both flux windows -- here the overlap window is **not** empty (`r_lo = 1.0032`, `r_hi = 1.0472 R_p`, 135 faces, `lower_profile_flux_state measured`), which on LHS 1140 b it is. `log10 Mdot = 9.58 g/s` at the cap; wall time about 11 minutes single-threaded. Golden snapshotted from that run and the matrix then verified case by case: **7/7 byte-identical**, 14 file comparisons, `wasp_full` re-run from scratch to check the snapshot (`count=14059`, reproduced exactly) |

T-E2 and T-E3 are the gates that let the work proceed incrementally; T-E5
through T-E7 are the plan's acceptance criteria restated as measurements.

## 8. Milestones

| # | content | gate |
|---|---|---|
| **E1** | schema, EXHALE reader, `kzz_cell`, resolved-config provenance, pair enforcement | T-E1, T-E2, T-E3, T-E4 -- **done 2026-08-27** (`docs/Update_EXHALE.md` section 76) |
| **E2** | Photochem adapter, one-way (profile produced from a prescribed T(p), no climate, no iteration) | a hand-checked LHS 1140 b profile drives a converged EXHALE run; T-E8 -- **adapters written and exercised on HD 209458 b, 2026-08-27** (`docs/Update_EXHALE.md` section 77); **both remaining gates closed with E3** |
| **E3** | clima step, tropopause and cold-trap water, `Kzz(p)` from the climate solution | cold-trap `f_H2O` documented and compared with the published estimate; the comparison stated either way -- **done 2026-08-27** (`docs/Update_EXHALE.md` section 78): `radiative_convective_column.py` + `--climate`; `f_H2O = 0.13-0.33 ppm` at a 0.5-1.0 bar tropopause against the published `~7 ppm` at 0.1 bar, stated as a departure and not tuned; the LHS 1140 b profile drives a JFNK `info = 0` run and T-E8 closes on C/N/O to 1.5e-14. `Kzz(p)` is NOT from the climate solution -- it provides none -- and the run states the adopted constant, recorded in the header |
| **E4** | the closure driver and the LHS 1140 b application | T-E5, T-E6, T-E7 -- **done 2026-08-27** (`docs/Update_EXHALE.md` section 79). `src/utils/element_flux_closure.py` + the trial fluxes actually imposed on the chemistry through `set_upper_bc`; the overlap window measured and rejected on physical grounds, the steady-flux window used in its place with the substitution recorded per iteration. T-E6 and T-E7 pass; T-E5 passes in the steady window and fails in the overlap window, which is the result rather than a defect. **The composition is now an output**, and on LHS 1140 b it comes back at `He/H = 2.09235` -- the value eddy mixing holds, not one the escape sets |
| **E5** | documentation and regression | `docs/input_schema.md` section 2c and `README.md` carry the key and the schema; `docs/Update_EXHALE.md` carries the change; T-E9, T-E10 -- **done 2026-08-27**. T-E9 was measured with E2 (`Update_EXHALE.md` section 77). T-E10 is `backup/regression/lower_profile`, now the seventh default case of `run_check.sh`, with the matrix paragraph of the project `CLAUDE.md` and the table in `README_HOWTO.md` updated beside it. The schema lives in `docs/input_schema.md` section 2d, and the reference manual now carries both the key and the file: `docs/EXHALE_user_manual.tex`, "The lower atmosphere as a profile". `docs/Update_EXHALE.tex` was brought level with the markdown through section 79 |

E1 and E2 are independent of Phase D's remaining items; E4 is not — it needs
diffusion on the direct steady route and the restart behavior of
`binary_diffusion_design.md` sections 7.3-7.4.

## 9. Decisions to confirm before implementation

**All six were confirmed as proposed on 2026-08-27 and are implemented as
written, except that the helium-flux arm of decision 5 belongs to E4.**

Expensive to reverse, because each fixes a file format, a user-visible key,
or a physical convention:

1. **File format.** A `#`-header plus fixed-column text table (proposed:
   consistent with every other EXHALE product, readable by
   `examples/exhale_io.py`'s schema-header convention, diffable). The
   alternatives are a YAML/JSON header with an attached table, or a binary
   product. Reversing this after the first LHS 1140 b run means rewriting
   every stored profile.
2. **Key names.** `Lower atmosphere profile:` for the input key,
   `lower_atmosphere_profile.dat` for the default file name,
   `photochem_to_lower_profile.py` / `vulcan_to_lower_profile.py` for the
   adapters, `element_flux_closure.py` for the driver. Also: whether
   `he_kzz` / `He_Kzz` keep their names once the coefficient can be a
   profile, or become `eddy_diffusion_*`. Renaming the key breaks existing
   input files, including the four regression cases' `base.inp` and the
   LHS 1140 b run directories.
3. **Matching-level rule.** Proposed: the match is stated by the producer in
   the header (`p_match_bar`), the file must extend above it, and EXHALE
   places its base there. The alternatives are letting EXHALE choose the
   match (e.g. where `tau = 1` or where `tau_chem = tau_adv`), or fixing it
   at 1 microbar. This choice determines what every stored profile has to
   contain.
4. **Single source versus enforced agreement.** Proposed: with a profile
   present, the scalar `base.inp` keys of the three physics categories are
   **refused**. The alternative is to accept both and enforce numerical
   agreement within a tolerance. The refusal is stricter and simpler; the
   tolerance version is friendlier to hand-editing. This is a user-visible
   error behavior.
5. **Flux update rule and tolerance.** Proposed: damped Picard with
   `omega = 0.5`, halving on failure to descend, undamped residual as the
   test, `tol = 0.05` justified by the 0.025 dex reproducibility floor, and
   `k_max = 8`. Whether the closure iterates the helium flux as well as
   hydrogen, or holds helium at zero and reports its measured value, is part
   of the same decision.
6. **Which code is the production chemistry.** P1 left this to the user
   (`vulcan_photochem_comparison.md` P1.7). This design assumes Photochem for
   production and VULCAN as the cross-check arm, because only Photochem has
   the climate model section 5 requires. If the user decides otherwise, the
   climate step needs an external radiative-convective model and section 5 is
   rewritten.
7. **What a restart does with a reservoir the handoff has moved.** Confirmed
   and implemented 2026-08-27, after the six above: the loaded column of an
   element the handoff states is renormalized onto the stated reservoir by one
   factor per element (section 3.3). It belongs on this list because it is a
   user-visible change of restart semantics -- a restart of such an element no
   longer reproduces its seed exactly. Restarts with no handoff, and elements
   the handoff is silent about, are bit-for-bit what they were.

## 10. Relation to the continuum IR coupling, item (G)

`docs/oxygen_chemistry_new_plan.md` section 4.1 ranks item (G) — the missing
continuum IR coupling — **above** the composition work this phase belongs to,
and for a measured reason: a converged molecular layer radiates itself down
to 190-400 K against `T_eq ~ 1100-1400 K`, and on the HD 209458 b molecular
example the collapse reaches the wind at `Mdot -0.34 dex`. Phase E does not
close it and must not be presented as closing it. What Phase E supplies is
`T(p)` at and *below* the matching level, from a climate model; EXHALE's own
thermal structure *above* the match, through the molecular layer up to the
H2 -> H front, is still set by line coolants with no continuum term, so the
converged temperature there remains a property of the model rather than a
prediction. The two are the same weakness seen from two sides — (G) is that
layer's energy not being constrained, (H) its composition — and closing
either alone does not make the layer a prediction. On LHS 1140 b the concern
is not the hot-Jupiter collapse but its cold-planet counterpart, which is
unmeasured: at `T_eq = 226 K` the layer between the base and the front is
where both items act, and no run yet says how much of the wind's energy
budget passes through it. That statement is a caveat on the Phase E result,
not a reason to defer it.
