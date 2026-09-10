# B1a: the active-equation inventory and the validity-state semantics (2026-09-06)

Step B1a of `docs/PLAN_20260906_rev2.md` section 3, written before A2 and
required by review items R4, R5 (`docs/PLAN_20260906_review.md`), F2 and 5.3
(`docs/PLAN_20260906_review2.md`). It answers three questions with one list:

1. which equations are *independent* in a given run, so that A2 evaluates the
   set the run actually solves rather than a hard-coded pair of species;
2. which of them can be evaluated today on a given state without changing it,
   so that A2's "evaluated / unavailable / not applicable" status has a source;
3. which validity states B6 has to report, who produces each, and what is
   recorded about it now.

The same list generates the D-physical prerequisites of a configuration by
active code path (F2), which is section 5.

## 0. Provenance

Everything here is READ from the source tree at
`/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/` on 2026-09-06,
HEAD `35d9dd5` with uncommitted work in the tree. **Nothing was MEASURED:** no
build, no run, no test was executed while writing this document. Every routine
named was located by grep in this tree.

**The tree is being edited while this was written.** `steady_residual.f90`,
`steady_newton.f90`, `ionization_equilibrium.f90` and
`diffusive_photochemistry.f90` were all modified between 08:59 and 09:01 on
2026-09-06, which is Step A1 landing; `diffusive_photochemistry.f90` moved by
about 256 lines during the reading. That is why a line number is not a citation
here: it is a pointer into a moving file and the routine name is the durable
one. The document was written with line numbers and they had gone stale within
the day; on 2026-09-07 (item DOCS-LINES) every citation below was re-anchored
to the file and the routine, or to the labelled block of a long routine, read
against the live tree. The same was done to
`docs/d0_governing_system_20260906.md`, whose numbers were older still. Where
D0 already states a fact this document cites the D0 section.

Two consequences of that same work are visible in the source and are recorded
here as state, not as a review: `carrier_returned_state_verdict`
(`diffusive_photochemistry.f90`) now exists as an acceptance decision on
the returned carrier state, and `steady_gates_met` has gained a fourth,
optional gate `n_no_chem_root` (`steady_residual.f90`, `steady_gates_met` and
`n_cells_without_chemical_root`). This inventory
describes the equation set A2 must evaluate; it does not describe or endorse
either decision.

**Updated 2026-09-06, after A2 step 4.** Sections 2.3, 2.4 and 2.5 said that
the elemental transport balances, the eliminated-species closure and the
H(n=2) system had no state-consistent evaluator. Each of them now has one, and
the three "State-consistent evaluation today" verdicts below are rewritten to
say what exists, with routine names re-read from the tree on that date.
Section 4.4 is corrected in the same pass: the two temperature floors it
listed as unbudgeted accepted corrections are no longer corrections at all.
Nothing else in this document was re-read on that date, so a statement outside
those four blocks is as old as the rest of the file.

Two source facts contradict the brief that commissioned this document and are
recorded rather than smoothed over:

- `backup/regression/wasp_he23off` is **not** a pure H/He case. It carries a
  `metals.inp` with C, N, O, Mg, Ca, Na, Fe, so it is helium plus ten-element
  trace metals with the triplet off. The molecule-free, metal-free
  configuration in the tree is `examples/tutorial_nometals` (no `metals.inp`,
  `Include He23S? True`, no `Stellar Teff`, so H(n=2) is off as well).
- `backup/regression/oxygen_chemistry` does not state
  `Molecular carrier transport`. The default resolves to
  `carrier_transport = thereis_oxychem` (`input_read.f90`, `input_read`), so that case
  runs four transported carriers, not zero.

---

## 1. Configuration flags

### 1.1 File presence

| file | read by | effect |
|---|---|---|
| `metals.inp` | `read_metals_input` (`metals_input_read.f90`), called from `input_read`; a missing file returns before any key is read | sets `X_C ... X_Fe`, and `cx_full`, `cx_O2p_H`, `cno_cool`, `eos_metals`, `pp_metals`. `thereis_metals` is any `X_* > 0` (`input_read.f90`, `input_read`) |
| `opacity.inp` | `read_opacity_input`, called from `set_energy_vectors` (`set_energy_vectors.f90`) | loads the opacity-model tables of the dispatcher (same routine; the dispatcher covers H I, He I, He II and He 2^3S only, D0 C14) |
| `base.inp` | `read_base_inp` (`input_read.f90`), called from `input_read`, with `apply_lower_atmosphere_profile` and `set_element_abundance` beside it; keys `T_base, r_base, HeH_base, Kzz_base, q_H2_base, p_base` are the labels `read_base_inp` matches | states the base level and the base composition of the lower-atmosphere handoff; a stated `p_base` fixes `n0` |
| `<lap_file>` from `Lower atmosphere profile:` | `read_lower_atmosphere_profile` (`lower_atmosphere_profile.f90`) | replaces the `base.inp` scalars by a table over pressure; carries the elemental reservoirs (metals on with no `metals.inp`) and `K_zz(p)` |

The composition scalars of `base.inp` are refused beside a profile
(`input_read.f90`, `read_base_inp` calling `refuse_scalar_key`), so the two
handoff forms never both state a base.

### 1.2 `input.inp` keys

Full recognized list: the `known_keys` parameter array in the head of
`input_read.f90`. The table below carries the keys that switch physics or a
solver; the purely numeric ones (luminosities, planet parameters, tolerances)
are omitted where they change no equation. "Word" is the token index
`get_word` reads. Every key is parsed in `input_read` (`input_read.f90`), one
long routine, so the citation below is the key label itself where nothing
narrower exists.

| key | variable (default; declared in `parameters.f90`) | what it activates |
|---|---|---|
| `He/H number ratio` | `HeH`, `thereis_He` (`input_read.f90`, `input_read`) | helium stages: `N_eq` 1 -> 3 (same routine, the `N_eq` assignment block) |
| `Include He23S?` | `thereis_HeITR` (`.true.`); parsed in `input_read` | the He 2^3S level row `tr_triplet_row` (`ion_residual_core.f90`), `System_HeH_TR*`, triplet cooling and the 10830 line |
| `metals.inp` presence | `thereis_metals` (`.false.`) | `metal_rows`, `metal_electron_sum` (`ion_residual_core.f90`), metal cooling, `cx_init` |
| `Molecular chemistry` | `thereis_mol` (`.false.`) | `System_HeH_mol*`, `mol_rates`, `molecular_reaction_heat`, `h3p_cooling`, `caloric_eos` H2 ladder |
| `Oxygen chemistry` | `thereis_oxychem` (`.false.`) | oxygen rows `oxygen_carrier_rows` (`System_HeH_mol.f90`), `water_photolysis_init` (called from `input_read`), CO handling, and by default the carrier operator |
| `Molecular carrier transport` | `carrier_transport` (`.false.` but conditional); resolved in `input_read` | `photochemical_transport_step` (`diffusive_photochemistry.f90`), the carrier set of `carrier_set_init` (same file) |
| `Ionization transport` | `ionization_transport` (`.false.`) | H+ as the fifth carrier; `x_hp_fixed` imposed on the sweep for `j >= 1` (`ionization_equilibrium.f90`, `ioniz_eq`). Refuses `Molecular chemistry: False`, `Molecular carrier transport: False`, `Coupled carrier solve`, `Solver: Newton` (`input_read.f90`, `input_read`) |
| `Coupled carrier solve` | `carrier_in_newton` (`.false.`) | `set_transported_species_rows` (`steady_newton.f90`): one row and one unknown per transported balance the configuration activates, and `nvar_jac` with them |
| `He_diffusion` | `he_diffusion` (`.false.`) | `element_diffusion_step` (`binary_element_diffusion.f90`), called from the attempted step of `EXHALE_main.f90` and from `steady_wind_with_element_diffusion` |
| `He_Kzz`, `He_alphaT`, `He_ambipolar`, `He_metal_diffusion` | `he_kzz` 0, `he_alphaT` 0, `he_ambipolar` `.true.`, `he_metal_diffusion` `.false.` | terms of the diffusive flux `J` (`settling_coefficient`, `drift_and_gradient_face_coefficients`) and the trace-metal arm (`solve_trace_element_in_hydrogen`), all `binary_element_diffusion.f90` |
| `Stellar Teff` **and** `Stellar radius` | `use_excited_H = (T_star_eff > 0) .and. (R_star > 0)` (`input_read.f90`, `input_read`) | `excited_H_update` and the whole H(n=2) system (`excited_hydrogen.f90`); the Balmer continuum band; `Spectrum type: Planck` is refused without both (same routine) |
| `Secondary_ionization` | `use_sec_ion` `.true.`, `sec_ion_immediate` `.false.` (`input_read.f90`, `input_read`) | the SvS85 partition inside the photoheating integrand (`util_ion_eq.f90`, `photoelectron_share` and the `photoionization_field_at_cell_H` / `_HHe` integrands); staged flip in the "Staged secondary ionization" block of `EXHALE_main.f90` |
| `Base IR field` | `base_ir_field` (`.false.`) | `fine_structure_line_transfer` incident field (`Cool_coeff.f90`) and `h3p_net_cooling_rate` (`h3p_cooling.f90`) |
| `Molecular IR bands` | `mol_ir_bands` (`.false.`) | `h2_line_emission_lte`, `h2o_band_emission_lte`, `co_band_emission_lte` and their net-cooling functions (`molecular_infrared_cooling.f90`) |
| `Molecular reaction heat` | `mol_reaction_heat` (`.true.`) | `molecular_chemical_heating` and `species_formation_energy` (`molecular_reaction_heat.f90`), summed in `ioniz_eq` (`ionization_equilibrium.f90`) |
| `Stellar LW flux` | `F_LW_star` (0); `lw_flux_stated`; from the SED in `input_read` | `lyman_werner.f90` H2 photodissociation and, with oxygen on, the H2O/OH photolysis of the same beam |
| `Stellar FUV B3/B4 flux`, `Stellar Lya flux` | `F_FUV_B3/B4`, `F_Lya_star` | `water_photolysis.f90`; the Ly-alpha stellar beam of `lya_rt`. `Stellar FUV B1 flux` was retired on 2026-09-06 when its 1110-1201 A band was merged into the Lyman-Werner interval |
| `H2 double ionization` | `h2_double_ionization` (`'chung80'`) | the H+ + H+ channel of `h2_photo_channels.f90` |
| `H2 neutral dissociation` | `h2_neutral_dissociation` (`.true.`) | the H + H channel over 33-41 eV |
| `Spectrum type` | `sp_type`, `do_read_sed`, `is_PL_sed`, `is_monochr` (`input_read.f90`, `input_read`) | which field `J_inc`/`sed_read` builds; `Monochromatic` below 24.6 eV zeroes helium outright (same routine) |
| `Use only EUV` | `thereis_Xray` (`input_read.f90`, `input_read`) | the X-ray band of the photon grid and the secondary-electron budget |
| `Atomic rate set`, `Legacy_HHe_rates`, `ATES_photoionization_rate` | `atomic_rate_set_k22`, `legacy_hhe_rates`, `ates_photoion_rate` | which recombination and collisional-ionization fits the H/He rows use (`Cool_coeff.f90`, the `atomic_rate_set_k22` branches of `alpha_rec_HII_B`, `alpha_rec_HeII_B` and the collisional-ionization coefficients) |
| `He_rec_coupling`, `He_H_charge_exchange` | `use_he_rec_coupling` `.true.`, `he_h_charge_exchange` `.true.` | the on-the-spot He recombination photon term; `he_h_cx_fvec` (`charge_exchange.f90`) |
| `Deexc heat` | `incl_deexc_heat` (`.true.`) | the 10.2 eV de-excitation return in the `heat_balmer` assembly of `excited_H_update` (`excited_hydrogen.f90`) |
| `Caloric EOS` | `caloric_eos_monatomic` (`.false.`) | the Roueff H2 ladder of `caloric_eos.f90` versus a constant index |
| `Photoelectron heating` | `photoheat_full_photon_energy` (`.false.`), `photoheat_photon_fraction` | `h nu - I` against the whole photon in the heating integrand |
| `Viscosity`, `Conduction` | `visc_on`, `cond_on` (both `.false.`) | `viscous_conduction_step` (marching) and `viscous_conduction_sources` (steady residual) |
| `Shapiro filter` | `shapiro_eps` (`-1`), `shapiro_every` (4) | `shapiro_filter`, called in the attempted step of `EXHALE_main.f90` after the energy stage and the transport stage |
| `Low-Mach damping` | `lowmach_damp_eps` (`-1`) | the JST fourth difference of `low_mach_dissipation.f90` |
| `Domain mode`, `Outer radius` | `spherical_domain` (`.false.`), `r_out_user` | spherical potential in place of the Roche potential (`grav_field.f90`, `phi` and `Dphi`) |
| `Base BC` | `base_bc_mode` (0 = density) | whether `n0` is stated or derived from `p_base` (`input_read.f90`, `input_read`) |
| `Energy solver` | `use_semi_implicit_energy` (`.true.`) | `solve_energy_semi_implicit` in the `coupled_source` loop of `EXHALE_main.f90` versus the explicit forward-Euler branch beside it |
| `Time stepping` | `use_local_dt` (`.false.`) | the cell-wise `dt_loc` of `eval_dt` (`eval_dt.f90`), which is then the step of the diffusion, carrier and energy operators |
| `Solver` | `use_newton_solver` (`.false.`), `newton_du_switch` | the JFNK hand-off of `EXHALE_main.f90`, the `newton_du_switch` block after the marching stop test |
| `IC mode` / `Transonic IC` / `Hot Parker IC` / `Load IC` | `ic_mode`, `transonic_ic`, `hot_parker_ic`, `do_load_IC` | which state the run starts from; `windae` runs the Wind-AE port |
| `Newton solver`, `Brent solver` | `use_newton_ieq`, `use_brent_tsolve` (both `.true.`) | the cell ionization solver and the post-process temperature root finder |
| `Molecular base`, `Lower atmosphere`, `Lower column` | `molecular_base`, `lower_atm_mode` (0), `lower_col_r1bar` | the analytic Koskinen column and the molecular base particle count |
| `Grid cells`, `Base grid`, `Grid type`, `Reconstruction scheme`, `Numerical flux`, `CFL`, `du_th`, `Reconstruction continuation` | grid and discretization | they change the discrete operator, not the equation set |

---

## 2. Independent equations per flag combination

Notation: an equation is **independent** if the run carries an unknown for it
that no other equation determines. An equation is a **closure** if it
determines a quantity that is then substituted, so that its residual is an
accuracy statement about an eliminated variable rather than a row of the
system.

### 2.1 Always present: the hydrodynamic rows

Three rows on the physical cells, `u = (rho, rho v, E)`
(`EXHALE_main.f90`, the state allocation of `allocate_state_vectors`;
`UW_conversions.f90`, `W_to_U` and `U_to_W`), D0 section 2.1.

| row | continuous form | marching operator | stationary residual |
|---|---|---|---|
| mass | `d_t rho + (1/r^2) d_r(r^2 rho v) = 0` | the RK3 stages of the attempted step (`EXHALE_main.f90`, operator tag `as_op_hydro`) | `R(:,1) = dF - S` (`steady_residual.f90`, module header and `assemble_residual`) |
| momentum | `+ d_r p + rho d_r Phi`, plus the viscous stress when armed | same, plus `viscous_conduction_step` in the same block | `R(:,2) = dF - S - F_mu` |
| energy | `d_t E + (1/r^2) d_r(r^2 v(E+p)) = heat - cool + viscous/conduction` | the RK3 flux carries the work; sources enter in the `coupled_source` loop and in `viscous_conduction_step` | `R(:,3) = dF_E - S_E - (heat - cool) - (w F_mu + q_mu + conduction)` |

Two facts about the mass row that A2 must not lose (D0 C2): the chemistry stage
overwrites `rho` and pushes it into the conserved row, and the stationary
residual has no such term. **Both were removed on 2026-09-06 (item B3c,
`docs/Update_EXHALE.md`):** `ioniz_eq` now takes the density `intent(in)`
(`ionization_equilibrium.f90`, the `n_io` declaration and its T2.1 note), the
mass sum of the returned composition comes back as the optional check
`rho_recon` built by `calc_rho`, and the pressure rebuild that pushed it into
the energy row is gone with the composition projection. The marching fixed
point and the residual zero remain different objects wherever `chem`,
`transport` or `filter` of `EXHALE_UPDATE_MAP` is non-zero
(`EXHALE_main.f90`, `update_map_end_step`).

**State-consistent evaluation today:** yes, with a caveat.
`assemble_residual(u, n_part, heat, cool, R)` (`steady_residual.f90`)
is a pure function of its arguments *given* `heat` and `cool`, which the caller
must have obtained by running `ioniz_eq` on the state first (the routine's own
header). `ioniz_eq` is not side-effect free: it rewrites `f_sp` (and, until
item B3c on 2026-09-06, `rho`). So a
certification of the hydrodynamic rows on an already-solved state is available,
but only because the caller supplies the source arrays; A2 cannot obtain them
without either accepting the existing arrays or re-running a mutating sweep.

### 2.2 Transported carrier balances

`carrier_set_init` (`diffusive_photochemistry.f90`) fixes the set once,
after every key is parsed:

| carrier | index | active when |
|---|---|---|
| H2 | `ic_H2` | always, when the operator runs at all |
| OH, H2O, CO | `ic_OH, ic_H2O, ic_CO` | `thereis_oxychem` |
| H+ | `ic_Hp` | `ionization_transport` |

Equation (`diffusive_photochemistry.f90`, module header), for each active
carrier:

```
d_t n_i + (1/r^2) d_r [ r^2 (n_i v + Phi_i) ] = P_i - L_i
Phi_i = -n_tot (D_i + K_zz) d_r f_i - n_i D_i [1/H_i - 1/H_atm]
```

Backward Euler, as the module header's discretization note states, assembled by
`carrier_residual` and solved by a block-tridiagonal Newton (`solve_carriers`),
all in `diffusive_photochemistry.f90`.
R4's point restated in the code's own terms: the steady form of one of these is
transport minus reaction, `div(n_i v + Phi_i) = P_i - L_i`, and
`carrier_steady_residual` builds exactly that by setting `dt_big = 1e30`
to kill the time term. Requiring `P_i - L_i = 0` would be a different
and wrong equation.

**CO was the exception inside the exception, and is no longer.** As audited it
was a transported reservoir with **no balance row and no destruction rate**
(D0 C26), the chemistry taking it as transported and
`limit_to_element_budget` applying a thermal equilibrium ceiling with no rate
behind it. **The ceiling was deleted on 2026-09-06 (item CEILING-DEL)** and the
CO row now carries published rates: `carrier_source`
(`diffusive_photochemistry.f90`) destroys CO by He+ + CO -> C+ + O + He
(UMIST RATE22 4068) and by photodissociation of the 912-1201 A beam with the
Visser et al. (2009) shielding function. The chemistry still takes the
transported CO rather than the local equilibrium value
(`ionization_equilibrium.f90`, `ioniz_eq`), and what
`limit_to_element_budget` (`diffusive_photochemistry.f90`) applies to
it is conservation only: a cell's carriers cannot hold more carbon or oxygen
nuclei than the cell has. The row is destruction only, so the region in which
that is legitimate is measured cell by cell and reported rather than assumed
(`carrier_co_domain_take`, same file).

**State-consistent evaluation today: no.** `carrier_steady_residual`
(`diffusive_photochemistry.f90`) is the nearest thing and it fails both halves
of R5:

- it sets `rcmax = 0`, and every optional output to zero, then returns when
  `.not. thereis_mol`, `.not. carrier_transport` or `.not. bg_ready`. The third
  is an active model whose background is not available, and it is
  indistinguishable at the call site from a converged zero;
- its arguments are `intent(in)` and it is not side-effect free. It writes the
  module-save arrays `row_terms`, `col_scale_car`, `row_scale_car`,
  `headroom_car` and the flag `headroom_set`, and it calls
  `carrier_photolysis` and `carrier_diffusivities`, which refresh
  further internal arrays. The headroom refresh is even conditional on
  `weno_mode .eq. 1`, so the value left behind depends on which
  evaluation context called it.

What else behaves like that (every routine in `diffusive_photochemistry.f90`):

| routine | initialized-zero early return | mutates |
|---|---|---|
| `carrier_steady_residual` | `thereis_mol`, `carrier_transport`, `bg_ready` | `row_terms`, scales, headroom, photolysis, advection correction |
| `carrier_transport_interval`, the body of `photochemical_transport_step` | `bg_ready` | the whole carrier state, the CO record, the limiter diagnostics |
| `relax_photochemical_composition` | `bg_ready` | (same module block) |
| `carrier_transport_diagnostics` | none, but reports the **last** solve, not the state handed in | nothing |
| `carrier_cell_is_constrained` | returns `.false.` before any solve | nothing |

`bg_ready` and `bg_cell` are owned by `ionization_equilibrium` and imported by
the `use` statement in the head of `diffusive_photochemistry.f90`, so "the
backgrounds of the state being
certified" is not something the caller can assert today: it is whatever the
last sweep left.

### 2.3 Elemental transport balances

`binary_element_diffusion.f90` (2357 lines when this was written; 3037 on
2026-09-07, after B4-1 and B4-1c), active under `He_diffusion`.

- **He/H partition**, unknown `X = rho_He/rho` (stated in the module header and
  set by `element_diffusion_step`), solved in
  advective form `rho (X^{n+1}-X^n)/dt + (1/r^2) d_r(r^2 J) + rho v d_r X = 0`
  by `solve_mass_fraction`, one backward-Euler step per call with frozen
  coefficients. Conservative form is in the module header and is not what
  is solved; helium mass is conserved to the hydro truncation in a transient
  (D0 C28). Since B4-1 (2026-09-06) that advective term is taken only when
  `solve_mass_fraction` is called with `advect` true, the fixed-wind
  relaxation: in the marching route the advection rides on the hydro's face
  mass fluxes inside the Runge-Kutta stages and the row is the diffusive half
  alone.
- **Trace elements**, under `He_metal_diffusion`
  (`solve_trace_element_in_hydrogen`): each metal mixing ratio moves against a
  fixed hydrogen background with **no counter-flux**, so zero net diffusive
  mass flux, which is structural for the H/He binary (stated at
  `helium_hydrogen_diffusion` and in the header), does not hold for them
  (D0 C29).
- The inner boundary is a Dirichlet reservoir, the outer a zero gradient
  (`element_diffusion_step`, the boundary block that sets `X_base` and `jlo`).
- In the steady route the operator is outside the residual entirely and is
  co-converged by a Picard alternation with under-relaxation
  (`EXHALE_main.f90`, `steady_wind_with_element_diffusion`), because
  `steady_residual` contains no diffusion (D0 C17).

**State-consistent evaluation: yes, since A2 step 4.**
`element_transport_residual` (`binary_element_diffusion.f90`,
public) returns the STATIONARY balance of both arms on a state handed in:

- the He/H partition, `div(r^2 J)/r^2 + rho v dX/dr = 0` (no source, no sink),
  obtained by evaluating the operator's own `composition_residual` at a step
  long enough that the time term is numerically absent, which is the device
  the carrier balance uses;
- each trace element's mixing-ratio balance, in the same form, with the
  element furthest out in units of its own terms reported at each cell so that
  an element whose rates are small beside another's is not hidden by it.

It returns the residual and the scale cell by cell, in `g cm^-3 s^-1` for
helium and `s^-1` for a trace element; the scale of a row is the sum of the
magnitudes of its own terms and the floor beneath it is a rate 1e-20 of the
base composition carried across the domain in one flow time `R0/v0`. Cell 1
carries no equation in either arm (it is the Dirichlet reservoir the operator
states) and reads zero against that floor.

`composition_residual` is still private but now returns its
residual and its scale cell by cell through optional arguments, instead of the
maximum alone. Two blocks were factored out of `element_diffusion_step` so
that the operator and the residual share one definition of the coefficients
rather than keeping two: `trace_element_transport_coefficients`,
one element's density, mean charge, stage-resolved diffusion
coefficient and settling coefficient, and `trace_face_coefficients`,
the Peclet-hybrid face coefficients of the trace flux.

**Side effects: none, and it is asserted.** The five diagnostics of the last
solve that the operator writes (`he_fraction_over_one`,
`he_fraction_under_zero`, `he_fraction_newton_steps`,
`he_fraction_newton_resid`, `trace_ratio_under_zero`) are held aside and put
back by `reinstate_diffusion_diagnostics`, which compares them
with the copy after the restore and prints a warning on any mismatch. Nothing
is re-solved and nothing is projected. `report_step` (only under
`EXHALE_DIFFUSION_CHECK=1`) is untouched and remains a diagnostic of the last
solve rather than of a supplied state.

### 2.4 Species eliminated by equilibrium, and the closure that eliminates them

`ioniz_eq` (`ionization_equilibrium.f90`) re-solves a local chemical
and ionization equilibrium in every cell every step; call sites in
`EXHALE_main.f90` (the startup equilibration, the `coupled_source` loop of the
attempted step, `equilibrate_loaded_composition`, `physical_handoff_check`,
`update_map_begin_step` and `steady_wind_with_element_diffusion`) and in
`eval_residual` (`steady_newton.f90`). The system is `N_eq` rows in `N_eq`
stage fractions (`input_read.f90`, `input_read`, the `N_eq` assignment block),
with the neutral stage of each element and the electron density eliminated
(D0 section 4.7):

| system | routine | `N_eq` | independent unknowns |
|---|---|---|---|
| H only | `ion_system_H` (`System_H.f90`) | 1 | `x_HII` |
| He + H | `ion_system_HeH` (`System_HeH.f90`) | 3 | `+ x_HeII, x_HeIII` |
| `+` triplet | `ion_system_HeH_TR` (`System_HeH_TR.f90`) | 4 | `+ x_HeITR` |
| `+` metals | `ion_system_HeH_metals` (`System_HeH_metals.f90`) / `ion_system_HeH_TR_metals` (`System_HeH_TR_metals.f90`) | `3+2n_melem` / `4+2n_melem` | two stages of each metal element |
| molecular | `ion_system_HeH_mol` (`System_HeH_mol.f90`) | 7 (+1 triplet) | `+ x_H2, x_H2p, x_H3p, x_HeHp` |
| `+` oxygen, `+` metals | `ion_system_HeH_mol_metals` (`System_HeH_mol_metals.f90`) | `+2 +2n_melem` | OH and H2O at `oxygen_row_base()`, metals at `metal_row_base()` |

The closures, and what makes each one a closure rather than a row:

- **element totals** eliminate the neutral stage: `n_HI = (1 - sum x) n_H` and
  `n_HeI = (1 - x2 - x3) n_He - n_HeHp` in `ion_system_HeH_mol`
  (`System_HeH_mol.f90`), `nm0(e) = (1 - x - x') n_X` in `metal_fractions`
  (`ion_residual_core.f90`);
- **charge neutrality** eliminates `n_e` (`ion_system_HeH_mol`;
  `calc_ne` in `utilities.f90`). It is never an unknown and never a row;
- **total mass is not an equation of the cell solve**: `rho` was recomputed
  after the solve by `calc_rho` (`ionization_equilibrium.f90`, `ioniz_eq`) and
  the fractions written back as `f_sp = n/n_io`. Since item B3c (2026-09-06)
  only the write-back remains: `calc_rho` builds the optional `rho_recon`,
  which is a check against the density handed in and is never applied;
- a **transported** quantity replaces its balance row by an identity row after
  the row scaling: `fvec(4) = x(4) - x_h2_fix` and
  `fvec(1) = x(1) - x_hp_fix`, both at the end of `ion_system_HeH_mol`
  (`System_HeH_mol.f90`), with the same block in `ion_system_HeH_mol_metals`.
  So a carrier is not eliminated by
  equilibrium; its row is imposed;
- when no root exists in a molecular cell the **constrained continuation** is
  promoted (`ioniz_eq` calling `equilibrium_from_molecular_limit`,
  `constrained_chemical_equilibrium.f90`): the unknowns become
  `u_k = ln n_k`, element conservation becomes explicit rows
  (`element_conservation_rows`) and charge neutrality is an identity at every
  iterate (`network_balance_rows`). The
  roles of row and constraint invert, which A2 must know, because "the
  eliminated species" is a different set in this branch.

Fast-process justification is stated in one place only: O(1D) is folded into
one channel because it has a single sink and a 1.6e-3 s lifetime
(`System_HeH_mol.f90`, the note above `oxygen_carrier_rows`). Every other
elimination rests on the local equilibrium assumption itself, which the code
violates knowingly across the H2 front when the carrier operator is off
(`input_read.f90`, `input_read`, the `carrier_transport` resolution block).

**State-consistent evaluation: yes, since A2 step 4.** What was missing was
never the residual: `normalized_reaction_residual`
(`ionization_equilibrium.f90`) is the measure the sweep judges a
root by, and the acceptance classes are recorded from it (section 4.1 below).
What was missing was the CELL STATE to evaluate it at outside the sweep. As
the sweep passes each cell it now records the rate state it built for that
cell (`ieq_rate_cell`, a copy of `ieq_cell`, declared in the module head), the
electron density the
turnover scale is referred to, the temperature, the particle count and the
metal coefficients in the order `set_metal_coeffs` takes them. Those arrays
are written by the sweep and never read by it, so no solve path depends on
them.

- `ionization_closure_residual_profile(rho, f_sp, res, res_triplet, ok, why)`
  reinstalls each cell in turn and returns the closure residual
  of every physical cell, plus the He 2^3S level row of the same evaluation
  (section 2.5). The composition is read from the CALLER'S `f_sp`, not from
  the sweep, so a state whose composition has moved since is measured as it
  stands; `ieq_state_vector` builds it in the fraction layout
  of the system the run solves, against the element totals the cell's rate
  state carries.
- `ionization_closure_residual_cell(j, x, nx, res, rowres)` is
  the same measurement for one cell at an ARBITRARY composition, which is what
  a perturbation test needs.
- `normalized_reaction_residual` gained an optional `row_out`, the same
  measure row by row, and an H-only branch: it served the He-branch acceptance
  only, so `System_H`'s single balance had no measure at all.

Because the measure is the sweep's own, a cell of acceptance class 1, 2, 3 or
5 reads at or below `ieq_res_tol` and a class-4 cell above it. Nothing is
solved, and `f_sp`, `rho` and every ledger and streak are untouched.

**Side effects: none, and it is asserted.** The scratch of one cell that the
installation writes -- `ieq_cell`, the `met_*` arrays of
`System_HeH_metals`, and the molecular and oxygen coefficients and row
turnover scales of `System_HeH_mol` -- is held aside in one value
(`ieq_closure_scratch`, declared in the module head) by
`save_ieq_closure_scratch`, put back by `restore_ieq_closure_scratch`, and
`ieq_closure_scratch_matches` compares the live scratch with
the copy element by element after the restore, printing a warning on any
mismatch. That check runs on every closure evaluation of every run. One item
cannot be asserted from outside its module and the code says so at the site:
`cx_kc`, the charge-exchange rate coefficients of a cell, is private to
`charge_exchange`, which exports `cx_set_cell` to write it and no accessor to
read it.

**Status semantics.** The closure is `unavailable`, never a number, when no
sweep has left a cell state (`ieq_rates_ready` false) or when the ELEMENT
TOTALS of the state handed in are not the ones those rates were built with --
a different gas, whose closure these rates cannot state.

### 2.5 Level equations

- **He 2^3S** is an unknown of the coupled system when the triplet is on:
  row `tr_triplet_row` (`ion_residual_core.f90`), which carries
  photoionization of the metastable, recombination into it, `q13`, `q31a`,
  `q31b`, electron-impact ionization, `A31` and the **total** He(2^3S)+H
  ionization `Q31 n_HI` (both Penning and associative branches quench it, as
  the note above that routine states). It is a level inside He I:
  `bsp_is_excited_level` (`species_table.f90`) is true for it alone, every
  budget sum skips that column, and the ground singlet is formed by difference
  in one place (`composition.f90`, `he_ground_singlet_density_column`).
- **H(n=2)** is **not an unknown of any solve**. `n2_populations`
  (`excited_hydrogen.f90`) solves a 2x2 2s/2p statistical equilibrium
  outside the ionization system, once each step, before it (the
  `excited_H_update` call that precedes `ioniz_eq` in the `coupled_source`
  loop of `EXHALE_main.f90`), and its products enter as a lagged added rate
  and a lagged added heat (`ionization_equilibrium.f90`, `ioniz_eq`).
  So its "equation" is a closure evaluated one outer pass behind the state it
  is a closure for. It is also not an absorber in `tauE`, so it removes no
  photons while its band neighbors do (D0 C6).

**State-consistent evaluation: yes for both, since A2 step 4.**

- **He 2^3S** is one row of the coupled system, so it comes out of the same
  evaluation as 2.4: `ionization_closure_residual_profile` returns it
  separately through the `row_out` vector of
  `normalized_reaction_residual`, at the row index
  `ieq_triplet_row()` (`ionization_equilibrium.f90`; x(4) in the
  atomic layouts, x(8) in the molecular one). It is reported beside the
  closure it belongs to instead of being folded into the maximum, so a refusal
  can say whether it is the level that is out.
- **H(n=2)** has an interface of its own:
  `excited_hydrogen_level_residual(T_in, n_in, f_sp_in, res, scale, ok, why)`
  (`excited_hydrogen.f90`) forms
  `r2p = L2p n2p - M12 n2s - S2p` and `r2s = L2s n2s - M21 n2p - S2s`
  [cm^-3 s^-1] for the populations the module holds, at the field it holds
  (`Jlya_arr` and the Balmer-continuum rates, which are the field this state's
  own lagged coupling used), and at the composition and temperature of the
  state handed in. `n2s_arr`, `n2p_arr`, `gph_balmer_HI` and `heat_balmer` are
  read and left as they are. The cell's value is the row of the two that is
  furthest out in units of its own terms; the scale is the sum of that row's
  term magnitudes, with a numerical floor of 1e-30 cm^-3 s^-1.

  The rate matrix and source vector are now `n2_rate_matrix`, the one
  definition of those coefficients: `n2_populations` is that call plus
  Cramer's rule, and the residual is that call plus the two rows, so the solve
  and its measure cannot drift apart.

  **This measures the lag that section 2.5 describes.** MEASURED on
  `wasp_full`: the row reads 1.0 -- the residual equal to the sum of its own
  terms -- at cell 435 at marching step 500, and 1.8e-16 at the final state,
  where `EXHALE_main` refreshes `excited_H_update` from the converged state
  before writing. The 1s excitation rates are exponential in T, so a step's
  temperature change moves the source vector by order unity while the
  populations still belong to the previous pass.

**Status semantics.** The H(n=2) row is `unavailable`, never a zero, until
`excited_H_update` has run once (`excited_H_field_ready`, declared in the head
of `excited_hydrogen.f90`): before that there is no field and no population,
and a residual would be a statement about zeros.

### 2.6 The energy and temperature equation, and where it is solved

Four different solves of the same physical balance exist (D0 section 5.6):

| context | equation | routine | acceptance |
|---|---|---|---|
| marching source stage, default | `e(T_new) - e(T_old) = dt (heat - cool(T_new))` divided by `(n_tot + n_e)` | `solve_energy_semi_implicit` (`energy_semi_implicit.f90`), called in the `coupled_source` loop of `EXHALE_main.f90` | **none when this was written**: exactly two Newton iterations, cooling re-evaluated once at the end of iteration 1, no residual test (D0 C3). Replaced on 2026-09-06 (item B2): the solve is bracketed and residual-controlled by `energy_balance_init` / `energy_balance_update` / `energy_balance_finish` in the same file, and reaching the lower bracket end is a `FLOOR` failure, not an adopted state |
| marching source stage, `Energy solver: Explicit` | `u(3) += dt (heat - cool)` | the explicit branch beside it in the same loop | none |
| stationary | `R(:,3) = 0` | `assemble_residual` (`steady_residual.f90`) | `steady_gates_met` (`steady_residual.f90`) |
| post-process `_adv` | steady advected balance in one cell, upwind | `T_equation` (`T_equation.f90`), called from `post_process_adv` | a bracketing scan then Brent (`solve_T_brent`, `brent_root`) |

The post-process solve uses a **different channel set** from the hydro: no H3+,
no molecular infrared bands, molecular ions omitted from `n_e`, and
on the heating side only photoheating, the He recombination coupling and the
He(2^3S)+H Penning heat (`post_process_adv.f90`, `post_process_adv`, its
heating and electron-density assembly). That is D0 C8, and it
matters to A2 because the `_adv` files are what the transit tools read.

Between the sweep and the energy stage sat the composition projection that F2
names: solve the composition at fixed `T`, recompute the particle count,
rebuild `p` from the **same** `T`, `W_to_U`. In the monatomic limit that moved
`u_th` by `(3/2) k_B T Delta n_part` with no source behind it (D0 C1), and
there was **no molecular or metal guard** on the sequence, which was F2's
point: a pure H/He run with changing ionization ran it every step. **Removed on
2026-09-06 (item B3c, `docs/Update_EXHALE.md`):** rows 7 to 9 are now one local
source step per cell at fixed volume, iterated to a fixed point, with the
energy row anchored on `u_th_old` and no `comp_p_from_T` and no `W_to_U` in the
loop (`EXHALE_main.f90`, the `coupled_source` loop and the note above it).

---

## 3. Constraints

| constraint | form | imposed or checked | where |
|---|---|---|---|
| element totals, cell solve | neutral stage by difference | **imposed**, by construction of the fraction layout | `ion_system_HeH_mol` (`System_HeH_mol.f90`); `metal_fractions` (`ion_residual_core.f90`) |
| element totals, constrained branch | explicit residual rows | **imposed as equations** | `element_conservation_rows` (`constrained_chemical_equilibrium.f90`) |
| charge neutrality | `n_e = sum Z n_ion` | **imposed** (identity, never a row) | `calc_ne` (`utilities.f90`); `ion_system_HeH_mol` (`System_HeH_mol.f90`); identically at every iterate in the constrained branch (`network_balance_rows`) |
| fractions inside the simplex | `x >= 0`, each element's tracked stages summing to at most 1 | **checked after the solve**, and a violation at or below 1e-6 is **clamped** | `ionization_fractions_physical`; `clamp_fractions_to_element_budget`, gated inside `ioniz_eq` (all `ionization_equilibrium.f90`) |
| carrier element budget | carriers rescaled onto the element simplex; conservation only, the CO ceiling having been deleted on 2026-09-06 (item CEILING-DEL) | **imposed after the solve**; above headroom a trial is **refused** rather than clamped in the coupled solve | `limit_to_element_budget` (`diffusive_photochemistry.f90`); refusal in `eval_residual` (`steady_newton.f90`, the `headroom_ok` block) |
| element ratio `n_El/n_H` and absolute nucleus density | two invariants | **checked, never corrected**; off unless `EXHALE_ELEMENT_ASSERT` is 1 or 2 | `element_census_verify` and its 1e-9 gate (`element_census.f90`), call sites in `EXHALE_main.f90`, `ionization_equilibrium.f90` and `diffusive_photochemistry.f90` |
| charge neutrality as a diagnostic | computed beside the element census | **reported, not gated** | `element_census.f90`, the census type declaration |
| total mass | not an equation of the cell solve; `rho` recomputed from the accepted composition | **neither**: it was an output of the chemistry until item B3c (2026-09-06) made the density `intent(in)` and the mass sum a check | `ionization_equilibrium.f90`, `ioniz_eq` and `calc_rho` |
| helium mass under diffusion | zero net diffusive flux structural for the binary; **not** for trace metals | measured and reported (`mass closure`), not gated | `binary_element_diffusion.f90`, `helium_hydrogen_diffusion`, `element_diffusion_step` and `solve_trace_element_in_hydrogen` |

### Boundary equations by cell class

| class | cells | equations |
|---|---|---|
| base ghosts (the external reservoir) | `1-Ng .. 0` | characteristic face condition at `r_edg(0)`: the LODI `C^-` compatibility relation plus the `(p, s)` reservoir, with the count set by the eigenvalue signs, `0 < v < c` giving two reservoir conditions and one interior relation (`base_boundary.f90`, module header and `characteristic_base_face_state`). The ghosts written are volume averages of the hydrostatic isentrope through the face state, not copies of it (`base_ghost_averages`, called from `Apply_BC_W` in `Apply_BC.f90`). Supersonic inflow over-specifies the boundary and is **diagnosed only** (`check_base_inflow_is_subsonic`, `Apply_BC.f90`, called from the attempted step of `EXHALE_main.f90`) - D0 C33 |
| base ghosts, composition | same | H2 gets a Dirichlet inflow composition applied as an advective inflow term when the base face flux is inflowing and a handoff states the partition (`diffusive_photochemistry.f90`, the `base_dirichlet` and `base_inflow` declarations and their note, applied in `carrier_source`); face 0 carries no diffusive flux (same note); the proton deliberately gets **no** Dirichlet value (`carrier_set_init`, "the proton is never stated"); the imposed proton partition is applied only for `j >= 1`, so the two ghosts re-solve their ionization locally (D0 C18, `ionization_equilibrium.f90`, `ioniz_eq`); the diffusion operator has its own inner Dirichlet reservoir (`binary_element_diffusion.f90`, `element_diffusion_step`) |
| physical cells | `1 .. N` | the rows of sections 2.1 to 2.5. In the stationary solve only these are unknowns; `set_base_fix` can additionally anchor the first `nfix` of them by replacing their rows with `Y - Yfix` (`steady_newton.f90`, `set_base_fix` and `apply_base_fix`) |
| outer ghosts | `N+1 .. N+Ng` | free outflow: zero-gradient copy, with linear extrapolation added under WENO3 as an accuracy device, dropped back to zero gradient when it leaves `rho > 0, p > 0` and **counted** (`free_outflow_ghost`, `Apply_BC.f90`) |

---

## 4. Validity states (B6 semantics)

Five states, per review 2 section 5.3. For each: what produces it, how it is
recorded today, what A2 needs.

### 4.1 Active unvalidated physics

| item | producer | recorded today | A2 needs |
|---|---|---|---|
| **CO ceiling active** (retired) | the thermodynamic ceiling that produced this state, in `limit_to_element_budget`, was **deleted on 2026-09-06 (item CEILING-DEL)**: CO is now destroyed by rates in the carrier row, and `limit_to_element_budget` applies conservation only. This state therefore has **no producer** | nothing. The destruction model's domain record (`carrier_co_domain_take` / `carrier_co_domain_record`, `diffusive_photochemistry.f90`) replaced it in every consumer: the end-of-run report, the seven `co_domain_*` keys of the resolved-config file (`write_setup_report.f90`, `write_resolved_config`) and the `dom_*` header block of `output/Oxygen_chemistry.txt`. It is reported under 4.2, not here | nothing; `certification.f90` reports `n_active_unvalidated_physics = -1`, the report's own "not produced" convention (the report type declaration and `read_validity_states`) |
| **oxygen energy ledger incomplete** | the oxygen cycle has photolysis heat but no collisional reaction-energy ledger; the retained O(1D) excitation energy is absent; the He(2^3S) associative branch heat (8.1 eV) is deposited nowhere (`ionization_equilibrium.f90`, the heating assembly of `ioniz_eq`; `molecular_reaction_heat.f90`, module header) - D0 C25 | **nothing**. There is no flag, counter or header field | a static "oxygen active" implication: if `thereis_oxychem` then the energy ledger is incomplete by construction, so the state is unconditional for that configuration |
| **H3+ model domain** | `h3p_cooling.f90`: the collider density is clamped to the first table column so the emission does not vanish as `n(H2) -> 0`, and the published piecewise fits join discontinuously - D0 C24. **Both were replaced on 2026-09-06 (item B3b-H3):** below the table's 1e6 cm^-3 `h3p_nonlte_factor` follows the exact collisional limit, linear in `n(H2)` and zero at zero, and `h3p_emission_lte` blends the four Table 5 segments over the 5 per cent below each join | **nothing** when this was written; the clamp was silent. Four domain counters (`h3p_n_below_fit_T`, `h3p_n_above_fit_T`, `h3p_n_outside_nonlte_T`, `h3p_n_below_collider`) now exist and are the record | a "clamped" signal by cell, or a static implication from `thereis_mol` plus the local `n(H2)` range |
| **H2 caloric/chemistry state mismatch** | the Roueff ladder of `caloric_eos.f90` versus the H2 partition function the chemistry used - D0 C22. **Closed on 2026-09-06 (item B3b-H2Q):** one partition function, `h2_partition_function` (`caloric_eos.f90`), gives `u_rv(T)` as its first moment and the chemistry's `K_eq` as its free energy (`mol_rates.f90`, the H2 dissociation equilibrium block) | nothing | static implication from `thereis_mol` |
| **metal Voronov transcription** | the code's own standing "VERIFY against Voronov 1997, Table 1" note (`Cool_coeff.f90`, above the metal collisional-ionization coefficients) | a source comment | static implication from `thereis_metals` |

### 4.2 Out-of-domain closure

Every closure with a stated domain, and where the domain is checked:

| closure | stated domain | check |
|---|---|---|
| He recombination into 2^3S / 1^1S (Oklopcic and Hirata 2018) | 5e3 to 2e4 K, the only H/He recombination fit with a stated range | **not enforced** (D0 section 4.3) |
| Penning He(2^3S)+H (Taylor 2025) | 200 to 10000 K, below 2 per cent | evaluated unclamped outside, accuracy claim withdrawn at the code site |
| CHIANTI metal-line fits | 1e3 to 1e5 K (`Cool_coeff.f90`, the validity note above the `cool_*_chianti` block) | a smooth coronal cutoff below 1000 K, `exp(-x^2)`, `w` from `Coronal cutoff width` (`coronal_excitation_cutoff`); the legacy branch is deliberately unguarded (stated in the same note) |
| Ca II two-level, C I fine structure | fitted range | `T` clamped to `[Tlo, Thi]` (`Cool_coeff.f90`, `fine_structure_upsilon` and `h_impact_deexcitation`) |
| H2O and CO infrared bands | 50 to 2000 K | clamped outside (`molecular_infrared_cooling.f90`, the module header's cross-section range note and `molecular_infrared_init`) |
| H2 line sum | unbounded in `T` at the code site | none |
| H3+ Miller 2013 fits | 30 to 5000 K in `T`, tabulated `n(H2)` | `T` clamped at 30 K and at 5000 K (`h3p_cooling.f90`, `h3p_emission_lte`), collider density held at the table edge above 1e14 cm^-3 and taken to the exact collisional limit below 1e6 (`h3p_nonlte_factor`, since item B3b-H3 on 2026-09-06); each of the four is counted |
| H2 self-shielding table | 700-3200 K, 1e12-1e14 cm^-3, 1e12-5e21 cm^-2 | clamped, stated in the module header of `h2_self_shielding_table.f90` and carried by `h2_shield_interp`, `h2_lw_pump_cross_section` and `h2_shield_locate` |
| Blanc's law mixture | exact for a trace species | stated in `stage_mixture_diffusion` (`binary_element_diffusion.f90`); not checked |
| Chapman-Enskog first approximation | the whole diffusion coefficient set | stated; not checked |
| secondary ionization (Dalgarno, Yan and Liu 1999 eq. 14) | stated for `x <= 0.1` | extrapolated to `x = 1`; He/H fixed at 0.1; no metals; 30 eV floor |
| static-slab escape probability (Neufeld 1990 eq. 3.27, Harrington 1973 eq. 40) | static slab damping wings, `(a tau)^(1/3) > 10` | boost is a tuned parameter (D0 C7) |
| bremsstrahlung Gaunt factor | table range | clamped to the table (`gbar_ff`, `Cool_coeff.f90`) |

**Recorded today: nothing.** A grep of the tree for `out_of_domain`,
`validity_state` or `domain_flag` returns no match. Every clamp above is
silent: the returned value carries no signal that it came from an edge. A2
needs at minimum an activation count for each closure, and preferably the worst
excursion, on the pattern the diffusion module already uses for its range clip
(`he_fraction_over_one` / `he_fraction_under_zero`, declared in the head of
`binary_element_diffusion.f90`, which exist precisely because "the clipped `X`
cannot tell an overshoot from an exact 1"). Two of the rows above have gained
that record since: the H3+ domain counters (item B3b-H3, 2026-09-06) and the CO
destruction domain ledgers (item B3b-CO, 2026-09-06).

### 4.3 Rejected numerical trial with no adopted contribution (does not invalidate)

| producer | recorded today |
|---|---|
| steady-solver line-search trials and Jacobian probes | the sweep state kind `ieq_state_steady_candidate` (declared in the head of `ionization_equilibrium.f90`) makes a probe neither raise nor clear the non-root streak (`nonroot_streak_update`), and the acceptance ledger prints one section for each state kind (`write_ioniz_eq_acceptance_report`, called at the end of `EXHALE_main.f90`) |
| RK stage positivity failures retaken at half `dt` | `n_steps_dt_halved`, `n_dt_halvings`, raised in the attempted step and printed by `write_run_counter_report` (`EXHALE_main.f90`) |
| first-order flux repair | three separate numbers, and the separation is exactly the one B6 needs: `calls` (stages offered the repair), `faces` (interfaces repaired), and `faces in accepted steps` - the repairs of an attempt discarded by the `dt` bisection are counted in the second and **not** the third (`EXHALE_main.f90`, `write_run_counter_report`) |
| carrier trials above the headroom | refused, not clamped (the element-headroom note in `carrier_steady_residual`, `diffusive_photochemistry.f90`; the refusal itself in `eval_residual`, `steady_newton.f90`) |

This state is the best served today. The flux-repair counter is the model to
follow: an attempted-step count and an accepted-step count as separate numbers.

### 4.4 Unbudgeted accepted correction (invalidates)

A correction that alters a conserved quantity or a physical field in the
adopted state with no source term behind it.

| correction | producer | recorded today |
|---|---|---|
| **composition projection** (removed 2026-09-06, item B3c) | the temperature-preserving pressure reset that stood in the marching source stage of `EXHALE_main.f90` (D0 C1) | **nothing**. `EXHALE_UPDATE_MAP` reports a `chem` column when armed (`update_map_end_step`), but it is a diagnostic mode, not a run record |
| **`rho` rewritten by the chemistry** (removed 2026-09-06, item B3c) | `ionization_equilibrium.f90`, `ioniz_eq`, which now takes the density `intent(in)` (D0 C2) | as above |
| **simplex clamp of the ionization fractions** | `clamp_fractions_to_element_budget` (`ionization_equilibrium.f90`), applied at or below 1e-6 | counted in the acceptance ledger as class 3 (a projected/handback state that rechecks as a root) |
| **class-4 acceptance** (a non-root held by the relaxation amnesty) | `nonroot_streak_update` (`ionization_equilibrium.f90`): reported loudly, counted, and allowed for at most `ieq_nonroot_streak_stop` consecutive sweeps, after which the run stops | `n_acc(4)`, `acc_resmax(4)`, longest non-root streak, in the acceptance ledger |
| **CO destruction domain** (the ceiling's successor) | `carrier_co_domain_take` (`diffusive_photochemistry.f90`): the destruction-only CO row is legitimate only where `tau_dest / tau_res` is at or below `co_domain_f_dom = 0.1`, and cells outside that, cells above the shielding table's excitation-temperature limit and cells in which He+ + CO leads the He+ loss are counted | **cumulative, by ledger family**, read by `carrier_co_domain_record`; entered into the class 4.2 out-of-domain closure count of the certification (`certification.f90`, `read_validity_states`) |
| **carrier element rescale** | `limit_to_element_budget` | `pct_limited`, `pct_worst_limit`, `pct_cell_constrained` by cell (last call only, `carrier_transport_diagnostics` and `carrier_cell_is_constrained`) |
| **ghost extrapolation dropped to zero gradient** | `free_outflow_ghost` (`Apply_BC.f90`) | `n_ghost_cells_positivity_limited`, printed by `write_run_counter_report` (`EXHALE_main.f90`) |
| **reconstruction face states scaled toward the cell average** | positivity limiter | `n_faces_positivity_limited`, printed in the same report |
| **helium fraction range clip** | `binary_element_diffusion` | `he_fraction_over_one`, `he_fraction_under_zero`, `trace_ratio_under_zero` (declared in the module head); printed only under `EXHALE_DIFFUSION_CHECK=1` |

Two of these, the composition projection and the `rho` rewrite, were the
largest unbudgeted corrections in the code and were the only ones with **no run
record at all**. Everything else in the table has at least a counter. Both were
removed on 2026-09-06 (item B3c), so the category now reads zero, as the
closing paragraph of this section says.

**The two temperature floors have left this table (2026-09-06, B2 and B2b).**
Neither stage clamps any more: the energy update (`energy_semi_implicit.f90`)
and the conduction stage (`viscous_conduction.f90`) are both bracketed and
residual-controlled, and reaching the lower bracket end is a FAILURE that
stops the run with a status of its own rather than a state the run adopts. A
rejected attempt contributes nothing to the state a certification judges, so
the two counters belong to section 4.3 and neither invalidates:

| attempt | producer | recorded today |
|---|---|---|
| **energy-update cells at the lower bracket end** | `energy_semi_implicit.f90`, status `FLOOR` | `n_energy_floor_hits`, `n_energy_floor_cells()`; read by `certification.f90` into `n_energy_floor_attempts` and printed in the report's attempts row, beside the validity states and never decisive |
| **conduction cells at the lower bracket end** | `viscous_conduction.f90`, status `CONDUCTION_FLOOR` | `n_conduction_floor_hits`, `n_conduction_floor_cells()`; read into `n_conduction_floor_attempts`, same row |

The unbudgeted-correction category itself stays in the report and in the
certification verdict, and reads zero today: it is where a producer of a
genuine unbudgeted correction (the composition projection and the `rho`
rewrite above, once B3c gives them a record) reports.

### 4.5 Specified external reservoir (informational)

| reservoir | producer | recorded today |
|---|---|---|
| base `(p, s)` reservoir | `set_base_reservoir(ntot_bc + dp_bc, 1, ...)` called from `input_read` (`input_read.f90`), continued to the face by `continue_hydrostatic_isentrope` (`base_boundary.f90`) | the resolved-config report; the base level source string `base_level_source` |
| base composition handoff | `base.inp` or the profile; H2 Dirichlet inflow at the base face | `write_resolved_config` (`write_setup_report.f90`); the profile's `solution_id` pairing |
| base infrared field | `Base IR field`: the coolants see a diluted `B_nu(T0)` from below | the key echo; `write_output` carries the flag |
| diffusion inner Dirichlet | `element_diffusion_step` (`binary_element_diffusion.f90`) | nothing |
| base inflow admitted by the characteristic BC | `check_base_inflow_is_subsonic` | `base_inflow_mach_max`, and the supersonic case with a step count and first step (`EXHALE_main.f90`, `write_run_counter_report`); **always reported** for a run with any base inflow, which is the right pattern for an informational state |

---

## 5. Prerequisite generator

Rules, from `PLAN_20260906_rev2.md` section 3 (Step D) with F2's correction:

- **always**: A0 and A0-impl in physical mode, A1, A2, the exact output state
  from B3c's assembly, B6 clear or the run excluded, C1 and C3 for the input
  field;
- **B2** wherever the semi-implicit thermal path runs. That is every
  configuration below: `use_semi_implicit_energy` defaults `.true.`
  (`parameters.f90`) and none of these cases sets `Energy solver:`;
- **B3c** wherever the composition reset of the marching source stage
  (`EXHALE_main.f90`) runs, which was unconditional, so again every
  configuration below. B3c landed on 2026-09-06 and removed that reset;
- **B3a** wherever recoverable rejection is claimed;
- **A3** with active carrier transport;
- **B3b** subitems only where the process is active;
- **B4** with species or element transport active;
- **B5** never mandatory: it is not required where a time-converged solution
  passes A2 certification.

`n/a` below means "recorded as not applicable", not "skipped".

| configuration | source | active equation set beyond hydro | A3 | B3b subitems | B4 | B5 |
|---|---|---|---|---|---|---|
| pure H/He, power law | `examples/tutorial_nometals` (no `metals.inp`, `Include He23S? True`, no `Stellar Teff` so H(n=2) off) | 4 local rows: H+, He+, He++, He 2^3S | n/a | photoevent partition only; H2, oxygen, metal, H3+ items **n/a** | n/a | not required if the marching solution certifies |
| He 2^3S + metals | `backup/regression/wasp_full` (and `wasp_he23off`, same metals, triplet off) | 4 + 2x7 local rows (C, N, O, Mg, Ca, Na, Fe stated); H(n=2) active (`Stellar Teff` and `Stellar radius` both given) | n/a | photoevent partition; metal reaction energies; H(n=2)/Lyman-alpha ownership. H2, oxygen, H3+ **n/a** | n/a | as above |
| hot Uranus molecular | `backup/regression/mol_base_handoff` (`Molecular chemistry: True`, `base.inp`, `Solver: Newton`, fixed-step snapshot with `maxsteps`) | 8 local rows (H+, He+, He++, H2, H2+, H3+, HeH+, He 2^3S); no transport | n/a | H2 photoevent channels; the common H2 state model; H3+ emission model and domain; molecular reaction heats. Oxygen **n/a** | n/a | it uses `Solver: Newton`, so B5's interfaces apply, but as a fixed-step snapshot it is a relaxation state under A0 and not certifiable |
| carriers | `backup/regression/mol_carrier` (adds `Molecular carrier transport: True`) | the above **plus** one transported H2 balance, whose row in the cell solve becomes the identity `x4 - x_h2_fix` | **required** | as above | **required** (species transport) | as above |
| oxygen | `backup/regression/oxygen_chemistry` (`Oxygen chemistry: True`, `Stellar LW/FUV/Lya` fluxes, metals C/N/O; carrier transport **on by default**) | four transported balances (H2, OH, H2O, CO), with CO carrying no kinetic row; oxygen rows OH and H2O in the cell solve | **required** | oxygen reaction energies **or exclusion** (4.1); CO model or exclusion; FUV photolysis; plus the molecular items | **required** | as above |
| diffusion | `backup/regression/mol_diffusion` (`He_diffusion: True`, `He_Kzz: 1.0e9`, molecular, `base.inp`) | one elemental transport balance for `X = rho_He/rho`, on top of the 8 local rows | n/a (no carrier transport stated, oxygen off, so `carrier_transport` stays false) | molecular items | **required** (element transport) | as above |
| profile | `backup/regression/lower_profile` (`Lower atmosphere profile:`, `He_diffusion: True`, `He_metal_diffusion: True`, metals from the profile, no molecules) | He/H elemental balance **and** the trace-element balances, which have no counter-flux (D0 C29); atomic + metal local rows; `K_zz(p)` from the profile | n/a | metal items | **required**, and it owns C29 | as above |
| H+ transport | `backup/regression/hp_front` (`Molecular chemistry`, `Molecular carrier transport`, `Ionization transport: True`, `Load IC? True`, `base.inp`) | two transported balances (H2 and H+), both replacing their cell rows by identities; the two lower ghosts still re-solve H/H+ locally (D0 C18) | **required** | molecular items | **required** | **impossible today**: this configuration refuses `Solver: Newton`, PTC and the coupled carrier solve (`input_read.f90`, `input_read`, the `Ionization transport` refusal block), and the steady packing has no proton unknown (D0 C19), so the only route to acceptance is a marching solution certified by A2 |

Every row above additionally carries: A0, A0-impl, A1, A2, B2, B3c, B6, C1,
C3. None of the eight is eligible for D-physical today, for the reason F2
gives: the temperature-preserving composition reset is part of the source
update in all of them.

---

## 6. Open questions

Only where the code admits two readings.

1. **Is `bg_ready` a property of the state or of the module?** It is owned by
   `ionization_equilibrium` and imported by the carrier module
   (the `use` statement in the head of `diffusive_photochemistry.f90`). A2 must decide whether "the
   backgrounds of the state being certified" means (a) `bg_ready` plus an
   assertion that the last sweep was on this state, or (b) a background
   snapshot passed in with the state. The code supports neither today.
2. **Which context refreshes the carrier headroom?** The refresh is gated on
   `weno_mode .eq. 1` (`diffusive_photochemistry.f90`, the element-headroom
   note in `carrier_steady_residual`), that is on being
   the outer iterate rather than a probe or a line-search trial. A2's evaluation
   is a fourth context that the gate does not name, and reading it as "the
   iterate" or as "a probe" gives different frozen constraints.
3. **Is the CO ceiling an unbudgeted accepted correction or active unvalidated
   physics? Closed, and not by an answer.** The ceiling was deleted on
   2026-09-06 (item CEILING-DEL): CO is destroyed by published rates in its
   own row, so nothing removes carbon and oxygen without a rate and nothing
   moves chemical energy without a source. What is reported in its place is
   the destruction model's domain, under 4.2.
4. **Does a class-3 acceptance (a clamped root) count as an unbudgeted
   correction?** The clamp is bounded by 1e-6 and is argued as round-off of a
   root on a face of the allowed region (`ionization_equilibrium.f90`, the note
   above the `clamp_fractions_to_element_budget` gate in `ioniz_eq`),
   which reads as a numerical repair; but it does move the composition of the
   adopted state.
5. **Whose equation is the H(n=2) system?** It is solved outside every solver
   and fed back lagged (2.5). A2 can treat it as a closure with an accuracy
   condition on the lag, or as an active equation that is simply never
   iterated. The code states the lag as a design choice
   (`excited_hydrogen.f90`, module header) and does not measure it.
