# Physics probe tests

Maintained acceptance tests derived from the read-only diagnostic programs in
`docs/audit_20260905/`. The originals stay where they are and are not moved,
edited or converted mechanically: they are the record of what each review
measured at `35d9dd5`. What is here instead asserts, against a stated
reference and a stated tolerance, the physics those probes exposed, and exits
nonzero when an assertion fails.

Every driver links the **production** sources under `src/modules`, never a copy
of them. `source_closure.py` follows the `use` statements of each driver
through `src/modules` (and the generated `build_stamp.f90`, which
`utilities.f90` reads) and prints the closure in an order in which each file
compiles after the modules it uses; `run.sh` compiles exactly that list.

## Running

```bash
src/tests/physics_probe/run.sh
```

Objects, module files and executables go to `build/tests/physics_probe/`,
which the script creates. It needs nothing else built first.

`EXHALE_TEST_OUT` sends them somewhere else instead, which is what lets two
builds of this suite in one tree run at the same time; the drivers read the
same directory as their scratch, under the name `EXHALE_TEST_WORK`, and
setting that name selects it too:

```bash
EXHALE_TEST_OUT=/tmp/pp src/tests/physics_probe/run.sh
# against a private build:
make OBJDIR=build_x EXE=EXHALE_x.x
EXHALE_OBJDIR=$PWD/build_x EXHALE_TEST_OUT=/tmp/pp src/tests/physics_probe/run.sh
```

`EXHALE_OBJDIR` is the build the suite reads its `build_stamp.f90` from, the
only file it takes from a build tree, since the drivers compile the
production sources themselves; when the named build has not generated one,
`run.sh` writes a stamp saying so and compiles that. Nothing tested here
depends on its contents. No test here runs `EXHALE.x`, so `EXHALE_EXE` has no
effect on this suite.

`FC` selects the compiler (default `gfortran`); the flags are GNU syntax
(`-O0 -g -fcheck=all -fbacktrace -fopenmp`).

Two further keys are read inside the suite, both optional and both without
effect when unset:

- `EXHALE_PROBE_STATE` names the directory of a written state, from which
  group 37 (`molecular_energy_recipients.f90`) reads the three cells its
  uncertainty record is quoted at. The pair read is
  `{Hydro_ioniz_IC.txt, Ion_species_IC.txt}`, the state a run hands back;
  the pair without `_IC` is the product written FROM a measurement of that
  state and is a different state. Unset, nothing is read and the driver's
  compiled cells stand, so the suite runs with no state in the tree.
- `EXHALE_REACTION_HEAT_RECIPIENTS=0` is a production key, read by
  `molecular_reaction_heat.f90` inside the closure of that same driver. It
  puts the four recipient corrections of group 37 back where they were;
  running the driver with it set is the control, and the rows that state the
  correction then fail by construction.

A full build and run takes about thirteen seconds (MEASURED 13.4 s, 34
Fortran drivers and 7 Python scripts, 1604 assertions, gfortran, one thread,
2026-09-18).

Output is one line per assertion,

```
PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
```

and the exit status is nonzero if any assertion failed or any build failed.
Lines beginning `note:`, `DIAGNOSTIC`, `MEASURED` or `#` are diagnostics, not
assertions: they carry a number, a table or a record that has no comparison
attached. A `SKIP` line reports a check whose external input is absent, with
the path that was looked for; it is not a failure (only
`metal_photoion_table_transcription.py` emits one).

## Files

| File | Contents |
|---|---|
| `run.sh` | Builds into `build/tests/physics_probe/`, runs every driver, exits nonzero on any failure. |
| `source_closure.py` | The dependency-ordered list of production sources a driver needs. |
| `assertion_report.f90` | The verdict line, the failure count, and the four comparison kinds (relative, absolute, within a factor, strictly positive). Every comparison is written so that a NaN fails. |
| `riemann_wave_speeds.f90` | Groups 1 and 2 below. |
| `h2_channel_detector_ratio.f90` | Group 3. |
| `h3p_cooling_limits.f90` | Group 4. |
| `h2_rovibrational_identity.f90` | Group 5. |
| `photoevent_energy_ledger.f90` | Group 6. |
| `radiation_exchange_ledger.py` | Group 7. |
| `wind_ae_bridge_composition.f90` | Group 8. |
| `opacity_model_p_ledger.f90` | Group 9. |
| `lya_band_transmission.f90` | Group 10. |
| `helium_level_cooling_ledger.f90` | Group 12. |
| `lyman_werner_cell_mean.f90` | Group 13. |
| `voronov_carbon_ionization.f90` | Group 14. |
| `atomic_mass_and_radius_constants.f90` | Group 15. |
| `jupiter_radius_python_tools.py` | Group 16. |
| `ionization_threshold_turn_on.f90`, `threshold_literal_uniqueness.py`, `threshold_literal_exceptions.txt` | Group 17. |
| `h2_infrared_line_populations.f90`, `h2_partition_sum_uniqueness.py` | Group 18. |
| `lya_escape_probability.f90` | Group 19. |
| `roe_low_mach_velocity_jump.f90` | Group 20. |
| `weno3_reconstruction_order.f90` | Group 21. |
| `positivity_limiter_scaling.f90` | Group 22. |
| `lya_beam_cell_mean.f90` | Group 23. |
| `fine_structure_escape_probability.f90` | Group 24. |
| `xuv_cell_mean_attenuation.f90` | Group 25. |
| `photoionization_field_substitution.f90` | Group 26. |
| `photoionization_field_self_consistency.f90` | Group 27. |
| `absorbed_fraction_switch.f90` | Group 28. |
| `species_formation_energy_table.f90`, `species_formation_energy_python_table.py` | Group 29. |
| `heating_channel_closure.f90`, `heating_sum_uniqueness.py` | Group 30. |
| `co_destruction.f90` | Group 31. |
| `co_helium_ion_sink.f90` | Group 32. |
| `h3p_recombination_branching.f90` | Group 33. |
| `water_photolysis_lyman_alpha_yields.f90` | Group 34. |
| `metal_photoionization_fits.f90`, `metal_photoion_table_transcription.py` | Group 35. |
| `recombination_coefficient_fits.f90` | Group 36. |
| `molecular_energy_recipients.f90` | Group 37. |
| `molecular_third_body_in_the_chemistry.f90` | Group 38. |

The group numbers run 1 to 38 with no 11: that number was never used. Groups
20 to 38 are in the order `run.sh` builds and runs them; the earlier numbers
are in the order they were written, which is not the same order.

## What each test asserts

### 1. Lax-Friedrichs signal speed (`riemann_wave_speeds.f90`)

Origin: the LLF block of `docs/audit_20260905/audit_probe.f90`.
Production routine: `lax_friedrichs_flux` (`src/modules/flux/Num_Fluxes.f90`).

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `llf_signal_speed` | Viscosity coefficient inferred from the returned mass flux for WL=(1,-2,0.6), WR=(2,-2,1.2), gamma=5/3 | `max(abs(v_L)+a_L, abs(v_R)+a_R)` = 3, the spectral radius bound the flux is monotone under | 1e-12 relative |
| `llf_mass_flux` | Mass flux of the same pair | -4.5, the value that bound gives | 1e-12 relative |
| `llf_mirror_mass_flux` | Mass flux of the pair mirrored by `(rho,v,p) -> (rho,-v,p)` with the sides exchanged | minus the unmirrored mass flux, an exact symmetry of the Euler equations | 1e-12 relative |

**RED at 35d9dd5.** `Num_Fluxes.f90` line 287 sets the coefficient to
`max(abs(v_L+a_L), abs(v_R+a_R))`, which cancels rather than adds when the
velocity is negative: it returns 1 instead of 3, the mass flux -3.5 instead of
-4.5, and the mirrored pair 4.5 instead of 3.5.

**GREEN on all three** with `a1 = max(abs(v_L)+a_L, abs(v_R)+a_R)`
(Phase 1 item A1), the spectral-radius bound the Perthame and Shu positivity
lemma is stated for.

### 2. ROE star state (`riemann_wave_speeds.f90`)

Origin: the ROE block of `audit_probe.f90`, `roe_equal_pressure.f90` and
`roe_flux_vacuum_probe.f90`. Production routines: `speed_estimate_ROE`
(`src/modules/flux/speed_estimate_ROE.f90`) and the ROE branch of `Num_flux`,
reached by setting `flux` of `global_parameters` to `'ROE'`, which is what the
input key `Numerical flux:` writes. The interface the estimator answers to,
including the three status values, is `docs/a2_roe_interface.md`.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `roe_rarefaction_scale_invariance_left`, `..._right` | Star sound speeds of WL=(1,-1,1), WR=(1,1,0.1) with rho and p both scaled by 10 | the unscaled values: the Euler equations are invariant under `(rho,p) -> (s rho, s p)` at fixed velocity and fixed p/rho | 1e-12 relative |
| `roe_shock_scale_invariance_left`, `..._right` | The same for WL=(2,5,2), WR=(8,-5,0.2), which selects the two-shock branch | the unscaled values | 1e-12 relative |
| `roe_expansion_status` | Status of the equal-pressure expansion WL=(1,-1,1), WR=(1,1,1) | `ROE_STAR_OK`: the two gas edges still meet, so a star state exists | exact |
| `roe_expansion_star_pressure` | Its star pressure | strictly positive | strict |
| `roe_expansion_star_sound_speed_positive_left`, `..._right` | Its star sound speeds | strictly positive and finite | strict |
| `roe_expansion_star_sound_speed_left`, `..._right` | The same values | `sqrt(gamma) - (gamma-1)/2` = 0.957661..., the exact two-rarefaction star sound speed | 5e-2 relative, the accuracy of an approximate estimator |
| `roe_contact_status`, `roe_contact_star_pressure`, `roe_contact_star_velocity`, `roe_contact_star_sound_speed_left`, `..._right` | Star state of the stationary contact WL=(1,0,1), WR=(0.125,0,1) | `ROE_STAR_OK`, p_* = 1, u_* = 0, and each star sound speed equal to its own data value: with no velocity jump and no pressure jump the linearized branch is exact, so the data are the exact solution | 1e-12, not the estimator's 5e-2 |
| `roe_shock_status`, `roe_shock_star_pressure`, `roe_shock_star_velocity`, `roe_shock_star_sound_speed_left`, `..._right` | Star state of the asymmetric two-shock collision WL=(1,0.5,1), WR=(3,-0.5,1.5) | the exact Riemann solution p_* = 2.3768663335, u_* = -0.2357045218, c_*L = 1.5505344803, c_*R = 1.0024512384 | 5e-2 relative |
| `roe_rarefaction_status`, `roe_rarefaction_star_pressure`, `roe_rarefaction_star_velocity`, `roe_rarefaction_star_sound_speed_left`, `..._right` | Star state of the asymmetric two-rarefaction expansion WL=(1,-2,1), WR=(0.5,1,0.2) | the exact Riemann solution p_* = 0.0201759337, u_* = 0.0987433798, c_*L = 0.5914133221, c_*R = 0.5160777075 | 5e-2 relative |
| `roe_vacuum_status` | Status of the vacuum-producing expansion WL=(1,-10,1), WR=(1,10,1) | `ROE_STAR_VACUUM`: the two rarefaction fans separate and no star state exists | exact |
| `roe_vacuum_edge_speed_left`, `..._right` | Its gas-edge speeds `u -/+ 2c/(gamma-1)` | -6.127017 and +6.127017 | 1e-6 relative |
| `roe_vacuum_update_density`, `roe_vacuum_update_thermal_energy` | Density and thermal energy of the left cell after one conservative Euler update at `dt/dx = 0.01` with the production ROE flux of WL=(1,-10,1), WR=(1,10,1) | strictly positive; the exact solution is a vacuum with zero interface flux, and no update may produce a negative density or a negative thermal energy | strict |

**RED at 35d9dd5** on nine of the ten assertions the suite then carried, the
exception being `roe_vacuum_update_density`. The two-rarefaction branch took a
square root where the exact estimate has the exponent `2 gamma/(gamma-1)`, and
the two-shock branch dropped the factor `rho` from the star densities, so
neither branch was invariant under the scaling: the rarefaction sound speeds
moved by a factor 0.66 and the shock ones by 3.2 when rho and p were
multiplied by 10. The equal-pressure expansion has `p_max/p_min = 1`, so the
PVRS guess was kept and it is negative (-0.291): both star sound speeds came
back NaN. The vacuum-producing expansion left the density positive (0.9) but
the thermal energy at -2.93.

**GREEN** on all 30 of them with the estimator of `docs/a2_roe_interface.md`
(Phase 1 item A2): the two branch formulas are the published ones, the vacuum
case is a status and not a square root of a negative number, and a face
without an admissible star state takes the HLLE flux, which leaves the left
cell at density 0.9 and thermal energy +2.53.

### 3. H2 channel split (`h2_channel_detector_ratio.f90`)

Origin: `docs/audit_20260905/h2_detector_ratio.f90`. Production routines:
`h2_channel_cross_sections` and `h2_channel_stoichiometry`
(`src/modules/functions/h2_photo_channels.f90`) and
`frac_H2_dissociative_ionization` (`src/modules/functions/cross_sec.f90`), at
80 eV with `h2_double_ionization_model = 'chung80'` and unit total cross
section.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `h2_detector_event_ratio` | `(S+D)/M` of the returned channels | `r = f/(1-f)` with `f = frac_H2_dissociative_ionization(80)`, the tabulated signal ratio of Chung et al. (1993) Table II, whose detector "responds with a single pulse" (p. 886) for the unresolved double event | 1e-12 relative |
| `h2_double_event_fraction` | `D/(S+D)` | 0.20, Chung's 20 per cent read as an event fraction of the same signal | 1e-12 relative |
| `h2_channel_sum` | Sum of the four channels | 1, the total absorption cross section passed in | 1e-12 relative |
| `h2_channel_nuclei_closure_1..4` | `2 dH2 + 2 dH2p + dH + dHp` of each channel | 0 | 1e-12 absolute |
| `h2_channel_charge_closure_1..4` | `dH2p + dHp - dele` of each channel | 0 | 1e-12 absolute |

Both tolerances on the first two are round-off, not the 4-5 per cent Chung
quotes for the data: what is tested is that the channels reproduce the numbers
the code itself holds, not the accuracy of the measurement.

**RED at 35d9dd5** on the two ratio assertions, **GREEN** on the sum and the
eight closures. The module reads the same two numbers as a proton-number
ratio, `(S+2D)/M = r` and `2D/(S+2D) = 0.20`, so it returns
`(S+D)/M = 0.2545` against `r = 0.2827` and `D/(S+D) = 0.1111` against 0.20.
The reading of the published text is the second review's
(`docs/development_plan_20260905_review2.md` Q1) against the published
version in `references/`.

### 4. H3+ cooling limits (`h3p_cooling_limits.f90`)

Origin: the H3+ lines of `audit_probe.f90` and of `roe_equal_pressure.f90`.
Production routines: `h3p_cooling_rate`, `h3p_emission_lte`,
`h3p_nonlte_factor` (`src/modules/lower_atmosphere/h3p_cooling.f90`).

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `h3p_cooling_vanishes_without_h2` | `Lambda(1000 K, n_H3+=1, n_H2=0) / Lambda(..., n_H2=1e2)` | 0: the routine models excitation by H2 collisions only, so with no collision partner there is no emission | 1e-12 absolute |
| `h3p_cooling_collisional_slope` | `Lambda(n_H2=1e2)/Lambda(n_H2=1)` at 1000 K | 100: two decades of a rate that is linear in the collision-partner density far below thermalization (Miller et al. 2013 Table 6 puts that near 1e10 cm^-3 at 1000 K; 1 and 100 cm^-3 are four to six decades below the lowest tabulated density) | factor 2 |
| `h3p_lte_emission_join_300`, `_800`, `_1800` | Relative jump of `h3p_emission_lte` across each fit join | 0: one physical function of temperature fitted piecewise | 1e-3 absolute (the fits are quoted to four figures) |
| `h3p_nonlte_factor_high_density_join` | Relative jump of `h3p_nonlte_factor(5000 K, n_H2)` at the last tabulated density 1e14 cm^-3 | 0: the factor reaches unity from below, continuously | 1e-3 absolute |

**RED at 35d9dd5**, all six. The density argument is clamped at the lowest
tabulated density, so the rate is the same constant at n_H2 = 0, 1 and 100
(ratio 1 where 100 is required, and nonzero without any H2 at all); the LTE
emission jumps by 43 per cent at 300 K, -0.11 per cent at 800 K and -2.35 per
cent at 1800 K; and the non-LTE factor steps from 0.9985 to 1 at 1e14 cm^-3,
a jump of 0.15 per cent. The 800 K join fails marginally, by a factor 1.13
over its tolerance.

### 5. H2 rovibrational energy identity (`h2_rovibrational_identity.f90`)

Origin: the H2 block of `audit_probe.f90`. Production routines:
`h2_partition_function` and `h2_rovibrational_sum`
(`src/modules/states/caloric_eos.f90`), the one Boltzmann sum over the one H2
level ladder; `h2_rovibrational_energy_and_heat_capacity` (same module), the
tabulated internal energy the energy equation carries; `keq_H_H_to_H2`
(`src/modules/lower_atmosphere/mol_rates.f90`), the equilibrium constant built
from that same partition function; and `equilibrium_constant_conc`
(`src/modules/lower_atmosphere/oxygen_rates.f90`), the same equilibrium
constant from the NIST-JANAF Shomate data, an independent source.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `h2_rovibrational_energy_at_300`, `_1000`, `_2000`, `_4000`, `_8000` | `T d ln Q/d ln T` of the chemistry partition function, by a centered difference of relative step 1e-5 | `u_rv/k` returned by the caloric EOS table at the same temperature: for one level set, `u/k = T^2 d ln Q/dT` exactly | 1e-3 relative |
| `h2_partition_function_identity_at_300` … `_8000` | the same logarithmic derivative | `u_rv/k` from the exact level sum rather than its table | 1e-8 relative |
| `h2_equilibrium_constant_table_at_1000`, `_2000`, `_4000` | `keq_H_H_to_H2`, the tabulated production path | the closed form `Lam(2m) Q exp(D0 hc/kT)/(4 Lam(m))^2` from the same `Q` | 5e-4 relative, the tabulation error `mol_rates` states |
| `h2_equilibrium_constant_janaf_at_1000`, `_2000`, `_4000` | the same `keq_H_H_to_H2` | `n(H2)/n(H)^2` from the NIST-JANAF Shomate enthalpies and entropies of H and H2 | 1e-3 relative |
| `h2_partition_function_movement_at_1000`, `_2000`, `_4000` | `Q(observed ladder)/Q(Huber-Herzberg model ladder)` | 1.017599, 1.036501, 1.077285, the movement B1 decision 9 adopts | 1e-4 relative |

The finite difference has a truncation error of order 1e-10, so the 1e-3
tolerance of the first group is a statement about the level sets and not about
the differencing.

**RED at 35d9dd5**, the first five: the chemistry energy was below the EOS
energy by 0.63 per cent at 300 K, 1.8 at 1000 K, 3.1 at 2000 K, 5.1 at 4000 K
and 8.3 at 8000 K. The two routines described different level sets, so the
equilibrium constant and the energy equation disagreed about how much energy
an H2 molecule holds.

**GREEN 2026-09-06.** The chemistry now reads the same partition function the
equation of state is a moment of (B1 T6.4, decision 9: the observed Roueff et
al. 2019 ladder). The retired model ladder placed the levels too high, not
too few: of the 8000 K gap only 3.1 per cent came from the 55 levels it did
not reach. The independent NIST-JANAF route, which deviated by 1.0 per cent at
600 K rising to 10.1 per cent at 5000 K against the model ladder, now agrees
to better than 0.1 per cent.

### 6. Photoevent energy ledger (`photoevent_energy_ledger.f90`)

Origin: `docs/audit_20260905/photoevent_energy_check.py`. That script
evaluates a **superseded** equation as a counterexample and is therefore not
an assertion about the code; it is not ported. What is asserted here is the
accepted equation of `docs/development_plan_20260905_rev3.md` section 4.2
item 3 and section 4.3 (A5): for one event in an isolated cell,

```
Delta u_th + Delta E_chem = absorbed photon energy.
```

Production quantities: `photoelectron_share`
(`src/modules/radiation/util_ion_eq.f90`), the heat share the H, He and H2
heating integrals all call; `e_th_HI` and `e_th_H2_dd`
(`src/modules/init/parameters.f90`), the thresholds those calls pass it;
`h2_dissociation_energy_eV` (`src/modules/lower_atmosphere/mol_rates.f90`),
the single definition of D0(H2). `photoelectron_share` **is** callable in
isolation: its module's closure is 20 production modules plus the generated
build stamp, needs no LAPACK and compiles in about 1.5 s, so the test calls
the production function and not a copy of its formula.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `h_photoionization_event_energy_sum` | `share(e_th_HI, 20 eV)*20 + e_th_HI` | 20 eV | 1e-12 relative |
| `h2_double_chemical_product_energy` | `D0(H2) + 2 e_th_HI` | 31.675 eV, the value the plan quotes for `Delta E_chem` of the double channel | 1e-3 relative |
| `h2_double_photoevent_energy_sum` | `share(e_th_H2_dd, 80 eV)*80 + D0(H2) + 2 e_th_HI` | 80 eV | 1e-12 relative |

**GREEN** on the first two, **RED** on the third at 35d9dd5. H
photoionization closes exactly (6.4 eV of heat plus 13.6 eV stored). The H2
double channel deposits `E - 51.4 = 28.6` eV as heat and stores 31.678 eV,
which accounts for 60.278 of the 80 eV: **19.722 eV has no recipient**
anywhere in the accounting. The driver prints that number on a `note:` line.

The 0.003 eV by which the production `D0 + 2 e_th_HI` exceeds the plan's
31.675 eV is `e_th_HI = 13.6` eV in `parameters.f90` against the
13.598434599 eV the audit probe used; it is four orders of magnitude below
the gap.

### 7. Radiation exchange ledger (`radiation_exchange_ledger.py`)

Origin: `docs/audit_20260905/review3_ledger_checks.py`. Pure arithmetic over
the convention itself; it holds no production symbol and is **GREEN by
construction**. It is kept so that the ownership rule the Fortran tests are
written against is pinned in one readable place.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `two_cell_radiation_exchange_material_sum` | Material energy change of a pair of cells when one emits a 10 eV photon that the other fully absorbs | 0: "escaping" means leaving the local material system of the emitting cell, not leaving the atmosphere (rev3 section 4.2 item 3) | 1e-12 absolute |
| `h_photoionization_event_energy_sum_chi_13.600000000`, `..._13.598434599` | `(E - chi) + chi` at E = 20 eV | 20 eV, for either value of the ionization potential | 1e-12 absolute |

The driver also prints, on a `note:` line, the 10 eV of material energy the
rejected rule (subtract only domain-escaping emission) would create out of
nothing.

### 8. Wind-AE bridge composition (`wind_ae_bridge_composition.f90`)

Production routine: `wae_read_exhale_input`
(`src/modules/wind_ae/wae_exhale_input.f90`), the front end that maps an
EXHALE `input.inp` onto the Wind-AE parameter list, reached by `IC mode:
windae` and by the standalone `wind_ae_ic.x`.

Wind-AE never stores a number ratio: it forms every number density from the
mass fractions and the atomic masses it was given, `n_j = rho
HX(j)/atomic_mass(j)` (`wae_soe.f90` 138, 166, 216; `wae_glq_rates.f90` 67).
The helium abundance the solver works at is therefore
`(HX(2)/m_He)/(HX(1)/m_H)`, and it must equal the `He/H number ratio` of the
input file: a run asks for a helium abundance and the solver has to relax at
that one. The driver writes a minimal `input.inp`, reads it back through the
production routine and reconstructs the ratio, for He/H = 0, 0.0793, 0.1 and
1.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `wind_ae_helium_hydrogen_number_ratio_at_<He/H>` | `(HX(2)/m_He)/(HX(1)/m_H)` after the mapping | the `He/H number ratio` of the file | 1e-12 relative |
| `wind_ae_mass_fraction_sum_at_<He/H>` | `HX(1) + HX(2)` | 1: hydrogen and helium are the only species Wind-AE carries | 1e-12 relative |

Both tolerances are round-off: these are exact identities of a unit
conversion, not fits.

**RED at 35d9dd5** on the three nonzero abundances, GREEN on the mass-fraction
sums and on He/H = 0. The mass fractions were built with helium weighing
exactly 4 hydrogen masses while the solver was handed the CODATA masses, so
the reconstructed ratio was 1.0071696 times the requested one: 0.0798685
against 0.0793, 0.1007170 against 0.1, 1.0071696 against 1.

### 9. Opacity model 'P' radiative ledger (`opacity_model_p_ledger.f90`)

Production routine: `PH_heat_HHe` (`src/modules/radiation/util_ion_eq.f90`),
called on a twelve-cell uniform slab of H I, He I, He II, He 2^3S, C I and
Mg II under the power-law spectrum of `backup/regression/wasp_full`, with
`opa_pf = 10` in cell 4 and 1 elsewhere. `opa_pf` is the cross-section factor
f(p) that opacity model 'P' leaves in place (`opacity_pT_factor`); the driver
sets it directly, so the assertions are about the ledger of that state and not
about the pressure parameterization producing it.

The reference is the beam itself, on the same energy grid, with the same
production cross sections and the same columns: dtau(E) = sigma(E) n dr
opa_pf, and the flux F(tau) = F_XUV exp(-tau) of a_tau = 0.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `energy_absorbed_equals_beam_loss[opa_pf_10]`, `[opa_pf_1]` | `q_abs*dr` of the cell, recovered from the returned `heat` and `q` as `heat/q` | The energy the beam loses across the cell, `INT F(tau_out) [1 - exp(-dtau)] dE` | 1e-5 relative |
| `HI_ionizations_equal_beam_photon_loss[opa_pf_10]`, `[opa_pf_1]` | `P_HI n_HI dr` | The photons the beam loses across the cell times the H I share `dtau_HI/dtau` | 1e-5 relative |

The identity is exact only for a thin cell (the code illuminates a cell with
the field at its inner face), so the slab is built with a cell optical depth of
6.3e-8 and the measured deviation is 3e-8, two orders below the tolerance.

Added in Phase 1 batch 2a, not at 35d9dd5: **RED before item 2a-OPAC** on the
two `opa_pf_10` assertions, both low by exactly the factor f = 10 (9.4287e-04
against 9.4287e-03 erg cm^-2 s^-1, and 2.7707e+07 against 2.7707e+08
cm^-2 s^-1), and GREEN on the two `opa_pf_1` assertions in both states.
**GREEN on all four** with the item's change in place.

### 10. Ly-alpha band transmission (`lya_band_transmission.f90`)

Origin: Phase 1 item 2a-LYA, extended by item B4-6c of
`docs/PLAN_20260906_rev2.md`. Production routines: `fuv_lw_photon_field`
(`src/modules/radiation/util_ion_eq.f90`, the shared 912-2304 A beam),
`lya_stellar_beam_transmission`,
`lya_stellar_beam_transmission_cell_mean` and
`lya_line_center_optical_depth` (`src/modules/radiation/lya_rt.f90`),
`water_photolysis_rate`
(`src/modules/lower_atmosphere/water_photolysis.f90`).

Band B2 of the photolysis set is the H I Ly-alpha resonance line, so the
stellar flux that drives it reaches a cell through the atomic hydrogen column
above that cell. The beam is resonantly scattered rather than destroyed, so its
transmission is not `exp(-tau)` but the fraction of the stellar profile whose
wings penetrate to that depth, `erfc(x1/(sqrt(2) X_s))`, the expression
`lya_rt.f90` already builds its penetrating stellar beam from (and whose
provenance is recorded there). A photolysis rate is what the molecules of a cell undergo
averaged over the cell, so the band must carry the MEAN of that transmission
over the cell and not its value at a face, exactly as the continuum side of
the same rate already does. The size of that correction is set by the change
of the erfc argument across the cell and not by the cell's line depth.

Three H I columns are used, all isothermal on a uniform grid with the H2O and
OH densities set to zero so that every continuum depth vanishes and each rate
is the free-streaming rate times the line factor alone: an extended
exponential layer (scale height 0.3 R_p, column 3.0e18 cm^-2), the
molecular base of a hot Uranus (scale height 0.05 R_p, base line-centre depth
1e7, so that one cell holds most of the depth above it), and a graded column
whose cell line-centre depths run from 1 at the top to 1e3 at the base. A
fourth configuration adds a flat H2O density to the extended layer for the
photon budget.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `lya_band_free_streams_without_hydrogen` | `j_H2O(B2)` at the base with `n_HI = 0` | `sigma_H2O(B2) N_ph(B2)`, the thin limit | 1e-12 relative |
| `lya_band_carries_beam_cell_mean_{thin_layer,molecular_base,graded_1_to_1e3}` | max over cells of `j_H2O(j,B2)`/free-stream against the cell mean of the beam over that cell's own line-centre depth | an independent fine quadrature of the production transmission (composite Simpson, 8192 intervals in `u = sqrt(tau)`) | 1e-6 relative |
| `lya_band_face_over_mean_{thin_layer,molecular_base,graded_1_to_1e3}` | largest ratio of the star-ward face value to that mean, reference-only | above 1 | strict |
| `lya_band_test_column_cm2` | the H I column the extended layer builds | strictly positive | strict |
| `lya_band_base_transmission_below_half` | the base ratio | 0 | 0.5 absolute, the physical statement that this column suppresses the band |
| `lya_band_leaves_other_band_unchanged_{LW,B3,B4}` | max relative change of each other band between `n_HI = 0` and the layer | 0: only B2 has H I as a line absorber | 1e-12 absolute |
| `lya_band_transmission_falls_inward` | max inward increase of `j_H2O(B2)` over the free-streaming rate | 0: the column only grows inward | 1e-12 absolute |
| `lya_band_photon_budget_cell_mean` | photons the continuum absorbers remove from B2, summed over the column, read back from the production rate | the exact loss, a fine integration of `T_lya exp(-tau_c)` over each cell | the covariance bound of the product-of-means form, `0.25 dtau_c (T_out - T_in)/<T>` maximized over the column, measured as 1.399e-4 |
| `lya_band_photon_budget_face_error` | the same budget with the face value of the line factor, reference-only | its error exceeds that bound | strict |

The band multiplies two cell means, the line factor and the continuum factor,
rather than taking the mean of their product, because the continuum mean is
formed inside `water_photolysis_rate`. Both factors fall monotonically across
a cell, so the product of the means is a lower bound on the mean of the
product, short of it by their covariance; the budget assertion is that bound.
MEASURED on the budget column: the product of means is 6.0e-6 relative below
the exact loss, against a bound of 1.4e-4, while the face rule is 1.2e-2
above it.

**RED before item 2a-LYA** on three of the ten assertions the driver then
had: the band free-streamed to the base, so the transmission ratio was 1.0
against the 0.2449 of the test column. **RED before item B4-6c** on four of
the sixteen it then had: the band read the beam at the cell's star-ward face,
which is 1.095 times the cell mean in the extended layer, 1.005 times it in
the graded column and **18.2 times** it in the base cell of the molecular
base, and the column photon budget stood 1.2 per cent above the beam's own
loss. **GREEN on all sixteen** with the cell mean in place; the base ratio of
the extended layer is then 0.2066. The driver carries fifteen assertions now
(MEASURED 2026-09-18): the FUV band list is LW, B2, B3, B4, so there is one
`lya_band_leaves_other_band_unchanged_*` row fewer than when it was written.

## Summary at 35d9dd5 (Phase 0)

Seven groups, 41 assertions: 15 pass, 26 fail. Failing: 3 of group 1, 9 of
group 2 (of 10), 2 of group 3 (of 11), 6 of group 4, 5 of group 5, 1 of
group 6 (of 3). Group 7 (3 assertions) was green by construction; groups 8, 9
and 10 did not exist. Every failure is a defect the Phase 0 audit identified;
the tests are expected to stay red until the corresponding Phase 1 to Phase 4
item lands.

## Summary after Phase 1 batch 2a (measured 2026-09-05)

Ten groups, 83 assertions: 69 pass, 14 fail, 5.1 s. Measured on the working
tree with batch 2a landed and the LLF item of batch 2b already in it (see the
note below the table). Batch 2b adds further drivers, so the group count above
ten belongs to that batch's record, not to this one.

| Group | Assertions | Pass | Fail | State |
|---|---|---|---|---|
| 1 LLF signal speed | 3 | 3 | 0 | turned green |
| 2 ROE star state | 30 | 30 | 0 | turned green, and grew from 10 assertions to 30 |
| 3 H2 channel split | 11 | 9 | 2 | unchanged; item A8 |
| 4 H3+ cooling limits | 6 | 0 | 6 | unchanged; Phase 3 |
| 5 H2 rovibrational identity | 5 | 0 | 5 | unchanged; Phase 3 |
| 6 Photoevent energy ledger | 3 | 2 | 1 | unchanged; item A5, Phase 3 |
| 7 Radiation exchange ledger | 3 | 3 | 0 | green by construction |
| 8 Wind-AE bridge composition | 8 | 8 | 0 | new with item 2a-WAE, green |
| 9 Opacity model P ledger | 4 | 4 | 0 | new with item 2a-OPAC, green |
| 10 Ly-alpha band transmission | 10 | 10 | 0 | new with item 2a-LYA, green |

What turned green and why:

- **Group 2** (ROE), by item A2 of batch 2a: the estimator of
  `docs/a2_roe_interface.md`, with the corrected two-rarefaction exponent and
  two-shock densities, the vacuum test taken first as a status, and the HLLE
  flux on any face without an admissible star state.
- **Groups 8, 9 and 10** are new drivers that came with items 2a-WAE, 2a-OPAC
  and 2a-LYA, each written RED against the defect its item removes and green
  afterwards. Their references and tolerances are sections 8, 9 and 10 above.
- **Group 1** (LLF) is green from item A1, which is a batch 2b item and was
  already in the working tree when this was measured; it is not a batch 2a
  result.

The 14 remaining failures are the same defects Phase 0 identified, each waiting
on its own item: the Chung detector inversion (group 3, item A8), the H3+
cooling limits and the H2 rovibrational energy identity (groups 4 and 5,
Phase 3), and the 19.72 eV an H2 double photoevent at 80 eV leaves unassigned
(group 6, item A5, Phase 3).

### 12. Helium level cooling ledger (`helium_level_cooling_ledger.f90`)

Origin: the two helium findings of
`docs/audit_20260905/source_audit/group_B_report.md`
("eval_cool applies the He I ground-state collisional coefficients to the total
neutral helium ...", "He II recombination cooling is charged at alpha_B ...").
Production routines exercised: `eval_cool`
(`src/modules/radiation/util_ion_eq.f90`) and the coefficient routines of
`Cool_coeff.f90` it assembles.

References. The state vector carries one neutral-helium column, n(He I), and it
CONTAINS the 2^3S metastable (`composition.f90`, where the difference is
taken). The two levels are 19.8 eV apart: the 24.6 eV ionization potential and
the Cen (1992) collisional-excitation coefficient are the ground singlet's, the
4.8 eV potential (Black 1981), the 10830 A excitation and the 2^3S -> 2^1S /
2^1P conversions are the metastable's. Each density must therefore appear once,
under its own level's coefficient. For the recombination channel the electron
gas loses kT per recombination the balance actually performs, so the cooling
coefficient is kT times the coefficient the balance removes He+ with: with the
metastable tracked that is `rcheiiB + rcheiTR` (`rec_HeII_11S + rec_HeII_23S`,
Oklopcic & Hirata 2018), without it the case-B coefficient of the
`Atomic rate set:` key.

Configuration: a three-cell slab at 1e4, 2e4 and 5e4 K with
n(1^1S) = 1e5 cm^-3 and n(2^3S) = 1e2 cm^-3 (the ratio 1e-3 of the finding; the
largest measured on `wasp_full` is 2.2e-4). The first six assertions run with
n(H I) = n(He II) = 0, so the two breakdown columns hold the helium terms alone
and the comparison carries no cancellation. The three temperatures are there
because the 24.6 eV coefficients are exponentially suppressed below about
2e4 K, so a double count of the metastable is invisible at 1e4 K (1.2e-11
relative) and is of order n(2^3S)/n(1^1S) at 5e4 K (6.0e-4 relative).

All tolerances are 1e-12 relative: every comparison is an identity of the
assembly, not a fit. The 24.6 eV prefactor is written in the test with the same
single-precision literal `3.940e-11` the assembly uses, so the comparison is of
the level the coefficient is charged to and not of the precision of the
constant (that duplicated literal is a separate open item of rev3 section 10.3).

| Assertion | Measures | Reference | Tolerance |
|---|---|---|---|
| `helium_collisional_ionization_by_level_T{10000,20000,50000}K` | `cool_chan(:,2)`, the collisional-ionization column | `n_e (24.6 eV a_ion_HeI n(1^1S) + 4.8 eV ci_HeI23S n(2^3S))` | 1e-12 relative |
| `helium_collisional_excitation_by_level_T{10000,20000,50000}K` | `cool_chan(:,4)`, the He I collisional-excitation column | `n_e (lambda_coex_HeI n(1^1S) + (10830 A + 0.80 q31a + 1.40 q31b) n(2^3S))` | 1e-12 relative |
| `he_ii_recombination_cooling_per_event_metastable_T{10000,20000,50000}K` | `rec_cool_HeII / (k T)` with the metastable tracked | `rec_HeII_11S + rec_HeII_23S`, the two channels the balance removes He+ with | 1e-12 relative |
| `recombination_cooling_column_T{10000,20000,50000}K` | `cool_chan(:,1)`, the recombination column | `n_e (lambda_rec_HII n_HII + lambda_rec_HeII n_HeII + lambda_rec_HeIII n_HeIII)` | 1e-12 relative |
| `he_ii_recombination_cooling_per_event_case_b_T{10000,20000,50000}K` | `rec_cool_HeII / (k T)` with the metastable off | `alpha_rec_HeII_B`, the coefficient the balance keeps there | 1e-12 relative |

**RED before item 2b-HEI** on nine of the fifteen: the six level assertions (the
ground-state coefficients were charged to the summed neutral helium, so the
metastable paid them as well as its own) and the three metastable-tracked
recombination assertions (2.807e-13 charged against 3.640e-13 removed at
1e4 K). The three `recombination_cooling_column` assertions and the three
case-B ones passed before and after; they are the controls. **GREEN on all
fifteen** after the item.

What this test does NOT cover: with `He_rec_coupling: True` the balance further
weights the ground-term capture by the fraction of its 24.6 eV photons that do
not go back into He I and adds the 0.25 alpha_B excited-singlet capture. That
weight is a function of the absorber densities, not of temperature, and is not
carried into the cooling (the reason is written at `alpha_rec_HeII_total` in
`Cool_coeff.f90`).

### 13. Lyman-Werner rate as a cell mean (`lyman_werner_cell_mean.f90`)

Origin: the third finding of `docs/audit_20260905/source_audit/group_C_report.md`
and rev 3 section 10.3, item 2b-KLW.
Production routines: `fuv_lw_photon_field` (`src/modules/radiation/util_ion_eq.f90`),
`lyman_werner_dissociation_rate` and `lyman_werner_dissociation_rate_cell_mean`
(`src/modules/lower_atmosphere/lyman_werner.f90`),
`h2_lw_dissociation_cross_section`
(`src/modules/lower_atmosphere/h2_self_shielding_table.f90`).

`calc_column_dens_one` accumulates from the top down, so `NH2col(j)` already
contains the whole of cell j: it is the column at that cell's inner face, and
`NH2col(j+1)` is the column at its star-ward face. The dissociation cross
section falls by four decades across the self-shielding transition and fastest
at the H2 front, where one cell spans a large fraction of a decade of column,
so the inner-face value applied to the whole cell is low in one direction. The
H2O and OH continua of the same beam already take the exact cell mean.

The configuration is an isothermal (1300 K) exponential H2 layer on a uniform
grid of 40 cells whose width is 2.5 scale heights, so the column grows by
e^2.5 across each cell and the cross section falls by 1.12 decades across the
deepest one. The oxygen chemistry is off, so the continuum depth of the band
vanishes and the line term is what is measured. The reference is the mean of
the production local rate over the cell's column interval by the midpoint rule
on 1000 steps.

| Assertion | Measures | Reference | Tolerance |
|---|---|---|---|
| `lw_deepest_cell_drop_above_0p8_dex_margin` | the configuration is in the steep regime | the cross section falls by more than 0.8 dex across the deepest cell | strictly positive margin |
| `lw_rate_is_the_cell_mean` | worst relative departure of `k_lw(j)` from the reference cell mean, over every cell but the outermost | 0 | 1e-3 absolute |
| `lw_inner_face_rate_below_the_mean` | the inner-face rate divided by the cell mean, deepest cell | 0.31406, MEASURED here | factor 1.01 |
| `lw_inner_face_rate_never_above_the_mean` | the same ratio, largest over all cells, less 1 | 0 | 1e-12 absolute |
| `lw_cell_mean_of_a_uniform_column` | the mean over a cell across which the column does not grow | the local rate at that column | 1e-12 relative |
| `lw_synthetic_cell_spans_one_decade` | the bisected synthetic cell | a cross-section ratio of exactly 10 | 1e-6 relative |
| `lw_synthetic_cell_mean_matches_reference` | the cell mean of that synthetic cell | the 1000-point reference integral | 1e-3 relative |

**RED before item 2b-KLW** on `lw_rate_is_the_cell_mean` (worst cell 0.686, i.e.
the rate was 69 per cent below the mean of the cell it was applied to) and on
`lw_inner_face_rate_below_the_mean` (1.000, the rate being the inner-face value
itself). **GREEN after**: worst departure 3.3e-6, inner-face ratio 0.31406, so
the deepest cell of this configuration was carrying a rate 3.18 times too
small. The other five assertions pass before and after; they are the controls
on the quadrature and on the configuration.

### 14. Carbon electron-impact ionization rows (`voronov_carbon_ionization.f90`)

Origin: the C I finding of `docs/audit_20260905/source_audit/group_B_report.md`
and rev 3 section 10.2, item 2c-CI.
Production routines: `ion_coeff_CI_range` and `ion_coeff_CII_range`
(`src/modules/radiation/Cool_coeff.f90`).

Voronov (1997, ADNDT 65, 1) fits every electron-impact ionization rate
coefficient by `k(T) = A (1 + P sqrt(U)) U^K exp(-U)/(X + U)` with
`U = dE/(k_B T)` and one row of `(dE, P, A, X, K)` per ion in his Table 1. The
row of neutral carbon is `(11.3 eV, 0, 6.85e-8, 0.193, 0.25)`: `P` is zero, so
the bracket is absent. The row of C+ is `(24.4 eV, 1, 1.86e-8, 0.286, 0.24)`
and does carry it. A transcription is an identity, not a fit, so the tolerance
is round-off; the reference constants are written with the same literal kinds
the production rows use, so what is compared is the algebraic form and the
constants and not the rounding of a default-real literal. `Cool_coeff.f90`
writes the ionization potentials to four figures throughout (C I 11.26 eV,
C II 24.38 eV, N I 14.53 eV, O I 13.62 eV) where the table prints three, and
the reference follows the code in this; the driver prints the ratio Voronov's
own 11.3 eV would give as a `note:` line (0.791 at 2000 K, 0.952 at 1e4 K,
0.988 at 5e4 K), which is a diagnostic and not an assertion.

C II is asserted alongside C I because the two rows share the one functional
form: it is the control showing that the form is right and that only the C I
constants were in question.

| Assertion | Measures | Reference | Tolerance |
|---|---|---|---|
| `c_i_electron_impact_ionization_T2000K` | `ion_coeff_CI_range` at 2000 K | `6.85e-8 U^0.25 exp(-U)/(0.193 + U)` = 1.256337e-37 cm^3 s^-1 | 1e-12 relative |
| `c_i_electron_impact_ionization_T10000K` | the same at 1e4 K | 2.076865e-14 cm^3 s^-1 | 1e-12 relative |
| `c_i_electron_impact_ionization_T50000K` | the same at 5e4 K | 2.274529e-09 cm^3 s^-1 | 1e-12 relative |
| `c_ii_electron_impact_ionization_T2000K` | `ion_coeff_CII_range` at 2000 K | `1.86e-8 (1 + sqrt(U)) U^0.24 exp(-U)/(0.286 + U)` = 2.040039e-70 cm^3 s^-1 | 1e-12 relative |
| `c_ii_electron_impact_ionization_T10000K` | the same at 1e4 K | 4.737535e-21 cm^3 s^-1 | 1e-12 relative |
| `c_ii_electron_impact_ionization_T50000K` | the same at 5e4 K | 5.589994e-11 cm^3 s^-1 | 1e-12 relative |

**RED before item 2c-CI** on the three C I assertions: the row read
`(P, X, K) = (0.193, 0.25, 0.25)`, X's published value taken as P and K's as X,
which adds a `(1 + 0.193 sqrt(U))` factor that neutral carbon does not have and
moves the pole of the denominator. Measured 3.213429e-37, 3.510707e-14 and
2.924777e-09 cm^3 s^-1, i.e. 2.558, 1.690 and 1.286 times the published rate.
**GREEN after.** The three C II assertions pass before and after.

### 15. Atomic mass and radius constants (`atomic_mass_and_radius_constants.f90`)

Origin: findings 1 and 3 of `docs/audit_20260905/source_audit/group_D_report.md`
and decision 14 of `docs/development_plan_20260905_rev3.md` section 10.5
(item 2c-MASS).
Production data: `bsp_mass` (`src/modules/init/species_table.f90`), the
astronomical constants of `global_parameters`
(`src/modules/init/parameters.f90`) and the copies the Wind-AE front end holds
(`src/modules/wind_ae/wae_exhale_input.f90`).

`bsp_mass` is stated in units of the density normalization
`mu = 1.67353284e-24 g`, the hydrogen ATOM and not the atomic mass unit, so
helium must carry `m_He/m_H = 6.6464731e-24/1.67353284e-24 = 3.97152` and
HeH+ that plus one hydrogen. Separately, `Planet radius [R_J]` is read twice
from the same file, by `input_read` and by the Wind-AE front end, and each
multiplies it by its own Jupiter radius: one input line must name one planet,
so the two constants are one number, the IAU 2015 nominal equatorial radius
7.1492e9 cm (Prsa et al. 2016, AJ 152, 41, Table 1). The same holds for the
two masses.

| Assertion | Measures | Reference | Tolerance |
|---|---|---|---|
| `helium_atom_mass_in_hydrogen_atom_units` | `bsp_mass(He I)/bsp_mass(H I)` | 3.9715224, the CODATA mass ratio | 1e-4 relative |
| `helium_ii_mass_equals_helium_i` | the He II row | the He I row | 1e-12 relative |
| `helium_iii_mass_equals_helium_i` | the He III row | the He I row | 1e-12 relative |
| `helium_metastable_mass_equals_helium_i` | the He 2^3S row | the He I row | 1e-12 relative |
| `heh_plus_mass_closes_on_its_two_nuclei` | the HeH+ row | `bsp_mass(He I) + bsp_mass(H I)` | 1e-15 relative |
| `jupiter_radius_iau_2015_nominal_equatorial` | `global_parameters` `RJ` | 7.1492e9 cm | exact |
| `jupiter_radius_of_wind_ae_front_end` | the front end's `RJ` | `global_parameters` `RJ` | exact |
| `jupiter_mass_of_wind_ae_front_end` | the front end's `MJ` | `global_parameters` `MJ` | exact |
| `solar_mass_of_wind_ae_front_end` | the front end's `MSUN` | `global_parameters` `Msun` | exact |
| `astronomical_unit_of_wind_ae_front_end` | the front end's `AU` | `global_parameters` `AU` | exact |
| `helium_mass_of_binary_element_diffusion` | `m_He_amu` of `binary_element_diffusion` | `bsp_mass(He I)` | exact |
| `hydrogen_mass_of_binary_element_diffusion` | `m_H_amu` of the same module | `bsp_mass(H I)` | exact |
| `helium_mass_of_wind_ae_ic_writer` | `ICW_mHe_over_mH` of `wae_ic_writer` | `bsp_mass(He I)/bsp_mass(H I)` | 1e-5 relative |
| `helium_mass_of_analytic_lower_column` | `m_He_over_m_H` of `lower_column` | `bsp_mass(He I)/bsp_mass(H I)` | exact |

**RED before item 2c-MASS** on five assertions: the helium ratio measured
4.0 against 3.9715224 (0.72 per cent high), `RJ` measured 6.9911e9 cm (the
volumetric mean radius, 2.26 per cent below the equatorial one), and the front
end's `RJ`, `MJ` and `MSUN` each differing from the global constant.
**GREEN after** those nine.

The last five rows are item 2c-MASS2, which closed the copies of the helium
mass that the table's move left behind, and the second copy of the
astronomical unit. MEASURED, RED before it on four of the five: the binary
element diffusion, the Wind-AE IC writer and the analytic lower column each
measured 4.0 against the table's 3.9715 (0.72 per cent high), and the front
end's `AU` measured 1.49598e13 cm against 1.495978707e13 (8.6e-8 relative
high); `hydrogen_mass_of_binary_element_diffusion` passed already, and is
there because the closure needs both masses and not only helium. GREEN after
all five. The one loose tolerance is the Wind-AE IC writer: the Wind-AE
tree is compiled into the standalone `wind_ae_ic.x` and cannot read the
species table, so it forms the ratio from the two CODATA atom masses it hands
the solver, 6.6464790722e-24/1.67353284e-24 = 3.9715259, which agrees with
the table's 3.9715 to 6.5e-6. The fourth copy of the helium mass, the
background collision partners of the carrier transport, is asserted in
`src/tests/species_masses/`: that module's closure reaches the ionization
sweep and therefore needs LAPACK, which the drivers here do not link.

### 16. The Jupiter radius of the Python tools (`jupiter_radius_python_tools.py`)

Origin: the same decision 14, item 2c-MASS2. `Planet radius [R_J]` is read by
the code and, separately, by the Python tools that prepare its input
(`src/utils/run_lower.py`, `vulcan_to_base.py`, `lower_profile_schema.py`,
`EXHALE_interface_state.py`, `vulcan_driver.py`) and by those that analyse
its output (`EXHALE_plots.py`, `eta_approx.py`, `exhale_transit_lib.py`,
`roche_recon.py`). Each held its own copy, so a regenerated `base.inp` or
lower-atmosphere profile would have been built for a planet 2.26 per cent
smaller than the one the code relaxes.

The one Python definition is now `RJ_CM` of `examples/exhale_io.py`, and this
test asserts that every listed tool resolves its own `R_J` to that value. The
assertion is made on the VALUE the tool ends up with, not on the presence of
an import line, so a tool that imported the constant and then overwrote it
would fail here. The three tools that cannot be imported in a test process
(`EXHALE_plots.py` reads `./input.inp` at import, `eta_approx.py` asks the
terminal for a radius, `roche_recon.py` binds its `R_J` only inside its
self-test) have their single module-level assignment read as text and
evaluated with `RJ_CM` bound.

| Assertion | Measures | Reference | Tolerance |
|---|---|---|---|
| `jupiter_radius_of_parameters_f90` | `RJ` parsed out of `parameters.f90` | 7.1492e9 cm, IAU 2015 nominal equatorial | exact |
| `jupiter_radius_of_the_python_definition_exists` | that `examples/exhale_io.py` defines `RJ_CM` | it does | exact |
| `jupiter_radius_of_the_python_definition` | `exhale_io.RJ_CM` | the Fortran value | exact |
| `jupiter_radius_of_<tool>` (eight rows) | the radius each of the other tools runs with | `RJ_CM`, or `RJ_CM/100` for the two that work in metres | exact, 1e-15 for those two |
| `jupiter_radius_of_vulcan_driver` | that it imports `RJ` from `run_lower` | an import, not a literal | exact |
| `volumetric_mean_jupiter_radius_absent_from_the_tools` | the nine sources scanned for 6.9911e9 or 6.9911e7 | no occurrence | exact |

MEASURED, RED before item 2c-MASS2 on ten of the twelve assertions it then
made: `jupiter_radius_of_parameters_f90` passed, because item 2c-MASS had
already moved the Fortran, and `jupiter_radius_of_vulcan_driver` passed
because that tool already imported its radius from `run_lower.py` (which held
the stale copy). Every other tool carried 6.9911e9 (or 6.9911e7 where it works
in SI), 2.26 per cent below the equatorial radius, and there was no `RJ_CM` to
import. GREEN after: all thirteen (the thirteenth is the existence assertion's
companion, which only runs once a definition exists).

### 17. Ionization thresholds written once (`ionization_threshold_turn_on.f90`, `threshold_literal_uniqueness.py`)

Origin: item 2c-THR of Phase 1 batch 2c. The photon grid puts every
ionization threshold on a bin edge, the photoelectron energy is `h nu - e_th`,
the secondary-ionization targets and the collisional-ionization cooling remove
`e_th` per event, and the photoionization cross sections turn on there. Those
are one physical energy each, and the source has to state it once. It did not:
`parameters.f90` carried `e_th_HeI = 24.6` while the He I cross section turned
on at the Verner et al. fit's own `E_th = 24.59`, so the band
[24.59, 24.60] eV, 0.115 per cent of the He I photoionization rate of the
default power law, was integrated as zero; the He 2^3S cross section turned on
at 4.78 eV, 0.02 eV below the grid floor `e_th_HeTR = 4.80`; and every
hydrogenic cross section turned on at `0.99999` of its threshold.

The values are now the measured ionization energies (NIST ASD; for H2 the NIST
Chemistry WebBook adiabatic value), the metastable being the He I energy minus
the 1s2s 3S1 excitation energy.

| Assertion | Measures | Reference | Tolerance |
|---|---|---|---|
| `<species>_threshold_is_the_measured_potential` (four rows) | `e_th_HI`, `e_th_HeI`, `e_th_HeII`, `e_th_H2` | 13.598434599, 24.587389, 54.417765, 15.425927 eV | 1e-9 relative |
| `he_metastable_threshold_is_singlet_minus_exc` | `e_th_HeTR` | `24.587389 - 19.819614` eV | 1e-12 relative |
| `<species>_threshold_in_erg` (four rows) | `e_th_*_erg * erg2eV` | the eV constant | 1e-15 relative |
| `h2_<channel>_channel_threshold` (three rows) | `h2_channel_threshold` | `e_th_H2`, `e_th_H2_di`, `e_th_H2_dd` | exact |
| `<species>_cross_section_zero_below_threshold` (five rows) | the cross section 1e-6 eV below its threshold, through `photoion_sigma` | 0 | exact |
| `<species>_cross_section_positive_above_threshold` (five rows) | the same 1e-6 eV above | strictly positive | -- |
| `he_i_ates_fit_*`, `constant_opacity_model_*` | the opt-in He I fit and the constant opacity model at the same two energies | as above | exact |
| `h_n2_edge_is_a_quarter_of_the_h_i_threshold` | `e_th_HI_n2` of `J_inc.f90`, the Balmer-continuum edge | `e_th_HI/4` | exact |
| `threshold_literals_written_once` | production `.f90` lines, comments and string bodies removed, carrying the numeric value of one of these potentials outside its declaration | none, beyond the listed exceptions | exact |

The exceptions of the last assertion are in `threshold_literal_exceptions.txt`
with a reason each: numbers that are coefficients of a published fit (the dE of
the Voronov 1997 electron-impact ionization rates, the Rydberg prefactor of the
Oklopcic/Lampon collision-strength rates, the scaling energy of the Yan et al.
1998 sum-rule tail) or rows of tabulated data. The file also supports an entry
marked `open`, for a duplicate in a file the item writing the test did not own;
such an entry is printed as a DIAGNOSTIC instead of failing, so that it stays
visible. There is none at present.

MEASURED, RED before item 2c-THR: 11 of the 22 assertions the driver could
then make failed (the four erg constants did not yet exist), and the literal
scan found 13 duplicate statements of these potentials in five files. GREEN
after: all 27 driver assertions and the scan.

### 18. H2 infrared line populations from the shared partition function (`h2_infrared_line_populations.f90`, `h2_partition_sum_uniqueness.py`)

The H2 rovibrational ladder of Roueff et al. (2019) fixes three things at once:
the internal energy the caloric equation of state carries, the free energy the
H + H <-> H2 equilibrium constant is built from, and the LTE level populations
that weight the H2 quadrupole lines of `molecular_infrared_cooling`. They are
moments of one Boltzmann sum, so the source forms that sum once, in
`h2_partition_function` of `caloric_eos.f90`.

| Assertion | Measures | Reference | Tolerance |
|---|---|---|---|
| `h2_line_emission_from_shared_Q_at_<T>` (five rows, 100 to 8000 K) | `h2_line_emission_lte`, the optically thin LTE power per H2 molecule | the same line sum rebuilt from `h2_partition_function` | 1e-12 relative |
| `h2_partition_sum_uniqueness` | production `.f90` lines under `src/modules`, comments and string bodies removed, naming `h2_lev_T` or `h2_lev_g` outside the module that declares them and the module that sums them | none | exact |

The temperatures are nodes of the module's own log-spaced evaluation grid, so
the table lookup returns the tabulated value and the assertion tests the
populations rather than an interpolation.

MEASURED, before item B3b-IR: the five emission assertions already passed, and
the level sum the cooling module then carried agreed with
`h2_partition_function` to the bit at every one of the five temperatures. The
`exp(-min(E/T, 700))` guard of that sum binds only where a level has
`E/T > 700`, i.e. below `51966/700 = 74.2` K, which is below the whole range
over which the layer radiates; where it does bind it replaces a term that
underflows to zero by one of order 1e-304, which cannot move a sum whose
ground term is 1. The uniqueness scan is what was RED: it found the second sum
at `molecular_infrared_cooling.f90:214`. GREEN after: the scan, and the five
emission assertions unchanged to every printed digit.

Two diagnostics are printed rather than asserted. The first is the ratio of
that retired guarded sum to `h2_partition_function`, the transition record of
what the de-duplication changed; it is 0 to the bit at all five temperatures.
The second is the table ceiling: `h2_line_emission_lte` is clamped above the
8000 K end of its grid, and at 10000 K it returns 1.181e-18 erg s^-1 against
the direct sum's 1.579e-18, 25 per cent low. The module header states that the
H2 line sum "has no such limit; it is evaluated directly from the level
ladder", which is true of the sum and not of the tabulation it is stored in.

### 19. Ly-alpha escape probability (`lya_escape_probability.f90`)

Ly-alpha scatters with partial redistribution (R_II-A, Hummer 1962): a photon
shifts by about one Doppler width per scattering and cannot be re-emitted
directly in the far wing, so escape from a thick static slab is a random walk
in frequency, and the mean number of scatterings before escape is proportional
to the line-centre optical depth and independent of the Voigt parameter. The
published solution is Neufeld (1990), ApJ 350, 216, eq. (3.27) with his
continuum destruction opacity set to zero, which for a source in the mid-plane
of a slab of line-centre optical half thickness `B` is Harrington (1973),
MNRAS 162, 43, eq. (40), `<N> = (4 sqrt(6)/pi^2) u_2 B = 0.909316 B` with
`u_2` Catalan's constant. `<N>` counts absorptions, so a generated photon is
emitted `1 + <N>` times and the escape chance per emission that closes the
two-level budget `n2p A beta = P` is `beta = 1/(1 + <N>)`, which is exact in
the thin limit as well and needs no interpolation.

| Assertion | Measures | Reference | Tolerance |
|---|---|---|---|
| `lya_slab_midplane_factor_is_catalan` | `lya_slab_source_position_factor(1/2)` | `u_2 = 0.9159656`, Harrington below his eq. (40) | 1e-12 relative |
| `lya_slab_factor_face_symmetry` | `Phi(0.3) - Phi(0.7)` | 0: the slab does not know which face is which | 1e-14 absolute |
| `lya_slab_factor_on_a_face_vanishes` | `Phi(0)` | 0 | exact |
| `lya_slab_factor_matches_its_series` | `Phi` at five arguments | `sum_{n odd} sin(n pi xi)/n^2` summed to 4e6 terms in the driver | 1e-9 relative |
| `lya_slab_scatterings_harrington_eq40_{1,2,3}` | `lya_slab_scatterings_before_escape(tau,tau)` at `tau = 1e8, 1e6, 1e4` | `(4 sqrt(6)/pi^2) u_2 tau`, written out in the driver | 1e-12 relative |
| `lya_wing_escape_is_one_over_one_plus_N_{1,2,3}` | `lya_wing_escape_probability` at the same depths | `1/(1 + 0.909316 tau)` | 1e-12 relative |
| `lya_slab_scatterings_face_symmetry` | `<N>(1e7, 3e6) - <N>(3e6, 1e7)` | 0 | 1e-8 absolute |
| `lya_wing_escape_thin_limit`, `_small_depth` | `beta` at zero and at `tau = 1e-9` | 1 and `1/(1 + 0.909316e-9)` | exact, 1e-12 relative |
| `lya_wing_escape_never_exceeds_one` | `beta` over fourteen decades of depth | at most 1 | strict |
| `lya_face_partition_sums_to_one` | the two shares of `lya_starward_escape_fraction` | 1, Neufeld's own statement under his eq. (2.25) | 1e-14 relative |
| `lya_face_partition_midplane_is_half`, `_source_on_the_face` | the same at a mid-plane source and at a source lying on the planet-ward face | 1/2 and 0 | 1e-14 relative, exact |
| `lya_closure_carries_the_published_escape` | `Jint n1s (A beta + D)` from `jlya_escape_prob` over an exponential H I column of base depth 1e8 | `(2 h nu^3/c^2)(g1s/g2p) P (1 - beta)` with `beta` the published escape probability | 1e-13 absolute |
| `lya_wing_domain_cells_{seen,out_of_domain}`, `lya_wing_domain_threshold` | `lya_wing_domain_record` over that column | the cells of the column with `(a tau)^(1/3) < 10` (Neufeld section Va), counted in the driver | 1e-14 relative |

The closure assertion is written as the equality of the two sides rather than
by recovering `n2p`, because `1 - beta` vanishes in the thin outer cells and
dividing by it would test round-off instead of the closure.

**RED before item LYA-BETA** on `lya_closure_carries_the_published_escape`: the
field was closed with the one-flight wing form `pi^(-1/4) sqrt(a/tau)`, the
complete-redistribution value, and the assertion measured 4.79e2 against a
tolerance of 1e-13. GREEN after, at 2.19e-16. The other assertions are
statements about routines that item introduced and have no earlier value.

Two further things the driver measured while it was written. The first is that
`beta_tot` was assembled as `1 - (1 - beta_esc)(1 - beta_sob)`, which forms
`1 - x` for `x` of order 1e-8 at the base and recovers it by subtraction,
losing eight of its digits; it is now the union of two independent chances,
`beta_esc + beta_sob(1 - beta_esc)`, and the closure assertion went from 1.4e-9
to 2.2e-16. The second is the domain record itself: in the test column, 49 of
64 cells fail `(a tau)^(1/3) > 10`, which is the same picture the stored
`wasp_full` column gives.

### 20. Low-Mach scaling of the ROE velocity jump (`roe_low_mach_velocity_jump.f90`)

Origin: item L25 step 2 of the LHS 1140 b stationary campaign
(`docs/lhs1140b_stationary_L25_20260916.md` section on the input key).
Production routines: `Num_flux` and `Phys_flux`
(`src/modules/flux/Num_Fluxes.f90`), taking the ROE branch because `flux` of
`global_parameters` is set to `'ROE'`, and the module variable
`low_mach_velocity_jump`, which the input key `Low Mach velocity jump:`
writes and which defaults to `.false.`.

The correction is Rieper (2011), J. Comput. Phys. 230, 5263, eq. 3.15 with
the local Mach number of the Roe averages, `dvel -> min(Ma_Roe, 1) dvel`,
which in one dimension is `|U_Roe|/a_Roe`. The reference of every assertion
is the analytic three-wave Roe flux (Toro 2009, eq. 11.29) written out in the
driver from the Roe averages of the pair, so that the returned flux is
compared with an independent evaluation and not with itself. The flux is
AFFINE in the factor and only the two acoustic coefficients carry the
velocity jump, which is what lets assertion group (3) below separate the
scaled part from the untouched one. The state pairs are chosen so that
neither eigenvalue is transonic and the rarefaction entropy correction does
not fire; that the analytic reference reproduces the returned flux to
rounding is itself the check of that choice.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `roe_key_off_row_{m,p,e}` | The three rows of the flux of WL=(1,0.30,1), WR=(1.2,0.34,1.1) with the key off | the analytic Roe flux at factor 1 | 1e-12 relative |
| `roe_key_off_repeatable_row_{m,p,e}` | The same flux taken again after a call with the key on | the first value: the branch keeps no state | exact |
| `roe_supersonic_face_mach` | `Ma_Roe` of WL=(1,5,1), WR=(1.2,5.2,1.1) | at least 1 | at least |
| `roe_supersonic_unchanged_row_{m,p,e}` | The same flux with the key on, less the flux with it off | 0: at a supersonic face the factor is exactly 1 | exact |
| `roe_low_mach_face_mach` | `Ma_Roe` of the pair built to sit at 1e-3 | 1e-3 | 5e-3 relative |
| `roe_low_mach_scaled_row_{m,p,e}` | Its three flux rows with the key on | the analytic flux with the velocity jump multiplied by that Mach number | 1e-12 relative |
| `roe_pressure_jump_part_unchanged_row_{m,p,e}` | The flux at factor 0 implied by the two calls, the flux affine in the factor | the analytic flux at factor 0: the pressure-jump and entropy parts of the dissipation did not move | 1e-9 relative |
| `roe_alternating_face_velocity` | `v_Roe` of WL=(1,-1e-4,1), WR=(1,+1e-4,1) | 0 exactly: equal densities put the Roe weight at 1/2 | exact |
| `roe_alternating_central_flux_row_{m,p,e}` | The corrected flux of that pair | the central flux `(F_L + F_R)/2` alone: with drho = dp = 0 and factor 0 no dissipation is left | 1e-15 absolute |
| `roe_alternating_mass_flux_uncorrected`, `..._corrected` | Its mass row, with the key off and on | 0, with and without the factor | 1e-18 absolute |
| `roe_alternating_momentum_damping_removed` | The momentum row the factor removes | `rho_avg a_avg dvel / 2`, the damping of the 2 dr alternating velocity mode | 1e-10 relative |

MEASURED 2026-09-18: 24 assertions, all pass with the key at its default.
The last group is the physical content of the option: where the two
velocities cancel in the Roe average the factor is 0, so the velocity-jump
dissipation is removed entirely and the alternating mode is left undamped.
The driver has no RED-at-entry statement: the branch and the test came in
together with item L25.

### 21. Convergence order of the face reconstruction (`weno3_reconstruction_order.f90`)

Origin: item WENO3-ORDER of `docs/Update_EXHALE_stage2.md` (2026-09-07),
measurement and fix in one item.
Production routines: `Reconstruct`, the vector reconstruction and
`weno3_geometry_coefficients` (`src/modules/states/Reconstruction.f90`).

The finite-volume update divides the difference of the two face fluxes by the
shell volume, so the state a cell carries is the volume average and the state
a face wants is a point value. Writing `V = r^3/3`, the volume average of a
cell is the plain average over an interval in V and a face in r is a face in
V, so the reconstruction is the Cartesian one in that frame. The driver
builds the volume averages by 8-point Gauss-Legendre quadrature (exact for
degree 15 in r, so the input carries no error of its own at the level
measured) and reports L_inf and L_1 of the face error with the rate between
successive resolutions.

The ESWENO3 smoothness factors of Yamaleev and Carpenter (2009), J. Comput.
Phys. 228, 3025, tend to 1 as O(h^2) on smooth data (their eq. 24), so in a
smooth region the scheme is its linear one and the order is the order of the
linear weights. Two of the four coefficients the code applied were the other
cell's volume share, `C1(j) = c/(b+c)` where the candidate needs `b/(b+c)`,
and likewise `C2(j-1)` where the candidate needs `b/(a+b)`. An O(h) error in
a coefficient multiplying an O(h) jump is an O(h^2) error in the face value,
so the scheme was second order. The swap was inherited from ATES v2.0.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `weno3_local_copy_matches_production` | The driver's copy of the formula, taking its coefficients as arguments, on the production coefficients | the production routine | 1e-14 absolute |
| `weno3_vector_matches_scalar_faces` | The vector reconstruction | the scalar one, face by face | 1e-13 absolute |
| `weno3_optimal_weights_exact_on_volume_quadratic` | The face value of the optimal-coefficient scheme on a quadratic in V | the exact face value: a third-order scheme reproduces a parabola in the volume coordinate exactly | 1e-12 absolute |
| `weno3_optimal_weights_third_order_{uniform,mixed_grid,gaussian}` | Convergence rate with the derived ideal coefficients, on a grid uniform in r, on the production stretched grid, and on a Gaussian profile | at least 2.8 | shortfall 0, exact |
| `weno3_corrected_stencil_coefficients_{uniform,mixed}` | The same with only the two stencil coefficients corrected and the production weight ratios kept | at least 2.8: the two coefficients are the whole fix | shortfall 0, exact |
| `weno3_production_third_order_on_uniform_volume` | The production scheme on a grid uniform in the VOLUME coordinate, where `C1 = C2` and the swap is invisible | at least 2.8: the control that identifies the cause rather than exhibiting a cure | shortfall 0, exact |
| `weno3_production_third_order`, `..._mixed` | The production scheme on the grid uniform in r and on the stretched one | at least 2.8 | shortfall 0, exact |

Each rate assertion is one-sided and is written as the SHORTFALL
`max(0, floor - rate)` against zero at zero tolerance, so a NaN rate fails.
The driver also prints, as diagnostics, the full L_inf and L_1 tables of six
schemes and a scan of the ESWENO floor `eps = (kappa dr_j)^2` over
`kappa = 1e-2` to 1e1, which is the measurement showing that the floor is not
what sets the order.

**RED before the fix of 2026-09-07** on `weno3_production_third_order` and
`..._mixed`: rate 1.997 on the grid uniform in r and 1.988 on the production
grid, against 3.000 and 2.991 with the volume shares corrected. MEASURED
2026-09-18: 11 assertions, all pass.

### 22. Positivity scaling of the reconstructed face states (`positivity_limiter_scaling.f90`)

Origin: item REF-ZHANGSHU of `docs/Update_EXHALE_stage2.md`, with
`docs/p49_front_energy_mode.md` section 6.1; the published version of the
paper was read.
Production routines: `positivity_limited_faces`, `positivity_scaling` and
`positive_variable_scaling` (`src/modules/states/Reconstruction_step.f90`).

The limiter cites Zhang and Shu (2010), J. Comput. Phys. 229, 8918, section
2.2, their eqs. (2.4), (2.5), (2.10), (2.11) and the quadratic (2.12). Their
limiter acts on the CONSERVED vector in two steps with an absolute floor;
this code scales the PRIMITIVE vector about the cell average in ONE step,
with `theta` the minimum over density and pressure and the relative floor
`eps = epsilon(1) q_avg`. The three consequences are what is asserted. The
reference state is a hypersonic layer at Mach 60, the state the limiter fires
in: the thermal pressure is a part in 3e3 of the total energy density, so a
reconstruction across a jump tips it negative.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `admissible_reconstruction_untouched`, `admissible_reconstruction_count` | An admissible reconstruction passed through the limiter, and the count of limited faces | unchanged to the bit, and zero | exact |
| `primitive_scaling_pressure_positive`, `..._density_positive`, `primitive_scaling_internal_energy` | The limited face of a reconstruction whose pressure is negative | strictly positive: pressure is itself a scaled variable, so `p/(gamma-1)` has the sign the scaling put there and the failure mode the quadratic of their eq. (2.12) exists to avoid has no path here | strict |
| `limited_face_count` | The count of limited faces of that state | 1 | exact |
| `theta_is_pressure_limited` | The `theta` the routine chose | `(p_a - eps)/(p_a - p_rec)`, pressure being the binding variable | 1e-14 relative |
| `limited_state_is_convex_combination` | The limited density | `rho_a + theta (rho_rec - rho_a)` | 1e-14 relative |
| `conserved_path_root_is_a_root` | The pressure at the root of their eq. (2.12), found by the driver | the comparison floor `eps = 1e-3 p_avg` | 1e-8 relative |
| `primitive_theta_at_most_conserved` | The code's `theta` | at most their `t_eps` at the same floor: along the conserved path `p(t)` is concave, so at equal eps this limiter is the more restrictive of the two and never the less safe one | at most |
| `bare_update_never_negative`, `bare_update_reaches_zero` | The sweep of the unclamped scaled value | never negative over the sweep, and exactly zero at least once: `eps` is one unit in the last place of the cell average, so the scaled value is a cancellation of the same size as the floor | at least |
| `clamped_scaling_delivers_the_floor` | The same sweep with the clamp | no state at zero | exact |
| `production_limiter_never_delivers_zero` | The production routine driven over the same reconstructed values | no face density or pressure at or below zero | exact |
| `endpoint_mean_departs_from_average` | The mean of a cell's two limited endpoints, one `theta` per face state | at least 1e-2 away from the cell average: conservativity, the third property their section 2.2 requires, is lost by a per-face `theta` | at least |
| `cell_theta_restores_average` | The same with their single `theta` per cell | at most 1e-14 from the average | at most |
| `nan_face_scales_to_the_average` | `positive_variable_scaling` of a NaN face value | 0, the cell average | exact |
| `ulp_floor_lost_in_conserved_round_trip` | The limited face mapped to conserved variables and back | a pressure differing from the limited one by at least the limited pressure itself: at Mach 60 a floor of one ulp of p is nineteen decades below E, so their admissible set of eq. (2.6) cannot be tested at it | at least |

What of their theorem carries over is stated in the driver header: positivity
of the face values is their Lemma 2.5 and needs no quadrature; their Theorem
2.1 does not carry over, its hypotheses being the N-point Legendre
Gauss-Lobatto set of their eq. (1.7) with `2N - 3 >= k` and the CFL bound of
their eq. (2.1), neither of which a finite-volume reconstruction to two face
values supplies for `k = 2`.

MEASURED 2026-09-18: 18 assertions, all pass. The clamp of the scaled density
and pressure to the floor is the change of 2026-09-07; before it the bare
update landed on exactly zero for about a part in 200 of the reconstructions
that fire the limiter, and a zero face density divides in `v = m/rho` and in
the sound speed of `Num_Fluxes.f90`.

### 23. The Ly-alpha beam a cell is pumped by (`lya_beam_cell_mean.f90`)

Origin: item B4-6b of `docs/Update_EXHALE_stage2.md`.
Production routines, all of `src/modules/radiation/lya_rt.f90`:
`lya_line_center_optical_depth`, `lya_stellar_beam_transmission`,
`lya_stellar_beam_transmission_cell_mean`,
`lya_photosphere_attenuation_cell_mean` and `jlya_escape_prob`.

Ly-alpha is resonantly SCATTERED, so the stellar beam is not attenuated by
`exp(-tau)`: the broad stellar line penetrates through its own wings and the
fraction reaching line-centre depth `tau` is `erfc(c sqrt(tau))` with
`c = sqrt(a/sqrt(pi))/(sqrt(2) Xs)`, the one-flight penetration condition
recorded at its definition site. The transmission is therefore not
exponential in the column and its cell mean has no closed form. The rate a
cell undergoes is what its atoms undergo averaged over the cell, and the sum
over the column of that mean times the cell depth is exactly the depth
integral of the beam, which is the number of scatterings the column performs.
The size of the correction is set by the change of the erfc ARGUMENT across
the cell, not by the cell depth: a cell of depth 10 deep in a column of depth
1e7 changes it by 1e-6, while the base cell of an exponential layer changes
it by order unity.

Configuration: 40 cells of width 0.05 R_p on a planet of 1e10 cm, an
isothermal 8000 K exponential H I layer of scale height 0.05 R_p whose base
line-centre depth is 1e7, no electrons and no protons so that the internal
source vanishes and the field is the stellar beam alone, and the trapping
boost set to 1 so that the field is the beam mean times a constant.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `lya_cell_depth_is_own_absorber` | Worst relative departure of `dtau(j)` from `n_HI C_lya dr / Dnu_D` | 0: a cell's line-centre depth is its own absorber and nothing else | 1e-12 absolute |
| `lya_top_cell_face_depth_is_zero` | `tau_out` of the outermost cell | 0 | 1e-99 absolute |
| `lya_top_cell_carries_its_own_depth` | `tau` of that cell | its own `dtau` | 1e-12 relative |
| `lya_faces_of_neighbouring_cells_agree` | Worst departure of `tau_out(j)` from `tau(j+1)` | 0: one face, one depth | 1e-12 absolute |
| `lya_pumping_field_is_cell_mean` | Worst departure of the field `jlya_escape_prob` hands the H(n=2) pumping, over the free-streaming value, from the cell mean of the beam over that cell's own depth | an independent fine quadrature (composite Simpson, 8192 intervals in `u = sqrt(tau)`) | 1e-6 relative |
| `lya_face_value_error_exceeds_half` | The largest relative error of the star-ward face value against that mean, less 0.5 | strictly positive, reference-only | strict |
| `lya_erfc_argument_step_of_base_cell` | The largest step of the erfc argument across a cell | strictly positive, reference-only | strict |
| `lya_column_scattering_budget_cell_mean` | `sum_j <T_s>_j dtau_j` over the column | the beam's own depth integral, from the same fine reference | 1e-6 relative |
| `lya_column_scattering_budget_face_error` | The face rule's error on that budget, less 1 per cent | strictly positive, reference-only | strict |
| `lya_beam_mean_dtau{0.30,1.00,3.00,10.00}_{top,photosphere}` | `lya_stellar_beam_transmission_cell_mean` at each of four cell depths, at the top of the column and at the depth where the beam is half absorbed | the fine reference | 1e-6 relative |
| `lya_beam_mean_square_dtau{...}_{top,photosphere}` | The mean of its square at the same eight places | the same reference for the second moment | 1e-6 relative |
| `lya_beam_face_error_dtau{...}_{top,photosphere}` | The face value's relative error at each | strictly positive, reference-only | strict |
| `lya_photosphere_factor_is_cell_mean` | `lya_photosphere_attenuation_cell_mean`, the mean of the parameterized `1/(1+tau)` attenuation the `jlya_mode = 0` field uses | its closed form `ln[(1+tau_in)/(1+tau_out)]/dtau` | 1e-12 absolute |
| `lya_photosphere_factor_column_identity` | The sum of that mean times the cell depth over the column | `ln(1 + tau_total)`, which it telescopes to exactly | 1e-12 relative |
| `lya_photosphere_face_value_error` | The inner-face value's error against that mean | strictly positive, reference-only | strict |
| `lya_photosphere_factor_thin_limit` | The same mean at zero depth | 1, the unattenuated field | 1e-12 relative |

The depth identities and the closed-form `1/(1+tau)` mean are exact and are
tested at round-off; the beam mean is a three-point Gauss-Legendre quadrature
and is tested against the independent fine reference, never against itself.

**RED before item B4-6b.** The stellar Ly-alpha depth was a face value, so
the pumping field of a cell was the beam at that cell's star-ward face.
MEASURED 2026-09-18: 37 assertions, all pass.

### 24. Escape probability of the fine-structure lines (`fine_structure_escape_probability.f90`)

Origin: item LYA-BETA-METALS of `docs/Update_EXHALE_stage2.md`, the
measurement taken by the worker and the fix applied by the advisor.
Production routine: `line_escape_probability_one_face`
(`src/modules/radiation/Cool_coeff.f90`), which
`fine_structure_line_transfer` evaluates on the outward and inward
line-centre columns of each cell and hands to the [C I], [C II], [N II] and
[O I] statistical equilibrium as `A_ul -> beta A_ul`.

These are forbidden lines of a static, thermally broadened layer. The Voigt
parameter is of order 1e-8, so a photon escapes out of the Doppler CORE and
the published solution is the complete-redistribution Doppler-slab one, not
the damping-wing solution the Ly-alpha closure uses. de Jong, Boland and
Dalgarno (1980), A&A 91, 68, appendix B give as their eq. (B-7)
`beta = [1 - exp(-2.34 tau)]/(4.68 tau)` for `tau < 7` and
`beta = 1/(4 tau [ln(tau/sqrt(pi))]^(1/2))` for `tau >= 7`, "accurate to
within 10% for small and intermediate tau and exact at very large tau". Since
`phi` is normalized, `phi(0) = 1/sqrt(pi)` and the argument of (B-7) is the
frequency-integrated depth `tau_dJ = sqrt(pi) tau_centre`. The independent
check of that convention is Hollenbach and McKee (1979), ApJS 41, 555,
eq. (5.10).

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `fs_escape_is_de_jong_b7_thin_branch` | Worst relative departure of the routine at `tau` = 1e-4 to 6 | (B-7)'s thin branch, the two constants written out in the driver | at most 1e-14 |
| `fs_escape_is_de_jong_b7_thick_branch` | The same at `tau` = 7 to 1e7 | (B-7)'s thick branch | at most 1e-12 |
| `fs_escape_transparent_layer_is_one_half` | `beta(0)` | 1/2: half the photons leave through the near face, which is why the cell closure adds the two faces | 1e-14 relative |
| `fs_escape_branch_point` | The depth at which the routine switches branch | `sqrt(pi) exp(2.34^2/4) = 6.967558963071725`, where the two expressions of (B-7) meet, (B-7) quoting it as 7 | 1e-12 relative |
| `fs_escape_branch_is_continuous` | The two values a part in 1e12 on either side of it | equal: the switch loses only the `exp(-2.34 tau_c)` term the thick branch drops, 8e-8 relative | 1e-6 relative |
| `de_jong_b7_agrees_with_hollenbach_mckee_5_10` | Worst departure of `2 beta_dJ(sqrt(pi) tau)` from `eps_HM(tau)` at `tau` = 3 to 1e6 | at most 5e-2, the published accuracy of two independent fits | at most |
| `de_jong_b7_is_hollenbach_mckee_at_large_tau` | The same at `tau = 1e5` | `eps_HM(1e5)`: asymptotically they are the same expression | 1e-4 relative |
| `fs_escape_line_centre_argument_departure_case1..8` | `beta(tau_centre)/beta(sqrt(pi) tau_centre)` at the eight comparison depths | at least 1.75 | at least |
| `fs_escape_departure_at_shipped_depths` | The same ratio at `tau_centre = 0.14` | at most 1.13 | at most |

**The production argument is the (B-7) one.** `fine_structure_line_transfer`
multiplies each cell's line-centre depth by `sqrt(pi)` before any argument is
formed (`Cool_coeff.f90`, the `dl = dr_j R0 sqrt(pi)` line of that routine,
since 2026-09-07; ISSUES_20260909 records the correction under REF-LODI). The
last two rows above measure what the departure WOULD be if the line-centre
column were passed uncorrected, as a property of the function
`beta(tau)/beta(sqrt(pi) tau)`: a factor rising from 1 in the thin limit to
1.77 asymptotically with a maximum of 2.11 at the branch point, and 12 per
cent at the largest line-centre depth any shipped case reaches (0.14, [O I]
63 um at the base of the hot-Uranus gate). They are the record of the size of
the defect that correction removed, not an open one (the driver's header said
"open defect, not fixed" until 2026-09-18, which was stale). MEASURED
2026-09-18: 16 assertions, all pass.

### 25. The XUV beam a cell's rate is evaluated with (`xuv_cell_mean_attenuation.f90`)

Origin: item B4-6 of `docs/Update_EXHALE_stage2.md`
(`utilities.f90` `cell_mean_attenuation`, `util_ion_eq.f90`).
Production routine: `PH_heat_H` (`src/modules/radiation/util_ion_eq.f90`),
through `calc_column_dens`.

`calc_column_dens` accumulates from the top down, so `N1(j)` holds the WHOLE
of cell j and the depth built from it is the depth at that cell's INNER face.
The rate that belongs to the cell is the mean over the cell of the field,
which for a uniform absorber density, the rectangle rule the column
integration itself uses, is `<exp(-tau)> = exp(-tau_out)(1 - exp(-dtau))/dtau`.

Configuration: a twelve-cell uniform slab (N = 8 plus two ghost cells at each
end) of H I alone, no helium, no metals, no secondary ionization,
`opa_pf = 1`, `a_tau = 0`, and a MONOCHROMATIC 20 eV spectrum, so that a
cell's optical depth is one number rather than a spectrum of them and
"dtau = 1" means dtau = 1. `n_HI` is set from the target depth; nothing else
changes between the four cases `dtau` = 0.01, 0.3, 1 and 3.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `rate_is_cell_mean[dtau=...]` | `P_HI(j)` of every cell, at each of the four depths | an INDEPENDENT quadrature of the same integral (composite three-point Gauss-Legendre in the depth across the cell, 128 segments, truncation below 1e-16 relative), not the closed form the code uses | 1e-12 relative |
| `column_ionizations_equal_beam_loss[dtau=...]` | `sum_j P_HI(j) n_HI dr` over the column | `Phi_in [1 - exp(-tau_total)]`, the photons the beam loses between its two ends: with one absorber every photon lost is one H I photoionization, and the cell mean telescopes over the column exactly | 1e-12 relative |
| `inner_face_rule_ratio[dtau=...]` | The inner-face rate divided by the cell mean, on the reference arithmetic alone | its closed form `dtau exp(-dtau)/(1 - exp(-dtau))`, which is 0.582 at `dtau = 1`, a rate low by 41.8 per cent | 1e-12 relative |

The photon-number tolerance is 1e-12 and not round-off because the columns
are accumulated cell by cell, so a last-digit perturbation of the accumulated
depth is amplified by `tau_total`, which reaches 36 in the `dtau = 3` case,
and twelve such terms are summed. The third row is reference-only and keeps
the test discriminating: it fails if the two discretizations ever coincide.
MEASURED 2026-09-18: 56 assertions, all pass.

### 26. The order the field and its opacity are brought into agreement (`photoionization_field_substitution.f90`)

Origin: item COST3 of `docs/Update_EXHALE_stage2.md`, a measurement item
(`util_ion_eq.f90`, `ionization_equilibrium.f90`).
Production routines: `advance_starward_columns` and `calc_column_dens`, and
the equilibrium sweep `xuv_block_H` / `xuv_block_HHe`
(`src/modules/radiation/ionization_equilibrium.f90`).

A cell's photoionization rate is set by the column of the cells OUTSIDE it
and by its own, and by nothing inside it. The sweep runs outside in, so it
can be given the columns of the composition it has ALREADY solved (a
substitution) or the columns of one whole composition held fixed for the
traversal (a lagged, Jacobi map). Both are iterated to the same statement,
composition in equilibrium with the field its own opacity makes.

Configuration: the slab of group 25, made thick enough to hold an ionization
front and resolved finely enough that one cell is a fraction of an optical
depth. H I alone, no helium, no metals, no secondary ionization,
`opa_pf = 1`, `a_tau = 0`, a monochromatic 20 eV spectrum, cells 2.5e7 cm
wide, case-B recombination 2.59e-13 cm^3 s^-1, and the density set from the
unattenuated rate so that the outermost cell sits at an ionization parameter
of 30. The chemistry is photoionization equilibrium at a fixed recombination
coefficient, `P_HI (1-x) n_H = alpha n_H^2 x^2`, solved in closed form. It is
not EXHALE's network; it is the simplest composition that responds to the
field the way the real one does, which is all the ordering statement needs.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `starward_column_matches_calc_column_dens` | `N1_face(j)` of every cell, from `advance_starward_columns` run over the whole grid in one block | `calc_column_dens` at `j+1`, BIT FOR BIT: same rectangle rule, same `opa_pf` weight, same order of accumulation | exact |
| `starward_column_of_the_outermost_cell_is_zero` | The same at the outermost cell | 0: nothing sits above it | exact |
| `starward_column_leaves_the_whole_column` | The column the block leaves behind | `calc_column_dens` at the innermost cell | exact |
| `the_two_orderings_have_one_fixed_point` | Largest difference of the two ionization profiles, each iterated to 1e-12 | 0: the reordering may not move the answer | 1e-9 absolute |
| `the_two_orderings_have_one_field` | Largest relative difference of the two fields at those fixed points, through the production routine | 0 | 1e-6 absolute |
| `substitution_passes_to_the_csm_tolerance` | Passes the substitution takes to reach 1e-8, the tolerance the coupled source step stops at | at most 12 | at most |
| `lagged_passes_to_the_csm_tolerance` | Passes the lagged map takes to the same tolerance | at least 18 | at least |

MEASURED on this slab (250 cells, neutral optical depth 90, front near cell
190): the lagged map takes 24 passes to reach 1e-8 and the substitution 8,
and 29 against 11 to reach round-off. What is left for the substitution to
iterate is each cell's OWN optical depth, which no ordering of the column can
remove; that is why the ratio is about three here and not the grid size. The
two bounds are set wide enough for the arithmetic to move and narrow enough
that a change putting the field back on the lagged column would fail the
first. MEASURED 2026-09-18: 259 assertions, all pass, of which 253 are the
first row, one per cell.

### 27. The field one cell is solved in (`photoionization_field_self_consistency.f90`)

Origin: item FIELD-SELF of `docs/Update_EXHALE_stage2.md`
(`ionization_equilibrium.f90`, the self-field loop, default one pass).
Production routines: `xuv_self_field_passes`
(`src/modules/radiation/ionization_equilibrium.f90`),
`photoionization_field_at_cell_H` / `_HHe`
(`src/modules/radiation/util_ion_eq.f90`) and `cell_mean_attenuation`
(`src/modules/functions/utilities.f90`).

A cell's attenuated field is `exp(-tau_out)` times the cell mean of its OWN
attenuation, and that own depth `dtau = sum_abs sigma_nu n_abs dr` is built
from the very densities the cell is being solved for. One field build and one
solve therefore return a composition whose own opacity is not the opacity the
field was built with.

Configuration: the slab of groups 25 and 26, with the cell made ONE optical
depth thick at the neutral density, which is where the cell mean departs from
the face value by tens of per cent and where the self-consistency defect is
largest. `n_H = 1e9` cm^-3, ionization parameter 1 at the unattenuated cell
so that the equilibrium fraction sits near a half where `dx/dP` is largest,
an outside column of 0.7 optical depths for assertion (2), and the same
closed-form photoionization equilibrium as group 26.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `unattenuated_rate_is_positive` | The rate of an empty cell with nothing above it | strictly positive | strict |
| `one_cell_field_is_the_cell_mean_of_its_own_column` | The rate of the same cell filled with neutral hydrogen, still nothing above it | the unattenuated rate times `(1 - exp(-dtau))/dtau`, the single mean of `cell_mean_attenuation` | 1e-14 relative |
| `field_is_one_cell_mean_over_both_depths` | The rate with a column outside as well | the unattenuated rate times `exp(-tau_out)(1 - exp(-dtau))/dtau`, the ONE expression the routine forms from the two depths together | 1e-14 relative |
| `a_second_attenuation_would_be_visible` | Relative distance from a face value multiplied by a cell mean twice | at least 1e-2: what the routine returns is not the attenuation composed twice | at least |
| `two_evaluations_at_one_state_return_one_rate`, `..._one_heating_rate`, `..._one_efficiency` | Two evaluations at one composition | equal bit for bit: an iterating caller is iterating a map and not a history | exact |
| `one_pass_leaves_a_self_consistency_defect` | The move when the field is re-formed from the composition one pass returned and the cell solved again | at least 1e-4, well below the measured value and well above round-off, so the row states "one pass is not self-consistent" and nothing finer | at least |
| `self_consistent_cell_defect_is_at_round_off` | The same defect after the iteration has run | at most 1e-12 | at most |
| `self_consistent_cell_passes_are_bounded` | The pass count it takes | at most 60 | at most |
| `self_field_mode_contracts` | The ratio of successive moves | at most 0.9 | at most |
| `one_pass_from_the_fixed_point_returns_it` | A cell handed the converged composition and solved ONCE | the same state: the two schemes differ in the path, not in the answer | 1e-12 absolute |

MEASURED on this cell, and printed as diagnostics: the defect after one pass
is of order a per cent of the ionized fraction itself. That defect is the
statement `xuv_self_field_passes` rests on, and the statement its default of
one leaves open. MEASURED 2026-09-18: 12 assertions, all pass.

### 28. One absorbed fraction, two beams (`absorbed_fraction_switch.f90`)

Origin: item HYG-absorbed of `docs/Update_EXHALE_stage2.md`
(`water_photolysis.f90`, and a comment of `utilities.f90`).
Production routines: `absorbed_fraction_per_unit_depth` and
`cell_mean_attenuation` (`src/modules/functions/utilities.f90`), the XUV cell
mean of `util_ion_eq.f90`, and `water_photolysis_rate` and
`hydroxyl_photolysis_rate`
(`src/modules/lower_atmosphere/water_photolysis.f90`).

The quantity is `(1 - exp(-d))/d`, the fraction of the photons entering a
cell that the cell absorbs per unit of its own optical depth, and the
statement is that the XUV beam and the FUV bands take it from ONE definition.
The grid is 1501 points logarithmically spaced over `d` = 1e-12 to 1e3, plus
the switch point itself and its two neighbours in the double-precision grid,
so that both branches and the crossing are covered.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `absorbed_fraction_error_vs_quad` | Largest relative error over the whole grid | a QUADRUPLE-PRECISION evaluation of `(1 - exp(-d))/d`; the closed form subtracts two numbers that approach each other as `d -> 0`, so its relative error is about `eps/d` and is 1e-8 just above the switch at `d = 1e-8`, while the series `1 - d/2 + d^2/6` truncates at 4e-25 there | 2e-8 absolute |
| `absorbed_fraction_error_above_dtau_1em6` | The same error restricted to `d >= 1e-6` | the same reference | 1e-9 absolute |
| `water_band_rate_is_xuv_cell_mean_bitwise` | Count of grid points at which `water_photolysis_rate(F,ib,0,d,1)` differs from `s_b N_b cell_mean_attenuation(0,d)` | 0: both sides reach the same function through the same use association, so a second copy of the expression with any other switch point or series fails somewhere in the 15 decades scanned | exact |
| `hydroxyl_band_rate_is_xuv_cell_mean_bitwise` | The same for OH | 0 | exact |
| `absorbed_fraction_at_zero_depth` | `absorbed_fraction_per_unit_depth(0)` | 1, the limit, returned exactly; a cell of no optical depth is what every empty cell of a run carries | exact |

The two bitwise rows and not the error bound are what pin the shared
definition: a switch point moved up far enough that the series truncation
shows would break neither bound. MEASURED 2026-09-18: 5 assertions, all pass.

### 29. The one formation-energy table (`species_formation_energy_table.f90`, `species_formation_energy_python_table.py`)

Origin: item B3b-ENTH of `docs/Update_EXHALE_stage2.md`, against
`docs/b1_target_system_20260906.md` T1.2, T1.3 and T1.9.
Production routines: `species_formation_energy` and
`molecular_reaction_energy_eV`
(`src/modules/lower_atmosphere/molecular_reaction_heat.f90`);
`oxygen_formation_energy_eV` and `photolysis_threshold_erg`
(`src/modules/lower_atmosphere/oxygen_rates.f90`);
`h2_dissociation_energy_eV` (`src/modules/lower_atmosphere/mol_rates.f90`);
the H2 level ladder of
`src/modules/lower_atmosphere/molecular_infrared_data.f90`; and the
composition metadata of `src/modules/init/species_table.f90`, against which
every reaction is balanced.

T1.2 fixes one reference state for the whole code, every element a neutral,
ground-state, free atom at rest and the free electron at zero, and T1.9
requires every reaction energy of every network to be a DIFFERENCE of that
one table rather than a number written beside the reaction. A table has
properties a list of reaction enthalpies does not: a forward and a reverse
channel are exact negatives, a closed cycle releases exactly zero, and two
routes between the same two states release the same energy. Each of them
fails if any energy is transcribed twice. The driver also carries its own
transcription of the 17-reaction stoichiometry, taken from the reaction
names, so the module's stoichiometry table is never the reference for itself.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `eps_zero_{HI,HeI,CI,OI,NaI,FeI,electron}` | The neutral ground-state entries and the free electron | 0 | exact |
| `eps_{HII,HeII}_is_IP_{H,He}`, `eps_HeIII_is_cumulative` | The H and He ion entries | `e_th_HI`, `e_th_HeI`, `e_th_HeI + e_th_HeII` | exact, exact, 1e-12 absolute |
| `eps_{CII,CIII,NaII,FeIII}_cumulative` | Four metal stages | the sum of the potentials below them, 11.26, 11.26 + 24.38, 5.139 and 7.902 + 16.199 eV, READ from `mion_ethr` | 1e-12 absolute |
| `eps_He23S_excitation`, `eps_H_n2_is_Lyman_alpha` | The two excited levels | `e_th_HeI - e_th_HeTR`, and `0.75 e_th_HI`, the Lyman-alpha energy | 1e-12 absolute |
| `eps_H2_is_minus_D0`, `eps_H2_uses_the_one_D0` | The H2 entry | `-36118.11 cm^-1` in eV (Huber and Herzberg 1979), and `-h2_dissociation_energy_eV()` | 1e-9 absolute, exact |
| `h2_ladder_zero_is_v0J0`, `h2_ladder_limit_is_D0` | The lowest and highest term energies of the H2 ladder | 0, and D0 above that zero: the equation of state holds the ladder and the reservoir holds the bond, with no overlap (T1.3) | exact, 5e-5 relative |
| `reaction_energy_<R>` (17 rows) | Each of the 17 network reaction energies | the difference of table entries built from the driver's own stoichiometry | 1e-12 absolute |
| `balance_<R>` (17 rows) | Nuclei and charge of each reaction | 0 | exact |
| `reverse_pair_R9_R10`, `reverse_pair_R12_R15` | A channel plus its reverse | 0 | exact |
| `same_change_R12_R14` | Two channels with the same reactants and products | equal | exact |
| `route_R8_R6_equals_R5` | R8 followed by R6 | R5 with an H2 spectator | 1e-12 absolute |
| `cycle_R8_header_value`, `cycle_R6_header_value`, `cycle_leaves_IH2_minus_D0` | The three energies the module header quotes for the cycle that makes this channel 81 per cent of the base heating | 1.70 eV, 9.25 eV, and `I(H2) - D0(H2)` | 0.02 eV, 0.02 eV, 1e-12 absolute |
| `eps_zero_{O,C,H}_shomate`, `shomate_H2_uses_the_one_D0` | The Shomate species | 0, and the one `eps(H2)`: not a second fit of the bond | exact |
| `eps_{OH,H2O,CO}_bond` | The three oxygen molecules | -4.39888, -9.51622 and -11.11569 eV, MEASURED from the F coefficients of the lowest Shomate intervals | 1e-4 absolute |
| `eps_O1D_excitation` | O(1D) | `e_excite_O1D_erg` in eV: an excited state carries its excitation energy only | 1e-12 absolute |
| `threshold_{H2O_OH_H,H2O_H2_O,H2O_O_H_H,OH_O_H}_is_difference` | The four photolysis thresholds | the differences of the same table entries | 1e-12 absolute |
| `H2O_two_step_equals_one_step`, `H2O_H2_branch_differs_by_D0_H2` | Taking H2O apart in one step or in two, and the H2-forming branch against the full atomization | 0, and exactly one H2 bond | 1e-12 absolute |
| `oxygen_O1_reverse_is_exact_negative`, `oxygen_O6_carries_O1D_excitation`, `oxygen_O2_then_O1_is_overall` | The collisional oxygen channels | 0: the ledger they will be read from is closed, though no consumer deposits them today | 1e-12 absolute |

The driver writes the table it asserted to
`species_formation_energy_table.dat` in `EXHALE_TEST_WORK`, and
`species_formation_energy_python_table.py` reads that file and compares it
entry by entry with the `EPS_EV` dictionary of
`src/utils/formation_energy_flux_diagnostic.py`, which carries its own copy
because it reads output files and never links the code. MEASURED 2026-09-18:
78 Fortran assertions and 37 Python ones, all pass.

### 30. One assembly of the heating rate (`heating_channel_closure.f90`, `heating_sum_uniqueness.py`)

Origin: item HEATBRK of `docs/Update_EXHALE_stage2.md` (2026-09-06).
Production routine: `heating_of_composition`
(`src/modules/radiation/util_ion_eq.f90`), the one place the heating of a
cell is put together: the ionization sweep takes the total it returns for the
energy equation, `output/Heating_breakdown.txt` writes the channel array it
fills, and the advection-corrected post-process calls it again on its own
composition.

The property that makes one assembly worth having is that the total is the
running sum of the channels, so a deposit cannot reach the energy equation
and be absent from the breakdown. The second assertion is what makes the
first mean anything: a closure test on a cell whose molecular, oxygen and
metastable deposits are all zero closes trivially, and it was exactly two
silent zeros in a second copy of this sum, the associative He(2^3S) branch
and the collisional oxygen channels, that the breakdown file was missing.

The cell is a molecular base: 1400 K, 1e13 cm^-3 of H2 over a partly ionized
H/He gas with the He 2^3S metastable, the oxygen carriers OH, H2O and CO,
trace metals, an attenuated XUV field, a Lyman-Werner band flux and FUV
photolysis rates. The numbers are order-of-magnitude values of the
hot-Uranus gate, not a solution of anything: nothing here is asserted against
a physical reference, only the closure of the sum and the fact that every
channel fires.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `heat_total_is_channel_sum` | Worst relative departure, over the cells, of the sum of the channel array from the returned total | 0 | 1e-14 absolute |
| `every_heat_channel_exercised` | Count of channels with no deposit on this cell | 0; the driver names any channel that is silent | exact |

`heating_sum_uniqueness.py` is the scan that keeps the assembly single. Every
`.f90` file under `src/modules` is read with comments and character strings
removed and searched for the production entities that FORM a heating deposit;
each may be named in the module that defines it and in
`src/modules/radiation/util_ion_eq.f90`, and nowhere else. Reference: zero
files outside those two per name, tolerance zero. The scan cannot catch an
assembly written out of literal arithmetic; what it catches is the way all
three copies were actually written, by calling these routines. MEASURED
2026-09-18: 2 Fortran assertions and 1 Python one, all pass.

### 31. CO destruction (`co_destruction.f90`)

Origin: item B3b-CO of `docs/Update_EXHALE_stage2.md`, with
`docs/b3b_co_destruction_design_20260906.md` and
`docs/co_destruction_rates_literature_20260906.md`.
Production routines: `co_self_shielding` and `co_self_shielding_texc5`
(`src/modules/lower_atmosphere/co_self_shielding_table.f90`);
`co_photodissociation_rate`, `co_photodissociation_rate_cell_mean` and
`heat_per_co_dissociation`
(`src/modules/lower_atmosphere/co_photodissociation.f90`);
`rk_D1_Hep_CO`, `rk_CO_radiative_association` and `photolysis_threshold_erg`
(`src/modules/lower_atmosphere/oxygen_rates.f90`);
`oxygen_reaction_energy_eV` and `species_formation_energy`
(`src/modules/lower_atmosphere/molecular_reaction_heat.f90`).

References, all READ from the publisher PDF
`references/Visser_2009A&A_503_323.pdf` unless said otherwise: Visser, van
Dishoeck and Black (2009), A&A 503, 323, Table 6 (the shipped shielding
function, b(CO) = 0.3 km/s, T_ex(CO) = 50 K, T_ex(H2) = 501.5 K) and Table 5
(the same block at T_ex(CO) = 5 K, the paper's reference set), their Table 1
(the 37 predissociating bands with lambda_0, f_v0 and eta), and UMIST RATE22
entries 4068 and 8597 (Millar, Walsh, Van de Sande and Markwick 2024, A&A
682, A109).

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `co_shield_table6_grid_points` | Worst departure of `co_self_shielding` at every grid point of Table 6 | the transcribed table | 1e-12 absolute |
| `co_shield_texc_spread_max`, `..._min` | Largest and smallest ratio of Table 5 to Table 6 over the grid | 22.0701 and 0.548675, a REGRESSION on the transcription of the two blocks and not a physical tolerance; it is the size of the excitation-temperature uncertainty and is larger than the factor two Visser's section 5.2 quotes for a whole cloud model, that sentence being about a model rate and this the worst single grid point | 1e-4 relative |
| `co_shield_interpolation_base_cell` | `co_self_shielding_texc5` at log N(CO) = 18.96, log N(H2) = 21.80 | 2.6e-5, the log-bilinear interpolation of Table 5 computed by hand in the literature document on the base cell of the `oxygen_chemistry` case | 2e-2 relative |
| `co_shield_clamp_above_co_axis`, `..._h2_axis` | The table above the top of either axis | its corner entry, Visser's section 5.1 stating why that is the right behavior and not a clamp of convenience | 1e-12 relative |
| `co_shield_zero_column_unshielded` | A column of zero | 1 | 1e-12 relative |
| `co_shield_tex_limit` | The excitation-temperature limit the run's domain record counts cells against | 512 K, Visser's section 4.4 | 1e-12 relative |
| `co_photon_energy_37_band_mean` | `e_co_photon_erg` | the f-weighted mean of `hc/lambda_0` over the 37 rows of Table 1, transcribed in the driver and recomputed, so the constant cannot drift from the data it came from | 1e-6 relative |
| `co_thin_rate_is_sigma_times_photon_flux` | The rate at Theta = 1 with no continuum | `sigma_CO F/<hv>` | 1e-12 relative |
| `co_band_photon_energy_matches_lw` | The `<hv>` that rate divides by | the Lyman-Werner band's own constant, not a copy that can drift from it | 1e-12 relative |
| `co_cross_section_against_visser_k0` | 1.5496e-17 cm^2, the shipped 912-1201 A cross section restated over 912-1110 A | `k0/F_Draine = 2.590e-10/1.232e7`, Visser's Table 6 header over the Draine 912-1110 A photon flux (Draine and Bertoldi 1996, Table 1) | within a factor 1.5, Visser's own 20 per cent accuracy loosened by the difference in spectral shape stated in `co_photodissociation.f90` section 2 |
| `co_hep_rate_{300,1500,3000}K` | `rk_D1_Hep_CO` at three temperatures | 1.60e-9 cm^3 s^-1: RATE22 4068 has beta = gamma = 0, so an accidental Arrhenius factor would break it | 1e-15 relative |
| `co_hep_reaction_energy_from_table` | `oxygen_reaction_energy_eV(D1)` | the same difference formed here from `species_formation_energy` with the runtime species keys | 1e-14 relative |
| `co_hep_reaction_energy_value` | The same | `e_th_HeI - IP(C I) + eps(CO)`, the helium and carbon ionization potentials less the C=O bond | 1e-12 relative |
| `co_hep_reaction_is_exothermic` | The same | strictly positive | strict |
| `co_bond_energy_from_table` | D0(CO) of the photolysis channel | `-eps(CO)`: the two channels of one molecule cannot state two bond energies | 1e-14 relative |
| `co_photodissociation_deposit`, `..._positive` | `heat_per_co_dissociation` | `<hv> - D0(CO)`, the rest of the photon, and strictly positive | 1e-6 relative, strict |
| `co_cell_mean_against_fine_quadrature` | The cell mean over a cell spanning three decades of CO column and two of H2 column with a continuum depth across it | a 4000-segment midpoint sum of the same integrand | 1e-3 relative |
| `co_cell_mean_head_segment` | The cell mean of a cell whose star-ward column is zero | the same fine sum over the head segment | 1e-3 relative |
| `co_cell_mean_uniform_cell_is_point_rate` | The cell mean of a cell with no column gradient | the point rate | 1e-12 relative |
| `co_rate_zero_without_band_flux` | The cell mean at zero band flux | 0 | exact |
| `co_radiative_association_rate_1440K` | `rk_CO_radiative_association(1440 K)` | `4.69e-19 (T/300)^1.52 exp(50.5/T)`, RATE22 8597 recomputed from its own four numbers | 1e-12 relative |
| `co_radiative_association_positive` | The same | strictly positive; it enters no row and only the record | strict |

MEASURED 2026-09-18: 27 assertions, all pass. The driver was written RED
against a mutated table entry and a mutated band row, which is how the two
transcriptions were checked.

### 32. He+ + CO as a species sink (`co_helium_ion_sink.f90`)

Origin: item B3b-CO2 of `docs/Update_EXHALE_stage2.md`, "the He+ + CO channel
now consumes the helium ion".
Production routines: `mol_heh_rows`, `set_mol_coeffs` and
`set_mol_turnover_rates`
(`src/modules/nonlinear_system_solver/System_HeH_mol.f90`); `rk_D1_Hep_CO`
(`src/modules/lower_atmosphere/oxygen_rates.f90`);
`oxygen_reaction_energy_eV`
(`src/modules/lower_atmosphere/molecular_reaction_heat.f90`).

UMIST RATE22 entry 4068, He+ + CO -> O + C+ + He, alpha = 1.60e-9, beta = 0,
gamma = 0, method M, accuracy A: one helium ion, one CO molecule and 2.2117
eV of product kinetic energy per event. The helium leaves neutral, so the
free neutral He that closes the helium budget receives it and no other helium
row moves. The cell is the shielded molecular base of the `oxygen_chemistry`
case in order of magnitude (T = 1440 K, n_H = 1e13, n_He = 1e12,
n(CO) = 2.4e10 cm^-3), so the numbers the assertions carry are the ones the
code meets there.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `hep_row_co_sink` | The He+ row of `mol_heh_rows` evaluated at n(CO) = 0 and at n(CO) > 0 with everything else held, differenced | `-k_D1 n(CO) n(He+)` | 1e-12 relative |
| `hep_row_co_sink_is_a_loss` | `sink + abs(sink)` | 0: a sink written with the wrong sign would pass a magnitude test and fail this one | exact |
| `hep_row_derivative_co_term` | The gain in `d(row 2)/d n(He+)` by central differences when the CO is switched on, step 1e-4 of n(He+) | `-k_D1 n(CO)`, the entry the numerical Jacobian of `hybrd1` builds | 1e-8 relative |
| `co_moves_only_the_helium_ion_row` | Count of the other six rows that move | 0: the reaction touches helium and carbon alone, and carbon is not a row of this system | exact |
| `heat_and_species_event_rate` | The event rate the heating assembly deposits channel 18 for | the event rate the row performs, one ion per event | 1e-14 relative |
| `reaction_energy_D1_exothermic` | `oxygen_reaction_energy_eV(D1)` | strictly positive, and the SAME number the row is paired with (group 31 asserts where it comes from) | strict |
| `hep_turnover_scale_co_term` | The rise in `1/mol_inv_turnover(2)`, the bound the He+ row is divided by before `hybrd1` reads it | `k_D1 n(CO) n_He`: the bound puts the whole helium element on the ion, as it does for every other term of that row | 1e-10 relative |
| `hep_row_unchanged_without_oxychem` | The row with the oxygen chemistry off and a stale `n_co` left in the cell state | its value at n(CO) = 0 | exact |

MEASURED 2026-09-18: 8 assertions, all pass.

### 33. The two product channels of H3+ recombination (`h3p_recombination_branching.f90`)

Origin: item REF-CHEM of `docs/Update_EXHALE_stage2.md` (comments and
documents only; `mol_rates.f90`).
Production routines: `rk_R6_H3p_dr_H2` and `rk_R7_H3p_dr_3H`
(`src/modules/lower_atmosphere/mol_rates.f90`). Together they are the H3+
sink of the molecular base, and so they set how much H3+ survives to radiate
through `h3p_cooling.f90`.

Larsson, McCall and Orel (2008), Chem. Phys. Lett. 462, 145, p. 149, give the
thermal rate constant of the Kokoouline and Greene calculation as
`alpha(300 K) = (7.2 +- 1.1)e-8` cm^3 s^-1, "in good agreement with the new
storage ring results", and the three-body branching ratio 0.70 +- 0.07, "in
very good agreement with the CRYRING storage ring results". Neither constant
in the module is printed in that paper: each is the total split by the
branching and carried in temperature by the index 0.65 of the earlier
storage-ring fits, which Yelle (2004), Icarus 170, 167, Table 1 prints as
R16a and R16b.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `h3p_dissociative_recombination_total_T{300,1000,3000}K` | The sum of the two channels | `7.2e-8 (300/T)^0.65` | 1e-12 relative |
| `h3p_recombination_to_h2_and_h_T{300,1000,3000}K` | The H2 + H channel | `0.30` of that total | 1e-12 relative |
| `h3p_recombination_to_three_atoms_T{300,1000,3000}K` | The 3H channel | `0.70` of that total | 1e-12 relative |
| `h3p_three_body_branching_fraction_T{300,1000,3000}K` | The 3H share | 0.70 at every temperature | 1e-12 relative |

The three temperatures are the base, the H3+ layer above it and the top of
the molecular region, the whole range over which the extrapolation is used.
The tolerance is round-off because these are transcriptions. The `note:`
lines carry the ratio to the older Sundstrom and Datz pair that Yelle's table
and Frelikh and Murray-Clay (2026) Table 1 both use; Larsson's own reading of
that pair (p. 149) is that "the early results obtained at CRYRING [23,24] and
ASTRID [27], which gave results just above or at 1e-7 cm^3 s^-1, were
slightly too high because of rotational excitations", so the ratio is a
superseded value over the one that supersedes it and not a disagreement
between equals. MEASURED 2026-09-18: 12 assertions, all pass.

### 34. H2O photolysis yields on the Ly-alpha band (`water_photolysis_lyman_alpha_yields.f90`)

Origin: item REF-CHEM of `docs/Update_EXHALE_stage2.md`.
Production data: `qy_H2O_OH_H`, `qy_H2O_H2_O1D` and `qy_H2O_O_H_H` on the
`fuv_band` table of `src/modules/lower_atmosphere/oxygen_rates.f90`. The
three-body branch is open in no other band, and it is the only photolytic
source of ground-state O in the option, so a shift of the triplet moves where
O and OH come from at the base.

Slanger and Black (1982), J. Chem. Phys. 77, 2432, measured H(2S) and
O(1D, 3P) yields from H2O photolysis at 1216 A. Their eqs. (1) to (4) are the
four channels, and p. 2435 reads "We may thus set yields of 78% for processes
1 + 2, 10% for process 3, and 12% for process 4". The 0.78 is therefore
Phi1 + Phi2 and not a single ground-state channel, and the three are not
independent: the 10% is adopted from Stief et al., the 12% is derived against
it, and the 78% is what normalization leaves.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `lyman_alpha_band_count` | That a band of the table contains 1215.67 A | it does: the yields are a single-wavelength measurement and mean nothing on a band that does not contain the line | 1e-12 relative |
| `bands_with_three_body_branch_open` | The number of bands with the three-body branch open | 1 | 1e-12 relative |
| `h2o_photolysis_yield_oh_plus_h_lyman_alpha` | `qy_H2O_OH_H` of that band | 0.78 | 1e-12 relative |
| `h2o_photolysis_yield_h2_plus_o1d_lyman_alpha` | `qy_H2O_H2_O1D` | 0.10 | 1e-12 relative |
| `h2o_photolysis_yield_o_plus_2h_lyman_alpha` | `qy_H2O_O_H_H` | 0.12 | 1e-12 relative |
| `h2o_photolysis_yield_sum_lyman_alpha` | Their sum | 1, which is what makes the merged OH channel a partition of the absorptions rather than a subset of them | 1e-12 relative |

Tolerance is round-off: these are transcribed constants, not fits. The driver
prints the band index and its edges on a `note:` line. MEASURED 2026-09-18: 6
assertions, all pass.

### 35. Every metal photoionization cross section (`metal_photoionization_fits.f90`, `metal_photoion_table_transcription.py`)

Origin: item REF-METALS of `docs/Update_EXHALE_stage2.md` (`cross_sec.f90`,
`Cool_coeff.f90`).
Production routine: `metal_photoion_sigma`
(`src/modules/functions/cross_sec.f90`), the one place the photo-table column
order is defined and the one entry through which `set_energy_vectors` fills
`sigma_tab`.

Two published fits, each where it is valid. Verner, Ferland, Korista and
Yakovlev (1996), ApJ 465, 487, Table 1 and eqs. (1) to (4) are the
OUTER-SHELL fit on `E_th <= E < E_max`; Verner and Yakovlev (1995), A&AS 109,
125, Table 1 and eq. (1) are the same expression with `y_0 = y_1 = 0` and the
orbital quantum number kept in Q, the fit of the same outer shell at
`E >= E_max` and the fit at every energy for K I, for which the 1996 Table 1
has no row. `E_max` is the threshold of the outermost inner shell: Verner's
own routine `phfit2` switches the outer shell there (its variable `einn`) and
sets it to zero for potassium. `metal_photoion_sigma` returns the TOTAL
photoabsorption, the outer shell plus every subshell open at that energy,
which is what a caller of `phfit2` obtains by summing it over `is = 1 ... 7`,
and the published expression the driver compares against is built the same
way. The parameter rows are typed into the driver independently of the table
installed in `cross_sec.f90`, from
`references/verner_photo/photo.dat` and `references/verner_photo/table1.dat`,
so this driver is a second transcription and not a copy of the first.

Seventeen ions: C I, C II, O I, O II, N I, N II, Mg I, Mg II, Si I, Si II,
Ca I, Ca II, Na I, K I, S I, Fe I, Fe II.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `sigma_<ion>_at_{1.05,2.00,8.00}_Eth` | The cross section at 1.05, 2 and 8 times the threshold, so that the onset rise `((x-1)^2 + y_w^2)`, the `y^(-Q)` fall and the `(1 + sqrt(y/y_a))^(-P)` tail are each exercised | the published expression evaluated with the published constants | 1e-10 relative |
| `sigma_<ion>_at_{120,350,1240}eV` | The same at three fixed XUV energies, the last the top of the photon grid | the same | 1e-10 relative |
| `sigma_below_threshold_<ion>` | The cross section at 0.999 of the threshold | 0: there is no bound-free absorption below it | exact |
| `sigma_below_Emax_is_1996_<ion>`, `sigma_above_Emax_is_1995_<ion>` | 1e-6 eV under and over `E_max` | the 1996 fit, and the 1995 fit plus every subshell whose turn-on is `E_max` itself | 1e-10 relative |
| `sigma_is_1995_everywhere_<ion>` | K I, whose `E_max` is zero | the 1995 fit from the threshold up | 1e-10 relative |
| `sigma_{under,over}_subshell_turn_on_<ion>_<m>` | Every subshell turn-on, 1e-6 eV on either side | the published sum: under it the shell contributes nothing, over it exactly its 1995 fit | 1e-10 relative |
| `sigma_rises_at_subshell_turn_on_<ion>_<m>` | The step across each turn-on | strictly positive: a shell that opens raises the total | strict |
| `sigma_column_zero`, `sigma_column_past_end` | A column index outside the table | 0 | exact |
| `sigma_zero_below_ethr_<ion>`, `sigma_positive_above_ethr_<ion>` | Each fit 1e-6 eV under and over the threshold `species_table` carries | 0, and strictly positive: ONE definition of each metal threshold, since the photon grid puts `mion_ethr` on a bin EDGE and charges `h nu - mion_ethr` against it, so a cross section turning on at its own copy mis-charges the band between the two | exact, strict |

Mg II is the one ion whose 1996 row IS its 1995 row, so its two sides of
`E_max` agree by construction and the switch is invisible in it. The
threshold window of 1e-6 eV is narrower than any gap this can catch (the one
it did catch was Mg II, 0.005 eV) and wider than the double-precision spacing
of these energies.

Nothing here is a statement about what the absorption is CHARGED to: the
model advances the ion one stage and hands the degradation cascade one
electron of `h nu - E_th(outer)` whichever shell absorbed, and that
under-count, with the neglected fluorescence, is stated at the table in
`cross_sec.f90`.

`metal_photoion_table_transcription.py` is the separate field-by-field
comparison of the INSTALLED table with the author's distribution. It runs
`src/utils/metal_photoion_table.py --check`, which reads the emitted block
back and compares every field, `E_max` and the handover included, with
`photo.dat`, `table1.dat` and the `PH1` DATA statements of `phfit2.f`. The
distribution lives outside the repository, so its absence is reported as a
`SKIP` line with the path that was looked for and not as a failure; what is
then untested is the transcription, and every other assertion about these
cross sections still runs. MEASURED 2026-09-18: 368 Fortran assertions and 23
Python ones, all pass.

### 36. Every recombination coefficient of the balance (`recombination_coefficient_fits.f90`)

Origin: item REF-METALS of `docs/Update_EXHALE_stage2.md`.
Production routines, all of `src/modules/radiation/Cool_coeff.f90`:
`alpha_rr_metal`, `alpha_dr_metal`, `alpha_rec_metal`, `alpha_rec_FeI_Huang`,
`alpha_rec_FeII_Huang`, `alphaB_HII_new`, `alphaB_HeII_new`,
`alphaB_HeIII_new` and `alpha1_HeII_mao`.

Which source each stage takes: Badnell radiative recombination (Badnell 2006,
ApJS 167, 334; the author's current tabulation, `references/badnell_rr`),
`alpha_RR = A/[sqrt(T/T0)(1+sqrt(T/T0))^(1-b)(1+sqrt(T/T1))^(1+b)]` with
`b = B + C exp(-T2/T)`, for fourteen metal daughters (C I, C II, N I, N II,
O I, O II, Mg I, Mg II, Si I, Si II, Ca II, Na I, K I, S I) and H I, He I,
He II; Badnell dielectronic recombination (`references/badnell_dr`),
`alpha_DR = T^(-3/2) sum_i c_i exp(-E_i/T)`, for the same fourteen and He I;
Mao and Kaastra (2016), A&A 587, A84, eq. (9) at n = 1, subtracted from the
total to make case B for H I, He I and He II and used on its own as the He I
ground-capture coefficient; Verner's `rrfit.f` (version 4, 1999) power law,
whose Ca rows are those of Shull and Van Steenberg (1982), ApJS 48, 95, for
Ca I alone, the Badnell tabulation reaching the Mg-like sequence and having
no K-like row and neither set having a dielectronic fit there; and Huang et
al. (2023), ApJ 951, 123, eqs. (5) and (6), their fits to the Nahar iron
data, for Fe I and Fe II, where the Badnell dielectronic project stops short
of the Mn-like and Cr-like sequences.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `alpha_rr_<daughter>_at_{3000,10000,30000}K` | `alpha_rr_metal` of each of the fourteen daughters | the Badnell radiative form evaluated here from the published row | 1e-10 relative |
| `alpha_dr_<daughter>_at_{3000,10000,30000}K` | `alpha_dr_metal` of the same | the Badnell dielectronic sum from the published coefficients | 1e-10 relative |
| `alpha_rec_<daughter>_at_{3000,10000,30000}K` | `alpha_rec_metal` | the sum of the two | 1e-10 relative |
| `alpha_rr_CaI_at_{...}K`, `alpha_dr_CaI_at_{...}K` | Ca I | `1.120e-13 (T/1e4)^(-0.900)`, which is `rrec(:,20,20)`, the row of the twenty-electron neutral calcium the recombination forms; and zero, there being no dielectronic term | 1e-10 relative, exact |
| `alpha_rec_FeI_at_{...}K`, `alpha_rec_FeII_at_{...}K` | Iron | Huang eqs. (5) and (6) written out in the driver | 1e-10 relative |
| `alphaB_HII_at_{...}K`, `alphaB_HeII_at_{...}K`, `alphaB_HeIII_at_{...}K` | The three case-B coefficients | the Badnell total, with the He I dielectronic series where it applies, less the Mao and Kaastra n = 1 direct capture | 1e-10 relative |
| `alpha1_HeII_at_{...}K` | The He I ground-capture coefficient | the Mao and Kaastra n = 1 expression alone | 1e-10 relative |
| `alpha_rec_unknown_daughter` | `alpha_rec_metal('CaIII', 1e4)` | 0: a daughter the tabulation does not carry is not a rate | exact |

Every row is a transcription, so the tolerance is round-off and nothing else.
The three temperatures 3e3, 1e4 and 3e4 K span the wind: the low one is where
the dielectronic term of a near-neutral stage is still switching on, the high
one where it dominates.

**RED at the time this was written** on the Ca I radiative row: it read
`(6.78e-13, 0.80)`, which is `rrec(:,20,19)`, the rate forming Ca II. The
rate was 6.05 times the published one at 1e4 K. MEASURED 2026-09-18: 151
assertions, all pass.

### 37. Who receives the energy of the molecular network (`molecular_energy_recipients.f90`)

Origin: item L31 of the LHS 1140 b stationary campaign
(`docs/lhs1140b_stationary_L31_energy_cycles_20260917.md`), which covers all
eight blocks A to H: its section 5.3 is the uncertainty bracket of block G
and its section 5.5 the collider sum of block H, and it is the item that took
the driver to its present 81 assertions. Block H was measured again on a
solved state by L34b section 2.4
(`docs/lhs1140b_stationary_L34b_20260918.md`), and the `EXHALE_PROBE_STATE`
reader that block G's three cells now come from is L37 section 6
(`docs/lhs1140b_stationary_L37_20260918.md`).
Production routines: `molecular_chemical_heating`,
`molecular_reaction_energy_eV`, `dissociative_recombination_n2_eV` and
`dissociative_recombination_n2_source`
(`src/modules/lower_atmosphere/molecular_reaction_heat.f90`);
`k3b_H_H_to_H2`, `k3b_H_H_to_H2_atomic_H`, `k3b_H_H_to_H2_monatomic`,
`h2_association_collider_density`, `rk_R12_H2_thdis`, `rk_R15_3body_H2`,
`keq_H_H_to_H2`, `rk_R5_H2p_dr`, `rk_R6_H3p_dr_H2`, `rk_R7_H3p_dr_3H` and
`rk_R16_HeHp_dr` (`src/modules/lower_atmosphere/mol_rates.f90`);
`h2_total_decay_rate_max`, `h2_vibrational_heat_fraction`, `gamma_10_H`,
`gamma_10_H2`, `gamma_10_He`, `e_vib_v1_eV`, `e_vib_v5_eV` and `e_vib_v6_eV`
(`src/modules/lower_atmosphere/h2_vibrational_relaxation.f90`); and
`mol_heh_rows`, `set_mol_coeffs` and `h2_third_body_density`
(`src/modules/nonlinear_system_solver/System_HeH_mol.f90`).

Three channels leave a product in an excited state instead of giving the gas
kinetic energy: the dissociative recombinations R5 (H2+ + e) and R16
(HeH+ + e) leave one hydrogen atom in n = 2 (Takagi 2002, Phys. Scr. T96, 52;
Guberman 1994, Phys. Rev. A 49, R4277), and R6 (H3+ + e) leaves its H2
fragment vibrationally hot with a distribution peaking at v = 5-6
(Kokoouline, Greene and Esry 2001, Nature 412, 891; Strasser et al. 2001,
Phys. Rev. Lett. 86, 779). A fourth correction is of a rate rather than a
recipient: the third body of the R12/R15 pair is resolved by collider (Cohen
and Westberg 1983, J. Phys. Chem. Ref. Data 12, 531, p. 559) instead of the
total density carrying the M = H2 coefficient. A fifth is of the data the
thermalized fraction is built from: the all-level radiative rate of the
code's own ladder in place of the v = 1 one, and three published quantum
calculations in place of two 1979 fits (Lique 2015 for H, Jozwiak et al. 2024
for He, Le Bourlot et al. 1999 for H2).

**Block A, the enthalpy is untouched by the change of recipient.**
`enthalpy_R5_unchanged`, `enthalpy_R16_unchanged` and `enthalpy_R6_unchanged`
against 10.9478, 11.7534 and 9.2500 eV (1e-3, 1e-3 and 2e-2 absolute), the
values the audit of the network READ from the same table.
`n2_excitation_is_the_table_entry` against
`species_formation_energy(H_n2)` exactly and
`n2_excitation_is_0.75_IP_H` against `0.75 e_th_HI` (1e-14 relative).
`R5_heat_plus_excitation_is_enthalpy`,
`R16_heat_plus_excitation_is_enthalpy` and
`R6_prompt_plus_internal_is_enthalpy`, each zero at 1e-13: heat plus radiated
or excited is the reaction energy of the one formation table, so a recipient
change that moved an enthalpy would fail here. `R5_kinetic_share` 0.749 eV,
`R16_kinetic_share` 1.554 eV (2e-3 absolute), and
`R6_internal_share_is_the_v5_v6_peak` 2.4795 eV (1e-3 absolute), the peak of
the H3+ product distribution between v = 5 and v = 6 of the code's own
ladder.

**Block B, R12 and R15 are an exact detailed-balance pair, collider by
collider.** At 300, 800, 1500 and 3000 K, for a gas of one third body:
`detailed_balance_M_{H2,H,He}_T<T>` against `keq_H_H_to_H2` (1e-13 relative),
and `collider_{H2_is_k1_H2,H_is_k1_H,He_is_k1_monatomic}_T<T>` against the
single coefficient times the density (1e-13 relative). Then the efficiency
ordering at 1000 K: `efficiency_H_over_H2_at_1000K` 1.983 and
`efficiency_monatomic_over_H2_at_1000K` 0.4281 (1e-3 relative).

**Block C, the collider sum reduces to what it replaced.**
`collider_sum_is_ntot_in_pure_H2`: with n(H) and n(He) zero and n(H2) the
whole density, the H2-equivalent third body is the total density exactly, so
the correction is a change of mixture and not of normalization.

**Block D, the radiative and collisional data are the published ones.**
`ladder_reduction_has_run` (exact), `all_level_A_max_matches_the_ladder`
against the driver's own inventory of the 302-level ladder (1e-3 relative)
and `all_level_A_max_over_v1_total` 6.556 (2e-3 relative), which is the
statement that the rate is the all-level maximum and not the v = 1 value.
`lique_H_thermal_v1_v0_at_300K` against Lique (2015) Table 2's 1.8e-13
cm^3 s^-1 (2e-2 relative); `jozwiak_He_thermal_v1_v0_at_{808,1000}K` against
2.64e-15 and 7.02e-15, obtained independently from the same VizieR tables
(1e-2 relative); `lebourlot_H2_thermal_v1_v0_at_{1000,1500}K` against
2.631e-15 and 2.302e-14 (2e-3 relative), from the state-to-state table the
Cloudy distribution carries for a para-H2 perturber, which is the ONE H2-H2
set the paper published. `jozwiak_He_detailed_balance_{0300,1000,1500,2000}K`
takes the excitation and de-excitation rows READ from
`references/Jozwiak_2024J_A+A_685_A113/ph2-rat.dat` and asserts their ratio
against the Boltzmann factor of the v = 1 term energy the module itself
carries, both levels having j = 0 and so equal degeneracies; Le Bourlot et
al.'s file distributes only the de-excitation direction, so the same check
cannot be made on the H2 table. The ordering of the three colliders is
`helium_share_of_the_collider_coefficient` 9.341e-4 (2e-3 relative),
`atomic_H_is_three_orders_above_H2` (at least 1e3) and
`H2_is_no_better_a_quencher_than_He` (at most 1).

**Block E, the ledger deposits what the recipients leave it.** Each cell
carries ONE channel, every density that would open another being zero, so the
assembly's whole output is that channel's. `ledger_R5_deposits_enthalpy_less_n2`,
`ledger_R16_deposits_enthalpy_less_n2`,
`ledger_R6_withholds_the_internal_share` and
`ledger_R15_third_body_is_atomic_H`, each against the production rate
coefficients times the densities (1e-12 relative);
`heat_fraction_vanishes_without_colliders` exactly, which is the row that
would catch an internal share deposited promptly.
`n2_source_is_k5_nH2p_ne_plus_k16_nHeHp_ne` (1e-14 relative),
`n2_source_is_present` and `n2_source_vanishes_without_molecular_ions`: the
H(n=2) source is the same two rates and the same densities the ledger
charged, so one event is not counted at two rates.
`ledger_R5_same_without_excited_H` holds the deposit fixed when
`use_excited_H` is off: with the level carried the 10.199 eV reaches it,
without it the energy leaves as radiation, and either way the gas is charged
the same.

**Block F, the complete cycles close.** Four routes turn one H2 into two H
atoms, driven by one photon, and each must liberate `IP(H2) - D0(H2)` between
the gas, the H(n=2) level, the 153 nm photon of R23 and the H2 internal
energy the ledger withholds. `cycle_photoionization_then_R5_closes`,
`cycle_photoionization_R8_R6_closes`,
`cycle_photoionization_HeHp_R16_closes`,
`cycle_R23_nascent_H2p_then_R5_closes`,
`cycle_R15_then_R12_liberates_nothing` and
`the_four_routes_liberate_the_same_energy`, at 1e-9 absolute. The gas term of
each step is MEASURED, not restated: it is the production assembly's own
output on a cell carrying that step, divided by the step's event rate. A
share deposited twice shows as a surplus, a share dropped as a deficit.

**Block G, the uncertainty record on three cells.** Five quantities of the
ledger are approximations rather than measured numbers, and each is varied
HERE ALONE, on the state the code meets, so that the heat bracket of each can
be read off separately: (a) the recipient fraction of the 10.199 eV that
leaves the gas at R5 and R16, taken to zero, which is the whole recipient
correction undone; (b) the R6 branching, the share of R6 events that leave
the fragment excited at all; (c) the R6 internal energy at the +/-15 per cent
published uncertainty of the PEAK's position, the distance from the peak to
the mean of the distribution being quantified by neither source and not
inside this bracket; (d) the quench fraction from its model value to 1, the
thermalized limit; (e) the helium efficiency, argon's coefficient scaled by
Cohen and Westberg's stated +/-0.3 in log, which moves both directions of the
R12/R15 pair. The nominal heat is the production assembly's and every
variation is an exact difference against it, formed from the same production
rate coefficients and reaction energies, so no term of the ledger is restated
in the driver. Asserted per cell:
`recipient_bracket_is_the_n2_source_<case>` against
`dissociative_recombination_n2_source` times the excitation energy (1e-13
relative), which is the statement that one event gives one excitation counted
in one place; `bracket_contains_the_nominal_<case>` and
`bracket_below_contains_the_nominal_<case>` (at least 0). Then
`quench_bracket_is_closed_at_the_base`, `1 - f` at most 1e-5: at the
molecular depth the ladder is thermalized, so the scalar fraction and a
level-population model cannot differ there by anything the run reports.

The three cells are a RULE, not a list of numbers: **molecular depth** is
physical cell 1, the base of the column; the **H2 front** is the first cell
where `2 n(H2)/n_H` has fallen to half its value at cell 1; the **dilute
upper column** is the cell nearest `r = 1.842260` R_p, and that radius IS the
choice, no threshold on the profile picking that cell, because `2 n(H2)/n_H`
is flat and rising outward there. `EXHALE_PROBE_STATE` names the directory
the three cells are read from, and the pair read is
`{Hydro_ioniz_IC.txt, Ion_species_IC.txt}`, the state a run hands back; the
pair without `_IC` is the product written FROM a measurement of that state
and is a different state. The species columns are taken BY NAME from the
file's own `# columns` header, and the electron density and the particle
count are formed as the code's mass policy forms them: a column whose name
ends in III carries two charges and one ending in II carries one, H2p, H3p
and HeHp carry one each, and HeITR is a level inside He I and enters neither
sum. Unset, nothing is read and the compiled cells stand, so the suite runs
with no state in the tree; the driver prints the distance of the read cells
from the compiled ones as one number. Pointed at an atomic state it stops
with `FAIL probe_state_has_no_column_H2`, and at a directory with no state
pair with `FAIL probe_state_file_unreadable`, rather than reading a zero.

**Block H, one collider sum evaluated as each of its three callers calls
it.** The chemistry forms it inside the rows
(`System_HeH_mol::h2_third_body_density`, at the composition the row is
evaluated at), the carrier balance reaches the same call through
`mol_heh_rows` and writes it into the R15 channel of its row record, and the
heat ledger forms it in `molecular_chemical_heating`.
`collider_sum_chemistry_equals_the_routine` and
`collider_sum_carrier_row_equals_the_routine` are exact, and
`collider_sum_heat_ledger_equals_the_routine` is read out of the assembly on
a cell whose only open channel is the thermal dissociation R12, where the
whole deposit is `k12 n_third n(H2) q(R12)`, so the third body the ledger
used is recovered without dividing anything (1e-14 relative). If any one of
them formed its own, the H2 abundance would be set at one rate and its heat
charged at another.

The four corrections are turned off together by
`EXHALE_REACTION_HEAT_RECIPIENTS=0`; running this driver with that key set is
the control, and the rows of A, C (the mixture), D and E that state the
correction then FAIL by construction. MEASURED 2026-09-18: 81 assertions, all
pass with the key unset.

### 38. One three-body rate for the composition and for the energy (`molecular_third_body_in_the_chemistry.f90`)

Origin: the item "one three-body rate for the H2 chemistry and for its heat"
of `docs/Update_EXHALE_stage2.md` (2026-09-16), with
`docs/lhs1140b_stationary_L7g_model_20260916.md`.
Production routines: `mol_heh_rows`, `set_mol_coeffs` and
`set_mol_turnover_rates`
(`src/modules/nonlinear_system_solver/System_HeH_mol.f90`);
`h2_association_collider_density`, `k3b_H_H_to_H2`,
`k3b_H_H_to_H2_atomic_H`, `k3b_H_H_to_H2_monatomic`, `rk_R15_3body_H2`,
`rk_R12_H2_thdis` and `keq_H_H_to_H2`
(`src/modules/lower_atmosphere/mol_rates.f90`).

R15, H + H + M -> H2 + M, and its reverse R12 run on a COLLIDER SUM,
`k1(H2) n(H2) + k1(H) n(H) + k1(Ar) n(He)`, with the three recommended
coefficients of Cohen and Westberg (1983), J. Phys. Chem. Ref. Data 12, 531,
p. 559. Two of those densities are unknowns of the molecular balance system,
so the rate is a function of the state the row is evaluated at and cannot be
hoisted once per cell from the total heavy-particle density. If it is, the H2
abundance is set by one rate and its heat charged at another, and nothing in
either answer says so.

The cell is the certified LHS 1140 b molecular base in order of magnitude
(T = 808 K, n(H) = 1.87e12, n(He) = 6.17e12, n(H2) = 5.15e11 cm^-3), so the
numbers the assertions carry are the ones the code meets there.

| Assertion | Quantity | Reference | Tolerance |
|---|---|---|---|
| `one_k15_for_the_row_and_the_ledger` | The R15 rate the H2 row forms, read from the row's own channel decomposition | the rate the energy ledger forms at the same state, `rk_R15_3body_H2` of `h2_association_collider_density` | 1e-14 relative |
| `collider_sum_is_not_the_total_density` | The total-density coefficient over the collider-sum one | at least 1.1, so that the row above has content; at this cell it stands 27 per cent above | at least |
| `R12_and_R15_share_the_third_body` | R12's channel divided by R15's | `n(H2)/(K_eq n(H)^2)`, with no collider density left in it, which is what makes the pair an exact detailed balance whatever the mixture is | 1e-12 relative |
| `dR15_dnH2_is_the_H2_collider` | `d(R15 rate)/d n(H2)` by central differences of the production rows | `k1(H2) n(H)^2`, the analytic derivative of the collider sum; with the third body frozen it is zero | 1e-8 relative |
| `dR15_dnHI_has_the_collider_term` | `d(R15 rate)/d n(H)` | `k1(H2)[eff_H n(H)^2 + 2 n_third n(H)]`: the atom is both a third body and a reactant, and a hoisted coefficient keeps only the second term | 1e-6 relative |
| `dR15_dnHI_collider_term_is_resolved` | The collider term's share of that derivative | at least 1e-2: it is 20 per cent of the whole derivative at this cell and not a rounding of the reactant term | at least |
| `dR12_dnHI_is_the_atomic_third_body` | `d(R12 rate)/d n(H)` | `k12 eff_H n(H2)`, which is zero with a hoisted coefficient | 1e-6 relative |
| `H2_turnover_bound_rises_with_the_helium` | The H2 row's turnover bound when the cell's helium is raised | strictly larger, helium being a third body of this channel | at least 1 + 1e-6 |

The four derivative rows are the entries the numerical Jacobian of `hybrd1`
builds, so the check is on the quantity the solver actually uses. MEASURED
2026-09-18: 8 assertions, all pass.
