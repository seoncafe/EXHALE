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
(`-O0 -g -fcheck=all -fbacktrace -fopenmp`). A full build and run takes about
ten seconds (measured 10.3 s, sixteen drivers and five Python groups,
2026-09-06).

Output is one line per assertion,

```
PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
```

and the exit status is nonzero if any assertion failed or any build failed.
Lines beginning `note:` are diagnostics, not assertions: they carry a number
that has no comparison attached.

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
| `h2_infrared_line_populations.f90`, `h2_partition_sum_uniqueness.py` | Group 18. |
| `lya_escape_probability.f90` | Group 19. |

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
| `lya_band_leaves_other_band_unchanged_{LW,B1,B3,B4}` | max relative change of each other band between `n_HI = 0` and the layer | 0: only B2 has H I as a line absorber | 1e-12 absolute |
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
the sixteen it has now: the band read the beam at the cell's star-ward face,
which is 1.095 times the cell mean in the extended layer, 1.005 times it in
the graded column and **18.2 times** it in the base cell of the molecular
base, and the column photon budget stood 1.2 per cent above the beam's own
loss. **GREEN on all sixteen** with the cell mean in place; the base ratio of
the extended layer is then 0.2066.

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
