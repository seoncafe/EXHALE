# D0: the governing system of EXHALE v1.00, on paper (2026-09-06)

Phase 2 deliverable of `docs/development_plan_20260905_rev3.md` section 4.2,
Step 3 of `docs/development_plan_20260905_execution.md`. This document records
**what the code solves today**, before discretization where the equations
permit it and at the discrete level where the code has no continuous
counterpart. It proposes no fix and adopts no new physics. Every statement
about behavior carries the file and line that produces it.

## 0. How to read this document, and what its labels mean

**Provenance.** Every statement below is READ from the source tree at
`/nfs/mocafe/kiseon/RT_Codes/ExoAtmosphere/EXHALE_v1.00/` as it stood on
2026-09-06 (HEAD `35d9dd5`, with the uncommitted Phase 1 batches 2a, 2b and 2c
in the working tree). **Nothing in this document was MEASURED here**: no build,
no run, no test, no regression case was executed while writing it. Where a
number appears it is one of

- a constant or expression read from the source (labeled by its file and line);
- a value quoted from a source comment, which records a measurement made when
  that code was written (marked "code site records");
- a value quoted from `docs/physics_numerics_audit_20260905.md`,
  `docs/code_status_20260905.md` or `docs/Update_EXHALE_stage2.md` (marked with the
  document).

No quoted number was re-measured. Where a document and the code disagree, the
code is authoritative and the disagreement is recorded as its own item.

**Citations are by file and routine, not by line (2026-09-07, item
DOCS-LINES).** This document was written with `file.f90:NNN` pointers, and the
series that followed moved every one of them within days. They were rewritten
to the file and, where the file holds more than one routine, the routine, the
declaration or the labelled comment block. Two conventions make the short forms
unambiguous: a bare `input_read.f90` is the key parser `input_read` unless
another routine is named, and a bare `EXHALE_main.f90` is the main program
`Hydro_ioniz`, whose blocks carry the banner comments this document quotes
("Main temporal loop", "THE ATTEMPTED STEP (B3a)", "Ionization Equilibrium",
"THE ADOPTION BOUNDARY", "Staged secondary ionization") and whose fourteen
step operations are numbered by the `as_op_*` tags of `attempted_step.f90`.
Statements about code that the series has since changed are kept as the dated
record they are, with the change and its `docs/Update_EXHALE_stage2.md` item named at
the site.

**Term entries.** Each term is given as (a) its physical meaning, (b) its
source or the statement that it is a model choice, (c) its validity range or
the approximation it makes, (d) the routine and line that computes it now, and
(e) a consistency verdict. Verdicts are OK, or a reference to a numbered
consistency item `C<n>` in section 8, which is where every inconsistency this
document found is collected.

**One-line verdict.** EXHALE solves the spherically symmetric single-fluid
Euler equations with a composition-dependent caloric equation of state, a
locally solved chemical state, an attenuated stellar radiation field, and up to
three operator-split composition transport operators. The composition, the
radiation field and the thermal state are advanced by a Lie splitting whose
pieces do not all see the same state at the same time level, and the steady
solver enforces a strictly smaller system than the marching loop advances.
Those two facts, and their consequences, are the subject of sections 6 and 7.

---

## 1. Unknowns and state

### 1.1 What is actually advanced

The conserved vector has **three components and no species**:

```
u = (rho, rho*v, E),   E = (1/2) rho v^2 + e_int(rho, p, composition)
```

declared `src/EXHALE_main.f90`, allocated `u(3, 1-Ng:N+Ng)` at
`src/EXHALE_main.f90`, with `Ng = 2` (`src/modules/init/parameters.f90`)
and `N` a runtime cell count (`parameters.f90`, key `Grid cells:` parsed at
`src/modules/files_IO/input_read.f90`). The three rows are fixed by the
conversions `W_to_U` / `U_to_W`,
`src/modules/functions/UW_conversions.f90`. The gravitational
potential is **not** inside `u(3)`; gravity enters as a momentum source and as a
flux-form work term (section 2.3).

Formation (chemical) energy is **not** stored in `u(3)` either. This is the
plan's decision 2 as it stands today: the evolved variable is kinetic plus
thermal internal energy, and chemical energy appears only through source terms
(section 5). Decision 2 is still open (section 9).

**Normalization** (`input_read.f90`, symbols declared
`parameters.f90`): lengths in `R0`, velocities in
`v0 = sqrt(k_B T0 / mu)` with `mu = m_H`, densities in `rho0 = n0 mu`, pressures
in `p0 = n0 mu v0^2 = n0 k_B T0`, temperatures in `T0`, times in `t_s = R0/v0`,
volumetric rates in `q0 = n0 mu v0^3 / R0`, and the Jeans parameter
`b0 = G M_p mu /(k_B T0 R0)`. `n0` is the base H+He nucleus density and must have
exactly one source, else the run stops (`input_read.f90`).

### 1.2 The composition array

All composition lives in one array,

```
f_sp(1-Ng:N+Ng, n_species),   n_species = 40      parameters.f90
```

holding **species density per unit mass** in code units. Its column map is the
canonical species table `src/modules/init/species_table.f90`:
`isp_HI=1 ... isp_HeTR=6`; 27 metal ion stages
(three stages, neutral / +1 / +2, for C, O, N, Mg, Si, Ca, Fe; two for Na, K, S,
`melem_top`); `isp_H2=34, isp_H2p=35, isp_H3p=36, isp_HeHp=37`;
`isp_OH=38, isp_H2O=39, isp_CO=40`. Species mass, charge
and nucleus counts are the `bsp_mass`, `bsp_charge` and `bsp_nH`/`bsp_nHe`/
`bsp_nO`/`bsp_nC` parameter arrays of the same file; element weights the
`melem_A_u` and `melem_A` arrays.

`bsp_is_excited_level` (`species_table.f90`) is `.true.` for He 2^3S
alone: the metastable is a **level inside He I**, and `n(He I)` is carried as
the total neutral helium population including it. Every budget sum (mass,
particle count, element nuclei, free electrons, collision partners) must skip
that column while the level population is still solved and, under diffusion,
transported. The rule is stated at `species_table.f90` and is the
reason the ground singlet is formed by difference in exactly one place,
`src/modules/functions/composition.f90`.

### 1.3 Transported, locally solved, or diagnostic

| Species | Class | Line establishing it |
|---|---|---|
| H I | local equilibrium; closes the hydrogen budget | `System_HeH_mol.f90`; write-back re-closure `diffusive_photochemistry.f90` |
| H II | local equilibrium by default; **transported** with `Ionization transport: True` | default `.false.` `parameters.f90`; carrier `diffusive_photochemistry.f90`; imposed on the sweep `ionization_equilibrium.f90` |
| He I, He II, He III | local equilibrium; the He **element total** moves under `He_diffusion` | stages `ionization_equilibrium.f90`; element `EXHALE_main.f90` |
| He 2^3S | local equilibrium, as a level of He I | `ion_residual_core.f90` |
| H(n=2) | **diagnostic plus lagged-explicit feedback**; never an unknown of any solve | `excited_hydrogen.f90`; populations; feedback only as an added rate `ionization_equilibrium.f90` and an added heat |
| H2 | **transported** whenever the carrier operator runs; else local equilibrium | carrier 1 `diffusive_photochemistry.f90`; row `System_HeH_mol.f90` |
| H2+, H3+, HeH+ | local equilibrium, explicitly excluded from the carriers | `diffusive_photochemistry.f90`; rows `System_HeH_mol.f90` |
| OH, H2O | transported with oxygen chemistry; else rows `iox`, `iox+1` | carriers; rows `System_HeH_mol.f90` |
| CO | **no balance row at all**: transported reservoir, or a closed-form equilibrium density | `ionization_equilibrium.f90` vs (`co_equilibrium_density`, `oxygen_rates.f90`) |
| 27 metal ion stages | local equilibrium; element totals move only under `he_metal_diffusion` (default `.false.`, `parameters.f90`) | rows `ion_residual_core.f90`; exclusion `diffusive_photochemistry.f90` |
| electrons | **derived**, never an unknown: the charge-neutrality sum | `calc_ne` `utilities.f90`; in the residuals `System_HeH_mol.f90` |

### 1.4 Where the state lives

`module global_parameters` (`parameters.f90`) does **not** hold the
hydrodynamic state. `u, u1, u2, u_old, W, WL, WR, dF, S, rho, v, E, p, T, cs,
f_sp, nhi ... nm, ne, n_tot, heat, cool, eta, mom, dt_loc, Rres` are local
allocatables of `program Hydro_ioniz`, declared `EXHALE_main.f90` and
allocated in `allocate_state_vectors`, `EXHALE_main.f90`.
`global_parameters` holds the grid (`r, r_edg, dr_j`, `parameters.f90`),
the gravitational potential samples (`Gphi_c, Gphi_i`), the opacity
multiplier `opa_pf`, the eddy coefficient `kzz_cell`, the
H(n=2) feedback arrays `gph_balmer_HI, heat_balmer` and the
photon grid and cross-section vectors.

Two mirror arrays live in `ionization_equilibrium.f90`: `nmol_eq(:,4)`
and `nox_eq(:,3)`, in cm^-3. **These, not `f_sp`, are what
`write_output` prints for the molecular and oxygen columns**
(`write_output.f90`), which is why every write is preceded by
`molecular_carrier_densities_from_state` (`ionization_equilibrium.f90`;
call sites `EXHALE_main.f90`).

### 1.5 The packed cell state

`src/modules/nonlinear_system_solver/ion_cell_state.f90` defines three derived
types, each with one saved module instance declared `!$omp threadprivate`,
because the ionization sweep is parallel over cells and a
threadprivate copy must need no allocation.

- `type ion_rates` / `ieq_cell`: the photoionization rates
  `P_HI, P_HeI, P_HeII, P_HeITR, P_H2` and the three disjoint subsets of
  `P_H2` (`P_H2_di`, `P_H2_dd`, `P_H2_nd`); `k_LW` with self-shielding
  already applied; recombination, collisional ionization, and the
  triplet channels `A31, q13, q31a, q31b, Q31`; `nh, nhe, T_K, ntot`, where
  `ntot` is the total gas particle count from `calc_ntot` and is the third body
  M of R12/R13/R15; the He/H charge exchange pair
  `kcx_He0_Hp, kcx_Hep_H0`; the oxygen family counts; and
  the three independent imposed-partition flags `x_h2_fixed`, `x_ox_fixed`,
  `x_hp_fixed` with their values.
- `type adv_rates` / `adv_cell`: the post-process advection cell
  state, carrying `c1 = dr/v`, the upstream populations and `xe_metal`. Note
  `xheiS_old` is the ground singlet alone.
- `type teq_state` / `teq_cell`: the post-process temperature
  residual state.

---

## 2. Hydrodynamics

### 2.1 The system

Spherically symmetric, single fluid, finite volume with area `A = r^2` and
volume element `dV = (r_+^3 - r_-^3)/3`:

```
d_t rho     + (1/r^2) d_r (r^2 rho v)           = 0
d_t (rho v) + (1/r^2) d_r (r^2 rho v^2) + d_r p = -rho dPhi/dr
d_t E       + (1/r^2) d_r (r^2 v (E+p))         = -rho v dPhi/dr
```

Flux differences are formed in `src/modules/time_step/RK_rhs.f90`: the
face loop calls `Num_flux` for each face and stores `face_flux`/`face_p`; the
cell loop then differences them with the cell geometry `A_p`, `A_m`, `dV`.

(a) Mass, momentum and energy conservation of a compressible inviscid flow.
(b) Textbook; the ATES v2.0 heritage (Caldiroli et al. 2021) is the origin of
the scheme. (c) Single fluid, one temperature, one velocity, spherical
symmetry; valid where the gas is collisional and the species are momentum- and
thermally coupled, which the code does not test at run time (audit section 5.7).
(d) `RK_rhs.f90`. (e) See C15 for the momentum row, C2 for the mass row.

### 2.2 Fluxes and reconstruction

Three Riemann solvers, selected by the required key `Numerical flux:`
(`input_read.f90`), dispatched at `Num_Fluxes.f90`:

| key | routine | notes |
|---|---|---|
| `LLF` | `lax_friedrichs_flux` (`Num_Fluxes.f90`) | Rusanov coefficient `alpha = max(|v_L|+a_L, |v_R|+a_R)`. Sources cited at the code site: Rusanov 1961; Perthame and Shu 1996; Zhang and Shu 2010. This is the Phase 1 batch 2b correction (A1); the previous `max(|v_L+a_L|, |v_R+a_R|)` omitted the negative acoustic characteristic (audit N2). |
| `HLLC` | `Num_flux` (`Num_Fluxes.f90`) | Davis/Einfeldt speeds, star states. Default choice for molecular production runs. |
| `ROE` | `Num_flux` (`Num_Fluxes.f90`) | Roe average, Harten-Hyman entropy fix, flux as Toro eq. 11.29; falls back to `hlle_flux` when `speed_estimate_ROE` returns no admissible star state, counted in `n_faces_roe_hlle`. |

`speed_estimate_ROE.f90` implements Toro (2009) ch. 9: PVRS eq. 9.20,
two-rarefaction eq. 9.32, two-shock eq. 9.42, switch `Q_user = 2` (eq. 9.43),
vacuum test eq. 4.76. **Stated validity** (module header): one ideal gas with a
constant adiabatic index on both sides. The consequence is enforced:
`input_read.f90` refuses `Numerical flux: ROE` together with
`Molecular chemistry: True` unless `Caloric EOS: monatomic`.

Across a face each side uses its own cell's `gamma_eff`
(`Num_Fluxes.f90`); only the Roe average and the star-state estimate take
the arithmetic mean `gam_face`.

**Low-Mach dissipation** (`src/modules/flux/low_mach_dissipation.f90`) is a
gated Jameson-Schmidt-Turkel (1981) fourth-difference stress **added to the
numerical flux** in `RK_rhs.f90`, so the marching right-hand side and
the steady residual solve the same equation:

```
D_p(j+1/2) = eps4 g(M) rho_f lambda_f (v_{j+2} - 3 v_{j+1} + 3 v_j - v_{j-1})
D_E(j+1/2) = v_f D_p(j+1/2),     g = [max(0, 1 - M_f^2/M_th^2)]^2
```

assembled in `low_mach_dissipation.f90`; the mass component is identically
zero. Key `Low-Mach damping: <eps4> [<M_th>]` (`input_read.f90`), default
`lowmach_damp_eps = -1.0`, off (`parameters.f90`). (c) The code site states
this is a numerical dissipation with no physical counterpart, admissible only
where it is negligible against the physical fluxes, and its size must be
reported (magnitude routine); stability bound
`eps4 < 1/(16 CFL)`.

**Reconstruction.** PLM with a generalized MC/minmod limiter, `theta = 2.0`
(`PLM_rec.f90`, limiter, cited as Kappeli 2016), or ESWENO3
(`Reconstruction.f90`). Key `Reconstruction scheme: PLM | WENO3 |
PLM+WENO3` (`input_read.f90`). The two-stage switch lives in the
marching loop, not in the reconstruction module: armed and fired there
(`EXHALE_main.f90`, the "PLM -> WENO3 continuation" block), with an optional homotopy
`R_lambda = (1-lambda) R_PLM + lambda R_WENO3` in
`steady_residual.f90` that reproduces each pure scheme bitwise at the
end points.

Positivity guards, which change the discrete flux but are not part of the
continuous system: face-state limiter `Reconstruction.f90` (Zhang and
Shu 2010) and the single-interface first-order flux repair
`RK_rhs.f90` (Hu, Adams and Shu 2013; Stone et al. 2020 section 4.6).

### 2.3 Gravity

(a) The potential of the planet, plus the stellar tidal and centrifugal terms in
the co-rotating frame. (b) Model choice, selected by input; the Roche form is
the ATES heritage form. (d) `src/modules/functions/grav_field.f90`:

```
spherical (grav_field.f90, phi and Dphi)
    Phi  = -b0/r ,  dPhi = +b0/r^2

Roche, the default (:22-25, 41-44)
    Phi  = -b0/r - b0 Mrapp/(atilde - r)
           - b0 (1+Mrapp)/(2 atilde^3) (atilde Mrapp/(1+Mrapp) - r)^2
```

with `Mrapp = M_star/M_p`, `atilde = a_orb/R0` (`input_read.f90`).
Key `Domain mode: Spherical` with a mandatory `Outer radius [R_p]:`
(`input_read.f90`); the default is Roche
(`spherical_domain = .false.`, `parameters.f90`), with `r_max` at the L1
radius (`input_read.f90`).

Samples are pre-tabulated as **potentials** at cell centers and edges,
`set_gravity_grid.f90`. Gravity enters twice:

- momentum source `S(2) = -(1/2)(rho_lowerface + rho_upperface)
  (Gphi_i(j) - Gphi_i(j-1))/dr`, `src/modules/states/Source.f90`, with the two
  densities the reconstructed face states passed at `RK_rhs.f90`;
- energy, in **flux form**: `dF3p = dAp Fp(1) (Gphi_i(j) - Gphi_c(j))
  - dAm Fm(1) (Gphi_i(j-1) - Gphi_c(j))`, `RK_rhs.f90`, with
  `S(3) = 0` (`Source.f90`). This is the well-balanced form: for a constant
  mass flux it reduces to `rho v dPhi/dr` with the same discrete potential the
  momentum source uses.

(c) A spherically symmetric wind in a Roche potential is not a
multidimensional Roche-overflow solution; the geometry, the irradiation
average and the outer boundary must be stated together (audit section 5.7).
(e) OK, with the caveat of C31.

**There is no rescaling of gravity below the escape radius anywhere in the
hydrodynamics.** `r_esc` (`Escape radius [R_p]:`, `input_read.f90`) is a
diagnostic and convergence window only: `define_grid.f90` sets `j_min`,
guarded against an empty window and clamped to `max(j_min,1)`; the marching
convergence measure is the mass-flux spread over `[j_min:N]`
(`EXHALE_main.f90`); the residual norm is reported for the wind and
the layer separately and combined by the larger
(`steady_residual.f90`), with no region switch in any row scale. The
only gravity scaling in the tree is `b0_eff` in `set_IC.f90`, which
builds the initial hydrostatic column and never enters the evolution
equations.

### 2.4 Equation of state

**Thermal.** `p = (n_tot + n_e) T` in code units, composition-exact:
`comp_p_from_T` and `comp_T_from_p` (`composition.f90`), with
`n_tot` and `n_e` from the single policy point `get_species_densities`.
Mass per hydrogen nucleus `comp_mass_per_H`; base
particle count `comp_ntot_bc`, which subtracts the hydrogen nuclei
bound into H2 when the molecular base is on. Metals enter the bulk budget when
`eos_include_metals` (default true, `metals.inp` key `eos_metals 0|1`).

**Caloric.** `src/modules/states/caloric_eos.f90`:

```
e         = (3/2)(n_tot + n_e) k T + n(H2) k u_rv(T)
C_V/k     = (3/2)(n_tot + n_e) + n(H2) c_rv(T)
gamma_eff = 1 + (n_tot + n_e)/(C_V/k)
```

(module header; `internal_energy_of_mixture`,
`heat_capacity_of_mixture`, `adiabatic_index_at_T`, the
four hydro maps, and the Newton inverse `temperature_of_mixture`).

- Atomic gas: `gamma = 5/3` exactly, written as the legacy expression so an
  atomic run is bitwise the constant-index code (the `thereis_mol` gate and
  the constant-index branch of `caloric_state_from_composition`).
- Molecular gas: **only H2** carries internal degrees of freedom, and it carries
  the complete measured rovibrational ladder, 302 bound levels of X^1 Sigma_g^+
  up to 51966 K, from **Roueff et al. 2019, A&A 630, A58, table 2**, reused from
  `molecular_infrared_data`, tabulated in `ln T` and interpolated by a
  C^1 cubic Hermite because the JFNK residual
  differentiates this map.
- (c) Stated approximations: H2+, H3+, HeH+, OH, H2O and CO are counted as
  monatomic, justified by measured ratios at or below 1.7e-7 of n(H2);
  the ladder is the electronic ground state only, so dissociation is a chemical
  source and not a heat capacity; the ortho/para **equilibrium**
  mixture is used, agreeing with a frozen 3:1 mixture to four decimals above
  300 K; outside the table the interpolant continues linearly with the
  end-point heat capacity.
- Key `Caloric EOS: ladder | monatomic`, default `ladder`
  (`input_read.f90`, `parameters.f90`); `monatomic` is labeled a
  comparison option and acts by returning `u_rv = c_rv = 0`.
- (e) C22: the chemistry module builds a **different** H2 internal-state model
  (`mol_rates::q_rovib_H2`), so one thermodynamic potential does not generate
  both the equation of state and the equilibrium constants.

The sound speed used by the Riemann solvers is the **frozen-composition** one.
That is the correct fast-wave limit if chemistry cannot relax during an acoustic
disturbance; the equilibrium (reacting) compressibility is different, and the
stiff relaxation limit of the split scheme has not been shown to recover it
(audit section 6.5). Recorded as C32.

### 2.5 Viscosity and conduction

`src/modules/time_step/viscous_conduction.f90`. **Both default off**, so an
unkeyed run is byte-identical to the inviscid code.

(a) Newtonian viscous stress and Fourier conduction in spherical symmetry.
(b) The module header explicitly **rejects Koskinen et al. (2022)
ApJ 929, 52, eqs. (B5)/(B6) as printed** and implements the Navier-Stokes result
instead, because the implemented dissipation is manifestly non-negative and
satisfies `w F_mu + q_mu = div(tau . w)`, which the printed pair does not:

```
tau_rr = (4/3) mu (w' - w/r),   tau_theta = tau_phi = -(1/2) tau_rr
F_mu   = (4/3)(1/r^2) d_r(r^2 mu d_r w) - (4/3)(d_r mu)(w/r) - (8/3) mu w/r^2
q_mu   = (4/3) mu (w' - w/r)^2   >= 0
S_E    = w F_mu + q_mu + (1/r^2) d_r(r^2 kappa d_r T)
```

(d) `viscous_momentum_coeffs`, `viscous_dissipation`,
`thermal_conduction_coeffs`, `thermal_conduction_source`,
and the single pair `viscous_conduction_sources` used by **both** the
marching step and the steady residual.

Coefficients: `kappa(T) = 4.45e4 (T/1000 K)^0.7` erg cm^-1 s^-1 K^-1, Watson,
Donahue and Walker (1981), also quoted as Erkaev et al. 2016 eq. (A6);
`mu(T) = (4/15)(m_H/k_B) kappa(T)` through the Eucken ratio 15/4
(`viscous_conduction.f90`, the conductivity and viscosity coefficients).

(c) Stated validity: these are the **neutral atomic hydrogen** values.
Above the ionization front the Spitzer electron conductivity and the Coulomb ion
viscosity are much larger and are not covered; the terms are included for the
dense, largely neutral base.

(d, continued) They are **operator split**, not in `RK_rhs`:
`viscous_conduction_step` is a Crank-Nicolson stage called once per
marching step at `EXHALE_main.f90`: momentum diffusion solved
tridiagonally at fixed pressure so the kinetic-energy change is exactly the
viscous work, then `C(T*-T)/dt = (1/2)[Q(T)+Q(T*)] + q_mu`.
Boundaries: Dirichlet base ghost, zero
diffusive flux at the outer face. A temperature floor of
`0.01 T0` is applied and counted, and a floored cell is flagged as
**not** a zero of the steady residual.

### 2.6 Boundary conditions

**Inner (base): one variant only.** `src/modules/states/base_boundary.f90`
imposes a characteristic condition **at the face** `r_edg(0)`, with the number
of conditions set by the sign and magnitude of the face Mach number:

| regime | reservoir conditions | interior conditions |
|---|---|---|
| `0 < v < c`, subsonic inflow | 2 | 1 |
| `-c < v < 0`, reversal | 1 | 2 |
| `v > c`, supersonic inflow | 3 | 0, and not well posed |
| `v < -c`, supersonic outflow | 0 | 3 |

The interior relation is the LODI `C^-` compatibility condition
`p_b - rho_i c_i v_b = p_i - rho_i c_i v_i`, cited to Thompson 1987,
Poinsot and Lele 1992, and Carlson 2011 (NASA/TM-2011-217181) section 2. The
closed-form invariant is deliberately not used because the caloric equation of
state gives each cell its own `gamma_eff`.

The reservoir is `(p, s)` at the base level `r_base_level`, carried to the face
along its own hydrostatic isentrope (`set_base_reservoir`).
Constraints actually imposed (`characteristic_base_face_state`):
interior state continued to the face along cell 1's hydrostatic isentrope with
`rho v r^2` constant; reservoir continued to the face;
`p_b = p_res + (1/2) w (p_i - rho_i c_i v_i - p_res)` with
`base_incoming_invariant_weight = 0.0`, that is `p_b = p_res`, the
alternative recorded as measured and rejected; density from the
reservoir on inflow and from the interior isentrope on reversal, blended by a
C^1 smoothstep in face Mach number over `base_face_mach_blend = 1.0e-6`;
`v_b = v_i + (p_b - p_i)/(rho_i c_i)`; a
supersonic-outflow blend to the interior; a Mach cap
`base_face_mach_max = 5.0`, counted; and the reservoir alone
when cell 1 is inadmissible (`reservoir_face_state`).

**Face state and cell values are different objects.** `base_boundary_states`
returns the face state at `r_edg(0)`, the ghost **cell averages**,
and the point state at `r_edg(-1)`. The ghosts are 8-point Gauss-Legendre volume
averages of the hydrostatic isentrope through the face state
(`base_ghost_averages`), with `rho v r^2` held constant through them; they are
explicitly not copies of the face value, and the module header
records the half-cell defect that this removes (a measurement of
0.2502 to 0.2514 of `rho g` at cell 1 over an eight-fold refinement, driving a
first-order velocity artifact from -247 to -30 cm/s). Consumers: `Apply_BC_W`
writes the ghosts (`Apply_BC.f90`); `Rec_BC` puts the face state into
`WL_out(:,0)` and deliberately leaves `WR_out(:,0)`, the interior
reconstruction down to the base face, to the Riemann solver.

Base keys: `Base BC: density` (default) or `Base BC: pressure [<p_ubar>]`
(`input_read.f90`), which state where the base level is, not a boundary
variant. Four former base keys are hard errors
(`input_read.f90`, message). Supersonic inflow is
over-specified and is only **diagnosed**:
`check_base_inflow_is_subsonic` (`Apply_BC.f90`), called
`EXHALE_main.f90`. Recorded as C33.

**Outer: free outflow.** `free_outflow_ghost`, called from `Apply_BC_W`:
zero-gradient copy for PLM, linear extrapolation for WENO3 with a positivity guard that
falls back to zero-gradient and counts the event. Reconstructed
face states at `Rec_BC`. No inflow and no fixed-pressure outer
condition exists.

**Shapiro filter** (not a boundary condition): a 1-2-1 low-pass on the
conservative variables, `Apply_BC.f90`, key
`Shapiro filter: <eps> [<every>]`, default off (`parameters.f90`), whose
comment records that with the filter on the HD 209458 b cold-start solution is
driven off the transonic saddle into the infall attractor.

---

## 3. Radiation

### 3.1 The incident field

Four spectrum types, key `Spectrum type:` (`input_read.f90`), anything else
stops the run:

| type | flag | routine |
|---|---|---|
| `Power-law` | `is_PL_sed` | `J_inc` `J_inc.f90` |
| `Planck` | `spectrum_is_planck` `J_inc.f90` | `planck_stellar_flux_nu/_eV` |
| `Load` | `do_read_sed` | `read_sed` `sed_read.f90` |
| `Monochromatic` | `is_monochr` | `set_energy_vectors.f90` |

(a) The stellar spectral flux at the planet. (b) Power law: the ATES
(Caldiroli et al. 2021) two-segment parameterization, `J = J_norm E^PLind`
normalized separately on the EUV and X-ray segments with the split set by
`Lrapp = 10^(LX-LEUV)` (`J_inc.f90`, `PLind = -1` handled).
Planck: `F_nu = pi B_nu(T_eff) (R_star/a)^2`, cited to Rybicki and Lightman 1979
eqs. 1.51 and 1.13, with a Wien-tail branch above `h nu/kT = 700`.
Loaded: two columns, wavelength in angstrom and `F_lambda` at the
planet (`sed_read.f90`), converted.
(c) The Planck header states the model choice and its limits: line
blanketing and the Balmer jump make it an **upper bound** over 4.8 to 13.6 eV,
while the chromospheric EUV of an active star exceeds it above 13.6 eV by orders
of magnitude.

**One spectrum type builds every band** (decision 13 of the plan, implemented in
Phase 1 batch 2a-PLANCK and 2c-BALMER). `stellar_flux_eV` (`J_inc.f90`)
is the single accessor at arbitrary energy, and `spectrum_covers_eV`
 the coverage predicate. **There is no separate sub-Lyman field**:
the 4.8 to 13.6 eV band is the same spectrum evaluated on the same grid. This
closes the audit finding that the sub-Lyman field was 500 times too faint.

Normalization: `J_XUV = (10^LX + 10^LEUV)/(4 pi a^2)`
(`set_energy_vectors.f90`). For `Planck` the absolute scale comes from
`T_eff, R_star`, not from LX/LEUV, so the integrated grid flux differs from the
nominal `J_XUV`; the setup report prints both. For `Load`, LX and
LEUV are recomputed from the file (`sed_read.f90`).

The dayside geometry factor `dayside_dilution()` (`parameters.f90`,
key `2D approximate method:`) is 1, 0.5 or 0.25 and is applied once: to the XUV
grid at `set_energy_vectors.f90`, to the five FUV bands at
`util_ion_eq.f90`, to the Balmer band at `excited_hydrogen.f90`, to the
stellar Lyman-alpha beam at `lya_rt.f90`. (c) The header
`parameters.f90` records that the shielded-band shell average (0.20 to
0.33 for the Lyman-Werner band over a molecular layer, `docs/e2_lw_geometry.md`
section 5.3) is **not** implemented. Recorded as C31.

### 3.2 The photon grid and its quadrature

For the two analytic types the grid is built from six bands with fixed bin
budgets (`set_energy_vectors.f90`): the H(n=2) band (20 bins), the
sub-Lyman band `[e_abs_low, 13.598]` (89 bins, `num_TR`, the user's decision of
2026-09-05), then 50 bins each over `[13.598, 24.587]`, `[24.587, 54.418]`,
`[54.418, e_mid]` and the X-ray band. Bins are geometrically spaced
(`geometric_band_edges`), `e_v` is the geometric mean of the two
edges, `de_v` the exact width, so `sum(de_v)` equals the band span
exactly.

The energy integral is a **midpoint rule** `sum(integrand * de_v)`:
`PH_heat_H` `util_ion_eq.f90`, `PH_heat_HHe`.

**Threshold handling, analytic path.** `ionization_thresholds_active`
(`set_energy_vectors.f90`) collects every threshold of an active
absorber; `band_sub_boundaries` cuts each band at every threshold
strictly inside it; `bins_by_log_width` shares the budget by
logarithmic width with at least one bin each. **Every threshold is a bin edge**,
so every bin center lies inside one interval of continuity of every cross
section. The header records the audit measurement of the retired
grid (P_HI 1.095, P_HeI 1.016, P_HeII 1.030 times the exact integral) and the
current claim that every metal rate is within 5e-4 of the exact
integral.

**The loaded-SED path is a different construction.** `read_sed` builds its own
grid, and `set_energy_vectors` only flips it. There

```
de_v(1)      = -0.5*(e_v(2) - e_v(1))
de_v(2:Nl-1) = -0.5*(e_v(3:Nl) - e_v(1:Nl-2))        sed_read.f90
de_v(Nl)     = -0.5*(e_v(Nl) - e_v(Nl-1))
```

that is, central differences around table nodes, with **no threshold cutting**.
This is C5, and it matters because the 2026-09-05 decision puts every planet
example on its literature SED.

**Coverage stop.** `photon_grid_floor_eV` (`sed_read.f90`) is the lowest
active absorber threshold: `e_th_HI/4 = 3.3996` eV when the excited-hydrogen
coupling is armed, else `sub_lyman_absorber_threshold_eV`, which
returns the He 2^3S threshold 4.768 eV or the lowest active low-ionization-
potential metal (K I 4.341, Na I 5.139, Ca I 6.113, Mg I 7.646, Fe I 7.902,
Si I 8.152, S I 10.36 eV), else 13.598 eV. A loaded table that stops above that
floor **stops the run**, naming the starved absorbers and the two remedies
(`sed_read.f90`). This is the user's decision that an SED short of the
triplet band always stops.

### 3.3 Columns and optical depths

`src/modules/functions/utilities.f90`, three routines with the identical
recurrence: `calc_column_dens` (H I, He I singlet, He II, He 2^3S),
`calc_column_dens_one` (H2, H2O, OH), `calc_column_dens_metals`
 (every photoionizable metal ion).

```
N(N+Ng) = dr_j(N+Ng) R0 n(N+Ng) opa_pf(N+Ng)
N(j)    = N(j+1) + n(j) dr_j(j) R0 opa_pf(j)
```

(a) The absorbing column between the star and the inner face of cell j.
(b) Model choice: a **radial, plane-parallel beam**. There is no `cos(theta)`,
no slant path and no spherical path-length factor anywhere in the column; the
three-dimensional geometry enters only as the scalar dilution on the flux.
`lyman_werner.f90` states this explicitly: the radial column is the
substellar ray exactly, so the rate is the substellar rate and the dilution
turns it into a shell average. (c) Valid for a beam whose optical depth changes
slowly across a cell; the cell-mean corrections that exist are listed below.
(d) As cited. (e) See C27 (the Lyman-alpha depth uses a different rule) and
C14 (`opa_pf` multiplies absorbers the opacity dispatcher does not cover).

`dr_j(j) = r_edg(j) - r_edg(j-1)` is the cell's own width
(`define_grid.f90`), which is the Phase 1 batch 2c-DRJ fix of
the audit finding that the stored width was the width of cell j+1.

Optical depths are formed from these columns:
`tauE = (s_hi N1 + s_hei N15 + s_heii N2 + s_heiTR NTR + s_h2 NH2col +
sum_k sigma_tab(:,k) Nm_col(j,k)) * 1e-18`, `util_ion_eq.f90`; and the
attenuated flux weight `int_f = F_XUV exp(-tauE)/(1 + a_tau tauE) opa_pf(j)`,
where the `1/(1+a_tau tau)` factor is the ATES `alpha`
correction, zero unless selected (`input_read.f90`).

**H(n=2) is not in `tauE`.** The Balmer continuum removes no photons from the
beam, while He 2^3S and the low-ionization-potential metals in the adjacent
band do. Recorded as C6.

### 3.4 Photoionization and photoheating, absorber by absorber

The generic form implemented in `PH_heat_HHe` (`util_ion_eq.f90`):

```
P_s(j)  = INT F_E e^{-tau}/(1+a tau) f(p) sigma_s(E) / E dE                [s^-1]
h1_s(j) = INT F_E e^{-tau}/(1+a tau) f(p) sigma_s(E) (1 - E_th/E) f_heat dE
heat(j) = sum_s n_s(j) h1_s(j)
```

The heating share is `photoelectron_share(e_th, e) = 1 - e_th/e`
(`util_ion_eq.f90`): **the photoelectron only**, `h nu - I`. The
comparison-only option `Photoelectron heating: full | <fraction>` charges the
whole photon (`parameters.f90`). `f_heat = fhv` is the
secondary-ionization heat fraction, applied only above `E_th + E_sec_ion` with
`E_sec_ion = 30.0` eV (`parameters.f90`). The routine returns
`heat_of_one_X` with the density divided out, and the contraction with
composition happens in exactly one routine,
`photoheating_of_composition`, whose six-column split
`heat_chan` is exact by construction.

| absorber | rate | heating | threshold | cross section |
|---|---|---|---|---|
| H I | `int_1` | `acc_HI` | `e_th_HI = 13.598434599` eV | hydrogenic Kramers/Gaunt `sigma(E,Z,E_th)` `cross_sec.f90` |
| He I (1^1S) | `int_15` | `acc_HeI` | `e_th_HeI = 24.587389` | Verner, Ferland, Korista and Yakovlev 1996 Table 1 (`cross_sec.f90`); the legacy ATES two-term fit under `ates_photoion_rate` |
| He II | `int_2` | `acc_HeII` | `e_th_HeII = 54.417765` | hydrogenic with the measured potential (`opacity_models.f90`) |
| He 2^3S | `int_TR` | `acc_HeTR` | `e_th_HeTR = 4.767775` | two VFKY96 forms plus a log-linear bridge across the Cooper minimum, `cross_sec.f90`, fitted to Norcross (1971) and TOPbase at about 4 per cent |
| 17 metal ions | `int_m` | `h1m_loc` | `mion_ethr(i)` | Verner et al. 1996, dispatcher `cross_sec.f90`, table filled once `set_energy_vectors.f90` |
| H2, total | `int_h2`, `P_H2` | see channels | `e_th_H2 = 15.425927` | `sigma_H2` `cross_sec.f90`: Backx et al. 1976 Table 1 over 15.4 to 18 eV, Samson and Haddad 1994 Table 1 over 18 to 300 eV, Yan, Sadeghpour and Dalgarno 1998 eq. 19 sum-rule tail above 300 eV |

The thresholds are single definitions in `parameters.f90` since Phase 1 batch
2c-THR, and every cross section turns on exactly there.

**H2 channels** (`src/modules/functions/h2_photo_channels.f90`), four mutually
exclusive final states summing identically to `sigma_H2`, built at
:

| channel | reaction | threshold | share |
|---|---|---|---|
| `ICH_M` | H2 + h nu -> H2+ + e | 15.4 eV | `(1-f_di) s_ion` |
| `ICH_S` | H2 + h nu -> H + H+ + e | 18.076 eV | `(1-q_D) f_di s_ion` |
| `ICH_D` | H2 + h nu -> H+ + H+ + 2e | 51.4 eV, a **vertical** threshold | `q_D f_di s_ion` |
| `ICH_N` | H2 + h nu -> H + H | window 33 to 41 eV | `f_n sigma_H2` |

`f_di` is from **Chung, Lee, Masuoka and Samson (1993), J. Chem. Phys. 99, 885,
Table II**, 69 rows over 18.076 to 124 eV, held at its 124 eV value above
(`cross_sec.f90`). `f_n` is from their Table 1 (`h2_photo_channels.f90`,
flagged as an assumption of the source). `q_D` is a **model, not a measurement**
: key `H2 double ionization`, default `chung80`; the detector-event
inversion, with the verbatim Chung et al. p. 886 quote justifying it, is at
 (Phase 1 batch 2b item A8). Each channel is charged its own threshold
in the heating integrand (`util_ion_eq.f90`), with the neutral channel
charged `1 - D0(H2)/E`. A documented approximation: the double
channel makes two electrons sharing the excess but reuses the single-electron
partition, biasing toward too little heating. (e) C23: the 51.4 eV vertical
threshold spent in the heating calculation exceeds the asymptotic chemical
energy of the represented products by 19.725 eV, with no identified recipient.

**H(n=2) Balmer continuum** (`src/modules/radiation/excited_hydrogen.f90`) is a
separate, **decoupled, lagged-explicit** channel (module header):

```
gamma_2 = INT_{3.3996}^{13.598} F_E/(h nu) sigma_2(E) dE   balmer_band_integrals
hpe_2   = the same integrand weighted by h(nu - nu_2)       balmer_band_integrals
sigma_2(E) = 1.4e-17 (e_th_HI_n2/E)^3                      module head, same file
```

hydrogenic continuum cited to Osterbrock and Ferland (2006); the
field is the run's own spectrum type, the 2c-BALMER fix. The 2s/2p
populations solve the two-by-two system of **Christie, Arras and Li (2013),
ApJ 772, 144, eqs. 12-13** (`n2_populations`) with the rate
coefficients of `hydrogen_n2_rates.f90`, a module created explicitly to end the
`G2s`/`G2p` shadowing incident (module header). The feedback is an effective
extra H I photoionization rate `gph_balmer_HI = gamma2 n2tot/n_HI`
injected in `ioniz_eq` (`ionization_equilibrium.f90`), plus `heat_balmer` added
in the same sweep. The quadrature is a 400-point uniform trapezoid on its own
grid (`balmer_band_integrals`, `excited_hydrogen.f90`), which is C6.

**He recombination photons** (`util_ion_eq.f90`, header):
an on-the-spot Draine (2011) eqs. 14.16/14.17 parameterization; the emitted
photon is shared among H I, H2, every metal ion below `E_c` and, above 24.6 eV,
He I, in proportion to `n_s sigma_s(E_c)`, with the escaping fraction dropped.
Four stated validity limits, including that these photoelectrons
bypass the Shull and van Steenberg partition and that `tau_c` uses `dr_j`, that
is a radially escaping photon, whereas the true mean chord is longer by a factor
of order 2. Default on (`parameters.f90`).

### 3.5 Secondary ionization

`src/modules/radiation/electron_energy_degradation.f90`, two sources
deliberately divided (module header):

- **Shull and van Steenberg (1985), ApJ 298, 268**, their eqs. (1)-(2) with
  Table 2 coefficients, supply the ionization budget:
  `svs85_fion_HI(x) = 0.3908 (1-x^0.4092)^1.7592`,
  `svs85_fion_HeI(x) = 0.0554 (1-x^0.4614)^1.6660`.
- **Dalgarno, Yan and Liu (1999), ApJS 125, 237** supply the heating fraction
  (their eq. 14, `eta = 1 + (eta0-1)/(1 + C x^a)`, Table 7), the mean
  energy per ion pair (eq. 13, Table 4), and everything molecular.

Targets receiving secondary ionizations: **H I, He I, H2**. The accumulators
`acc_secHI`, `acc_secHeI`, `acc_secH2` are built for every absorber whose
photoelectron exceeds threshold (`photoionization_field_at_cell_HHe`,
`util_ion_eq.f90`), and the rates are added there; `P_H2` additionally receives
`Psec_H2_di = dissoc_ion_per_H2p Psec_H2`, one proton per 22 H2+, from the
paragraph after Dalgarno et al. eq. (10)
(`electron_energy_degradation.f90`). The branchings return per target
particle, so no rate is divided by a vanishing neutral density.

The heat fraction (`photoelectron_energy_partition`, final line):

```
f_h = 1 - (1 - eta_mix)(1-x) + f_diss + f_vib + f_fluor
f_h = min(max(f_h,0), 1 - f_ion)
```

with `eta_mix` the Dalgarno eq. (12) H/H2 mixture and three molecular
heat returns absent from Dalgarno's own `eta`: dissociation kinetic energy,
direct vibrational excitation of v = 1,2 quenched by
`f_vib_heat`, and B/C fluorescence return.

(c) Stated validity: Dalgarno eq. (14) is stated by its authors for
`x <= 0.1` and is deliberately extrapolated to `x -> 1` in the wind; SvS85 is
fitted over `1e-4 < x < 1` with an `E0 > 100` eV footnote; both fix
`n(HeI)/n(HI) = 0.1`, so a helium-to-hydrogen ratio of order 1 to 10 is an order
of magnitude outside both; neither carries metals, while the ionized fraction
handed in includes metal electrons, which is an explicit extension; Dalgarno
requires `T` below about 2000 K for ground-vibrational H2. The composition
renormalization from the fixed SvS85 ratio onto the cell's own densities is in
the same routine; the H2 collision weight `dal_h2_weight = 1.89` is read off
Dalgarno eqs. (9)/(10) and the identification is explicitly not made in the
paper. The 30 eV floor discards photoelectrons between 13.6 and 30 eV
that can ionize hydrogen, which in the hot-Uranus band carry 18 per cent of the
incident XUV energy, because no consulted source resolves the partition there
(`util_ion_eq.f90`).

**Staging.** Key `Secondary_ionization:` (`input_read.f90`): `False` off;
absent or `True` **staged** (the default, `parameters.f90`); `Immediate`
active from step 0, labeled for comparison only (`parameters.f90`). The
coupling is gated everywhere by
`sec_on = use_sec_ion .and. sec_ion_active` (`util_ion_eq.f90`).
Staging in the marching loop: initial arming `EXHALE_main.f90`; the flip
at once the wind has converged without it, which then **re-arms
every convergence gate**; a hold of `N_stall` steps; and a flip
before the JFNK finish, because GMRES stagnates if handed a
state relaxed without the coupling. The state is written into and read from the
profile header (`utilities.f90`) so that a restart does not re-stage a
coupling the state was converged under.

### 3.6 Lyman-alpha

`src/modules/radiation/lya_rt.f90`. **Closure: the escape probability of the
static plane-parallel damping-wing slab solution** (Neufeld 1990, ApJ 350, 216,
eq. 3.27 at zero continuum destruction; Harrington 1973, MNRAS 162, 43, eq. 40
for the mid-plane source), not Monte Carlo; the justification is
that WASP-121b has `tau0 = 1.2e8`, so an acceleration-free Monte Carlo is
intractable and core skipping is excluded.

```
DnuD     = nu_lya sqrt(2kT/mu)/c                          jlya_escape_prob
a        = (A_2p1s/4 pi)/DnuD                             jlya_escape_prob
tau      = top-down rectangle sum of n_HI C_lya/DnuD dr   lya_line_center_optical_depth
<N>      = (4 sqrt(6)/pi^2)[(tau_up+tau_dn)/2] Phi(xi)    jlya_escape_prob
beta_esc = 1/(1 + <N>)                                    jlya_escape_prob
beta_sob = (1 - e^{-tau_S})/tau_S                         jlya_escape_prob
beta_tot = beta_esc + beta_sob (1 - beta_esc)             jlya_escape_prob
J_int    = (2h nu^3/c^2)(g1s/g2p) P (1-beta_tot)
           / [ (A_2p1s beta_tot + D_2p) n_HI ]            jlya_escape_prob
J_star   = xi F_Lya_star/(4 pi Dnu_star) T_star [...]     jlya_escape_prob
J_lya    = J_int + J_star                                 jlya_escape_prob
```

(a) The Voigt-profile-averaged mean intensity at Lyman-alpha, cited as Huang
et al. 2023 eq. 11. It is **not** a cooling suppression factor and not
an in-situ emissivity. (b) Modes: 0, the default, is the parameterized Huang
et al. (2017) eq. 6 form in `excited_hydrogen.f90`; 1 loads an external
profile (key `Jlya RT file:`); 2 runs `jlya_escape_prob` in line
(key `Jlya escape-prob: True`). `lya_star_boost` is described at the code site as
tuned to Huang et al. 2023 Fig. 11, that is, a model parameter and not a
measured quantity. (c) The slab solution holds in the damping wings of a
static, homogeneous slab, `(a tau)^(1/3) > 10` (Neufeld section Va), which
`lya_wing_domain_record` counts cell by cell; the two faces of the slab share
one escaping population (Neufeld eq. 2.25) and are not two probabilities
multiplied, while the Sobolev channel is a separate closure joined to it as a
union of independent chances. (d) As cited. (e) **C7**: the field feeds only the
2s/2p pump (`excited_hydrogen.f90`), and **there is no Lyman-alpha
cooling suppression anywhere**: `lambda_coex_HI` (`Cool_coeff.f90`) is the
optically thin form, unmodified by `beta_esc` or `tau_lya`, in the same cells
where `beta_tot` is far below one.

The one field genuinely shared with another consumer is
`lya_stellar_beam_transmission`, used both by the n=2 pumping and
by FUV band B2 (`util_ion_eq.f90`), a deliberate single definition
(Phase 1 batch 2a-LYA).

### 3.7 Lyman-Werner and the FUV bands

`src/modules/lower_atmosphere/lyman_werner.f90`. Band 912 to 1110 angstrom, the
Draine and Bertoldi (1996), ApJ 468, 269 interval (module header). **Changed on
2026-09-06 (item LW-NORM-B):** the band is 912 to 1201 angstrom, the interval
the self-shielding table's line list actually covers, and `e_lw_photon_erg` is
`1.88021e-11 erg = 11.7354 eV` accordingly (READ 2026-09-07).

```
k_LW,thin = sigma_LW F_LW / <h nu>
sigma_LW  = zeta_diss(0)/F = 4.17e-11/1.208e7 = 3.452e-18 cm^2   module head
<h nu>    = 2hc/(912+1110 A) = 12.2635 eV               e_lw_photon_erg, same file
```

from DB96 Table 1, Table 2 and the annotation inside their Fig. 1 panel, which
the code site notes is in neither caption. (c) Stated validity:
`sigma_LW` depends on the shape inside the band, and repeating the arithmetic
for DB96's Draine-1978 field gives 2.685e-18, 22 per cent lower, so the flat
`F_lambda` value is good to about 25 per cent for a band that is not strongly
tilted; a band dominated by a single emission line is outside the range.

The live rate is structurally DB96 eq. (40), line self-shielding times
`exp(-tau_cont)`, with no dust, and the version actually called is
the **cell mean** `lyman_werner_dissociation_rate_cell_mean`, a
composite three-point Gauss-Legendre rule on segments of at most 0.05 dex of
column and 0.5 of continuum depth, with an accuracy of 2.5e-4 recorded at the
code site over 1500 random cells; the star-ward face carries
`NH2col(j+1)`. This is Phase 1 batch 2b.

**Self-shielding is a table, not a fit.**
`src/modules/lower_atmosphere/h2_self_shielding_table.f90` is a generated file
(rebuilt by `src/utils/h2_shielding_table_line_by_line.py`) tabulating
`sigma_diss`, `p_single` and `p_eff` on a 39 x 9 x 3 grid in
`(N_H2, T, n_H)`. Line data: **Abgrall, Roueff and Drira (2000), A&AS 141, 297**,
on a 4e5-point frequency grid with Voigt profiles and LTE populations, a
plane-parallel slab, flat `F_lambda` normalized to 343 erg cm^-2 s^-1;
checked against the Meudon PDR code (Le Petit et al. 2006) to 5 to 22
per cent over 1e19 to 1e21 cm^-2. The trapping ratio `p_eff/p_single`
comes from CLOUDY (Ferland et al. 2017) and is plane-parallel, provenance name
`SLAB_SURROGATE`. (c) Stated validity: T from 700 to
3200 K, n_H from 1e12 to 1e14 cm^-3, N_H2 from 1e12 to 5e21 cm^-2; outside,
the edge value is clamped, "a statement that nothing was calculated there, not
that the value holds". The two published fits are retained but **not called by
the rate**, with the reason. DB96 eq. (39) then remained live for a different
quantity, the band share removed from the shared beam. **That use was removed
on 2026-09-06 (item LW-NORM-B):** the share the beam loses is now the table's
own pump absorption, `h2_lw_band_photon_fraction_absorbed` and
`h2_lw_dissociation_per_absorbed_photon`
(`h2_self_shielding_table.f90`), read by `fuv_lw_photon_field`
(`util_ion_eq.f90`) into `p_lw_absorbed`.

The Doppler parameter is thermal only, with no turbulent term.
The fragment heat is 0.4 eV to the H+H pair, cited to Black and Dalgarno (1977),
p. 418.

`src/modules/lower_atmosphere/water_photolysis.f90` adds four more bands (B1
1110 to 1201, B2 the Lyman-alpha line, B3 1231 to 1450, B4 1451 to 2304
angstrom) for reactions O3, O4, O5 (H2O) and O7 (OH), with
`j_b = s_b N_b exp(-tau_b)`, `s_b` the wavelength-weighted mean cross section.
Cross sections are the `photodissociation` datasets of Photochem's
`H2O.h5` and `OH.h5`, concatenated by Photochem from Huebner and Mukherjee
(2015), Heays, Bosman and van Dishoeck (2017), A&A 602, A105, and Ranjan et al.
(2020). (c) Stated approximations: band-averaged `s_b` with
`exp(-tau)` makes deep rates a **lower bound**; the shape
sensitivity is large, and on a real HD 189733 b spectrum `s_b/s_b(flat)` is
1.03, 1.03, 1.35 and 0.22 for H2O on LW, B1, B3, B4, so **B4 is off by a factor
4.6 to 6.1**; a line-dominated M-dwarf B1/B3 is the failure case.
The absorbers deliberately omitted are listed in the module header.

All three LW absorbers share one beam through `fuv_lw_photon_field`
(`util_ion_eq.f90`) using the exact identity
`T_out - T_in = e^{-tau_out}[(1-A_out)(1-e^{-dtau}) + dA e^{-dtau}]`.

### 3.8 The opacity dispatcher

`src/modules/radiation/opacity_models.f90` selects the photoionization cross
section itself, configured by the presence of `opacity.inp`. Four models
(`photoion_sigma`): `A` analytic (the default, reproducing
`cross_sec.f90` exactly), `C` constant, `P` constant plus a Robinson and Catling
(2012) pressure broadening `opacity_pT_factor(p) = 1 + a (p/p_pivot)^n`, and
`T` tabulated with log-log interpolation.

It feeds the **same** optical depth as the photoionization opacity; there is no
second `tau`. The model is applied once, at grid construction
(`set_energy_vectors.f90`), and the `P` multiplier is carried in **both**
the columns (`utilities.f90`) and the local flux weight
(`util_ion_eq.f90`), so a cell removes from the beam exactly the
photons it absorbs. That invariant was restored in Phase 1 batch 2a-OPAC.

(e) **C14**: only H I, He I, He II and He 2^3S go through the dispatcher.
`s_h2` (`set_energy_vectors.f90`) and `sigma_tab` for all 17 metal ions
 call their analytic curves directly, while `opa_pf` still multiplies
their columns and their local absorption.

---

## 4. Species balance

### 4.1 The closure

The default closure is **local chemical and ionization equilibrium, re-solved
every step in every cell**, by `ioniz_eq`
(`ionization_equilibrium.f90`), called from `EXHALE_main.f90` (the startup
equilibration, the `coupled_source` loop of the attempted step,
`equilibrate_loaded_composition`, `physical_handoff_check`,
`update_map_begin_step` and `steady_wind_with_element_diffusion`) and from
`eval_residual` (`steady_newton.f90`). The cell sweep is
parallel over cells, serial on the first step because
that warm start reads the neighboring cell.

### 4.2 The `System_HeH_*` family

**Sizing** (`input_read.f90`):

```
no He                      N_eq = 1
He, no triplet             N_eq = 3
He + triplet               N_eq = 4
metals, no triplet         N_eq = 3 + 2*n_melem
metals + triplet           N_eq = 4 + 2*n_melem
molecular                  N_eq = 7
molecular + triplet        N_eq = 8
 + oxygen chemistry        N_eq = N_eq + 2
 + metals                  N_eq = N_eq + 2*n_melem
```

**Dispatch**, in two places that must agree: the solve
(`ionization_equilibrium.f90`) and
the residual re-evaluation used for acceptance. The row bases
`mbase = metal_row_base()` and `iox = oxygen_row_base()` are computed once and
shared, with the single definitions at
`System_HeH_mol_metals.f90`.

The solver differs by branch: the H-only, HeH and HeH_metals branches use
`solve_ieq`, an analytic-Jacobian damped Newton with a MINPACK `hybrd1` fallback
(`newton_solver.f90`); the TR, TR_metals, mol and mol_metals branches call
`hybrd1` directly, because those merged residuals carry no written Jacobian
(`System_HeH_TR_metals.f90`, `System_HeH_mol_metals.f90`).

**Unknowns.** Stage fractions of one element each:

| routine | N_eq | unknowns |
|---|---|---|
| `ion_system_H` `System_H.f90` | 1 | `x1 = n_HII/n_H` |
| `ion_system_HeH` `System_HeH.f90` | 3 | `+ x2 = n_HeII/n_He`, `x3 = n_HeIII/n_He` |
| `ion_system_HeH_TR` `System_HeH_TR.f90` | 4 | `+ x4 = n_HeITR/n_He` |
| `ion_system_HeH_metals` `System_HeH_metals.f90` | 3+2 n_melem | metals from `x4` |
| `ion_system_HeH_TR_metals` `System_HeH_TR_metals.f90` | 4+2 n_melem | layout |
| `ion_system_HeH_mol` `System_HeH_mol.f90` | 7 (+1) | `x4 = 2n_H2/n_H, x5 = 2n_H2p/n_H, x6 = 3n_H3p/n_H, x7 = n_HeHp/n_H`, triplet `x8`; layout |
| `ion_system_HeH_mol_metals` `System_HeH_mol_metals.f90` | above `+2 +2 n_melem` | oxygen at `oxygen_row_base`, metals at `metal_row_base`; layout |

**Rows.** The shared atomic rows are in `ion_residual_core.f90`: `heh_rows`
 (each of the form `n_X Gamma_photo + (n_X beta_coll - alpha_rec n_X+)
n_e`), the analytic Jacobian pieces, the triplet row `tr_triplet_row`
 (photoionization of the metastable, recombination into it, collisional
excitation `q13` and de-excitation `q31a`, `q31b`, electron-impact ionization,
radiative decay `A31`, and the **total** He(2^3S)+H ionization `Q31 n_HI`,
because both the Penning and the associative branch quench the metastable,
stated), the TR four-row form `heh_tr_rows` including the
Penning proton source, `metal_electron_sum` and
`metal_rows`.

The molecular rows are `mol_heh_rows` (`System_HeH_mol.f90`): H+, He+, He++,
H2, H2+, H3+, HeH+ and He 2^3S. The four H2 photoabsorption channels are
carried as three disjoint subsets of `P_H2` with the remainder leaving H2+,
documented at the routine, citing Yan, Sadeghpour and Dalgarno (1998) section 4 and Dalgarno,
Yan and Liu (1999) after their eq. (10). The oxygen rows are
`oxygen_carrier_rows` in the same file, with
O4 and O6 entering as one channel because O(1D) has exactly one sink in the
audited set and a 1.6e-3 s lifetime, argued in the comment above it.

**Charge exchange** enters two ways that never overlap: the generic Huang
Table-4 assembly `cx_add_to_fvec` pointed at `cx_metal_base`
(`ionization_equilibrium.f90`), and the dedicated He/H
group-B pair `he_h_cx_fvec` (`charge_exchange.f90`), present in every
system containing helium, gated by `he_h_charge_exchange` (default `.true.`).
Group B is excluded from `cx_act` by `cx_init` precisely so it is never
counted twice, and a sign flag reconciles the two row orientations
(`System_HeH.f90` and `System_HeH_TR.f90`).

**Row scaling and imposed rows.** Molecular rows are divided by their own
turnover rate so the helium and molecular blocks reach `hybrd1` with the same
weight (`System_HeH_mol.f90`, scales). A
transported quantity then **replaces** a balance row with an identity row,
applied after the scaling so the row is exactly `x - x_fix`:
`fvec(4) = x(4) - x_h2_fix` and `fvec(1) = x(1) - x_hp_fix`.

**The constrained fallback.** When no starting point yields a root in a
molecular cell, `equilibrium_from_molecular_limit` is promoted
(`ionization_equilibrium.f90`,
`constrained_chemical_equilibrium.f90`). It uses `u_k = ln(n_k)` so every
iterate is positive, the CEA method of Gordon and McBride 1994, NASA RP-1311;
element conservation becomes explicit residual rows;
charge neutrality is identical rather than approximate; and the
radiation field is walked from six decades below the cell's own up to full
strength by natural-parameter continuation, Allgower and Georg 2003.
Only the molecular branch is promoted, because the constrained form is aimed at
species carrying two elements at once, HeH+ and the oxygen carriers, whose two
budgets a fraction layout can satisfy only one at a time
(`ionization_equilibrium.f90`).

### 4.3 Where the rates are defined

| family | routine | source stated at the code site |
|---|---|---|
| H/He radiative and dielectronic recombination | `alphaB_HII_new` etc. `Cool_coeff.f90` | Badnell RR (2023 update of Badnell 2006, ApJS 167, 334) plus Badnell adf09 DR, minus Mao and Kaastra (2016), A&A 587, A84, alpha_1 |
| alternates | `Cool_coeff.f90` | `legacy_hhe_rates` gives Hui and Gnedin 1997; `atomic_rate_set_k22` gives Koskinen et al. 2022 Table 1 R1/R2 |
| He triplet recombination | `alpha_rec_HeII_11S`, `_23S` `Cool_coeff.f90`, summed | Oklopcic and Hirata 2018, ApJ 855, 11; the only H/He recombination fits with a stated fitted range, 5e3 to 2e4 K, **not enforced** |
| metal recombination | `alpha_rr_metal`, `alpha_dr_metal`, 14 Badnell blocks | Ca II from Shull and Van Steenberg 1982 with no DR term; Fe I and Fe II from Huang et al. 2023 eqs. (5)-(6), because the Badnell DR project stops at the P-like sequence |
| collisional ionization | `Cool_coeff.f90` and 17 metal routines | Voronov 1997, ADNDT 65, 1. The metal block carries a standing "VERIFY against Voronov 1997, Table 1" note: the code does not consider that transcription verified |
| C I Voronov row | `Cool_coeff.f90` | corrected in Phase 1 batch 2c-CI (X and K had been in the wrong slots) |
| charge exchange | `charge_exchange.f90` | groups by Huang et al. 2023 Table 4; group B1 Kimura et al. 1993 (tabulated 6e3 to 1e5 K), B2 Stancil et al. 1998 from Zygelman et al. 1989 (**published range 1 to 1000 K**, evaluated from the 1140 K base upward), the non-radiative companion Kingdon and Ferland 1996 |
| Penning He(2^3S)+H | `Cool_coeff.f90` region, `ioniz_HeI23S_H` | Taylor 2025 temperature-dependent fit, errors below 2 per cent over 200 to 10000 K per the code site, evaluated unclamped outside with the accuracy claim withdrawn |
| H(n=2) | `hydrogen_n2_rates.f90` | Christie et al. 2013 Table 2 collisional; Draine 2011 Table 2 R2/R8/R9 |
| molecular network | `mol_rates.f90` | Koskinen et al. 2022 Table 1, R1 to R23 |
| oxygen | `oxygen_rates.f90` | Shomate reverse rates; `co_equilibrium_density` |

(e) Duplicated definitions are collected as C10 through C13.
**He(2^3S)+He(2^3S) is not implemented at all**: no metastable-metastable
channel exists in the tree.

### 4.4 Transport operators

Three optional operators break the local closure. All three are operator split
against the hydrodynamics, and none is folded into `u`.

| operator | key | default | default at | implementation | call site |
|---|---|---|---|---|---|
| He/H element diffusion | `He_diffusion: True` | `.false.` | `parameters.f90` | `element_diffusion_step` (`binary_element_diffusion.f90`) | the attempted step of `EXHALE_main.f90` and `steady_wind_with_element_diffusion` |
| eddy coefficient | `He_Kzz:` | 0.0 | `parameters.f90` | array `kzz_cell` | filled in `init.f90` |
| thermal diffusion | `He_alphaT` | 0.0, a no-op | `parameters.f90` | `settling_coefficient` | |
| ambipolar field | `He_ambipolar: False` | `.true.`, on | `parameters.f90` | `settling_coefficient` | |
| metal trace diffusion | `He_metal_diffusion` | `.false.` | `parameters.f90` | `solve_trace_element_in_hydrogen` | `element_diffusion_step` |
| molecular carrier transport | `Molecular carrier transport` | `.false.`, **but conditional** | `parameters.f90` | `photochemical_transport_step` (`diffusive_photochemistry.f90`) | the attempted step of `EXHALE_main.f90` |
| ionization transport (H+) | `Ionization transport: True` | `.false.` | `parameters.f90` | fifth carrier of the same operator | same |
| coupled carrier solve | `Coupled carrier solve` | `.false.` | `parameters.f90` | `set_transported_species_rows` (`steady_newton.f90`) | the JFNK hand-off of `EXHALE_main.f90` |

**The carrier default is conditional.** If the key is absent, `input_read`
(`input_read.f90`) sets `carrier_transport = thereis_oxychem`, and the reason
is recorded beside it: the oxygen cycle exists to compute a base partition
that a local steady state cannot produce, so it turns transport on; a molecular
run without it keeps a closure that is violated on its own solution across the
H2 front.

**The carrier equation** (`diffusive_photochemistry.f90`):

```
d_t n_i + (1/r^2) d_r [ r^2 (n_i v + Phi_i) ] = P_i - L_i
Phi_i = -n_tot (D_i + K_zz) d_r f_i - n_i D_i [1/H_i - 1/H_atm]
```

Backward Euler, assembled in `carrier_residual`: time term, diffusive and
drift face divergence, advection, chemistry source. **The advective term
changed on 2026-09-07 (item B4-1c):** in the marching loop the carriers are
advected as `F_rho Y_c/m_c` inside the Runge-Kutta stages
(`advected_carrier_register` and `species_advection_stage`,
`binary_element_diffusion.f90`) and the carrier row carries the diffusive half
alone. **Item B5b made the remaining advective term the same one**: where a
row carries it -- the stationary balance and the fixed-wind relaxation -- it
is the divergence of those same face species fluxes
(`carrier_advective_divergence`, `species_flux_divergence`), so one
discretization of the term exists in the code. The cell-velocity upwind
difference and its deferred van Leer correction are gone. Solved by an inner
block-tridiagonal Newton (`solve_carriers`, `block_thomas`) whose advective
entries are the donor-cell linearization of the face-flux term.

Base boundary: the direction is the sign of the **face mass flux at the base
face**, never `v(1)`. With an inflowing face the face composition is
reconstructed from the inner ghosts, which carry the handoff partition where
a handoff states one (`q_H2_base`) and the base cell's own where none does;
face 0 carries no diffusive flux. The proton deliberately gets no stated
partition.

Headroom: `carrier_headroom` is the largest n(H2) the cell's
hydrogen budget admits. `limit_to_element_budget` rescales
carriers back onto the element simplex, including the outer ghosts. It applied
the CO thermal ceiling until that was deleted on 2026-09-06 (item
CEILING-DEL); what is left in the routine is conservation only. A trial state above the headroom is refused by the
coupled solve rather than clamped.

**How transported values reach the chemistry**: they are imposed, not re-solved.
`ioniz_eq` (`ionization_equilibrium.f90`) sets `x_h2_fixed`/`x_ox_fixed` and
`x_hp_fixed`, the last **only for interior cells j >= 1**, because the
two lower ghosts are the inflow reservoir for which the operator states no
ionization fraction. This is C18. CO is taken as transported in the same
routine rather than re-formed from the local state, with the reason stated
there.

**Gates.** `Ionization transport: True` stops the run unless
`Molecular chemistry` (`input_read.f90`) and
`Molecular carrier transport` are on, and is incompatible with
`Coupled carrier solve` and with `Solver: Newton` (`input_read.f90`,
`input_read`).
That last refusal is C19.

### 4.5 Element census

`src/modules/functions/element_census.f90`, 475 lines. **It checks; it does not
correct**: `element_census_verify` declares `rho` and `f_sp` as `intent(in)`.

It computes, for hydrogen, helium and the ten metal elements, the **nucleus**
number density in each cell with every carrier at its exact stoichiometric
multiplicity, and beside it the positive charge density. The multiplicities are
**not written here**: they are read from `bsp_nH`, `bsp_nHe`, `bsp_nO`, `bsp_nC`,
`bsp_charge`, `mion_elem`, `mion_stage` of the species table, so
adding a carrier means adding a row there. Excited levels are skipped.

Two distinct invariants: the **ratio** `n_El/n_H` is gated across every
composition operator including `ioniz_eq`, since that routine rewrites `rho`
from its own `calc_rho` and may move all densities by one common factor;
the **absolute** nucleus density is gated additionally when the
caller states that its operator holds `rho` fixed, which the carrier
write-back and the transport step are. The mass closure is reported always.
Charge neutrality is computed but is a reported diagnostic, not a gated
invariant.

Default gate 1.0e-9, justified as three decades above
the round-off measured on the marching states and a decade below the cell
solver's own `xtol = sqrt(eps) = 1.5e-8`. The module is **off** unless
`EXHALE_ELEMENT_ASSERT` is 1 (report) or 2 (report and stop).
Call sites: `EXHALE_main.f90`;
`ionization_equilibrium.f90`;
`diffusive_photochemistry.f90`.

### 4.6 Binary element diffusion

`src/modules/functions/binary_element_diffusion.f90`, 2357 lines.

**Unknown: the helium mass fraction** `X = rho_He/rho`.
Continuity form in the header, but **solved in advective form**
(residual):

```
rho (X^{n+1} - X^n)/dt + (1/r^2) d_r (r^2 J) + rho v d_r X = 0
```

The stated reason is that the routine receives only the new `rho, v,
T, dt` and not the hydro face mass fluxes, so an independent conservative
advection of `rho X` could not keep a uniform `X` uniform. The consequence
recorded in the code: helium mass is conserved only to the order of the hydro
truncation in a transient, the two forms agreeing exactly at steady state. This
is C28. **Changed on 2026-09-06 (item B4-1):** the marching route carries the
element advection on the hydro face mass fluxes inside the Runge-Kutta stages,
and `solve_mass_fraction` is called with `advect` false, so the advective form
below is what the fixed-wind relaxation runs.

One backward-Euler step per call with coefficients frozen within the step;
finite volume with faces at `r_edg(0:N)`; advection one-sided
upwind on the **cell** velocity, not face-averaged, because face-averaged
upwinding froze a cell whose two face velocities straddle zero and produced the
helium hole at the breathing HD 209458 b base. Nonlinear
in `X`, solved by a Newton iteration with a tridiagonal Jacobian and step
halving. Inner Dirichlet reservoir, outer zero gradient.

**The flux** (coded, face coefficients):

```
J = -rho (D12 + K_zz) dX/dr - rho D12 G X(1-X)
G = [(m_He - m_c1) m_H g - (Zbar_He - Zbar_c1) eE]/kT
    + alpha_T dlnT/dr - dln(psi)/dr        [cm^-1]
```

| term | present | source at the code site |
|---|---|---|
| Fickian gradient | yes | module header |
| gravity, barodiffusion | yes | module header |
| **ambipolar electric field** | **yes, and computed** in `settling_coefficient` | **Koskinen et al. (2013) section 2.1**, quoted |
| thermal diffusion | present, no-op by default | **no source given** |
| eddy `K_zz` | gradient term only | module header |
| `-dln(psi)/dr` | yes | derived in place |

The ambipolar term is not a hard-wired constant:
`eE = -(1/n_e) d(n_e kT)/dr` is evaluated from the solved electron density and
temperature by central difference, and the same array is handed
out (`eEout`) so the trace-metal solve uses one definition;
the comment records that the metal loop previously assumed the
proton-plasma constant `eE = m_H g/2`, "which is wrong wherever helium or the
metals carry a significant share of the electrons". Stated limits.

Peclet hybrid; the two factors of `X(1-X)` are taken from
**opposite sides of the face** so the drift flux vanishes at both
ends of the composition axis, which is what bounds `X` in [0,1].

**Coefficients**, all Chapman-Enskog first approximation
`D_st = 3kT/(16 n mu_st Omega^{1,1})` with `n` the collision
partner density:

| limit | routine | source at the code site | note |
|---|---|---|---|
| hard sphere | `hard_sphere_pair_diffusion` | Banks and Kockarts (1973), Aeronomy Part B | no table or equation number given |
| polarization | `polarization_pair_diffusion` | Schunk and Nagy, Ionospheres, eq. 4.88 | the book is **not** in `references/`, so the constant was re-derived from the collision integral |
| ion-neutral | `ion_neutral_pair_diffusion` | the two collision integrals add, so frictions add | |
| Coulomb | `coulomb_pair_diffusion` | Paquette et al. 1986, ApJS 61, 177, section II; Schunk and Nagy eq. 4.142 | neither in `references/`, constant re-derived |
| Coulomb logarithm | `coulomb_logarithm` | Spitzer 1962, section 5.2 | floored at 1 |

Mixture rule: **Blanc's law**, stated valid exactly for a trace
species and as the standard approximation otherwise.
Polarizabilities from Schwerdtfeger and Nagle (2019), Mol. Phys. 117, 1200.

**Zero net diffusive mass flux.** For the H/He binary it is enforced
**structurally, by construction**, not by an algebraic subtraction: the gas is
split into exactly two components whose masses sum to `rho` exactly
(module header), so `X_He = 1 - X_1` identically and only one
diffusive flux is independent; the code never forms a second flux.
`project_elements` then meets both element
totals exactly, depositing each element's shortfall
into its neutral ground stage. The closure is measured and
reported, not assumed. **It is not enforced for the trace metals**
(stated), which is C29. HeH+ is in
neither carrier list and neither mean charge,
and metal nuclei are outside `n` and the reduced mass.

**Coulomb suppression is not a separate factor**: it is produced by the
stage-resolved pair sum, in which the He++/H+ Coulomb pair dominates the Blanc
harmonic sum in the ionized region. The physical statement and its size are at
: the Coulomb momentum-transfer cross section at about 1e4 K is
roughly 1e-12 cm^2 against 1e-15 cm^2 for a hard sphere, so a single neutral
coefficient would let helium and the metals settle two to three orders of
magnitude too fast. It is exposed as the diagnostic ratio `D_eff/D_neutral`.

`K_zz` enters the **gradient coefficient only, never the drift**,
because eddy mixing transports the mixture as a whole and has no preferred
species. The operator reads only the array `kzz_cell`, filled once
by `eddy_diffusion_on_grid` (`lower_atmosphere_profile.f90`) from
`init.f90`: the scalar `he_kzz` in every cell with no profile, or the
profile's `Kzz` column interpolated onto the grid and
clamped outside.

**The transported pair is H and He as elements**, not species: the chemistry
redistributes hydrogen among H, H+, H2, H2+, H3+ and helium among its stages
without changing either element's mass, so the same equation holds above and
below the molecular front. Fifteen stage pairs are dispatched to
hard sphere, ion-neutral or Coulomb and combined by Blanc's law
(`stage_mixture_diffusion`, `carrier_fractions`). Three documented exclusions:
He 2^3S is not a separate carrier, HeH+ is in neither list, and resonant
H+ + H is not a term because
both partners carry the same element.

The steady route runs this operator in an outer **Picard alternation** because
the steady residual contains no diffusion (`EXHALE_main.f90`), with
under-relaxation starting at 0.5 and halved to a floor of 0.125 on any pass
whose drift failed to fall. `relax_element_composition`
 does not step at the hydro CFL: it builds a composition time scale
per cell, grows it geometrically, and advects on the **face mass fluxes**
`F_rho` of the state, as `div(F_rho X)` (item B5c). It used to carry the
non-conservative `rho v dX/dr` on a smoothed steady mass flux
`mdot/(4 pi r^2)`, because below about 1.02 R_p the converged states do not
satisfy continuity and relaxing that form on `rho_j v_j` converged to "the
composition of a flow that neither conserves mass nor exists"; the
conservative term needs no such repair, and its fixed point is the zero of the
stationary elemental row that judges the pass.

### 4.7 Independent variables and constraints

**In the fraction systems, the constraints eliminate unknowns.**

- **Element totals eliminate the neutral stage**: hydrogen
  `n_HI = (1 - x1 - x4 - x5 - x6 - x7) n_H` (`System_HeH_mol.f90`), atomic
  form `System_HeH.f90`; helium `n_HeI = (1 - x2 - x3) n_He - n_HeHp`
  (`System_HeH_mol.f90`), with the ground singlet a further difference;
  metals `nm0(e) = (1 - x(ix) - x(ix+1)) n_X`
  (`ion_residual_core.f90`); free atomic O closes the oxygen budget the way
  atomic H closes the hydrogen one (`System_HeH_mol.f90`).
- **Charge neutrality eliminates `n_e`**, which is never an unknown and never a
  row: `System_HeH_mol.f90`, atomic forms `System_HeH.f90`,
  `System_HeH_TR.f90`, `System_HeH_TR_metals.f90`, with the metal charges
  added by `metal_electron_sum`. The triplet is neutral and contributes nothing
  (`System_HeH_TR_metals.f90`).
- **Total mass is not an equation of the cell solve.** `rho` is recomputed after
  the solve from the accepted composition (`calc_rho`,
  `ionization_equilibrium.f90`), and the fractions are written back as
  `f_sp = n_species/n_io`. This is exactly why the element census
  gates the ratio always and the absolute densities only under a stated fixed
  density.

So the count is `N_eq` rows against `N_eq` stage fractions, with the totals and
neutrality absorbed. The physical constraints re-enter as **bounds checked after
the solve**: `ionization_fractions_physical` requires every
fraction non-negative and each element's tracked stages summing to at most one,
with the reason that the residuals are polynomials admitting
roots outside the simplex, and such a root gives a negative neutral density that
poisons the photoionization integrals however small its residual;
`clamp_fractions_to_element_budget` is applied only when the
violation is at or below 1e-6, that is round-off of a root on a face of the
allowed region.

**In the constrained formulation the roles invert**
(`constrained_chemical_equilibrium.f90`): the unknowns are `ln(n_k)`,
element conservation becomes explicit residual rows
(`element_conservation_rows`), and charge neutrality stays an
identity holding at every iterate. The motivation is that HeH+ and
the oxygen carriers each appear in two conservation rows at once, which a
post-solve rescaling of stage fractions cannot satisfy.

---

## 5. Energy sources and sinks

### 5.1 The three slots, up front

| slot | code location | what enters |
|---|---|---|
| **(A) RK3 right-hand side** | `EXHALE_main.f90`; source `Source.f90` | **nothing radiative**. `S(1) = 0`, `S(2)` gravity plus the geometric pressure term for PLM, **`S(3) = 0`** (`Source.f90`). Adiabatic expansion cooling is carried by the flux, not by a channel |
| **(B) operator-split source stage** | `EXHALE_main.f90`, then | `solve_energy_semi_implicit(u,W,dt_loc,heat,cool,f_sp,count)` (default, `parameters.f90`) or the explicit `u(3,:) += dt_loc*(heat-cool)`; then `viscous_conduction_step` |
| **(C) steady residual** | `steady_residual.f90` | `R(3,:) = dF(3,:) - S(3,:) - (heat - cool) - Sene`, the same arrays and the same `viscous_conduction_sources` |
| **(D) diagnostic only** | `util_ion_eq.f90`; `write_output.f90` | `Cooling_breakdown.txt`, `Heating_breakdown.txt`, `Lyman_Werner.txt`, the FUV ledger |
| **(E) post-process only** | `T_equation.f90` from `post_process_adv.f90` | the `*_adv` profiles |

### 5.2 Heating channels

All are assembled into one array in `ioniz_eq` (`ionization_equilibrium.f90`),
exported as `heat_out = heat/q0`. Every row below is a term of that one sum.

| channel | note |
|---|---|
| photoheating of H I, He I, He II, He 2^3S, H2, metals | section 3.4; the secondary-electron partition is a **reduction of the integrand** (`util_ion_eq.f90`), not a separate additive term |
| H(n=2) Balmer photoelectric plus Lyman-alpha de-excitation | `heat_balmer = Hpe + Hdx`, `excited_hydrogen.f90`; `Hdx` returns 10.2 eV per collisional de-excitation, gated by `incl_deexc_heat`, default on (`parameters.f90`). **Lagged one outer pass** |
| He recombination coupling | `util_ion_eq.f90` |
| Penning He(2^3S)+H | 6.2 eV per event; `f_penning_HeI23S = 0.9` (`Cool_coeff.f90`) |
| Penning He(2^3S)+H2 | 4.394 eV per event |
| Lyman-Werner fragment kinetic energy | 0.4 eV, Black and Dalgarno 1977 p. 418; the 4.48 eV bond energy is paid by the photon |
| Lyman-Werner fluorescence | landing energy 2.060 eV at 700 K to 2.128 eV at 3200 K from Abgrall et al. (2000) A-values (`h2_vibrational_relaxation.f90`, clamped), replacing the 2.0 eV adopted by Burton, Hollenbach and Tielens (1990); heat fraction `q/(q+A10)` from Hollenbach and McKee (1979) p. 586, with He omitted because no He collider rate exists in either source |
| molecular reaction heat, R5 to R23 | see below |
| FUV photolysis of H2O and OH | the `heat_fuv` sum, section 3.7 |

**Molecular reaction heat**
(`src/modules/lower_atmosphere/molecular_reaction_heat.f90`) is built
from **one** species enthalpy table (`species_enthalpies`) with zero
at H + He + e at rest, so the heat released by `A+B -> C+D` is
`h(A)+h(B)-h(C)-h(D)`, positive meaning heating. Reactions summed:
R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15 (with the Hollenbach and McKee
1979 collisional branching), R16, R17, R18, R19, R20, R23. Sources:
IP(H) and IP(He) from NIST ASD; IP(H2) from the NIST WebBook, 15.42593 eV;
D0(H2) from Huber and Herzberg 1979 through `mol_rates`; D0(H3+) from Mizus
et al. 2019, Mol. Phys. 117, 1663; D0(HeH+) from Coxon and
Hajigeorgiou 1999, Tables 9 and 8. **Deliberately excluded** to
avoid double counting: all photon-driven reactions, radiative
recombination R1/R2, collisional ionization R3/R4, Penning, and H/He charge
exchange. The He(2^3S) associative branch to HeH+ (8.1 eV) is **not deposited
anywhere**, stated at the code site. Default on (`parameters.f90`).

Viscous and conduction heating are their own operator-split stage and their own
residual term, not part of `heat` (section 2.5).

### 5.3 Cooling channels

Assembled by `eval_cool` (`util_ion_eq.f90`, grid-level line transfer
and parallel dispatch) and `eval_cool_cells`, which forms the total:

```
cool = ne*(brem + coex + reco + coio) + cool_M + cool_H3p + cool_H2 + cool_H2O + cool_CO
```

| channel | formula and source at the code site | stated validity |
|---|---|---|
| recombination | `lambda_rec_HII = 3.435e-30 T x^1.970/(1+(x/2.250)^0.376)^3.720`, Hui and Gnedin 1997 (`Cool_coeff.f90`); `lambda_rec_HeII = kT alpha_rec_HeII_total(T)`, deliberately tied to the coefficient the balance uses; `lambda_rec_HeIII` | the He II tie is documented; C11 |
| collisional ionization | `e_th_X_erg * a_ion_X * n_X`; Voronov 1997 coefficients, with the energy removed taken as the **measured** potential and not Voronov's fit `dE` (`Cool_coeff.f90`) | |
| bremsstrahlung | `1.426e-27 sqrt(T) sum n_ion Z^2 gbar(Z)`; Gaunt factor from van Hoof et al. 2014, MNRAS 444, 420, **clamped to the table** (`Cool_coeff.f90`) | molecular ions excluded from the charge sum with the reason; C30 |
| collisional excitation, Lyman-alpha | `lambda_coex_HI = 7.5e-19/(1+sqrt(T/1e5)) exp(-118348/T)` (`Cool_coeff.f90`), a Cen 1992-type form, **no citation on that line** | optically thin, no escape probability; C7 |
| He I, He II excitation | `lambda_coex_HeI` and `lambda_coex_HeII` (`Cool_coeff.f90`), assembled in the same place as the Lyman-alpha row | |
| He 2^3S channels | 10830 line `1.16e-20 sqrt(T) exp(-13179/T)` (`Cool_coeff.f90`), plus 2^3S to 2^1S (0.80 eV) and 2^3S to 2^1P (1.40 eV) with Bray (2000) collision strengths | ground-to-triplet `q13` excluded to avoid double counting |
| metal lines, C/N/O | `Cool_coeff.f90`: `Lambda = T^{-1/2} sum A_i exp(-T_i/T)`, fitted to CHIANTI v11.0.2, script `cooling_data/fit_cno_formulas.py`, maximum error 0.2 to 2.1 per cent | **1e3 to 1e5 K** |
| Mg I, Mg II, Ca II, Na I | exact two-level skeleton with Burgess-Tully type-1 upsilon; maximum error 2.8, 1.1, 1.2, 0.01 per cent | 1e3 to 1e5 K |
| Fe II | `cool_FeII_func` (`Cool_coeff.f90`), **overridden** in `eval_cool_cells` (`util_ion_eq.f90`) by `cool_FeII_ne`: the coronal form overestimates the dense base by about 1e4, so the density-dependent multilevel coefficient is used instead | |
| ground-term fine structure C I, C II, N II, O I | `util_ion_eq.f90`: exact statistical-equilibrium solution with electron **and** hydrogen colliders, line trapping `beta_fs` and the incident base infrared field `nbar_fs` (`Cool_coeff.f90`) | the line transfer covers **only** those eight lines; every other metal ion keeps `beta = 1` |
| coronal cutoff | `Cool_coeff.f90`: every CHIANTI-derived coefficient multiplied by `exp(-x^2)`, `x = (1000 K/T - 1)/w`, exactly 1 above 1000 K, default `w = 0.1` | rationale; the legacy branch is deliberately unguarded |
| H3+ infrared | `util_ion_eq.f90`: `n_H3+ 4 pi E_LTE(T) s(T, n_H2) 1e7`, Miller, Stallard, Tennyson and Melin (2013), J. Phys. Chem. A 117, 9770, Table 5 piecewise over 30-300, 300-800, 800-1800, 1800-5000 K and Table 6 for the non-LTE factor (`h3p_cooling.f90`); net form with the base field | C24 |
| H2, H2O, CO infrared bands | net `n_X (emission(T) - W absorption(T))`, `molecular_infrared_cooling.f90`, equations; H2O/CO band means from the Photochem correlated-k coefficients (Wogan et al. 2025, PSJ 6, 256) through HELIOS-K from HITEMP; H2 lines from Roueff et al. 2019 | stated: optically thin with no escape probability, LTE populations, H2O/CO tables span **50 to 2000 K** and are clamped outside, the H2 line sum unbounded in T, the lower atmosphere taken black. Default off (`parameters.f90`) |

Adiabatic expansion is **not** a channel here: it is the pressure-work part of
`dF(3,:)`.

### 5.4 The temperature equation

**`T_equation.f90` is not used by the marching loop or by the Newton solver.**
Its only callers are `post_process_adv.f90` (`solve_T_brent`)
(`hybrd1`).

It solves the **steady advected energy balance in one cell**, upwind
differenced: in the caloric branch,
`mu_- rho v e(x_H2, x) - mu_+ rho v e(x_H2, T_old) - (coeff x + mu_+ mu_- dr
(heat_old - cool(x))) = 0`, and correspondingly in the monatomic branch. It is a
residual **including** the advection and pressure-work terms, not a pure
heating-equals-cooling balance and not a time-implicit update. The unknown is
the single scalar `x(1) = T/T0`. The root finder is an 80-point
log-spaced upward scan for the first sign change, then Brent,
because the metal-cooled residual is non-monotone and carries a spurious hot
root that `hybrd1` can land on.

`n_e` there is `n_HII + n_HeII + 2 n_HeIII` plus the metal stages, with
**molecular ions deliberately omitted**, and the cooling assembled at
 covers recombination, collisional ionization, bremsstrahlung,
excitation and the metals with the same Fe II and C I/C II/N II/O I overrides,
but **not** H3+ and **not** the molecular infrared bands. `heat_old` is `theat`
from `post_process_adv.f90`, which is photoheating plus the He
recombination coupling plus the He(2^3S)+H Penning heat only: **no Balmer, no
Lyman-Werner, no molecular chemical, no FUV heating**. All of that is C8.

### 5.5 The semi-implicit energy update

`src/modules/time_step/energy_semi_implicit.f90`, `solve_energy_semi_implicit`,
called once each marching step after the ionization solve. The equation
is `e(T_new) - e(T_old) = dt (heat - cool(T_new))` divided by
`(n_tot + n_e)`.

**Replaced on 2026-09-06 (item B2, `docs/Update_EXHALE_stage2.md`).** What this
section describes is the update as it stood on 2026-09-06. The solve is now
bracketed and residual-controlled (`energy_balance_init`,
`energy_balance_update`, `energy_balance_finish` in the same file), cooling is
evaluated at the returned temperature, and reaching the lower bracket end is a
`FLOOR` failure that stops the run rather than a clamped state the run adopts.
The paragraphs below are kept as the dated record C3 rests on.

- `T_old = p/(n_tot + n_e)`.
- `cool_trial = cool`: the incoming cooling, which `ioniz_eq` evaluated
  at the **old** temperature and at the **post-sweep** composition
  (`ionization_equilibrium.f90`).
- `dC/dT` from one perturbed call at `T_old + delta`.
- Exactly **two** Newton iterations, with `dC/dT` held constant and
  the damping derivative `dF_dT = 1 + c_factor |dC_dT|` using the
  absolute value rather than `max(0, dC_dT)`, whose reason (an undamped explicit
  step on the falling cooling branch, the source of a base and front limit cycle)
  is recorded at the statement.
- Cooling is re-evaluated **once**, at the end of iteration 1.
- The composition `f_sp` is frozen for the whole routine and
  `heat` is held fixed across the update.
- A temperature floor `T_trial < 0.01` is clamped and **counted**; a
  floored cell has not converged to a physical temperature.
- On exit `cool = cool_trial`, `p = (n_tot+n_e) T_trial`, `u(3)` rebuilt.

(e) **C3**: the second iteration changes the temperature again after the last
cooling evaluation, so the returned `cool` belongs to the previous iterate unless
the final correction vanishes, and there is no residual acceptance test (audit
N5). `heat` is never re-evaluated, although in molecular gas it contains
temperature-dependent collisional reaction terms, including endothermic
dissociation.

### 5.6 The three-way split, stated precisely

**In the marching source stage (B) and in the steady residual (C), identically**:
every heating row of section 5.2 and every cooling row of section 5.3, plus the
viscous and conduction sources through one routine.

**Only inside the semi-implicit update**: nothing physically new; it re-evaluates
the same `eval_cool` at a perturbed and at an updated temperature.

**Only inside `T_equation` (the `*_adv` post-process)**: the same H/He and metal
channels, **missing** H3+, the molecular infrared bands and the molecular
electron donors on the cooling side, and **missing** Balmer, Lyman-Werner,
molecular chemical and FUV heat on the heating side. The two omissions are
internally consistent, because `post_process_adv.f90` calls
`eval_cool` without `nmol` and `nox`.

**Diagnostic only, never fed back**: `Cooling_breakdown.txt`
(`util_ion_eq.f90`, which calls the same `eval_cool` with `cool_chan`
and prints `max|sum(channels)/cool - 1|`, and adds the
Planck-mean optical depths of H2O and CO measured nowhere else);
`Heating_breakdown.txt` (which **reconstructs** the heating channel
by channel); `Lyman_Werner.txt` (`write_output.f90`); the FUV ledger.

---

## 6. The operator split and the time levels

This section is the heart of D0's item 8 and is where the audit's "residual is
not the update map" finding lives.

### 6.1 The marching step, in order

Loop head `EXHALE_main.f90` ("Main temporal loop"). Let **n** be the state at
step entry. The tag column is the operation index of `attempted_step.f90`,
which is how b1 section 7.1 and the rejection reports name the same operations.
**Rows 15 to 17 were removed on 2026-09-06 (item B3c):** rows 14 to 18 are now
one local source step per cell at fixed volume, iterated to a fixed point, with
no `comp_p_from_T` and no `W_to_U` between the sweep and the energy update.

| # | operation tag | operator | reads | writes |
|---|---|---|---|---|
| 0 | `as_op_dt` | `eval_dt(W,dt,dt_loc)` | `W` at level **n** (written of the previous step) | `dt, dt_loc` |
| 1 | `as_op_checkpoint` | `u_old = u` | `u^n` | `u_old` |
| 2 | `as_op_hydro` | `reconstruction_continuation_rhs(u,...)` = Reconstruct + `RK_rhs` | `u^n` | `WL, WR, dF, S`, face fluxes |
| 3 | `as_op_hydro` | RK stage 1 `u1 = u - dt_loc (dF - S)` | `u^n` | `u1` |
| 4 | `as_op_hydro` | `Apply_BC(u1)` | `u1` interior | `u1` ghosts |
| 5 | `as_op_hydro` | positivity test, repair, else halve `dt` and retry | `u1` | `u1, dt` |
| 6 | `as_op_hydro` | RK stage 2 `u2 = (3u + u1 - dt_loc(dF-S))/4`, BC, guards | `u1`, `u^n` | `u2` |
| 7 | `as_op_hydro` | RK stage 3 `u = (u + 2(u2 - dt_loc(dF-S)))/3`, BC, guards | `u2`, `u^n` | **post-hydro `u`** |
| 8 | `as_op_primitives` | `U_to_W(u,W)`; `rho, v, p` | post-hydro `u` | `W, rho, v, p` |
| 9 | `as_op_primitives` | `get_species_densities` | **post-hydro `rho`**, `f_sp^n` | `nhi ... n_tot` |
| 10 | `as_op_primitives` | `comp_T_from_p` | **post-hydro `p`**, pre-sweep `n_tot, ne` | `T` |
| 11 | `as_op_diffusion` | `element_diffusion_step` | post-hydro `rho, v, T`; `f_sp^n` | `f_sp`, element split |
| 12 | `as_op_carriers` | `photochemical_transport_step` | post-hydro `rho, v`; `f_sp` after 11 | `f_sp` carrier columns |
| 13 | `as_op_excited_H` | `excited_H_update` | post-hydro `T, rho, v`; post-transport `f_sp` | `gph_balmer_HI, heat_balmer`, the n=2 diagnostics |
| **14** | `as_op_ioniz_eq` | **`ioniz_eq(T,rho,f_sp,heat,cool,eta)`** | **post-hydro `T` and `rho`, post-transport `f_sp`** | `f_sp` (all stages), **`rho` in place**, `heat, cool, eta`, `nmol_eq, nox_eq` |
| 15 | `as_op_composition` | `get_species_densities` | post-sweep `rho, f_sp` | `nhi ... n_tot` |
| 16 | `as_op_composition` | `comp_p_from_T` | **the same `T` the sweep was given**, post-sweep `n_tot, ne` | `p` |
| 17 | `as_op_composition` | `W(1,:)=rho; W(3,:)=p; W_to_U(W,u)` | post-sweep `rho, p` | **`u`: the chemistry stage changes the conserved mass and energy** |
| 18 | `as_op_energy` | `solve_energy_semi_implicit` (or the explicit update beside it) | `u,W` after 17; `heat` frozen at the sweep; **post-sweep `f_sp`** | `u(3,:)`, `W(3,:)`, `cool` |
| 19 | `as_op_bc` | `Apply_BC(u)` | | ghosts |
| 20 | `as_op_conduction` | `viscous_conduction_step` with its own `U_to_W` and `comp_T_from_p` at the **pre-energy-step** `n_tot, ne` | post-energy `u` | `u` |
| 21 | `as_op_shapiro` | Shapiro filter if armed, then BC | `u` | `u` |
| 22 | after `as_op_boundary` | `U_to_W(u,W)`; `rho, v, p, E` | **`u^{n+1}`** | `W, rho, v, p, E` |
| 23 | after `as_op_boundary` | `check_base_inflow_is_subsonic` | `W` | nothing |
| 24 | after `as_op_boundary` | `get_species_densities`, `comp_T_from_p` | level n+1 | `T^{n+1}` |
| 25 | after `as_op_boundary` | `mom = rho v r^2`; `du = mass_flux_spread(mom)` | level n+1 | `du` |
| 26 | after `as_op_boundary` | `dtu` over rows and cells `[j_min:N]` | `u^{n+1}, u_old` | `dtu` |
| 27 | after `as_op_boundary` | every `N_resid` steps: `assemble_residual`, `residual_norms`, `steady_gates_met` | `u^{n+1}`, this step's `heat, cool` | `Rres`, gates |
| 28 | after `as_op_boundary` | continuation, PLM to WENO3 hand-off, stop flags | | |
| 29 | after `as_op_boundary` | staged secondary-ionization flip; every gate re-armed | | `sec_ion_active` |
| 30 | after `as_op_boundary` | JFNK hand-off when armed | `u^{n+1}, f_sp` | `u, f_sp, W, ..., heat, cool` |

So the composition is the Lie splitting

```
hydro (SSP-RK3, BC after every stage)
  -> element diffusion -> carrier transport -> excited hydrogen
  -> ionization equilibrium and the radiative rates
  -> energy update -> BC -> viscosity and conduction (Crank-Nicolson)
  -> BC -> Shapiro -> BC -> diagnostics
```

### 6.2 The two time-level facts that matter

**The ionization sweep sees the post-hydro density and temperature**, and the
post-transport composition. `U_to_W` converts the state the three RK
stages left; `rho` comes from that; `T` from `comp_T_from_p` at
 on the post-hydro pressure; `ioniz_eq` is handed exactly
those. This is the stated design: the sweep solves the remaining
stages against the transported partition rather than recomputing it.

**The energy update sees the post-ionization composition** but a cooling
evaluated at the old temperature and a heating frozen at the sweep
(section 5.5).

### 6.3 The composition projection changes energy with no source (C1)

Between the ionization sweep and the energy update, the code did the
following: it solved the composition **at fixed temperature**; it recomputed
the particle count; it rebuilt the pressure from the **same** temperature and
the **new** particle count; and it called `W_to_U`, which turned that
pressure into `u(3)`. **Removed 2026-09-06 (item B3c).** In the monatomic
limit,

```
u_th = (3/2) n_part k_B T ,   Delta u_projection = (3/2) k_B T (n_part,new - n_part,old)
```

and this contribution is not obtained by integrating any physical heating term.
The molecular equation of state adds a change in the H2 rovibrational reservoir
as well. Nuclear conservation removes neither. The audit
(`docs/physics_numerics_audit_20260905.md` section N1, priority P0) demonstrates
the arithmetic through the production conversion routine and states the physical
point: creating a thermal particle does not provide its equilibrium kinetic
energy for free, and dissociation changes particle count and binding energy
together, so it cannot be represented by a temperature-preserving composition
reset plus an independently applied heat source without a derivation.

### 6.4 The chemistry stage is not conservative in `u` (C2)

`ioniz_eq` declared its density argument `intent(inout)`
(`ionization_equilibrium.f90`) and overwrote it with `calc_rho` of the
equilibrium composition; `W(1,:) = rho` in the marching source stage of
`EXHALE_main.f90` and `W_to_U` pushed that into the **conserved mass row**.
**Removed 2026-09-06 (item B3c, `docs/Update_EXHALE_stage2.md`):** the density is
`intent(in)` and `calc_rho` builds the optional `rho_recon`, which is a check
and is never applied. The code's own
diagnostic names this: `EXHALE_UPDATE_MAP` reports a `chem` column and states
explicitly that a non-zero `chem`, `transport` or `filter` column is a term of
the production map that the steady residual does not contain
(`EXHALE_main.f90`).

### 6.5 The written state is not one state (C4)

`src/tests/grid_and_gates/output_state_consistency.sh` compares
`Hydro_ioniz.txt` column 6 (`heat`) against `Heating_breakdown.txt` column 4
(`heat_total`) over the physical cells, reference 0, tolerance 1e-6.
`docs/Update_EXHALE_stage2.md` section 6 (item 2c-STATE)
records the test as RED at 1.2e-5, with the cause stated there: the heat column
is assembled from the pre-sweep rates while the breakdown re-evaluates on the
post-sweep composition, together with the `excited_H_update` call between the
loop and the final write. That number is quoted from the update log and was not
re-measured here.

**A second, separate finding:** the test's own source citation is stale. Its
header and `src/tests/grid_and_gates/README.md:120-145` describe a
block at `EXHALE_main.f90` that forced `sec_ion_active = .true.` after
the marching loop. At the current tree there is no such block: the only two
`sec_ion_active = .true.` assignments are at `EXHALE_main.f90`,
both **inside** the loop, and the post-loop comment states the
opposite policy ("the coupling is left exactly as the loop left it"), which is
the 2c-STATE change. The code the test cites no longer exists; whether the test
now passes was not measured here.

The in-code counterpart is `assert_written_state_is_the_accepted_one`
(`EXHALE_main.f90`, called), which re-measures
the flux spread on the state as written and warns if it differs from the last
steady solve's accepted value by more than 1e-10 relative; it returns
immediately when no steady solve ran.

### 6.6 Temporal order of the split

RK3 on the homogeneous hydrodynamic stages does not make the Lie-split
composition and energy scheme third order in time; backward Euler on the sources
is first order. The audit states this (section 6.4) and notes the actual temporal
order has not been established with a coupled reference solution. Nothing in the
code claims otherwise, and no measurement of the split's order exists in the
tree. Recorded as an open item, not a defect.

---

## 7. The steady solver

### 7.1 What the residual is

`src/modules/time_step/steady_residual.f90`:

```
R(:,1) = dF - S                                                   (mass)
R(:,2) = dF - S - F_mu                                            (momentum)
R(:,3) = dF_E - S_E - (heat - cool) - (w F_mu + q_mu + conduction) (energy)
```

`assemble_residual` calls the **same**
`reconstruction_continuation_rhs` the marching stages use, then, only
when `transport_active()`, subtracts the **same** `viscous_conduction_sources`
the marching loop relaxes, so that the marching fixed point is the
residual zero for those terms.

**Three hydrodynamic rows, plus one H2 row when the carrier transport is coupled (`nvar_jac >= 4`, sections 7.2 and 7.3). There are no other species rows and no ionization rows.**
Ionization is a nested inner equilibrium: the caller supplies `heat` and `cool`
by calling `ioniz_eq` on the state first.

### 7.2 What the unknown vector carries

The conserved hydrodynamic variables on the **physical cells only**, `j = 1..N`;
the ghosts are not unknowns and are refilled by `Apply_BC` on every residual
evaluation. Flat ordering `Y(nvar (j-1) + k) = u(k,j)` (`pack_U`
`steady_newton.f90`).

- Default `nvar_jac = 3`, that is **3N unknowns**.
- With `Molecular carrier transport` and `carrier_in_newton`,
  `set_transported_species_rows(.true.)` makes `nvar_jac = 4` and the fourth
  unknown is **n(H2)** per cell, packed in code density units and not cm^-3
  (`pack_carrier`, with the reason: the Krylov step size
  would otherwise be set by a component thirteen decades away).

**What is dropped or held relative to the marching state** (this is C17):

- the whole composition `f_sp` is **eliminated locally**, re-solved by `ioniz_eq`
  inside each residual evaluation (the routine header, "method A");
- element (He/H) diffusion is **not in the residual at all** and is co-converged
  by an outer damped Picard loop (`EXHALE_main.f90`);
- with `nvar_jac = 3` the carrier partition is **frozen**: `x_h2_fixed` pins
  `x_h2` at whatever the last transport step wrote
  (`ionization_equilibrium.f90`);
- `x_hp_fixed` pins the H/H+ partition from the transported proton, but **only
  for `j >= 1`**, the ghosts re-solving locally (C18); the H+
  balance row is then replaced by `fvec(1) = x(1) - x_hp_fix`
  (`System_HeH_mol.f90`, `System_HeH_mol_metals.f90`);
- `set_base_fix(Y, nfix)` (`steady_newton.f90`) can anchor the first
  `nfix` physical cells, replacing their rows by `Y - Yfix`.

**Consequence, and the code says so**: the marching production map contains
`chem`, `transport` and `filter` terms the residual does not
(`EXHALE_main.f90`). The marching fixed point and the residual zero
are therefore not the same object except where those terms vanish.

### 7.3 Residual evaluation order

`eval_residual` (`steady_newton.f90`) reproduces the marching order
:

1. `f_sp = f_sp_seed`: the seed is an explicit argument, so `F` is a
   function of `Y` alone. The header records the measurement that forced this
  : with a shared workspace, 1406 of 1500 entries of `F` at one state
   differed by 3.0e-3 of the row scale against a 1e-5 tolerance; with the seed
   held, three evaluations are bitwise identical.
2. `unpack_U(Y,u)`, with the carrier unknown written into
   `f_sp(:,isp_H2)` and checked against `carrier_headroom`.
3. `get_species_densities` on the interior **before** the ghosts are filled,
   because `Apply_BC` reads the global `n_part_cell1`.
4. `Apply_BC` -> `U_to_W` -> `get_species_densities` -> `comp_T_from_p` ->
   `excited_H_update` -> `ioniz_eq` -> `get_species_densities` -> **a second**
   `Apply_BC`, so the ghosts carry the composition the fluxes will
   use.
5. `assemble_residual`, `pack_R`.
6. the carrier row when present, converted by `tscale_code = R0/(v0 n0)`.
7. admissibility gating.

`frozen_residual` is the same with `heat`, `cool` and the particle
count held fixed, hence strictly local and exactly banded.

### 7.4 The solver

**Preconditioner**: banded, built by colored finite differences from the
**frozen-radiation** residual (`build_banded_jac`); geometry
`kl = ku = 8, ncolor = 17` for `nvar = 3` and `kl = ku = 11, ncolor = 23` for
`nvar = 4`, the bandwidth following from the WENO3 stencil.
Factored with LAPACK `dgbtrf`/`dgbtrs` and used as the **right** preconditioner;
the weakly non-local radiation response is left to the outer iteration.

**Krylov**: `pgmres`, right-preconditioned GMRES(m), a **single
cycle**, `m = 40` at the marching hand-off, relative tolerance 1e-1.
Matrix-free products at `jv_product`. A probe
point outside the describable set is answered by halving the finite-difference
step, at most six times; if even the smallest is inadmissible the
cycle stops at the directions it built.

**Scaling**: `cell_state_scales` is the column scaling;
`cell_row_scales` copies it for the three hydrodynamic rows and
differs only for the carrier row. Since section 143 of the changelog the row
scale of the measure is separate from the column scaling, because scaling the
system by the measure's scales left the hot-Uranus solve with no descent
direction after 179 iterations against 9.

**PTC shift**: the system solved is `(I/dtau + J) dY = -F`, with an
SER-style ramp on acceptance and `dtau = max(0.25 dtau, dtau0)` on rejection.

**Trust region**: `EXHALE_TRUST_REGION=1`, and **only on the four-unknown
branch**; dogleg with `tr_delta`, Newton leg at
`tr_gm_rtol` default 1e-1. Off by default, and when on it replaces the line
search and the damped Gauss-Newton escape.

**Line search**: backtracking on the scaled merit, at most 20 trials, with a
non-monotone Grippo test over 5 iterates and a switch to monotone Armijo after a
stagnation restart. Acceptance requires positivity of the density
and pressure on the physical cells and non-negativity of the carrier unknown, and
a trial the code cannot describe is rejected **before** the merit comparison.
Trials are evaluated with WENO weights **recomputed**, while the
inner Newton model uses frozen weights, because the frozen test
accepted 7 of 14 steps that raised the true residual.

**Escapes**: a damped Gauss-Newton when `dtau` is at its floor;
`n_no_descent_max = 12`; `n_stall_best = 20`, after which the solve returns to
the best iterate and switches to monotone Armijo, and after another 20 stops as
STAGNATED, distinct from "no descent direction exists".

**Acceptance**: `steady_gates_met` (`steady_residual.f90`), the single
definition used by JFNK, PTC and the marching stop, is three gates: the residual
gate `rnorm < resid_tol`; the **flux gate** `flux_spread_of_state(u) <
flux_spread_th`, the spread of the Riemann **face** mass flux over faces bounding
cells at or above `r_flux = 1.2` with the sign kept; and the carrier
gate, applied only when the caller says the carrier is an unknown,
which the marching loop does not (`EXHALE_main.f90`). `rnorm` is
`resid_relnorm`, which norms the wind and the layer **separately
and combines them by the larger**, because a sum would let the outer volume
average the inner column away; the carrier row is deliberately not
folded in. The returned state is the best iterate by residual, but
never one that fails the flux gate over one that passes it.

Defaults: `resid_th = -1.0` with `N_resid = 500` (`parameters.f90`),
and `1e-5` when `Resid tol:` was not given (`EXHALE_main.f90`);
`flux_spread_th = 2.0e-5`, `r_flux = 1.2` (`parameters.f90`);
`carrier_resid_th = 1.0e-3`; `newton_du_switch = 1.0e-2`.

**Outer loop**: `steady_wind_with_element_diffusion`
(`EXHALE_main.f90`), a damped Picard alternation of the steady solve at
fixed composition, `relax_element_composition` for He/H, and
`relax_photochemical_composition` for the carriers when they are not Newton
unknowns; `omega` from 0.5 down to a floor of 0.125; at most 20 passes. After
each solve the state is refreshed from `u` through the full chain including a
second `comp_T_from_p` so that `(p, T, f_sp)` is consistent.

**Composition consistency at the end**: `steady_newton.f90` repeatedly
feeds the output composition into a new residual evaluation at fixed
hydrodynamic unknowns and stops when the scalar residual norm changes little.
The audit (N7) notes that a stable norm can coexist with a moving weak species,
and that during a single evaluation the rates use the temperature derived from
the seed composition while the post-sweep composition can change the equation of
state and the particle count before the hydrodynamic residual is assembled.

### 7.5 The consequence for a transported proton

`Ionization transport: True` refuses `Solver: Newton`, `EXHALE_PTC` and
`Coupled carrier solve` (`input_read.f90`, `input_read`), and the steady unknown packing
supports three hydrodynamic quantities plus optionally one H2 density and no
proton continuity unknown. So a molecular layer with a transported H+ has
**no steady solver**: marching cannot reach its 1e7 to 1e8 s relaxation time and
Newton cannot hold it. This is C19, and it is the first entry of section 4 of
`docs/code_status_20260905.md`.

### 7.6 The post-process the transit tools consume

`post_process_adv.f90` replaces the local ionization balance of each cell by the
steady advection-ionization ordinary differential equation integrated upwind
across the cell, then re-solves the stationary energy equation on the
corrected composition, and writes `Hydro_ioniz_adv.txt` and
`Ion_species_adv.txt`. Ten outer passes. Density and velocity are **not**
recomputed.

Stated assumptions, which matter because these are the profiles the transit
tools read:

1. **Molecule-free reconstruction** (module header): H2, H2+, H3+, HeH+, OH,
   H2O and CO are excluded; `calc_ne` and `calc_ntot` are called without `nmol`.
   The molecular and oxygen columns of the `_adv` files are the **equilibrium**
   values, and `write_output.f90` prints that limitation into the file.
   This is C20.
2. Three conditions under which a cell keeps the equilibrium ionization:
   inflow (`v <= 0`), because the upwind stencil takes the wrong
   upstream cell; a Damkohler number above `Da_local_equilibrium = 1e2` formed
   from the **slowest** relaxation rate among all solved populations; and
   `x_HII,eq < 1e-6`, because the residuals carry `x_HI` and
   the extracted ion fraction is then worse than 1 per cent relative.
3. In the energy equation only the inflow condition is applied,
   because the other two are statements about the ionization balance and not
   about the discretization.
4. Metals follow `pp_metals`: 0 metal-free, 1 frozen equilibrium (the default),
   2 re-solved.

---

## 8. The term table and the consistency items

### 8.1 Consistency items

| id | statement | lines |
|---|---|---|
| **C1** | The composition update changes thermal energy at fixed temperature, separately from the explicit energy sources. Not obtained by integrating any heating term. | `EXHALE_main.f90`; audit N1  **Closed 2026-09-06 (B3c):** rows 7 to 9 are one local source step at fixed volume, the energy row anchored on `u_th_old`; no `comp_p_from_T` and no `W_to_U` in the loop. |
| **C2** | The chemistry stage overwrites `rho` and pushes it into the conserved mass row, so it is not conservative in `u`; the steady residual contains no such term. | `ionization_equilibrium.f90`; `EXHALE_main.f90`  **Closed 2026-09-06 (B3c):** the density is `intent(in)` to `ioniz_eq` and the mass sum of the returned composition is a checked diagnostic. |
| **C3** | The semi-implicit energy update performs two iterations but re-evaluates cooling once, so the returned `cool` belongs to the previous iterate; `heat` is frozen and contains temperature-dependent reaction terms; no residual acceptance test. | `energy_semi_implicit.f90`; audit N5  **Closed 2026-09-06 (B2):** the update is a bracketed, residual-controlled solve; cooling at the returned temperature; a floor is a failure, not a clamp. |
| **C4** | `Hydro_ioniz.txt` `heat` and `Heating_breakdown.txt` `heat_total` are not one state (quoted RED at 1.2e-5). Separately, the test's source citation is stale: the block it names no longer exists. | test; `EXHALE_main.f90`; `docs/Update_EXHALE_stage2.md` section 6  **Citation half closed 2026-09-06** (test header rewritten to the rate lag; measured RED at 1.19e-5, header `sec_ion=F`); the rate-lag half stays for Phase 3. |
| **C5** | Two photon-grid constructions: threshold-cut geometric bins for the analytic spectra, central differences around table nodes with no threshold cutting for a loaded SED. | `set_energy_vectors.f90` vs `sed_read.f90`  **Closed 2026-09-06** (item 2c-SEDGRID, `Update_EXHALE_stage2.md` section 6: the loaded grid is edge-based with every active threshold inserted; `stellar_flux_eV` reads the table rows `e_sed_node`/`F_sed_node`). |
| **C6** | The Balmer continuum is integrated on its own 400-point trapezoid over a band the photon grid also covers, and H(n=2) is not an absorber in `tauE`, so it removes no photons from the beam while its band neighbors do. | `excited_hydrogen.f90`; `util_ion_eq.f90` |
| **C7** | Lyman-alpha: the cooling function charges the full 10.2 eV per collisional excitation as escaping radiation with no escape probability, in the same cells where `lya_rt` computes an escape probability far below one, while the de-excitation heat returns 10.2 eV. | `Cool_coeff.f90` and `util_ion_eq.f90` vs `lya_rt.f90` and `excited_hydrogen.f90` |
| **C8** | The same H/He and metal cooling sum is written twice, and the post-process temperature is solved with a different set of channels from the hydro temperature: no H3+, no molecular infrared bands, no molecular electrons in `n_e`, no Balmer, Lyman-Werner, molecular chemical or FUV heat; `post_process_adv.f90` calls `eval_cool` without `nheiTR`. | `util_ion_eq.f90` vs `T_equation.f90`; `post_process_adv.f90` |
| **C9** | The `Lyman_Werner.txt` `heat_LW` column carries only the fragment half of the Lyman-Werner heat; the energy equation applies fragment plus fluorescence, the latter the larger. | `write_output.f90` vs `ionization_equilibrium.f90` |
| **C10** | Two case-B hydrogen recombination fits: one drives the ionization balance, the other the n=2 cascade and the Lyman-alpha source. `Atomic rate set:` moves only the first. | `Cool_coeff.f90` vs `hydrogen_n2_rates.f90`, used `excited_hydrogen.f90`, `lya_rt.f90` |
| **C11** | `lambda_rec_HII` is an independent Hui and Gnedin fit and is **not** switched by `Atomic rate set:`, while its balance partner is. | `Cool_coeff.f90` |
| **C12** | Byte-equal duplicates with no compiler-visible link: Voronov H I and He I coefficients; the Mao and Kaastra He II alpha_1 twelve lines apart; the bremsstrahlung prefactor; and the entire H(n=2) rate set restated in the Python transit tool. | `Cool_coeff.f90` vs `mol_rates.f90`; `Cool_coeff.f90`; `util_ion_eq.f90` vs `T_equation.f90`; `exhale_transit_lib.py:136-175` vs `hydrogen_n2_rates.f90` |
| **C13** | A dead duplicate of the He/H charge exchange pair, about 4 per cent and 500 K away from the live one, with no callers. | `mol_rates.f90` vs `charge_exchange.f90` |
| **C14** | The opacity dispatcher covers H I, He I, He II and He 2^3S only; H2 and the 17 metal ions bypass it while `opa_pf` still multiplies their columns and their local absorption. | `set_energy_vectors.f90` vs |
| **C15** | The pressure enters the momentum row two different ways: inside the flux plus a geometric source for PLM, as a direct gradient for WENO3. | `Num_Fluxes.f90` and `Source.f90` vs `RK_rhs.f90` |
| **C16** | The dummy `alpha` is threaded through the whole flux path, never assigned and never read; the code reads an uninitialized variable on every flux call. | `EXHALE_main.f90`; `steady_residual.f90`; `RK_rhs.f90`; `Num_Fluxes.f90`  **Closed 2026-09-06**: the dummy removed from the three signatures and ten call sites; `hydrostatic_column` byte-identical against a same-compiler control. |
| **C17** | The steady residual carries three (or four) rows and no species rows; the composition is eliminated by a nested equilibrium and element diffusion is outside the residual entirely, co-converged by Picard. | `steady_residual.f90`; `steady_newton.f90`; `EXHALE_main.f90` |
| **C18** | The imposed proton partition is applied only for `j >= 1`; the two lower ghosts re-solve their ionization locally against the same rates. | `ionization_equilibrium.f90` |
| **C19** | `Ionization transport: True` refuses the Newton solver, PTC and the coupled carrier solve, and the steady packing has no proton unknown, so that configuration has no steady solver at all. | `input_read.f90`, `input_read`; `steady_newton.f90` |
| **C20** | The `_adv` profiles the transit tools read are reconstructed on a molecule-free H/He background; their molecular columns are the equilibrium sweep's values. | `post_process_adv.f90`; `write_output.f90` |
| **C21** | Species advection does not use the hydro face Riemann mass flux or the spherical cell volume: the carrier residual uses a cell-centered derivative of a mass-normalized fraction and `r_center^2 dr` for the diffusion divergence, and the base direction comes from a wind-region average while the magnitude is local. | audit section 5.4; `diffusive_photochemistry.f90`  **Transport half closed 2026-09-07 (B4-1c):** in the marching loop the carriers are advected as `F_rho Y_c/m_c` on the hydro face mass flux inside the Runge-Kutta stages; the cell-centered form survives only in the fixed-wind relaxation. |
| **C22** | Chemistry and the caloric equation of state use different H2 internal-state models, so one thermodynamic potential does not generate both. | `caloric_eos.f90` (Roueff ladder) vs `mol_rates::q_rovib_H2`; audit P4  **Closed 2026-09-06 (B3b-H2Q, B3b-IR):** one `h2_partition_function` over the Roueff ladder serves the EOS, the chemistry's `K_eq` and the infrared populations. |
| **C23** | H2 double photoionization spends the 51.4 eV vertical threshold in the heating calculation while the represented products hold 31.675 eV; the difference has no identified recipient. | `util_ion_eq.f90`; `h2_photo_channels.f90`; audit P1  **Closed 2026-09-06 (B3b-PE, wired by B4-6):** the 19.725 eV fragment kinetic energy of the double channel and the H(2p)+H(2s) prompt radiation of the neutral window are assigned by `h2_channel_energy_recipients`. |
| **C24** | H3+ cooling clamps the collider density to the first table column, so it does not vanish as `n(H2) -> 0`, and the published piecewise fits join discontinuously. | `h3p_cooling.f90`; audit P2, P3  **Closed 2026-09-06 (B3b-H3):** below the table the non-LTE factor follows the exact collisional limit, linear in `n(H2)`, zero at zero. |
| **C25** | The oxygen cycle has photolysis heat but no collisional reaction-energy ledger; the retained O(1D) excitation energy is stated as absent from the heat sum; the He(2^3S) associative branch heat is not deposited anywhere. | `ionization_equilibrium.f90`; `molecular_reaction_heat.f90`; audit P5 |
| **C26** | CO has no balance row and no kinetic destruction rate: it is an inert transported reservoir until an algebraic equilibrium ceiling removes it. | `ionization_equilibrium.f90`; `diffusive_photochemistry.f90` `limit_to_element_budget`; audit P6  **Closed 2026-09-06 (B3b-CO, B3b-CO2, CEILING-DEL):** the CO carrier row carries He+ + CO -> C+ + O + He at the UMIST RATE22 4068 rate and shielded photodissociation on the Lyman-Werner beam (Visser et al. 2009), with the He+ sink in the molecular He+ row and the domain of the destruction-only row measured in each cell and reported; the equilibrium ceiling and its four counters are deleted. MEASURED at the deletion: the ceiling never fired in the current tree, so every output column was byte-identical across it. |
| **C27** | The absorber columns are a rectangle rule on cell widths carrying `opa_pf`; the Lyman-alpha line-center optical depth is a trapezoid on cell centers without it. | `utilities.f90` vs `lya_rt.f90` |
| **C28** | Element diffusion is solved in advective, not conservative, form, so helium mass is conserved only to the order of the hydro truncation in a transient. | `binary_element_diffusion.f90`  **Closed 2026-09-06 (B4-1):** the marching route carries the element advection on the hydro face mass fluxes inside the Runge-Kutta stages and `solve_mass_fraction` runs with `advect` false; the advective form is left to the fixed-wind relaxation. |
| **C29** | Zero net diffusive mass flux is structural for the H/He binary but is **not** enforced for the trace metals, which move against a fixed hydrogen background with no counter-flux. | `binary_element_diffusion.f90` |
| **C30** | Molecular ions contribute to `n_e` but are excluded from the bremsstrahlung charge sum. | `util_ion_eq.f90` |
| **C31** | The three-dimensional geometry is one scalar dilution on the flux; columns are radial and plane-parallel; the shielded-band shell average is not implemented, and the fluid closure is never replaced when the collisional diagnostic fails. | `parameters.f90`; `utilities.f90`; audit sections 5.5, 5.7 |
| **C32** | The Riemann solvers use a frozen-composition sound speed; the equilibrium reacting compressibility differs and the stiff relaxation limit has not been shown to recover it. | `caloric_eos.f90`; audit section 6.5 |
| **C33** | Supersonic inflow at the base is over-specified by the characteristic count and is only diagnosed, never prevented. | `base_boundary.f90`; `Apply_BC.f90`; `EXHALE_main.f90` |

### 8.2 Term table

| term | equation | routine:line | source | validity | consistency |
|---|---|---|---|---|---|
| mass flux | `(1/r^2) d_r(r^2 rho v)` | `RK_rhs.f90`; `Num_Fluxes.f90` | Godunov, ATES heritage | collisional single fluid | C2 |
| momentum flux and pressure | `d_r(rho v^2) + d_r p` | `RK_rhs.f90`; `Source.f90` | as above | PLM and WENO3 differ in form | C15 |
| energy flux | `(1/r^2) d_r(r^2 v(E+p))` | `RK_rhs.f90` | as above | | OK |
| gravity, momentum | `-rho d_r Phi` | `Source.f90`; `grav_field.f90` | Roche or spherical, model choice | 1-D Roche is not a 3-D overflow solution | C31 |
| gravity, energy | flux form `dF3p` | `RK_rhs.f90` | well-balanced form | | OK |
| numerical dissipation, LLF | `alpha = max(|v_L|+a_L, |v_R|+a_R)` | `Num_Fluxes.f90` | Rusanov 1961; Zhang and Shu 2010 | | OK, batch 2b |
| ROE star state | Toro 2009 eqs. 9.20, 9.32, 9.42, 4.76 | `speed_estimate_ROE.f90` | Toro (2009) ch. 9 | constant index only; refused with molecules | OK, batch 2a |
| low-Mach damping | JST fourth difference | `low_mach_dissipation.f90` | Jameson, Schmidt and Turkel 1981 | numerical only; `eps4 < 1/(16 CFL)` | OK, default off |
| caloric EOS | `e, C_V, gamma_eff` | `caloric_eos.f90` | Roueff et al. 2019 A&A 630 A58 table 2 | H2 ladder only; trace molecules monatomic; ortho/para equilibrium | C22 |
| viscosity | `(4/3) mu (w' - w/r)` and `q_mu` | `viscous_conduction.f90` | Navier-Stokes; the code rejects Koskinen 2022 (B5)/(B6) as printed | neutral atomic H coefficient; no Coulomb viscosity | OK, default off |
| conduction | `(1/r^2) d_r(r^2 kappa d_r T)` | `viscous_conduction.f90` | Watson, Donahue and Walker 1981 | no Spitzer electron conduction | OK, default off |
| base boundary | LODI `C^-` at the face plus a (p,s) reservoir | `base_boundary.f90` | Thompson 1987; Poinsot and Lele 1992; Carlson 2011 section 2 | supersonic inflow not well posed | C33 |
| outer boundary | zero gradient or linear extrapolation | `Apply_BC.f90` | model choice | | OK |
| incident field | `J_inc`, Planck, or a loaded table | `J_inc.f90`; `sed_read.f90` | ATES parameterization; Rybicki and Lightman 1979 | Planck is an upper bound below 13.6 eV and a lower bound above | OK |
| photon grid | geometric bins cut at every active threshold | `set_energy_vectors.f90` | quadrature design, batch 2c | analytic spectra only | **C5** |
| column density | rectangle rule to the inner face | `utilities.f90` | model choice, radial plane-parallel beam | no slant path, no shell average | C27, C31 |
| photoionization rate | `INT F_E e^{-tau} sigma/E dE` | `util_ion_eq.f90` | Verner et al. 1996; hydrogenic | | OK |
| photoheating | `(1 - E_th/E)` times the above | `util_ion_eq.f90` | standard `h nu - I` | comparison option charges the whole photon | C23 |
| H2 channels | four disjoint final states | `h2_photo_channels.f90` | Chung et al. 1993 Tables I and II; Yan et al. 1998 section 4 | `q_D` is a model; data uncertainty 4 to 5 per cent | C23 |
| H(n=2) rate and heat | `gamma_2`, `hpe_2`, 2s/2p system | `excited_hydrogen.f90` | Christie, Arras and Li 2013 eqs. 12-13; Osterbrock and Ferland 2006 | lagged explicit; own quadrature; not an absorber | C6, C10 |
| secondary ionization | SvS85 fractions, DYL99 heating | `electron_energy_degradation.f90` | Shull and van Steenberg 1985 eqs. 1-2 Table 2; Dalgarno, Yan and Liu 1999 eqs. 12-14 | eq. 14 stated for `x <= 0.1`, extrapolated to 1; He/H fixed at 0.1; no metals; 30 eV floor | OK, limits stated |
| Lyman-alpha field | static-slab escape probability | `lya_rt.f90` | Neufeld 1990 eqs. 2.25, 3.27; Harrington 1973 eq. 40; Huang et al. 2023 eq. 11 | static slab damping wings, `(a tau)^(1/3) > 10`; boost is a tuned parameter | C7 |
| Lyman-Werner rate | `sigma_LW F/<h nu>` times shielding | `lyman_werner.f90`, `lyman_werner_dissociation_rate` and `lyman_werner_dissociation_rate_cell_mean` | Draine and Bertoldi 1996 Table 1, Table 2, Fig. 1 annotation | 25 per cent on band shape; cell mean to 2.5e-4 | OK, batch 2b |
| H2 self-shielding | table in `(N_H2, T, n_H)` | `h2_self_shielding_table.f90` | Abgrall, Roueff and Drira 2000; CLOUDY trapping ratio | 700-3200 K, 1e12-1e14 cm^-3, 1e12-5e21 cm^-2, clamped outside | OK, limits stated |
| FUV photolysis | `s_b N_b exp(-tau_b)` | `water_photolysis.f90` | Huebner and Mukherjee 2015; Heays et al. 2017; Ranjan et al. 2020 | band-mean is a lower bound; **B4 off by 4.6 to 6.1** | OK, stated |
| H/He/metal equilibrium | `N_eq` rows, fractions | `ion_residual_core.f90`; `System_HeH_mol.f90` | Badnell; Voronov 1997; Huang et al. 2023 Table 4; Koskinen et al. 2022 Table 1 | Voronov metal transcription not verified at the code site | C10-C13 |
| Penning | Taylor 2025 fit, `f_penning = 0.9` | `Cool_coeff.f90`; `ion_residual_core.f90` | Taylor 2025 | 200-10000 K, evaluated outside unclamped | OK, stated |
| carrier transport | `d_t n + div(n v + Phi) = P - L` | `diffusive_photochemistry.f90` | Blanc's law mixture; hard-sphere pairs | trace-solute reading of H2 when H2 is most of the gas | C21 |
| element diffusion | `rho d_t X + div(r^2 J)/r^2 + rho v d_r X = 0` | `binary_element_diffusion.f90` | Chapman-Enskog; Koskinen et al. 2013 section 2.1 for the field; Paquette et al. 1986; Spitzer 1962 | advective form; two sources not held locally and re-derived | C28, C29 |
| element census | ratio and absolute nucleus invariants | `element_census.f90` | exact bookkeeping identity | checks only; off unless the variable is set | OK |
| heating sum | 9 additive channels | `ionization_equilibrium.f90` | see section 5.2 | | C9, C25 |
| cooling sum | `ne(brem+coex+reco+coio) + metals + H3+ + IR` | `util_ion_eq.f90` | Hui and Gnedin 1997; Cen-type; van Hoof et al. 2014; CHIANTI v11.0.2; Miller et al. 2013; Roueff et al. 2019 | CHIANTI fits 1e3-1e5 K with a 1000 K cutoff; IR tables 50-2000 K clamped | C7, C8, C24, C30 |
| marching energy update | `e(T_new) - e(T_old) = dt(heat - cool(T_new))` | `energy_semi_implicit.f90` | model choice | two iterations, no residual test | C1, C3 |
| post-process temperature | steady advected balance in one cell | `T_equation.f90` | model choice | different channel set from the hydro | C8 |
| steady residual | `L(U) + S = 0` on three or four rows | `steady_residual.f90` | model choice | no species rows, no diffusion | C17, C19 |
| steady acceptance | residual, flux spread, carrier | `steady_residual.f90` | model choice | the carrier gate is off unless the carrier is an unknown | OK |

---

## 9. Open questions for the user

Each is a place where the code admits two readings, or where the plan records a
decision that is still open. One sentence each.

1. **Energy convention (plan decision 2).** Does `U(3)` stay kinetic plus
   thermal with the external-exchange identity and the local ownership rule, or
   does it become a formation-inclusive total energy?
2. **The composition projection (C1).** Is the projection to be removed by
   deriving a conservative source update, or is the thermal-only convention to
   be kept with matched chemical source terms?
3. **The mass row (C2).** Should the chemistry stage stop writing `rho`, or
   should the steady residual gain the corresponding production term so that the
   two maps have the same fixed point?
4. **Charged transport closure (plan decision 5).** Derived zero net current
   with the ambipolar field, or a common-velocity approximation justified over an
   explicit domain, given that today the proton has zero molecular diffusion and
   is neither?
5. **CO (plan decision 6). Answered 2026-09-06:** the one-sided destruction
   model was derived and validated with its domain argument (items B3b-CO,
   B3b-CO2), and the ceiling was then deleted rather than its runs excluded
   (item CEILING-DEL). Original question: derive and validate a one-sided
   destruction model with a timescale and domain argument, or exclude
   ceiling-active runs?
6. **H3+ emission model (plan decision 7).** Scaled LTE with the Table 6
   departures, or unscaled LTE with the factors as a separate model, and over
   what temperature range?
7. **Oxygen ledger (plan decision 9).** If the oxygen reaction-energy ledger is
   still incomplete at the end of Phase 3, label and exclude those
   configurations, or hold them?
8. **Atomic-line ionization transport (plan decision 10).** Should the H+ carrier
   be selectable without molecular chemistry, and if so is it reported as an
   approximation?
9. **The carrier retry design (plan decision 4).** A carrier substep retry inside
   the transport step, or an outer attempted-step controller re-entering at the
   hydrodynamic stage?
10. **The loaded-SED quadrature (C5).** Should the threshold-cut grid
    construction be extended to a loaded table, given that the planet examples
    now use literature SEDs?
11. **The `_adv` profiles (C8, C20).** Should the post-process temperature be
    solved with the same cooling function as the hydro, and should the `_adv`
    files stop carrying molecular columns that were never corrected?
12. **The `T_equation` and `eval_cool` duplication (C8).** One assembly with two
    entry points, or two assemblies with a test that pins them together?
13. **Trace-metal diffusion (C29).** Should the trace-metal solve gain a counter-flux in
    the hydrogen component, or should it be labeled a trace approximation with a
    stated validity domain?
14. **`alpha` (C16).** Removed 2026-09-06 (no statement read it); no question remains.
    (Its removal changes four interfaces and no numbers.)
15. **The stale test citation (C4).** Fixed 2026-09-06 with the rate-lag description; the mismatch itself (1.19e-5) is the Phase 3 item.
    `output_state_consistency.sh` and its README be corrected now, or after the
    Phase 3 item that owns the remaining 1.2e-5 mismatch?

---

## 10. What this document does not establish

- It does not measure anything. Every number is read from the source, from a
  source comment, or from a cited document, and none was re-measured.
- It does not verify any citation against its publication. Every attribution in
  sections 2 to 5 is what the code comment claims; several carry no table or
  equation number, and the metal collisional-ionization block carries the code's
  own unresolved "VERIFY" note (`Cool_coeff.f90`).
- It does not cover the Wind-AE reimplementation (`src/modules/wind_ae/`), which
  is a separate solver with its own recombination and cooling data, nor the
  transit post-processing beyond the point where it consumes the `_adv` profiles.
- It does not derive the charged transport closure, the CO model, or the
  independent species space for the stationary system. Those are items 5, 9 and 6
  of the plan's section 4.2 and are the subject of the decisions in section 9.
