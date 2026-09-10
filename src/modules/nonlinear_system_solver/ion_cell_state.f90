	module ion_cell_state
	! Named-field cell state for the equilibrium ionization residuals. It
	! replaces reading the position-packed params() inside the residuals
	! (§5.2 Inc 4, docs/refactor_plan_system_composition_parser.md): params
	! stays as the MINPACK transport argument (hybrd1/fdjac1/solve_ieq), but
	! the converted equilibrium systems no longer read it -- they read the
	! named fields below instead. The field names follow the physical
	! quantities of the equilibrium params layout.
	!
	! The ionization cell sweep (ioniz_eq) is OpenMP-parallel over cells, so
	! the single module instance is threadprivate: each thread fills its own
	! copy field by field, exactly like the threadprivate met_* arrays in
	! System_HeH_metals. All components are fixed-size scalars (no allocatable
	! components), so the threadprivate copy needs no allocation of its own on each thread.
	!
	! The advection-correction residuals (System_implicit_adv_H/HeH/HeH_TR) and
	! the post-process temperature residual (T_equation / solve_T_brent) read
	! their own named-field types (adv_rates/teq_state) below. Their loop in
	! post_process_adv is serial (master thread), but the instances are kept
	! threadprivate for uniformity with ieq_cell.

	implicit none
	public

	! EVERY FIELD CARRIES ITS ZERO AS ITS DEFAULT, so that an ion_rates that
	! has not been filled yet is the state of a cell with no gas and no
	! rates: no photoionization, no recombination, no collisional or
	! charge-exchange channel, no third body, and no imposed carrier
	! partition (the three x_*_fixed flags false). That is the only value
	! this type has an unambiguous meaning for. The invariant matters most
	! for the grid-sized arrays of ionization_equilibrium (bg_cell,
	! ieq_rate_cell and the two solver copies), which are allocated before
	! the first sweep and hashed, saved and restored by the attempted-step
	! checkpoint from that moment on: without the defaults, an allocation
	! served from a recycled block leaves them holding heap contents, and a
	! true in x_h2_fixed read that way would impose an H2 partition no
	! caller set.
	type ion_rates
		real*8 :: P_HI = 0.0d0
		real*8 :: P_HeI = 0.0d0
		real*8 :: P_HeII = 0.0d0
		real*8 :: rchiiB = 0.0d0
		real*8 :: rcheiiB = 0.0d0
		real*8 :: rcheiiiB = 0.0d0
		real*8 :: nh = 0.0d0
		real*8 :: nhe = 0.0d0
		real*8 :: a_ion_HI = 0.0d0
		real*8 :: a_ion_HeI = 0.0d0
		real*8 :: a_ion_HeII = 0.0d0
		real*8 :: a_ion_HeITR = 0.0d0
		real*8 :: rcheiTR = 0.0d0
		real*8 :: A31 = 0.0d0
		real*8 :: P_HeITR = 0.0d0
		real*8 :: q13 = 0.0d0
		real*8 :: q31a = 0.0d0
		real*8 :: q31b = 0.0d0
		real*8 :: Q31 = 0.0d0
		real*8 :: P_H2 = 0.0d0
		! The dissociative part of P_H2 [s^-1]: the sub-rate of it that
		! leaves H + H+ + e- instead of H2+ + e- (threshold 18.08 eV against
		! 15.4). A SUBSET of P_H2, so P_H2 is still the whole H2 destruction
		! rate and P_H2 - P_H2_di is what makes H2+.
		real*8 :: P_H2_di = 0.0d0
		! The double-ionization part of P_H2 [s^-1]: the sub-rate that
		! leaves H+ + H+ + 2e- (threshold 51.4 eV). A SUBSET of P_H2, like
		! P_H2_di and disjoint from it, so P_H2 is still the whole H2
		! destruction rate. One event makes TWO protons.
		! Zero unless a double-ionization model is selected.
		real*8 :: P_H2_dd = 0.0d0
		! The neutral-dissociation part of P_H2 [s^-1]: the sub-rate that
		! leaves H + H and no ion at all (the sub-unity photoionization
		! yield over 33-41 eV). A SUBSET of P_H2, like the two above and
		! disjoint from both. Zero unless h2_neutral_dissociation.
		real*8 :: P_H2_nd = 0.0d0
		! Lyman-Werner photodissociation rate of H2 [s^-1], already carrying
		! the self-shielding of the star-ward H2 column (lyman_werner.f90).
		! Zero unless the run supplies a Lyman-Werner band flux. Assigned by
		! ioniz_eq for every molecular cell.
		real*8 :: k_LW = 0.0d0
		real*8 :: T_K = 0.0d0
		! Total gas-particle density [cm^-3], electrons excluded: the third
		! body M of the three-body molecular reactions R12/R13/R15. It is
		! calc_ntot's sum (one particle per species), NOT rho/m_H.
		real*8 :: ntot = 0.0d0
		! He <-> H charge-exchange rate coefficients (Huang Table 4 group B):
		! kcx_He0_Hp = He0+H+ -> He++H0, kcx_Hep_H0 = He++H0 -> He0+H+.
		real*8 :: kcx_He0_Hp = 0.0d0
		real*8 :: kcx_Hep_H0 = 0.0d0
		! Oxygen chemistry (the A2 option). n_ofam is the FREE OXYGEN FAMILY
		! of the cell [cm^-3]: every oxygen nucleus except the one locked in
		! CO, i.e. the reservoir shared by O I, O II, O III, OH and H2O. The
		! oxygen unknowns are fractions of it, and its counterpart n_co is
		! the carbon monoxide density, a background of the solve: CO has no
		! formation channel in the one-sided model (B3b-CO), so it is not a
		! balance row; it is the carrier's value in the cell, destroyed by
		! He+ + CO and by photodissociation in the carrier row, and read by
		! the He+ row of mol_heh_rows for the same He+ + CO loss (B3b-CO2).
		! Assigned by ioniz_eq for every cell of an oxygen-chemistry run.
		real*8 :: n_ofam = 0.0d0
		real*8 :: n_co = 0.0d0
		! Imposed carrier partitions. Where a partition is not a local root,
		! the balance row that would have computed it is replaced by the
		! value, and the other rows keep their balances and are solved
		! against it -- which is what keeps the ionization stages consistent
		! with the imposed partition. The fractions are in the same units as
		! the unknowns they replace:
		!   x_h2_fix  = 2 n_H2 / n_H(nuclei)      -> x(4)
		!   x_oh_fix  = n_OH  / n_O(free family)  -> x(iox)
		!   x_h2o_fix = n_H2O / n_O(free family)  -> x(iox+1)
		!
		! Hydrogen and oxygen carry SEPARATE flags because they are imposed
		! for different physical reasons and do not always occur together:
		!
		!   x_h2_fixed   the H2 partition is owned by something other than
		!                this cell's local balance -- the transported
		!                carriers (milestone M3, diffusive_photochemistry),
		!                or the lower-boundary reservoir, whose composition
		!                is incoming data the shielded base cell cannot
		!                derive for itself (base_h2_composition_imposed).
		!   x_ox_fixed   the OH/H2O partition is owned by the transport
		!                solve. Only carrier transport sets this one.
		!
		! One flag for both would tie them together: imposing H2 at the base
		! ghost would then silently impose OH and H2O as well, at whatever
		! the oxygen fractions happened to hold, and zero the oxygen
		! carriers of a run that never asked for it.
		!   x_hp_fixed   the H/H+ partition is owned by the transported
		!                proton (Ionization transport: True), not by this
		!                cell's local photoionization/recombination
		!                balance. A THIRD flag, for the same reason the
		!                first two are separate: it is imposed for a
		!                different physical reason (the ionization time
		!                exceeds the flow time) and in a different set of
		!                runs, and tying it to x_h2_fixed would impose an
		!                ionization fraction on every run that carries H2.
		!                x_hp_fix = n_H+ / n_H(nuclei)        -> x(1)
		logical :: x_h2_fixed = .false.
		logical :: x_ox_fixed = .false.
		logical :: x_hp_fixed = .false.
		real*8 :: x_h2_fix = 0.0d0
		real*8 :: x_oh_fix = 0.0d0
		real*8 :: x_h2o_fix = 0.0d0
		real*8 :: x_hp_fix = 0.0d0
	end type ion_rates

	type(ion_rates), save :: ieq_cell
	!$omp threadprivate(ieq_cell)

	! Named-field cell state for the advection-correction ionization residuals
	! (System_implicit_adv_H/HeH/HeH_TR). One type covers BOTH advection
	! layouts: the H-only and HeH systems read the leading fields, the HeH_TR
	! system additionally reads the triplet channels (rcheiTR..xheiTR_old).
	! heh_loc is the effective He/H for the electron density (global HeH
	! normally, the local nhe/nh with He_diffusion). Filled field by field by
	! post_process_adv before each hybrd1 advection solve.
	!
	! xheiS_old is the upstream population of the GROUND SINGLET He(1^1S)
	! alone, not of the summed He I: the advection systems carry the singlet
	! and the He(2^3S) metastable as separate unknowns, so that neither is
	! ever formed as the difference of the other two (a difference that
	! collapses to zero once the metastable holds most of the neutral He).
	! Without the triplet the two coincide, all He I being in the singlet.
	type adv_rates
		real*8 :: c1
		real*8 :: xhi_old
		real*8 :: xheiS_old
		real*8 :: xheiii_old
		real*8 :: nh
		real*8 :: P_HI
		real*8 :: P_HeI
		real*8 :: P_HeII
		real*8 :: rchiiB
		real*8 :: rcheiiB
		real*8 :: rcheiiiB
		real*8 :: a_ion_HI
		real*8 :: a_ion_HeI
		real*8 :: a_ion_HeII
		real*8 :: a_ion_HeITR
		real*8 :: heh_loc
		real*8 :: rcheiTR
		real*8 :: A31
		real*8 :: P_HeITR
		real*8 :: q13
		real*8 :: q31a
		real*8 :: q31b
		real*8 :: Q31
		real*8 :: xheiTR_old
		! He <-> H charge-exchange rate coefficients (Huang Table 4 group B),
		! set per cell at the advection call site (see he_h_cx_rates).
		real*8 :: kcx_He0_Hp
		real*8 :: kcx_Hep_H0
		! Metal electrons of the cell, counted per H nucleus (n_e,metal/n_h).
		! The recombination terms of the advection residuals scale with the
		! TOTAL free electron density, and in the shielded base the metals are
		! the dominant electron donors (n_e,metal/n_e = 0.4-1.1 there), so the
		! residuals add this to the H/He electrons. Same definition as the
		! equilibrium residual's metal_electron_sum (X+ once, X++ twice) and as
		! calc_ne; zero when metals are absent.
		real*8 :: xe_metal = 0.0d0
	end type adv_rates

	type(adv_rates), save :: adv_cell
	!$omp threadprivate(adv_cell)

	! Named-field cell state for the post-process temperature residual
	! (T_equation / solve_T_brent). Filled by post_process_adv before each
	! cell temperature solve.
	type teq_state
		real*8 :: nhi
		real*8 :: nhii
		real*8 :: nhei
		real*8 :: nheii
		real*8 :: nheiii
		real*8 :: mup
		real*8 :: mum
		real*8 :: rhov
		real*8 :: coeff
		real*8 :: dr
		real*8 :: Told
		real*8 :: heaold
		! n(H2)/(n_tot + n_e) of the cell, for the caloric EOS of the energy
		! residual. Zero for an atomic gas, which is the monatomic case.
		real*8 :: x_h2
		! The same fraction in the UPSTREAM cell, the one the upwind
		! difference takes its energy from. It selects the branch of the
		! residual together with x_h2 above: a cell is solved with the
		! caloric energy whenever either side of the interface holds H2,
		! because an atomic cell fed by a molecular neighbor receives
		! molecular energy. Zero for an atomic upstream.
		real*8 :: x_h2_up = 0.0d0
		! SPECIFIC internal energy of that upstream cell,
		! e_up = E(x_h2_up, T_up)/mu_up [code units, energy per unit mass],
		! with E the energy per particle of the caloric EOS. This is the
		! quantity the flow transports into the cell, so it is evaluated at
		! the UPSTREAM composition and the upstream temperature: the
		! rovibrational heat capacity of H2 belongs to the gas that holds
		! the molecules, and across a dissociation front the two cells store
		! different energy at the same temperature. Formed by the caller,
		! and held fixed while the root finder varies this cell's T.
		real*8 :: e_up = 0.0d0
		! Divergence of the mass flux of the cell, div(rho v) [code units],
		! the coefficient of the enthalpy flux term of the steady
		! internal-energy equation (T_equation, where the algebra is written
		! out). It is a property of the state the residual is handed, not of
		! the trial temperature, so it is held fixed while the root finder
		! varies T, exactly as the densities are.
		!
		! Zero is the stationary mass flux rho v r^2 = const, for which the
		! term is absent from the equation; the field then contributes an
		! exact zero and the residual is the advected balance without it.
		real*8 :: div_rhov = 0.0d0
	end type teq_state

	type(teq_state), save :: teq_cell
	!$omp threadprivate(teq_cell)

	end module ion_cell_state
