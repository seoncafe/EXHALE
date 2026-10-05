	module ion_cell_state
	! Named-field cell state for the equilibrium ionization residuals. It
	! replaces reading the position-packed params() inside the residuals.
	! params stays as the MINPACK transport argument
	! (hybrd1/fdjac1/solve_ieq), but
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
		! 1^1S -> 2^3S by electron impact [cm^3 s^-1]: the direct excitation
		! plus the excitation of the higher triplets, which cascade into
		! 2^3S (utils_ion_eq: HeITR_coeffs)
		real*8 :: q13 = 0.0d0
		! 2^3S -> 1^1S electron-impact de-excitation [cm^3 s^-1], the
		! detailed-balance reverse of the direct excitation (Cool_coeff.f90)
		real*8 :: q31g = 0.0d0
		real*8 :: q31a = 0.0d0
		! 2^3S -> 2^1P and every singlet above it (utils_ion_eq: HeITR_coeffs)
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
		! kcx_He0_Hp = He0+H+ -> He++H0, kcx_Hep_H0 = He++H0 -> He0+H+,
		! and kcx_Hepp_H0 = He2+ + H0 -> He+ + H+ (radiative).
		real*8 :: kcx_He0_Hp = 0.0d0
		real*8 :: kcx_Hep_H0 = 0.0d0
		real*8 :: kcx_Hepp_H0 = 0.0d0
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
		! The cell this context describes (1..N; 0 = not a grid cell), for
		! diagnostics written from inside a parallel sweep.
		integer :: jcell = 0
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
		!                carriers (diffusive_photochemistry),
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
		!   x_heii_fixed, x_heiii_fixed
		!                the two ionized stages of helium, owned by the
		!                same transported partition and measured per
		!                helium NUCLEUS:
		!                x_heii_fix  = n_He+  / n_He(nuclei)  -> x(2)
		!                x_heiii_fix = n_He++ / n_He(nuclei)  -> x(3)
		!                They carry their own flags rather than riding on
		!                x_hp_fixed because a run may transport the proton
		!                while the helium stages stay local: which stages
		!                are carried is the set the transport operator
		!                registered, and a flag per stage is what states
		!                that set to the residual.
		!   x_hetr_fixed the He 2^3S population, owned by the transported
		!                level (He 2^3S transport: True), per helium
		!                NUCLEUS:
		!                x_hetr_fix = n(He 2^3S) / n_He(nuclei) -> x(x_hetr_row)
		!                x_hetr_row is the row of the level balance in the
		!                system the run solves (4 atomic, 8 molecular;
		!                ieq_triplet_row), set with the flag, because this
		!                module and ion_residual_core cannot ask the sweep.
		!
		!                Rows 2 and 3 exist in every system solved where
		!                helium is in the mixture; the one system that
		!                offers row 1 alone (System_H) is solved only where
		!                there is no helium, which "Ionization transport"
		!                refuses.
		logical :: x_h2_fixed = .false.
		logical :: x_ox_fixed = .false.
		logical :: x_hp_fixed = .false.
		logical :: x_heii_fixed = .false.
		logical :: x_heiii_fixed = .false.
		real*8 :: x_h2_fix = 0.0d0
		real*8 :: x_oh_fix = 0.0d0
		real*8 :: x_h2o_fix = 0.0d0
		real*8 :: x_hp_fix = 0.0d0
		real*8 :: x_heii_fix = 0.0d0
		real*8 :: x_heiii_fix = 0.0d0
		logical :: x_hetr_fixed = .false.
		integer :: x_hetr_row = 0
		real*8 :: x_hetr_fix = 0.0d0
		! THE H2 PARTITION OF THE NON-IONIZED HYDROGEN, stated by the
		! lower-boundary reservoir (the base handoff) for a lower ghost.
		! The handoff gives x2 = 2 n(H2)/(2 n(H2) + n(H I)), the fraction
		! of the NEUTRAL hydrogen nuclei bound in H2 (the lower-atmosphere
		! model never saw the wind's ionizing field); the ionized fraction
		! is the ghost's own. The row the H2 balance (4) is replaced by is
		! therefore
		!     x(4) - x2 (1 - x_ion(x)) = 0,
		!     x_ion = x(1) + x(5) + x(6) + x(7)
		! (h2_fraction_of_reservoir), solved TOGETHER with the ionization
		! rows: x_ion is an unknown of the same solve, not the value of a
		! previous pass. Pinning x(4) at x2 (1 - x_ion) of a previous
		! state leaves 1 - x(4) of the hydrogen for H I and the ions, and
		! where the ionization balance asks for more than that room no
		! composition satisfies the rows (md/Update_EXHALE_stage3.md
		! section 105: 0.0315 of ions against 0.0302 of room in the
		! LHS 1140 b He/H 2e4 ghost). In the implicit row the neutral
		! hydrogen left is (1 - x2)(1 - x_ion) >= 0 at every x_ion, so the
		! row never empties the simplex. Exclusive with x_h2_fixed.
		logical :: x_h2_neutral_partition_fixed = .false.
		real*8 :: x_h2_neutral_partition = 0.0d0
	end type ion_rates

	type(ion_rates), save :: ieq_cell
	!$omp threadprivate(ieq_cell)

	! Named-field cell state for the advection-correction ionization residuals
	! (System_implicit_adv_H/HeH/HeH_TR). One type covers BOTH advection
	! layouts: the H-only and HeH systems read the leading fields, the HeH_TR
	! system additionally reads the triplet channels (rcheiTR..xheiTR_hist).
	! heh_loc is the effective He/H for the electron density (global HeH
	! normally, the local nhe/nh with He_diffusion). Filled field by field by
	! post_process_adv before each hybrd1 advection solve.
	!
	! THE STEP. Each residual is one step of the backward differentiation
	! formula for dx/dr = R(x)/v in the fractions x (post_process_adv,
	! variable_step_bdf2_weights):
	!
	!     x_j - x_hist = c1 R_j(x_j) ,   c1 = g h_j / v_el,j ,
	!
	! with x_hist = a1 x_{j-1} - a2 x_{j-2} and h_j = r_j - r_{j-1}. The
	! second-order step has the variable-step BDF2 weights (a1, a2, g); the
	! first step of the recursion, and a step retaken for positivity, is
	! backward Euler, a1 = g = 1 and a2 = 0, for which x_hist is the
	! upstream fraction x_{j-1}. The rates R_j and the velocity are both
	! those of the cell the step lands on.
	!
	! THE VELOCITY IS THAT OF THE ELEMENT'S NUCLEI. x is a fraction of the
	! nuclei of one element, and those nuclei move with their own velocity
	! v_el = v + w_el, w_el the diffusive (gradient, eddy, settling) drift
	! of the element against the mass-weighted velocity v: with the element
	! continuity div(n_el v_el) = 0 the steady stage balance
	! div(x n_el v_el) = S becomes v_el dx/dr = S/n_el. c1 weights the
	! hydrogen row with v_H and c1_he the helium rows with v_He; without
	! element diffusion both velocities are v.
	!
	! xheiS_hist is the history of the GROUND SINGLET He(1^1S) alone, not of
	! the summed He I: the advection systems carry the singlet and the
	! He(2^3S) metastable as separate unknowns, so that neither is ever formed
	! as the difference of the other two (a difference that collapses to zero
	! once the metastable holds most of the neutral He). Without the triplet
	! the two coincide, all He I being in the singlet.
	type adv_rates
		real*8 :: c1      ! g h_j/v_H, the rate weight of the hydrogen row
		real*8 :: c1_he   ! g h_j/v_He, the rate weight of the helium rows
		real*8 :: xhi_hist
		real*8 :: xheiS_hist
		real*8 :: xheiii_hist
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
		real*8 :: q31g
		real*8 :: q31a
		real*8 :: q31b
		real*8 :: Q31
		real*8 :: xheiTR_hist
		! He <-> H charge-exchange rate coefficients (Huang Table 4 group B),
		! set per cell at the advection call site (see he_h_cx_rates).
		real*8 :: kcx_He0_Hp
		real*8 :: kcx_Hep_H0
		real*8 :: kcx_Hepp_H0
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
		! Neutral helium by level: the ground singlet n(1^1S) and the 2^3S
		! metastable (zero where it is not tracked). The cooling acts on each
		! level with its own coefficients (radiative_cooling_of_cell).
		real*8 :: nheiS
		real*8 :: nheiTR = 0.0d0
		real*8 :: nheii
		real*8 :: nheiii
		! Ground-capture escape weights of the H II, He II and He III
		! recombinations of the cell (utils_ion_eq:
		! ground_capture_escape_weights), functions of the densities, held
		! fixed while the root finder varies T.
		real*8 :: y_HI   = 0.0d0
		real*8 :: y_gnd  = 0.0d0
		real*8 :: y_HeII = 0.0d0
		! Mean molecular weight of the cell (mup) and of its upstream
		! neighbor (mum); mup*mum is the factor the residual is multiplied by
		! (T_equation), nothing else.
		real*8 :: mup
		real*8 :: mum
		! rho_j v_j of the cell [code units]
		real*8 :: rhov
		! Pressure-work coefficient mum*v_j*(rho_j - rho_hist), with rho_hist
		! the history of the density formed with the weights of the step
		! (the adv_rates note on THE STEP).
		real*8 :: coeff
		! The width the right-hand side of the step is multiplied by: g*h_j
		! for the variable-step BDF2 step, h_j for backward Euler, with
		! h_j = r_j - r_{j-1} [code units].
		real*8 :: dr_step
		! Heating of the cell [code units]
		real*8 :: heaold
		! n(H2)/(n_tot + n_e) of the cell, for the caloric EOS of the energy
		! residual. Zero for an atomic gas.
		real*8 :: x_h2
		! THE HISTORY OF THE SPECIFIC INTERNAL ENERGY the flow carries in,
		! e_hist = a1 e_{j-1} - a2 e_{j-2} with the weights of the step,
		! each e_k = E(x_H2,k, T_k)/mu_k [code units, energy per unit mass]
		! evaluated at the composition and temperature of ITS OWN cell: the
		! rovibrational heat capacity of H2 belongs to the gas that holds the
		! molecules, so across a dissociation front two cells at the same
		! temperature store different energy, and an atomic cell fed by a
		! molecular neighbor receives molecular energy. Backward Euler is the
		! upstream cell's energy alone. Formed by the caller and held fixed
		! while the root finder varies this cell's T.
		real*8 :: e_hist = 0.0d0
		! Divergence of the mass flux of the cell, div(rho v) [code units],
		! the coefficient of the enthalpy flux term of the steady
		! internal-energy equation (T_equation, where the algebra is written
		! out). post_process_adv takes it from the mass row of the state, the
		! face mass fluxes of the Riemann solve differenced over the cell
		! (the statement in its energy block). It is a property of the state
		! the residual is handed, not of the trial temperature, so it is held
		! fixed while the root finder varies T, exactly as the densities are.
		!
		! Zero is a stationary mass flux, for which the term is absent from
		! the equation; the field then contributes an exact zero and the
		! residual is the advected balance without it.
		real*8 :: div_rhov = 0.0d0
		! THE TRANSPORT TERMS OF THE ENERGY EQUATION the run solved and the
		! advected balance must keep (post_process_adv, 2026-09-30): heat
		! conduction Q_c = b_lo T_{j-1} + b_di T_j + b_up T_{j+1} (the
		! triplets of viscous_conduction) and the divergence of the enthalpy
		! flux of the element and carrier fluxes, entered as heating [code
		! units]. heat_extra is the part held fixed while the root finder
		! varies this cell's T (the neighbors' conduction and the enthalpy
		! flux divergence), cond_diag the coefficient b_di of this cell's own
		! T. Zero for a run without these terms, and the residual is then the
		! one it always was, to the bit.
		real*8 :: heat_extra = 0.0d0
		real*8 :: cond_diag  = 0.0d0
	end type teq_state

	type(teq_state), save :: teq_cell
	!$omp threadprivate(teq_cell)

	contains

	! Fraction of a cell's hydrogen nuclei that is ionized, in the
	! molecular layout: H+ (x(1)) and the hydrogen nuclei held in the
	! molecular ions, two in H2+ (x(5)), three in H3+ (x(6)) and one in
	! HeH+ (x(7)), the multiplicities those unknowns already carry.
	pure double precision function ionized_hydrogen_nuclei_fraction(x)     &
	                                                          result(x_ion)
	real*8, intent(in) :: x(:)
	x_ion = x(1) + x(5) + x(6) + x(7)
	end function ionized_hydrogen_nuclei_fraction

	! The H2 row of a reservoir cell whose neutral hydrogen partition x2
	! is stated (x_h2_neutral_partition_fixed): x(4) = x2 (1 - x_ion),
	! the H2 nuclei fraction of ALL the hydrogen when x2 of the
	! non-ionized hydrogen is bound in H2 and x_ion is ionized.
	pure double precision function h2_fraction_of_reservoir(x2, x_ion)     &
	                                                         result(x_h2)
	real*8, intent(in) :: x2, x_ion
	x_h2 = x2*(1.0d0 - x_ion)
	end function h2_fraction_of_reservoir

	end module ion_cell_state
