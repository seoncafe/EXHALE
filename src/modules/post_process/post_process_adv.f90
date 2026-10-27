	module post_processing
	! Advective correction of H/He and trace-metal ionization profiles.
	! The _adv state omits H2, H2+, H3+, and HeH+: those densities are not
	! inputs or re-solved here. calc_ne, calc_ntot, and
	! hydrogen_helium_nuclei_density therefore omit molecular particles
	! and electrons, and nh=nhi+nhii counts only free H nuclei. Results
	! describe atomic/ionized escape flow; molecular-base and metal profiles
	! inside a molecular layer inherit this approximation.

	use global_parameters
	use ion_cell_state, only: ieq_cell, adv_cell, teq_cell
	use species_table, only: n_mion, n_melem, melem_i0, melem_top, mion_stage
	use utils
	use System_implicit_adv_H
	use System_implicit_adv_HeH
	use System_implicit_adv_HeH_TR
	use System_HeH_metals, only: ion_system_HeH_metals, set_metal_coeffs
	use charge_exchange,   only: cx_set_cell, he_h_cx_rates
	use Cooling_Coefficients
	use utils_ion_eq
	use composition, only: he_ground_singlet_density
	use output_write
	use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
	use caloric_eos, only: h2_particle_fraction, internal_energy_of_mixture
	use equation_T
	use opacity_models           ! opacity_pT_factor for the 'P' model
	! The mass row of the state and the scale that makes it dimensionless:
	! the SAME operator the stationary certification measures that row by,
	! called here on the state the post-process is handed rather than
	! re-derived from it.
	use Conversion,          only: W_to_U
	use Reconstruction_step, only: n_faces_positivity_limited
	use steady_residual_mod, only: assemble_residual, residual_row_scale,  &
	                               row_terms_describe_state,               &
	                               carrier_enthalpy_divergence
	use certification,       only: cert_scale_floor, cert_tol_energy
	! The species the reconstruction below omits, for the cell class it does
	! not model (see adv_unsupported).
	use ionization_equilibrium, only: nmol_eq, nox_eq
	! The molecular channels that destroy the H/He populations this
	! post-process solves, for the composition class (see where it is formed).
	use mol_rates, only: rk_R10_Hp_H2v4, rk_R13_Hp_H2_M,                &
	                     rk_R17_Hep_H2_diss, rk_R23_H2_Hep_cx
	! The transport terms of the run's energy equation (heat conduction and
	! the enthalpy flux of the element and carrier fluxes), which the
	! advected balance keeps with the corrected temperature (the energy
	! block below, "THE TRANSPORT TERMS").
	use viscous_conduction, only: conduction_active,                      &
	                              thermal_conduction_coeffs,              &
	                              conduction_base_level_T
	use binary_element_diffusion, only: interdiffusion_enthalpy_active,   &
	                              interdiffusion_enthalpy_divergence_of_state, &
	                              helium_diffusive_face_flux,             &
	                              element_mass_fractions, mixture_mass_sum

	implicit none

	! Validity limits of the advection correction (see post_process_adv, where
	! both are used and their derivation is written out).
	!
	! Da_local_equilibrium uses the smallest listed population-rate estimate.
	! Fast relaxation permits the local-equilibrium approximation when that
	! equilibrium varies slowly over the relaxation length. This estimate is
	! not an eigenvalue of the coupled reaction Jacobian.
	real*8, parameter :: Da_local_equilibrium = 1.0d2
	! xHII_adv_min avoids extracting trace ions by subtracting x_HI from 1.
	! sqrt(eps)=1.5e-8 is the MINPACK step tolerance and relative Jacobian
	! perturbation, not a bound on root accuracy; the measured loss of trace
	! ions motivating this guard is recorded in condition (iii) below.
	real*8, parameter :: xHII_adv_min = 1.0d-6
	! Correct a cell only if its production mass residual satisfies
	! |R_1(j)|/s_1(j) <= adv_conditional_tol (steady_residual.f90,
	! mass_flux_row_scale; output_write). This is the fractional face
	! mass-flux change across the cell and bounds the conditional correction
	! error. Whole-state certification uses cert_tol_mass separately and is
	! recorded in the file header. enthalpy_ratio_report_level=1 reports
	! when mass-divergence enthalpy exceeds the retained energy terms; this
	! is a sensitivity measure, not a stationarity certificate.
	real*8, parameter :: enthalpy_ratio_report_level = 1.0d0

	! Integrate ionization and internal-energy ODEs along r with
	! variable-step BDF2, using all rates and velocities at r_j:
	! y_j-a1 y_{j-1}+a2 y_{j-2}=g h_j F(y_j,r_j),
	! h_j=r_j-r_{j-1}, w=h_j/h_{j-1},
	! a1=(1+w)^2/(1+2w), a2=w^2/(1+2w), g=(1+w)/(1+2w).
	! These are the derivative weights of the quadratic through three
	! nodes (Li & Liao 2022, arXiv:2201.00527, eq. 1.4; Hairer,
	! Norsett & Wanner, Solving Ordinary Differential Equations I,
	! 2nd ed., sect. III.5, not re-read here). At w=1 they give
	! (3/2)y_j-2y_{j-1}+(1/2)y_{j-2}=h F_j. The local error is O(h^3).
	!
	! Start or restart after a refused cell with backward Euler. Retake a
	! BDF2 step with inadmissible populations or no temperature root by
	! backward Euler. Its local error is O(h^2); both methods settle stiff
	! relaxation, whereas trapezoid has amplification
	! (1+z/2)/(1-z/2) -> -1 for z=-Da -> -infinity.
	! Variable-step BDF2 uses w <= 1+sqrt(2) for zero stability
	! (Grigorieff 1983, Numer. Math. 42, 359, cited by Li & Liao;
	! the 1983 paper was not read);
	! larger ratios use backward Euler. Measured LHS 1140 b Mixed grids
	! have w <= 1.018 (2026-09-27).
	real*8, parameter :: bdf2_step_ratio_limit = 2.414213562373095d0

	! Advection-corrected H/He ionized fractions for the current cell, pinned
	! while the metal re-solve (pp_metals=2) adjusts only the metal stages.
	! Set per cell before each ion_system_metals_pp / hybrd1 call.
	real*8, save :: pp_xHII_fix   = 0.0d0
	real*8, save :: pp_xHeII_fix  = 0.0d0
	real*8, save :: pp_xHeIII_fix = 0.0d0

	! MEASUREMENT ONLY, default off; read once, at the first call of
	! post_process_adv (adv_column_controls_read).
	! EXHALE_ADV_COLUMN_FLOOR_IT=<n>: after the column iteration has met its
	! step test, take n further Newton steps before the final rows are
	! formed, so that the residual the log reports is the stagnation level
	! of the iteration (the rounding floor the margin of atol_E is anchored
	! on) and not the level the step test happened to stop at.
	! EXHALE_ADV_COLUMN_EXCLUDE_UNROOTED=1: the cells whose cell-by-cell
	! energy step found no admissible root in the marching sweep are fixed
	! values of the column solve (at the run's own temperature, written
	! adv_failed) instead of unknowns of it, which is how the column was
	! built until 2026-10-01 (energy_column_newton says why they are
	! unknowns by default).
	integer, save :: adv_column_floor_it = 0
	logical, save :: adv_column_exclude_unrooted = .false.
	logical, save :: adv_column_controls_read = .false.

	contains

	! Thermal Damkohler number Da=t_cross |heat-cool|/u_th, with
	! t_cross in s, net rate in erg cm^-3 s^-1, and u_th in erg cm^-3.
	! Da compares the net thermal-energy change over one crossing with
	! the cell's current energy. The code retains the run's temperature
	! for Da > 1 and permits an advected correction for Da <= 1.
	! Near heat=cool, |heat-cool| may be small even when the local
	! relaxation rate is large; Da does not measure that rate. As v -> 0
	! with nonzero net rate, t_cross grows and the correction switches
	! off before its outward-flow ODE loses meaning.
	! For u_th <= 0, return a large Da because no thermal state can be
	! advanced from that cell.
	pure real*8 function thermal_damkohler_number(t_cross, q_rad, u_th)
	real*8, intent(in) :: t_cross, q_rad, u_th

	! No internal energy to rewrite is the Da -> infinity end of the same
	! statement, and it is not a temperature the correction can start from.
	if (.not. (u_th > 0.0d0)) then
		thermal_damkohler_number = huge(1.0d0)
	else
		thermal_damkohler_number = t_cross*q_rad/u_th
	endif

	end function thermal_damkohler_number

	! Diagnostic ratio of the mass-divergence enthalpy term to retained
	! terms in rho v de/dr-p v dln(rho)/dr+h div(rho v)=heat-cool.
	! With derivatives from the BDF2 or backward-Euler history divided by
	! dr_step, use
	! q_adv=|rho v(e_j-e_hist)|,
	! q_prs=|(p/rho) v(rho_j-rho_hist)|,
	! q_enth=|h div(rho v) dr_step|, h=e+p/rho.
	! The ratio q_enth/max(q_adv,q_prs,1e-99) exceeds one when the omitted mass
	! divergence can set the corrected temperature. It is a sensitivity
	! diagnostic; 1e-99 makes the all-zero case zero, not 0/0. A small
	! ratio alone cannot certify a stationary mass row;
	! that decision uses |R_1|/s_1 in post_process_adv. The report uses
	! enthalpy_ratio_report_level=1. Arguments are the cell and its step
	! history.
	pure real*8 function enthalpy_flux_term_ratio(h_up, div_rhov, dr_step, &
	                     rho_up, v_up, e_up, e_hist, w_up, rho_hist)
	real*8, intent(in) :: h_up      ! enthalpy per unit mass of the cell
	real*8, intent(in) :: div_rhov  ! div(rho v), in code density/time units
	real*8, intent(in) :: dr_step   ! width the step multiplies F by
	real*8, intent(in) :: rho_up    ! density of the cell
	real*8, intent(in) :: v_up      ! velocity of the cell
	real*8, intent(in) :: e_up      ! internal energy per unit mass, the cell
	real*8, intent(in) :: e_hist    ! ... its history
	real*8, intent(in) :: w_up      ! p/rho of the cell
	real*8, intent(in) :: rho_hist  ! history of the density
	real*8 :: q_enth,q_adv,q_prs

	q_adv  = abs(rho_up*v_up*(e_up - e_hist))
	q_prs  = abs(w_up*v_up*(rho_up - rho_hist))
	q_enth = abs(h_up*div_rhov*dr_step)
	! A cell in which all three terms vanish carries no equation at all, and
	! the floor makes that ratio zero rather than 0/0: nothing to refuse.
	enthalpy_flux_term_ratio = q_enth/max(q_adv, q_prs, 1.0d-99)

	end function enthalpy_flux_term_ratio

	! The weights of one variable-step BDF2 step (see bdf2_step_ratio_limit
	! for the formula and its derivation): the step into r_j from r_{j-1} and
	! r_{j-2}, with h_step = r_j - r_{j-1} and h_prev = r_{j-1} - r_{j-2} in
	! any common unit.
	pure subroutine variable_step_bdf2_weights(h_step, h_prev,            &
	                                          w_prev, w_prev2, w_rate)
	real*8, intent(in)  :: h_step, h_prev
	real*8, intent(out) :: w_prev   ! a1, the weight of y_{j-1}
	real*8, intent(out) :: w_prev2  ! a2, the weight of y_{j-2}
	real*8, intent(out) :: w_rate   ! g, F_j is multiplied by g h_step
	real*8 :: step_ratio

	step_ratio = h_step/h_prev
	w_prev  = (1.0d0 + step_ratio)**2/(1.0d0 + 2.0d0*step_ratio)
	w_prev2 = step_ratio**2/(1.0d0 + 2.0d0*step_ratio)
	w_rate  = (1.0d0 + step_ratio)/(1.0d0 + 2.0d0*step_ratio)

	end subroutine variable_step_bdf2_weights

	subroutine post_process_adv(rho,v,p,T_in,heat,cool,eta,   &
                                  nhi_in,nhii_in,		    &
                                  nhei_in,nheii_in,nheiii_in,   &
                                  nheiTR_in, nm_in, f_sp_in,       &
                                  derived_state)


	real*8, dimension(1-Ng:N+Ng), intent(in) :: rho,v,p,T_in
   real*8, dimension(1-Ng:N+Ng), intent(in) :: heat,cool
   real*8, dimension(1-Ng:N+Ng), intent(in) :: eta
   real*8, dimension(1-Ng:N+Ng), intent(in) :: nhi_in,nhii_in
   real*8, dimension(1-Ng:N+Ng), intent(in) :: nhei_in,nheii_in,   &
      							  				        nheiii_in,nheiTR_in
   ! Converged equilibrium metal densities (dimensionless, n0 units), used
   ! by the metal-aware post-process modes (pp_metal_mode = 1 frozen, 2 re-solve).
   real*8, dimension(1-Ng:N+Ng,n_mion), intent(in) :: nm_in
   ! The composition of the recorded state, handed to the residual assembly
   ! below: its energy row reads the composition (the interdiffusion
   ! enthalpy flux of a diffusing mixture). Only the mass row is read here.
   real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp_in
   ! How the solves that produced the written rows ended, as a whole
   ! (adv_derived_state, output_write): handed to the caller and written
   ! into the header of both _adv files.
   type(adv_derived_state), intent(out) :: derived_state
	
	integer j,k
	integer :: n_pp_reject       ! cell-by-cell T solves rejected as non-physical
	! Cell solves that did not converge (hybrd1 info /= 1) and therefore kept
	! the incoming equilibrium state: the advection ionization system, the
	! metal stage re-solve, and the energy equation. Summed over the passes,
	! not over the cells of the last one: a cell reverted in an early pass
	! feeds the upwind cascade of that pass and has to be reported even if the
	! last pass converges everywhere.
	integer :: n_adv_noconv, n_metal_noconv, n_T_noconv
	! Energy cell solves whose iterate was kept although the iteration
	! reported no further progress, because its residual was already at
	! the cancellation floor of the terms of the equation (see the solve).
	integer :: n_T_res_root
	integer :: Neq_adv,lwa_adv   ! advection system size (metal-independent)
	integer :: Neq_mpp,lwa_mpp   ! metal re-solve system size (pp_metals=2)
	! Passes of the post-process (the outer loop below).
	integer, parameter :: n_pp_passes = 10
	 
	real*8, dimension(1-Ng:N+Ng) ::  T_K,p_out,T_out     ! T_K is in K; p_out and T_out use code units.
	real*8, dimension(1-Ng:N+Ng) ::  nh,nhe,ne,n_tot
	! Metal electron density [cgs]: the metal part of ne (X+ once, X++ twice),
	! handed to the advection residuals as adv_cell%xe_metal. Zero when metals
	! are absent or the post-process runs metal-free (nm_w = 0).
	real*8, dimension(1-Ng:N+Ng) ::  ne_metal
	! Ionized fraction of the H+He nuclei, for the SvS85 secondary ionization.
	real*8, dimension(1-Ng:N+Ng) ::  xion
	
	! Discard scratch: distinct locals for the discarded intent(out) slots so no
	! two out-arguments in the same call alias one another.
	real*8, dimension(1-Ng:N+Ng) ::  dum_v1,dum_v2,dum_v3,dum_v4,dum_v5,dum_v6

   ! Photo ionization rates
   real*8, dimension(1-Ng:N+Ng) ::  P_HI,P_HeI,P_HeII,P_HeITR
   ! Metal photoionization rates for each ion (filled by PH_heat_HHe; used for the
   ! metal re-solve in mode 2, otherwise discarded).
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  P_m
   ! The part of P_m that takes a neutral straight to X++ (an autoionizing
   ! inner-shell vacancy), for the same re-solve.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  P_m2
   ! Working metal densities driving the post-process heating/cooling [cgs].
   ! Built from nm_in per pp_metal_mode: 0 -> zero (metal-free, legacy _adv),
   ! 1 -> frozen eq metals, 2 -> re-solved. nm_out is its dimensionless (n0)
   ! copy written to the _adv ion-species file.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  nm_w, nm_out
   ! Line transfer of the ground-term fine-structure lines (escape
   ! probabilities and the incident lower-atmosphere field), frozen at the
   ! profile the temperature solve starts from (see equation_T pp_beta_fs).
   real*8, dimension(1-Ng:N+Ng,n_fsline) ::  beta_fs_pp, nbar_fs_pp

   ! Recombination coefficients
   real*8, dimension(1-Ng:N+Ng) ::  rchiiB,rcheiiB,rcheiiiB,rcheiTR
   ! Recombination radiation absorbed on the spot (use_h_rec_escape,
   ! use_he_rec_coupling; zero-effect when both are off).
   real*8, dimension(1-Ng:N+Ng) ::  rchiiB_hrc,rcheiiB_hrc,rcheiiiB_hrc
   real*8, dimension(1-Ng:N+Ng) ::  dP_HI_hrc,dP_HeI_hrc,dheat_hrc
   ! H2 density seen by the He-recombination coupling. The _adv reconstruction
   ! is molecule-free (module-header composition note), so it is identically
   ! zero here and the coupling reduces to the H I / He I competition; the H2
   ! rate it returns is discarded for the same reason.
   real*8, dimension(1-Ng:N+Ng) ::  nh2_pp, dP_H2_hrc
   ! Metal share of the recombination photons, added to P_m, and its
   ! part that ends in X++, added to P_m2
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  dP_m_hrc, dP_m2_hrc
   ! Metal recombination/ionization rates for each ion returned by eval_cool.
   ! In the re-solve mode (pp_metals=2) they feed the cell-by-cell metal
   ! ionization-balance solve; in the frozen mode they are discarded.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  rec_m_pp,aion_m_pp
   ! Each element's total metal density (sum of stages, held fixed across the
   ! re-solve since ionization only redistributes among an element's stages).
   real*8, dimension(1-Ng:N+Ng,n_melem) ::  nm_tot_pp
   ! Each element's metal coefficients handed to set_metal_coeffs for the
   ! re-solve (canonical element order); built per cell from the 2-D rates.
   real*8, dimension(n_melem) ::  meg_ntot,meg_g0,meg_g1,meg_b0,meg_b1,  &
                                  meg_a1,meg_a2,meg_g02
   integer, dimension(n_melem) ::  meg_top
   integer :: i0,top,im
      
   real*8, dimension(1-Ng:N+Ng) ::  q13,q31g,q31a,q31b,Q31
   ! Ground-capture escape weights of the H II, He II and He III
   ! recombinations of each cell (ground_capture_escape_weights), for the
   ! cooling the temperature root balances
   real*8, dimension(1-Ng:N+Ng) ::  y_HI_pp, y_gnd_pp, y_HeII_pp
	real*8 :: A31
 	
 	! Ionization coefficients
   real*8, dimension(1-Ng:N+Ng) ::  a_ion_HI,a_ion_HeI,a_ion_HeII,a_ion_HeITR
      
   ! Heating, cooling
   real*8, dimension(1-Ng:N+Ng) ::  theat,tcool
   ! Cooling of the composition each pass STARTS from, at that pass's
   ! temperature [erg cm^-3 s^-1], from the first eval_cool of the pass. It is
   ! the radiative side of the thermal Damkohler number of the energy solve
   ! (see the loop), which has to be formed before the temperature is solved
   ! for and therefore cannot use tcool.
   real*8, dimension(1-Ng:N+Ng) ::  tcool_in
   ! Photoheating of ONE particle of each absorber [erg s^-1] from the
   ! attenuated field, the rate side of the heating assembly.
   real*8, dimension(1-Ng:N+Ng) ::  h1_HI_pp,h1_HeI_pp,h1_HeII_pp,       &
                                    h1_HeTR_pp,h1_H2_pp
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  h1_m_pp
   ! Molecular and oxygen carrier densities of the _adv composition, and the
   ! Lyman-Werner and FUV band rates: all identically zero or, for the
   ! dissociations per pump, unity, because the reconstruction carries no
   ! molecules (module header). The two .false. flags of the heating call
   ! are what actually keeps those channels out; these arrays exist because
   ! the assembly takes the composition as arguments.
   real*8, dimension(1-Ng:N+Ng,4) ::  nmol_pp
   real*8, dimension(1-Ng:N+Ng,3) ::  nox_pp
   real*8, dimension(1-Ng:N+Ng) ::  k_lw_pp, p_lw_pp, k_co_pp
   real*8, dimension(1-Ng:N+Ng,n_fuv_band) ::  j_fuv_pp
   ! The heating deposits of the _adv composition, channel by channel. The
   ! _adv files carry the total; the columns are kept because the assembly
   ! forms the total as their sum.
   real*8, dimension(1-Ng:N+Ng,n_heat_channel) ::  heat_chan_pp
      
 	! Updated species densities
   real*8, dimension(1-Ng:N+Ng) :: nhi,nhii
   ! Neutral helium is carried as its two populations, the ground singlet
   ! nheiS = n(1^1S) and the metastable nheiTR = n(2^3S); the summed He I
   ! density nhei is their SUM, formed wherever a routine wants the total.
   ! The advection systems solve for the two populations directly, so the
   ! singlet is never obtained as nhei - nheiTR, a difference that loses all
   ! its digits once the metastable holds most of the neutral helium.
   real*8, dimension(1-Ng:N+Ng) :: nhei,nheii,nheiii,nheiTR,nheiS
   real*8, dimension(1-Ng:N+Ng) :: mmw
   real*8, dimension(1-Ng:N+Ng) :: nhi_w,nhii_w
   real*8, dimension(1-Ng:N+Ng) :: nhei_w,nheii_w,nheiii_w,nheiTR_w
      							  
	      
	      
   real*8 :: dr	                           ! Grid spacing
   real*8 :: tol,dpmpar                      ! Equilibrium system setup
      
      
	! Validity of the advection correction, cell by cell (filled by the block
	! just before the ionization loop, where the three conditions are stated).
	logical, dimension(1-Ng:N+Ng) :: adv_correction_valid
	real*8  :: t_cross        ! residence time of the gas in the cell [s]
	real*8  :: t_cross_He     ! the same for the helium nuclei [s]
	real*8  :: nu_relax_He    ! slowest relaxation rate of the helium rows
	! THE VELOCITY OF EACH ELEMENT'S NUCLEI [code units], v_el = v + w_el
	! (see where they are formed): the hydrogen row of the advection
	! systems is stepped on v_nuc_H, the helium rows on v_nuc_He.
	real*8, dimension(1-Ng:N+Ng) :: v_nuc_H, v_nuc_He
	real*8, dimension(0:N) :: J_he_face
	real*8, dimension(1-Ng:N+Ng,1+n_melem) :: Y_elem
	real*8, dimension(1-Ng:N+Ng) :: msum_pp
	real*8  :: J_he_cell, rho_cgs_cell
	real*8  :: nu_relax       ! relaxation rate of the slowest species [1/s]
	real*8  :: rec_HeII_tot   ! total He II -> He I recombination coefficient
	real*8  :: Da_slowest     ! Damkohler number of that species
	real*8  :: xHII_eq        ! equilibrium H ionized fraction of the cell
	integer :: n_adv_eq       ! cells left at the equilibrium ionization

	! Validity of the energy correction, cell by cell (see the energy loop,
	! where both quantities are stated and formed).
	real*8  :: t_cross_E      ! residence time of the gas in the cell [s]
	real*8  :: q_rad          ! |heating - cooling| [erg cm^-3 s^-1]
	real*8  :: u_th           ! internal energy density [erg cm^-3]
	real*8  :: Da_thermal     ! thermal Damkohler number of the cell
	real*8  :: Da_thermal_max ! largest of those over the tested cells
	real*8  :: r_Da_thermal_max    ! radius at which it is largest
	integer :: n_T_local      ! cells left at the run's own temperature
	! The enthalpy flux of the mass-flux divergence, the third term of the
	! energy equation (see the statement before the energy loop), against
	! the other two terms: the specific internal energy of the cell, its
	! p/rho and its enthalpy per unit mass, for the sensitivity diagnostic.
	real*8  :: e_up           ! internal energy per unit mass of the cell
	real*8  :: w_up,h_up      ! p/rho and the enthalpy per unit mass, this cell
	real*8  :: enthalpy_flux_ratio, enthalpy_flux_ratio_max
	real*8  :: r_enthalpy_flux_max ! radius at which that ratio is largest
	integer :: n_enthalpy_term_dominant ! cells whose ratio is above the level

	! Whether the mass row of the cell is within the fraction a corrected
	! row is accurate to: .false. where it stands above
	! adv_conditional_tol, so the terms the steady equations drop are
	! larger than that fraction of the ones they keep (see the block before
	! the ionization validity conditions). Formed once, on the state the
	! post-process was handed, and read by every pass.
	logical, dimension(1-Ng:N+Ng) :: mass_row_within_tol
	! The measure itself, row by row, written as a column of
	! Hydro_ioniz_adv.txt: the CONDITION of each row travels with its
	! verdict. Zero in a row the loop below does not reach.
	real*8, dimension(1-Ng:N+Ng)   :: adv_mass_row
	! The state itself in conservative variables, its assembled residual, and
	! the particle count the residual is given: the three arguments of the
	! production mass operator.
	real*8, dimension(3,1-Ng:N+Ng) :: W_state, u_state, R_state
	real*8, dimension(1-Ng:N+Ng)   :: n_part_state
	real*8  :: mass_row          ! |R_1(j)|/s_1(j) of one cell
	real*8  :: mass_row_max      ! largest of those over the tested cells
	real*8  :: r_mass_row_max    ! radius at which it is largest
	integer :: n_limited_before  ! face-limiter ledger of the RUN, restored
	logical :: terms_are_run_state ! the stored row terms are this state's
	! Mean molecular weight of that same state, for the term sizes above.
	real*8, dimension(1-Ng:N+Ng) :: mmw_in
	integer :: n_not_stationary ! cells the stationarity condition refuses
	integer :: n_stationarity_only ! ... that the ionization conditions do not

	! Particle count of the species the reconstruction omits, and of those it
	! carries, for the cell class it does not model (adv_unsupported).
	real*8  :: n_omitted, n_carried
	logical, dimension(1-Ng:N+Ng) :: cell_class_modelled
	! The loss rates [s^-1] of one solved population by the molecular
	! channels the reconstruction omits and by the channels it carries, for
	! the composition class (formed with the cell class).
	real*8  :: loss_omitted, loss_carried, n_h2_cell
	! Cells whose particle count or the balance of one solved population is
	! carried more by the omitted species than by those the reconstruction
	! holds: neither the composition nor the temperature of such a row is a
	! statement about that gas, and both fields are unsupported.
	integer :: n_species_class_refused
	! The same statement about the ENERGY equation (see where it is formed):
	! the heating and cooling of the run's own state that this post-process
	! does not carry, against those it does, and the verdict of each cell.
	real*8  :: q_omitted, q_carried
	real*8, dimension(1-Ng:N+Ng)   :: theat_of_run_composition
	logical, dimension(1-Ng:N+Ng) :: energy_class_modelled
	integer :: n_energy_class_refused
	real*8  :: r_energy_class_top

	! Validity of the temperature and of the composition of each row, in the
	! five values of the two-field schema (output_write); written as the last
	! two columns of Hydro_ioniz_adv.txt and Ion_species_adv.txt. The values
	! are assigned where each decision is taken, so no reader of this routine
	! has to reconstruct which condition a row carries.
	integer, dimension(1-Ng:N+Ng) :: adv_T_status, adv_comp_status


   ! Legacy MINPACK argument sized for ion_system_metals_pp. The residuals
   ! read named cell-state fields; params is not consumed.
   real*8, dimension(60) :: params
   real*8, dimension(40) :: paramsT
      
	real*8 :: vp,vm
	real*8 :: sys_x_T(1)
   real*8 :: wa_T(8)
   ! The cell energy solve kept the run's own temperature: it did not
   ! converge, or its root was outside the band a temperature of this gas
   ! can occupy.
   logical :: T_solve_fell_back
   ! How the cell energy solve ended (energy_step): a root, a root at the
   ! cancellation floor of its own terms, or one of the four ways it keeps
   ! the run's own temperature.
   integer :: T_outcome
   integer, parameter :: T_root_found = 0, T_root_at_floor = 1,           &
                         T_no_bracket = 2, T_not_converged = 3,           &
                         T_non_positive = 4, T_out_of_band = 5

   ! THE STEP (bdf2_step_ratio_limit, this module). The variable-step BDF2
   ! weights of the step INTO each cell, from the two cells below it, and
   ! whether the grid admits that step there (zero-stability).
   real*8,  dimension(1-Ng:N+Ng) :: bdf_a1, bdf_a2, bdf_g
   logical, dimension(1-Ng:N+Ng) :: bdf2_admissible
   ! Whether the step in progress is the BDF2 step (else backward Euler).
   logical :: use_bdf2
   ! THE TRANSPORT TERMS OF THE ENERGY EQUATION (the statement at the
   ! energy loop): whether the run solved any, the conduction triplets and
   ! the enthalpy flux divergence of the current temperature profile, the
   ! neighbor temperatures the marching sweep reads, the step each corrected
   ! cell took, and whether the column Newton solve is the one evaluating.
   logical :: pp_transport_on, pp_newton_mode
   real*8, dimension(N) :: pp_cond_lo, pp_cond_di, pp_cond_up
   real*8, dimension(1-Ng:N+Ng) :: pp_heat_rel, pp_T_nbr, pp_T_march
   ! The two parts of pp_heat_rel, entered as heating: the enthalpy flux
   ! divergence of the element fluxes and of the carrier fluxes.
   real*8, dimension(1-Ng:N+Ng) :: pp_heat_elem, pp_heat_carr
   real*8  :: pp_T_bath
   logical, dimension(1-Ng:N+Ng) :: pp_bdf2_used
   integer :: n_col_newton_it
   real*8  :: col_newton_step
   ! How the column energy solve of the pass ended (adv_col_* of
   ! output_write), why it was rejected, the cell that triggered the
   ! rejection, the number of its unknowns, and the passes whose column was
   ! rejected.
   integer :: col_outcome, col_reason, col_cell, col_unknowns
   integer :: n_col_passes_rejected
   ! The column solve in progress has failed a test (reject_column), and
   ! the temperature of the triggering cell at that test [code units].
   logical :: col_bad
   real*8  :: col_T_trigger
   ! The terms of one cell's energy row (energy_row_terms_of_cell), and
   ! the seven physically grouped terms the acceptance of the column is
   ! measured against (energy_row_grouped_terms).
   integer, parameter :: n_row_terms = 14, n_row_groups = 7
   ! The worst energy row of the column of the pass (energy_column_newton):
   ! its cell, |R| [erg cm^-3 s^-1], |R|/S_E, and the tolerance over S_E.
   integer :: col_worst_cell
   real*8  :: col_worst_R, col_worst_rel, col_worst_tol_rel
   ! Physical cells of the column of the pass whose marching energy step
   ! found no root (energy_column_newton).
   integer :: n_col_unrooted
   ! Width the right-hand side of a step is multiplied by, g h_j or h_j
   ! [code units], and the histories of the density and of the specific
   ! internal energy that step differences against.
   real*8  :: dr_step, rho_hist, e_hist
   ! The divergence of the mass flux of each cell, div(rho v), from the
   ! mass row of the state (see the statement before the energy loop).
   real*8, dimension(1-Ng:N+Ng) :: div_rhov_state
   ! Ledger of the steps of the last pass: BDF2 and backward Euler steps,
   ! and the BDF2 steps retaken by backward Euler.
   integer :: n_bdf2_comp, n_be_comp, n_bdf2_comp_retaken
   integer :: n_bdf2_T, n_be_T, n_bdf2_T_retaken
      
      
      
      
   ! Numerical tolerance for system solution
   tol = sqrt(dpmpar(1))

   ! The two measurement keys of the column energy solve (see the
   ! declaration of adv_column_floor_it), read once.
   if (.not. adv_column_controls_read) call read_column_energy_solve_controls()

   ! Advection system size. The advection-correction systems
   ! (adv_implicit_H/HeH/HeH_TR) solve only H/He fractions and never
   ! the metal stages, so the count is independent of the (possibly
   ! larger) global N_eq used by the equilibrium solver. Using N_eq
   ! here would feed hybrd1 uninitialized fvec/x entries for the metal
   ! slots. The matching MINPACK workspace size is lwa_adv.
   if (.not.thereis_He) then
      Neq_adv = 1
   else if (thereis_HeITR) then
      Neq_adv = 4
   else
      Neq_adv = 3
   endif
   lwa_adv = (Neq_adv*(3*Neq_adv + 13))/2

   ! Metal re-solve system size (pp_metals=2 below). ion_system_metals_pp
   ! pins the three H/He rows and solves the metal stages from row 4, so the
   ! system is 3 + 2*n_melem rows -- NOT the global N_eq, which is larger
   ! whenever the equilibrium layout carries extra unknowns (the He 2^3S
   ! metastable, or the molecular H2/H2+/H3+/HeH+ block). Passing N_eq there
   ! leaves those trailing rows of fvec unwritten, i.e. hybrd1 iterating on
   ! uninitialized residuals; the same reason lwa_adv exists above.
   Neq_mpp = 3 + 2*n_melem
   lwa_mpp = (Neq_mpp*(3*Neq_mpp + 13))/2

   ! THE VELOCITY OF EACH ELEMENT'S NUCLEI. The fractions the advection
   ! systems step are fractions of the nuclei of ONE element, and with
   ! He_diffusion those nuclei do not move with the mass-weighted velocity
   ! v: helium drifts against the rest of the gas with the diffusive mass
   ! flux J (gradient, eddy and settling; helium_diffusive_face_flux, the
   ! flux the element row of the run is the divergence of), so
   !
   !    v_He = v + J/(rho Y) ,   v_H = v - J/(rho (1 - Y)) ,
   !
   ! Y the helium mass fraction and component 1 (hydrogen with the metals
   ! and the heavy nuclei of the molecules) carrying -J. With the element
   ! continuity div(n_el v_el) = 0 of a steady state, the stage balance
   ! div(x n_el v_el) = S of ionization_stage_transport.f90 (its
   ! equation (1) without the stage eddy term) is v_el dx/dr = S/n_el,
   ! which is the step the advection systems take with c1 = g h_j/v_el.
   ! Stepping every element on v instead gives each the residence time of
   ! the mixture: MEASURED on the LHS 1140 b He/H 2.09 reference
   ! (Update_EXHALE_stage3 section 82), v_He/v = 0.82-0.97 and
   ! v_H/v = 2.1-1.08 at 1.3-4 R_p, and the He+ fraction, frozen in beyond
   ! 2 R_p, came out 3.5 % low and H+ 6 % high against the transported
   ! stages of the solution.
   ! J is a face flux; the cell value is the mean of the two faces, and
   ! faces 0 and N carry none (the operator's own boundary: the element
   ! crosses either end only with the gas). Ghost cells keep v.
   !
   ! NOT CARRIED: the stage's own eddy term -n_el K dx/dr of (1), which is
   ! second order in r and has no place in a marching step. Its size
   ! against the advective stage flux, K |dln x/dr|/v_el, MEASURED on the
   ! same state (K = 1e9 cm^2 s^-1): He+ 0.035 at 1.3 R_p, 0.005 at
   ! 1.5 R_p, below 0.001 beyond 1.75 R_p; H+ 0.007 and 0.001 there. The
   ! step is therefore valid where v_el >> K |dln x/dr|, which holds over
   ! the corrected rows of that state; a run with a larger K or a slower
   ! wind has to be checked against it.
   v_nuc_H  = v
   v_nuc_He = v
   if (thereis_He .and. he_diffusion) then
      call helium_diffusive_face_flux(rho, T_in, f_sp_in, J_he_face)
      call element_mass_fractions(f_sp_in, Y_elem)
      call mixture_mass_sum(f_sp_in, msum_pp)
      do j = 1, N
         J_he_cell    = 0.5d0*(J_he_face(j-1) + J_he_face(j))
         rho_cgs_cell = rho(j)*n0*mu*msum_pp(j)
         if (Y_elem(j,1) .gt. 0.0d0)                                    &
            v_nuc_He(j) = v(j) + J_he_cell/(rho_cgs_cell*Y_elem(j,1))/v0
         if (Y_elem(j,1) .lt. 1.0d0)                                    &
            v_nuc_H(j)  = v(j) - J_he_cell/                              &
                          (rho_cgs_cell*(1.0d0 - Y_elem(j,1)))/v0
      enddo
   endif

	
	
	! Dimensional total number density profile and temperature
	T_K = T_in*T0
	         
   nhi    = nhi_in*n0
	nhii   = nhii_in*n0
   if (thereis_He) then
		nheii  = nheii_in*n0
		nheiii = nheiii_in*n0
		! Zero-init the triplet so the innermost ghost cell (1-Ng), which
		! the advection loop below never assigns when HeITR is off, does
		! not write uninitialized memory to the output column.
		nheiTR = 0.0
		if (thereis_HeITR) nheiTR = nheiTR_in*n0
		! Split the incoming He I into its two populations. The equilibrium
		! solution reports the summed He I and the metastable, so the singlet
		! is formed here once and then carried as a population of its own.
		! The subtraction is well conditioned wherever it is made: the
		! equilibrium metastable is bounded by its own balance, whose loss
		! side is dominated by the 2^3S -> 1^1S decay A31, and stays orders of
		! magnitude below the summed He I. It is taken through
		! he_ground_singlet_density all the same, so that a restart whose two
		! helium columns disagree cannot start the advection from a negative
		! singlet, and so that the run reports it if one did.
		nheiS  = nhei_in*n0
		if (thereis_HeITR)                                              &
			nheiS = he_ground_singlet_density(nhei_in, nheiTR_in)*n0
		nhei   = nheiS + nheiTR
	else
		! Helium off. The helium-free branch of the ionization loop below
		! zeroes these at its END, but the first pass reads them before that
		! -- in nhe, in the electron sum and in eval_cool -- so they are
		! defined here.
		nheiS  = 0.0
		nhei   = 0.0
		nheii  = 0.0
		nheiii = 0.0
		nheiTR = 0.0
	endif


	! Select the metal treatment for the post-process (set via 'pp_metals'
	! in metals.inp -> pp_metal_mode). Metals are off => zero either way.
	!   0 metal-free (legacy), 1 frozen eq metals, 2 re-solve.
	! pp_metal_on (module flag in equation_T) tells the cell-by-cell temperature
	! solve to include the metal cooling/brem/n_e terms.
	if (thereis_metals .and. pp_metal_mode >= 1) then
		nm_w        = nm_in*n0          ! eq metal densities [cgs] (frozen or seed)
		pp_metal_on = .true.
	else
		nm_w        = 0.0d0
		pp_metal_on = .false.
	endif

	! Re-solve mode (pp_metals=2): the metal ionization balance is recomputed
	! per cell inside the post-process loop at the advection-corrected H/He and
	! the post-process temperature (see the block after the H/He advection
	! solve). The element totals are conserved by ionization (only the stage
	! split changes), so freeze each element's total from the eq densities and
	! reuse it as a constant throughout. nm_w starts at the eq split and is
	! overwritten with the re-solved split each iteration.
	if (pp_metal_mode == 2 .and. thereis_metals .and. thereis_He) then
		do im = 1,n_melem
			i0 = melem_i0(im)
			if (melem_top(im) >= 2) then
				nm_tot_pp(:,im) = nm_w(:,i0) + nm_w(:,i0+1) + nm_w(:,i0+2)
			else
				nm_tot_pp(:,im) = nm_w(:,i0) + nm_w(:,i0+1)
			endif
		enddo
	endif


	! Non-converged cell solves, accumulated over every pass (see declaration)
	n_adv_noconv   = 0
	n_metal_noconv = 0
	n_T_noconv     = 0
	n_T_res_root   = 0
	! The column energy solve: nothing solved yet, no pass rejected.
	col_outcome  = adv_col_not_run
	col_reason   = adv_col_reason_none
	col_cell     = 0
	col_unknowns = 0
	n_col_passes_rejected = 0
	col_worst_cell  = 0
	col_worst_R     = 0.0d0
	col_worst_rel   = 0.0d0
	col_worst_tol_rel = 0.0d0

	! The weights of the step into each cell. The two lowest rows have no
	! two cells below them and are never the target of a BDF2 step.
	bdf_a1 = 1.0d0
	bdf_a2 = 0.0d0
	bdf_g  = 1.0d0
	bdf2_admissible = .false.
	do j = 3-Ng,N+Ng
		call variable_step_bdf2_weights(r(j) - r(j-1), r(j-1) - r(j-2),     &
		                                bdf_a1(j), bdf_a2(j), bdf_g(j))
		bdf2_admissible(j) = (r(j) - r(j-1) <=                              &
		                      bdf2_step_ratio_limit*(r(j-1) - r(j-2)))
	enddo

	! Iterate the post processing: a fixed number of passes, with no test of
	! convergence between them, so the outer iteration of the product is
	! never verified (adv_outer_unverified).
	do k = 1,n_pp_passes

	n_bdf2_comp = 0;  n_be_comp = 0;  n_bdf2_comp_retaken = 0
	n_bdf2_T    = 0;  n_be_T    = 0;  n_bdf2_T_retaken    = 0
	
	! Summed He I from the two populations carried through the pass, then the
	! H and He NUCLEI totals from the one shared definition (utils).  The
	! molecular and oxygen carriers are omitted because this pass does not
	! carry them (see the module-header composition note); the He 2^3S
	! nuclei enter through nhei, which already sums the two populations.
	nhei = nheiS + nheiTR
	call hydrogen_helium_nuclei_density(nhi,nhii,nhei,nheii,nheiii,nh,nhe)

	! Free electron density (assuming overall neutrality; nm_w adds the
	! metal electrons under the eos_metals policy). Molecular-ion electrons are
	! excluded here -- the post-process does not carry the molecular densities
	! (see the module-header composition note).
	call calc_ne(nhii,nheii,nheiii,ne,nm_w)

	! Ionized fraction handed to the photoelectron partition: the TOTAL free
	! electron density over the H and He nuclei, the quantity Dalgarno, Yan &
	! Liu (1999) section 7 define ("the number density ratio of the electrons
	! to the hydrogen and helium nuclei").  This one is of the composition
	! the pass STARTS from, and serves the photoionization rates that drive
	! the composition sweep; the sweep and the metal re-solve change the
	! densities, so the heating of the composition they return is formed
	! with the ratio of that composition, recomputed after them (THE COUNTS
	! OF THE COMPOSITION, below).  On this path the
	! molecular-ion electrons are missing from ne for the reason just given,
	! so a molecular run's post-process carries a slightly low ratio; the
	! equilibrium pass, which does have them, does not.
	xion = min(max(ne/max(nh + nhe, 1.0d-99), 0.0d0), 1.0d0)

	! Metal electrons alone, for the electron density inside the advection
	! residuals. Same definition as calc_ne above and as the equilibrium
	! residual's metal_electron_sum (stage 1 counts one electron, stage 2 two).
	! Unlike calc_ne this is NOT conditioned on eos_include_metals: that switch
	! governs whether the metals enter the gas mass/particle budget, while the
	! recombination terms of the residuals need the true free electron density.
	! nm_w is already zero when metals are off or the post-process runs
	! metal-free, so this is zero there.
	ne_metal = 0.0d0
	do im = 1,n_mion
		if (mion_stage(im) .gt. 0)                                        &
			ne_metal = ne_metal + dble(mion_stage(im))*nm_w(:,im)
	enddo

	! Cell-by-cell opacity pressure factor ('P' model; =1 otherwise)
	do j = 1-Ng,N+Ng
		opa_pf(j) = opacity_pT_factor((nh(j)+nhe(j)+ne(j))*kb_erg*T_K(j))
	enddo

   ! Calculate the photoionization rates

	! nh2_solved: the solved state's H2 decides the molecular gate of the
	! low-energy photoelectron partition (low_energy_electron_partition),
	! as it does in the equilibrium solve; this gas carries no H2 itself.
	if (thereis_He) then
		call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm_w, xion,     &
					 P_HI,P_HeI,P_HeII,P_HeITR, P_m,        &
					 dum_v1,dum_v2, P_m2=P_m2,              &
					 nh2_solved=nmol_eq(:,1))
  	else
	  	call PH_heat_H(nhi, xion, P_HI,dum_v1,dum_v2)
  	endif

   !---- Recombination rates ----!

	! nmol is not passed: the _adv reconstruction is molecule-free, so
	! eval_cool builds the atomic electron sum and leaves out the H3+ infrared
	! cooling (module-header composition note). Both omissions are the same
	! approximation and both end at the molecular layer.
	call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm_w,            &
	  			   rchiiB,rcheiiB,rcheiiiB, rec_m_pp,             &
				   a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m_pp,          &
				   tcool_in, a_ion_HeITR=a_ion_HeITR)

 	
 	if (thereis_HeITR) then
		! rcheiiB becomes the case-B capture into the singlets; the escaping
		! ground captures are added below.
		call HeITR_coeffs(T_K,rcheiTR,rcheiiB,A31,q13,q31g,q31a,q31b,Q31)

	endif

	! Recombination radiation absorbed on the spot (recombination_radiation_
	! absorbed). Replace the recombination coefficients by the net ones and
	! add the photoionization rates driving the advection ODE (the heating
	! of the same photons is part of the one heating assembly below,
	! heating_of_composition, which calls the same routine; dheat_hrc is not
	! read). Mirrors ionization_equilibrium; molecule-free (nh2_pp = 0).
	if (use_h_rec_escape .or. (use_he_rec_coupling .and. thereis_He)) then
		nh2_pp = 0.0d0
		call recombination_radiation_absorbed(T_K, nhi, nhii, nh2_pp,     &
		                     nhei, nheii, nheiii, nheiTR, ne, nm_w,        &
		                     A31, q31a, q31b,                              &
		                     rchiiB_hrc, rcheiiB_hrc, rcheiiiB_hrc,        &
		                     dP_HI_hrc, dP_HeI_hrc, dP_H2_hrc, dP_m_hrc,   &
		                     dheat_hrc, dP_m2=dP_m2_hrc, xion=xion,        &
		                     nh2_solved=nmol_eq(:,1))
		rchiiB = rchiiB_hrc
		if (use_he_rec_coupling .and. thereis_He) then
			rcheiiB  = rcheiiB_hrc
			rcheiiiB = rcheiiiB_hrc
		endif
		P_HI    = P_HI  + dP_HI_hrc
		P_HeI   = P_HeI + dP_HeI_hrc
		! The metal share of the same photons feeds the metal re-solve
		! (pp_metals=2) exactly as it feeds the equilibrium solve.
		if (thereis_metals) then
			P_m  = P_m  + dP_m_hrc
			P_m2 = P_m2 + dP_m2_hrc
		endif
	endif


	! The steady advection correction is conditional on the recorded
	! flow's mass row: |R_1|/s_1 <= adv_conditional_tol. The measure is the
	! fractional change of numerical face mass flux across a cell, obtained
	! from assemble_residual and residual_row_scale, as in stationary
	! certification. The file records it as adv_mass_row for every row;
	! whole-state certification with cert_tol_mass is separate.
	!
	! A larger mass residual makes omitted divergence terms comparable to
	! the retained steady terms and refuses both thermal and ionization
	! corrections. The cell retains the run's temperature and equilibrium
	! composition. In backup/regression/hydrostatic_column, about 8% mass
	! flux change per cell drove the full equation to 1.7e7 K from 1124 K;
	! omitting enthalpy divergence yielded a 0.78 K adiabat. Outflow sign
	! and local Damkohler tests alone do not detect this error.
	!
	! A small mass row is necessary but does not establish momentum or
	! energy balance. Evaluate the refusal once on the input state at the
	! first pass, using its equilibrium composition, so the refused set
	! does not depend on iteration history (write_adv_validity_header).
	if (k == 1) then
		call calc_mmw(nh,nhe,ne,mmw_in,nm_w)

		! The state in conservative variables. Whether it is the very state
		! whose row terms the run last stored is asked BEFORE assembling, so
		! the answer is about the run and not about this call, and it is
		! reported: the primitive-to-conservative round trip through the
		! caloric EOS is not the identity to the last bit.
		W_state(1,:) = rho
		W_state(2,:) = v
		W_state(3,:) = p
		call W_to_U(W_state, u_state)
		terms_are_run_state = row_terms_describe_state(u_state)

		! The particle count the residual is given. Only the MASS row is read
		! below, and that row is the flux difference of rho v alone: it does
		! not contain the particle count, the heating or the cooling, which
		! enter the momentum and energy rows and the operator-split transport
		! sources. The count handed over is the one of the reconstruction this
		! post-process solves with, so nothing here reads a composition the
		! rest of the routine does not.
		call calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_part_state,nm_w)
		n_part_state = n_part_state + ne

		! The face-limiter ledger counts the UPDATES OF THE RUN, and the run
		! reports it after this routine returns; assembling a residual is a
		! measurement of a state and not a step, so the count is put back.
		n_limited_before = n_faces_positivity_limited
		call assemble_residual(u_state, n_part_state, f_sp_in, heat, cool,  &
		                       R_state)
		n_faces_positivity_limited = n_limited_before

		mass_row_within_tol = .true.
		adv_mass_row        = 0.0d0
		n_not_stationary    = 0
		mass_row_max        = 0.0d0
		r_mass_row_max      = 0.0d0
		do j = 2-Ng,N+Ng
			! The cell's own mass row over the largest term that row holds.
			! R_state(1,1-Ng) is not read: the right-hand side is formed over
			! 2-Ng..N+Ng, which is also the range of every loop below, and
			! adv_mass_row stays zero in that row, which the status fields
			! report as not_evaluated.
			mass_row = abs(R_state(1,j))                                    &
			           /max(residual_row_scale(1,j,u_state), cert_scale_floor)
			adv_mass_row(j) = mass_row
			if (mass_row > mass_row_max) then
				mass_row_max   = mass_row
				r_mass_row_max = r(j)
			endif
			! Written so that a row that is not a number is refused rather
			! than accepted by a comparison that is false either way.
			if (.not. (mass_row <= adv_conditional_tol)) then
				mass_row_within_tol(j) = .false.
				n_not_stationary        = n_not_stationary + 1
			endif
		enddo

		! THE DIVERGENCE OF THE MASS FLUX the energy equation holds, cell by
		! cell: the mass row just assembled, which is the face mass fluxes of
		! the Riemann solve of this state differenced over the cell's own
		! control volume (the mass row has no source, so R_1 is that
		! divergence). See the statement before the energy loop for why this
		! and not a difference of the cell-centred rho v r^2.
		div_rhov_state        = 0.0d0
		div_rhov_state(2-Ng:) = R_state(1,2-Ng:)

		! THE SENSITIVITY OF THE ANSWER, next to the decision above. How much
		! of the energy equation the non-stationarity of the recorded flow
		! carries is the size of its enthalpy flux term against the two terms
		! the balance keeps (enthalpy_flux_term_ratio, this module), formed
		! with the step the energy loop takes into the cell (bdf2_step_ratio_
		! limit) on the run's own temperature. It is reported and refuses
		! nothing: the same ratio is reached with a mass divergence of any
		! size once the kept terms are large.
		enthalpy_flux_ratio_max  = 0.0d0
		r_enthalpy_flux_max      = 0.0d0
		n_enthalpy_term_dominant = 0
		do j = 2-Ng,N+Ng
			! The three terms, with the common 1/dr_step divided out, in the
			! variables of the residual: e = E(x_H2,T)/mu with E the energy
			! per particle of the caloric EOS, p/rho = T/mu, h = e + p/rho.
			! Each point of the history is evaluated at the composition of
			! ITS OWN cell, as the residual does: the rovibrational energy of
			! H2 travels with the gas that holds the molecules.
			use_bdf2 = bdf2_admissible(j)
			dr_step  = step_width(j, use_bdf2)
			if (use_bdf2) then
				e_hist   = bdf_a1(j)*specific_internal_energy(j-1, T_in(j-1), mmw_in(j-1)) &
				         - bdf_a2(j)*specific_internal_energy(j-2, T_in(j-2), mmw_in(j-2))
				rho_hist = bdf_a1(j)*rho(j-1) - bdf_a2(j)*rho(j-2)
			else
				e_hist   = specific_internal_energy(j-1, T_in(j-1), mmw_in(j-1))
				rho_hist = rho(j-1)
			endif
			e_up = specific_internal_energy(j, T_in(j), mmw_in(j))
			w_up = T_in(j)/mmw_in(j)
			h_up = e_up + w_up
			enthalpy_flux_ratio = enthalpy_flux_term_ratio(h_up,        &
			     div_rhov_state(j), dr_step, rho(j), v(j), e_up, e_hist,   &
			     w_up, rho_hist)
			if (enthalpy_flux_ratio > enthalpy_flux_ratio_max) then
				enthalpy_flux_ratio_max = enthalpy_flux_ratio
				r_enthalpy_flux_max     = r(j)
			endif
			if (enthalpy_flux_ratio > enthalpy_ratio_report_level)          &
				n_enthalpy_term_dominant = n_enthalpy_term_dominant + 1
		enddo

		! THE CELL CLASS THE RECONSTRUCTION MODELS. This post-process solves
		! an H/He + trace-metal gas: calc_ne and calc_ntot are called without
		! the molecular and oxygen carriers, so in a cell where those omitted
		! species carry more of the particle count than the ones it does
		! carry, the counts every equation here is solved with are not the
		! counts of that gas, and neither the temperature nor the composition
		! of the row is a statement about it. The comparison is between two
		! counts of the same cell and carries no threshold; a run with no
		! molecules and no oxygen chemistry has n_omitted = 0 in every cell.
		cell_class_modelled = .true.
		if (thereis_mol .or. thereis_oxychem) then
			do j = 1-Ng,N+Ng
				n_omitted = 0.0d0
				if (thereis_mol)                                             &
					n_omitted = n_omitted + nmol_eq(j,1) + nmol_eq(j,2)       &
					          + nmol_eq(j,3) + nmol_eq(j,4)
				if (thereis_oxychem)                                         &
					n_omitted = n_omitted + nox_eq(j,1) + nox_eq(j,2)         &
					          + nox_eq(j,3)
				n_carried = n_part_state(j)
				cell_class_modelled(j) = (n_omitted <= n_carried)
			enddo
		endif

		! THE COMPOSITION CLASS: the same statement about the BALANCE of each
		! population the advection systems solve. They are molecule-free, so
		! the destruction of H+, He+ and He(2^3S) by H2 (R10 and R13; R17 and
		! R23; the He(2^3S) + H2 ionization) is absent from them, and in a
		! cell where that omitted loss of one population exceeds the loss the
		! systems carry for it, the corrected population answers another
		! balance than the gas has. MEASURED 2026-10-01 on the molecular
		! photochemical conduction state: at 1.005 to 1.08 R_p, where x(H2)
		! falls from 0.38 to 4e-3, the corrected n(He 2^3S) was 70 to 5600
		! times the run's, whose balance holds the He(2^3S) + H2 channel;
		! those rows passed the particle-count class above. The comparison is
		! between two rates of the same population of the same cell, on the
		! state handed in, and carries no threshold. A row this class
		! refuses carries the run's own composition and temperature
		! (adv_unsupported in both fields).
		n_species_class_refused = 0
		if (thereis_mol) then
			do j = 1-Ng,N+Ng
				n_h2_cell = nmol_eq(j,1)
				if (.not. (n_h2_cell > 0.0d0)) cycle
				if (.not. cell_class_modelled(j)) cycle
				! H+: R10 (H+ + H2(v>=4)) and R13 (H+ + H2 + M) against
				! recombination.
				loss_omitted = (rk_R10_Hp_H2v4(T_K(j))                      &
				             + rk_R13_Hp_H2_M(n_part_state(j)))*n_h2_cell
				loss_carried = rchiiB(j)*ne(j)
				if (.not. (loss_omitted <= loss_carried))                  &
					cell_class_modelled(j) = .false.
				if (thereis_He) then
					! He+: R17 (dissociative) and R23 (charge transfer) with
					! H2 against recombination and ionization to He++.
					loss_omitted = (rk_R17_Hep_H2_diss(T_K(j))                &
					             + rk_R23_H2_Hep_cx())*n_h2_cell
					rec_HeII_tot = rcheiiB(j)
					if (thereis_HeITR) rec_HeII_tot = rec_HeII_tot + rcheiTR(j)
					loss_carried = P_HeII(j)                                  &
					             + (rec_HeII_tot + a_ion_HeII(j))*ne(j)
					if (.not. (loss_omitted <= loss_carried))              &
						cell_class_modelled(j) = .false.
				endif
				if (thereis_HeITR) then
					! He(2^3S): its ionization by H2 against every loss the
					! He 2^3S row of the advection system holds (the
					! He(2^3S) relaxation rate of the validity block).
					loss_omitted = ioniz_HeI23S_H2(T_K(j))*n_h2_cell
					loss_carried = A31 + P_HeITR(j)                           &
					     + (q31g(j) + q31a(j) + q31b(j) + a_ion_HeITR(j))*ne(j) &
					     + Q31(j)*nhi(j)
					if (.not. (loss_omitted <= loss_carried))              &
						cell_class_modelled(j) = .false.
				endif
				if (.not. cell_class_modelled(j))                          &
					n_species_class_refused = n_species_class_refused + 1
			enddo
		endif

		! THE ENERGY CLASS THE RECONSTRUCTION MODELS. The energy equation
		! solved below balances the heating and cooling this post-process
		! assembles for an H/He + trace-metal gas; the run's own energy
		! equation balanced every channel it solved, the molecular and oxygen
		! ones included (the heat and cool handed in). On the state handed
		! in, at the run's temperature and composition, the difference of the
		! two is exactly what the post-process omits, and where it is larger
		! than what the post-process carries, the temperature the equation
		! returns is set by the absence of the channels that set the gas's
		! temperature, and the row is not a statement about that gas.
		! MEASURED on the LHS 1140 b molecular states (He/H = 1.6, 500
		! cells): below 1.08 R_p the molecular chemistry heat is 1.4e-7 of
		! the 1.5e-7 erg cm^-3 s^-1 deposited at 1.001 R_p and the H3+
		! infrared band carries the cooling, while the carried channels hold
		! 8.4e-9 heating; the energy equation there then returned 1280 /
		! 1051 / 600 K at 1.002 R_p on 500 / 1000 / 2000 cells against a
		! run temperature of 1500 K. Like the particle count above, the
		! comparison is between two numbers of the same cell and carries no
		! threshold; a run in which every channel is carried has
		! q_omitted = 0 up to the arithmetic of the two assemblies.
		call calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm_w)
		call heating_of_current_composition(theat_of_run_composition)
		energy_class_modelled  = .true.
		n_energy_class_refused = 0
		r_energy_class_top     = 0.0d0
		do j = 1-Ng,N+Ng
			q_omitted = abs(heat(j) - theat_of_run_composition(j))           &
			          + abs(cool(j) - tcool_in(j)/q0)
			q_carried = abs(theat_of_run_composition(j))                     &
			          + abs(tcool_in(j)/q0)
			! Written so that a comparison that is not a number refuses.
			energy_class_modelled(j) = (q_omitted <= q_carried)
			if (.not. energy_class_modelled(j)) then
				n_energy_class_refused = n_energy_class_refused + 1
				r_energy_class_top     = max(r_energy_class_top, r(j))
			endif
		enddo
	endif


	!---- Where is the advection correction valid? ----!
	!
	! The correction replaces the local ionization balance of a cell by the
	! steady advection-ionization ODE, integrated along the flow by the step
	! of the recursion (bdf2_step_ratio_limit). Three
	! conditions make that replacement carry no information; where any of them
	! holds the cell keeps the converged equilibrium ionization instead.
	!
	!  (i)   Inflow, v <= 0 in the cell or the one below it, or the same for
	!        the nuclei of hydrogen or helium (v_nuc_H, v_nuc_He) -- PHYSICAL.
	!        The step takes the upstream state from the cells below, which are
	!        not upstream when the gas (or one element) moves inward (the
	!        breathing base, or helium settling faster than the wind lifts
	!        it). The residence time h/v is then negative as well.
	!
	!  (ii)  Da = (dr/v_el)*nu_relax > Da_local_equilibrium -- PHYSICAL. The
	!        Damkohler number compares the time the gas spends in the cell with
	!        the relaxation time of the level populations. Da >> 1 means the
	!        populations approach local equilibrium rapidly. Neglecting
	!        advection is an approximation when that equilibrium varies slowly
	!        over the relaxation length; it is not an exact ODE solution.
	!
	!        nu_relax is the minimum of these frozen population-rate estimates,
	!        without diagonalizing the coupled reaction Jacobian:
	!           H I/H II      : P_HI + (a_ion_HI + alpha_HII)*n_e
	!           He I/He II    : P_HeI
	!                           + (a_ion_HeI + alpha_HeII + alpha_HeI23S)*n_e
	!           He II/He III  : P_HeII + (a_ion_HeII + alpha_HeIII)*n_e
	!           He(2^3S)      : A31 + P_HeITR
	!                           + (q31g + q31a + q31b + a_ion_HeITR)*n_e
	!                           + Q31*n_HI
	!        (the He(2^3S) row is exactly the loss side of fvec(4) of
	!        adv_implicit_HeH_TR, with the same rate coefficients from
	!        HeITR_coeffs / eval_cool -- no rate is redefined here.)
	!
	!        Each element's rates are taken with the residence time of that
	!        element's nuclei, h/v_H for the hydrogen row and h/v_He for the
	!        helium rows, the times their steps use.
	!
	!        Taking the minimum is what makes the gate a statement about the
	!        cell rather than about one species: the systems solve the whole
	!        H/He vector at once, so a cell may be pinned to equilibrium only
	!        if EVERY solved population is equilibrated. He(2^3S) relaxes
	!        orders of magnitude more slowly than H (A31 = 1.27e-4 s^-1 sets
	!        the floor), so gating the vector on the H rate alone froze the
	!        metastable at its equilibrium value in cells where it is in fact
	!        advected -- an order-of-magnitude step in the _adv 2^3S profile
	!        wherever a sharp H ionization front crossed the threshold.
	!        n_e here is the metal-inclusive electron density, the same one the
	!        equilibrium solve used.
	!
	!  (iii) x_HII,eq < xHII_adv_min -- NUMERICAL. The residuals carry the
	!        neutral fraction x_HI, with ions extracted as (1-x_HI)*n_h.
	!        MINPACK uses xtol=sqrt(eps)=1.5e-8 and a relative finite-difference
	!        perturbation of that size; neither is an absolute root-error bound.
	!        In the shielded-base case x_HII,eq~1e-13, extracted fractions were
	!        observed at about 1e-8 with either sign, contaminating the upwind
	!        history. The 1e-6 cutoff avoids that measured trace-ion failure.
	!
	! (i) and (ii) are statements about the flow and (iii) about the
	! representation of the unknown, so none of them depends on whether metal
	! cooling is switched on.

	! (iv) The mass row of the cell stands above adv_conditional_tol, so the
	!      terms the steady equations drop are larger than that fraction of
	!      the ones they keep and neither the ionization nor the energy
	!      correction is accurate to the stated fraction there -- PHYSICAL,
	!      stated and measured in the block above, which also records why it
	!      is not covered by (i) or (ii).
	!
	! The two status fields are built here and in the energy loop, each where
	! its own decision is taken. The composition field is complete when the
	! ionization loop below has run; the temperature field is completed in
	! the energy loop, which has conditions of its own. Both start at
	! not_evaluated, so a row no loop reaches says that rather than claiming
	! a correction: the ionization loop runs from 2-Ng and the energy loop
	! from 3-Ng, and the rows below those are ghosts with no upstream state
	! for the step to read.

	adv_correction_valid = .true.
	adv_correction_valid(1-Ng) = .false.   ! inner boundary: never corrected
	adv_T_status    = adv_not_evaluated
	adv_comp_status = adv_not_evaluated
	n_adv_eq = 0
	n_stationarity_only = 0
	do j = 2-Ng,N+Ng
		if (v(j) <= 0.0d0 .or. v(j-1) <= 0.0d0 .or.                       &
		    v_nuc_H(j) <= 0.0d0 .or. v_nuc_H(j-1) <= 0.0d0 .or.           &
		    v_nuc_He(j) <= 0.0d0 .or. v_nuc_He(j-1) <= 0.0d0) then
			! The gas, or the nuclei of one element, enter the cell from
			! above: the outward step has no upstream state.
			adv_correction_valid(j) = .false.
		else
			! The residence time of the nuclei of each element in the cell,
			! h_j/v_el: the velocity of the cell the rates are evaluated in,
			! as in the step (bdf2_step_ratio_limit).
			t_cross    = (r(j) - r(j-1))*R0/(v_nuc_H(j)*v0)
			t_cross_He = (r(j) - r(j-1))*R0/(v_nuc_He(j)*v0)
			nu_relax = P_HI(j) + (a_ion_HI(j) + rchiiB(j))*ne(j)
			Da_slowest = t_cross*nu_relax
			if (thereis_He) then
				! He II -> He I recombination: with the triplet on, rcheiiB is
				! the singlet channel alone (HeITR_coeffs overwrites it) and
				! rcheiTR is the triplet one, exactly as the He I row of the
				! residuals adds them.
				rec_HeII_tot = rcheiiB(j)
				if (thereis_HeITR) rec_HeII_tot = rec_HeII_tot + rcheiTR(j)
				nu_relax_He = min(                                          &
				     P_HeI(j)  + (a_ion_HeI(j)  + rec_HeII_tot)*ne(j),      &
				     P_HeII(j) + (a_ion_HeII(j) + rcheiiiB(j) )*ne(j))
				if (thereis_HeITR)                                          &
					nu_relax_He = min(nu_relax_He, A31 + P_HeITR(j)          &
					     + (q31g(j) + q31a(j) + q31b(j) + a_ion_HeITR(j))*ne(j) &
					     + Q31(j)*nhi(j))
				Da_slowest = min(Da_slowest, t_cross_He*nu_relax_He)
			endif
			xHII_eq = nhii_in(j)/max(nhi_in(j) + nhii_in(j), 1.0d-300)
			if (Da_slowest > Da_local_equilibrium .or.                     &
			    xHII_eq < xHII_adv_min)                                    &
				adv_correction_valid(j) = .false.
		endif
		! Condition (iv) last, so that the count says how many cells it adds
		! beyond (i) to (iii) rather than how many it holds on.
		if (.not. mass_row_within_tol(j)) then
			if (adv_correction_valid(j)) n_stationarity_only = n_stationarity_only + 1
			adv_correction_valid(j) = .false.
		endif
		! The cell class the reconstruction does not model (particle count
		! or the balance of a solved population, formed with the stationarity
		! condition): the composition of the row is the run's own, and the
		! upwind cascade takes the run's composition as its history there.
		if (.not. cell_class_modelled(j)) adv_correction_valid(j) = .false.
		if (.not. adv_correction_valid(j)) n_adv_eq = n_adv_eq + 1
		! The composition of a cell any of the four conditions refuses is the
		! equilibrium one the run converged to: the run's own value, kept.
		! A cell that is corrected here can still fail its own solve, which
		! the loop below records. The temperature field gets what conditions
		! (i) and (iv) say about it -- the energy step has no
		! upstream state under inflow either, and a non-stationary cell has
		! no steady energy equation -- and the energy loop adds the rest.
		if (adv_correction_valid(j)) then
			adv_comp_status(j) = adv_corrected
		else if (.not. cell_class_modelled(j)) then
			adv_comp_status(j) = adv_unsupported
		else
			adv_comp_status(j) = adv_retained
		endif
	enddo


   ! Evolve the species with the advection term kept, in the ODE form
   ! v dx/dr = R(x), cell by cell upward by the step of the recursion
   ! (bdf2_step_ratio_limit): BDF2 where the cell below was itself
   ! integrated in this pass, backward Euler at the start of the recursion
   ! and after every cell it does not integrate. Point values are used for
   ! cell averages; they agree to O(dr^2). The recorded velocity profile is
   ! held fixed.

   if (.not.thereis_He) then

		do j = 2-Ng,N+Ng
			! Outside its validity range the advection correction carries no
			! information; the cell keeps the converged equilibrium ionization.
			if (.not. adv_correction_valid(j)) then
				nhi(j)  = nhi_in(j)*n0
				nhii(j) = nhii_in(j)*n0
				cycle
			endif

			! The rates of the cell
			adv_cell%nh = nh(j)
			adv_cell%P_HI = P_HI(j)
			adv_cell%rchiiB = rchiiB(j)
			adv_cell%a_ion_HI = a_ion_HI(j)
			! Metal electrons of this cell, per H nucleus, for the electron
			! density the recombination terms of the residual see.
			adv_cell%xe_metal = ne_metal(j)/max(nh(j),1.0d-30)

			use_bdf2 = bdf2_admissible(j) .and.                              &
			           adv_comp_status(j-1) == adv_corrected
			call composition_step(j, use_bdf2)
			if (use_bdf2 .and. info == 1 .and.                              &
			    .not. populations_are_admissible()) then
				n_bdf2_comp_retaken = n_bdf2_comp_retaken + 1
				use_bdf2 = .false.
				call composition_step(j, use_bdf2)
			endif
			if (use_bdf2) then
				n_bdf2_comp = n_bdf2_comp + 1
			else
				n_be_comp = n_be_comp + 1
			endif

			! Non-converged cell: keep the equilibrium ionization (see the
			! H/He branch below for why the returned iterate is discarded).
			if (info /= 1) then
				nhi(j)  = nhi_in(j)*n0
				nhii(j) = nhii_in(j)*n0
				n_adv_noconv = n_adv_noconv + 1
				adv_comp_status(j) = adv_failed
				cycle
			endif

			! Extract solution profiles
			nhi(j)    = sys_x(1)*nh(j)
			nhii(j)   = (1.0 - sys_x(1))*nh(j)

		enddo

		   nheiS  = 0.0
		   nhei   = 0.0
		   nheii  = 0.0
		   nheiii = 0.0
		   nheiTR = 0.0

	else

		do j = 2-Ng,N+Ng
			! Outside its validity range the advection correction carries no
			! information; the cell keeps the converged equilibrium H/He
			! ionization (see the no-He branch).
			if (.not. adv_correction_valid(j)) then
				call pin_cell_to_equilibrium(j)
				cycle
			endif

			! The rates of the cell
			adv_cell%nh  = nh(j)
			adv_cell%P_HI  = P_HI(j)
			adv_cell%P_HeI  = P_HeI(j)
			adv_cell%P_HeII  = P_HeII(j)
			adv_cell%rchiiB  = rchiiB(j)
			adv_cell%rcheiiB = rcheiiB(j)
			adv_cell%rcheiiiB = rcheiiiB(j)
			adv_cell%a_ion_HI = a_ion_HI(j)
			adv_cell%a_ion_HeI = a_ion_HeI(j)
			adv_cell%a_ion_HeII = a_ion_HeII(j)
			! He <-> H charge-exchange rate coefficients (Huang Table 4 group
			! B, and He2+ + H0), read by he_h_cx_fvec_adv in the H/He adv
			! systems. T-only, so evaluate once per cell; the adv residual
			! adds nothing when he_h_charge_exchange is off.
			call he_h_cx_rates(T_K(j), adv_cell%kcx_He0_Hp,                &
			                           adv_cell%kcx_Hep_H0,                &
			                           adv_cell%kcx_Hepp_H0)
			! Effective He/H for the electron density inside the adv system:
			! the global HeH normally, the local (diffused) nhe/nh when
			! He_diffusion is on.
			if (he_diffusion) then
				adv_cell%heh_loc = nhe(j)/max(nh(j),1.0d-30)
			else
				adv_cell%heh_loc = HeH
			endif
			! Metal electrons of this cell, per H nucleus, for the electron
			! density the recombination terms of the residual see.
			adv_cell%xe_metal = ne_metal(j)/max(nh(j),1.0d-30)

			if (thereis_HeITR) then
				adv_cell%rcheiTR = rcheiTR(j)
				adv_cell%A31 = A31
				adv_cell%P_HeITR = P_HeITR(j)
				adv_cell%q13 = q13(j)
				adv_cell%q31g = q31g(j)
				adv_cell%q31a = q31a(j)
				adv_cell%q31b = q31b(j)
				adv_cell%Q31 = Q31(j)
				adv_cell%a_ion_HeITR = a_ion_HeITR(j)
			endif

			use_bdf2 = bdf2_admissible(j) .and.                              &
			           adv_comp_status(j-1) == adv_corrected
			call composition_step(j, use_bdf2)
			if (use_bdf2 .and. info == 1 .and.                              &
			    .not. populations_are_admissible()) then
				n_bdf2_comp_retaken = n_bdf2_comp_retaken + 1
				use_bdf2 = .false.
				call composition_step(j, use_bdf2)
			endif
			if (use_bdf2) then
				n_bdf2_comp = n_bdf2_comp + 1
			else
				n_be_comp = n_be_comp + 1
			endif

			! A cell whose advection system did not converge carries no
			! correction: the returned iterate satisfies neither the
			! advection balance it was asked to solve nor the equilibrium
			! balance it started from, so it is not a state of the gas. The
			! equilibrium solution of that same cell is, and it is what the
			! validity conditions above already fall back to.
			if (info /= 1) then
				call pin_cell_to_equilibrium(j)
				n_adv_noconv = n_adv_noconv + 1
				adv_comp_status(j) = adv_failed
				cycle
			endif

			! Extract solution profiles
			nhi(j)    = sys_x(1)*nh(j)
			nhii(j)   = (1.0 - sys_x(1))*nh(j)
			nheiS(j)  = sys_x(2)*nhe(j)
			nheiii(j) = sys_x(3)*nhe(j)
			if (thereis_HeITR) then
				nheiTR(j) = sys_x(4)*nhe(j)
				nheii(j)  = (1.0 - sys_x(2) - sys_x(3) - sys_x(4))*nhe(j)
			else
				nheiTR(j) = 0.0
				nheii(j)  = (1.0 - sys_x(2) - sys_x(3))*nhe(j)
			endif
			nhei(j)   = nheiS(j) + nheiTR(j)

		enddo

	endif ! End if thereis_He
	      
      
      
   !------- Fix stationarity of new pressure profile -------!
      
   !---- Update densities and temperature ----!
	
	! Number densities of the H/He the sweep returned, and the free electron
	! density of that H/He with the metal split the pass started from: the
	! electron density the charge-exchange coefficients of the metal
	! re-solve below are set at. The counts every later step reads are
	! formed after that re-solve (THE COUNTS OF THE COMPOSITION, below).
   nhei = nheiS + nheiTR
	call hydrogen_helium_nuclei_density(nhi,nhii,nhei,nheii,nheiii,nh,nhe)
   call calc_ne(nhii,nheii,nheiii,ne,nm_w)


	!---- Re-solve metal ionization (pp_metals=2) ----!
	! Recompute the metal ionization balance per cell at the advection-
	! corrected H/He and the post-process temperature T_K. The full coupled
	! H/He+metal residual (ion_system_HeH_metals) is reused, but rows 1-3 are
	! pinned by the ion_system_metals_pp wrapper so only the metal stages move;
	! the electron density inside the residual then combines the fixed H/He
	! with the re-solved metals (matching the equilibrium solver's n_e). The
	! photo/collisional/recombination rates (P_m, aion_m_pp, rec_m_pp) and the
	! charge-exchange couplings to H+/He+ come from the current pass, so the
	! k-loop drives metals, H/He and T to a joint fixed point. Element totals
	! are conserved (nm_tot_pp), so only the stage split is updated.
	if (pp_metal_mode == 2 .and. thereis_metals .and. thereis_He) then
		do j = 1-Ng,N+Ng
			if (nh(j) <= 0.0d0) cycle
			! A row whose composition is refused carries the run's own
			! composition, the metal stages included.
			if (adv_comp_status(j) .ne. adv_corrected) then
				nm_w(j,:) = nm_in(j,:)*n0
				cycle
			endif

			! Charge-exchange rate coefficients for this cell's temperature
			! and electron density.
			call cx_set_cell(T_K(j), ne(j))

			! Each element's metal coefficients (canonical order) from the 2-D
			! photo/collisional/recombination rate arrays.
			do im = 1,n_melem
				i0  = melem_i0(im)
				top = melem_top(im)
				meg_top(im)  = top
				meg_ntot(im) = nm_tot_pp(j,im)
				meg_g0(im)   = P_m(j,i0)
				meg_b0(im)   = aion_m_pp(j,i0)
				meg_a1(im)   = rec_m_pp(j,i0+1)
				if (top >= 2) then
					meg_g1(im)  = P_m(j,i0+1)
					meg_g02(im) = P_m2(j,i0)
					meg_b1(im)  = aion_m_pp(j,i0+1)
					meg_a2(im)  = rec_m_pp(j,i0+2)
				else
					meg_g1(im)  = 0.0d0
					meg_g02(im) = 0.0d0
					meg_b1(im)  = 0.0d0
					meg_a2(im)  = 0.0d0
				endif
			enddo
			call set_metal_coeffs(n_melem, meg_ntot, meg_g0, meg_g1,   &
			                      meg_b0, meg_b1, meg_a1, meg_a2, meg_top, &
			                      meg_g02)

			! Pin the advection-corrected H/He fractions for the wrapper.
			pp_xHII_fix   = nhii(j)/nh(j)
			pp_xHeII_fix  = nheii(j)/max(nhe(j),1.0d-300)
			pp_xHeIII_fix = nheiii(j)/max(nhe(j),1.0d-300)

			! Named-field cell state for the residual: once rows 1-3 are
			! pinned the H/He rates drop out, and the metal rows consume n_h,
			! n_he and the charge-exchange temperature T_K. pp
			! runs serially on the master thread, so this master threadprivate
			! copy is the one read by ion_system_HeH_metals inside the wrapper.
			! params stays only the MINPACK transport argument (unread).
			params    = 0.0d0
			ieq_cell%nh  = nh(j)
			ieq_cell%nhe = nhe(j)
			! The temperature the residual's charge exchange is evaluated at
			! (cx_add_to_fvec) is ieq_cell%T_K, which has to be the
			! temperature cx_set_cell filled the coefficients at, this cell's.
			! Left unset it was the last cell of the equilibrium sweep, and
			! the charge-exchange guard aborted every pp_metals=2 run (found
			! 2026-10-01 on backup/regression/wasp_full with pp_metals 2).
			ieq_cell%T_K = T_K(j)

			! Initial guess: pinned H/He fractions + current metal split.
			sys_x(1) = pp_xHII_fix
			sys_x(2) = pp_xHeII_fix
			sys_x(3) = pp_xHeIII_fix
			do im = 1,n_melem
				i0 = melem_i0(im)
				sys_x(4+2*(im-1)) = nm_w(j,i0+1)/max(nm_tot_pp(j,im),1.0d-30)
				if (melem_top(im) >= 2) then
					sys_x(5+2*(im-1)) = nm_w(j,i0+2)/max(nm_tot_pp(j,im),1.0d-30)
				else
					sys_x(5+2*(im-1)) = 0.0d0
				endif
			enddo

			call hybrd1(ion_system_metals_pp,Neq_mpp,sys_x,sys_sol,   &
			            tol,info,wa,lwa_mpp,params)

			! A stage split that did not converge is not an ionization
			! balance of the cell; restore the equilibrium split, which is.
			if (info /= 1) then
				nm_w(j,:)      = nm_in(j,:)*n0
				n_metal_noconv = n_metal_noconv + 1
				! The metal stages of this row are the equilibrium split
				! while its H/He is corrected; the row's composition solve
				! did not converge and the field says so.
				adv_comp_status(j) = adv_failed
				cycle
			endif

			! Extract the re-solved stage split (element totals conserved).
			do im = 1,n_melem
				i0 = melem_i0(im)
				nm_w(j,i0+1) = nm_tot_pp(j,im)*sys_x(4+2*(im-1))
				if (melem_top(im) >= 2) then
					nm_w(j,i0+2) = nm_tot_pp(j,im)*sys_x(5+2*(im-1))
					nm_w(j,i0)   = nm_tot_pp(j,im)                          &
					             *(1.0d0 - sys_x(4+2*(im-1)) - sys_x(5+2*(im-1)))
				else
					nm_w(j,i0)   = nm_tot_pp(j,im)*(1.0d0 - sys_x(4+2*(im-1)))
				endif
			enddo
		enddo
	endif


	!---- THE COUNTS OF THE COMPOSITION THE PASS NOW HOLDS ----!
	! Every composition update of the pass is final here: the H/He sweep,
	! with its cells restored to equilibrium, and under pp_metals=2 the
	! metal re-solve, with its cells restored to the equilibrium split. The
	! metal stages carry electrons, so the free electron density, the
	! particle count and the mean molecular weight are formed again from
	! that composition before anything reads them: the ionized fraction of
	! the photoelectron partition, the heating, the energy solve (its
	! caloric energy, its thermal Damkohler number, the mean molecular
	! weight of its histories) and the output pressure. Total number
	! density and electrons include the metal nuclei and electrons under
	! eos_metals; the molecular species are excluded, the post-process does
	! not carry them (see the module-header composition note).
	call calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm_w)
	call calc_ne(nhii,nheii,nheiii,ne,nm_w)
	call calc_mmw(nh,nhe,ne,mmw,nm_w)
	! The ionized fraction of the photoelectron partition, of the same
	! composition (the definition at the top of the pass).
	xion = min(max(ne/max(nh + nhe, 1.0d-99), 0.0d0), 1.0d0)

	!---- Heating of the advection-corrected composition ----!
	! (heating_of_current_composition, below)
	call heating_of_current_composition(theat)

	
	!---- Solve stationary energy equation ----!
	! This procedure uses the same velocity profile and 
	!	the ionization profile after the advection correction
	
	! Initialize temperature at ghost cells
	T_out = T_K/T0
	
	! The mean molecular weight the solve reads is the one formed above
	! from the final composition of the pass: nm_w adds the metal mass and
	! nuclei under the eos_metals policy, so the _adv temperature solve
	! uses the same composition as the main loop.

	! Line transfer of the ground-term fine-structure lines, from the
	! incoming profile, so the cell-by-cell energy solve balances the same
	! metal cooling eval_cool reports.
	call fine_structure_line_transfer(T_K, nm_w, beta_fs_pp, nbar_fs_pp)

	! The ground-capture escape weights of the three recombinations at the
	! composition the temperature is solved at, from the one definition
	! eval_cool also uses (molecule-free: no H2 absorber on this path).
	y_HI_pp   = 0.0d0
	y_gnd_pp  = 0.0d0
	y_HeII_pp = 0.0d0
	if (use_h_rec_escape .or. (use_he_rec_coupling .and. thereis_He)) then
		nh2_pp = 0.0d0
		call ground_capture_escape_weights(nhi, nh2_pp, nhei, nheii,     &
		                                   nm_w, y_HI_pp, y_gnd_pp,       &
		                                   y_HeII_pp)
	endif

	! Count cell-by-cell temperature solves rejected as non-physical (metal modes).
	n_pp_reject = 0
	n_T_local   = 0
	Da_thermal_max   = 0.0d0
	r_Da_thermal_max = 0.0d0

	! The temperature correction solves the supplied profile's steady
	! internal-energy equation (T_equation):
	! div(rho*e*v) + p*div(v) = heating - cooling
	! = rho*v*de/dr - p*v*dln(rho)/dr + h*div(rho*v), h=e+p/rho.
	! The last term vanishes only where the state's mass flux is stationary.
	! teq_cell%div_rhov receives the signed divergence of the Riemann face
	! mass flux over the cell's spherical control volume, the same mass row
	! used by RK_rhs and stationary certification. This avoids a spurious
	! cell-centered odd-even difference at a nearly hydrostatic base.
	! A correction on a state that fails continuity is not a physical gas
	! temperature; the mass-row validity gate keeps that cell's run T.
	! The correction to the divergence and measured profile effects are in
	! ../EXHALE_v1.00/md/Update_EXHALE_stage3.md section 29.10.

	! THE TRANSPORT TERMS. A run with heat conduction, or with the enthalpy
	! flux of diffusing elements or transported carriers, solved an energy
	! equation that holds them, and in a conduction-dominated layer they
	! carry the energy (md/Update_EXHALE_stage3.md sections 40-45: 0.5 to
	! 100 times the heating on the LHS 1140 b states). The advected balance
	! (E) keeps them, evaluated with the CORRECTED temperature: conduction
	! couples each cell to both neighbors, so the column is not a marching
	! recursion any more. The marching sweep below reads the cell below at
	! its corrected value and the cell above at the profile the pass starts
	! from (a first approximation), and the column Newton solve after it
	! solves the corrected cells together, with the conductivity and the
	! enthalpy flux divergence re-evaluated at every iterate
	! (energy_column_newton). Cells that keep the run's own temperature are
	! fixed values of that solve.
	pp_transport_on = conduction_active() .or.                            &
	                  interdiffusion_enthalpy_active() .or.               &
	                  associated(carrier_enthalpy_divergence)
	pp_newton_mode  = .false.
	pp_bdf2_used    = .false.
	! The column outcome is of THIS pass: a pass without a column to solve
	! says so, whatever the pass before it did.
	col_outcome     = adv_col_not_run
	col_reason      = adv_col_reason_none
	col_cell        = 0
	col_unknowns    = 0
	n_col_newton_it = 0
	col_worst_cell  = 0
	col_worst_R     = 0.0d0
	col_worst_rel   = 0.0d0
	col_worst_tol_rel = 0.0d0
	n_col_unrooted  = 0
	if (pp_transport_on) then
		pp_T_nbr = T_out
		call transport_terms_of_profile(T_out)
	endif

	do j = 3-Ng,N+Ng ! Start from first computational cell

		vm = v(j-1)
		vp = v(j)
		! The mass row is above the fraction a corrected row is accurate to
		! (condition (iv) of the validity block above): the steady equation
		! solved here would drop terms larger than that fraction of the ones
		! it keeps, so the cell keeps the run's own temperature. Tested
		! before the inflow sign, matching the order the status is built in:
		! it refuses more of the row than the conditions that follow.
		if (.not. mass_row_within_tol(j)) then
			T_out(j)        = T_in(j)
			adv_T_status(j) = adv_retained
			cycle
		endif
		! The energy balance of the cell is carried by channels this
		! post-process omits (the energy class, formed with the stationarity
		! condition above): the row keeps the run's own temperature, and the
		! field says the closure does not cover it.
		if (.not. energy_class_modelled(j) .or. .not. cell_class_modelled(j)) then
			T_out(j)        = T_in(j)
			adv_T_status(j) = adv_unsupported
			cycle
		endif
		! Inflow (condition (i) of the ionization validity block above): the
		! cell keeps the converged eq temperature. The energy step takes its
		! history from the cells below, just like the ionization step, so it
		! is invalid wherever the gas moves inward, and it would otherwise
		! land on the spurious hot root that then cascades up. This is a
		! property of the discretization, so it does not depend on the metal
		! switch. Conditions (ii) and (iii) of that block are statements about
		! the ionization balance and its representation and are deliberately
		! NOT applied here; the energy solve has the Damkohler condition of its
		! own that follows, formed from the terms of ITS equation.
		if (vp <= 0.0d0 .or. vm <= 0.0d0) then
			T_out(j)        = T_in(j)
			adv_T_status(j) = adv_retained
			cycle
		endif
		dr = r(j) - r(j-1)

		! Evaluate Da with updated-composition heating, start-of-pass cooling,
		! and the current temperature and particle count.
		! For Da>1 the code retains T_in instead of solving the advected
		! energy equation. This net-rate test does not establish a local
		! heat=cool root; near balance, it may miss fast thermal relaxation.
		! Decide before solving the corrected temperature.
		! To Be Checked: comment and code disagree, see md/comment_revision_20261027.md
		t_cross_E = dr*R0/(vp*v0)
		q_rad     = abs(theat(j)*q0 - tcool_in(j))
		u_th      = (n_tot(j) + ne(j))*kb_erg*T0                          &
		            *internal_energy_of_mixture(h2_particle_fraction(j),  &
		                                        T_out(j))
		Da_thermal = thermal_damkohler_number(t_cross_E, q_rad, u_th)
		if (Da_thermal > Da_thermal_max) then
			Da_thermal_max   = Da_thermal
			r_Da_thermal_max = r(j)
		endif
		if (Da_thermal > 1.0d0) then
			T_out(j)        = T_in(j)
			n_T_local       = n_T_local + 1
			adv_T_status(j) = adv_retained
			cycle
		endif

		! The step (bdf2_step_ratio_limit): BDF2 where the temperature of the
		! cell below was itself integrated in this pass, backward Euler at the
		! start of the recursion and after a cell that kept the run's own
		! temperature. A BDF2 step that finds no admissible root is retaken
		! by backward Euler before the cell falls back to the run's own
		! temperature.
		use_bdf2 = bdf2_admissible(j) .and. adv_T_status(j-1) == adv_corrected
		call energy_step(j, use_bdf2, sys_x_T(1), T_outcome)
		if (use_bdf2 .and. T_outcome >= T_no_bracket) then
			n_bdf2_T_retaken = n_bdf2_T_retaken + 1
			use_bdf2 = .false.
			call energy_step(j, use_bdf2, sys_x_T(1), T_outcome)
		endif
		if (use_bdf2) then
			n_bdf2_T = n_bdf2_T + 1
		else
			n_be_T = n_be_T + 1
		endif
		pp_bdf2_used(j) = use_bdf2

		! The ledger of how the cell solve ended (energy_step), and the run's
		! own temperature wherever it gave no root of this gas.
		select case (T_outcome)
		case (T_root_at_floor)
			n_T_res_root = n_T_res_root + 1
		case (T_not_converged)
			n_T_noconv   = n_T_noconv + 1
		case (T_no_bracket, T_non_positive, T_out_of_band)
			n_pp_reject  = n_pp_reject + 1
		end select
		T_solve_fell_back = (T_outcome >= T_no_bracket)
		if (T_solve_fell_back) sys_x_T(1) = T_in(j)

		! The cell solve did not converge, or its root was not a temperature
		! of this gas, so the row carries the run's own temperature.
		if (T_solve_fell_back) then
			adv_T_status(j) = adv_failed
		else
			adv_T_status(j) = adv_corrected
		endif

		! Extract solution profiles
		T_out(j) = sys_x_T(1)

	enddo

	! The corrected cells together, with the transport terms of the
	! corrected temperature (the statement before the marching sweep).
	if (pp_transport_on) call energy_column_newton()

	! Update pressure and temperature
	p_out = (n_tot + ne)/n0*T_out
	T_K = T_out*T0
      
	enddo ! End loop on post processing

	! THE CELL CLASS THE CLOSURE DOES NOT COVER, over both fields. Formed on
	! the state handed in (the block before the validity conditions) and
	! applied last, because it is a statement about what this post-process
	! models and stands above the conditions that name the flow and the
	! solve: in such a cell the particle and electron counts every equation
	! here was solved with are not the counts of that gas, so neither the
	! temperature nor the composition of the row describes it.
	do j = 1-Ng,N+Ng
		if (.not. cell_class_modelled(j)) then
			adv_T_status(j)    = adv_unsupported
			adv_comp_status(j) = adv_unsupported
		endif
	enddo

	! Report how many cells the advection correction was not applied to
	! (inflow, local ionization equilibrium, an unrepresentable ion fraction,
	! or a mass row above the fraction a corrected row is accurate to; see
	! the validity block above).
	! Counted on the last pass.
	if (n_adv_eq > 0) then
		write(*,'(a,i0,a,i0,a)') ' (post_process_adv) advection correction: ', &
		   n_adv_eq, ' of ', N+2*Ng-1,                                          &
		   ' cells kept at the equilibrium ionization.'
	endif

	! Report cells retained by the net-rate thermal Da>1 test and its maximum.
	! This test does not establish a local thermal-equilibrium root.
	write(*,'(a,i0,a,i0,a,es9.2,a,f7.4,a)')                                  &
	   ' (post_process_adv) energy correction: ', n_T_local, ' of ',          &
	   N+2*Ng-2, ' cells kept the run temperature (thermal Damkohler > 1);'   &
	   //' largest Damkohler ', Da_thermal_max, ' at r = ',                   &
	   r_Da_thermal_max, ' Rp.'

	! Report how far from stationary the state handed in is, on the mass row
	! the stationary certification measures, together with the fraction a
	! corrected row is accurate to and the count of cells above it. The
	! largest measure is reported whether or not any cell crossed the
	! fraction, so the margin of the condition on this state is in the log,
	! and the second count says how many cells the condition adds beyond the
	! ionization conditions (i) to (iii), which is what it is worth as a
	! condition of its own. Every row's own measure is in the adv_mass_row
	! column of the file, so a reader is not left with the maximum alone.
	write(*,'(a,es9.2,a,f7.4,a,es9.2,a,i0,a,i0,a,i0,a)')                     &
	   ' (post_process_adv) conditional correction: mass row |R_1|/s_1 of'    &
	   //' the input state up to ', mass_row_max, ' at r = ',                 &
	   r_mass_row_max, ' Rp; above ', adv_conditional_tol,                    &
	   ', the fraction a corrected row is accurate to, in ',                  &
	   n_not_stationary, ' of ', N+2*Ng-1, ' cells, of which ',               &
	   n_stationarity_only, ' the ionization conditions do not already refuse.'
	! Whether that row was measured on the very state the run last assembled
	! one for: the primitive-to-conservative round trip through the caloric
	! EOS is not the identity to the last bit, and a reader of the log should
	! not have to assume it is.
	if (.not. terms_are_run_state)                                           &
		write(*,'(a)') ' (post_process_adv) stationarity: the mass row was'   &
		   //' assembled from the primitive state handed in, which is the'    &
		   //' run''s state re-formed and not its conserved variables.'

	! Report the SENSITIVITY of the energy correction to the same
	! non-stationarity: the size of the enthalpy flux of the mass-flux
	! divergence against the two terms the steady balance keeps. It refuses
	! nothing (enthalpy_flux_term_ratio says why); the count is of cells in
	! which that term is the largest in the equation, so the temperature
	! returned there is set by it.
	write(*,'(a,es9.2,a,f7.4,a,i0,a,i0,a,es9.2,a)')                          &
	   ' (post_process_adv) energy correction: enthalpy flux of the'          &
	   //' mass-flux divergence / other terms up to ',                        &
	   enthalpy_flux_ratio_max,                                              &
	   ' at r = ', r_enthalpy_flux_max, ' Rp; dominant in ',                  &
	   n_enthalpy_term_dominant, ' of ', N+2*Ng-1,                            &
	   ' cells (reported, not a refusal; level ',                             &
	   enthalpy_ratio_report_level, ').'

	! Report the cells in which the balance of a solved population is
	! carried more by the omitted molecular channels than by those the
	! advection systems hold (the composition class).
	if (n_species_class_refused > 0)                                        &
		write(*,'(a,i0,a,i0,a)') ' (post_process_adv) composition class: ',  &
		   n_species_class_refused, ' of ', N+2*Ng,                         &
		   ' cells lose H+, He+ or He(2^3S) faster to H2 than to the'//     &
		   ' channels the advection systems carry; they keep the run'//     &
		   ' composition and temperature.'

	! Report the cells whose energy balance the post-process does not carry
	! (the energy class, formed on the state handed in).
	if (n_energy_class_refused > 0)                                         &
		write(*,'(a,i0,a,i0,a,f8.4,a)') ' (post_process_adv) energy class: ',&
		   n_energy_class_refused, ' of ', N+2*Ng,                            &
		   ' cells have more heating and cooling in channels the'//         &
		   ' post-process omits than in those it carries (up to r = ',       &
		   r_energy_class_top, ' Rp); they keep the run temperature.'

	! Report the steps the two recursions took on the last pass
	! (bdf2_step_ratio_limit): BDF2 steps, backward Euler steps (the start of
	! the recursion and every restart after a cell it did not integrate),
	! and the BDF2 steps retaken by backward Euler for a negative population
	! or no admissible temperature root.
	write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)')                              &
	   ' (post_process_adv) steps of the last pass: composition ',          &
	   n_bdf2_comp, ' BDF2 and ', n_be_comp, ' backward Euler (',           &
	   n_bdf2_comp_retaken, ' retaken from BDF2); energy ', n_bdf2_T,       &
	   ' BDF2 and ', n_be_T, ' backward Euler (', n_bdf2_T_retaken,         &
	   ' retaken from BDF2).'

	! Report no-bracket, nonpositive, and out-of-band cell outcomes.
	if (n_pp_reject > 0) then
		write(*,'(a,i0,a,i0,a)') ' (post_process_adv) T solve: ',            &
		   n_pp_reject, ' of ', N+2*Ng-2,                                     &
		   ' cells fell back to eq T (no bracket, non-positive or out-of-band root).'
	endif

	! Report the cell solves that did not converge and were therefore left at
	! the equilibrium state, summed over the passes (see the declaration).
	if (n_adv_noconv > 0)                                                   &
		write(*,'(a,i0,a)') ' (post_process_adv) advection system: ',         &
		   n_adv_noconv, ' cell solves did not converge and kept the'         &
		   //' equilibrium ionization.'
	if (n_metal_noconv > 0)                                                 &
		write(*,'(a,i0,a)') ' (post_process_adv) metal re-solve: ',           &
		   n_metal_noconv, ' cell solves did not converge and kept the'       &
		   //' equilibrium stage split.'
	if (n_T_noconv > 0)                                                     &
		write(*,'(a,i0,a)') ' (post_process_adv) energy equation: ',          &
		   n_T_noconv, ' cell solves did not converge and kept the'           &
		   //' equilibrium temperature.'
	if (n_T_res_root > 0)                                                   &
		write(*,'(a,i0,a)') ' (post_process_adv) energy equation: ',          &
		   n_T_res_root, ' cell solves stopped improving at a residual'       &
		   //' already at the cancellation floor of their own terms; the'     &
		   //' iterate is the root and was kept.'

      
	!---- Update cooling rates ----!

	! Molecule-free electron sum and no H3+ cooling, as at the first eval_cool
	! call above. T_equation, which solved for this T_K, balanced the same
	! assembly (Cool_coeff: radiative_cooling_of_cell) at the same densities;
	! the fine-structure escape probabilities it held fixed are those of the
	! profile the pass started from, which eval_cool re-forms at this T_K.
	call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm_w,            &
	  			   dum_v1,dum_v2,dum_v3, rec_m_pp,                   &
				   dum_v4,dum_v5,dum_v6, aion_m_pp,                      &
				   tcool, nheiTR=nheiTR)

	! Adimensionalize
	tcool = tcool/q0

 
   ! Adimensionalize ion densities before writing
   nhi_w    = nhi/n0
   nhii_w   = nhii/n0
   nhei_w   = nhei/n0
   nheii_w  = nheii/n0
   nheiii_w = nheiii/n0
   nheiTR_w = nheiTR/n0

      
   ! THE DERIVED STATE (adv_derived_state, output_write): how the solves
   ! that produced the rows about to be written ended, as a whole. The
   ! column fields are those of the last pass, which is the pass the rows
   ! are; the row counts are of the two fields as written. The verdict is
   ! a rejection when the column of the last pass was rejected or any row
   ! is failed; otherwise the product is unverified, never accepted.
   ! WHAT IS TESTED: the energy row of every column unknown at the final
   ! profile of the last pass (energy_column_newton), with the transport
   ! coefficients formed there -- but the conductivity and the enthalpy
   ! divergences of the INPUT composition f_sp_in, and the fine-structure
   ! escape probabilities of the profile the pass started from: lagged
   ! coefficients, not those of the derived composition. WHAT IS NOT: the
   ! full energy residual of the state on one derived composition with
   ! every coefficient its own, the convergence of the outer iteration
   ! (a fixed number of passes, outer=unverified), and the momentum
   ! balance of the corrected pressure (rho and v are written back
   ! unchanged).
   derived_state%recorded             = .true.
   derived_state%n_chem_solves_failed = n_adv_noconv + n_metal_noconv
   if (any(adv_comp_status .eq. adv_failed)) then
      derived_state%chemistry = adv_chem_cells_failed
   else
      derived_state%chemistry = adv_chem_complete
   endif
   derived_state%column                 = col_outcome
   derived_state%column_reason          = col_reason
   derived_state%column_iterations      = n_col_newton_it
   derived_state%column_cell            = col_cell
   derived_state%column_unknowns        = col_unknowns
   derived_state%column_worst_cell         = col_worst_cell
   derived_state%column_worst_residual     = col_worst_R
   derived_state%column_worst_residual_rel = col_worst_rel
   derived_state%cells_marching_unrooted   = n_col_unrooted
   derived_state%passes                 = n_pp_passes
   derived_state%column_passes_rejected = n_col_passes_rejected
   ! Physical cells 1..N only: the ghost rows are written but are not
   ! cells of the column, and the census does not count them.
   derived_state%cells_failed      = count(adv_T_status(1:N) .eq. adv_failed &
                                     .or. adv_comp_status(1:N) .eq. adv_failed)
   derived_state%cells_retained    = count(adv_T_status(1:N) .eq. adv_retained &
                                     .or. adv_comp_status(1:N) .eq. adv_retained)
   derived_state%cells_unsupported = count(adv_T_status(1:N) .eq. adv_unsupported &
                                     .or. adv_comp_status(1:N) .eq. adv_unsupported)
   derived_state%outer             = adv_outer_unverified
   ! A failed ghost row is not a failed cell, but it can only be failed as
   ! an unknown of a rejected column, which the first clause already holds.
   derived_state%rejected = (col_outcome .eq. adv_col_rejected) .or.        &
                            (derived_state%cells_failed .gt. 0)

   ! Write updated thermodynamic and ionization profiles. Metals are written
   ! in dimensionless (n0) units, matching the other species; zero in the
   ! metal-free mode, frozen/re-solved eq densities otherwise.
   nm_out = nm_w/n0

   ! A REFUSED ROW CARRIES THE RUN'S OWN STATE, as the file legend states,
   ! and it carries the very numbers the run wrote, not a round trip of
   ! them through the cgs arrays of this routine: a composition field that
   ! is not corrected (retained, failed, unsupported, not evaluated) is the
   ! handed-in composition in every species column, the metal stages
   ! included whatever pp_metals is; a temperature field that is retained,
   ! unsupported or not evaluated is the handed-in temperature; and a row
   ! whose two fields are both refused carries the handed-in pressure,
   ! heating and cooling as well, since nothing of it was solved here. The
   ! one exception is the temperature of a failed (2) row of a rejected
   ! column energy solve, which is the marching profile of that pass
   ! (the legend says so).
   do j = 1-Ng,N+Ng
      if (adv_comp_status(j) .ne. adv_corrected) then
         nhi_w(j)    = nhi_in(j)
         nhii_w(j)   = nhii_in(j)
         nhei_w(j)   = 0.0d0
         nheii_w(j)  = 0.0d0
         nheiii_w(j) = 0.0d0
         nheiTR_w(j) = 0.0d0
         if (thereis_He) then
            nhei_w(j)   = nhei_in(j)
            nheii_w(j)  = nheii_in(j)
            nheiii_w(j) = nheiii_in(j)
            if (thereis_HeITR) nheiTR_w(j) = nheiTR_in(j)
         endif
         if (thereis_metals) nm_out(j,:) = nm_in(j,:)
      endif
      if (adv_T_status(j) .eq. adv_retained .or.                        &
          adv_T_status(j) .eq. adv_unsupported .or.                     &
          adv_T_status(j) .eq. adv_not_evaluated) T_out(j) = T_in(j)
      if (adv_T_status(j) .ne. adv_corrected .and.                      &
          adv_T_status(j) .ne. adv_failed .and.                         &
          adv_comp_status(j) .ne. adv_corrected) then
         p_out(j) = p(j)
         theat(j) = heat(j)
         tcool(j) = cool(j)
      endif
   enddo
   ! The two status fields say, row by row, whether the temperature and the
   ! composition of that row are the steady correction or the run's own
   ! state, so a reader of the file (and of a spectrum built from it) can
   ! tell which rows carry which.
   call write_output(rho,v,p_out,T_out,theat,tcool,eta,                &
                     nhi_w,nhii_w,nhei_w,nheii_w,nheiii_w,              &
                     nheiTR_w,nm_out,'ad', adv_T_status, adv_comp_status,  &
                     adv_mass_row, derived_state)

	contains

	! Ionization state of one cell as the equilibrium solver left it: the
	! densities this post-process was handed. Both the validity conditions and
	! a non-converged cell fall back to it, so it is written once.
	subroutine pin_cell_to_equilibrium(jc)
	integer, intent(in) :: jc

	nhi(jc)    = nhi_in(jc)*n0
	nhii(jc)   = nhii_in(jc)*n0
	nheii(jc)  = nheii_in(jc)*n0
	nheiii(jc) = nheiii_in(jc)*n0
	nheiTR(jc) = 0.0
	if (thereis_HeITR) nheiTR(jc) = nheiTR_in(jc)*n0
	nheiS(jc)  = he_ground_singlet_density(nhei_in(jc)*n0, nheiTR(jc))
	nhei(jc)   = nheiS(jc) + nheiTR(jc)

	end subroutine pin_cell_to_equilibrium

	! Heating of the current composition [code units]. Use the shared
	! heating_of_composition assembly with the attenuated field recomputed
	! from this pass, subject to the molecular and oxygen omissions below.
	!
	! The composition reconstructed here carries no molecular and no oxygen
	! carriers (see the header of this module), which is what the two
	! composition flags below say; the molecular and oxygen deposits are
	! therefore absent from the _adv heating, as are the molecular carriers
	! from its n_e and its n_tot.
	subroutine heating_of_current_composition(heat_out)
	real*8, dimension(1-Ng:N+Ng), intent(out) :: heat_out
	! The metal photoionization rates the field call also returns; the
	! pass's own P_m (with the recombination photons added) is not touched.
	real*8, dimension(1-Ng:N+Ng,n_mion) :: P_m_of_field

	if (thereis_He) then
		call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm_w, xion,     &
		                 dum_v1,dum_v2,dum_v3,dum_v4, P_m_of_field, &
		                 dum_v6,dum_v5,                             &
		                 heat_of_one_HI   = h1_HI_pp,               &
		                 heat_of_one_HeI  = h1_HeI_pp,              &
		                 heat_of_one_HeII = h1_HeII_pp,             &
		                 heat_of_one_HeTR = h1_HeTR_pp,             &
		                 heat_of_one_H2   = h1_H2_pp,               &
		                 heat_of_one_mion = h1_m_pp,                &
		                 nh2_solved = nmol_eq(:,1))
  	else
	  	call PH_heat_H(nhi, xion, dum_v1,dum_v6,dum_v2,             &
	  	               heat_of_one_HI = h1_HI_pp)
		h1_HeI_pp  = 0.0d0
		h1_HeII_pp = 0.0d0
		h1_HeTR_pp = 0.0d0
		h1_H2_pp   = 0.0d0
		h1_m_pp    = 0.0d0
  	endif

	nmol_pp = 0.0d0
	nox_pp  = 0.0d0
	k_lw_pp = 0.0d0
	! No CO in the reconstructed composition, so no CO photodissociation
	! rate: the two CO channels of the assembly are gated on with_oxygen,
	! which is .false. here, and this array is what they would contract.
	k_co_pp = 0.0d0
	p_lw_pp = 1.0d0
	j_fuv_pp = 0.0d0
	call heating_of_composition(T_K,                                      &
	         nhi,nhii,nhei,nheii,nheiii,nheiTR, nm_w, nmol_pp, nox_pp,    &
	         ne, n_tot,                                                   &
	         h1_HI_pp,h1_HeI_pp,h1_HeII_pp,h1_HeTR_pp,h1_H2_pp,h1_m_pp,   &
	         A31,q31a,q31b,Q31,                                           &
	         k_lw_pp, p_lw_pp, k_co_pp, j_fuv_pp, j_fuv_pp,               &
	         .false., .false., heat_out, heat_chan_pp,                    &
	         nh2_solved=nmol_eq(:,1))

	! Adimensionalize
	heat_out = heat_out/q0

	end subroutine heating_of_current_composition

	! THE WIDTH THE RIGHT-HAND SIDE OF THE STEP INTO CELL jc IS MULTIPLIED
	! BY [code units]: g h_j for BDF2, h_j for backward Euler
	! (bdf2_step_ratio_limit).
	real*8 function step_width(jc, bdf2)
	integer, intent(in) :: jc
	logical, intent(in) :: bdf2
	step_width = r(jc) - r(jc-1)
	if (bdf2) step_width = bdf_g(jc)*step_width
	end function step_width

	! The history of the fraction num/den for the step into cell jc,
	! a1 x_{j-1} - a2 x_{j-2} (BDF2) or x_{j-1} (backward Euler). Cell jc-2
	! is read only by the BDF2 step.
	real*8 function history_of_fraction(num, den, jc, bdf2)
	real*8, dimension(1-Ng:N+Ng), intent(in) :: num, den
	integer, intent(in) :: jc
	logical, intent(in) :: bdf2
	if (bdf2) then
		history_of_fraction = bdf_a1(jc)*num(jc-1)/den(jc-1)                &
		                    - bdf_a2(jc)*num(jc-2)/den(jc-2)
	else
		history_of_fraction = num(jc-1)/den(jc-1)
	endif
	end function history_of_fraction

	! SPECIFIC internal energy of cell kc at the code temperature T_code and
	! mean molecular weight mu, E(x_H2,kc, T)/mu, evaluated at the H2 share of
	! ITS OWN cell (see the e_hist field of ion_cell_state).
	real*8 function specific_internal_energy(kc, T_code, mu)
	integer, intent(in) :: kc
	real*8,  intent(in) :: T_code, mu
	specific_internal_energy =                                          &
	     internal_energy_of_mixture(h2_particle_fraction(kc), T_code)/mu
	end function specific_internal_energy

	! ONE STEP OF THE IONIZATION RECURSION INTO CELL jc: the history of the
	! step and its rate weights c1 = (g) h_j/v_H and c1_he = (g) h_j/v_He,
	! each element's rows on the velocity of its own nuclei, then the
	! cell's advection system solved from the fractions of the last pass.
	! The rates of the cell are loaded by the caller; the result is left in
	! sys_x and info.
	subroutine composition_step(jc, bdf2)
	integer, intent(in) :: jc
	logical, intent(in) :: bdf2

	adv_cell%c1       = step_width(jc, bdf2)*R0/(v_nuc_H(jc)*v0)
	adv_cell%c1_he    = step_width(jc, bdf2)*R0/(v_nuc_He(jc)*v0)
	adv_cell%xhi_hist = history_of_fraction(nhi, nh, jc, bdf2)
	sys_x(1) = nhi(jc)/nh(jc)
	if (.not. thereis_He) then
		call hybrd1(adv_implicit_H,Neq_adv,sys_x,sys_sol,               &
		            tol,info,wa,lwa_adv,params)
		return
	endif
	adv_cell%xheiS_hist  = history_of_fraction(nheiS,  nhe, jc, bdf2)
	adv_cell%xheiii_hist = history_of_fraction(nheiii, nhe, jc, bdf2)
	sys_x(2) = nheiS(jc)/nhe(jc)
	sys_x(3) = nheiii(jc)/nhe(jc)
	if (thereis_HeITR) then
		adv_cell%xheiTR_hist = history_of_fraction(nheiTR, nhe, jc, bdf2)
		sys_x(4) = nheiTR(jc)/nhe(jc)
		call hybrd1(adv_implicit_HeH_TR,Neq_adv,sys_x,sys_sol,          &
		            tol,info,wa,lwa_adv,params)
	else
		call hybrd1(adv_implicit_HeH,Neq_adv,sys_x,sys_sol,             &
		            tol,info,wa,lwa_adv,params)
	endif
	end subroutine composition_step

	! Whether the fractions a composition step returned in sys_x are
	! populations: none negative, the H fractions summing to one and the He
	! fractions (the once-ionized one being the rest) to one.
	logical function populations_are_admissible()
	real*8 :: x_heii_rest
	populations_are_admissible = (sys_x(1) >= 0.0d0 .and. sys_x(1) <= 1.0d0)
	if (.not. thereis_He) return
	x_heii_rest = 1.0d0 - sys_x(2) - sys_x(3)
	if (thereis_HeITR) x_heii_rest = x_heii_rest - sys_x(4)
	populations_are_admissible = populations_are_admissible .and.       &
	     sys_x(2) >= 0.0d0 .and. sys_x(3) >= 0.0d0 .and.                 &
	     x_heii_rest >= 0.0d0
	if (thereis_HeITR) populations_are_admissible =                     &
	     populations_are_admissible .and. sys_x(4) >= 0.0d0
	end function populations_are_admissible

	! ONE STEP OF THE ENERGY RECURSION INTO CELL jc (the equation is stated
	! before the energy loop, the residual in T_equation). Returns the root
	! in code units and how the solve ended; T_root is the run's own
	! temperature whenever outcome is T_no_bracket or above.
	subroutine energy_step(jc, bdf2, T_root, outcome)
	integer, intent(in)  :: jc
	logical, intent(in)  :: bdf2
	real*8,  intent(out) :: T_root
	integer, intent(out) :: outcome
	real*8  :: xT(1), fT(1)
	real*8  :: e_cell_root
	logical :: ok_brent

	! The width of the step and the histories it differences against. The
	! specific energy of each history cell is formed at ITS OWN composition
	! and the temperature this loop left it with, which is the quantity the
	! flow carries in (the e_hist field of ion_cell_state).
	dr_step = step_width(jc, bdf2)
	if (bdf2) then
		rho_hist = bdf_a1(jc)*rho(jc-1) - bdf_a2(jc)*rho(jc-2)
		e_hist   = bdf_a1(jc)*specific_internal_energy(jc-1, T_out(jc-1),   &
		                                               mmw(jc-1))           &
		         - bdf_a2(jc)*specific_internal_energy(jc-2, T_out(jc-2),   &
		                                               mmw(jc-2))
	else
		rho_hist = rho(jc-1)
		e_hist   = specific_internal_energy(jc-1, T_out(jc-1), mmw(jc-1))
	endif

	teq_cell%nhi    = nhi(jc)
	teq_cell%nhii   = nhii(jc)
	teq_cell%nheiS  = nheiS(jc)
	teq_cell%nheiTR = nheiTR(jc)
	teq_cell%y_HI   = y_HI_pp(jc)
	teq_cell%y_gnd  = y_gnd_pp(jc)
	teq_cell%y_HeII = y_HeII_pp(jc)
	teq_cell%nheii  = nheii(jc)
	teq_cell%nheiii = nheiii(jc)
	teq_cell%mup    = mmw(jc)
	teq_cell%mum    = mmw(jc-1)
	teq_cell%rhov   = rho(jc)*v(jc)
	teq_cell%coeff  = mmw(jc-1)*v(jc)*(rho(jc) - rho_hist)
	teq_cell%dr_step = dr_step
	teq_cell%heaold = theat(jc)
	! Composition entry of the caloric EOS: the H2 share of the
	! particle-plus-electron count of this cell, as the equilibrium solve
	! left it.
	teq_cell%x_h2   = h2_particle_fraction(jc)
	teq_cell%e_hist = e_hist
	! Coefficient of the enthalpy flux term, the mass row of the state
	! (see the statement before the energy loop). It does not change while
	! the root finder varies the temperature.
	teq_cell%div_rhov = div_rhov_state(jc)
	call set_transport_terms_of_cell(jc)
	! Metal densities for this cell [cgs] go through the equation_T module
	! array (the 27-ion vector does not fit params). pp_metal_on gates
	! whether T_equation adds the metal cooling/brem/n_e terms.
	pp_nm_cell(:)  = nm_w(jc,:)
	pp_beta_fs(:)  = beta_fs_pp(jc,:)
	pp_nbar_fs(:)  = nbar_fs_pp(jc,:)

	outcome = T_root_found
	T_root  = T_in(jc)

	! With metals and use_brent_tsolve, solve_T_brent searches upward for
	! the first bracket in [max(0.05*T_in,1 K),4*T_in] to avoid a hot root
	! of the non-monotone residual. This is the first root in that search
	! range, not a proof of physical uniqueness. "Brent solver: False"
	! selects MINPACK and the [0.5,2]*T_in rejection band. Metals-off uses
	! MINPACK without that temperature band.
	if (pp_metal_on .and. use_brent_tsolve) then
		call solve_T_brent(paramsT, T_in(jc), xT(1), ok_brent)
		if (.not. ok_brent) then
			outcome = T_no_bracket
			return
		endif
		T_root = xT(1)
		return
	endif

	xT(1) = T_out(jc)
	call hybrd1(T_equation,1,xT,fT,tol,info,wa_T,8,paramsT)
	! WHAT info MEANS AND WHAT THE ROOT IS. MINPACK's info states how its
	! ITERATION ended, not whether the iterate is a root: info = 4 and 5 are
	! returned when the steps stop improving the residual, which is what a
	! stalled search and an ARRIVED one look like alike. The root of a scalar
	! equation is defined by its residual, so the iterate is kept whenever
	! that residual is negligible against the largest term the equation
	! holds, and only an iterate that is not a root is refused.
	!
	! THE SCALE. The terms of the residual as T_equation assembles them: the
	! advected energy of this cell and its history, the enthalpy flux of the
	! mass divergence, the compression work, and the photoheating. A sum of
	! terms of size s cannot be formed to better than a few machine epsilons
	! of s, so 1e2*epsilon(s) is the level at which the equation is an
	! identity in double precision. MEASURED on the three cells of the
	! LHS 1140 b 45 Rp wind that reached this branch under the upwind step
	! (before 2026-09-27): |R|/s = 1.5e-17, 2.8e-17 and 5.2e-17.
	!
	! WHY IT MATTERS MORE THAN ONE CELL. The corrected temperature of a cell
	! is the history of the cells above it, so a discarded root is not a
	! local blemish: every row above it is integrated from a different
	! starting value.
	if (info /= 1) then
		e_cell_root = internal_energy_of_mixture(teq_cell%x_h2, xT(1))
		if (abs(fT(1)) <= 1.0d2*epsilon(1.0d0)*max(                      &
		       abs(teq_cell%mum*teq_cell%rhov*e_cell_root),                &
		       abs(teq_cell%mup*teq_cell%mum*teq_cell%rhov*e_hist),        &
		       abs(teq_cell%mum*dr_step*teq_cell%div_rhov                  &
		           *(e_cell_root + xT(1))),                                &
		       abs(teq_cell%coeff*xT(1)),                                  &
		       abs(teq_cell%mup*teq_cell%mum*dr_step*theat(jc)),            &
		       abs(teq_cell%mup*teq_cell%mum*dr_step*(teq_cell%heat_extra     &
		           + teq_cell%cond_diag*xT(1))))) then
			outcome = T_root_at_floor
		else
			outcome = T_not_converged
			return
		endif
	endif
	! A non-positive root is not a temperature, whatever else is in the
	! gas, so this test is not conditional on the metals. It matters
	! because T_out feeds the NEXT post-process pass: eval_cool takes
	! sqrt(T/T0) in the Badnell recombination fit (rr_badnell), so a
	! negative T there is a NaN cooling rate in an ordinary build and an
	! abort under -ffpe-trap=invalid (MEASURED on the He/H = 1 molecular
	! case, metals-off: cells 278-280 at -42, -640 and -2474 K on the second
	! pass under the upwind step).
	if (.not. (xT(1) > 0.0d0)) then
		outcome = T_non_positive
		return
	endif
	! With metal cooling, additionally refuse an out-of-band root: the
	! metal-cooled residual is non-monotone and carries a second, spurious
	! HOT root that hybrd1 can land on.
	if (pp_metal_on) then
		if (xT(1) > 2.0d0*T_in(jc) .or. xT(1) < 0.5d0*T_in(jc)) then
			outcome = T_out_of_band
			return
		endif
	endif
	T_root = xT(1)
	end subroutine energy_step


	subroutine transport_terms_of_profile(Tprof)
	! The conduction triplets of the temperature profile Tprof [code units]
	! (thermal_conduction_coeffs, the operator of the run's own energy
	! equation, conductivity of the composition the post-process was handed)
	! and the enthalpy flux divergence of the element and carrier fluxes of
	! the same profile, entered as heating (the run's energy row ADDS the
	! divergence, steady_residual). pp_T_bath is the base-level temperature
	! the operator holds in the ghost (conduction_base_level_T).
	real*8, dimension(1-Ng:N+Ng), intent(in) :: Tprof
	real*8, dimension(1-Ng:N+Ng) :: s_elem, s_carr
	pp_heat_rel  = 0.0d0
	pp_heat_elem = 0.0d0
	pp_heat_carr = 0.0d0
	if (conduction_active()) then
		call thermal_conduction_coeffs(Tprof, f_sp_in, pp_cond_lo,        &
		                               pp_cond_di, pp_cond_up)
		pp_T_bath = conduction_base_level_T(Tprof)
	else
		pp_cond_lo = 0.0d0;  pp_cond_di = 0.0d0;  pp_cond_up = 0.0d0
		pp_T_bath  = 0.0d0
	endif
	if (interdiffusion_enthalpy_active()) then
		call interdiffusion_enthalpy_divergence_of_state(rho, Tprof,      &
		                                                 f_sp_in, s_elem)
		pp_heat_elem = -s_elem
		pp_heat_rel  = pp_heat_rel - s_elem
	endif
	if (associated(carrier_enthalpy_divergence)) then
		call carrier_enthalpy_divergence(rho, Tprof, f_sp_in, s_carr)
		pp_heat_carr = -s_carr
		pp_heat_rel  = pp_heat_rel - s_carr
	endif
	end subroutine transport_terms_of_profile


	subroutine set_transport_terms_of_cell(jc)
	! The transport heating of cell jc for the energy residual: the
	! conduction of its two neighbors and the enthalpy flux divergence as
	! the fixed part, its own conduction coefficient as the part that
	! varies with its temperature (teq_state). The neighbor above is the
	! pass's profile in the marching sweep and the current iterate in the
	! column solve; the neighbor below is the current value either way (the
	! base bath for cell 1).
	integer, intent(in) :: jc
	real*8 :: T_below, T_above
	teq_cell%heat_extra = 0.0d0
	teq_cell%cond_diag  = 0.0d0
	if (.not. pp_transport_on) return
	if (jc .lt. 1 .or. jc .gt. N) return
	if (jc .eq. 1) then
		T_below = pp_T_bath
	else
		T_below = T_out(jc-1)
	endif
	if (pp_newton_mode) then
		T_above = T_out(jc+1)
	else
		T_above = pp_T_nbr(jc+1)
	endif
	teq_cell%heat_extra = pp_cond_lo(jc)*T_below + pp_cond_up(jc)*T_above &
	                    + pp_heat_rel(jc)
	teq_cell%cond_diag  = pp_cond_di(jc)
	end subroutine set_transport_terms_of_cell


	real*8 function energy_residual_of_cell(jc)
	! The residual of the energy equation of corrected cell jc at the
	! current profile T_out: the setup energy_step makes (histories, step,
	! densities, transport terms) and T_equation at T_out(jc).
	integer, intent(in) :: jc
	real*8  :: xv(1), fv(1)
	integer :: iflag_r
	dr_step = step_width(jc, pp_bdf2_used(jc))
	if (pp_bdf2_used(jc)) then
		rho_hist = bdf_a1(jc)*rho(jc-1) - bdf_a2(jc)*rho(jc-2)
		e_hist   = bdf_a1(jc)*specific_internal_energy(jc-1, T_out(jc-1),   &
		                                               mmw(jc-1))           &
		         - bdf_a2(jc)*specific_internal_energy(jc-2, T_out(jc-2),   &
		                                               mmw(jc-2))
	else
		rho_hist = rho(jc-1)
		e_hist   = specific_internal_energy(jc-1, T_out(jc-1), mmw(jc-1))
	endif
	teq_cell%nhi    = nhi(jc);    teq_cell%nhii   = nhii(jc)
	teq_cell%nheiS  = nheiS(jc);  teq_cell%nheiTR = nheiTR(jc)
	teq_cell%y_HI   = y_HI_pp(jc); teq_cell%y_gnd = y_gnd_pp(jc)
	teq_cell%y_HeII = y_HeII_pp(jc)
	teq_cell%nheii  = nheii(jc);  teq_cell%nheiii = nheiii(jc)
	teq_cell%mup    = mmw(jc);    teq_cell%mum    = mmw(jc-1)
	teq_cell%rhov   = rho(jc)*v(jc)
	teq_cell%coeff  = mmw(jc-1)*v(jc)*(rho(jc) - rho_hist)
	teq_cell%dr_step = dr_step
	teq_cell%heaold = theat(jc)
	teq_cell%x_h2   = h2_particle_fraction(jc)
	teq_cell%e_hist = e_hist
	teq_cell%div_rhov = div_rhov_state(jc)
	call set_transport_terms_of_cell(jc)
	pp_nm_cell(:) = nm_w(jc,:)
	pp_beta_fs(:) = beta_fs_pp(jc,:)
	pp_nbar_fs(:) = nbar_fs_pp(jc,:)
	xv(1) = T_out(jc)
	iflag_r = 1
	call T_equation(1, xv, fv, iflag_r, paramsT)
	energy_residual_of_cell = fv(1)
	end function energy_residual_of_cell


	subroutine energy_column_newton()
	! Solve corrected-cell temperatures and cells without a marching root
	! together; other cells are fixed. Each energy residual depends on two
	! temperatures below, its own, and one above through conduction, so the
	! finite-difference Jacobian has two sub- and one super-diagonal.
	! Recompute conductivity and enthalpy fluxes at every iterate, limiting
	! each step to 30 percent of the current temperature.
	! The iteration stops at relative step 1e-10 or 100 iterations. Accept
	! only after evaluating all final energy rows at their final profile:
	! |R| <= atol_E + cert_tol_energy*S_E, S_E=max|grouped physical term|.
	! atol_E = 10*eps*(sum|grouped terms| + |rho*v*e_j| +
	! |rho*v*e_hist| + |c_lo*T_below| + |c_di*T| + |c_up*T_above|).
	! The floor covers cancellation in advection and conduction before
	! grouping. Reject LU failure, nonfinite terms, nonpositive T, or a
	! residual above this tolerance; report the first or worst unknown.
	! The conduction-coupled root has no 0.5-2*T_in band. Cells that the
	! one-cell marching solve could not root also enter the coupled solve;
	! EXHALE_ADV_COLUMN_EXCLUDE_UNROOTED=1 fixes them for measurement only.
	! A rejected column returns the marching profile and labels every
	! corrected unknown adv_failed. The measured rounding floor, temperature
	! range and unrooted cells are in
	! ../EXHALE_v1.00/md/Update_EXHALE_stage3.md section 55.
	! The returned profile is then judged on its RESIDUAL. With the transport coefficients and
	! the enthalpy divergences formed again at the final profile (the state
	! the rows describe), the energy row of every unknown is evaluated in
	! erg cm^-3 s^-1, R = q0 fvec/(mup mum dr_step), together with its seven
	! physically grouped terms (energy_row_grouped_terms: the advective
	! energy divergence div(rho e v) as one term, the p div(v) work, the
	! heating, the cooling, the conduction divergence and the two enthalpy
	! divergences).
	integer, parameter :: kl = 2, ku = 1, ldab = 2*kl + ku + 1
	integer, parameter :: it_max = 100
	integer, allocatable :: cells(:), pos(:), ipiv_c(:)
	real*8,  allocatable :: ab_c(:,:), g0(:), dT(:)
	integer :: nc, i, m, jr, it, info_c
	real*8  :: h, gp, lam_c, dmax, step_rel
	! The energy row of the worst unknown (energy_row_terms_of_cell), at the
	! column's temperature and with that cell at the run's own, printed
	! when the residual test rejects the column.
	real*8  :: row_terms_rejected(n_row_terms), row_terms_run_T(n_row_terms)
	real*8  :: T_keep
	! The residual test of the final rows: the terms of one row, its
	! grouped terms, S_E, sum |terms|, |R| and its tolerance; the measure
	! of the rounding floor, |R|/(eps sum|terms|), of every unknown (with
	! and without the two halves of the advective difference in the sum).
	real*8  :: row_t(n_row_terms), row_g(n_row_groups)
	real*8  :: S_E, sum_abs, R_abs, tol_cell, q_cell, q_worst
	real*8, allocatable :: floor_units(:), floor_units_raw(:)
	integer :: i_worst
	! Newton steps taken after the step test was first met.
	integer :: n_floor_steps
	nc = 0
	allocate(pos(1-Ng:N+Ng));  pos = 0
	! The unknowns: the cells the marching sweep corrected and, unless
	! adv_column_exclude_unrooted, the cells whose marching energy step
	! found no admissible root (adv_failed at this point; see the
	! statement at the top of the routine), counted apart.
	n_col_unrooted = 0
	do jr = 3-Ng, N+Ng
		if (adv_T_status(jr) .eq. adv_corrected .or.                     &
		    (.not. adv_column_exclude_unrooted .and.                      &
		     adv_T_status(jr) .eq. adv_failed)) then
			nc = nc + 1;  pos(jr) = nc
			if (adv_T_status(jr) .eq. adv_failed .and. jr .ge. 1 .and.    &
			    jr .le. N) n_col_unrooted = n_col_unrooted + 1
		endif
	enddo
	n_col_newton_it = 0
	col_newton_step = 0.0d0
	col_outcome     = adv_col_not_run
	col_reason      = adv_col_reason_none
	col_cell        = 0
	! The census counts physical cells only: the upper ghost rows the sweep
	! reaches are unknowns of the solve but not cells of the column.
	col_unknowns    = count(pos(1:N) .gt. 0)
	if (nc .eq. 0) then
		deallocate(pos);  return
	endif
	allocate(cells(nc), ipiv_c(nc), ab_c(ldab,nc), g0(nc), dT(nc))
	do jr = 3-Ng, N+Ng
		if (pos(jr) .gt. 0) cells(pos(jr)) = jr
	enddo
	pp_T_march = T_out
	pp_newton_mode = .true.
	col_bad = .false.
	n_floor_steps = 0
	do it = 1, it_max
		call transport_terms_of_profile(T_out)
		do i = 1, nc
			g0(i) = energy_residual_of_cell(cells(i))
		enddo
		do i = 1, nc
			if (.not. ieee_is_finite(g0(i))) then
				call reject_column(adv_col_reason_nonfinite, cells(i));  exit
			endif
		enddo
		if (col_bad) exit
		ab_c = 0.0d0
		do i = 1, nc
			m = cells(i)
			h = 1.0d-7*max(abs(T_out(m)), 1.0d-30)
			T_out(m) = T_out(m) + h
			do jr = m-1, m+2
				if (jr .lt. 3-Ng .or. jr .gt. N+Ng) cycle
				if (pos(jr) .eq. 0) cycle
				gp = energy_residual_of_cell(jr)
				! AB(kl+ku+1+row-col, col) of the LAPACK band storage
				ab_c(kl+ku+1+pos(jr)-i, i) = (gp - g0(pos(jr)))/h
			enddo
			T_out(m) = T_out(m) - h
		enddo
		dT = -g0
		call dgbsv(nc, kl, ku, 1, ab_c, ldab, ipiv_c, dT, nc, info_c)
		if (info_c .ne. 0) then
			! info_c > 0 is the zero pivot U(info_c,info_c); info_c < 0 an
			! illegal argument, which names no cell.
			if (info_c .gt. 0) then
				call reject_column(adv_col_reason_linear_solve, cells(info_c))
			else
				call reject_column(adv_col_reason_linear_solve, 0)
			endif
			exit
		endif
		do i = 1, nc
			if (.not. ieee_is_finite(dT(i))) then
				call reject_column(adv_col_reason_nonfinite, cells(i));  exit
			endif
		enddo
		if (col_bad) exit
		lam_c = 1.0d0
		do i = 1, nc
			m = cells(i)
			if (abs(dT(i)) .gt. 0.3d0*T_out(m))                            &
				lam_c = min(lam_c, 0.3d0*T_out(m)/abs(dT(i)))
		enddo
		dmax = 0.0d0
		do i = 1, nc
			m = cells(i)
			T_out(m) = T_out(m) + lam_c*dT(i)
			step_rel = abs(lam_c*dT(i))/max(T_out(m), 1.0d-30)
			if (step_rel .gt. dmax) then
				dmax = step_rel
			endif
		enddo
		n_col_newton_it = it
		col_newton_step = dmax
		if (dmax .lt. 1.0d-10) then
			! The measurement-only continuation (adv_column_floor_it).
			n_floor_steps = n_floor_steps + 1
			if (n_floor_steps .gt. adv_column_floor_it) exit
		endif
	enddo
	if (.not. col_bad) then
		do i = 1, nc
			if (.not. ieee_is_finite(T_out(cells(i)))) then
				call reject_column(adv_col_reason_nonfinite, cells(i))
				exit
			endif
			if (.not. (T_out(cells(i)) .gt. 0.0d0)) then
				call reject_column(adv_col_reason_non_positive_T, cells(i))
				exit
			endif
		enddo
	endif
	! THE RESIDUAL OF THE FINAL ROWS (the statement at the top of the
	! routine), with the transport terms of the final profile.
	if (.not. col_bad) then
		call transport_terms_of_profile(T_out)
		allocate(floor_units(nc), floor_units_raw(nc))
		q_worst = -1.0d0
		i_worst = 1
		do i = 1, nc
			call energy_row_terms_of_cell(cells(i), row_t)
			call energy_row_grouped_terms(row_t, row_g, S_E, sum_abs)
			R_abs = abs(row_t(11))
			if (.not. (ieee_is_finite(R_abs) .and. ieee_is_finite(sum_abs))) then
				call reject_column(adv_col_reason_nonfinite, cells(i))
				exit
			endif
			tol_cell = 10.0d0*epsilon(1.0d0)*(sum_abs + abs(row_t(12))       &
			           + abs(row_t(13)) + row_t(14)) + cert_tol_energy*S_E
			q_cell   = R_abs/max(tol_cell, tiny(1.0d0))
			floor_units(i)     = R_abs/max(epsilon(1.0d0)*sum_abs, tiny(1.0d0))
			floor_units_raw(i) = R_abs/max(epsilon(1.0d0)*(sum_abs           &
			                     + abs(row_t(12)) + abs(row_t(13)) + row_t(14)),  &
			                     tiny(1.0d0))
			if (q_cell .gt. q_worst) then
				q_worst = q_cell;  i_worst = i
				col_worst_cell    = cells(i)
				col_worst_R       = R_abs
				col_worst_rel     = R_abs/max(S_E, tiny(1.0d0))
				col_worst_tol_rel = tol_cell/max(S_E, tiny(1.0d0))
			endif
		enddo
		if (.not. col_bad) then
			call write_rounding_floor_measure(floor_units, floor_units_raw)
			if (q_worst .gt. 1.0d0)                                         &
				call reject_column(adv_col_reason_residual_above_tolerance,   &
				                   cells(i_worst))
		endif
		deallocate(floor_units, floor_units_raw)
	endif
	! THE ENERGY ROW OF THE CELL THE RESIDUAL TEST REJECTS, term by term, at
	! the temperature the column solve returned and with that cell at the
	! run's own temperature (the rest of the profile at the solve's), so a
	! reader can see which terms fail to balance. Printed only.
	if (col_bad .and. col_reason .eq. adv_col_reason_residual_above_tolerance &
	    .and. col_cell .ge. 1 .and. col_cell .le. N+Ng) then
		call transport_terms_of_profile(T_out)
		call energy_row_terms_of_cell(col_cell, row_terms_rejected)
		T_keep = T_out(col_cell)
		T_out(col_cell) = T_in(col_cell)
		call transport_terms_of_profile(T_out)
		call energy_row_terms_of_cell(col_cell, row_terms_run_T)
		T_out(col_cell) = T_keep
		call write_energy_row_terms(col_cell, row_terms_rejected,        &
		                            row_terms_run_T)
	endif
	! Every unknown of a column that converges is a corrected cell, the
	! cells whose marching step found no root included: the column
	! supplied their root.
	if (.not. col_bad) then
		do i = 1, nc
			adv_T_status(cells(i)) = adv_corrected
		enddo
	endif
	if (col_bad) then
		T_out = pp_T_march
		do i = 1, nc
			adv_T_status(cells(i)) = adv_failed
		enddo
		col_outcome = adv_col_rejected
		n_col_passes_rejected = n_col_passes_rejected + 1
	else
		col_outcome = adv_col_converged_on_residual
	endif
	pp_newton_mode = .false.
	write(*,'(A,I0,A,I0,A,I0,A,I0,A,ES10.3,A,L1,A,A,A,A,A,I0,A)')          &
	     ' (post_process_adv) transport terms: column solve of ',          &
	     col_unknowns, ' cells (', n_col_unrooted, ' with no marching'//   &
	     ' root; and ', nc - col_unknowns, ' ghost rows), ',              &
	     n_col_newton_it, ' iterations, last relative step ',              &
	     col_newton_step, '; marching profile kept: ', col_bad, ' (',          &
	     trim(adv_column_outcome_name(col_outcome)), ', reason ',          &
	     trim(adv_column_reason_name(col_reason)), ', cell ', col_cell, ')'
	if (col_worst_cell .ne. 0)                                             &
		write(*,'(A,I0,A,F9.5,A,ES10.3,A,ES10.3,A,ES10.3,A)')              &
		     '   worst energy row: cell ', col_worst_cell, ' (r = ',        &
		     r(col_worst_cell), ' Rp), |R| = ', col_worst_R,                &
		     ' erg cm^-3 s^-1, |R|/S_E = ', col_worst_rel,                  &
		     ' against the tolerance ', col_worst_tol_rel,                  &
		     ' (atol_E/S_E + cert_tol_energy)'
	if (col_bad .and. col_cell .ge. 1-Ng .and. col_cell .le. N+Ng)         &
		write(*,'(A,I0,A,F9.5,A,ES12.5,A,ES12.5,A)') '   cell ', col_cell,   &
		     ' at r = ', r(col_cell), ' Rp: T = ', col_T_trigger*T0,       &
		     ' K at the rejection, the run''s own T = ', T_in(col_cell)*T0, ' K'
	deallocate(pos, cells, ipiv_c, ab_c, g0, dT)
	end subroutine energy_column_newton

	subroutine energy_row_terms_of_cell(jc, terms)
	! THE ENERGY ROW OF CELL jc TERM BY TERM [erg cm^-3 s^-1], at the
	! current profile T_out and the transport coefficients of the last
	! transport_terms_of_profile. The residual of T_equation is the row
	! multiplied by mup*mum*dr_step in code units, so every term below is
	! q0 * (its part of fvec)/(mup*mum*dr_step), the sign of each as it
	! stands in the equation
	!     rho v de/dr - p v dln(rho)/dr + h div(rho v)
	!       = heating + conduction + enthalpy divergences - cooling .
	!   terms(1)  T [K]
	!   terms(2)  rho v de/dr, the advective term
	!   terms(3)  -p v dln(rho)/dr, the compression term of the residual
	!   terms(4)  h div(rho v), the enthalpy flux of the mass divergence
	!   terms(5)  p div(v) = terms(3) + (p/rho) div(rho v)
	!   terms(6)  photoheating
	!   terms(7)  radiative cooling
	!   terms(8)  heat conduction divergence, entered as heating
	!   terms(9)  enthalpy flux divergence of the element fluxes (heating)
	!   terms(10) enthalpy flux divergence of the carrier fluxes (heating)
	!   terms(11) the residual of the row, left side minus right side
	!   terms(12) rho v e_j/dr_step, the first half of the advective difference
	!   terms(13) rho v e_hist/dr_step, its second half; terms(2)=12-13.
	!             Both already include division by the step width and q0.
	!   terms(14) |c_lo T_below| + |c_di T| + |c_up T_above|, the magnitudes
	!             of the three parts the conduction divergence is the sum of
	!             (the face fluxes whose difference sets its rounding)
	! The cooling is the remainder of the residual T_equation returns once
	! the other terms are taken out, which is the cooling that equation
	! holds; its rounding is that of the largest term.
	integer, intent(in) :: jc
	real*8,  intent(out) :: terms(n_row_terms)
	real*8 :: fv, Mw, Tc, e_c, T_below, T_above, cond
	fv = energy_residual_of_cell(jc)
	Mw = teq_cell%mup*teq_cell%mum*teq_cell%dr_step
	Tc = T_out(jc)
	e_c = internal_energy_of_mixture(teq_cell%x_h2, Tc)
	if (jc .eq. 1) then
		T_below = pp_T_bath
	else
		T_below = T_out(jc-1)
	endif
	T_above = 0.0d0
	if (jc .lt. N+Ng) T_above = T_out(jc+1)
	cond = 0.0d0
	terms(14) = 0.0d0
	if (jc .ge. 1 .and. jc .le. N) then
		cond = pp_cond_lo(jc)*T_below + pp_cond_up(jc)*T_above             &
		     + pp_cond_di(jc)*Tc
		terms(14) = abs(pp_cond_lo(jc)*T_below) + abs(pp_cond_up(jc)*T_above) &
		          + abs(pp_cond_di(jc)*Tc)
	endif
	terms(1)  = Tc*T0
	terms(2)  = (teq_cell%mum*teq_cell%rhov*e_c                          &
	           - teq_cell%mup*teq_cell%mum*teq_cell%rhov*teq_cell%e_hist)/Mw
	terms(3)  = -teq_cell%coeff*Tc/Mw
	terms(4)  = teq_cell%mum*teq_cell%dr_step*teq_cell%div_rhov*(e_c + Tc)/Mw
	terms(5)  = terms(3) + teq_cell%div_rhov*Tc/teq_cell%mup
	terms(6)  = teq_cell%heaold
	terms(8)  = cond
	terms(9)  = pp_heat_elem(jc)
	terms(10) = pp_heat_carr(jc)
	terms(11) = fv/Mw
	terms(7)  = terms(11) - terms(2) - terms(3) - terms(4)                &
	          + terms(6) + terms(8) + terms(9) + terms(10)
	terms(12) = teq_cell%mum*teq_cell%rhov*e_c/Mw
	terms(13) = teq_cell%mup*teq_cell%mum*teq_cell%rhov*teq_cell%e_hist/Mw
	terms(2:14) = terms(2:14)*q0
	end subroutine energy_row_terms_of_cell

	subroutine energy_row_grouped_terms(terms, grouped, S_E, sum_abs)
	! THE PHYSICALLY GROUPED TERMS OF ONE ENERGY ROW (terms of
	! energy_row_terms_of_cell), whose signed sum is the residual:
	!     div(rho e v) + p div(v) - heating + cooling - conduction
	!       - element enthalpy - carrier enthalpy = R ,
	! with div(rho e v) = rho v de/dr + e div(rho v) the advective energy
	! divergence as ONE term (the three terms of the residual regrouped:
	! h div(rho v) = e div(rho v) + (p/rho) div(rho v), and p div(v) =
	! -p v dln(rho)/dr + (p/rho) div(rho v)). S_E is the largest magnitude
	! among them and sum_abs the sum of their magnitudes.
	real*8, intent(in)  :: terms(n_row_terms)
	real*8, intent(out) :: grouped(n_row_groups), S_E, sum_abs
	grouped(1) = terms(2) + terms(3) + terms(4) - terms(5)
	grouped(2) = terms(5)
	grouped(3) = terms(6)
	grouped(4) = terms(7)
	grouped(5) = terms(8)
	grouped(6) = terms(9)
	grouped(7) = terms(10)
	S_E     = maxval(abs(grouped))
	sum_abs = sum(abs(grouped))
	end subroutine energy_row_grouped_terms

	subroutine write_rounding_floor_measure(fu, fu_raw)
	! THE ROUNDING FLOOR OF THE FINAL ROWS, in units of eps sum|terms|: the
	! measurement that anchors the margin 10 of atol_E (energy_column_newton).
	! Printed for every column: the largest, the median and how many of
	! the unknowns exceed 1, 3 and 10; and the largest and the median when
	! the parts the row is assembled from are added to the sum: the two
	! halves of the advective difference and the three parts of the
	! conduction divergence (the floor of a row whose terms are small
	! differences of large ones).
	real*8, intent(in) :: fu(:), fu_raw(:)
	real*8, allocatable :: srt(:), srt_raw(:)
	integer :: n
	n = size(fu)
	allocate(srt(n), srt_raw(n))
	srt = fu;  srt_raw = fu_raw
	call sort_ascending(srt)
	call sort_ascending(srt_raw)
	write(*,'(A,I0,A,ES10.3,A,ES10.3,A,3(I0,A))')                          &
	     '   rounding floor |R|/(eps sum|terms|) over the ', n,             &
	     ' unknowns: max ', srt(n), ', median ', srt((n+1)/2),              &
	     '; above 1 / 3 / 10: ', count(fu .gt. 1.0d0), ' / ',              &
	     count(fu .gt. 3.0d0), ' / ', count(fu .gt. 1.0d1), ''
	write(*,'(A,ES10.3,A,ES10.3,A,3(I0,A))')                              &
	     '   the same with the advective halves and the conduction parts'//  &
	     ' in the sum: max ', srt_raw(n), ', median ', srt_raw((n+1)/2),    &
	     '; above 1 / 3 / 10: ', count(fu_raw .gt. 1.0d0), ' / ',          &
	     count(fu_raw .gt. 3.0d0), ' / ', count(fu_raw .gt. 1.0d1), ''
	deallocate(srt, srt_raw)
	end subroutine write_rounding_floor_measure

	subroutine sort_ascending(x)
	! Insertion sort in place (the columns hold a few hundred unknowns).
	real*8, intent(inout) :: x(:)
	real*8  :: v
	integer :: a, b
	do a = 2, size(x)
		v = x(a);  b = a - 1
		do while (b .ge. 1)
			if (x(b) .le. v) exit
			x(b+1) = x(b);  b = b - 1
		enddo
		x(b+1) = v
	enddo
	end subroutine sort_ascending

	subroutine write_energy_row_terms(jc, t_rej, t_run)
	! The two sets of terms of energy_row_terms_of_cell, side by side, and
	! each residual against S_E, the largest of its grouped terms
	! (energy_row_grouped_terms).
	integer, intent(in) :: jc
	real*8,  intent(in) :: t_rej(n_row_terms), t_run(n_row_terms)
	character(len=44), parameter :: lbl(2:11) = [character(len=44) ::    &
	     'rho v de/dr (advective)', '-p v dln(rho)/dr (compression)',     &
	     'h div(rho v) (mass divergence)', 'p div(v) (= compression + p/rho div(rho v))', &
	     'photoheating', 'radiative cooling, net (enters with -)',       &
	     'conduction divergence (as heating)',                            &
	     'element enthalpy divergence (as heating)',                      &
	     'carrier enthalpy divergence (as heating)', 'residual of the row']
	integer :: it
	real*8  :: g_rej(n_row_groups), g_run(n_row_groups)
	real*8  :: S_rej, S_run, sa_rej, sa_run
	write(*,'(A,I0,A,F9.5,A)') '   energy row of cell ', jc, ' (r = ',     &
	     r(jc), ' Rp) [erg cm^-3 s^-1], R = q0 fvec/(mup mum dr_step):'
	write(*,'(5X,A44,2(1X,A13))') ' ', 'column T', 'run''s own T'
	write(*,'(5X,A44,2(1X,ES13.5))') 'T [K]', t_rej(1), t_run(1)
	do it = 2, 11
		write(*,'(5X,A44,2(1X,ES13.5))') lbl(it), t_rej(it), t_run(it)
	enddo
	call energy_row_grouped_terms(t_rej, g_rej, S_rej, sa_rej)
	call energy_row_grouped_terms(t_run, g_run, S_run, sa_run)
	write(*,'(5X,A44,2(1X,ES13.5))') 'div(rho e v) (advective, one term)',   &
	     g_rej(1), g_run(1)
	write(*,'(5X,A44,2(1X,ES13.5))') '|residual| / S_E (largest grouped term)', &
	     abs(t_rej(11))/max(S_rej, 1.0d-300),                             &
	     abs(t_run(11))/max(S_run, 1.0d-300)
	end subroutine write_energy_row_terms

	subroutine reject_column(reason, cell)
	! The first test the column solve fails names the reason and the cell
	! (energy_column_newton); the solve reads col_bad.
	integer, intent(in) :: reason, cell
	col_bad    = .true.
	col_reason = reason
	col_cell   = cell
	! The temperature of that cell at the rejection, for the log.
	col_T_trigger = 0.0d0
	if (cell .ge. 1-Ng .and. cell .le. N+Ng) col_T_trigger = T_out(cell)
	end subroutine reject_column

	subroutine read_column_energy_solve_controls()
	! The two measurement keys of the column energy solve
	! (adv_column_floor_it, adv_column_exclude_unrooted). A value of the
	! first that cannot be read stops the run, since a measurement that
	! silently did not take place is worse than none.
	character(len=64) :: env_value
	integer :: ios_key
	adv_column_controls_read = .true.
	! EXHALE_ADV_COLUMN_FLOOR_IT: Overrides the number of extra column Newton iterations after the step
	! test. Unset retains the default; measurement only.
	call get_environment_variable('EXHALE_ADV_COLUMN_FLOOR_IT', env_value)
	if (len_trim(env_value) .gt. 0) then
		read(env_value, *, iostat=ios_key) adv_column_floor_it
		if (ios_key .ne. 0 .or. adv_column_floor_it .lt. 0) then
			write(*,'(A,A)') ' (post_process_adv) EXHALE_ADV_COLUMN_FLOOR_IT=', &
			     trim(env_value)//' is not a nonnegative integer.'
			error stop 1
		endif
		write(*,'(A,I0,A)') ' (post_process_adv) EXHALE_ADV_COLUMN_FLOOR_IT=', &
		     adv_column_floor_it, ': further Newton steps after the step'//   &
		     ' test (measurement only).'
	endif
	! EXHALE_ADV_COLUMN_EXCLUDE_UNROOTED: Set to 1 to hold column cells where the marching energy solve
	! found no root. Unset includes them; measurement only.
	call get_environment_variable('EXHALE_ADV_COLUMN_EXCLUDE_UNROOTED', env_value)
	adv_column_exclude_unrooted = (trim(env_value) .eq. '1')
	if (adv_column_exclude_unrooted)                                      &
		write(*,'(A)') ' (post_process_adv) EXHALE_ADV_COLUMN_EXCLUDE_UNROOTED=1:'//&
		     ' cells the marching energy step found no root for are fixed'//  &
		     ' values of the column (measurement only).'
	end subroutine read_column_energy_solve_controls

	end subroutine post_process_adv


	! Residual for the metal-only re-solve (pp_metals=2). Reuses the full
	! coupled H/He+metal balance (ion_system_HeH_metals, including the Huang
	! Table-4 charge exchange), then overwrites the three H/He rows with
	! identity equations that pin x(1..3) to the advection-corrected fractions
	! stored in pp_xHII_fix / pp_xHeII_fix / pp_xHeIII_fix. hybrd1 therefore
	! leaves H/He fixed and moves only the metal stages, while the electron
	! density and charge exchange inside the coupled residual still see the
	! correct (fixed) H/He densities. The cell-by-cell metal coefficients must be
	! loaded via set_metal_coeffs and the fractions pinned before each call.
	subroutine ion_system_metals_pp(N_in,x,fvec,iflag,params)
	integer :: N_in,iflag
	real*8  :: x(N_in),fvec(N_in)
	real*8  :: params(60)

	call ion_system_HeH_metals(N_in,x,fvec,iflag,params)

	! Pin the H/He fractions (rows 1-3) to the advection-corrected values.
	fvec(1) = x(1) - pp_xHII_fix
	fvec(2) = x(2) - pp_xHeII_fix
	fvec(3) = x(3) - pp_xHeIII_fix

	end subroutine ion_system_metals_pp

	end module post_processing
