	module post_processing
	! Subroutine to correct the output ionization profiles
	!	taking into account the ionization term
	!
	! Composition approximation (documented): the _adv reconstruction treats the
	! gas as H/He + trace metals and EXCLUDES the molecular species (H2, H2+,
	! H3+, HeH+). The advection correction receives only the H/He (nhi..nheiTR)
	! and metal (nm_in) densities; the molecular densities are not passed in and
	! are not re-solved here, so calc_ne, calc_ntot and
	! hydrogen_helium_nuclei_density below are called WITHOUT the nmol argument
	! (and, for the last, without the oxygen carriers either). In a molecular
	! run this omits the neutral-H2 particle count and the molecular-ion
	! electrons from the _adv n_tot/ne, and nh = nhi+nhii counts only the free
	! H nuclei. This is acceptable where the _adv
	! post-process is used (atomic/ionized escape flow); a molecular base needs
	! a molecular-aware post-process instead. Trace metals may now be solved
	! together with the molecular network, and pp_metal_mode carries the metal
	! stages here as usual -- but on that same molecule-free H/He background, so
	! _adv metal profiles inside the molecular layer inherit the approximation.

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
	                               row_terms_describe_state
	use certification,       only: cert_scale_floor
	! The species the reconstruction below omits, for the cell class it does
	! not model (see adv_unsupported).
	use ionization_equilibrium, only: nmol_eq, nox_eq

	implicit none

	! Validity limits of the advection correction (see post_process_adv, where
	! both are used and their derivation is written out).
	!
	! Da_local_equilibrium: above this Damkohler number the gas relaxes to the
	! local ionization equilibrium many times over while it crosses the cell,
	! so the equilibrium solution already solves the advection ODE. The number
	! is formed with the SLOWEST-relaxing species of the solved system (see
	! post_process_adv): the advection systems solve all H/He populations at
	! once, so pinning the cell to equilibrium is legitimate only when every
	! one of them is locally equilibrated.
	real*8, parameter :: Da_local_equilibrium = 1.0d2
	! xHII_adv_min: the residuals carry the NEUTRAL fraction and the ion
	! density is extracted as (1-x_HI)*n_h, so an equilibrium ion fraction
	! below this cannot be represented to better than ~1% by a solver whose
	! absolute resolution on x_HI is sqrt(eps) = 1.5e-8.
	real*8, parameter :: xHII_adv_min = 1.0d-6
	! THE MEASURE AND THE FRACTION A CORRECTED ROW IS ACCURATE TO. The
	! measure is the mass row of the state the post-process is handed,
	! divided by the largest term the row itself contains, |R_1(j)|/s_1(j),
	! which is the fractional change of the FACE mass flux across cell j
	! (steady_residual.f90, mass_flux_row_scale) and the same operator the
	! stationary certification measures that row by. The fraction it is
	! compared with is adv_conditional_tol (output_write), which is where
	! the first-order argument is written out: a row at or below it is a
	! conditional correction accurate to that fraction of itself, a row
	! above it keeps the run's own value. The certification's own verdict
	! on the WHOLE state, made with cert_tol_mass, is a separate statement
	! and travels in the file header beside this one.
	! Where the enthalpy flux of the mass-flux divergence is reported as
	! larger than both terms the steady balance keeps. It is a SENSITIVITY
	! statement about the answer and no longer a refusal (the block in
	! post_process_adv says why): unity is the point at which the term the
	! non-stationarity of the recorded flow puts into the energy equation
	! exceeds the terms the equation keeps, so above it the temperature the
	! equation returns is set by that term. It is not a stationarity
	! certificate: a ratio below one is reached with a mass divergence of
	! any size whenever the kept terms are large.
	real*8, parameter :: enthalpy_ratio_report_level = 1.0d0

	! Advection-corrected H/He ionized fractions for the current cell, pinned
	! while the metal re-solve (pp_metals=2) adjusts only the metal stages.
	! Set per cell before each ion_system_metals_pp / hybrd1 call.
	real*8, save :: pp_xHII_fix   = 0.0d0
	real*8, save :: pp_xHeII_fix  = 0.0d0
	real*8, save :: pp_xHeIII_fix = 0.0d0

	contains

	! Thermal Damkohler number of a cell: how many times over the local net
	! radiative rate could rewrite the internal energy of the gas while the
	! gas crosses the cell.
	!   t_cross  residence time of the gas in the cell [s]
	!   q_rad    magnitude of the local net radiative rate,
	!            |heating - cooling| [erg cm^-3 s^-1]
	!   u_th     internal energy density of the cell [erg cm^-3]
	!
	! Da >> 1: the gas reaches the local radiative balance while it is in the
	! cell, so the temperature there is set by heating = cooling and what the
	! gas carried in has been forgotten; the advected energy balance can only
	! add integration error, and the temperature to keep is the one the run's
	! own energy equation converged to, which carries every channel the run
	! solved (the post-process carries fewer; see the heating block below).
	! Da << 1: the temperature is carried by the flow and the advected balance
	! is what sets it.
	!
	! This is the energy counterpart of the Damkohler condition (ii) of the
	! ionization validity block below, and it is what makes the v -> 0 limit
	! of the temperature correction a statement about the flow rather than
	! about the sign of v: t_cross = dr/v grows without bound as v -> 0, so at
	! any nonzero radiative rate the correction switches ITSELF off before the
	! upwind difference it is built on loses its meaning. Unity is not a
	! tunable threshold: it is the statement that one of the two terms of the
	! equation is larger than the other.
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

	! Size of the enthalpy flux of the mass-flux divergence against the two
	! terms the steady internal-energy balance keeps, on one cell. Written
	! out, the equation is
	!
	!     rho v de/dr  -  p v dln(rho)/dr  +  h div(rho v)  =  heat - cool
	!
	! upwind-differenced, and with the common 1/dr divided out its three terms
	! are, in the variables of the residual,
	!
	!     q_adv  = |rho v (e_j - e_{j-1})|
	!     q_prs  = |w v (rho_j - rho_{j-1})|          w = p/rho
	!     q_enth = |h div(rho v) dr|                  h = e + w
	!
	! The third vanishes for a stationary mass flux rho v r^2 and for nothing
	! else, so this ratio is how much of the equation the non-stationarity of
	! the recorded flow carries, and above one the temperature the equation
	! returns is set by that term rather than by the balance.
	!
	! IT IS A SENSITIVITY DIAGNOSTIC AND NOT A STATIONARITY CERTIFICATE.
	! The ratio is the size of one term against two others, so it falls below
	! one at any size of the mass divergence once the kept terms are large
	! enough: a cell can carry a large continuity residual and a small ratio
	! at the same time. Stationarity of the cell is decided in
	! post_process_adv by the mass row of the state itself, the operator the
	! stationary certification uses; this ratio is reported next to it, with
	! enthalpy_ratio_report_level as the level at which the term dominates.
	!
	! The arguments are the cell (subscript up) and its upwind neighbor
	! (subscript lo): the two points at which the residual evaluates the
	! state.
	pure real*8 function enthalpy_flux_term_ratio(h_up, div_rhov, dr,      &
	                     rho_up, v_up, e_up, e_lo, w_up, rho_lo)
	real*8, intent(in) :: h_up      ! enthalpy per unit mass of the cell
	real*8, intent(in) :: div_rhov  ! divergence of the mass flux [1/length]
	real*8, intent(in) :: dr        ! width of the cell
	real*8, intent(in) :: rho_up    ! density of the cell
	real*8, intent(in) :: v_up      ! velocity of the cell
	real*8, intent(in) :: e_up,e_lo ! internal energy per unit mass, two points
	real*8, intent(in) :: w_up      ! p/rho of the cell
	real*8, intent(in) :: rho_lo    ! density of the upwind point
	real*8 :: q_enth,q_adv,q_prs

	q_adv  = abs(rho_up*v_up*(e_up - e_lo))
	q_prs  = abs(w_up*v_up*(rho_up - rho_lo))
	q_enth = abs(h_up*div_rhov*dr)
	! A cell in which all three terms vanish carries no equation at all, and
	! the floor makes that ratio zero rather than 0/0: nothing to refuse.
	enthalpy_flux_term_ratio = q_enth/max(q_adv, q_prs, 1.0d-99)

	end function enthalpy_flux_term_ratio

	subroutine post_process_adv(rho,v,p,T_in,heat,cool,eta,   &
                                  nhi_in,nhii_in,		    &
                                  nhei_in,nheii_in,nheiii_in,   &
                                  nheiTR_in, nm_in)


	real*8, dimension(1-Ng:N+Ng), intent(in) :: rho,v,p,T_in
   real*8, dimension(1-Ng:N+Ng), intent(in) :: heat,cool
   real*8, dimension(1-Ng:N+Ng), intent(in) :: eta
   real*8, dimension(1-Ng:N+Ng), intent(in) :: nhi_in,nhii_in
   real*8, dimension(1-Ng:N+Ng), intent(in) :: nhei_in,nheii_in,   &
      							  				        nheiii_in,nheiTR_in
   ! Converged equilibrium metal densities (dimensionless, n0 units), used
   ! by the metal-aware post-process modes (pp_metal_mode = 1 frozen, 2 re-solve).
   real*8, dimension(1-Ng:N+Ng,n_mion), intent(in) :: nm_in
	
	integer i,j,k
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
	 
	real*8, dimension(1-Ng:N+Ng) ::  T_K,p_out,T_out     ! Dimensional temperature
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
   ! He recombination radiation -> H ionization coupling scratch
   ! (use_he_rec_coupling; zero-effect when off).
   real*8, dimension(1-Ng:N+Ng) ::  rcheiiB_hrc,dP_HI_hrc,dheat_hrc
   ! H2 density seen by the He-recombination coupling. The _adv reconstruction
   ! is molecule-free (module-header composition note), so it is identically
   ! zero here and the coupling reduces to the H I / He I competition; the H2
   ! rate it returns is discarded for the same reason.
   real*8, dimension(1-Ng:N+Ng) ::  nh2_pp, dP_H2_hrc
   ! Metal share of the He recombination photons (item P34), added to P_m
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  dP_m_hrc
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
                                  meg_a1,meg_a2
   integer, dimension(n_melem) ::  meg_top
   integer :: i0,top,im
      
   real*8, dimension(1-Ng:N+Ng) ::  q13,q31a,q31b,Q31
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
      							  
	      
	      
   real*8 :: TT                              ! Temperature component
   real*8 :: PIR_1,PIR_15,PIR_2,PIR_TR       ! Photoionization rates
   real*8 :: deltal                          ! Optical depth
   real*8 :: dr	                           ! Grid spacing
   real*8 :: iup_1,ilo_1,        &           ! Photoheating integral variables
             iup_15,ilo_15,      &
             iup_2,ilo_2,        &
             iup_TR,ilo_TR,	   &
             iup_f,ilo_f 
   real*8 :: elo,eup                         ! Energy parameters
   real*8 :: tol,dpmpar                      ! Equilibrium system setup
   real*8 :: Hea_1 		          	         ! Heating rates
   real*8 :: brem,coex,coio,reco             ! Cooling rates
   real*8 :: iup_H,ilo_H                     ! Heating rate integral variables         
      
      
	! Substitution in the ODE solution
	real*8 :: As

	! Validity of the advection correction, cell by cell (filled by the block
	! just before the ionization loop, where the three conditions are stated).
	logical, dimension(1-Ng:N+Ng) :: adv_correction_valid
	real*8  :: t_cross        ! residence time of the gas in the cell [s]
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
	! energy equation (see the statement in the energy loop): its coefficient
	! div(rho v), and its size against the other two terms.
	real*8  :: x_h2_cell      ! H2 share of the particle count of the cell
	real*8  :: x_h2_upwind    ! ... of the upwind cell the energy comes from
	! e_up and e_lo are the two ENDS OF THE INTERFACE, the cell at r(j) and
	! the cell at r(j-1), for the term sizes of the sensitivity diagnostic;
	! e_upwind is the same specific energy of the upwind cell handed to the
	! residual, which names it for the direction the flow comes from.
	real*8  :: e_up,e_lo      ! internal energy per unit mass at the two ends
	real*8  :: e_upwind       ! that of the upwind cell, for the residual
	real*8  :: w_up,h_up      ! p/rho and the enthalpy per unit mass, this cell
	real*8  :: Fmass_up,Fmass_lo   ! mass flux rho v r^2 at the two points
	real*8  :: dV_cv          ! volume between the two points, d(r^3)/3
	real*8  :: div_rhov       ! divergence of the mass flux of the cell
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

	! Validity of the temperature and of the composition of each row, in the
	! five values of the two-field schema (output_write); written as the last
	! two columns of Hydro_ioniz_adv.txt and Ion_species_adv.txt. The values
	! are assigned where each decision is taken, so no reader of this routine
	! has to reconstruct which condition a row carries.
	integer, dimension(1-Ng:N+Ng) :: adv_T_status, adv_comp_status


   ! 60 entries to match the ion_system_HeH_metals / ion_system_metals_pp
   ! dummy length used by the pp_metals=2 re-solve (only 1-11 are consumed;
   ! the advection systems use 1-22).
   real*8, dimension(60) :: params
   real*8, dimension(40) :: paramsT
      
	real*8 :: rhop,rhom,vp,mum,mup,vm
	real*8 :: sys_sol_T(1), sys_x_T(1)
	! Largest term the energy equation holds in the current cell, the
	! scale its residual is read against (see the solve below).
	real*8 :: energy_residual_scale
   real*8 :: wa_T(8)
   logical :: brent_ok                       ! Task 1: Brent T-solve status
   ! The cell energy solve kept the run's own temperature: it did not
   ! converge, or its root was outside the band a temperature of this gas
   ! can occupy.
   logical :: T_solve_fell_back
      
      
   !----------------------------------------------------------!      
      
   ! Global parameters
      
   ! Numerical tolerance for system solution
   tol = sqrt(dpmpar(1))

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

   !----------------------------------!
	
	! Preliminary profiles extraction
	
	! Dimensional total number density profile and temperature
	T_K = T_in*T0
	         
   ! Initialize vectors
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

	!----------------------------------!

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

	!----------------------------------!

	! Non-converged cell solves, accumulated over every pass (see declaration)
	n_adv_noconv   = 0
	n_metal_noconv = 0
	n_T_noconv     = 0
	n_T_res_root   = 0

	! Iterate the post processing
	do k = 1,10	! Usually 10 gives a good convergence
	
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
	! to the hydrogen and helium nuclei").  Reused for both PH_heat calls
	! below (nhi/nheii/nheiii unchanged between them).  On this path the
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

	if (thereis_He) then
		call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm_w, xion,     &
					 P_HI,P_HeI,P_HeII,P_HeITR, P_m,        &
					 dum_v1,dum_v2)
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
		call HeITR_coeffs(T_K,rcheiTR,rcheiiB,A31,q13,q31a,q31b,Q31)
		! NOTE: rcheiiB is alpha1 from Oklopcic - being overwritten

	endif

	! He recombination radiation ionizing H I (Draine 2011; default off).
	! Correct the He II recombination coefficient and the H I photoionization
	! rate driving the advection ODE (the heating correction is applied later,
	! to theat). Mirrors ionization_equilibrium.
	if (use_he_rec_coupling .and. thereis_He) then
		nh2_pp = 0.0d0
		call he_rec_coupling(T_K, nhi, nh2_pp, nhei, nheii, nheiTR,       &
		                     ne, nm_w, A31, q31a, q31b,                    &
		                     rcheiiB_hrc, dP_HI_hrc, dP_H2_hrc,            &
		                     dP_m_hrc, dheat_hrc)
		rcheiiB = rcheiiB_hrc
		P_HI    = P_HI + dP_HI_hrc
		! The metal share of the same photons feeds the metal re-solve
		! (pp_metals=2) exactly as it feeds the equilibrium solve.
		if (thereis_metals) P_m = P_m + dP_m_hrc
	endif

   !----------------------------------!

	!---- How far from stationary is the flow of this cell? ----!
	!
	! Every equation this post-process solves is a STEADY equation integrated
	! along the recorded flow, and the flow it integrates along carries a
	! fractional change of the face mass flux across the cell. Without a mass
	! source a physically stationary atmosphere satisfies div(rho v) = 0, so
	! that change measures how far the state departs from the one the steady
	! equations describe, and the terms those equations drop are the ones the
	! mass divergence puts into them: each is the measure times a term they
	! keep. THE CORRECTION IS THEREFORE FIRST ORDER IN THE MEASURE, and a
	! cell is corrected exactly where the measure is small enough for that
	! error to be the stated one (adv_conditional_tol, output_write: one
	! percent of the correction itself). Where it is not, the equation
	! returns a temperature the gas does not have and the cell is REFUSED
	! (MEASURED on backup/regression/hydrostatic_column, whose mass flux
	! falls about 8 percent from one cell to the next: the exact equation
	! reads that as a compression at the same rate and reaches 1.7e7 K in a
	! 1124 K column, while the same equation with the enthalpy flux term
	! dropped gives the adiabat of the supplied density, 0.78 K). A refused
	! cell keeps the run's own temperature and the equilibrium composition at
	! that temperature, exactly as the thermal Damkohler condition below
	! keeps them where the local radiative balance rather than the flow sets
	! the state.
	!
	! THE OPERATOR IS THE PRODUCTION MASS ROW, AND ITS VALUE TRAVELS WITH
	! EVERY ROW. The mass balance of the run is a finite-volume one:
	! the numerical FACE fluxes of the Riemann solve, differenced over the
	! cell's own control volume. A center-to-center difference of rho v r^2 is
	! the same quantity in the continuum and a different discrete test, and a
	! finite-volume steady state can carry a small mass row while the
	! center-sampled rho v r^2 is not constant. So the state is handed to
	! assemble_residual and its mass row read back through
	! residual_row_scale, the two routines the stationary certification calls
	! (certification.f90, hydro_row_entry); the ratio |R_1|/s_1 is the
	! fractional change of the face mass flux across the cell. Nothing of
	! that arithmetic is restated here, and the ratio of every row is written
	! out as the adv_mass_row column, so a reader has the condition of each
	! row and not only its verdict. The certification's own verdict on the
	! whole state, taken with cert_tol_mass, is a different statement and is
	! reported beside this one in the file header.
	!
	! THE REFUSAL COVERS THE COMPOSITION TOO. The ionization correction of a
	! cell is the same kind of object: the steady advection-ionization ODE
	! along the same flow. A cell whose mass row is above the fraction has no
	! accurate advective ionization correction either, so the refusal is applied to the
	! ionization validity below as well and the cell keeps the equilibrium
	! composition at the run's temperature. Its own conditions (i) and (ii) do
	! not cover these cells: (i) tests the SIGN of v and (ii) the ratio of the
	! residence time to the ionization relaxation time, and both are satisfied
	! by a fast, well-directed outflow whose mass flux changes by a large
	! fraction of itself across every cell.
	!
	! CONTINUITY IS NECESSARY AND NOT SUFFICIENT. A stationary mass row says
	! nothing about the momentum and energy rows of the same state, and the
	! certification of the input state, reported in the header of the files
	! written at the end, is what carries those; a row of this file is a
	! conditional correction on the recorded wind either way, never a
	! self-consistent solution (write_adv_validity_header).
	!
	! FORMED ONCE, ON THE STATE HANDED IN. rho, v and r are the recorded
	! profile and the post-process never changes them, so the refusal is a
	! property of that state; it is evaluated on the first pass, at the run's
	! own temperature and the mean molecular weight of the composition the
	! equilibrium solve returned. Recomputing it per pass would let the
	! refused set move with the iterate, and every pass feeds its upwind
	! cascade from that set, so the converged _adv product would depend on the
	! iteration history.
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
		call assemble_residual(u_state, n_part_state, heat, cool, R_state)
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

		! THE SENSITIVITY OF THE ANSWER, next to the decision above. How much
		! of the energy equation the non-stationarity of the recorded flow
		! carries is the size of its enthalpy flux term against the two terms
		! the balance keeps (enthalpy_flux_term_ratio, this module), formed on
		! the two points the upwind difference is taken between. It is
		! reported and refuses nothing: the same ratio is reached with a mass
		! divergence of any size once the kept terms are large.
		enthalpy_flux_ratio_max  = 0.0d0
		r_enthalpy_flux_max      = 0.0d0
		n_enthalpy_term_dominant = 0
		do j = 2-Ng,N+Ng
			! The divergence of the mass flux over the control volume whose
			! two bounding points are the two the upwind difference is taken
			! between: (A_p F_p - A_m F_m)/dV with F = rho v, A = r^2,
			! dV = d(r^3)/3, an exact zero wherever the two points carry the
			! same rho v r^2. This is the coefficient the energy residual is
			! given, so the diagnostic and the equation see one number.
			Fmass_up = rho(j)  *v(j)  *r(j)**2
			Fmass_lo = rho(j-1)*v(j-1)*r(j-1)**2
			dV_cv    = (r(j)**3 - r(j-1)**3)/3.0d0
			div_rhov = (Fmass_up - Fmass_lo)/dV_cv
			! The three terms, with the common 1/dr divided out, in the
			! variables of the residual: e = E(x_H2,T)/mu with E the energy
			! per particle of the caloric EOS, p/rho = T/mu, h = e + p/rho.
			! Each end of the upwind energy difference is evaluated at the
			! composition of ITS OWN cell, as the residual does: the
			! rovibrational energy of H2 travels with the gas that holds the
			! molecules, so across a dissociation front the two ends store
			! different energy at the same temperature.
			dr        = r(j) - r(j-1)
			x_h2_cell   = h2_particle_fraction(j)
			x_h2_upwind = h2_particle_fraction(j-1)
			e_up = internal_energy_of_mixture(x_h2_cell,   T_in(j))         &
			       /mmw_in(j)
			e_lo = internal_energy_of_mixture(x_h2_upwind, T_in(j-1))       &
			       /mmw_in(j-1)
			w_up = T_in(j)/mmw_in(j)
			h_up = e_up + w_up
			enthalpy_flux_ratio = enthalpy_flux_term_ratio(h_up,        &
			     div_rhov, dr, rho(j), v(j), e_up, e_lo, w_up, rho(j-1))
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
	endif

   !----------------------------------!

	!---- Where is the advection correction valid? ----!
	!
	! The correction replaces the local ionization balance of a cell by the
	! steady advection-ionization ODE, integrated upwind across the cell. Three
	! conditions make that replacement carry no information; where any of them
	! holds the cell keeps the converged equilibrium ionization instead.
	!
	!  (i)   Inflow, v <= 0 on either face -- PHYSICAL. The upwind
	!        discretization takes the upstream state from the cell below, which
	!        is not the upstream cell when the gas moves inward (the breathing
	!        base). The residence time dr/v is then negative as well.
	!
	!  (ii)  Da = (dr/v)*nu_relax > Da_local_equilibrium -- PHYSICAL. The
	!        Damkohler number compares the time the gas spends in the cell with
	!        the relaxation time of the level populations. Da >> 1 means the
	!        populations relax to local equilibrium many times over while the
	!        gas crosses the cell, so the equilibrium solution IS the solution
	!        of the ODE and the correction can only add integration error.
	!
	!        nu_relax is the SLOWEST relaxation rate among the species the
	!        advection system actually solves, each one being the total rate at
	!        which its own population is destroyed and re-formed:
	!           H I/H II      : P_HI + (a_ion_HI + alpha_HII)*n_e
	!           He I/He II    : P_HeI
	!                           + (a_ion_HeI + alpha_HeII + alpha_HeI23S)*n_e
	!           He II/He III  : P_HeII + (a_ion_HeII + alpha_HeIII)*n_e
	!           He(2^3S)      : A31 + P_HeITR
	!                           + (q31a + q31b + a_ion_HeITR)*n_e + Q31*n_HI
	!        (the He(2^3S) row is exactly the loss side of fvec(4) of
	!        adv_implicit_HeH_TR, with the same rate coefficients from
	!        HeITR_coeffs / eval_cool -- no rate is redefined here.)
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
	!        neutral fraction x_HI and the ion density is extracted as
	!        (1-x_HI)*n_h, so the ion fraction inherits the solver's ABSOLUTE
	!        resolution on x_HI: xtol = sqrt(eps) = 1.5e-8, which is also the
	!        forward-difference step of the MINPACK Jacobian. Below 1e-6 the
	!        extracted ion fraction is worse than 1% relative, and in a
	!        shielded base with x_HII,eq ~ 1e-13 it is quantized at 1e-8 with
	!        an arbitrary sign -- a negative ion density that the upwind
	!        cascade then carries into the cells above.
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
	! for the upwind difference to read.

	adv_correction_valid = .true.
	adv_correction_valid(1-Ng) = .false.   ! inner boundary: never corrected
	adv_T_status    = adv_not_evaluated
	adv_comp_status = adv_not_evaluated
	n_adv_eq = 0
	n_stationarity_only = 0
	do j = 2-Ng,N+Ng
		if (v(j) <= 0.0d0 .or. v(j-1) <= 0.0d0) then
			adv_correction_valid(j) = .false.
		else
			t_cross  = (r(j) - r(j-1))*R0/(v(j-1)*v0)
			nu_relax = P_HI(j) + (a_ion_HI(j) + rchiiB(j))*ne(j)
			if (thereis_He) then
				! He II -> He I recombination: with the triplet on, rcheiiB is
				! the singlet channel alone (HeITR_coeffs overwrites it) and
				! rcheiTR is the triplet one, exactly as the He I row of the
				! residuals adds them.
				rec_HeII_tot = rcheiiB(j)
				if (thereis_HeITR) rec_HeII_tot = rec_HeII_tot + rcheiTR(j)
				nu_relax = min(nu_relax,                                    &
				     P_HeI(j)  + (a_ion_HeI(j)  + rec_HeII_tot)*ne(j),      &
				     P_HeII(j) + (a_ion_HeII(j) + rcheiiiB(j) )*ne(j))
				if (thereis_HeITR)                                          &
					nu_relax = min(nu_relax, A31 + P_HeITR(j)                &
					     + (q31a(j) + q31b(j) + a_ion_HeITR(j))*ne(j)        &
					     + Q31(j)*nhi(j))
			endif
			Da_slowest = t_cross*nu_relax
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
		if (.not. adv_correction_valid(j)) n_adv_eq = n_adv_eq + 1
		! The composition of a cell any of the four conditions refuses is the
		! equilibrium one the run converged to: the run's own value, kept.
		! A cell that is corrected here can still fail its own solve, which
		! the loop below records. The temperature field gets what conditions
		! (i) and (iv) say about it -- the upwind energy difference has no
		! upstream state under inflow either, and a non-stationary cell has
		! no steady energy equation -- and the energy loop adds the rest.
		if (adv_correction_valid(j)) then
			adv_comp_status(j) = adv_corrected
		else
			adv_comp_status(j) = adv_retained
		endif
	enddo

   !----------------------------------!

   ! Evolve species including the advection term in the
   ! 	ODE form
   ! Note: we are using point values here instead of 
   ! 	volume averages; they agree up to O(dr^2)
      
   ! The ionization fraction at the inner boundary are taken 
	!	from the input vectors (completely neutral atmosphere)
     

   ! Loop to solve the differential equation
   ! It is implicitly assumed that the velocity fields does 
   !	not change by including the advection term
      
      
	      
   if (.not.thereis_He) then
	      
		do j = 2-Ng,N+Ng
			! Outside its validity range the advection correction carries no
			! information; the cell keeps the converged equilibrium ionization.
			if (.not. adv_correction_valid(j)) then
				nhi(j)  = nhi_in(j)*n0
				nhii(j) = nhii_in(j)*n0
				cycle
			endif

			! Substitutions
			dr  = (r(j) - r(j-1))*R0
			As  = dr/(v(j-1)*v0)

			! Advection coeff.
			adv_cell%c1 = As
			adv_cell%xhi_old = nhi(j-1)/nh(j-1)
			adv_cell%nh = nh(j)
			adv_cell%P_HI = P_HI(j)
			adv_cell%rchiiB = rchiiB(j)
			adv_cell%a_ion_HI = a_ion_HI(j)
			! Metal electrons of this cell, per H nucleus, for the electron
			! density the recombination terms of the residual see.
			adv_cell%xe_metal = ne_metal(j)/max(nh(j),1.0d-30)

			! Initial guess of solution
			sys_x(1) = nhi(j)/nh(j)     
			
			! Call hybrd1 routine (from minpack)
			call hybrd1(adv_implicit_H,Neq_adv,sys_x,sys_sol,   &
						tol,info,wa,lwa_adv,params)
			
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
      	
		   ! Force condition of zero helium
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

			! Substitutions
			dr  = (r(j) - r(j-1))*R0
			As  = dr/(v(j-1)*v0)

			! Advection coeff.
			adv_cell%c1  = As
			adv_cell%xhi_old  = nhi(j-1)/nh(j-1)
			adv_cell%xheiS_old  = nheiS(j-1)/nhe(j-1)
			adv_cell%xheiii_old  = nheiii(j-1)/nhe(j-1)
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
			! B), read by he_h_cx_fvec_adv in the H/He adv systems. T-only, so
			! evaluate once per cell; the adv residual adds nothing when
			! he_h_charge_exchange is off (bit-identical).
			call he_h_cx_rates(T_K(j), adv_cell%kcx_He0_Hp,                &
			                           adv_cell%kcx_Hep_H0)
			! Effective He/H for the electron density inside the adv system:
			! the global HeH normally (byte-identical legacy), the local
			! (diffused) nhe/nh when He_diffusion is on.
			if (he_diffusion) then
				adv_cell%heh_loc = nhe(j)/max(nh(j),1.0d-30)
			else
				adv_cell%heh_loc = HeH
			endif
			! Metal electrons of this cell, per H nucleus, for the electron
			! density the recombination terms of the residual see.
			adv_cell%xe_metal = ne_metal(j)/max(nh(j),1.0d-30)

			! Add more if HeITR is present
			if (thereis_HeITR) then
				adv_cell%rcheiTR = rcheiTR(j)
				adv_cell%A31 = A31
				adv_cell%P_HeITR = P_HeITR(j)
				adv_cell%q13 = q13(j)
				adv_cell%q31a = q31a(j)
				adv_cell%q31b = q31b(j)
				adv_cell%Q31 = Q31(j)
				adv_cell%a_ion_HeITR = a_ion_HeITR(j)
				adv_cell%xheiTR_old = nheiTR(j-1)/nhe(j-1)
				! Effective He/H for the electron density (see non-TR block).
				if (he_diffusion) then
					adv_cell%heh_loc = nhe(j)/max(nh(j),1.0d-30)
				else
					adv_cell%heh_loc = HeH
				endif
			endif
			
			! Initial guess of solution
			sys_x(1) = nhi(j)/nh(j) 
			sys_x(2) = nheiS(j)/nhe(j)
			sys_x(3) = nheiii(j)/nhe(j)
			if (thereis_HeITR) sys_x(4) = nheiTR(j)/nhe(j) 
			
			! Call hybrd1 routine (from minpack)
			if (thereis_HeITR) then 
				call hybrd1(adv_implicit_HeH_TR,Neq_adv,sys_x,sys_sol,   &
							tol,info,wa,lwa_adv,params)
			else
				call hybrd1(adv_implicit_HeH,Neq_adv,sys_x,sys_sol,   &
							tol,info,wa,lwa_adv,params)
			endif
				
			! A cell whose advection system did not converge carries no
			! correction: the returned iterate satisfies neither the
			! advection balance it was asked to solve nor the equilibrium
			! balance it started from, so it is not a state of the gas. The
			! equilibrium solution of that same cell is, and it is what the
			! three validity conditions above already fall back to.
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
	      
      
   !---------------------------------------------------------
      
   !------- Fix stationarity of new pressure profile -------!
      
   !---- Update densities and temperature ----!
	
	! Number densities      
   nhei = nheiS + nheiTR
	call hydrogen_helium_nuclei_density(nhi,nhii,nhei,nheii,nheiii,nh,nhe)

   ! Total number density (incl. metal nuclei under eos_metals). Molecular
   ! species are excluded -- the post-process does not carry them (see the
   ! module-header composition note).
   call calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm_w)

   ! Free electron density (assuming overall neutrality; incl. metal
   ! electrons under eos_metals; molecular-ion electrons excluded, as above)
   call calc_ne(nhii,nheii,nheiii,ne,nm_w)

	!----------------------------------!

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

			! Charge-exchange rate coefficients for this cell temperature.
			call cx_set_cell(T_K(j))

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
					meg_g1(im) = P_m(j,i0+1)
					meg_b1(im) = aion_m_pp(j,i0+1)
					meg_a2(im) = rec_m_pp(j,i0+2)
				else
					meg_g1(im) = 0.0d0
					meg_b1(im) = 0.0d0
					meg_a2(im) = 0.0d0
				endif
			enddo
			call set_metal_coeffs(n_melem, meg_ntot, meg_g0, meg_g1,   &
			                      meg_b0, meg_b1, meg_a1, meg_a2, meg_top)

			! Pin the advection-corrected H/He fractions for the wrapper.
			pp_xHII_fix   = nhii(j)/nh(j)
			pp_xHeII_fix  = nheii(j)/max(nhe(j),1.0d-300)
			pp_xHeIII_fix = nheiii(j)/max(nhe(j),1.0d-300)

			! Named-field cell state for the residual: only n_h, n_he are
			! consumed once rows 1-3 are pinned (the H/He rates drop out). pp
			! runs serially on the master thread, so this master threadprivate
			! copy is the one read by ion_system_HeH_metals inside the wrapper.
			! params stays only the MINPACK transport argument (unread).
			params    = 0.0d0
			ieq_cell%nh  = nh(j)
			ieq_cell%nhe = nhe(j)

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

	!----------------------------------!

	!---- Heating of the advection-corrected composition ----!

	! The photoheating of ONE particle of each absorber, from the same
	! attenuated field the equilibrium pass used, and then the ONE heating
	! assembly (utils_ion_eq) contracted with the advection-corrected
	! densities. This is the same routine the ionization sweep and the
	! heating breakdown call, so the _adv energy solve balances the heating
	! the run's own energy equation deposits and cannot drift from it.
	!
	! The composition reconstructed here carries no molecular and no oxygen
	! carriers (see the header of this module), which is what the two
	! composition flags below say; the molecular and oxygen deposits are
	! therefore absent from the _adv heating, as are the molecular carriers
	! from its n_e and its n_tot.
	if (thereis_He) then
		call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm_w, xion,     &
		                 dum_v1,dum_v2,dum_v3,dum_v4, P_m,          &
		                 dum_v6,dum_v5,                             &
		                 heat_of_one_HI   = h1_HI_pp,               &
		                 heat_of_one_HeI  = h1_HeI_pp,              &
		                 heat_of_one_HeII = h1_HeII_pp,             &
		                 heat_of_one_HeTR = h1_HeTR_pp,             &
		                 heat_of_one_H2   = h1_H2_pp,               &
		                 heat_of_one_mion = h1_m_pp)
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
	         .false., .false., theat, heat_chan_pp)

	! Adimensionalize
	theat = theat/q0

	!----------------------------------!
	
	!---- Solve stationary energy equation ----!
	! This procedure uses the same velocity profile and 
	!	the ionization profile after the advection correction
	
	! Initialize temperature at ghost cells
	T_out = T_K/T0
	
	! Calculate mean molecular weight. nm_w adds the metal mass/nuclei under
	! the eos_metals policy, so the _adv temperature solve uses the same
	! composition as the main loop (ne above already carries the metal
	! electrons via calc_ne).
	call calc_mmw(nh,nhe,ne,mmw,nm_w)

	! Line transfer of the ground-term fine-structure lines, from the
	! incoming profile, so the cell-by-cell energy solve balances the same
	! metal cooling eval_cool reports.
	call fine_structure_line_transfer(T_K, nm_w, beta_fs_pp, nbar_fs_pp)

	! Count cell-by-cell temperature solves rejected as non-physical (metal modes).
	n_pp_reject = 0
	n_T_local   = 0
	Da_thermal_max   = 0.0d0
	r_Da_thermal_max = 0.0d0

	!---- What the energy correction solves ----!
	!
	! The equation solved cell by cell below (T_equation) is the steady
	! internal-energy equation of the supplied profile,
	!
	!     div(u v) + p div(v)  =  heating - cooling
	!   = rho v de/dr  -  p v dln(rho)/dr  +  h div(rho v)             (E)
	!
	! upwind-differenced, with u the internal energy density, e the internal
	! energy per unit mass and h = e + p/rho the enthalpy per unit mass. The
	! two forms are one equation: div(rho e v) = e div(rho v) + rho v de/dr,
	! and p div(v) = (p/rho) div(rho v) - p v dln(rho)/dr.
	!
	! The third term is the enthalpy flux carried by the divergence of the
	! mass flux, h div(rho v) = h (1/r^2) d(rho v r^2)/dr. It vanishes for a
	! stationary mass flux rho v r^2 = const and for nothing else, so a
	! correction built without it describes a converged wind and no other
	! state. It is handed to the residual through teq_cell%div_rhov, whose
	! two branches T_equation states the algebra of.
	!
	! WHICH DIVERGENCE, AND WHY. div(rho v) is formed with the operator the
	! mass row of the state uses, (A_p F_p - A_m F_m)/dV with F = rho v,
	! A = r^2 and dV = d(r^3)/3 (RK_rhs), over the control volume whose two
	! bounding points are the two points the upwind energy difference is taken
	! between, r(j-1) and r(j). Those are the only two points at which this
	! residual evaluates the state, so the mass flux entering the term is the
	! state's own at the same two points, and the term is an exact zero
	! wherever the two carry the same rho v r^2. A centered difference of the
	! neighboring cells would not have that property and would measure a
	! divergence at a point the energy difference never visits. The FACE mass
	! flux of the state is a different object and is not read here: its faces
	! are not the two points of this difference. It is what decides WHETHER
	! this cell has a steady equation at all -- the mass row of the state,
	! measured in the stationarity block above by the operator the stationary
	! certification uses -- and the two answer two questions: that block
	! how far the state departs from stationary in this cell, this term what
	! the equation carries at that departure.
	!
	! WHAT THE TERM IS WORTH. MEASURED 2026-09-08 on
	! backup/regression/hydrostatic_column, a 300-step mechanical column whose
	! mass flux rho v r^2 runs over seven orders of magnitude across the
	! domain: without the third term the temperature is the upwind recursion
	! of the first two, whose continuum form is the adiabat T proportional to
	! rho^(gamma-1), and it falls from 1138 K at the base to 0.78 K at
	! 1.396 R_p while the run's own stays between 1084 and 1182 K. With the
	! third term the same column rises instead, to 1.7e7 K at 1.40 R_p and
	! 1.6e9 K at 2.96 R_p, because its mass flux falls by about 8 percent
	! from one cell to the next and (E) reads that as a compression at the
	! same rate. Both numbers are solutions of an equation on a state that
	! does not satisfy continuity, and neither is a temperature the gas has.
	! The term is therefore not a correction to that state, it is the state.
	! Such a cell is refused by condition (iv) of the validity block above,
	! which asks the mass row of the state itself; the size of this term
	! against the other two is reported next to that verdict as the
	! sensitivity of the answer, and the cells the condition refuses keep the
	! run's own temperature below.
	!
	! The check on the restored term is the case with a solution in closed
	! form: with heating and cooling negligible (E) integrates to
	! w = p/rho proportional to rho^(gamma-1) F^(-gamma), F = rho v r^2, and
	! src/tests/adv_static_limit marches a column of constant density and
	! F = r against that solution (first order in the cell width, MEASURED).

	do j = 3-Ng,N+Ng ! Start from first computational cell

		! Substitutions
		rhop = rho(j)
		rhom = rho(j-1)
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
		! Inflow (condition (i) of the ionization validity block above): the
		! cell keeps the converged eq temperature. The advection-corrected
		! energy solve is upwind-differenced just like the ionization solve, so
		! it is invalid wherever the gas moves inward, and it would otherwise
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
		mum = mmw(j-1)
		mup = mmw(j)

		! Local radiative balance (thermal_damkohler_number, this module): the
		! gas that spends longer in the cell than the local net radiative rate
		! needs to rewrite its internal energy has forgotten what it carried
		! in, so its temperature there is the root of heating = cooling and
		! not of the advected balance. The temperature to keep is then the
		! run's own, which is the root of that same local balance with every
		! channel the run solved. The pair is formed at the state the pass
		! starts from, since the correction has to be decided before the
		! temperature is solved for.
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

		! The divergence of the mass flux of this cell, on the two points the
		! upwind energy difference is taken between (see the statement above):
		! the mass row's flux-difference operator over the control volume
		! bounded by r(j-1) and r(j). The same operator, on the state handed
		! in, forms the stationarity condition (iv) of the validity block.
		Fmass_up = rho(j)  *vp*r(j)**2
		Fmass_lo = rho(j-1)*vm*r(j-1)**2
		dV_cv    = (r(j)**3 - r(j-1)**3)/3.0d0
		div_rhov = (Fmass_up - Fmass_lo)/dV_cv

		! The caloric state of the two cells the upwind energy difference is
		! taken between. e_upwind is the SPECIFIC internal energy of the
		! upwind cell, E(x_H2,up, T_up)/mu_up, which is the quantity the flow
		! carries into this cell: it is evaluated at the upwind composition
		! and the upwind temperature, so a cell with no H2 below a molecular
		! neighbor still receives the rovibrational energy of that gas. The
		! upwind temperature is the one the cell below was left with by this
		! same loop, T_out(j-1), the value the residual differences against.
		x_h2_cell   = h2_particle_fraction(j)
		x_h2_upwind = h2_particle_fraction(j-1)
		e_upwind    = internal_energy_of_mixture(x_h2_upwind, T_out(j-1))   &
		              /mum

	 	!--- Solve equation for temperature implicitly ---!
		
		! Parameters
		teq_cell%nhi  = nhi(j)
	 	teq_cell%nhii  = nhii(j)
	 	teq_cell%nhei  = nhei(j)
	 	teq_cell%nheii  = nheii(j)
	 	teq_cell%nheiii  = nheiii(j)
	 	teq_cell%mup  = mmw(j)
	 	teq_cell%mum  = mmw(j-1)
	 	teq_cell%rhov  = rhop*vp
	 	teq_cell%coeff  = mum*vp*(rhop-rhom)
	 	teq_cell%dr = dr
	 	teq_cell%Told = T_out(j-1)
	 	teq_cell%heaold = theat(j)
	 	! Composition entries of the caloric EOS: the H2 share of the
	 	! particle-plus-electron count of this cell and of the upwind one, as
	 	! the equilibrium solve left them, and the specific internal energy
	 	! the flow carries in with it.
	 	teq_cell%x_h2    = x_h2_cell
	 	teq_cell%x_h2_up = x_h2_upwind
	 	teq_cell%e_up    = e_upwind
	 	! Coefficient of the enthalpy flux term of (E), formed above from the
	 	! state this pass was handed. It does not change while the root finder
	 	! varies the temperature.
	 	teq_cell%div_rhov = div_rhov
	 	! Metal densities for this cell [cgs] go through the equation_T module
	 	! array (the 27-ion vector does not fit params). pp_metal_on gates
	 	! whether T_equation adds the metal cooling/brem/n_e terms.
	 	pp_nm_cell(:)  = nm_w(j,:)
	 	pp_beta_fs(:)  = beta_fs_pp(j,:)
	 	pp_nbar_fs(:)  = nbar_fs_pp(j,:)

	 	! Initial guess of solution
		sys_x_T(1) = T_out(j)

	 	! Task 1: solve the scalar energy equation by bracketing the physical
	 	! (lowest) root + Brent when metal cooling is on (the default). The
	 	! metal-cooled residual is non-monotone and has a second, spurious *hot*
	 	! root that a Newton/Powell solve (hybrd1) could land on; bracketing from
	 	! below selects the physical root structurally. Fall back to the converged
	 	! eq T if no bracket is found. With "Brent solver: False" (use_brent_tsolve
	 	! = .false.) the legacy MINPACK solve + 2x-band reject is used instead.
	 	! Metals-off always keeps the original MINPACK solve (monotone residual,
	 	! byte-identical).
	 	T_solve_fell_back = .false.
	 	if (pp_metal_on .and. use_brent_tsolve) then
	 		call solve_T_brent(paramsT, T_in(j), sys_x_T(1), brent_ok)
	 		if (.not. brent_ok) then
	 			sys_x_T(1)  = T_in(j)
	 			n_pp_reject = n_pp_reject + 1
	 			T_solve_fell_back = .true.
	 		endif
	 	else
	 		! Legacy MINPACK solve.
	 		call hybrd1(T_equation,1,sys_x_T,sys_sol_T,   &
	 		            tol,info,wa_T,8,paramsT)
	 		! A non-converged solve leaves an iterate that balances neither
	 		! the advected energy equation nor the equilibrium one; the
	 		! converged equilibrium temperature is the state to keep.
	 		!
	 		! WHAT info MEANS AND WHAT THE ROOT IS. MINPACK's info states how
	 		! its ITERATION ended, not whether the iterate is a root: info = 4
	 		! and 5 are returned when the steps stop improving the residual,
	 		! which is what a stalled search and an ARRIVED one look like
	 		! alike. At the root the residual sits at the cancellation floor
	 		! of the terms it is assembled from and no step can lower it, so
	 		! a solve that arrives to full precision is reported exactly as
	 		! one that never got there. The root of a scalar equation is
	 		! defined by its residual, so the iterate is kept whenever that
	 		! residual is negligible against the largest term the equation
	 		! holds, and only an iterate that is not a root falls back.
	 		!
	 		! THE SCALE. The terms of (E) as T_equation assembles them: the
	 		! advected internal energy of this cell and of its upwind
	 		! neighbor, the compression work, and the photoheating. A sum of
	 		! terms of size s cannot be formed to better than a few machine
	 		! epsilons of s, so 1e2*epsilon(s) is the level at which the
	 		! equation is an identity in double precision; an iterate that is
	 		! not a root stands orders of magnitude above it. MEASURED on the
	 		! three cells of the LHS 1140 b 45 Rp wind that reach this branch:
	 		! |R|/s = 1.5e-17, 2.8e-17 and 5.2e-17, against a fallback that
	 		! discarded roots good to every digit and, because the correction
	 		! is an upwind recursion, restarted the profile above them
	 		! (docs/lhs1140b_stationary_L10_20260913.md).
	 		!
	 		! WHY IT MATTERS MORE THAN ONE CELL. The corrected temperature of
	 		! a cell is differenced against the corrected temperature of the
	 		! cell below, so a discarded root is not a local blemish: every
	 		! row above it is integrated from a different starting value.
	 		if (info /= 1) then
	 			energy_residual_scale =                                        &
	 			   max(abs(mum*teq_cell%rhov*sys_x_T(1)),                      &
	 			       abs(mup*teq_cell%rhov*teq_cell%Told),                   &
	 			       abs((gamma_ad - 1.0d0)*teq_cell%coeff*sys_x_T(1)),      &
	 			       abs((gamma_ad - 1.0d0)*mup*mum*dr*theat(j)))
	 			if (abs(sys_sol_T(1)) <=                                       &
	 			    1.0d2*epsilon(1.0d0)*energy_residual_scale) then
	 				n_T_res_root = n_T_res_root + 1
	 			else
	 				sys_x_T(1) = T_in(j)
	 				n_T_noconv = n_T_noconv + 1
	 				T_solve_fell_back = .true.
	 			endif
	 		endif
	 		! A non-positive root is not a temperature, whatever else is in the
	 		! gas, so this test is not conditional on the metals. It matters
	 		! because T_out feeds the NEXT post-process pass: eval_cool takes
	 		! sqrt(T/T0) in the Badnell recombination fit (rr_badnell), so a
	 		! negative T there is a NaN cooling rate in an ordinary build and an
	 		! abort under -ffpe-trap=invalid. Measured on the He/H = 1 molecular
	 		! case, which is metals-off and so had no test at all: cells 278-280
	 		! come back at -42, -640 and -2474 K on the second pass. Keep the
	 		! converged equilibrium temperature -- the same state the
	 		! non-converged branch above keeps, and for the same reason.
	 		if (.not. (sys_x_T(1) > 0.0d0)) then
	 			sys_x_T(1)  = T_in(j)
	 			n_pp_reject = n_pp_reject + 1
	 			T_solve_fell_back = .true.
	 		endif
	 		! With metal cooling, additionally reject an out-of-band root: the
	 		! metal-cooled residual is non-monotone and carries a second,
	 		! spurious HOT root that hybrd1 can land on.
	 		if (pp_metal_on) then
	 			if (sys_x_T(1) > 2.0d0*T_in(j)  .or.   &
	 			    sys_x_T(1) < 0.5d0*T_in(j)) then
	 				sys_x_T(1)  = T_in(j)
	 				n_pp_reject = n_pp_reject + 1
	 				T_solve_fell_back = .true.
	 			endif
	 		endif
	 	endif

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

	! Report how many cells kept the run's own temperature because the local
	! radiative balance, not the flow, sets the temperature there, and how
	! close to that condition the rest of the column came. The largest
	! Damkohler number is reported whether or not any cell crossed unity: it
	! is the margin of the condition on this state, and a column that stays
	! three orders of magnitude below it and one that sits at 0.95 are
	! different statements about the same verdict.
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

	! Report how many cells fell back to the eq temperature: a non-positive
	! root (any run) or, with metals on, one outside the 0.5-2x band.
	if (n_pp_reject > 0) then
		write(*,'(a,i0,a,i0,a)') ' (post_process_adv) T solve: ',            &
		   n_pp_reject, ' of ', N+2*Ng-2,                                     &
		   ' cells fell back to eq T (non-positive or out-of-band root).'
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

   ! ---------------------------- !
      
	!---- Update cooling rates ----!

	! Molecule-free electron sum and no H3+ cooling, as at the first eval_cool
	! call above; T_equation, which solved for this T_K, assembles the same
	! channels.
	call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm_w,            &
	  			   dum_v1,dum_v2,dum_v3, rec_m_pp,                   &
				   dum_v4,dum_v5,dum_v6, aion_m_pp,                      &
				   tcool, nheiTR=nheiTR)

	! Adimensionalize
	tcool = tcool/q0

	!----------------------------------!
 
   ! Adimensionalize ion densities before writing
   nhi_w    = nhi/n0
   nhii_w   = nhii/n0
   nhei_w   = nhei/n0
   nheii_w  = nheii/n0
   nheiii_w = nheiii/n0
   nheiTR_w = nheiTR/n0

   !----------------------------------!
      
   ! Write updated thermodynamic and ionization profiles. Metals are written
   ! in dimensionless (n0) units, matching the other species; zero in the
   ! metal-free mode, frozen/re-solved eq densities otherwise.
   nm_out = nm_w/n0
   ! The two status fields say, row by row, whether the temperature and the
   ! composition of that row are the steady correction or the run's own
   ! state, so a reader of the file (and of a spectrum built from it) can
   ! tell which rows carry which.
   call write_output(rho,v,p_out,T_out,theat,tcool,eta,                &
                     nhi_w,nhii_w,nhei_w,nheii_w,nheiii_w,              &
                     nheiTR_w,nm_out,'ad', adv_T_status, adv_comp_status,  &
                     adv_mass_row)

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

	! End of subroutine
	end subroutine post_process_adv

	! ----------------------------------------------------------------- !

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

	! End of module
	end module post_processing
