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

	! THE STEP OF THE RECURSION. Every equation this post-process integrates
	! along the flow is an ODE in r, dy/dr = F(y, r): the ionization balance
	! v dx/dr = R(x) of the H, HeH and HeH_TR advection systems, and the steady
	! internal-energy equation (the statement before the energy loop). Both are
	! stiff wherever the chemistry or the radiative balance relaxes within a
	! cell (a Damkohler number above one). Each is advanced cell by cell with
	! ONE rule, the backward differentiation formula of order two on the
	! nonuniform grid (variable-step BDF2):
	!
	!     y_j - a1 y_{j-1} + a2 y_{j-2} = g h_j F(y_j, r_j) ,
	!     h_j = r_j - r_{j-1},  w = h_j/h_{j-1},
	!     a1 = (1+w)^2/(1+2w),  a2 = w^2/(1+2w),  g = (1+w)/(1+2w) .
	!
	! Derivation: the left side over g h_j is dP/dr at r_j for P the quadratic
	! through (r_{j-2},y_{j-2}), (r_{j-1},y_{j-1}), (r_j,y_j), since
	! dP/dr(r_j) = y_j (2h_j + h_{j-1})/(h_j (h_j + h_{j-1}))
	!            - y_{j-1} (h_j + h_{j-1})/(h_j h_{j-1})
	!            + y_{j-2} h_j/(h_{j-1} (h_j + h_{j-1})) ;
	! multiplying by g h_j = h_j (h_j + h_{j-1})/(2h_j + h_{j-1}) gives the
	! weights above. This is the variable-step BDF2 formula in the form
	! D2 v^n = d0(r_n,0) dv^n + d1(r_n,0) dv^{n-1} of Li & Liao, "Stability
	! of variable-step BDF2 and BDF3 methods", arXiv:2201.00527 (2022),
	! eq. (1.4), which was read for this; the textbook source they cite is
	! Hairer, Norsett & Wanner, Solving Ordinary Differential Equations I,
	! 2nd ed., sect. III.5 (not re-read here). At w = 1 they are the constant-step BDF2, (3/2) y_j - 2 y_{j-1}
	! + (1/2) y_{j-2} = h F_j. The formula is exact for quadratics, so its
	! local error is O(h^3) and the recursion is second order in the cell
	! width. The same weights differentiate the known profiles the energy
	! equation holds (the density), so every r-derivative of a step is taken
	! at r_j by one rule, and every other factor of the step (the rates, the
	! heating and cooling, the velocity that makes the residence time) is the
	! value at r_j.
	!
	! The recursion is started, and restarted after every cell it does not
	! integrate (a refused or failed cell), by backward Euler,
	! y_j - y_{j-1} = h_j F(y_j, r_j): one step of local error O(h^2), which
	! leaves the global order two. A BDF2 step whose populations or root are
	! not admissible (a negative population, no temperature root) is retaken
	! by backward Euler, which keeps the populations of a linear rate system
	! nonnegative; the log counts both.
	!
	! WHY NOT THE TRAPEZOID RULE. It is second order and A-stable but not
	! L-stable: its amplification factor (1 + z/2)/(1 - z/2) tends to -1 as
	! z = -Da -> -infinity, so a population that relaxes within the cell does
	! not settle on its local equilibrium but alternates about it from cell to
	! cell. BDF2, like backward Euler, is L-stable: both amplification roots
	! of its step tend to zero there, so in the stiff limit it returns the
	! local equilibrium.
	!
	! ZERO-STABILITY. Variable-step BDF2 is zero-stable for step ratios
	! w <= 1 + sqrt(2) (R. D. Grigorieff 1983, Numer. Math. 42, 359, as
	! cited by Li & Liao above; the 1983 paper was not read); a cell whose
	! ratio exceeds that bound is integrated by backward Euler. The Mixed
	! grids of the LHS 1140 b states have w <= 1.018 (MEASURED 2026-09-27).
	real*8, parameter :: bdf2_step_ratio_limit = 2.414213562373095d0

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
	! step it is built on loses its meaning. Unity is not a
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
	! differenced by the step of the recursion (bdf2_step_ratio_limit): every
	! r-derivative at r_j is (f_j - f_hist)/dr_step, with f_hist the history
	! a1 f_{j-1} - a2 f_{j-2} (backward Euler: f_{j-1}) and dr_step = g h_j
	! (backward Euler: h_j). With the common 1/dr_step divided out its three
	! terms are, in the variables of the residual,
	!
	!     q_adv  = |rho v (e_j - e_hist)|
	!     q_prs  = |w v (rho_j - rho_hist)|          w = p/rho
	!     q_enth = |h div(rho v) dr_step|            h = e + w
	!
	! The third vanishes for a stationary mass flux and for nothing else, so
	! this ratio is how much of the equation the non-stationarity of the
	! recorded flow carries, and above one the temperature the equation
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
	! The arguments are the cell (subscript up) and the history of the step
	! into it (subscript hist), the values the residual evaluates the state at.
	pure real*8 function enthalpy_flux_term_ratio(h_up, div_rhov, dr_step, &
	                     rho_up, v_up, e_up, e_hist, w_up, rho_hist)
	real*8, intent(in) :: h_up      ! enthalpy per unit mass of the cell
	real*8, intent(in) :: div_rhov  ! divergence of the mass flux [1/length]
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
                                  nheiTR_in, nm_in, f_sp_in)


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


   ! 60 entries to match the ion_system_HeH_metals / ion_system_metals_pp
   ! dummy length used by the pp_metals=2 re-solve (only 1-11 are consumed;
   ! the advection systems use 1-22).
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

	! Iterate the post processing
	do k = 1,10	! Usually 10 gives a good convergence

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
					 dum_v1,dum_v2, P_m2=P_m2)
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
		                     dheat_hrc, dP_m2=dP_m2_hrc)
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

   !----------------------------------!

	!---- Where is the advection correction valid? ----!
	!
	! The correction replaces the local ionization balance of a cell by the
	! steady advection-ionization ODE, integrated along the flow by the step
	! of the recursion (bdf2_step_ratio_limit). Three
	! conditions make that replacement carry no information; where any of them
	! holds the cell keeps the converged equilibrium ionization instead.
	!
	!  (i)   Inflow, v <= 0 in the cell or the one below it -- PHYSICAL. The
	!        step takes the upstream state from the cells below, which are not
	!        upstream when the gas moves inward (the breathing base). The
	!        residence time h/v is then negative as well.
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
	!                           + (q31g + q31a + q31b + a_ion_HeITR)*n_e
	!                           + Q31*n_HI
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
	! for the step to read.

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
			! The residence time of the gas in the cell, h_j/v_j: the
			! velocity of the cell the rates are evaluated in, as in the step
			! (bdf2_step_ratio_limit).
			t_cross  = (r(j) - r(j-1))*R0/(v(j)*v0)
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
					     + (q31g(j) + q31a(j) + q31b(j) + a_ion_HeITR(j))*ne(j) &
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
		! (i) and (iv) say about it -- the energy step has no
		! upstream state under inflow either, and a non-stationary cell has
		! no steady energy equation -- and the energy loop adds the rest.
		if (adv_correction_valid(j)) then
			adv_comp_status(j) = adv_corrected
		else
			adv_comp_status(j) = adv_retained
		endif
	enddo

   !----------------------------------!

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

			! Add more if HeITR is present
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
	! (heating_of_current_composition, below)
	call heating_of_current_composition(theat)

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

	!---- What the energy correction solves ----!
	!
	! The equation solved cell by cell below (T_equation) is the steady
	! internal-energy equation of the supplied profile,
	!
	!     div(u v) + p div(v)  =  heating - cooling
	!   = rho v de/dr  -  p v dln(rho)/dr  +  h div(rho v)             (E)
	!
	! differenced by the step of the recursion (bdf2_step_ratio_limit;
	! the residual term by term in T_equation), with u the internal
	! energy density, e the internal
	! energy per unit mass and h = e + p/rho the enthalpy per unit mass. The
	! two forms are one equation: div(rho e v) = e div(rho v) + rho v de/dr,
	! and p div(v) = (p/rho) div(rho v) - p v dln(rho)/dr.
	!
	! The third term is the enthalpy flux carried by the divergence of the
	! mass flux, h div(rho v) = h (1/r^2) d(rho v r^2)/dr. It vanishes for a
	! stationary mass flux rho v r^2 = const and for nothing else, so a
	! correction built without it describes a converged wind and no other
	! state. It is handed to the residual through teq_cell%div_rhov, whose
	! algebra T_equation states.
	!
	! WHICH DIVERGENCE, AND WHY. div(rho v) is the MASS ROW OF THE STATE: the
	! face mass fluxes of the Riemann solve of the state differenced over the
	! cell's own control volume, (A_{j+1/2} F_{j+1/2} - A_{j-1/2} F_{j-1/2})/dV_j
	! (RK_rhs), the operator condition (iv) above and the stationary
	! certification read. It is the divergence at the cell centre to second
	! order, the point at which the step evaluates every other factor, and
	! it is an exact zero wherever the two faces carry the same flux.
	!
	! Until 2026-09-27 the same operator was applied to the CELL-CENTRED
	! rho v r^2 of r(j-1) and r(j), the two points of the upwind difference
	! used then. That quantity is not the mass flux of a finite-volume state,
	! and at the base of a quasi-hydrostatic layer it carries the odd-even
	! pattern of the cell-centred velocity. MEASURED on the certified
	! LHS 1140 b state of He/H = 1.6 on 500 cells: at 1.0010-1.0014 R_p the
	! cell-centred rho v r^2 steps between 1.8e12 and 3.4e12 (v = 0.13,
	! 0.25, 0.21 cm/s) while the face mass flux changes by 1e-10 of itself
	! across each cell (the adv_mass_row column). Its difference put an
	! enthalpy flux 77 times the kept terms into the equation at 1.0010 R_p,
	! and the temperature returned at the first corrected cell, 805 / 4008 /
	! 579 K on 500 / 1000 / 2000 cells against a run temperature of 1465 /
	! 1272 / 1039 K there, was set by it and carried by the recursion up to
	! 1.1 R_p.
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
	! marching a column of constant density and F = r with this step
	! reproduces that solution at second order in the cell width: departure
	! 4.4e-4, 1.1e-4, 2.8e-5, 7.0e-6 on 50 to 400 cells of a stretched grid,
	! observed order 1.99, against 0.999 for backward Euler at every step
	! (MEASURED 2026-09-27, src/tests/adv_static_limit, rows N4 and N5).

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
		if (.not. energy_class_modelled(j)) then
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

	! THE HEATING OF THE COMPOSITION THE PASS HOLDS NOW [code units]. The
	! photoheating of ONE particle of each absorber, from the same attenuated
	! field the equilibrium pass used, and then the ONE heating assembly
	! (utils_ion_eq) contracted with the densities of the pass. This is the
	! same routine the ionization sweep and the heating breakdown call, so
	! the _adv energy solve balances the heating the run's own energy
	! equation deposits and cannot drift from it.
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
	         .false., .false., heat_out, heat_chan_pp)

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
	! step and its rate weight c1 = (g) h_j/v_j, then the cell's advection
	! system solved from the fractions of the last pass. The rates of the
	! cell are loaded by the caller; the result is left in sys_x and info.
	subroutine composition_step(jc, bdf2)
	integer, intent(in) :: jc
	logical, intent(in) :: bdf2

	adv_cell%c1       = step_width(jc, bdf2)*R0/(v(jc)*v0)
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
	! Metal densities for this cell [cgs] go through the equation_T module
	! array (the 27-ion vector does not fit params). pp_metal_on gates
	! whether T_equation adds the metal cooling/brem/n_e terms.
	pp_nm_cell(:)  = nm_w(jc,:)
	pp_beta_fs(:)  = beta_fs_pp(jc,:)
	pp_nbar_fs(:)  = nbar_fs_pp(jc,:)

	outcome = T_root_found
	T_root  = T_in(jc)

	! Solve the scalar energy equation by bracketing the physical (lowest)
	! root + Brent when metal cooling is on (the default). The metal-cooled
	! residual is non-monotone and has a second, spurious hot root that a
	! Newton/Powell solve (hybrd1) could land on; bracketing from below
	! selects the physical root structurally. With "Brent solver: False"
	! (use_brent_tsolve = .false.) the MINPACK solve + 2x-band reject is used
	! instead. Metals-off always takes the MINPACK solve (monotone residual).
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
		       abs(teq_cell%mup*teq_cell%mum*dr_step*theat(jc)))) then
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
