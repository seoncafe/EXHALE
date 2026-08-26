	module post_processing
	! Subroutine to correct the output ionization profiles
	!	taking into account the ionization term
	!
	! Composition approximation (documented): the _adv reconstruction treats the
	! gas as H/He + trace metals and EXCLUDES the molecular species (H2, H2+,
	! H3+, HeH+). The advection correction receives only the H/He (nhi..nheiTR)
	! and metal (nm_in) densities; the molecular densities are not passed in and
	! are not re-solved here, so calc_ne / calc_ntot below are called WITHOUT the
	! nmol argument. In a molecular run this omits the neutral-H2 particle count
	! and the molecular-ion electrons from the _adv n_tot/ne, and nh = nhi+nhii
	! counts only the free H nuclei. This is acceptable where the _adv
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
	use output_write
	use equation_T
	use opacity_models           ! opacity_pT_factor for the 'P' model

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

	! Advection-corrected H/He ionized fractions for the current cell, pinned
	! while the metal re-solve (pp_metals=2) adjusts only the metal stages.
	! Set per cell before each ion_system_metals_pp / hybrd1 call.
	real*8, save :: pp_xHII_fix   = 0.0d0
	real*8, save :: pp_xHeII_fix  = 0.0d0
	real*8, save :: pp_xHeIII_fix = 0.0d0

	contains
	
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
      
 	! Updated species densities
   real*8, dimension(1-Ng:N+Ng) :: nhi,nhii
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


   ! 60 entries to match the ion_system_HeH_metals / ion_system_metals_pp
   ! dummy length used by the pp_metals=2 re-solve (only 1-11 are consumed;
   ! the advection systems use 1-22).
   real*8, dimension(60) :: params
   real*8, dimension(40) :: paramsT
      
	real*8 :: rhop,rhom,vp,mum,mup,vm
	real*8 :: sys_sol_T(1), sys_x_T(1)
   real*8 :: wa_T(8)
   logical :: brent_ok                       ! Task 1: Brent T-solve status
      
      
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
		nhei   = nhei_in*n0
		nheii  = nheii_in*n0
		nheiii = nheiii_in*n0
		! Zero-init the triplet so the innermost ghost cell (1-Ng), which
		! the advection loop below never assigns when HeITR is off, does
		! not write uninitialized memory to the output column.
		nheiTR = 0.0
		if (thereis_HeITR) nheiTR = nheiTR_in*n0
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

	! Iterate the post processing
	do k = 1,10	! Usually 10 gives a good convergence
	
	! Use singlet if included
	nheiS = nhei
	if (thereis_HeITR) nheiS = nhei - nheiTR
	nh  = nhi  + nhii
	nhe = nheiS + nheii + nheiii
	if (thereis_HeITR) nhe = nhe + nheiTR

	! Ionized fraction of the H+He nuclei, for the SvS85 secondary-ionization
	! partition (metals excluded; nheiTR is neutral and not in the numerator).
	! Reused for both PH_heat calls below (nhi/nheii/nheiii unchanged between them).
	xion = min(max((nhii + nheii + nheiii)/max(nh + nhe, 1.0d-99), 0.0d0), 1.0d0)

	! Free electron density (assuming overall neutrality; nm_w adds the
	! metal electrons under the eos_metals policy). Molecular-ion electrons are
	! excluded here -- the post-process does not carry the molecular densities
	! (see the module-header composition note).
	call calc_ne(nhii,nheii,nheiii,ne,nm_w)

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
				   dum_v1, a_ion_HeITR=a_ion_HeITR)

 	
 	if (thereis_HeITR) then
		call HeITR_coeffs(T_K,rcheiTR,rcheiiB,A31,q13,q31a,q31b,Q31)
		! NOTE: rcheiiB is alpha1 from Oklopcic - being overwritten

	endif

	! He recombination radiation ionizing H I (Draine 2011; default off).
	! Correct the He II recombination coefficient and the H I photoionization
	! rate driving the advection ODE (the heating correction is applied later,
	! to theat). Mirrors ionization_equilibrium.
	if (use_he_rec_coupling .and. thereis_He) then
		call he_rec_coupling(T_K, nhi, nhei, nheii, nheiTR, ne,           &
		                     A31, q31a, q31b,                             &
		                     rcheiiB_hrc, dP_HI_hrc, dheat_hrc)
		rcheiiB = rcheiiB_hrc
		P_HI    = P_HI + dP_HI_hrc
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

	adv_correction_valid = .true.
	adv_correction_valid(1-Ng) = .false.   ! inner boundary: never corrected
	n_adv_eq = 0
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
		if (.not. adv_correction_valid(j)) n_adv_eq = n_adv_eq + 1
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
			
			! Extract solution profiles	
			nhi(j)    = sys_x(1)*nh(j)
			nhii(j)   = (1.0 - sys_x(1))*nh(j)

      	enddo
      	
		   ! Force condition of zero helium
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
				nhi(j)    = nhi_in(j)*n0
				nhii(j)   = nhii_in(j)*n0
				nhei(j)   = nhei_in(j)*n0
				nheii(j)  = nheii_in(j)*n0
				nheiii(j) = nheiii_in(j)*n0
				nheiTR(j) = 0.0
				if (thereis_HeITR) nheiTR(j) = nheiTR_in(j)*n0
				cycle
			endif

			! Substitutions
			dr  = (r(j) - r(j-1))*R0
			As  = dr/(v(j-1)*v0)

			! Advection coeff.
			adv_cell%c1  = As
			adv_cell%xhi_old  = nhi(j-1)/nh(j-1)
			adv_cell%xhei_old  = nhei(j-1)/nhe(j-1)
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
			sys_x(2) = nhei(j)/nhe(j)
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
				
			! Extract solution profiles	
			nhi(j)    = sys_x(1)*nh(j)
			nhii(j)   = (1.0 - sys_x(1))*nh(j)
			nhei(j)   = sys_x(2)*nhe(j)
			nheii(j)  = (1.0 - sys_x(2) - sys_x(3))*nhe(j)
			nheiii(j) = sys_x(3)*nhe(j)
			if (thereis_HeITR) then
				nheiTR(j) = sys_x(4)*nhe(j) 
			else
				nheiTR(j) = 0.0
			endif
			
		enddo
		
	endif ! End if thereis_He
	      
      
   !---------------------------------------------------------
      
   !------- Fix stationarity of new pressure profile -------!
      
   !---- Update densities and temperature ----!
	
	! Number densities      
   nheiS = nhei
	if (thereis_HeITR) nheiS = nhei - nheiTR
	nh  = nhi  + nhii 
	nhe = nheiS + nheii + nheiii
	if (thereis_HeITR) nhe = nhe + nheiTR

   ! Total number density (incl. metal nuclei under eos_metals). Molecular
   ! species are excluded -- the post-process does not carry them (see the
   ! module-header composition note).
   call calc_ntot(nhi,nhii,nhei,nheii,nheiii,nheiTR,n_tot,nm_w)

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

	!---- Update photoheating rate ----!
	
	if (thereis_He) then
		call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm_w, xion,     &
		                 dum_v1,dum_v2,dum_v3,dum_v4, P_m,          &
		                 theat,dum_v5)
  	else
	  	call PH_heat_H(nhi, xion, dum_v1,theat,dum_v2)
  	endif

	! He recombination radiation ionizing H I (Draine 2011; default off):
	! photoelectron heating from the coupled H ionizations, evaluated at the
	! advection-corrected densities (the rate/coefficient corrections were
	! applied to P_HI/rcheiiB before the advection solve above).
	if (use_he_rec_coupling .and. thereis_He) then
		call he_rec_coupling(T_K, nhi, nhei, nheii, nheiTR, ne,           &
		                     A31, q31a, q31b,                             &
		                     rcheiiB_hrc, dP_HI_hrc, dheat_hrc)
		theat = theat + dheat_hrc
	endif

	! Penning ionization heating (mirrors ionization_equilibrium): He(2^3S)+H0
	! -> He(1^1S)+H+ + e- releases e_th_HeI - e_th_HeTR - e_th_HI (= 6.2 eV).
	! theat here is the freshly recomputed photoheating (+ he_rec_coupling), so
	! this term is added once and is not double-counted. Advection-corrected
	! densities.
	if (thereis_HeITR) theat = theat                                 &
	     + nheiTR*nhi*Q31*(e_th_HeI - e_th_HeTR - e_th_HI)/erg2eV

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

	do j = 3-Ng,N+Ng ! Start from first computational cell
		
		! Substitutions
		rhop = rho(j)
		rhom = rho(j-1)
		vm = v(j-1) 
		vp = v(j)
		! Inflow (condition (i) of the validity block above): the cell keeps the
		! converged eq temperature. The advection-corrected energy solve is
		! upwind-differenced just like the ionization solve, so it is invalid
		! wherever the gas moves inward, and it would otherwise land on the
		! spurious hot root that then cascades up. This is a property of the
		! discretization, so it does not depend on the metal switch. Conditions
		! (ii) and (iii) are statements about the ionization balance and its
		! representation and are deliberately NOT applied here.
		if (vp <= 0.0d0 .or. vm <= 0.0d0) then
			T_out(j) = T_in(j)
			cycle
		endif
		dr = r(j) - r(j-1)
		mum = mmw(j-1)
		mup = mmw(j)

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
	 	if (pp_metal_on .and. use_brent_tsolve) then
	 		call solve_T_brent(paramsT, T_in(j), sys_x_T(1), brent_ok)
	 		if (.not. brent_ok) then
	 			sys_x_T(1)  = T_in(j)
	 			n_pp_reject = n_pp_reject + 1
	 		endif
	 	else
	 		! Legacy MINPACK solve.
	 		call hybrd1(T_equation,1,sys_x_T,sys_sol_T,   &
	 		            tol,info,wa_T,8,paramsT)
	 		! With metal cooling, reject a non-physical / out-of-band root (the
	 		! spurious hot root) and fall back to eq T (the legacy guard).
	 		if (pp_metal_on) then
	 			if (.not. (sys_x_T(1) > 0.0d0)  .or.   &
	 			    sys_x_T(1) > 2.0d0*T_in(j)  .or.   &
	 			    sys_x_T(1) < 0.5d0*T_in(j)) then
	 				sys_x_T(1)  = T_in(j)
	 				n_pp_reject = n_pp_reject + 1
	 			endif
	 		endif
	 	endif

		! Extract solution profiles
		T_out(j) = sys_x_T(1)
	 	
	enddo
	
	! Update pressure and temperature
	p_out = (n_tot + ne)/n0*T_out
	T_K = T_out*T0
      
	enddo ! End loop on post processing

	! Report how many cells the advection correction was not applied to
	! (inflow, local ionization equilibrium, or an unrepresentable ion
	! fraction; see the validity block above). Counted on the last pass.
	if (n_adv_eq > 0) then
		write(*,'(a,i0,a,i0,a)') ' (post_process_adv) advection correction: ', &
		   n_adv_eq, ' of ', N+2*Ng-1,                                          &
		   ' cells kept at the equilibrium ionization.'
	endif

	! Report how many base cells fell back to the eq temperature (metal modes).
	if (pp_metal_on .and. n_pp_reject > 0) then
		write(*,'(a,i0,a,i0,a)') ' (post_process_adv) metal-mode T solve: ', &
		   n_pp_reject, ' of ', N+2*Ng-2,                                     &
		   ' cells fell back to eq T (stiff base band).'
	endif

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
   call write_output(rho,v,p_out,T_out,theat,tcool,eta,                &
                     nhi_w,nhii_w,nhei_w,nheii_w,nheiii_w,              &
                     nheiTR_w,nm_out,'ad')

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
