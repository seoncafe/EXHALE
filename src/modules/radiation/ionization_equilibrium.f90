   module ionization_equilibrium
	! Evaluate the ionization structure and the heating and cooling functions for a given temperature

	use global_parameters
   use ion_cell_state, only: ieq_cell
   use species_table, only: n_mion, mion_fsp, n_melem, melem_i0,        &
                            melem_top, mion_stage,                       &
                            isp_H2, isp_H2p, isp_H3p, isp_HeHp
   use utils
   use utils_ion_eq
   use Cooling_Coefficients      ! eval_cool, recombination/ionization rates
   use System_HeH                ! Equilibrium equations
	use System_HeH_TR
	use System_HeH_mol            ! molecular network
	use System_HeH_mol_metals     ! merged molecular network + metals
	use lower_column, only: q_h2_equilibrium
	use lyman_werner_photodissociation, only:                            &
	                           lyman_werner_dissociation_rate,            &
	                           h2_self_shielding_factor,                  &
	                           h2_doppler_parameter, e_lw_fragment_erg
   use System_HeH_metals
   use System_HeH_TR_metals      ! merged He-triplet + metals system
   use charge_exchange, only: cx_set_cell, cx_metal_base,  &   ! Huang Table 4 charge exchange
                              he_h_cx_rates                     ! He <-> H pair (group B)
   use System_H
   use newton_solver, only: solve_ieq   ! Task 2: analytic-Jacobian Newton (+ hybrd1 fallback)
   use opacity_models            ! opacity_pT_factor for the 'P' model

   implicit none

	! molecular species densities (cols 1 H2, 2 H2+, 3 H3+, 4 HeH+;
	! zero unless thereis_mol).  Module state: written by the equilibrium
	! solve, read by write_output for the extra output columns.
	real*8, dimension(:,:), allocatable :: nmol_eq

	! Lyman-Werner photodissociation diagnostics, filled by the equilibrium
	! solve when a band flux is supplied and read by write_lyman_werner:
	! star-ward H2 column [cm^-2], the DB96 self-shielding factor, and the
	! resulting dissociation rate [s^-1]. Untouched (shielding 1, rate 0)
	! when the run supplies no band flux.
	real*8, dimension(:), allocatable :: NH2_col_lw
	real*8, dimension(:), allocatable :: f_shield_lw
	real*8, dimension(:), allocatable :: k_lw_diss

	! Run-wide totals of the atomic ionization root validation, reported once
	! at the end of the run (EXHALE_main) next to the Newton usage counters:
	! stored states rejected as a starting point, cell solves that needed a
	! second or third starting point, first roots that lay outside the
	! physical simplex, and cells where no starting point produced an
	! admissible root. All zero for a run that never leaves the simplex.
	integer, save :: ieq_n_reseed = 0, ieq_n_retry  = 0
	integer, save :: ieq_n_unphys = 0, ieq_n_noroot = 0

	! Run-wide total of molecular cells whose equilibrium roots all left the
	! physical simplex, so the closest one was clamped back onto the element
	! budget (ioniz_eq). Zero for a run whose molecular solve stays inside the
	! simplex everywhere, and zero for an atomic run.
	integer, save :: ieq_n_mol_clamped = 0

	! Run-wide histogram of the hybrd1 exit code of the molecular cell solves,
	! indexed by info (0 = improper input or iflag < 0, 1 = converged to tol,
	! 2 = iteration limit, 3 = xtol too small, 4/5 = no progress). Every
	! attempt of every molecular cell is counted, so the total exceeds the
	! cell count whenever a cell needs its second or third starting point.
	! Zero for an atomic run.
	integer, save :: ieq_n_mol_info(0:5) = 0

	contains

	subroutine ioniz_eq_allocate_arrays
	! Allocate the grid-sized module arrays once the number of cells N is
	! known; called from EXHALE_main right after input_read. The values are
	! the initializers the declarations used to carry.

	allocate(nmol_eq(1-Ng:N+Ng,4))
	allocate(NH2_col_lw(1-Ng:N+Ng), f_shield_lw(1-Ng:N+Ng),               &
	         k_lw_diss(1-Ng:N+Ng))

	nmol_eq     = 0.0d0
	NH2_col_lw  = 0.0d0
	f_shield_lw = 1.0d0
	k_lw_diss   = 0.0d0

	end subroutine ioniz_eq_allocate_arrays

	subroutine ioniz_eq(T_in,n_io,f_sp_io,heat_out,cool_out,q)
      	 		  
	integer :: j,im
	logical :: usednt                     ! Task 2: Newton-vs-fallback flag

	real*8, dimension(1-Ng:N+Ng),   intent(in) :: T_in
	! Density and species fractions are read on entry and overwritten with the
	! equilibrium result on exit (in place); callers must not alias them.
	real*8, dimension(1-Ng:N+Ng),   intent(inout) :: n_io
	real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp_io

	real*8, dimension(1-Ng:N+Ng) ::  T_K      ! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng) ::  nh,nhi,nhii,                   & ! Species densities
	                                 nhe,nhei,nheii,nheiii,nheiTR,  &
	                                 ne,n_in_dim,n_tot
	! Ionized fraction of the H+He nuclei, for the SvS85 secondary ionization.
	real*8, dimension(1-Ng:N+Ng) ::  xion
   ! Metal ion densities in canonical species_table order (col im maps
   ! to f_sp column mion_fsp(im)); used throughout in place of named
   ! scalars for each ion so the driver scales with the number of metals.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  nm
   ! Total density of each metal element (canonical element order),
   ! held constant across the cell sweep (sum of its three stages).
   real*8, dimension(1-Ng:N+Ng,n_melem) ::  nm_tot

   ! Metal recombination and collisional ionization rates for each ion from
   ! eval_cool (canonical order); bridged to the named rc*/a_ion_*
   ! scalars below for the params packing.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  rec_m,aion_m

   ! Photo ionization rates
   real*8, dimension(1-Ng:N+Ng) ::  P_HI,P_HeI,P_HeII,P_HeITR
   real*8, dimension(1-Ng:N+Ng) ::  P_H2      ! (molecular; zero unless mol)
   ! Doppler parameter of H2 in a cell, for the Lyman-Werner shielding.
   real*8 :: b_h2_lw
   ! Metal photoionization rates for each ion (canonical order) from PH_heat.
   real*8, dimension(1-Ng:N+Ng,n_mion) ::  P_m
                       	
   ! Heating, cooling
   real*8, dimension(1-Ng:N+Ng) ::  heat,cool    
                 
   ! Recombination coefficients
   real*8, dimension(1-Ng:N+Ng) ::  rchiiB,rcheiiB,rcheiiiB,rcheiTR

   ! He recombination radiation -> H ionization coupling scratch
   ! (use_he_rec_coupling; zero-effect when off).
   real*8, dimension(1-Ng:N+Ng) ::  rcheiiB_hrc,dP_HI_hrc,dheat_hrc

	real*8, dimension(1-Ng:N+Ng) :: q13,q31a,q31b,Q31
	real*8 :: A31

   ! Ionization coefficients
   real*8, dimension(1-Ng:N+Ng) ::  a_ion_HI,a_ion_HeI,a_ion_HeII,a_ion_HeITR

	! Equilibrium system setup
   real*8 :: tol,dpmpar
   real*8, dimension(60) :: params
   ! Each element's metal coefficients for ion_system_HeH_metals (canonical order),
   ! built per cell from the 2D rate arrays and handed to set_metal_coeffs.
   real*8, dimension(n_melem) :: meg_ntot,meg_g0,meg_g1,meg_b0,meg_b1, &
                                 meg_a1,meg_a2
   ! Highest stage per element handed to set_metal_coeffs (2 = three-stage,
   ! 1 = two-stage); see species_table::melem_top.
   integer, dimension(n_melem) :: meg_top
   integer :: i0,top
   ! Base index of the first metal element's X+ unknown in sys_x: 4 normally,
   ! 5 when the He triplet occupies x(4) (merged HeITR+metals system), and
   ! metal_row_base() = 8 or 9 in the molecular layout, where x(4..7) are the
   ! H nuclei bound into H2/H2+/H3+/HeH+ and x(8) the triplet.
   integer :: mbase
   ! Molecular solve (thereis_mol): count of cells whose roots all left the
   ! physical simplex, so the closest one was clamped back onto the element
   ! budget, and the local gas pressure [bar] the H2 dissociation
   ! equilibrium of the retry seed is evaluated at.
   integer :: n_mol_clamped
   ! Per-sweep histogram of the molecular hybrd1 exit code (see
   ! ieq_n_mol_info above); combined over threads by the reduction below.
   integer :: n_mol_info(0:5)
   real*8  :: pbar_loc

   ! Ionization solve: validation of the returned root.
   ! The equilibrium systems are polynomial and possess roots outside the
   ! physical simplex (negative stage fractions, or ionized stages of one
   ! element summing above its nucleus total); such a root makes the neutral
   ! density negative and is not a solution of the physical problem even when
   ! the algebraic residual vanishes. Each cell therefore tries up to three
   ! starting points and keeps the best admissible root; x_root_best holds it,
   ! ok_rank grades the attempt (0 none, 1 physical, 2 physical + converged),
   ! and the counters feed the one-line sweep summary.
   ! Largest N_eq over the layouts: molecular (7) + triplet (1) + metals.
   integer, parameter :: n_x_max = 8 + 2*n_melem
   integer :: iatt, info_ieq, ok_rank, best_rank
   integer :: n_ieq_reseed, n_ieq_retry, n_ieq_unphys, n_ieq_fail
   logical :: conv_ieq, phys_ieq
   real*8, dimension(n_x_max) :: x_root_best
   ! Molecular cell whose roots all left the physical simplex: how far the
   ! least-offending one lies outside it, so the closest to a state can be
   ! selected and clamped back onto the element budget.
   real*8  :: viol, viol_best

   ! Output heating,cooling and absorbed energy
   real*8, dimension(1-Ng:N+Ng),intent(out) :: heat_out,cool_out,q

   !----------------------------------------------------------!      
   ! Global parameters
      
   ! Numerical tolerance for system solution
   tol = sqrt(dpmpar(1))

	!----------------------------------!
	
	! Preliminary profiles exctraction
	
	! Dimensional total number density profile and temperature
	n_in_dim = n_io*n0
	T_K      = T_in*T0
		
	! Extract species profiles
	nhi    = f_sp_io(:,1)*n_in_dim    ! HI
	nhii   = f_sp_io(:,2)*n_in_dim    ! HII
	if (thereis_He) then
		nhei   = f_sp_io(:,3)*n_in_dim    ! HeI
		nheii  = f_sp_io(:,4)*n_in_dim    ! HeII
		nheiii = f_sp_io(:,5)*n_in_dim    ! HeIII
		nheiTR = f_sp_io(:,6)*n_in_dim    ! HeITR
	else
		! Enforce condition of zero helium
		nhei   = 0.0
		nheii  = 0.0
		nheiii = 0.0
		nheiTR = 0.0
	endif

    ! Metal ion densities in canonical species_table order (col im maps to
    ! f_sp column mion_fsp(im)).
    do im = 1,n_mion
       nm(:,im) = f_sp_io(:,mion_fsp(im))*n_in_dim
    enddo

	! molecular species (zero when thereis_mol is off)
	if (thereis_mol) then
		nmol_eq(:,1) = f_sp_io(:,isp_H2)  *n_in_dim
		nmol_eq(:,2) = f_sp_io(:,isp_H2p) *n_in_dim
		nmol_eq(:,3) = f_sp_io(:,isp_H3p) *n_in_dim
		nmol_eq(:,4) = f_sp_io(:,isp_HeHp)*n_in_dim
	endif

    ! Total density of each metal element (sum of its three stages, in the
    ! same neutral+singly+doubly order as the original nC=nci+ncii+nciii).
    do im = 1,n_melem
       i0 = melem_i0(im)
       if (melem_top(im) .ge. 2) then
          nm_tot(:,im) = nm(:,i0) + nm(:,i0+1) + nm(:,i0+2)
       else
          ! Two-stage element: neutral + singly ionized only.
          nm_tot(:,im) = nm(:,i0) + nm(:,i0+1)
       endif
    enddo
	
	! Total number densities (with molecules: H and He NUCLEI totals --
	! the molecular system conserves elements, and nmol_eq is zero otherwise)
	nh  = nhi  + nhii
	nhe = nhei + nheii + nheiii
	if (thereis_mol) then
		nh  = nh  + 2.0d0*(nmol_eq(:,1) + nmol_eq(:,2))                  &
		          + 3.0d0*nmol_eq(:,3) + nmol_eq(:,4)
		nhe = nhe + nmol_eq(:,4)
	endif

	! Ionized fraction of the H+He nuclei, for the SvS85 secondary-ionization
	! partition (metals excluded; nheiTR is neutral and not in the numerator).
	xion = min(max((nhii + nheii + nheiii)/max(nh + nhe, 1.0d-99), 0.0d0), 1.0d0)

	! Free electron density (assuming overall neutrality; nm adds the
	! metal electrons under the eos_metals policy, nmol_eq the molecular-ion
	! electrons -- it is zero for an atomic run, so the sum is unchanged there)
	call calc_ne(nhii,nheii,nheiii,ne,nm,nmol_eq)

	! Total gas-particle density (electrons excluded), i.e. the density of
	! third bodies M for the three-body molecular reactions R12/R13/R15, and
	! -- with the electrons added back -- the gas pressure p = (n_tot+n_e)kT
	! of the chemical-equilibrium retry seed below.  It is NOT n_in_dim:
	! n_in_dim = rho/m_H is a MASS density in m_H units (calc_rho weights
	! each species by bsp_mass), so it over-counts the particles by the mean
	! particle mass -- 2.3x at an H2-rich base, up to 4x in a He-dominated
	! one.  Only the molecular path reads it, so an atomic run is unchanged.
	if (thereis_mol) then
		call calc_ntot(nhi,nhii,nhei,nheii,nheiii,n_tot,nm,nmol_eq)
	else
		n_tot = 0.0d0
	endif

	! Cell-by-cell pressure-broadening factor for the opacity ('P' model).
	! opacity_pT_factor returns 1.0 for all other models, so opa_pf=1
	! and the column densities are unchanged (bit-identical).
	do j = 1-Ng,N+Ng
		opa_pf(j) = opacity_pT_factor((nh(j)+nhe(j)+ne(j))*kb_erg*T_K(j))
	enddo

	!----------------------------------!

	! Lyman-Werner photodissociation of H2 (Draine & Bertoldi 1996; see
	! src/modules/lower_atmosphere/lyman_werner.f90). The star-ward H2
	! column uses the same radial integration and opa_pf weighting as every
	! other absorber column, and the incoming (pre-solve) H2 density, like
	! the photoionization columns built inside PH_heat_HHe. The band lies
	! below the H I edge, so nothing else in the model absorbs it: there is
	! no dust, and the trace-metal continuum is negligible against the H2
	! line self-shielding (module header).
	if (thereis_mol .and. F_LW_star .gt. 0.0d0) then
		call calc_column_dens_one(nmol_eq(:,1), NH2_col_lw)
		do j = 1-Ng,N+Ng
			b_h2_lw = h2_doppler_parameter(T_K(j))
			f_shield_lw(j) = h2_self_shielding_factor(NH2_col_lw(j), b_h2_lw)
			k_lw_diss(j)   = lyman_werner_dissociation_rate(F_LW_star,     &
			                                       NH2_col_lw(j), T_K(j))
		enddo
	endif

	!----------------------------------!

    !---- Photoionization and photoheating ----!
      
	if (thereis_He) then
		if (thereis_mol) then
			call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm, xion,    &
			         P_HI,P_HeI,P_HeII,P_HeITR, P_m,             &
			         heat,q, nmol_eq(:,1),P_H2)
		else
      	call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm, xion,      &
      			     P_HI,P_HeI,P_HeII,P_HeITR, P_m,        &
      			     heat,q)
		endif
	else
		call PH_heat_H(nhi, xion, P_HI,heat,q)
		P_m = 0.0
	endif

	! excited-H H(n=2) feedback (zero unless use_excited_H). The Balmer
	! photoionization of H(n=2) adds an effective HI photoionization rate
	! [s^-1] (proton source), and the Balmer photoelectric (+ optional
	! collisional de-excitation) heating adds to the photoheating rate
	! [erg cm^-3 s^-1]. Both global arrays are filled from the previous
	! converged outer pass by excited_hydrogen::excited_H_update, so the
	! coupling is decoupled from the hydro sub-step. Default-off => no change.
	if (use_excited_H) then
		gph_ground_HI = P_HI         ! capture pure ground-state rate (pre-n=2)
		P_HI = P_HI + gph_balmer_HI
		heat = heat + heat_balmer
	endif


   !----------------------------------!
      
   ! Evaluate cooling rates and recombination/collisional 
   ! 	ionization rates
      
    ! nmol_eq gives the cooling the same electron density this routine
    ! balances the ionization against, and the H3+/H2 densities its infrared
    ! cooling channel needs (zero for an atomic run). Nothing is added to
    ! `cool` after this call: the marching temperature update rebuilds the
    ! cooling from eval_cool alone, so any term added here would be a source
    ! the marching relaxes without and the steady residual demands.
    call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii, nm,           &
			   	   	rchiiB,rcheiiB,rcheiiiB, rec_m,             &
			    	a_ion_HI,a_ion_HeI,a_ion_HeII, aion_m,         &
			    	cool, nheiTR=nheiTR, a_ion_HeITR=a_ion_HeITR,  &
			    	nmol=nmol_eq)

	! Capture the ground-state H proton-budget coefficients on
	! every pass (the converged pass is the one read out by write_excited_H).
	if (use_excited_H) then
		cion_HI  = a_ion_HI          ! collisional ionization [cm^3 s^-1]
		arec_HII = rchiiB            ! case-B recombination  [cm^3 s^-1]
	endif

	! Charge-exchange rate coefficients are evaluated per cell below
	! (cx_set_cell) before each metal ionization solve.

	if (thereis_HeITR) then
		call HeITR_coeffs(T_K,rcheiTR,rcheiiB,A31,q13,q31a,q31b,Q31)
		! NOTE: rcheiiB is alpha1 from Oklopcic - being overwritten
	endif

	! He recombination radiation ionizing H I (Draine 2011 on-the-spot y/z;
	! default off => no change). Uses the pre-solve (lagged) densities, like
	! the other lagged rate terms; corrects the He II recombination coefficient
	! rcheiiB, adds an H I photoionization rate to P_HI, and adds photoelectron
	! heating to heat. The He II recombination *cooling* (rec_cool_HeII =
	! kT alpha_B) is left unchanged; the mismatch is <= y alpha_1 kT, negligible.
	if (use_he_rec_coupling .and. thereis_He) then
		call he_rec_coupling(T_K, nhi, nhei, nheii, nheiTR, ne,           &
		                     A31, q31a, q31b,                             &
		                     rcheiiB_hrc, dP_HI_hrc, dheat_hrc)
		rcheiiB = rcheiiB_hrc
		P_HI    = P_HI + dP_HI_hrc
		heat    = heat + dheat_hrc
	endif

	! Penning ionization heating: He(2^3S)+H0 -> He(1^1S)+H+ + e- releases the
	! electron kinetic energy e_th_HeI - e_th_HeTR - e_th_HI (= 6.2 eV) into the
	! gas. Q31 is the total ionization rate, so only its Penning branch
	! (f_penning_HeI23S) carries this exothermicity; the associative branch
	! ends in HeH+ and has a different one. Lagged (pre-solve) densities, like
	! every other channel above; nheiTR is the same array he_rec_coupling
	! already consumes.
	if (thereis_HeITR) heat = heat                                    &
	     + f_penning_HeI23S*nheiTR*nhi*Q31                             &
	       *(e_th_HeI - e_th_HeTR - e_th_HI)/erg2eV

	! Molecular Penning ionization heating: He(2^3S)+H2 -> He(1^1S)+H2+ + e-
	! releases the electron kinetic energy (e_th_HeI - e_th_HeTR) - e_th_H2
	! (= 24.6 - 4.80 - 15.4 = 4.4 eV) into the gas. Lagged (pre-solve)
	! densities; nmol_eq(:,1) is the neutral-H2 number density. ioniz_HeI23S_H2
	! (Cool_coeff.f90) is the total, scaled here to the Penning branch. Zero
	! unless a molecular run also tracks the triplet.
	if (thereis_mol .and. thereis_HeITR) heat = heat                  &
	     + f_penning_HeI23S*nheiTR*nmol_eq(:,1)*ioniz_HeI23S_H2(T_K)   &
	       *((e_th_HeI - e_th_HeTR) - e_th_H2)/erg2eV

	! Lyman-Werner photodissociation heating: H2 + hv -> H + H leaves the
	! fragment pair with about 0.4 eV of kinetic energy (Black & Dalgarno
	! 1977, ApJS 34, 405, p. 418). The 4.48 eV bond energy is paid by the
	! absorbed photon, not by the gas, so it is NOT a thermal sink of this
	! channel. Lagged (pre-solve) H2 density, like every channel above.
	if (thereis_mol .and. F_LW_star .gt. 0.0d0)                       &
	     heat = heat + k_lw_diss*nmol_eq(:,1)*e_lw_fragment_erg

   !----------------------------------!

   ! Ionization equilibrium system solution

	if (.not.thereis_He) then ! If no helium

		! Ionization solves in each cell are independent (the count>0 warm-start uses
		! this cell's own previous-step value), so the sweep is OpenMP-parallel
		! over cells. sys_x/wa/info are threadprivate (global_parameters); only
		! the subroutine-local scratch is private. count==0 runs serial (the if
		! clause) because its first-step warm-start reads the neighbour cell.
		!$omp parallel do default(shared) schedule(dynamic,8)                  &
		!$omp   private(params, usednt) if(count > 0)
		do j = N+Ng,1-Ng,-1

			! Lazily allocate this thread's threadprivate NL scratch.
			if (.not. allocated(sys_x)) allocate(sys_x(N_eq))
			if (.not. allocated(wa))    allocate(wa(lwa))

			! Ionization equilibrium system setup: named-field cell state
			! (Inc 4). System_H reads these; params is now only the MINPACK
			! transport argument (unread by the converted system).
			ieq_cell%P_HI     = P_HI(j)
			ieq_cell%rchiiB   = rchiiB(j)
			ieq_cell%nh       = nh(j)
			ieq_cell%a_ion_HI = a_ion_HI(j)

			 ! Initial guess
			if (count.le.0) then
				if(r(j).le.(1.5))then
					sys_x(1) = r(j)-0.5
				else
					sys_x(1) = 1.0
				endif     
			else

				sys_x(1) = nhii(j)/nh(j)

			endif

		 	! Analytic-Jacobian Newton (Task 2); hybrd1 fallback inside solve_ieq.
			call solve_ieq(ion_system_H,jac_system_H,N_eq,sys_x,    &
				      params,tol,wa,lwa,usednt)

			! Extract solution profiles
			nhi(j)    = nh(j)*(1.0 - sys_x(1))
			nhii(j)   = nh(j)*sys_x(1)

		enddo
		!$omp end parallel do

		nhei   = 0.0
		nheii  = 0.0
		nheiii = 0.0
		nheiTR = 0.0

	else  ! If there is helium

		! Metal unknowns start at sys_x(4) normally, but shift to sys_x(5)
		! when the He triplet occupies sys_x(4) (merged HeITR+metals system).
		! In the molecular layout the molecular unknowns own sys_x(4..7) and
		! the metals follow at metal_row_base() (single definition, shared
		! with the merged molecular+metals residual).
		if (thereis_mol) then
			mbase = metal_row_base()
		else
			mbase = 4
			if (thereis_HeITR) mbase = 5
		endif

		! OpenMP-parallel cell sweep (see the no-He branch above). The cell-by-cell
		! metal coefficients (met_*, System_HeH_metals) and charge-exchange rates
		! (cx_kc, cx_metal_base) are threadprivate, so each thread keeps its own;
		! cx_metal_base is broadcast (copyin) and toggled 4<->5 per cell. All the
		! subroutine-local scratch is private. count==0 stays serial (neighbour
		! warm-start).
		! Failed-cell counters (combined via the reduction below); reset once
		! per equilibrium sweep. n_mol_clamped counts molecular cells whose roots
		! were all outside the physical simplex; the n_ieq_* counters are the
		! atomic root validation (see the declarations above).
		n_mol_clamped = 0
		n_mol_info(:) = 0
		n_ieq_reseed = 0
		n_ieq_retry  = 0
		n_ieq_unphys = 0
		n_ieq_fail   = 0
		!$omp parallel do default(shared) schedule(dynamic,8) copyin(cx_metal_base) &
		!$omp   private(params, usednt, i0, top, im, meg_ntot, meg_g0, meg_g1,      &
		!$omp           meg_b0, meg_b1, meg_a1, meg_a2, meg_top,                     &
		!$omp           pbar_loc, viol, viol_best,                                  &
		!$omp           iatt, info_ieq, ok_rank, best_rank, conv_ieq, phys_ieq,      &
		!$omp           x_root_best)                                                &
		!$omp   reduction(+:n_mol_clamped,n_mol_info,n_ieq_reseed,n_ieq_retry, &
		!$omp               n_ieq_unphys,n_ieq_fail) if(count > 0)
		do j = N+Ng,1-Ng,-1

			! Lazily allocate this thread's threadprivate NL scratch.
			if (.not. allocated(sys_x))   allocate(sys_x(N_eq))
			if (.not. allocated(sys_sol)) allocate(sys_sol(N_eq))
			if (.not. allocated(wa))      allocate(wa(lwa))

			! System coefficients: named-field cell state (Inc 4). The He
			! equilibrium systems read these; params is now only the MINPACK
			! transport argument (unread by the converted systems).
			ieq_cell%P_HI       = P_HI(j)
			ieq_cell%P_HeI      = P_HeI(j)
			ieq_cell%P_HeII     = P_HeII(j)
			ieq_cell%rchiiB     = rchiiB(j)
			ieq_cell%rcheiiB    = rcheiiB(j)
			ieq_cell%rcheiiiB   = rcheiiiB(j)
			ieq_cell%nh         = nh(j)
			ieq_cell%nhe        = nhe(j)
			ieq_cell%a_ion_HI   = a_ion_HI(j)
			ieq_cell%a_ion_HeI  = a_ion_HeI(j)
			ieq_cell%a_ion_HeII = a_ion_HeII(j)

			! He <-> H charge-exchange rate coefficients (Huang Table 4 group
			! B): read by he_h_cx_fvec/he_h_cx_jac in every He system. Depends
			! only on T, so evaluate once per cell here (cheap). The residual
			! routines add nothing when he_h_charge_exchange is off, so this is
			! harmless (and bit-identical) in that case.
			call he_h_cx_rates(T_K(j), ieq_cell%kcx_He0_Hp,               &
			                           ieq_cell%kcx_Hep_H0)

			! Add more if HeITR is present
			if (thereis_HeITR) then
				ieq_cell%rcheiTR = rcheiTR(j)
				ieq_cell%A31     = A31
				ieq_cell%P_HeITR = P_HeITR(j)
				ieq_cell%q13     = q13(j)
				ieq_cell%q31a    = q31a(j)
				ieq_cell%q31b    = q31b(j)
				ieq_cell%Q31     = Q31(j)
				ieq_cell%a_ion_HeITR = a_ion_HeITR(j)
			endif

			! molecular cell state (System_HeH_mol layout)
			if (thereis_mol) then
				if (.not. thereis_HeITR) then
					ieq_cell%rcheiTR = 0.0d0     ! no triplet channels
					ieq_cell%A31     = 0.0d0
					ieq_cell%P_HeITR = 0.0d0
					ieq_cell%q13     = 0.0d0
					ieq_cell%q31a    = 0.0d0
					ieq_cell%q31b    = 0.0d0
					ieq_cell%Q31     = 0.0d0
					ieq_cell%a_ion_HeITR = 0.0d0
				endif
				ieq_cell%P_H2 = P_H2(j)
				ieq_cell%k_LW = k_lw_diss(j)   ! 0 without a LW band flux
				ieq_cell%T_K  = T_K(j)
				ieq_cell%ntot = n_tot(j)        ! M for the 3-body rates
				! Compute the molecular rate coefficients that are invariant
				! across this cell's Newton solve (they depend only on T and
				! n_tot); the residual then reads them, like set_metal_coeffs.
				! n_tot, not n_in_dim: R13 and R15 fold the third body M into
				! their coefficient, so they need the same particle density
				! the R12 term of the residual reads (see the n_tot comment
				! above).
				call set_mol_coeffs(T_K(j), n_tot(j))
				! Turnover scale of each balance row, from the coefficients
				! just built and this cell's densities; the residual divides
				! by it so the helium and molecular blocks reach hybrd1 with
				! the same weight.
				call set_mol_turnover_rates(ne(j))
			endif

			! Each element's metal coefficients are handed to
			! ion_system_HeH_metals via set_metal_coeffs; the charge-
			! exchange rate coefficients are stored for this cell by
			! cx_set_cell (used inside the residual by cx_add_to_fvec).
			if (thereis_metals) then
				call cx_set_cell(T_K(j))

				! Build the metal coefficients for each element in canonical order
				! from the 2D rate arrays (col i0 = neutral, i0+1 = +,
				! i0+2 = ++) and store them for the residual.
				do im = 1,n_melem
					i0  = melem_i0(im)
					top = melem_top(im)
					meg_top(im)  = top
					meg_ntot(im) = nm_tot(j,im)
					meg_g0(im)   = P_m(j,i0)
					meg_b0(im)   = aion_m(j,i0)
					meg_a1(im)   = rec_m(j,i0+1)
					if (top .ge. 2) then
						meg_g1(im) = P_m(j,i0+1)
						meg_b1(im) = aion_m(j,i0+1)
						meg_a2(im) = rec_m(j,i0+2)
					else
						! No second ionization stage for this element.
						meg_g1(im) = 0.0d0
						meg_b1(im) = 0.0d0
						meg_a2(im) = 0.0d0
					endif
				enddo
				call set_metal_coeffs(n_melem, meg_ntot, meg_g0, meg_g1, &
				                    meg_b0, meg_b1, meg_a1, meg_a2,     &
				                    meg_top)

				! Metal rows of the molecular system get their turnover
				! scale too, after set_mol_turnover_rates has reset the
				! array and set_metal_coeffs has filled met_*.
				if (thereis_mol)                                     &
					call set_mol_metal_turnover_rates(ne(j))
			endif

			! Initial guess
			if (count .eq. 0) then
				if (j .eq. N+Ng) then
					sys_x(1) = 1.0
					sys_x(2) = 1.0
					sys_x(3) = 1.0
					if (thereis_HeITR) sys_x(4) = 0.01
					if (thereis_mol) then
						! start fully ionized aloft; molecules negligible
						sys_x(4:7) = 1.0d-10
						if (thereis_HeITR) sys_x(8) = 0.01
					endif
					if (thereis_metals) then
						! Each metal starts fully singly ionized.
						do im = 1,n_melem
							sys_x(mbase + 2*(im-1))   = 1.0  ! X+  frac
							sys_x(mbase+1 + 2*(im-1)) = 0.0  ! X++ frac
						enddo
					endif
				else
					sys_x(1) = nhii(j+1)/nh(j+1)
					sys_x(2) = nheii(j+1)/nhe(j+1)
					sys_x(3) = nheiii(j+1)/nhe(j+1)
					if (thereis_HeITR) &
						sys_x(4) = nheiTR(j+1)/nhe(j+1)
					if (thereis_mol) then
						sys_x(4) = max(2.0d0*nmol_eq(j+1,1)/nh(j+1), 1.0d-10)
						sys_x(5) = max(2.0d0*nmol_eq(j+1,2)/nh(j+1), 1.0d-12)
						sys_x(6) = max(3.0d0*nmol_eq(j+1,3)/nh(j+1), 1.0d-12)
						sys_x(7) = max(nmol_eq(j+1,4)/nh(j+1), 1.0d-14)
						if (thereis_HeITR) sys_x(8) = nheiTR(j+1)/nhe(j+1)
					endif
					if (thereis_metals) then
						do im = 1,n_melem
							i0 = melem_i0(im)
							sys_x(mbase+2*(im-1)) = nm(j+1,i0+1)        &
							                  /max(nm_tot(j+1,im),1.0d-30)
							if (melem_top(im) .ge. 2) then
								sys_x(mbase+1+2*(im-1)) = nm(j+1,i0+2)        &
								                  /max(nm_tot(j+1,im),1.0d-30)
							else
								sys_x(mbase+1+2*(im-1)) = 0.0d0
							endif
						enddo
					endif
				endif
			else
				sys_x(1) = nhii(j)/nh(j)
				sys_x(2) = nheii(j)/nhe(j)
				sys_x(3) = nheiii(j)/nhe(j)
				if (thereis_HeITR) sys_x(4) = nheiTR(j)/nhe(j)
				if (thereis_mol) then
					sys_x(4) = 2.0d0*nmol_eq(j,1)/nh(j)
					sys_x(5) = 2.0d0*nmol_eq(j,2)/nh(j)
					sys_x(6) = 3.0d0*nmol_eq(j,3)/nh(j)
					sys_x(7) = nmol_eq(j,4)/nh(j)
					if (thereis_HeITR) sys_x(8) = nheiTR(j)/nhe(j)
				endif
				if (thereis_metals) then
					do im = 1,n_melem
						i0 = melem_i0(im)
						sys_x(mbase+2*(im-1)) = nm(j,i0+1)              &
						                  /max(nm_tot(j,im),1.0d-30)
						if (melem_top(im) .ge. 2) then
							sys_x(mbase+1+2*(im-1)) = nm(j,i0+2)              &
							                  /max(nm_tot(j,im),1.0d-30)
						else
							sys_x(mbase+1+2*(im-1)) = 0.0d0
						endif
					enddo
				endif

				! The guess just built is this cell's previous state. If that
				! state is not physical it must not seed the solve again: an
				! unphysical root fed back as the next step's guess reproduces
				! itself indefinitely (the self-sticking negative H II / Fe I
				! state at the base). Start from the local ionization balance
				! instead, which lies inside the simplex by construction. The
				! molecular network keeps its own chemical-equilibrium retry.
				if (.not.thereis_mol) then
					if (.not.ionization_fractions_physical(sys_x,N_eq,mbase)) then
						call ionization_balance_at_fixed_ne(sys_x,N_eq,   &
						                                    mbase,ne(j))
						n_ieq_reseed = n_ieq_reseed + 1
					endif
				endif
			endif

		 	! Analytic-Jacobian Newton (Task 2); hybrd1 fallback inside solve_ieq.
			! The He metastable-triplet systems keep the MINPACK solve (no
			! analytic Jacobian written for the triplet kinetics).
			if (thereis_mol) then
				! Molecular network, with the trace metals solved in the same
				! system when they are present (System_HeH_mol_metals: the two
				! blocks share the free electron density, which the metals
				! dominate in the shielded molecular base).
				!
				! hybrd1 is bistable in its initial guess, because the cell
				! itself is: a molecular basin (the dense, optically-thick
				! base, where the true root is strongly molecular) and an
				! atomic basin (the wind above the H2 -> H front). From a
				! zero-molecular warm start the solve fails at the base; at
				! the front cells themselves neither warm start need land in
				! the right basin. One starting point in each basin is
				! therefore tried in turn:
				!   1  the guess built above (previous state / neighbour),
				!   2  molecular basin -- the H2 dissociation equilibrium of
				!      the local (p, T) with every element in its own
				!      ionization balance at the incoming n_e,
				!   3  atomic basin -- the same balance with no H2,
				! and a root is kept only if it is physical (every stage
				! fraction >= 0, each element's tracked stages summing to at
				! most its nuclei).
				!
				! A cell that reaches the solver tolerance on a physical root
				! takes that root and stops. A cell that does not keeps the
				! best root the three attempts produced, and when NONE of them
				! was admissible it keeps the one lying closest to the simplex
				! and clamps it back onto the element budget. Every outcome is
				! therefore a root of the network computed at THIS cell's
				! state; none of them is the cell's previous composition,
				! which would make the equilibrium -- and the steady residual
				! built on it -- a function of the sequence of evaluations
				! rather than of the state being evaluated.
				best_rank = 0
				viol_best = 0.0d0
				! Gas pressure p = (n_tot + n_e) kB T [bar], the same ideal-gas
				! law the EOS uses (comp_p_from_T).
				pbar_loc  = (n_tot(j) + ne(j))*kb_erg*T_K(j)/1.0d6
				do iatt = 1,3
					if (iatt .eq. 2) then
						call dissociation_ionization_balance_at_fixed_ne(  &
						                 sys_x,N_eq,mbase,ne(j),pbar_loc,T_K(j))
						n_ieq_retry = n_ieq_retry + 1
					else if (iatt .eq. 3) then
						call ionization_balance_at_fixed_ne(sys_x,N_eq,   &
						                                    mbase,ne(j))
						n_ieq_retry = n_ieq_retry + 1
					endif

					if (thereis_metals) then
						! Metals appended above the molecular unknowns; point
						! charge exchange at their rows for this solve.
						cx_metal_base = mbase
						call hybrd1(ion_system_HeH_mol_metals,N_eq,sys_x,  &
						            sys_sol,tol,info,wa,lwa,params)
						cx_metal_base = 4
					else
						call hybrd1(ion_system_HeH_mol,N_eq,sys_x,sys_sol, &
						            tol,info,wa,lwa,params)
					endif
					conv_ieq = (info .eq. 1)
					n_mol_info(max(0,min(5,info))) =                     &
						n_mol_info(max(0,min(5,info))) + 1

					phys_ieq = ionization_fractions_physical(sys_x,N_eq,mbase)
					if (iatt .eq. 1 .and. .not.phys_ieq)                &
						n_ieq_unphys = n_ieq_unphys + 1
					ok_rank = 0
					if (phys_ieq) ok_rank = 1
					if (phys_ieq .and. conv_ieq) ok_rank = 2
					! Strict improvement only, so attempt 1 wins any tie.
					if (ok_rank .gt. best_rank) then
						best_rank = ok_rank
						x_root_best(1:N_eq) = sys_x(1:N_eq)
					else if (best_rank .eq. 0) then
						! Nothing admissible yet: remember the root that lies
						! closest to the simplex, to be clamped onto it below.
						viol = element_budget_violation(sys_x,N_eq,mbase)
						if (iatt .eq. 1 .or. viol .lt. viol_best) then
							viol_best = viol
							x_root_best(1:N_eq) = sys_x(1:N_eq)
						endif
					endif
					if (ok_rank .eq. 2) exit
				enddo

				if (best_rank .lt. 2) then
					! No attempt met the solver tolerance; keep the best root
					! found. If none of them was admissible, keep the one that
					! lies closest to the simplex and clamp it onto the element
					! budget: the violations are the solver's own resolution of a
					! root sitting on a face (a fully dissociated or fully
					! ionized species), so the clamped state is a root of the
					! network to the accuracy the solve reached. What it is NOT
					! is the cell's previous composition -- inheriting that would
					! make the equilibrium, and the steady residual built on it,
					! a function of the sequence of evaluations rather than of
					! the state being evaluated.
					sys_x(1:N_eq) = x_root_best(1:N_eq)
					if (best_rank .eq. 0) then
						call clamp_fractions_to_element_budget(sys_x,N_eq,mbase)
						n_mol_clamped = n_mol_clamped + 1
					endif
				endif
			else
				! Atomic H/He (+ He 2^3S) (+ metals). A vanishing residual is
				! not sufficient: the root must also be physical, i.e. every
				! stage fraction non-negative and the ionized stages of each
				! element summing to at most its nucleus total. The cell is
				! solved from up to three starting points and the best
				! admissible root is kept:
				!   1  the guess built above (previous state / neighbour),
				!   2  the uncoupled ionization balance at the incoming n_e,
				!   3  the optically thick limit, all nuclei neutral.
				! Attempt 1 is accepted as it stands whenever it converges to
				! a physical root, so a healthy cell takes exactly the same
				! solver path as before.
				best_rank = 0
				do iatt = 1,3
					if (iatt .eq. 2) then
						call ionization_balance_at_fixed_ne(sys_x,N_eq,   &
						                                    mbase,ne(j))
						n_ieq_retry = n_ieq_retry + 1
					else if (iatt .eq. 3) then
						sys_x(1:N_eq) = 0.0d0
					endif

					if (thereis_HeITR .and. thereis_metals) then
						! Merged He-triplet + metals: triplet at sys_x(4),
						! metals at sys_x(5..). The cell-by-cell metal
						! coefficients (set_metal_coeffs) and charge-exchange
						! rates (cx_set_cell) were already loaded above in the
						! thereis_metals block; here we only point charge
						! exchange at the shifted metal rows for the merged
						! solve.
						cx_metal_base = 5
						call hybrd1(ion_system_HeH_TR_metals,N_eq,sys_x,  &
						            sys_sol,tol,info,wa,lwa,params)
						cx_metal_base = 4
						conv_ieq = (info .eq. 1)
					else if (thereis_HeITR) then
						call hybrd1(ion_system_HeH_TR,N_eq,sys_x,sys_sol, &
						            tol,info,wa,lwa,params)
						conv_ieq = (info .eq. 1)
					else if (thereis_metals) then
						call solve_ieq(ion_system_HeH_metals,             &
						            jac_system_HeH_metals,                &
						            N_eq,sys_x,params,tol,wa,lwa,usednt,  &
						            info_ieq)
						conv_ieq = (info_ieq .eq. 1)
					else
						call solve_ieq(ion_system_HeH,jac_system_HeH,     &
						            N_eq,sys_x,params,tol,wa,lwa,usednt,  &
						            info_ieq)
						conv_ieq = (info_ieq .eq. 1)
					endif

					phys_ieq = ionization_fractions_physical(sys_x,N_eq,mbase)
					if (iatt .eq. 1 .and. .not.phys_ieq)                &
						n_ieq_unphys = n_ieq_unphys + 1
					ok_rank = 0
					if (phys_ieq) ok_rank = 1
					if (phys_ieq .and. conv_ieq) ok_rank = 2
					! Strict improvement only, so attempt 1 wins any tie and
					! an already healthy cell is untouched.
					if (ok_rank .gt. best_rank) then
						best_rank = ok_rank
						x_root_best(1:N_eq) = sys_x(1:N_eq)
					endif
					if (ok_rank .eq. 2) exit
				enddo

				if (best_rank .lt. 2) then
					if (best_rank .eq. 1) then
						! No attempt met the solver tolerance; keep the best
						! physical root found.
						sys_x(1:N_eq) = x_root_best(1:N_eq)
					else
						! Every attempt left the physical simplex. Rather than
						! propagate negative densities, hand back the uncoupled
						! ionization balance, which is admissible by
						! construction, and report the cell.
						call ionization_balance_at_fixed_ne(sys_x,N_eq,   &
						                                    mbase,ne(j))
						n_ieq_fail = n_ieq_fail + 1
					endif
				endif
			endif

			! Extract solution profiles
			if (thereis_mol) then
				! guard tiny negatives from the NL solve
				sys_x(1:N_eq) = max(sys_x(1:N_eq), 0.0d0)
				nhii(j)   = nh(j)*sys_x(1)
				nmol_eq(j,1) = 0.5d0*sys_x(4)*nh(j)
				nmol_eq(j,2) = 0.5d0*sys_x(5)*nh(j)
				nmol_eq(j,3) = sys_x(6)*nh(j)/3.0d0
				nmol_eq(j,4) = sys_x(7)*nh(j)
				nhi(j)    = nh(j)*max(1.0d0 - sys_x(1) - sys_x(4)      &
				              - sys_x(5) - sys_x(6) - sys_x(7), 0.0d0)
				nheii(j)  = nhe(j)*sys_x(2)
				nheiii(j) = nhe(j)*sys_x(3)
				nhei(j)   = max(nhe(j)*(1.0d0 - sys_x(2) - sys_x(3))   &
				              - nmol_eq(j,4), 0.0d0)
				if (thereis_HeITR) nheiTR(j) = nhe(j)*sys_x(8)
			else
			nhi(j)    = nh(j)*(1.0 - sys_x(1))
			nhii(j)   = nh(j)*sys_x(1)
			nhei(j)   = nhe(j)*(1.0 - sys_x(2) - sys_x(3))
			nheii(j)  = nhe(j)*sys_x(2)
			nheiii(j) = nhe(j)*sys_x(3)
			if (thereis_HeITR) nheiTR(j) = nhe(j)*sys_x(4)
			endif
			if (thereis_metals) then
				do im = 1,n_melem
					i0 = melem_i0(im)
					nm(j,i0+1) = nm_tot(j,im)*sys_x(mbase+2*(im-1))
					if (melem_top(im) .ge. 2) then
						nm(j,i0+2) = nm_tot(j,im)*sys_x(mbase+1+2*(im-1))
						nm(j,i0)   = nm_tot(j,im)                            &
						           *(1.0 - sys_x(mbase+2*(im-1)) - sys_x(mbase+1+2*(im-1)))
					else
						! Two-stage element: neutral = total - singly ionized.
						nm(j,i0)   = nm_tot(j,im)*(1.0 - sys_x(mbase+2*(im-1)))
					endif
				enddo
			endif

		enddo
		!$omp end parallel do

		! One summary line per sweep when a molecular cell's roots were all
		! outside the physical simplex, so the closest one was clamped onto the
		! element budget.
		if (thereis_mol .and. n_mol_clamped .gt. 0) then
			write(*,'(A,I0,A)') ' (ioniz_eq) WARNING: every molecular '//     &
				'equilibrium root left the physical simplex at ',              &
				n_mol_clamped, ' cell(s), clamped onto the element budget'
		endif

		! One summary line per sweep whenever a state left the physical
		! simplex: how many cells had their stored state rejected as a
		! starting point, how many first roots were outside the simplex, and
		! how many ended on the ionization balance because no starting point
		! produced an admissible root. A cell that is merely retried because
		! the solver did not reach its tolerance is routine and stays silent
		! here; it is counted in the run-wide totals reported at the end.
		if (n_ieq_reseed + n_ieq_unphys + n_ieq_fail .gt. 0) then
			write(*,'(A,I0,A,I0,A,I0,A,I0,A)')                             &
				' (ioniz_eq) step ', count,                                &
				': ionization roots - ', n_ieq_reseed,                     &
				' stored state(s) rejected, ', n_ieq_unphys,               &
				' root(s) outside the simplex, no admissible root at ',    &
				n_ieq_fail, ' cell(s)'
		endif

		ieq_n_reseed = ieq_n_reseed + n_ieq_reseed
		ieq_n_retry  = ieq_n_retry  + n_ieq_retry
		ieq_n_unphys = ieq_n_unphys + n_ieq_unphys
		ieq_n_noroot = ieq_n_noroot + n_ieq_fail
		ieq_n_mol_clamped = ieq_n_mol_clamped + n_mol_clamped
		ieq_n_mol_info(:) = ieq_n_mol_info(:) + n_mol_info(:)

	endif

	
	! Density with atomic numbers (nm adds the metal mass under the
	! eos_metals policy)
   if (thereis_mol) then
      call calc_rho(nhi,nhii,nhei,nheii,nheiii,n_io,nm,nmol_eq)
   else
      call calc_rho(nhi,nhii,nhei,nheii,nheiii,n_io,nm)
   endif

   ! Abundancies profiles
   f_sp_io(:,1) = nhi/n_io
   f_sp_io(:,2) = nhii/n_io
   f_sp_io(:,3) = nhei/n_io
   f_sp_io(:,4) = nheii/n_io
   f_sp_io(:,5) = nheiii/n_io
   f_sp_io(:,6) = nheiTR/n_io
   ! molecular abundances
   if (thereis_mol) then
      f_sp_io(:,isp_H2)   = nmol_eq(:,1)/n_io
      f_sp_io(:,isp_H2p)  = nmol_eq(:,2)/n_io
      f_sp_io(:,isp_H3p)  = nmol_eq(:,3)/n_io
      f_sp_io(:,isp_HeHp) = nmol_eq(:,4)/n_io
   else
      f_sp_io(:,isp_H2:isp_HeHp) = 0.0d0
   endif

   ! Metal abundances in canonical order (col mion_fsp(im) of f_sp_io).
   do im = 1,n_mion
      f_sp_io(:,mion_fsp(im)) = nm(:,im)/n_io
   enddo

   ! Adimensional number density profile
	n_io = n_io/n0

   ! Adimensional heating and cooling rates
   heat_out = heat/q0
   cool_out = cool/q0
      
   ! Adjust value of pressure boundary condition (the base electron
   ! density in units of n0; with eos_metals the metal electrons are
   ! included, consistently with calc_ne). Molecular ions are deliberately
   ! omitted as trace electron donors: the base is nearly neutral, so the
   ! molecular-ion electrons are negligible in dp_bc.
   dp_bc = (nhii(1-Ng) + nheii(1-Ng) + 2.0*nheiii(1-Ng))/n0
   if (eos_include_metals .and. thereis_metals) then
      do im = 1,n_mion
         if (mion_stage(im) .gt. 0)                                     &
            dp_bc = dp_bc + dble(mion_stage(im))*nm(1-Ng,im)/n0
      enddo
   endif

	! End of subroutine
	end subroutine ioniz_eq

	!----------------------------------!

	logical function ionization_fractions_physical(x,n,mbase) result(ok)
	! Is a root of the equilibrium system a physically admissible state?
	!
	! Every unknown of these systems is the fraction of one element's nuclei
	! found in one ionization stage, so a physical state has
	!   (i)  every fraction >= 0, and
	!   (ii) for each element, the tracked ionized stages summing to <= 1,
	!        the neutral stage being the remainder 1 - sum.
	! The residuals are polynomials in the fractions and do admit roots
	! outside this simplex; such a root gives a negative neutral density and
	! is not a solution of the physical problem, however small the residual.
	! Downstream it poisons the photoionization integrals (a negative
	! absorber column gives negative photoheating) and the element budget.
	!
	! The comparisons carry a small tolerance so that a root sitting on a
	! face of the simplex (a fully neutral or fully ionized element) is not
	! rejected for round-off; it is far below the violations this test is
	! meant to catch (fractions of ~1e-6 to ~1e-1 outside the simplex).
	!
	! Widening it to the solver's own xtol = sqrt(machine epsilon) was tried
	! and measured: it lets a cell stop on a root that sits just outside a
	! face instead of retrying from another starting point, and the retried
	! root is the better one. On examples/15 (HD 209458 b, molecular, cold)
	! the wider band turns info = 0 at ||R|| = 6.105e-6 in 89 outer
	! iterations into info = 2 at 1.112e-4 in 133, so the band stays at
	! round-off.
	!
	! Layout: x(1) = H II / H, x(2) = He II / He, x(3) = He III / He,
	! x(4) = He 2^3S / He when the triplet is tracked (x(8) in the molecular
	! layout, where x(4..7) are the H nuclei bound in H2, H2+, H3+, HeH+),
	! and the metal stages X+ / X++ from x(mbase) upwards, two per element.

	integer, intent(in) :: n, mbase
	real*8,  intent(in) :: x(n)
	real*8, parameter   :: ftol = 1.0d-10
	integer :: im, ix
	real*8  :: s

	ok = .false.

	! Hydrogen nuclei: H II, plus the H bound in molecules where tracked.
	if (x(1) .lt. -ftol) return
	s = x(1)
	if (thereis_mol) then
		do ix = 4,7
			if (x(ix) .lt. -ftol) return
			s = s + x(ix)
		enddo
	endif
	if (s .gt. 1.0d0 + ftol) return

	! Helium nuclei: He II, He III, the 2^3S metastable (which the systems
	! carry as a separate level inside the neutral stage) and, in the
	! molecular layout, the He nucleus bound in HeH+. The molecular systems
	! close the neutral He as n_He(1 - x2 - x3) - n_HeH+, so a root leaving
	! that negative is inadmissible even when x2 + x3 <= 1; x(7) is HeH+ per
	! H nucleus, hence the n_H/n_He conversion from this cell's state.
	if (thereis_He) then
		if (x(2) .lt. -ftol) return
		if (x(3) .lt. -ftol) return
		s = x(2) + x(3)
		if (thereis_HeITR) then
			ix = 4
			if (thereis_mol) ix = 8
			if (x(ix) .lt. -ftol) return
			s = s + x(ix)
		endif
		if (thereis_mol .and. ieq_cell%nhe .gt. 0.0d0)                &
			s = s + x(7)*ieq_cell%nh/ieq_cell%nhe
		if (s .gt. 1.0d0 + ftol) return
	endif

	! Each metal element separately: X+ (+ X++ for the three-stage elements).
	if (thereis_metals) then
		do im = 1,n_melem
			ix = mbase + 2*(im-1)
			if (x(ix) .lt. -ftol) return
			s = x(ix)
			if (melem_top(im) .ge. 2) then
				if (x(ix+1) .lt. -ftol) return
				s = s + x(ix+1)
			endif
			if (s .gt. 1.0d0 + ftol) return
		enddo
	endif

	ok = .true.

	end function ionization_fractions_physical

	!----------------------------------!

	subroutine ionization_balance_at_fixed_ne(x,n,mbase,n_e)
	! Stage fractions of each element in its OWN ionization balance at a
	! given electron density: photoionization plus electron-impact ionization
	! against radiative recombination,
	!
	!    n_k (gamma_k + beta_k n_e) = alpha_{k+1} n_e n_{k+1},
	!
	! with the couplings between elements (charge exchange, and the
	! dependence of n_e on the unknowns themselves) dropped and n_e taken
	! from the incoming state. Writing u_k = gamma_k + beta_k n_e for the
	! rate out of stage k and d_k = alpha_k n_e for the rate back into stage
	! k-1, the three-stage solution is
	!
	!    (n_0, n_1, n_2) proportional to (d_1 d_2, u_0 d_2, u_0 u_1),
	!
	! non-negative and normalized to one, so the result always lies inside
	! the physical simplex whatever state it is asked to replace. It is used
	! as a starting point for the coupled solve, and as the state of last
	! resort for a cell where no starting point produced an admissible root.
	!
	! The He 2^3S fraction is the steady state of the metastable level at
	! those populations: recombination and collisional excitation feed it,
	! radiative decay, collisional de-excitation, photoionization and Penning
	! ionization on H0 drain it. It is capped by the neutral He fraction.
	!
	! Rates come from the same cell state the residuals read (ieq_cell and
	! the met_* coefficients of System_HeH_metals), so this is the incoming
	! cell's own physics, not a generic guess. Molecules are not part of this
	! balance: it is the molecule-free limit of the state, so in the molecular
	! layout the H2/H2+/H3+/HeH+ fractions x(4..7) are left at zero and the
	! metastable sits at x(8).

	integer, intent(in)  :: n, mbase
	real*8,  intent(in)  :: n_e
	real*8,  intent(out) :: x(n)
	integer :: itr
	real*8  :: u0,u1,d1,d2,w1,w2,s
	real*8  :: xneu, n_hi_loc

	x(1:n) = 0.0d0

	! Hydrogen
	u0 = ieq_cell%P_HI + ieq_cell%a_ion_HI*n_e
	d1 = ieq_cell%rchiiB*n_e
	s  = u0 + d1
	if (s .gt. 0.0d0) x(1) = u0/s

	if (thereis_He) then
		! Helium. Both He+ recombination channels return to neutral He: the
		! metastable capture rcheiTR is a branch of the recombination, so it
		! adds to the rate back into He I when the triplet is tracked.
		u0 = ieq_cell%P_HeI + ieq_cell%a_ion_HeI*n_e
		d1 = ieq_cell%rcheiiB*n_e
		if (thereis_HeITR) d1 = d1 + ieq_cell%rcheiTR*n_e
		u1 = ieq_cell%P_HeII + ieq_cell%a_ion_HeII*n_e
		d2 = ieq_cell%rcheiiiB*n_e
		w1 = u0*d2
		w2 = u0*u1
		s  = d1*d2 + w1 + w2
		if (s .gt. 0.0d0) then
			x(2) = w1/s
			x(3) = w2/s
		endif

		if (thereis_HeITR) then
			itr = 4
			if (thereis_mol) itr = 8
			xneu     = max(1.0d0 - x(2) - x(3), 0.0d0)
			n_hi_loc = max(1.0d0 - x(1), 0.0d0)*ieq_cell%nh
			s = ieq_cell%P_HeITR + ieq_cell%A31                          &
			  + n_hi_loc*ieq_cell%Q31                                    &
			  + (ieq_cell%q31a + ieq_cell%q31b                           &
			     + ieq_cell%a_ion_HeITR)*n_e
			if (s .gt. 0.0d0) x(itr) = min(                              &
			      n_e*(x(2)*ieq_cell%rcheiTR + xneu*ieq_cell%q13)/s, xneu)
		endif
	endif

	if (thereis_metals)                                                  &
		call metal_ionization_balance_at_fixed_ne(x,n,mbase,n_e)

	end subroutine ionization_balance_at_fixed_ne

	!----------------------------------!

	real*8 function element_budget_violation(x,n,mbase) result(viol)
	! How far a root of an equilibrium system lies outside the physically
	! allowed states: the largest of the negative stage fractions and of the
	! amounts by which one element's tracked stages exceed its nuclei, both
	! measured as fractions of the element. Zero for an admissible state.
	! Used to pick, among roots that are all inadmissible, the one closest to
	! a state (see clamp_fractions_to_element_budget). Same layout as
	! ionization_fractions_physical.

	integer, intent(in) :: n, mbase
	real*8,  intent(in) :: x(n)
	integer :: im, ix
	real*8  :: s

	viol = 0.0d0
	do ix = 1,n
		viol = max(viol, -x(ix))
	enddo

	s = x(1)
	if (thereis_mol) s = s + x(4) + x(5) + x(6) + x(7)
	viol = max(viol, s - 1.0d0)

	if (thereis_He) then
		s = x(2) + x(3)
		if (thereis_HeITR) then
			ix = 4
			if (thereis_mol) ix = 8
			s = s + x(ix)
		endif
		! The He nucleus bound in HeH+ (x(7) is per H nucleus), as in
		! ionization_fractions_physical.
		if (thereis_mol .and. ieq_cell%nhe .gt. 0.0d0)                &
			s = s + x(7)*ieq_cell%nh/ieq_cell%nhe
		viol = max(viol, s - 1.0d0)
	endif

	if (thereis_metals) then
		do im = 1,n_melem
			ix = mbase + 2*(im-1)
			s  = x(ix)
			if (melem_top(im) .ge. 2) s = s + x(ix+1)
			viol = max(viol, s - 1.0d0)
		enddo
	endif

	end function element_budget_violation

	!----------------------------------!

	subroutine clamp_fractions_to_element_budget(x,n,mbase)
	! Move a root that lies just outside the physically allowed states onto the
	! nearest allowed one: no negative populations, and no element with more
	! nuclei in its tracked stages than it has. Negative fractions are set to
	! zero and, where an element's tracked stages still sum above its nuclei,
	! they are scaled down to sum exactly to them, leaving the neutral stage
	! empty.
	!
	! This is the treatment of a root the solver has resolved to a face of the
	! allowed region -- a fully dissociated H2, a fully neutral or fully
	! ionized element -- and delivered as a small number of either sign. The
	! clamped state is that root to the accuracy the solve reached, so it
	! remains a solution of the equilibrium at THIS cell's state; what it must
	! not be replaced by is a composition carried over from an earlier
	! evaluation.

	integer, intent(in)    :: n, mbase
	real*8,  intent(inout) :: x(n)
	integer :: im, ix
	real*8  :: s, hehp_he

	do ix = 1,n
		if (x(ix) .lt. 0.0d0) x(ix) = 0.0d0
	enddo

	! Hydrogen nuclei: H II and the H bound in molecules where tracked.
	s = x(1)
	if (thereis_mol) s = s + x(4) + x(5) + x(6) + x(7)
	if (s .gt. 1.0d0) then
		x(1) = x(1)/s
		if (thereis_mol) then
			do ix = 4,7
				x(ix) = x(ix)/s
			enddo
		endif
	endif

	! Helium nuclei: He II, He III, the 2^3S metastable and, in the molecular
	! layout, the He nucleus bound in HeH+. x(7) is HeH+ per H nucleus and
	! belongs to BOTH element budgets; it has just been capped against the H
	! nuclei, so it is capped against the He nuclei here as well (which can
	! only relax the H budget), and the free-He stages are scaled into what
	! is left. hehp_he is exactly zero without molecules, where the scaling
	! reduces to the plain x/s of every non-molecular run.
	if (thereis_He) then
		hehp_he = 0.0d0
		if (thereis_mol .and. ieq_cell%nhe .gt. 0.0d0) then
			hehp_he = x(7)*ieq_cell%nh/ieq_cell%nhe
			if (hehp_he .gt. 1.0d0) then
				x(7)    = x(7)/hehp_he
				hehp_he = 1.0d0
			endif
		endif
		s = x(2) + x(3)
		ix = 0
		if (thereis_HeITR) then
			ix = 4
			if (thereis_mol) ix = 8
			s = s + x(ix)
		endif
		if (s .gt. 1.0d0 - hehp_he) then
			x(2) = x(2)*(1.0d0 - hehp_he)/s
			x(3) = x(3)*(1.0d0 - hehp_he)/s
			if (ix .gt. 0) x(ix) = x(ix)*(1.0d0 - hehp_he)/s
		endif
	endif

	! Each metal element separately.
	if (thereis_metals) then
		do im = 1,n_melem
			ix = mbase + 2*(im-1)
			s  = x(ix)
			if (melem_top(im) .ge. 2) s = s + x(ix+1)
			if (s .gt. 1.0d0) then
				x(ix) = x(ix)/s
				if (melem_top(im) .ge. 2) x(ix+1) = x(ix+1)/s
			endif
		enddo
	endif

	end subroutine clamp_fractions_to_element_budget

	!----------------------------------!

	subroutine dissociation_ionization_balance_at_fixed_ne(x,n,mbase,n_e,   &
	                                                       p_bar,T_gas)
	! H2 dissociation equilibrium at the local (p, T), with every element in
	! its own ionization balance at a given electron density: the limit of the
	! molecular network in which the couplings between the molecular ions and
	! the ionization of H and He are dropped.
	!
	!   - The H nuclei are partitioned between H2 and atomic H by the
	!     chemical-equilibrium mixing ratio q_H2(p, T) of Koskinen et al.
	!     (2022) Eq. 11, the same fit the molecular base boundary condition
	!     uses; per H nucleus the fraction bound in H2 is
	!     2 q (1 + He/H)/(1 + q) (lower_column::mu_mixture).
	!   - The atomic remainder is ionized in the H ionization balance of
	!     ionization_balance_at_fixed_ne, which also sets He, the 2^3S
	!     metastable and the metal stages.
	!   - H2+, H3+ and HeH+ are left at zero: each is a trace intermediate
	!     whose abundance is set by the very couplings this limit drops, and
	!     each holds far fewer H nuclei than H2 wherever the gas is molecular.
	!
	! Every fraction returned is non-negative and each element's stages sum to
	! at most its nuclei, so the starting point is a physically allowed state,
	! and it depends only on the local (p, T, n_e) and on the rate
	! coefficients of the cell -- not on the cell's previous composition.
	!
	! Molecular layout only: x(4) is the H-nucleus fraction bound in H2 there,
	! and the He 2^3S metastable sits at x(8).

	integer, intent(in)  :: n, mbase
	real*8,  intent(in)  :: n_e, p_bar, T_gas
	real*8,  intent(out) :: x(n)
	real*8 :: qh2, x_h2

	call ionization_balance_at_fixed_ne(x,n,mbase,n_e)

	qh2  = q_h2_equilibrium(p_bar, T_gas)
	x_h2 = 2.0d0*qh2*(1.0d0 + HeH)/(1.0d0 + qh2)
	if (x_h2 .gt. 1.0d0) x_h2 = 1.0d0
	x(4) = x_h2
	! Only the H nuclei left atomic are available to ionize.
	x(1) = x(1)*(1.0d0 - x_h2)

	end subroutine dissociation_ionization_balance_at_fixed_ne

	!----------------------------------!

	subroutine metal_ionization_balance_at_fixed_ne(x,n,mbase,n_e)
	! Metal stage fractions in each element's OWN ionization balance at a
	! given electron density, element by element in canonical order (the
	! metal part of ionization_balance_at_fixed_ne above; see there for the
	! balance itself). Written separately because the molecular retry needs
	! exactly this part: its H/He/molecular unknowns come from the
	! chemical-equilibrium H2 fit, its metal unknowns from here. Only
	! x(mbase..) is touched, so the caller's other seeds are preserved. An
	! element that is absent keeps its stages at zero, as its residual rows do.

	use System_HeH_metals, only: met_nelem, met_ntot, met_g0, met_g1,     &
	                             met_b0, met_b1, met_a1, met_a2, met_top

	integer, intent(in)    :: n, mbase
	real*8,  intent(in)    :: n_e
	real*8,  intent(inout) :: x(n)
	integer :: im, ix, top
	real*8  :: u0,u1,d1,d2,w1,w2,s

	do im = 1,met_nelem
		ix = mbase + 2*(im-1)
		x(ix)   = 0.0d0
		x(ix+1) = 0.0d0
		if (met_ntot(im) .le. 1.0d-30) cycle
		top = met_top(im)
		u0  = met_g0(im) + met_b0(im)*n_e
		d1  = met_a1(im)*n_e
		if (top .ge. 2) then
			u1 = met_g1(im) + met_b1(im)*n_e
			d2 = met_a2(im)*n_e
			w1 = u0*d2
			w2 = u0*u1
			s  = d1*d2 + w1 + w2
			if (s .gt. 0.0d0) then
				x(ix)   = w1/s
				x(ix+1) = w2/s
			endif
		else
			s = d1 + u0
			if (s .gt. 0.0d0) x(ix) = u0/s
		endif
	enddo

	end subroutine metal_ionization_balance_at_fixed_ne

	! End of module
	end module ionization_equilibrium
