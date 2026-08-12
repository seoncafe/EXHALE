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
	use h3p_cooling,  only: h3p_cooling_rate
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
	real*8, dimension(1-Ng:N+Ng,4), save :: nmol_eq = 0.0d0

	! Run-wide totals of the atomic ionization root validation, reported once
	! at the end of the run (EXHALE_main) next to the Newton usage counters:
	! stored states rejected as a starting point, cell solves that needed a
	! second or third starting point, first roots that lay outside the
	! physical simplex, and cells where no starting point produced an
	! admissible root. All zero for a run that never leaves the simplex.
	integer, save :: ieq_n_reseed = 0, ieq_n_retry  = 0
	integer, save :: ieq_n_unphys = 0, ieq_n_noroot = 0

	contains 
	
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
	                                 ne,n_in_dim
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
   ! Molecular solve (thereis_mol): count of cells whose equilibrium solve
   ! failed (kept previous state), and the physically-informed guess scalars
   ! for the retry (chemical-equilibrium H2 fit at the local p, T).
   integer :: n_molfail
   real*8  :: pbar_loc, qh2_loc, x2_loc

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

	! Cell-by-cell pressure-broadening factor for the opacity ('P' model).
	! opacity_pT_factor returns 1.0 for all other models, so opa_pf=1
	! and the column densities are unchanged (bit-identical).
	do j = 1-Ng,N+Ng
		opa_pf(j) = opacity_pT_factor((nh(j)+nhe(j)+ne(j))*kb_erg*T_K(j))
	enddo

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
    ! balances the ionization against (zero for an atomic run).
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

	! optically-thin H3+ infrared cooling (Miller+2013 fits with
	! the Table-6 non-LTE factor), evaluated at the pre-solve state like
	! every other channel.  Zero when molecules are off/absent.
	if (thereis_mol) then
		do j = 1-Ng,N+Ng
			cool(j) = cool(j) + h3p_cooling_rate(T_K(j),               &
			                     nmol_eq(j,3), nmol_eq(j,1))
		enddo
	endif

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
	! gas. Lagged (pre-solve) densities, like every other channel above; nheiTR
	! is the same array he_rec_coupling already consumes.
	if (thereis_HeITR) heat = heat                                    &
	     + nheiTR*nhi*Q31*(e_th_HeI - e_th_HeTR - e_th_HI)/erg2eV

	! Molecular Penning ionization heating: He(2^3S)+H2 -> He(1^1S)+H2+ + e-
	! releases the electron kinetic energy (e_th_HeI - e_th_HeTR) - e_th_H2
	! (= 24.6 - 4.80 - 15.4 = 4.4 eV) into the gas. Lagged (pre-solve)
	! densities; nmol_eq(:,1) is the neutral-H2 number density. The rate
	! coefficient is the Garcia Munoz (2025) Table A.5 fit penning_HeI23S_H2
	! (Cool_coeff.f90). Zero unless a molecular run also tracks the triplet.
	if (thereis_mol .and. thereis_HeITR) heat = heat                  &
	     + nheiTR*nmol_eq(:,1)*penning_HeI23S_H2(T_K)                  &
	       *((e_th_HeI - e_th_HeTR) - e_th_H2)/erg2eV

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
		! per equilibrium sweep. n_molfail counts molecular cells that kept
		! their previous state; the n_ieq_* counters are the atomic root
		! validation (see the declarations above).
		n_molfail    = 0
		n_ieq_reseed = 0
		n_ieq_retry  = 0
		n_ieq_unphys = 0
		n_ieq_fail   = 0
		!$omp parallel do default(shared) schedule(dynamic,8) copyin(cx_metal_base) &
		!$omp   private(params, usednt, i0, top, im, meg_ntot, meg_g0, meg_g1,      &
		!$omp           meg_b0, meg_b1, meg_a1, meg_a2, meg_top,                     &
		!$omp           pbar_loc, qh2_loc, x2_loc,                                   &
		!$omp           iatt, info_ieq, ok_rank, best_rank, conv_ieq, phys_ieq,      &
		!$omp           x_root_best)                                                 &
		!$omp   reduction(+:n_molfail,n_ieq_reseed,n_ieq_retry,n_ieq_unphys,   &
		!$omp               n_ieq_fail) if(count > 0)
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
				ieq_cell%T_K  = T_K(j)
				ieq_cell%ntot = n_in_dim(j)      ! M for the 3-body rates
				! Compute the molecular rate coefficients that are invariant
				! across this cell's Newton solve (they depend only on T and
				! n_tot); the residual then reads them, like set_metal_coeffs.
				call set_mol_coeffs(T_K(j), n_in_dim(j))
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
				!   2  molecular basin -- the chemical-equilibrium H2 fit at
				!      the local (p, T); the H-nucleus fraction bound in H2
				!      seeds sys_x(4), with the metal stages from their own
				!      ionization balance at the incoming n_e,
				!   3  atomic basin -- the molecule-free ionization balance of
				!      every element at the incoming n_e (H2/H2+/H3+/HeH+ set
				!      to zero, the metastable and the metals to their own
				!      steady states),
				! and, as in the atomic branch below, a converged root is kept
				! only if it is also physical (every stage fraction >= 0, each
				! element's tracked stages summing to at most its nuclei).
				best_rank = 0
				do iatt = 1,3
					if (iatt .eq. 3) then
						call ionization_balance_at_fixed_ne(sys_x,N_eq,   &
						                                    mbase,ne(j))
						n_ieq_retry = n_ieq_retry + 1
					else if (iatt .eq. 2) then
						pbar_loc = n_in_dim(j)*kb_erg*T_K(j)/1.0d6   ! gas pressure [bar]
						qh2_loc  = q_h2_equilibrium(pbar_loc, T_K(j))
						x2_loc   = 2.0d0*qh2_loc*(1.0d0 + HeH)/(1.0d0 + qh2_loc)
						if (x2_loc .gt. 1.0d0) x2_loc = 1.0d0
						sys_x(1) = nhii(j)/nh(j)      ! near-neutral H (dark base)
						sys_x(2) = nheii(j)/nhe(j)
						sys_x(3) = nheiii(j)/nhe(j)
						sys_x(4) = x2_loc             ! 2 n_H2 / n_H from the fit
						sys_x(5) = 1.0d-10            ! tiny H2+ seed
						sys_x(6) = 1.0d-10            ! tiny H3+ seed
						sys_x(7) = 1.0d-10            ! tiny HeH+ seed
						if (thereis_HeITR) sys_x(8) = nheiTR(j)/nhe(j)
						if (thereis_metals)                                &
							call metal_ionization_balance_at_fixed_ne(     &
							                     sys_x,N_eq,mbase,ne(j))
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
					endif
					if (ok_rank .eq. 2) exit
				enddo

				if (best_rank .lt. 2) then
					if (best_rank .eq. 1) then
						! No attempt met the solver tolerance; keep the best
						! physical root found.
						sys_x(1:N_eq) = x_root_best(1:N_eq)
					else
						! No starting point produced an admissible root: keep
						! this cell's previous state by restoring sys_x to the
						! pre-solve fractions (the extraction below then
						! reproduces it unchanged), and count the cell for the
						! one-line summary warning.
						sys_x(1) = nhii(j)/nh(j)
						sys_x(2) = nheii(j)/nhe(j)
						sys_x(3) = nheiii(j)/nhe(j)
						sys_x(4) = 2.0d0*nmol_eq(j,1)/nh(j)
						sys_x(5) = 2.0d0*nmol_eq(j,2)/nh(j)
						sys_x(6) = 3.0d0*nmol_eq(j,3)/nh(j)
						sys_x(7) = nmol_eq(j,4)/nh(j)
						if (thereis_HeITR) sys_x(8) = nheiTR(j)/nhe(j)
						if (thereis_metals) then
							do im = 1,n_melem
								i0 = melem_i0(im)
								sys_x(mbase+2*(im-1)) = nm(j,i0+1)         &
								                  /max(nm_tot(j,im),1.0d-30)
								if (melem_top(im) .ge. 2) then
									sys_x(mbase+1+2*(im-1)) = nm(j,i0+2)       &
									                  /max(nm_tot(j,im),1.0d-30)
								else
									sys_x(mbase+1+2*(im-1)) = 0.0d0
								endif
							enddo
						endif
						n_molfail = n_molfail + 1
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

		! One summary line per sweep if any molecular cell failed to converge
		! twice (kept its previous state).
		if (thereis_mol .and. n_molfail .gt. 0) then
			write(*,'(A,I0,A)') ' (ioniz_eq) WARNING: molecular '//        &
				'equilibrium failed at ', n_molfail,                       &
				' cells (kept previous state)'
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

	endif

	
	! Density with atomic numbers (nm adds the metal mass under the
	! eos_metals policy)
   if (thereis_mol) then
      call calc_rho(nhi,nhii,nhei,nheii,nheiii,nheiTR,n_io,nm,nmol_eq)
   else
      call calc_rho(nhi,nhii,nhei,nheii,nheiii,nheiTR,n_io,nm)
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

	! Helium nuclei: He II, He III and the 2^3S metastable, which the systems
	! carry as a separate level inside the neutral stage.
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
