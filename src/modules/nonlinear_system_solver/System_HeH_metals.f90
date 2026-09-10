	module System_HeH_metals
	! Ionization equilibrium system with H, He, and an arbitrary number of
	! metal elements. Each element carries up to two ionized stages; the
	! number actually solved is set per element by met_top (2 = neutral/+/++,
	! 1 = neutral/+ only).
	! Unknowns:
	!   x(1) = n_HII/n_H
	!   x(2) = n_HeII/n_He
	!   x(3) = n_HeIII/n_He
	! and, for each metal element e = 1..met_nelem (canonical order),
	!   x(4 + 2*(e-1)) = n_(Xe)II /n_Xe
	!   x(5 + 2*(e-1)) = n_(Xe)III/n_Xe   (pinned to 0 when met_top(e) < 2)
	!
	! The metal data for each element (total density + photoionization,
	! collisional-ionization and recombination coefficients) is supplied
	! per cell through the module-level arrays below, set by the driver via
	! set_metal_coeffs before each hybrd1 call. Those arrays are OpenMP
	! THREADPRIVATE (see the directive at their declaration), so each thread
	! of the parallel cell loop in ioniz_eq carries the coefficients of the
	! cell it is solving; the arrangement lets the system grow with the
	! number of metals without changing the argument list.
	!
	! Every photoionization rate this system reads -- ieq_cell%P_HI, P_HeI,
	! P_HeII and the met_g0/met_g1 of each element -- is a CONSTANT of the
	! solve although the field that sets it is a function of these unknowns.
	! The reason, and where that dependence is closed instead, are stated
	! once in System_HeH.
	!
	! Charge exchange with H/He (and, in full mode, metal-metal) is added
	! generically by cx_add_to_fvec from the charge_exchange module after
	! the photo/collisional/recombination balance is built (Huang et al.
	! 2023, Table 4). Absent elements (n_X <= 0) are force-zeroed, so the
	! same system handles C-only, C+O, C+N+O, C+N+O+Mg, etc.; charge
	! exchange vanishes automatically for any absent reactant.

	use global_parameters
	use ion_cell_state, only: ieq_cell
	use ion_residual_core, only: heh_rows, heh_crow, heh_jac_local,       &
	                             metal_fractions, metal_electron_sum,     &
	                             metal_rows
	use charge_exchange, only: cx_add_to_fvec, cx_add_to_jac,             &
	                           he_h_cx_fvec, he_h_cx_jac

	implicit none

	! Cell-by-cell metal element data, set by set_metal_coeffs.
	integer, save :: met_nelem = 0
	real*8, allocatable, save :: met_ntot(:)   ! total density of the element
	real*8, allocatable, save :: met_g0(:)     ! photoionization of neutral
	real*8, allocatable, save :: met_g1(:)     ! photoionization of singly ionized
	real*8, allocatable, save :: met_b0(:)     ! collisional ionization of neutral
	real*8, allocatable, save :: met_b1(:)     ! collisional ionization of singly ionized
	real*8, allocatable, save :: met_a1(:)     ! recombination of singly ionized
	real*8, allocatable, save :: met_a2(:)     ! recombination of doubly ionized
	! Highest ionization stage tracked per element (2 = up to X++, i.e. the
	! full three-stage system; 1 = only neutral + X+). For a two-stage
	! element the upper unknown x(ix+1) carries no physics and is pinned to
	! zero in the residual to keep the Jacobian non-singular.
	integer, allocatable, save :: met_top(:)

	! These cell-by-cell coefficients are set (set_metal_coeffs) and read inside the
	! ionization-equilibrium cell sweep, which is now OpenMP-parallel over cells.
	! Make each thread keep its own copy so concurrent cells do not clobber one
	! another. The allocatables are lazily allocated per thread on first use in
	! set_metal_coeffs (its `if (.not.allocated)` guard now runs per thread).
	!$omp threadprivate(met_nelem, met_ntot, met_g0, met_g1, met_b0, met_b1,  &
	!$omp                met_a1, met_a2, met_top)

	contains

	! Store the cell-by-cell metal element coefficients for ion_system_HeH_metals.
	subroutine set_metal_coeffs(nelem, ntot, g0, g1, b0, b1, a1, a2, top)
		integer, intent(in) :: nelem
		real*8, dimension(nelem), intent(in) :: ntot, g0, g1, b0, b1, a1, a2
		integer, dimension(nelem), intent(in) :: top
		if (.not. allocated(met_ntot)) then
			allocate(met_ntot(nelem), met_g0(nelem), met_g1(nelem),   &
			         met_b0(nelem), met_b1(nelem), met_a1(nelem),     &
			         met_a2(nelem), met_top(nelem))
		endif
		met_nelem = nelem
		met_ntot  = ntot
		met_g0    = g0
		met_g1    = g1
		met_b0    = b0
		met_b1    = b1
		met_a1    = a1
		met_a2    = a2
		met_top   = top
	end subroutine set_metal_coeffs

	subroutine ion_system_HeH_metals(N_eq,x,fvec,iflag,params)

	integer :: N_eq,iflag
	real*8  :: x(N_eq),fvec(N_eq)
	real*8  :: params(60)

	! H/He coefficients (params 1-11)
	real*8  :: g_hi,g_hei,g_heii        ! photoionization
	real*8  :: a_hii,a_heii,a_heiii     ! recombination
	real*8  :: b_hi,b_hei,b_heii        ! collisional ionization
	real*8  :: n_h,n_he,n_e
	! H/He densities
	real*8  :: n_hi,n_hii,n_hei,n_heii,n_heiii
	! Each element's metal densities (neutral/+/++). Charge exchange is added
	! generically afterwards by cx_add_to_fvec, so no separate CX arrays for each element
	! are needed here.
	real*8  :: nm0(met_nelem),nm1(met_nelem),nm2(met_nelem)
	real*8  :: n_X
	integer :: e,ix

	! Unpack H/He coefficients
	g_hi    = ieq_cell%P_HI
	g_hei   = ieq_cell%P_HeI
	g_heii  = ieq_cell%P_HeII
	a_hii   = ieq_cell%rchiiB
	a_heii  = ieq_cell%rcheiiB
	a_heiii = ieq_cell%rcheiiiB
	n_h     = ieq_cell%nh
	n_he    = ieq_cell%nhe
	b_hi    = ieq_cell%a_ion_HI
	b_hei   = ieq_cell%a_ion_HeI
	b_heii  = ieq_cell%a_ion_HeII

	! H/He densities from fractions
	n_hii   = x(1)*n_h
	n_hi    = (1.0 - x(1))*n_h
	n_heii  = x(2)*n_he
	n_heiii = x(3)*n_he
	n_hei   = (1.0 - x(2) - x(3))*n_he

	! Metal ion densities (neutral/+/++) from fractions, canonical order
	call metal_fractions(x, 4, met_nelem, met_ntot, nm0, nm1, nm2)

	! Electron density (X+ counts once, X++ twice). H/He first, then the
	! metals in canonical element order, reproducing the original sum.
	n_e = n_hii + n_heii + 2.0*n_heiii
	call metal_electron_sum(n_e, met_nelem, nm1, nm2)

	! Steady-state balance equations: fvec = 0 when
	! (production from lower stage) = (loss to upper stage). Charge
	! exchange is added to the relevant rows by cx_add_to_fvec below.

	! HI <-> HII, HeI <-> HeII, HeII <-> HeIII (standard H/He rows, shared
	! helper; n_e here already includes the metal electrons).
	call heh_rows(fvec, n_hi, n_hii, n_hei, n_heii, n_heiii, n_e,  &
	              g_hi, g_hei, g_heii, a_hii, a_heii, a_heiii,      &
	              b_hi, b_hei, b_heii)

	! Metal ionization balance, one element at a time (force-zero if the
	! element is absent; otherwise normal balance). Three-stage elements
	! solve both X0<->X+ and X+<->X++; two-stage elements solve only
	! X0<->X+ and pin the unused upper unknown to zero.
	call metal_rows(fvec, x, 4, met_nelem, met_ntot, met_g0, met_g1,   &
	                met_b0, met_b1, met_a1, met_a2, met_top,           &
	                nm0, nm1, nm2, n_e)

	! Add charge exchange (Huang Table 4) to the H, He and metal balance
	! rows. Rates were stored per cell by cx_set_cell; absent reactants
	! contribute zero, preserving the zero-abundance identity rows above.
	call cx_add_to_fvec(N_eq, fvec, nm0, nm1, nm2,                  &
	                    n_hi, n_hii, n_hei, n_heii, n_heiii)

	! He <-> H charge exchange (Huang Table 4 group B). Standard He row
	! (HeI->HeII positive), so he_row_sign = +1.
	call he_h_cx_fvec(fvec, ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0,  &
	                  n_hi, n_hii, n_hei, n_heii, 1.0d0)

	return

	end subroutine ion_system_HeH_metals

   ! Analytic Jacobian of ion_system_HeH_metals (Task 2). Each row has the form
   !   f_i = P_i(x_local) + C_i(x_local) * n_e
   ! so  J(i,k) = [photoion + dC_i/dx_local * n_e]_local  +  C_i * dn_e/dx_k.
   ! The second term is a rank-1 "n_e coupling" common to every row; the first
   ! is local to each element's own unknowns. Charge exchange is added by
   ! cx_add_to_jac (mirror of cx_add_to_fvec). Absent elements and the unused
   ! upper stage of two-stage elements are pinned to identity rows, matching the
   ! residual's fvec(ix)=x(ix) identities.
   subroutine jac_system_HeH_metals(N_eq,x,fjac,params)
   integer :: N_eq
   real*8  :: x(N_eq), fjac(N_eq,N_eq), params(60)
   real*8  :: g_hi,g_hei,g_heii, a_hii,a_heii,a_heiii, b_hi,b_hei,b_heii
   real*8  :: n_h,n_he,n_e
   real*8  :: n_hi,n_hii,n_hei,n_heii,n_heiii
   real*8  :: nm0(met_nelem),nm1(met_nelem),nm2(met_nelem)
   real*8  :: dne(N_eq), Crow(N_eq), n_X
   integer :: e,ix,i,k

   ! Unpack H/He coefficients (named cell state, same as the residual).
   g_hi=ieq_cell%P_HI; g_hei=ieq_cell%P_HeI; g_heii=ieq_cell%P_HeII
   a_hii=ieq_cell%rchiiB; a_heii=ieq_cell%rcheiiB; a_heiii=ieq_cell%rcheiiiB
   n_h=ieq_cell%nh; n_he=ieq_cell%nhe
   b_hi=ieq_cell%a_ion_HI; b_hei=ieq_cell%a_ion_HeI; b_heii=ieq_cell%a_ion_HeII

   n_hii=x(1)*n_h;  n_hi=(1.0-x(1))*n_h
   n_heii=x(2)*n_he; n_heiii=x(3)*n_he; n_hei=(1.0-x(2)-x(3))*n_he
   call metal_fractions(x, 4, met_nelem, met_ntot, nm0, nm1, nm2)
   n_e = n_hii + n_heii + 2.0*n_heiii
   call metal_electron_sum(n_e, met_nelem, nm1, nm2)

   ! Electron-density gradient dne(k) = d n_e / d x_k.
   dne = 0.0d0
   dne(1)=n_h; dne(2)=n_he; dne(3)=2.0*n_he
   do e=1,met_nelem
      ix=4+2*(e-1)
      dne(ix)   = met_ntot(e)
      dne(ix+1) = 2.0*met_ntot(e)
   enddo

   ! n_e-coefficient C_i of each row (the factor multiplying n_e in fvec_i).
   Crow = 0.0d0
   call heh_crow(Crow(1), Crow(2), Crow(3), n_hi, n_hii, n_hei, n_heii,  &
                 n_heiii, a_hii, a_heii, a_heiii, b_hi, b_hei, b_heii)
   do e=1,met_nelem
      ix=4+2*(e-1)
      if (met_ntot(e) .gt. 1.0d-30) then
         Crow(ix) = nm0(e)*met_b0(e) - met_a1(e)*nm1(e)
         if (met_top(e) .ge. 2) Crow(ix+1) = nm1(e)*met_b1(e) - met_a2(e)*nm2(e)
      endif
   enddo

   ! Rank-1 n_e coupling for every row: fjac(i,k) = C_i * dne(k).
   do i=1,N_eq
      do k=1,N_eq
         fjac(i,k) = Crow(i)*dne(k)
      enddo
   enddo

   ! Local "direct" terms (photoionization + dC_i/dx_local * n_e). H/He rows:
   call heh_jac_local(N_eq, fjac, n_h, n_he, n_e, g_hi, g_hei, g_heii,  &
                      a_hii, a_heii, a_heiii, b_hi, b_hei, b_heii)
   ! Metal rows (present elements only):
   do e=1,met_nelem
      ix=4+2*(e-1); n_X=met_ntot(e)
      if (n_X .le. 1.0d-30) cycle
      ! Row ix: X0<->X+ , f = nm0*g0 + (nm0*b0 - a1*nm1)*n_e
      fjac(ix,ix)   = fjac(ix,ix)   - n_X*met_g0(e)                    &
                    + (-n_X*met_b0(e) - met_a1(e)*n_X)*n_e
      fjac(ix,ix+1) = fjac(ix,ix+1) - n_X*met_g0(e) + (-n_X*met_b0(e))*n_e
      if (met_top(e) .ge. 2) then
         ! Row ix+1: X+<->X++ , f = nm1*g1 + (nm1*b1 - a2*nm2)*n_e
         fjac(ix+1,ix)   = fjac(ix+1,ix)   + n_X*met_g1(e) + (n_X*met_b1(e))*n_e
         fjac(ix+1,ix+1) = fjac(ix+1,ix+1) + (-met_a2(e)*n_X)*n_e
      endif
   enddo

   ! Charge exchange (Huang Table 4) -- mirror of cx_add_to_fvec.
   call cx_add_to_jac(N_eq, fjac, nm0, nm1, nm2,                       &
                      n_hi, n_hii, n_hei, n_heii, n_heiii,             &
                      met_ntot, n_h, n_he)

   ! He <-> H charge-exchange Jacobian (rows 1,2; cols 1,2,3), mirror of
   ! the he_h_cx_fvec call in the residual (standard He row, +1).
   call he_h_cx_jac(N_eq, fjac, ieq_cell%kcx_He0_Hp, ieq_cell%kcx_Hep_H0,  &
                    n_h, n_he, n_hi, n_hii, n_hei, n_heii)

   ! Pinned/identity rows LAST (CX never targets them): absent elements pin
   ! both stages; two-stage elements pin the unused X++ unknown.
   do e=1,met_nelem
      ix=4+2*(e-1)
      if (met_ntot(e) .le. 1.0d-30) then
         fjac(ix,:)   = 0.0d0;  fjac(ix,ix)     = 1.0d0
         fjac(ix+1,:) = 0.0d0;  fjac(ix+1,ix+1) = 1.0d0
      else if (met_top(e) .lt. 2) then
         fjac(ix+1,:) = 0.0d0;  fjac(ix+1,ix+1) = 1.0d0
      endif
   enddo
   end subroutine jac_system_HeH_metals

	end module System_HeH_metals
