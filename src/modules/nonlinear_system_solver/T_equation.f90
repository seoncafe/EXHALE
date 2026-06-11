	module equation_T
	! Equation for temperature at the steady state

	use global_parameters
	use utils, only : calc_ne
	use Cooling_Coefficients
	use species_table, only : n_mion, mion_iscool, mion_name,           &
	                          mion_z2, mion_elem, mion_stage, melem_Z,  &
	                          im_FeII

	implicit none

	! Per-cell metal state for the advection post-process temperature solve.
	! post_process_adv sets these (cgs densities for the current cell, and
	! the on/off switch) immediately before each hybrd1(T_equation,...) call,
	! so the implicit T it converges to balances the SAME metal line cooling
	! that eval_cool reports. Module-scope (not via params(40)) because the
	! 27-ion metal vector does not fit the legacy params array. Safe because
	! the post-process temperature loop is serial.
	real*8  :: pp_nm_cell(n_mion) = 0.0d0   ! current-cell metal densities [cm^-3]
	logical :: pp_metal_on        = .false. ! add metal cooling/brem/ne in T_equation

	contains
	
	subroutine T_equation(N_T_eq,x,fvec,iflag,params)
	
	integer :: N_T_eq,iflag
	real*8  :: x(N_T_eq),fvec(N_T_eq)
   real*8  :: params(40)
	real*8  :: nhi,nhii
	real*8  :: nhei,nheii,nheiii
	real*8  :: ne
	real*8  :: mum,mup
	real*8  :: rhov
	real*8  :: coeff
	real*8  :: dr
	real*8  :: Told,heaold
	real*8  :: reco,coio,brem,coex,cool,cool_M
	real*8  :: TT
	integer :: im

	! Parameters
	nhi    = params(1)
	nhii   = params(2)
	nhei   = params(3)
	nheii  = params(4)
	nheiii = params(5)
  	mup    = params(6)
   mum    = params(7)
   rhov   = params(8)
   coeff  = params(9)
   dr	   = params(10)
   Told   = params(11)
   heaold = params(12)
   ! Metal densities are supplied per-cell through the module array
   ! pp_nm_cell (cgs), set by post_process_adv; pp_metal_on gates whether
   ! metals contribute. params(13-18) are no longer used.

   ! Free electron density (includes metal ions)
   if (thereis_He) then
		ne = nhii + nheii + 2.0*nheiii
	else
		ne = nhii
	endif
	if (pp_metal_on) then
		do im = 1,n_mion
			ne = ne + dble(mion_stage(im))*pp_nm_cell(im)
		enddo
	endif
      
	! Substitutions
	TT = x(1)*T0
	! With metal cooling active the energy residual is much stiffer, so the
	! hybrd1 search can transiently overshoot to a negative trial T. Evaluating
	! brem (sqrt(TT)) or the metal tables there yields NaN, which then poisons
	! the Newton step. Floor TT to a small positive value while metals are on so
	! every cooling term stays finite; the physical root sits far above the
	! floor, so the converged T is unaffected. Mode 0 (pp_metal_on=.false.)
	! keeps the exact original behaviour.
	if (pp_metal_on) TT = max(TT, 1.0d0)

	!--- Evaluate cooling rates ---!
			
    ! Cooling rate
   reco  =  rec_cool_HII_func(TT)*nhii     & ! HII
         +  rec_cool_HeII_func(TT)*nheii   & ! HeII
         +  rec_cool_HeIII_func(TT)*nheiii   ! HeIII

   !-- Collisional ionization --!
      
   ! Cooling rate
   coio =  2.179e-11*ion_coeff_HI_func(TT)*nhi  	         & ! HI
           + 3.940e-11*ion_coeff_HeI_func(TT)*nhei 		   & ! HeI
	  		  + kb_erg*631515.0*ion_coeff_HeII_func(TT)*nheii    ! HeII

   !-- Bremsstrahlung --!

   ! Cooling rate (H/He plus every charged metal ion, Z^2-weighted, matching
   ! the brem accumulator in eval_cool). Neutral stages carry z2 = 0.
   brem = ih**2.0*GF_func(TT,ih)*nhii                       &  ! HII
        + ihe**2.0*GF_func(TT,ihe)*(nheii + nheiii)            ! He
   if (pp_metal_on) then
      do im = 1,n_mion
         if (mion_z2(im) == 0) cycle
         brem = brem + dble(mion_z2(im))                          &
                       *GF_func(TT, dble(melem_Z(mion_elem(im)))) &
                       *pp_nm_cell(im)
      enddo
   endif
   brem = 1.426e-27*sqrt(TT)*brem

   !-- Collisional excitation --!

   ! Collisional excitation
   coex = coex_rate_HI_func(TT)*nhi       &    ! HI
        + coex_rate_HeI_func(TT)*nhei     &    ! HeI
        + coex_rate_HeII_func(TT)*nheii        ! HeII

   !-- Metal line cooling (optically thin, beta = 1, as in eval_cool).
   ! Sum the same mion_iscool coolants eval_cool sums, using the scalar
   ! coefficient dispatcher so the converged T balances the reported cooling.
   cool_M = 0.0d0
   if (pp_metal_on) then
      do im = 1,n_mion
         if (.not. mion_iscool(im)) cycle
         if (im .eq. im_FeII) then
            ! density-dependent Fe II (matches eval_cool's c_metal override:
            ! the local ne selects the coronal->LTE-saturated coefficient)
            cool_M = cool_M + pp_nm_cell(im)*cool_FeII_ne_scalar(TT, ne)
         else
            cool_M = cool_M + pp_nm_cell(im)                          &
                              *cool_coeff_by_ion_scalar(im, TT)
         endif
      enddo
   endif

   ! Total cooling rate in erg/(s cm^3)
 	cool = (ne*(brem + coex + reco + coio) + ne*cool_M)/q0

	! Equation
   fvec(1) = mum*rhov*x(1) - mup*rhov*Told 		&
           - (g-1.0)*(coeff*x(1) + mup*mum*dr*(heaold - cool))
      
   ! End of subroutine
	end subroutine T_equation

   !-------------------------------------------------------------------!
   ! Scalar root-finder for the energy equation (Task 1: Brent).
   !
   ! T_equation is a single nonlinear equation in x = T/T0. With metal
   ! cooling its residual is NON-monotone and admits a second, spurious
   ! *hot* root; a Newton/Powell solver (hybrd1) can land on it. Bracketing
   ! the first sign change scanning upward from a low floor structurally
   ! selects the lowest (physical) root, then Brent polishes it. Used only
   ! by the post-processor (one pass), so the scan cost is negligible.
   !-------------------------------------------------------------------!

   ! Scalar residual of the energy equation at x = T/T0 (wraps T_equation).
   double precision function Tres(xx, params)
   real*8, intent(in) :: xx, params(40)
   real*8  :: xv(1), fv(1)
   integer :: iflag
   xv(1) = xx
   iflag = 1
   call T_equation(1, xv, fv, iflag, params)
   Tres = fv(1)
   end function Tres

   ! Bracket the lowest physical root and polish it with Brent. x_guess is
   ! the eq value (T_in/T0) used to anchor the search window. On return
   ! ok=.false. means no sign change was found in [0.05,4]*x_guess, and the
   ! caller should fall back to the eq temperature.
   subroutine solve_T_brent(params, x_guess, x_out, ok)
   real*8,  intent(in)  :: params(40), x_guess
   real*8,  intent(out) :: x_out
   logical, intent(out) :: ok
   integer, parameter :: nscan = 80
   real*8  :: xlo, xhi, xprev, fprev, xx, ff, a, b, fa, fb
   integer :: i

   ok    = .false.
   x_out = x_guess
   if (x_guess .le. 0.0d0) return
   xlo = max(0.05d0*x_guess, 1.0d0/T0)   ! floor ~1 K, below any physical root
   xhi = 4.0d0*x_guess                   ! covers the physical root
   xprev = xlo
   fprev = Tres(xlo, params)
   if (fprev .eq. 0.0d0) then
      x_out = xlo; ok = .true.; return
   endif
   do i = 1, nscan
      xx = xlo*(xhi/xlo)**(dble(i)/dble(nscan))   ! log-spaced scan
      ff = Tres(xx, params)
      if (ff .eq. 0.0d0) then
         x_out = xx; ok = .true.; return
      endif
      if (fprev*ff .lt. 0.0d0) then       ! first sign change = lowest root
         a = xprev; fa = fprev; b = xx; fb = ff
         call brent_root(a, b, fa, fb, params, x_out)
         ok = .true.
         return
      endif
      xprev = xx; fprev = ff
   enddo
   end subroutine solve_T_brent

   ! Brent's method on a sign-changing bracket [a,b] of Tres(.,params).
   subroutine brent_root(a_in, b_in, fa_in, fb_in, params, root)
   real*8, intent(in)  :: a_in, b_in, fa_in, fb_in, params(40)
   real*8, intent(out) :: root
   real*8, parameter   :: tol = 1.0d-10
   integer, parameter  :: itmax = 100
   real*8  :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm, eps
   integer :: it

   eps = epsilon(1.0d0)
   a = a_in; b = b_in; fa = fa_in; fb = fb_in
   c = b; fc = fb; d = b - a; e = d
   do it = 1, itmax
      if (fb*fc .gt. 0.0d0) then
         c = a; fc = fa; d = b - a; e = d
      endif
      if (abs(fc) .lt. abs(fb)) then
         a = b; b = c; c = a
         fa = fb; fb = fc; fc = fa
      endif
      tol1 = 2.0d0*eps*abs(b) + 0.5d0*tol
      xm = 0.5d0*(c - b)
      if (abs(xm) .le. tol1 .or. fb .eq. 0.0d0) exit
      if (abs(e) .ge. tol1 .and. abs(fa) .gt. abs(fb)) then
         s = fb/fa
         if (a .eq. c) then
            p = 2.0d0*xm*s
            q = 1.0d0 - s
         else
            q = fa/fc
            r = fb/fc
            p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
            q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
         endif
         if (p .gt. 0.0d0) q = -q
         p = abs(p)
         if (2.0d0*p .lt. min(3.0d0*xm*q - abs(tol1*q), abs(e*q))) then
            e = d; d = p/q
         else
            d = xm; e = d
         endif
      else
         d = xm; e = d
      endif
      a = b; fa = fb
      if (abs(d) .gt. tol1) then
         b = b + d
      else
         b = b + sign(tol1, xm)
      endif
      fb = Tres(b, params)
   enddo
   root = b
   end subroutine brent_root

	! End of module
	end module equation_T
