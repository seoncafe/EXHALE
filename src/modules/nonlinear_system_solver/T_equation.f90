	module equation_T
	! Equation for temperature at the steady state

	use global_parameters
	use ion_cell_state, only: teq_cell
	use utils, only : calc_ne
	use Cooling_Coefficients
	use species_table, only : n_mion, mion_iscool, mion_name,           &
	                          mion_z2, mion_elem, mion_stage, melem_Z,  &
	                          im_CI, im_CII, im_NII, im_OI, im_FeII

	implicit none

	! Cell-by-cell metal state for the advection post-process temperature solve.
	! post_process_adv sets these (cgs densities for the current cell, and
	! the on/off switch) immediately before each hybrd1(T_equation,...) call,
	! so the implicit T it converges to balances the SAME metal line cooling
	! that eval_cool reports. Module-scope (not via params(40)) because the
	! 27-ion metal vector does not fit the legacy params array. Safe because
	! the post-process temperature loop is serial.
	real*8  :: pp_nm_cell(n_mion) = 0.0d0   ! current-cell metal densities [cm^-3]
	logical :: pp_metal_on        = .false. ! add metal cooling/brem/ne in T_equation
	! Line-center escape probabilities of the ground-term fine-structure
	! lines for the current cell, set alongside pp_nm_cell from the profile
	! the post-process starts from. They are HELD FIXED while the root finder
	! varies T, exactly as the metal densities are: the optical depth is a
	! column over the whole atmosphere above the cell, so it is not a function
	! of this cell's trial temperature alone. 1 = optically thin, the value
	! used when metals are off.
	real*8  :: pp_beta_fs(n_fsline) = 1.0d0

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
	real*8  :: GF_z1,GF_z2   ! free-free Gaunt at ion net charge Z_ion = 1, 2
	integer :: im

	! Parameters. The energy-equation coefficients are read from the named-field
	! teq_cell state; params stays only the MINPACK transport argument (unread),
	! passed through by hybrd1 / solve_T_brent / Tres.
	nhi    = teq_cell%nhi
	nhii   = teq_cell%nhii
	nhei   = teq_cell%nhei
	nheii  = teq_cell%nheii
	nheiii = teq_cell%nheiii
  	mup    = teq_cell%mup
   mum    = teq_cell%mum
   rhov   = teq_cell%rhov
   coeff  = teq_cell%coeff
   dr	   = teq_cell%dr
   Told   = teq_cell%Told
   heaold = teq_cell%heaold
   ! Metal densities are supplied cell-by-cell through the module array
   ! pp_nm_cell (cgs), set by post_process_adv; pp_metal_on gates whether
   ! metals contribute.

   ! Free electron density (includes metal ions). Molecular ions are
   ! deliberately omitted as trace electron donors (negligible in the atomic
   ! post-process gas where this energy residual is solved).
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

   ! Cooling rate (H/He plus every charged metal ion). Free-free scales with
   ! the ion NET charge Z_ion: H II, He II and singly-ionized metals are
   ! Z_ion = 1; He III and doubly-ionized metals are Z_ion = 2. The Gaunt
   ! factor is evaluated at Z_ion. Matches the brem accumulator in eval_cool.
   GF_z1 = GF_func(TT, 1.0d0)
   GF_z2 = GF_func(TT, 2.0d0)
   brem = GF_z1*nhii + GF_z1*nheii + 4.0d0*GF_z2*nheiii   ! HII, HeII, HeIII
   if (pp_metal_on) then
      do im = 1,n_mion
         if (mion_z2(im) == 0) cycle
         if (mion_stage(im) == 2) then
            brem = brem + dble(mion_z2(im))*GF_z2*pp_nm_cell(im)
         else
            brem = brem + dble(mion_z2(im))*GF_z1*pp_nm_cell(im)
         endif
      enddo
   endif
   brem = 1.426e-27*sqrt(TT)*brem

   !-- Collisional excitation --!

   ! Collisional excitation
   coex = coex_rate_HI_func(TT)*nhi       &    ! HI
        + coex_rate_HeI_func(TT)*nhei     &    ! HeI
        + coex_rate_HeII_func(TT)*nheii        ! HeII

   !-- Metal line cooling. Sum the same mion_iscool coolants eval_cool
   ! sums, using the scalar coefficient dispatcher so the converged T
   ! balances the reported cooling. Trapping of the ground-term
   ! fine-structure lines enters through pp_beta_fs, everything else is
   ! optically thin.
   cool_M = 0.0d0
   if (pp_metal_on) then
      do im = 1,n_mion
         if (.not. mion_iscool(im)) cycle
         if (im .eq. im_FeII) then
            ! density-dependent Fe II (matches eval_cool's c_metal override:
            ! the local ne selects the coronal->LTE-saturated coefficient)
            cool_M = cool_M + pp_nm_cell(im)*cool_FeII_ne_scalar(TT, ne)
         else if (cno_chianti .and. im .eq. im_CI) then
            ! saturated [C I] 609/370um ground term (matches eval_cool)
            cool_M = cool_M + pp_nm_cell(im)                          &
                              *cool_CI_ne_func(TT, ne, nhi,           &
                                 pp_beta_fs(ifs_CI609),               &
                                 pp_beta_fs(ifs_CI370))
         else if (cno_chianti .and. im .eq. im_CII) then
            ! saturated [C II] 158um ground term (matches eval_cool)
            cool_M = cool_M + pp_nm_cell(im)                          &
                              *cool_CII_ne_func(TT, ne, nhi,          &
                                 pp_beta_fs(ifs_CII158))
         else if (cno_chianti .and. im .eq. im_NII) then
            ! saturated [N II] 205/122um ground term (matches eval_cool)
            cool_M = cool_M + pp_nm_cell(im)                          &
                              *cool_NII_ne_func(TT, ne, nhi,          &
                                 pp_beta_fs(ifs_NII205),              &
                                 pp_beta_fs(ifs_NII122))
         else if (cno_chianti .and. im .eq. im_OI) then
            ! saturated [O I] 63/145/44um ground term (matches eval_cool)
            cool_M = cool_M + pp_nm_cell(im)                          &
                              *cool_OI_ne_func(TT, ne, nhi,           &
                                 pp_beta_fs(ifs_OI63),                &
                                 pp_beta_fs(ifs_OI145),               &
                                 pp_beta_fs(ifs_OI44))
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
