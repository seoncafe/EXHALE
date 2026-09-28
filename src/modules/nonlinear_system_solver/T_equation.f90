	module equation_T
	! Equation for temperature at the steady state

	use global_parameters
	use ion_cell_state, only: teq_cell
	use caloric_eos, only: internal_energy_of_mixture
	use utils, only : calc_ne
	use Cooling_Coefficients
	use species_table, only : n_mion, mion_stage

	implicit none

	! Cell-by-cell metal state for the advection post-process temperature solve.
	! post_process_adv sets these (cgs densities for the current cell, and
	! the on/off switch) immediately before each hybrd1(T_equation,...) call.
	! The cooling T_equation balances is radiative_cooling_of_cell
	! (Cool_coeff), the one assembly eval_cool also evaluates, so the T it
	! converges to balances the cooling eval_cool reports at it.
	! Module-scope (not via params(40)) because the 27-ion metal vector does
	! not fit the legacy params array. Safe because the post-process
	! temperature loop is serial.
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
	! Photon occupation number of the lower-atmosphere infrared field at each
	! of those lines, set from the same profile. 0 = no incident field, the
	! value used when metals are off or "Base IR field" is off.
	real*8  :: pp_nbar_fs(n_fsline) = 0.0d0

	! THE BRACKET TOLERANCE OF THE TEMPERATURE ROOT (brent_root), and the
	! hook that replaces it for a whole run (EXHALE_TEQ_TOL=<x>, default
	! off). Brent stops when the bracket has shrunk to
	! 2*epsilon*|b| + tol/2 in x = T/T0, so 1e-10 pins the root to about
	! 1e-10*T0 ~ 1e-6 K and the root is a piecewise map of the state at that
	! size. A non-positive or unreadable value leaves it at 1e-10.
	real*8,  save :: teq_bracket_tol      = 1.0d-10
	logical, save :: teq_bracket_tol_read = .false.

	contains

	! THE BRACKET TOLERANCE THE TEMPERATURE ROOT IS ASKED FOR. Returns
	! what EXHALE_TEQ_TOL names, or the argument when the variable is
	! unset, unreadable or non-positive. Pure of module state, so a caller
	! may cache the answer and a test may ask it twice with two
	! environments.
	real*8 function temperature_bracket_tolerance(tol_default)
	real*8, intent(in) :: tol_default
	character(len=32)  :: env_teq
	real*8             :: tol_asked
	temperature_bracket_tolerance = tol_default
	call get_environment_variable('EXHALE_TEQ_TOL', env_teq)
	if (len_trim(env_teq) .eq. 0) return
	tol_asked = 0.0d0
	read(env_teq,*,err=311,end=311) tol_asked
  311	continue
	if (tol_asked .gt. 0.0d0) temperature_bracket_tolerance = tol_asked
	end function temperature_bracket_tolerance

	
	subroutine T_equation(N_T_eq,x,fvec,iflag,params)
	
	integer :: N_T_eq,iflag
	real*8  :: x(N_T_eq),fvec(N_T_eq)
   real*8  :: params(40)
	real*8  :: nhi,nhii
	real*8  :: nheiS,nheiTR,nheii,nheiii
	real*8  :: ne
	real*8  :: mum,mup
	real*8  :: rhov
	real*8  :: coeff
	real*8  :: dr_step
	real*8  :: heaold
	real*8  :: x_h2
	real*8  :: e_hist
	real*8  :: e_cell
	real*8  :: div_rhov
	real*8  :: cool
	real*8  :: TT
	! The cooling channels of the cell, and its metal densities as the cell
	! assembly takes them (zero with the metals off)
	real*8  :: chan(n_cool_chan), nm_c(n_mion)
	integer :: im

	! Parameters. The energy-equation coefficients are read from the named-field
	! teq_cell state; params stays only the MINPACK transport argument (unread),
	! passed through by hybrd1 / solve_T_brent / Tres.
	nhi    = teq_cell%nhi
	nhii   = teq_cell%nhii
	nheiS  = teq_cell%nheiS
	nheiTR = teq_cell%nheiTR
	nheii  = teq_cell%nheii
	nheiii = teq_cell%nheiii
  	mup    = teq_cell%mup
   mum    = teq_cell%mum
   rhov   = teq_cell%rhov
   coeff  = teq_cell%coeff
   dr_step = teq_cell%dr_step
   heaold = teq_cell%heaold
   x_h2   = teq_cell%x_h2
   e_hist = teq_cell%e_hist
   div_rhov = teq_cell%div_rhov
   ! Metal densities are supplied cell-by-cell through the module array
   ! pp_nm_cell (cgs), set by post_process_adv; pp_metal_on gates whether
   ! metals contribute.

   ! Free electron density, the charge sum of calc_ne (utilities) for the
   ! molecule-free post-process gas, which is what eval_cool forms from the
   ! same densities there: H II, He II and He III, and the metal stages
   ! under the eos_metals policy.
   if (thereis_He) then
		ne = nhii + nheii + 2.0d0*nheiii
	else
		ne = nhii
	endif
	nm_c = 0.0d0
	if (pp_metal_on) nm_c = pp_nm_cell
	if (pp_metal_on .and. eos_include_metals .and. thereis_metals) then
		do im = 1,n_mion
			if (mion_stage(im) .gt. 0)                                  &
				ne = ne + dble(mion_stage(im))*nm_c(im)
		enddo
	endif
      
	! Substitutions
	TT = x(1)*T0
	! The hybrd1 search can transiently overshoot to a negative trial T. Every
	! cooling rate below is a fit in T that has no value there: the free-free
	! term takes sqrt(TT), and the recombination coefficients are a NaN for
	! TT < 0 and diverge as TT -> 0+. The NaN then poisons the Newton step,
	! and under -ffpe-trap=invalid it aborts the run (measured in the
	! post-process of the He/H = 1 molecular case, which is metals-off). So
	! floor the argument of the rate functions at a small positive temperature,
	! unconditionally: it is the domain of the fits, not an option. The
	! physical root sits far above the floor -- the coldest base in the
	! repository is ~200 K -- so the converged T is unaffected, and the residual
	! stays real and finite where the solver is only passing through.
	TT = max(TT, 1.0d0)

	!--- Cooling rate ---!
	! The atomic, ionic and metal channels of the cell from the ONE assembly
	! (Cool_coeff: radiative_cooling_of_cell), at this trial temperature and
	! the cell's densities, frozen ground-capture escape weights and
	! fine-structure line transfer. The post-process gas carries no
	! molecules, so eval_cool adds nothing to these channels there either.
	call radiative_cooling_of_cell(TT, ne, nhi, nhii, nheiS, nheiTR,     &
	                               nheii, nheiii, teq_cell%y_HI,          &
	                               teq_cell%y_gnd, teq_cell%y_HeII, nm_c, &
	                               pp_beta_fs, pp_nbar_fs, chan)

	! Total cooling rate in code units
	cool = sum(chan)/q0

	! Equation
   ! Steady internal-energy equation of the profile,
   !
   !     div(u v) + p div(v)  =  heating - cooling
   !   = rho v de/dr  -  p v dln(rho)/dr  +  h div(rho v)
   !
   ! with u the internal energy density, e the internal energy per unit mass
   ! and h = e + p/rho the enthalpy per unit mass. The second form follows
   ! from div(rho e v) = e div(rho v) + rho v de/dr and
   ! p div(v) = (p/rho) div(rho v) - p v dln(rho)/dr, so the three terms
   ! below are that equation and not an approximation to it.
   !
   ! THE STEP. Every r-derivative of the equation is taken at r_j with the
   ! one rule of the step (ion_cell_state, adv_rates): for any profile f,
   !
   !     df/dr (r_j)  ~  (f_j - f_hist)/dr_step ,
   !     f_hist = a1 f_{j-1} - a2 f_{j-2} ,  dr_step = g h_j ,
   !
   ! the variable-step BDF2 differentiation formula, or backward Euler
   ! (a1 = g = 1, a2 = 0). The unknown enters through e_j alone, the
   ! density through rho_j - rho_hist (coeff), and every other factor is
   ! the value at r_j. Multiplying the equation by mup*mum*dr_step, with
   ! e = E(x_H2,T)/mu and p/rho = T/mu (E the energy per particle of
   ! internal_energy_of_mixture, mu the mean molecular weight, mup for this
   ! cell and mum for the upstream one), gives the residual below term by
   ! term:
   !
   !     mum*rhov*E(x_H2,T)               rho v e_j
   !   - mup*mum*rhov*e_hist              rho v e_hist
   !   + mum*dr_step*div_rhov*(E + T)     dr_step h_j div(rho v)
   !   - coeff*T                          (p/rho) v (rho_j - rho_hist)
   !   - mup*mum*dr_step*(heat - cool)    dr_step (heating - cooling)
   !
   ! THE HISTORY IS THE ENERGY THE FLOW CARRIES IN, so each e_k of e_hist is
   ! the SPECIFIC internal energy of its own cell, E(x_H2,k, T_k)/mu_k,
   ! formed by the caller at that cell's composition and temperature: the
   ! rovibrational heat capacity of H2 belongs to the gas that holds the
   ! molecules, so across a dissociation front two cells at the same
   ! temperature store different energy, and a cell with no H2 fed by a
   ! molecular neighbor still receives molecular energy.
   !
   ! ONE FORM FOR EVERY GAS. internal_energy_of_mixture returns
   ! T/(gamma_ad - 1) exactly for an atomic cell, so the monatomic equation
   ! is this residual with that E; it is not a separate branch. (Until
   ! 2026-09-27 the atomic interface was written out separately, multiplied
   ! by gamma_ad - 1, to keep an older arithmetic to the bit.)
   !
   ! div_rhov = 0 is a stationary mass flux, for which the enthalpy flux
   ! term is absent from the equation; the product is then an exact zero
   ! and the residual is the advected balance of the other terms alone.
   e_cell  = internal_energy_of_mixture(x_h2, x(1))
   fvec(1) = mum*rhov*e_cell                                            &
           - mup*mum*rhov*e_hist                                        &
           + mum*dr_step*div_rhov*(e_cell + x(1))                       &
           - (coeff*x(1) + mup*mum*dr_step*(heaold - cool))

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
   integer, parameter  :: itmax = 100
   real*8  :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm, eps, tol
   integer :: it

   ! Read once for the run, and announced only when it is armed, so that a
   ! run without the hook is the run without the code.
   if (.not. teq_bracket_tol_read) then
      teq_bracket_tol      = temperature_bracket_tolerance(teq_bracket_tol)
      teq_bracket_tol_read = .true.
      if (teq_bracket_tol .ne. 1.0d-10)                                   &
         write(*,'(A,ES11.3)') ' (T_equation) EXHALE_TEQ_TOL: the'//      &
              ' temperature bracket closes at ', teq_bracket_tol
   endif
   tol = teq_bracket_tol

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
