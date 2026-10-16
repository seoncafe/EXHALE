      module newton_solver
      use global_parameters,   only: use_newton_ieq
      use chemical_row_layout, only: fraction_resolution
      ! THE CELL SOLVE OF THE IONIZATION AND CHEMICAL SYSTEMS, AND WHEN IT
      ! HAS CONVERGED.
      !
      ! Every unknown of these systems is the fraction of the nuclei of one
      ! element held by one species (chemical_row_layout). A trace species
      ! (H3+ at 1e-6 of the hydrogen, He++ at 1e-15 of the helium, H+ at
      ! 1e-14 of the hydrogen in molecular gas) is as much a part of the
      ! answer as an O(1) fraction: its density enters its own cooling, the
      ! electron budget of a shielded layer, and the carried stages of the
      ! transport solve. The solve has therefore converged when the Newton
      ! correction at the returned state resolves EVERY unknown relative to
      ! its own magnitude,
      !     |dx_i| <= rtol |x_i| + fraction_resolution   for every i,
      ! with rtol the stopping tolerance of the sweep (sqrt(epsilon) unless
      ! EXHALE_IEQ_TOL names another) and fraction_resolution = epsilon(1),
      ! the spacing of the numbers at the element total, below which no
      ! change of a fraction moves any sum of its element (the resolution
      ! floor of the acceptance measure, normalized_reaction_residual). The
      ! correction dx = -J^{-1} f is the Newton estimate of the remaining
      ! error of each unknown, so after a step that meets the test the
      ! state is resolved to about rtol^2 of each unknown: far inside the
      ! acceptance of every kinetic row (|f_i| <= 1e-6 of its group's
      ! turnover, 1e-3 of its own gross rate for H3+).
      !
      ! Neither test that stood before is a statement about each unknown:
      ! an absolute step of 1e-11 (1 + max |x|) leaves a stage whose root
      ! is 1e-11 at a relative error of order one, and MINPACK's relative
      ! step on the scaled norm of all unknowns (hybrd1, xtol) is dominated
      ! by the O(1) fractions, so a trace unknown is left unresolved
      ! (md/Update_EXHALE_stage3.md sections 116.2, 116.6 and 118.1).
      !
      ! Interfaces (matching the existing residuals so they can be passed as
      ! actual arguments):
      !   fcn(n, x, fvec, iflag, params)     -- residual (MINPACK form)
      !   jac(n, x, fjac, params)            -- analytic Jacobian d fvec_i/d x_k

      implicit none
      private
      public :: solve_ieq, newton_dense, newton_fd,                        &
                hybrd1_resolving_every_unknown,                           &
                correction_resolves_every_unknown,                        &
                newton_solver_read_environment
      ! Run-wide usage counters (reported once at the end of the run): how many
      ! ionization-equilibrium solves the analytic-Jacobian Newton handled vs.
      ! how many fell back to MINPACK hybrd1.
      integer, save, public :: nt_calls = 0, nt_fallback = 0
      ! Reproducible-validation switch: if the environment variable
      ! EXHALE_FORCE_HYBRD1 is set (to a non-empty value), solve_ieq skips Newton
      ! and always uses MINPACK hybrd1. Lets the same binary produce both the
      ! Newton and the reference hybrd1 outputs for an A/B comparison.
      logical, save :: nt_force = .false.

      contains

      subroutine newton_solver_read_environment()
      ! Serial, once: the EXHALE_FORCE_HYBRD1 validation switch (any value
      ! set = every ionization cell solve goes to hybrd1 instead of Newton).
      character(len=8) :: envval
      call get_environment_variable('EXHALE_FORCE_HYBRD1', envval)
      nt_force = (len_trim(envval) .gt. 0)
      end subroutine newton_solver_read_environment


      ! The stopping test of the cell solve (module header): the Newton
      ! correction dx at x resolves every unknown relative to its own
      ! magnitude, or to the resolution of a fraction where the unknown is
      ! smaller than that. A NaN component fails it.
      logical function correction_resolves_every_unknown(n, dx, x, rtol)  &
                       result(ok)
      integer, intent(in) :: n
      real*8,  intent(in) :: dx(n), x(n), rtol
      integer :: i
      ok = .false.
      do i = 1, n
         if (.not. (abs(dx(i)) .le. rtol*abs(x(i)) + fraction_resolution))&
            return
      enddo
      ok = .true.
      end function correction_resolves_every_unknown


      ! Analytic-Jacobian Newton first; on failure hybrd1 from the same
      ! guess, and the Newton again from hybrd1's result, whose own stop
      ! does not resolve a trace unknown (module header).
      ! used_newton returns .true. iff the first Newton converged (for the
      ! fallback count). The optional converged flag reports whether the
      ! returned x met the stopping test of the module header (1) or hybrd1
      ! alone reported convergence (1), or neither (0), so a caller can tell
      ! an accepted solution from a solver that gave up.
      subroutine solve_ieq(fcn, jac, n, x, params, tol, wa, lwa, used_newton,  &
                           converged)
      external :: fcn, jac, hybrd1
      integer, intent(in)    :: n, lwa
      real*8,  intent(inout) :: x(n)
      real*8,  intent(in)    :: params(*), tol
      real*8,  intent(inout) :: wa(lwa)
      logical, intent(out)   :: used_newton
      integer, intent(out), optional :: converged
      integer :: info, info_m, info_r
      real*8  :: xsave(n), fvec(n), xm(n)

      ! The EXHALE_FORCE_HYBRD1 validation switch is read ONCE, serially,
      ! by newton_solver_read_environment (input_read calls it after the
      ! keys); a caller that never called it runs with the switch off. The
      ! lazy read that stood here was removed on 2026-09-13 (review P2): it
      ! wrote two shared variables without synchronization and relied on
      ! the main program's serial first sweep to be safe.

      xsave = x
      ! solve_ieq is called from the OpenMP-parallel ionization cell sweep; the
      ! run-wide usage counters must be updated atomically.
      !$omp atomic
      nt_calls = nt_calls + 1
      if (nt_force .or. .not. use_newton_ieq) then
         info = 0                      ! legacy / validation: use the hybrd1 path
      else
         call newton_dense(fcn, jac, n, x, params, tol, info)
      endif
      if (info .eq. 1) then
         used_newton = .true.
         if (present(converged)) converged = 1
      else
         ! Newton failed -> MINPACK solve from the same guess, then the
         ! Newton from its result (kept only when it meets the test).
         x = xsave
         call hybrd1(fcn, n, x, fvec, tol, info_m, wa, lwa, params)
         used_newton = .false.
         info_r = 0
         if (.not. (nt_force .or. .not. use_newton_ieq)) then
            xm = x
            call newton_dense(fcn, jac, n, x, params, tol, info_r)
            if (info_r .ne. 1) x = xm
         endif
         if (present(converged))                                          &
            converged = merge(1, 0, info_m .eq. 1 .or. info_r .eq. 1)
         !$omp atomic
         nt_fallback = nt_fallback + 1
      endif
      end subroutine solve_ieq


      ! hybrd1, then the forward-difference Newton from its result: hybrd1
      ! stops on a step relative to the scaled norm of all the unknowns,
      ! which the O(1) fractions dominate; the Newton continues until the
      ! correction resolves every unknown (module header). The Newton's
      ! result is kept only when it meets that test, and info is then 1;
      ! otherwise x and info are hybrd1's.
      subroutine hybrd1_resolving_every_unknown(fcn, n, x, fvec, rtol,     &
                                                info, wa, lwa, params)
      external :: fcn, hybrd1
      integer, intent(in)    :: n, lwa
      real*8,  intent(inout) :: x(n), fvec(n), wa(lwa)
      real*8,  intent(in)    :: params(*), rtol
      integer, intent(out)   :: info
      real*8  :: xh(n)
      integer :: info_r, iflag

      call hybrd1(fcn, n, x, fvec, rtol, info, wa, lwa, params)
      ! The validation switch EXHALE_FORCE_HYBRD1 asks for hybrd1 alone.
      if (nt_force) return
      xh = x
      call newton_fd(fcn, n, x, params, rtol, info_r)
      if (info_r .eq. 1) then
         info = 1
      else
         x = xh
      endif
      ! fvec is the residual at the x returned, as hybrd1's contract says.
      iflag = 1
      call fcn(n, x, fvec, iflag, params)
      end subroutine hybrd1_resolving_every_unknown


      ! Newton with the analytic Jacobian. info = 1 when the correction at
      ! the returned x resolves every unknown (module header), 0 otherwise
      ! (singular Jacobian, no admissible step, NaN, or maxit reached); the
      ! caller then falls back to hybrd1.
      subroutine newton_dense(fcn, jac, n, x, params, rtol, info)
      external :: fcn, jac
      integer, intent(in)    :: n
      real*8,  intent(inout) :: x(n)
      real*8,  intent(in)    :: params(*), rtol
      integer, intent(out)   :: info
      call damped_newton_to_every_unknown(fcn, jac, .true., n, x, params,  &
                                          rtol, info)
      end subroutine newton_dense


      ! The same Newton with a forward-difference Jacobian (the systems
      ! without an analytic one: the molecular and the He 2^3S layouts).
      subroutine newton_fd(fcn, n, x, params, rtol, info)
      external :: fcn
      integer, intent(in)    :: n
      real*8,  intent(inout) :: x(n)
      real*8,  intent(in)    :: params(*), rtol
      integer, intent(out)   :: info
      call damped_newton_to_every_unknown(fcn, fcn, .false., n, x, params, &
                                          rtol, info)
      end subroutine newton_fd


      ! Damped Newton whose stop is the correction test of the module
      ! header. A step x + lam dx (lam halved from 1) is taken when it
      ! lowers the residual norm (Armijo, as the solver always did) OR the
      ! simplified Newton correction at the trial point, J^{-1} f(x + lam dx)
      ! with the Jacobian of x, measured relative to each unknown
      ! (D_i = 1/(|x_i| + fraction_resolution)), falls below (1 - lam/4)
      ! of that of the step (a monotonicity test on the error estimate
      ! itself, invariant under any scaling of the rows). The second
      ! test is what lets a trace unknown be resolved once the large rows
      ! sit at their rounding and the residual norm, which they dominate,
      ! can no longer decrease.
      subroutine damped_newton_to_every_unknown(fcn, jac, analytic, n, x,  &
                                                params, rtol, info)
      external :: fcn, jac
      logical, intent(in)    :: analytic
      integer, intent(in)    :: n
      real*8,  intent(inout) :: x(n)
      real*8,  intent(in)    :: params(*), rtol
      integer, intent(out)   :: info
      integer, parameter :: maxit = 60, maxls = 20
      real*8  :: fvec(n), fjac(n,n), dx(n), dxbar(n), xnew(n), fnew(n)
      real*8  :: dsc(n), fnorm, fnorm_new, lam, theta0, theta
      integer :: ipiv(n), it, ils, iflag, sinfo
      logical :: accepted

      info  = 0
      iflag = 1
      call fcn(n, x, fvec, iflag, params)
      fnorm = sqrt(sum(fvec*fvec))
      if (.not. (fnorm .eq. fnorm)) return          ! NaN at the guess
      if (fnorm .eq. 0.0d0) then
         info = 1; return
      endif

      do it = 1, maxit
         if (analytic) then
            call jac(n, x, fjac, params)
         else
            call forward_difference_jacobian(fcn, n, x, fvec, fjac, params)
         endif
         call lu_factor(n, fjac, ipiv, sinfo)
         if (sinfo .ne. 0) return                   ! singular -> fall back
         dx = -fvec
         call lu_solve(n, fjac, ipiv, dx)
         if (.not. all(dx .eq. dx)) return          ! NaN correction

         ! Converged: the correction resolves every unknown. It is applied
         ! when the corrected state is finite (it lies within the
         ! resolution of x either way).
         if (correction_resolves_every_unknown(n, dx, x, rtol)) then
            xnew = x + dx
            call fcn(n, xnew, fnew, iflag, params)
            if (all(fnew .eq. fnew)) x = xnew
            info = 1
            return
         endif

         dsc    = 1.0d0/(abs(x) + fraction_resolution)
         theta0 = sqrt(sum((dsc*dx)**2))
         lam = 1.0d0
         accepted = .false.
         do ils = 1, maxls
            xnew = x + lam*dx
            call fcn(n, xnew, fnew, iflag, params)
            fnorm_new = sqrt(sum(fnew*fnew))
            if (fnorm_new .eq. fnorm_new) then
               if (fnorm_new .lt. (1.0d0 - 1.0d-4*lam)*fnorm) then
                  accepted = .true.
               else
                  dxbar = -fnew
                  call lu_solve(n, fjac, ipiv, dxbar)
                  theta = sqrt(sum((dsc*dxbar)**2))
                  if (theta .eq. theta .and.                              &
                      theta .lt. (1.0d0 - 0.25d0*lam)*theta0)             &
                     accepted = .true.
               endif
            endif
            if (accepted) exit
            lam = 0.5d0*lam
         enddo
         if (.not. accepted) return                 ! no progress -> fall back

         x     = xnew
         fvec  = fnew
         fnorm = fnorm_new
         if (fnorm .eq. 0.0d0) then
            info = 1; return
         endif
      enddo
      ! maxit exhausted without meeting the test -> fall back.
      end subroutine damped_newton_to_every_unknown


      ! Forward-difference Jacobian of fcn at x (fvec = fcn(x)), with the
      ! step rule of MINPACK fdjac1: h_j = sqrt(epsilon) |x_j|, or
      ! sqrt(epsilon) at x_j = 0, so that a trace unknown is differenced at
      ! its own size; the step used is the representable (x_j + h) - x_j.
      subroutine forward_difference_jacobian(fcn, n, x, fvec, fjac, params)
      external :: fcn
      integer, intent(in)    :: n
      real*8,  intent(in)    :: x(n), fvec(n), params(*)
      real*8,  intent(out)   :: fjac(n,n)
      real*8  :: xp(n), fp(n), h
      integer :: j, iflag
      iflag = 1
      do j = 1, n
         xp = x
         h  = sqrt(epsilon(1.0d0))*abs(x(j))
         if (h .eq. 0.0d0) h = sqrt(epsilon(1.0d0))
         xp(j) = x(j) + h
         h = xp(j) - x(j)
         call fcn(n, xp, fp, iflag, params)
         fjac(:,j) = (fp - fvec)/h
      enddo
      end subroutine forward_difference_jacobian


      ! LU factorization with partial pivoting, a = P L U in place (L unit
      ! lower, below the diagonal; U on and above it). info = 1 if a is
      ! singular.
      subroutine lu_factor(n, a, ipiv, info)
      integer, intent(in)    :: n
      real*8,  intent(inout) :: a(n,n)
      integer, intent(out)   :: ipiv(n), info
      integer :: i, k, p
      real*8  :: amax, t(n)

      info = 0
      do k = 1, n
         p = k; amax = abs(a(k,k))
         do i = k+1, n
            if (abs(a(i,k)) .gt. amax) then
               amax = abs(a(i,k)); p = i
            endif
         enddo
         ipiv(k) = p
         if (.not. (amax .gt. 0.0d0)) then
            info = 1; return
         endif
         if (p .ne. k) then
            t = a(k,:); a(k,:) = a(p,:); a(p,:) = t
         endif
         do i = k+1, n
            a(i,k) = a(i,k)/a(k,k)
            a(i,k+1:n) = a(i,k+1:n) - a(i,k)*a(k,k+1:n)
         enddo
      enddo
      end subroutine lu_factor


      ! Solve a y = b with the factors of lu_factor (y returned in b).
      subroutine lu_solve(n, a, ipiv, b)
      integer, intent(in)    :: n, ipiv(n)
      real*8,  intent(in)    :: a(n,n)
      real*8,  intent(inout) :: b(n)
      integer :: i, k
      real*8  :: t
      do k = 1, n
         if (ipiv(k) .ne. k) then
            t = b(k); b(k) = b(ipiv(k)); b(ipiv(k)) = t
         endif
      enddo
      do i = 2, n
         b(i) = b(i) - sum(a(i,1:i-1)*b(1:i-1))
      enddo
      do i = n, 1, -1
         if (i .lt. n) b(i) = b(i) - sum(a(i,i+1:n)*b(i+1:n))
         b(i) = b(i)/a(i,i)
      enddo
      end subroutine lu_solve

      end module newton_solver
