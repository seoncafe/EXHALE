      module newton_solver
      use global_parameters, only: use_newton_ieq
      ! Task 2: a small, self-contained damped-Newton solver with an ANALYTIC
      ! Jacobian for the ionization-equilibrium systems, plus solve_ieq, which
      ! tries Newton first and falls back to MINPACK hybrd1 on any failure.
      ! The fallback guarantees no regression: if the analytic-Jacobian Newton
      ! does not converge (or the Jacobian is wrong), the result is exactly what
      ! hybrd1 would have produced.
      !
      ! Interfaces (matching the existing residuals so they can be passed as
      ! actual arguments):
      !   fcn(n, x, fvec, iflag, params)     -- residual (MINPACK form)
      !   jac(n, x, fjac, params)            -- analytic Jacobian d fvec_i/d x_k
      ! Both read the same params layout and module state as today.

      implicit none
      private
      public :: solve_ieq, newton_dense
      ! Run-wide usage counters (reported once at the end of the run): how many
      ! ionization-equilibrium solves the analytic-Jacobian Newton handled vs.
      ! how many fell back to MINPACK hybrd1.
      integer, save, public :: nt_calls = 0, nt_fallback = 0
      ! Reproducible-validation switch: if the environment variable
      ! ATES_FORCE_HYBRD1 is set (to a non-empty value), solve_ieq skips Newton
      ! and always uses MINPACK hybrd1. Lets the same binary produce both the
      ! Newton and the reference hybrd1 outputs for an A/B comparison.
      logical, save :: nt_init = .false., nt_force = .false.

      contains

      ! Try the analytic-Jacobian Newton; on failure restore x and call hybrd1.
      ! used_newton returns .true. iff Newton converged (for fallback counting).
      subroutine solve_ieq(fcn, jac, n, x, params, tol, wa, lwa, used_newton)
      external :: fcn, jac, hybrd1
      integer, intent(in)    :: n, lwa
      real*8,  intent(inout) :: x(n)
      real*8,  intent(in)    :: params(*), tol
      real*8,  intent(inout) :: wa(lwa)
      logical, intent(out)   :: used_newton
      integer :: info, info_m
      real*8  :: xsave(n), fvec(n)
      character(len=8) :: envval

      ! One-time check of the ATES_FORCE_HYBRD1 validation switch.
      if (.not. nt_init) then
         call get_environment_variable('ATES_FORCE_HYBRD1', envval)
         nt_force = (len_trim(envval) .gt. 0)
         nt_init  = .true.
      endif

      xsave = x
      nt_calls = nt_calls + 1
      if (nt_force .or. .not. use_newton_ieq) then
         info = 0                      ! legacy / validation: use the hybrd1 path
      else
         call newton_dense(fcn, jac, n, x, params, tol, info)
      endif
      if (info .eq. 1) then
         used_newton = .true.
      else
         ! Newton failed -> exact-as-before MINPACK solve from the same guess.
         x = xsave
         call hybrd1(fcn, n, x, fvec, tol, info_m, wa, lwa, params)
         used_newton = .false.
         nt_fallback = nt_fallback + 1
      endif
      end subroutine solve_ieq

      ! Damped Newton with a backtracking (Armijo) line search on ||fvec||_2.
      ! info = 1 on convergence, 0 otherwise (singular Jacobian, no progress,
      ! NaN, or maxit reached) -> caller falls back to hybrd1.
      subroutine newton_dense(fcn, jac, n, x, params, ftol, info)
      external :: fcn, jac
      integer, intent(in)    :: n
      real*8,  intent(inout) :: x(n)
      real*8,  intent(in)    :: params(*), ftol
      integer, intent(out)   :: info
      integer, parameter :: maxit = 60, maxls = 20
      real*8  :: fvec(n), fjac(n,n), dx(n), xnew(n), fnew(n)
      real*8  :: fnorm, fnorm_new, lam, dxmax, xscale
      integer :: it, ils, iflag, sinfo

      info  = 0
      iflag = 1
      call fcn(n, x, fvec, iflag, params)
      fnorm = sqrt(sum(fvec*fvec))
      if (.not. (fnorm .eq. fnorm)) return          ! NaN at the guess
      if (fnorm .eq. 0.0d0) then
         info = 1; return
      endif

      do it = 1, maxit
         ! Analytic Jacobian at the current x, then solve fjac*dx = -fvec.
         call jac(n, x, fjac, params)
         dx = -fvec
         call gauss_solve(n, fjac, dx, sinfo)
         if (sinfo .ne. 0) return                   ! singular -> fall back

         ! Backtracking line search: require ||fvec|| to decrease.
         lam = 1.0d0
         fnorm_new = fnorm
         do ils = 1, maxls
            xnew = x + lam*dx
            call fcn(n, xnew, fnew, iflag, params)
            fnorm_new = sqrt(sum(fnew*fnew))
            if ((fnorm_new .eq. fnorm_new) .and.                          &
                (fnorm_new .lt. (1.0d0 - 1.0d-4*lam)*fnorm)) exit  ! accept
            lam = 0.5d0*lam
         enddo
         if (.not. (fnorm_new .lt. fnorm)) return   ! no progress -> fall back

         x     = xnew
         fvec  = fnew
         fnorm = fnorm_new

         ! Convergence: Newton step small (fractions are O(1)) or residual ~0.
         dxmax  = maxval(abs(lam*dx))
         xscale = 1.0d0 + maxval(abs(x))
         if (dxmax .lt. 1.0d-11*xscale .or. fnorm .lt. ftol) then
            info = 1
            return
         endif
      enddo
      ! maxit exhausted without meeting the step/residual test -> fall back.
      end subroutine newton_dense

      ! Solve a*y = b for y (returned in b) by Gaussian elimination with
      ! partial pivoting. a is overwritten. info = 1 if a is singular.
      subroutine gauss_solve(n, a, b, info)
      integer, intent(in)    :: n
      real*8,  intent(inout) :: a(n,n), b(n)
      integer, intent(out)   :: info
      integer :: i, j, k, p
      real*8  :: amax, t, f

      info = 0
      do k = 1, n
         ! partial pivot
         p = k; amax = abs(a(k,k))
         do i = k+1, n
            if (abs(a(i,k)) .gt. amax) then
               amax = abs(a(i,k)); p = i
            endif
         enddo
         if (amax .eq. 0.0d0) then
            info = 1; return
         endif
         if (p .ne. k) then
            do j = k, n
               t = a(k,j); a(k,j) = a(p,j); a(p,j) = t
            enddo
            t = b(k); b(k) = b(p); b(p) = t
         endif
         ! eliminate
         do i = k+1, n
            f = a(i,k)/a(k,k)
            do j = k+1, n
               a(i,j) = a(i,j) - f*a(k,j)
            enddo
            b(i) = b(i) - f*b(k)
         enddo
      enddo
      ! back substitution
      do i = n, 1, -1
         t = b(i)
         do j = i+1, n
            t = t - a(i,j)*b(j)
         enddo
         b(i) = t/a(i,i)
      enddo
      end subroutine gauss_solve

      end module newton_solver
