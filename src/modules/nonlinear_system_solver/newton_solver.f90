      module newton_solver
      use global_parameters,   only: use_newton_ieq
      use chemical_row_layout, only: fraction_resolution, row_identity,    &
                                     row_pinned_fraction,                  &
                                     row_species_balance,                  &
                                     row_stage_boundary_flow,              &
                                     row_summed_he_i_balance,              &
                                     row_imposed_reservoir
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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
      ! A norm-based stopping test can leave trace fractions unresolved.
      ! See md/Update_EXHALE_stage3.md sections 116.2, 116.6, and 118.1.
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
                newton_solver_read_environment,                           &
                composition_with_prescribed_fractions
      ! The cell solve of composition_with_prescribed_fractions solves the
      ! free fractions only: the identity and pinned rows of a cell are data
      ! of its chemical equations, not unknowns (that routine's header).
      logical, parameter :: composition_fraction_reduction_on = .true.
      ! Run-wide usage counters (reported once at the end of the run): how many
      ! ionization-equilibrium solves the analytic-Jacobian Newton handled vs.
      ! how many fell back to MINPACK hybrd1.
      integer, save, public :: nt_calls = 0, nt_fallback = 0
      ! EXHALE_FORCE_HYBRD1: any nonempty value selects MINPACK hybrd1
      ! for every cell solve; otherwise use_newton_ieq selects the solver.
      logical, save :: nt_force = .false.

      abstract interface
         subroutine chemical_residual_function(n, x, fvec, flag, params)
         integer :: n, flag
         real*8  :: x(n), fvec(n), params(*)
         end subroutine chemical_residual_function
      end interface
      ! The prescribed composition and free rows of the active cell solve.
      ! Pointers refer to that invocation's local arrays only until its
      ! reduced solve returns. Each OpenMP thread owns its mapping and
      ! residual procedure; no callback uses another thread's cell state.
      logical, save :: composition_mapping_active = .false.
      integer, save :: chemical_composition_size = 0
      integer, pointer, save :: chemical_free_rows(:) => null()
      real*8, pointer, save :: chemical_prescribed_fractions(:) => null()
      real*8, pointer, save :: chemical_trial_fractions(:) => null()
      real*8, pointer, save :: chemical_trial_residual(:) => null()
      procedure(chemical_residual_function), pointer, save ::             &
                                     chemical_full_residual => null()
      !$omp threadprivate(composition_mapping_active,                    &
      !$omp chemical_composition_size, chemical_free_rows,                &
      !$omp chemical_prescribed_fractions, chemical_trial_fractions,       &
      !$omp chemical_trial_residual, chemical_full_residual)

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

      ! input_read calls newton_solver_read_environment once, before the
      ! OpenMP cell sweep. Without that call, EXHALE_FORCE_HYBRD1 is off.

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
      real*8, external :: enorm
      logical :: accepted

      info  = 0
      iflag = 1
      if (.not. all(ieee_is_finite(x))) return
      if (.not. ieee_is_finite(rtol) .or. rtol .lt. 0.0d0) return
      call fcn(n, x, fvec, iflag, params)
      if (.not. all(ieee_is_finite(fvec))) return
      ! A squared norm may underflow for a nonzero residual.
      fnorm = enorm(n, fvec)
      if (.not. ieee_is_finite(fnorm)) return
      if (all(fvec .eq. 0.0d0)) then
         info = 1; return
      endif

      do it = 1, maxit
         if (analytic) then
            call jac(n, x, fjac, params)
         else
            call forward_difference_jacobian(fcn, n, x, fvec, fjac, params)
         endif
         if (.not. all(ieee_is_finite(fjac))) return
         call lu_factor(n, fjac, ipiv, sinfo)
         if (sinfo .ne. 0) return                   ! singular -> fall back
         dx = -fvec
         call lu_solve(n, fjac, ipiv, dx)
         if (.not. all(ieee_is_finite(dx))) return   ! nonfinite correction

         ! Converged: the correction resolves every unknown. It is applied
         ! when the corrected state is finite (it lies within the
         ! resolution of x either way).
         if (correction_resolves_every_unknown(n, dx, x, rtol)) then
            xnew = x + dx
            if (.not. all(ieee_is_finite(xnew))) return
            call fcn(n, xnew, fnew, iflag, params)
            if (.not. all(ieee_is_finite(fnew))) return
            x = xnew
            info = 1
            return
         endif

         dsc    = 1.0d0/(abs(x) + fraction_resolution)
         theta0 = sqrt(sum((dsc*dx)**2))
         lam = 1.0d0
         accepted = .false.
         do ils = 1, maxls
            xnew = x + lam*dx
            if (.not. all(ieee_is_finite(xnew))) then
               lam = 0.5d0*lam
               cycle
            endif
            call fcn(n, xnew, fnew, iflag, params)
            if (.not. all(ieee_is_finite(fnew))) then
               lam = 0.5d0*lam
               cycle
            endif
            fnorm_new = enorm(n, fnew)
            if (ieee_is_finite(fnorm_new)) then
               if (fnorm_new .lt. (1.0d0 - 1.0d-4*lam)*fnorm) then
                  accepted = .true.
               else
                  dxbar = -fnew
                  call lu_solve(n, fjac, ipiv, dxbar)
                  theta = sqrt(sum((dsc*dxbar)**2))
                  if (all(ieee_is_finite(dxbar)) .and.                   &
                      ieee_is_finite(theta) .and.                         &
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
         if (all(fvec .eq. 0.0d0)) then
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


      ! THE CELL SOLVE OF A COMPOSITION WITH PRESCRIBED FRACTIONS: only the
      ! fractions that are unknowns of this cell are solved for.
      !
      ! row_kinds is the cell's row layout (chemical_row_layout,
      ! ionization_row_layout and mark_constraint_rows_of_cell). Two kinds
      ! of row state a value, not an equation of the cell:
      !   row_identity          x = 0 of an absent element or of the unused
      !                         X++ of a two-stage element; its datum is 0;
      !   row_pinned_fraction   x - x_fix of a fraction the flow carries;
      !                         its datum is x_fix, prescribed_fractions(i)
      !                         (impose_transported_fractions).
      ! Every other kind is an unknown: the kinetic rows and the imposed
      ! reservoir row of a lower ghost, x(H2) - x2 (1 - x_ion) + x_O, which
      ! couples H2 to the ionized fraction of the same solve.
      !
      ! WHY IT IS EXACT. Let C be the prescribed rows with data a and U the
      ! free rows; X(y) is the full composition with X_C = a and X_U = y.
      ! The residual fcn writes each pinned row as x_i - a_i after every
      ! other term, and each identity row as x_i, to which no reaction adds
      ! anything where the element density is zero or the stage is unused
      ! (metal_rows, cx_add_to_fvec), so F_C(X(y)) = 0 identically, and
      ! F(X(y)) = 0 is the same system as G(y) = F_U(X(y)) = 0. G is solved
      ! with the complete residual: the full composition is rebuilt before
      ! EVERY evaluation, so every reaction, the electron density and every
      ! coupling to a prescribed species stay in the free rows. The Newton
      ! correction of the full system has d_C = 0, and its free part solves
      ! J_UU d_U = -F_U, the Jacobian of G; the stopping test of the module
      ! header is applied to every free fraction, and is met by the
      ! prescribed ones identically. What changes is the iteration path:
      ! hybrd1's unit-scaled norm, trust region and evaluation budget
      ! 200 (n_free + 1) are those of the smaller system.
      !
      ! THE RETURNED STATE is the full composition, prescribed rows at
      ! their data, and residual_full is fcn at it. A prescribed row that
      ! is not exactly zero there (a layout that does not describe the
      ! residual), a value that is not finite, or a negative callback flag
      ! marks the solve unconverged (solve_status 0, or the flag), as does
      ! a datum that is not finite, which returns the entry composition
      ! with residual_full = huge. The molecular caller (ioniz_eq) counts
      ! such a status 0 in the info = 0 bin of its hybrd1 ledger, beside
      ! MINPACK's own code 0, and judges the returned composition like any
      ! unconverged attempt. A row kind outside the list above, a
      ! layout with no prescribed row, and the EXHALE_FORCE_HYBRD1
      ! validation switch take the full solve of
      ! hybrd1_resolving_every_unknown. The mapping, the full state and the
      ! work arrays are local to the call. The module residual reads them
      ! through a threadprivate mapping installed for the reduced solve and
      ! cleared on return. Its procedure address needs no stack trampoline.
      ! The cell sweep calls synchronously on the owning OpenMP thread.
      ! Nested calls on that same thread are not supported: the entry guard
      ! stops before an active mapping can be overwritten, including calls
      ! from the caller's residual. Sequential calls start with no mapping.
      recursive subroutine composition_with_prescribed_fractions(fcn,     &
                                       n_full,                            &
                                       fractions, residual_full,           &
                                       tolerance, solve_status, workspace, &
                                       workspace_size, parameters,         &
                                       row_kinds, prescribed_fractions)
      procedure(chemical_residual_function) :: fcn
      integer, intent(in)    :: n_full, workspace_size, row_kinds(n_full)
      real*8,  intent(inout) :: fractions(n_full), residual_full(n_full)
      real*8,  intent(inout) :: workspace(workspace_size)
      real*8,  intent(in)    :: tolerance, parameters(*)
      real*8,  intent(in)    :: prescribed_fractions(n_full)
      integer, intent(out)   :: solve_status
      integer, target :: free_rows(n_full)
      integer :: n_free, ir_full, ir_free, callback_flag
      real*8, target :: composition_fixed(n_full), composition_trial(n_full)
      real*8, target :: residual_trial(n_full)
      real*8  :: fractions_free(n_full), residual_free(n_full)
      ! hybrd1's work array, of the size it needs for n_full unknowns.
      real*8  :: workspace_free(n_full*(3*n_full+13)/2)
      logical :: classified

      if (composition_mapping_active)                                    &
         error stop 'COMPOSITION SOLVE: nested invocation on one thread'
      composition_mapping_active = .true.

      ! The prescribed rows take their data; the free rows are listed.
      classified = composition_fraction_reduction_on .and. .not. nt_force
      n_free = 0
      composition_fixed = prescribed_fractions
      do ir_full = 1, n_full
         select case (row_kinds(ir_full))
         case (row_identity)
            composition_fixed(ir_full) = 0.0d0
         case (row_pinned_fraction)
            continue
         case (row_species_balance, row_stage_boundary_flow,               &
               row_summed_he_i_balance, row_imposed_reservoir)
            n_free = n_free + 1
            free_rows(n_free) = ir_full
            fractions_free(n_free) = fractions(ir_full)
         case default
            classified = .false.
         end select
      enddo
      ! Nothing to prescribe, or a layout this routine does not know: the
      ! full solve.
      if (.not. classified .or. n_free .eq. n_full) then
         call hybrd1_resolving_every_unknown(fcn, n_full, fractions,       &
                                             residual_full, tolerance,     &
                                             solve_status, workspace,      &
                                             workspace_size, parameters)
         composition_mapping_active = .false.
         return
      endif
      if (.not. all(ieee_is_finite(composition_fixed))) then
         solve_status = 0
         residual_full = huge(1.0d0)
         composition_mapping_active = .false.
         return
      endif

      ! The free fractions; with none, the prescribed state is the answer
      ! and the residual below is its test.
      if (n_free .gt. 0) then
         chemical_composition_size = n_full
         chemical_free_rows => free_rows(1:n_free)
         chemical_prescribed_fractions => composition_fixed
         chemical_trial_fractions => composition_trial
         chemical_trial_residual => residual_trial
         chemical_full_residual => fcn
         call hybrd1_resolving_every_unknown(                              &
                    chemical_residual_of_free_fractions, n_free,           &
                    fractions_free(1:n_free), residual_free(1:n_free),     &
                    tolerance, solve_status, workspace_free,               &
                    size(workspace_free), parameters)
         nullify(chemical_free_rows, chemical_prescribed_fractions,       &
                 chemical_trial_fractions, chemical_trial_residual,       &
                 chemical_full_residual)
         chemical_composition_size = 0
      else
         solve_status = 1
      endif

      ! The full composition and its residual.
      fractions = composition_fixed
      do ir_free = 1, n_free
         fractions(free_rows(ir_free)) = fractions_free(ir_free)
      enddo
      callback_flag = 1
      call fcn(n_full, fractions, residual_full, callback_flag, parameters)
      if (callback_flag .lt. 0) solve_status = callback_flag
      if (.not. all(ieee_is_finite(fractions)) .or.                        &
          .not. all(ieee_is_finite(residual_full))) solve_status = 0
      ! Every prescribed row must be exactly zero at its datum: this checks
      ! the layout handed in against the rows the residual writes.
      do ir_full = 1, n_full
         if (row_kinds(ir_full) .eq. row_identity .or.                     &
             row_kinds(ir_full) .eq. row_pinned_fraction) then
            if (residual_full(ir_full) .ne. 0.0d0) solve_status = 0
         endif
      enddo

      composition_mapping_active = .false.
      end subroutine composition_with_prescribed_fractions


      ! The residual of the free fractions: the full composition at the
      ! prescribed data and the free values, evaluated by the calling
      ! thread's full residual, and its free rows.
      subroutine chemical_residual_of_free_fractions(n_unknown, unknowns, &
                                                     residual_unknown,    &
                                                     flag, parameter_data)
      integer :: n_unknown, flag
      real*8  :: unknowns(n_unknown), residual_unknown(n_unknown)
      real*8  :: parameter_data(*)
      integer :: ir_unknown
      if (.not. composition_mapping_active .or.                          &
          .not. associated(chemical_full_residual))                       &
         error stop 'COMPOSITION RESIDUAL: no active cell mapping'
      if (n_unknown .ne. size(chemical_free_rows))                         &
         error stop 'COMPOSITION RESIDUAL: free-row count mismatch'
      chemical_trial_fractions = chemical_prescribed_fractions
      do ir_unknown = 1, n_unknown
         chemical_trial_fractions(chemical_free_rows(ir_unknown)) =       &
                                                   unknowns(ir_unknown)
      enddo
      call chemical_full_residual(chemical_composition_size,             &
                            chemical_trial_fractions,                     &
                            chemical_trial_residual, flag, parameter_data)
      do ir_unknown = 1, n_unknown
         residual_unknown(ir_unknown) =                                   &
                   chemical_trial_residual(chemical_free_rows(ir_unknown))
      enddo
      end subroutine chemical_residual_of_free_fractions

      end module newton_solver
