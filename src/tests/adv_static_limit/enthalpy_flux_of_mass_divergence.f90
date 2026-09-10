      program enthalpy_flux_of_mass_divergence
      ! THE THIRD TERM OF THE ADVECTED ENERGY EQUATION (T_equation.f90).
      !
      ! The steady internal-energy equation of a supplied profile is
      !
      !     div(u v) + p div(v)  =  heating - cooling
      !   = rho v de/dr - p v dln(rho)/dr + h div(rho v) ,   h = e + p/rho
      !
      ! and the third term, the enthalpy flux carried by the divergence of
      ! the mass flux, reaches the residual through teq_cell%div_rhov. It
      ! vanishes for a stationary mass flux rho v r^2 = const and for nothing
      ! else, so the rows below pin two things: the algebra of the term in
      ! each of the two branches of the residual, and what the equation gives
      ! on a column whose mass flux is NOT stationary, where the term is the
      ! whole difference between the exact solution and the adiabat.
      !
      ! THE COLUMN WITH A KNOWN EXACT SOLUTION. Take rho constant, mu = 1,
      ! no heating and no cooling, and a mass flux F = rho v r^2 = r, that is
      ! v = 1/r. With Q = 0 the equation reduces to
      !
      !     dln(w)/dr = (gamma - 1) dln(rho)/dr - gamma dln(F)/dr ,  w = p/rho
      !
      ! so w, and with mu = 1 the temperature, follows
      ! T proportional to rho^(gamma-1) F^(-gamma), here T proportional to
      ! r^(-gamma). Without the third term the same column gives
      ! dln(w)/dr = (gamma - 1) dln(rho)/dr, the adiabat, which for a constant
      ! density is T = const: the term is not a correction to that answer, it
      ! is the answer.
      !
      ! One line per assertion:
      !   PASS|FAIL <name> measured= reference= tol=
      use global_parameters, only: T0, q0, n0, gamma_ad, thereis_He
      use ion_cell_state, only: teq_cell
      use caloric_eos, only: internal_energy_of_mixture
      use equation_T, only: Tres, solve_T_brent
      implicit none

      real*8  :: params(40)
      real*8  :: xj, res_off, res_on, added, expected
      real*8  :: div0, err(4), rate
      integer :: nfail, i
      integer, parameter :: ncell(4) = (/ 50, 100, 200, 400 /)

      nfail  = 0
      params = 0.0d0

      ! Scales of the residual. q0 divides the cooling sum, so it has to be
      ! nonzero even where every density is; T0 sets the temperature unit.
      T0 = 1.0d3
      q0 = 1.0d0
      n0 = 1.0d0
      thereis_He = .false.

      ! ---------------------------------------------------------------!
      ! Rows 1 and 2: the algebra of the term, one row per branch.
      !
      ! Dividing the residual by mup*mum*dr puts it in the form above, with
      ! e = E(x_H2,T)/mu and p/rho = T/mu, E being the energy per particle of
      ! internal_energy_of_mixture. So the term the residual has to carry is
      ! mup*mum*dr*h_j*div(rho v) = mum*dr*(E(x_H2,x) + x)*div_rhov, and in
      ! the monatomic branch, which carries the extra factor (gamma_ad - 1)
      ! of its own scaling, that collapses to gamma_ad*mum*dr*x*div_rhov.
      ! Each row forms the residual twice on the same cell, once with the
      ! field zero and once with it set, and measures the difference against
      ! that expression.
      ! ---------------------------------------------------------------!
      call one_cell(0.0d0)
      xj = 1.1d0

      res_off = Tres(xj, params)
      teq_cell%div_rhov = 0.37d0
      res_on  = Tres(xj, params)
      added    = res_on - res_off
      expected = gamma_ad*teq_cell%mum*teq_cell%dr*teq_cell%div_rhov*xj
      call chk('enthalpy_flux_term_algebra_monatomic', added, expected,     &
               1.0d-12, nfail)

      call one_cell(0.4d0)
      res_off = Tres(xj, params)
      teq_cell%div_rhov = 0.37d0
      res_on  = Tres(xj, params)
      added    = res_on - res_off
      expected = teq_cell%mum*teq_cell%dr*teq_cell%div_rhov                 &
                 *(internal_energy_of_mixture(teq_cell%x_h2, xj) + xj)
      call chk('enthalpy_flux_term_algebra_caloric', added, expected,       &
               1.0d-12, nfail)

      ! ---------------------------------------------------------------!
      ! Row 3: a stationary mass flux carries no term at all.
      !
      ! div(rho v) is formed in post_process_adv as the mass row's
      ! flux-difference operator over the control volume between the two
      ! points the upwind difference is taken between, so two points with the
      ! same rho v r^2 give exactly zero, and the residual is then the
      ! advected balance of the first two terms alone, to the bit. This is
      ! the property that leaves a converged wind, and every path that never
      ! sets the field, where it was.
      ! ---------------------------------------------------------------!
      call one_cell(0.0d0)
      div0 = stationary_divergence()
      teq_cell%div_rhov = div0
      res_on = Tres(xj, params)
      teq_cell%div_rhov = 0.0d0
      res_off = Tres(xj, params)
      call chk_exact('enthalpy_flux_absent_at_stationary_mass_flux',        &
                     res_on - res_off, nfail)
      call chk_exact('stationary_mass_flux_divergence_is_zero', div0, nfail)

      ! ---------------------------------------------------------------!
      ! Rows 4 to 6: the column with the known exact solution.
      !
      ! Row 4 marches it with the term and measures the departure from
      ! T proportional to r^(-gamma) on the finest grid; row 5 measures that
      ! the departure falls with the cell width, which is what a first-order
      ! upwind difference of this equation has to do; row 6 marches the same
      ! column with the term removed and measures the departure of the
      ! adiabat from the same exact solution, which is the size of what the
      ! term carries.
      ! ---------------------------------------------------------------!
      do i = 1, 4
         err(i) = column_departure(ncell(i), .true.)
      enddo
      call chk_lt('exact_equation_solution_on_a_diverging_column',          &
                  err(4), 1.0d-2, nfail)
      ! First order halves the departure when the cell width is halved, so
      ! two refinements quarter it; the bound admits anything at least first
      ! order and refuses a scheme that has stopped converging.
      rate = err(4)/err(2)
      call chk_lt('exact_equation_solution_converges_with_cell_width',      &
                  rate, 0.6d0, nfail)
      write(*,'(a,4(es11.4,1x))') '  departure at 50, 100, 200, 400'//     &
         ' cells: ', err(1), err(2), err(3), err(4)
      ! The adiabat of a constant density is T = const, and the exact
      ! solution falls as r^(-gamma), so over r = 1 to 2 the two part by
      ! 2^gamma - 1 = 2.17.
      call chk_gt('adiabat_departs_from_the_exact_solution',                &
                  column_departure(ncell(4), .false.), 1.0d0, nfail)

      if (nfail .gt. 0) stop 1

      contains

      ! One cell of a wind, in the units of the residual. The numbers are a
      ! cell of no particular run: the rows that use it measure a difference
      ! of the residual with itself, so only their scale matters.
      subroutine one_cell(x_h2)
      real*8, intent(in) :: x_h2
      teq_cell%nhi    = 0.0d0
      teq_cell%nhii   = 0.0d0
      teq_cell%nhei   = 0.0d0
      teq_cell%nheii  = 0.0d0
      teq_cell%nheiii = 0.0d0
      teq_cell%mup    = 1.2d0
      teq_cell%mum    = 1.3d0
      teq_cell%rhov   = 2.5d0
      teq_cell%coeff  = -0.4d0
      teq_cell%dr     = 3.0d-3
      teq_cell%Told   = 1.0d0
      teq_cell%heaold = 0.0d0
      teq_cell%x_h2   = x_h2
      teq_cell%div_rhov = 0.0d0
      end subroutine one_cell

      ! The divergence post_process_adv forms, on two points that carry the
      ! same mass flux: (A_p F_p - A_m F_m)/dV with F = rho v, A = r^2 and
      ! dV = d(r^3)/3. The two states are chosen so that rho v r^2 is the
      ! same BINARY number at both points and not merely the same to
      ! round-off: r_hi/r_lo = 5/4 and rho_lo = 25/16 are exact in binary, so
      ! both products evaluate to 25/16 with no rounding, and the difference
      ! is an exact zero. That is the statement being made, that stationarity
      ! removes the term exactly rather than nearly.
      real*8 function stationary_divergence()
      real*8 :: rlo, rhi, rho_lo, rho_hi, v_lo, v_hi
      rlo    = 1.0d0
      rhi    = 1.25d0
      rho_lo = 1.5625d0        ! = (5/4)^2
      v_lo   = 1.0d0
      rho_hi = 1.0d0
      v_hi   = 1.0d0
      stationary_divergence = (rho_hi*v_hi*rhi*rhi - rho_lo*v_lo*rlo*rlo)  &
                              /((rhi**3 - rlo**3)/3.0d0)
      end function stationary_divergence

      ! March the column rho = 1, v = 1/r, mu = 1, no heating and no cooling
      ! from r = 1 to r = 2 on n cells and return the largest relative
      ! departure of the marched temperature from the exact solution
      ! T proportional to r^(-gamma). with_term = .false. removes the
      ! enthalpy flux and marches the adiabat instead.
      real*8 function column_departure(n, with_term)
      integer, intent(in) :: n
      logical, intent(in) :: with_term
      real*8  :: rlo, rhi, dr, xprev, xnew, xexact, d
      integer :: j
      logical :: ok

      dr = 1.0d0/dble(n)
      xprev = 1.0d0                     ! T at r = 1, in units of T0
      column_departure = 0.0d0
      do j = 1, n
         rlo = 1.0d0 + dble(j-1)*dr
         rhi = 1.0d0 + dble(j)*dr
         teq_cell%nhi    = 0.0d0
         teq_cell%nhii   = 0.0d0
         teq_cell%nhei   = 0.0d0
         teq_cell%nheii  = 0.0d0
         teq_cell%nheiii = 0.0d0
         teq_cell%mup    = 1.0d0
         teq_cell%mum    = 1.0d0
         teq_cell%rhov   = 1.0d0/rhi    ! rho v with rho = 1, v = 1/r
         teq_cell%coeff  = 0.0d0        ! constant density
         teq_cell%dr     = rhi - rlo
         teq_cell%Told   = xprev
         teq_cell%heaold = 0.0d0
         teq_cell%x_h2   = 0.0d0
         if (with_term) then
            ! F = rho v r^2 = r at both points, divided by the volume between
            ! them, exactly as post_process_adv forms it.
            teq_cell%div_rhov = (rhi - rlo)/((rhi**3 - rlo**3)/3.0d0)
         else
            teq_cell%div_rhov = 0.0d0
         endif
         call solve_T_brent(params, xprev, xnew, ok)
         if (.not. ok) then
            column_departure = huge(1.0d0)
            return
         endif
         xprev  = xnew
         xexact = rhi**(-gamma_ad)
         d = abs(xnew - xexact)/xexact
         if (d .gt. column_departure) column_departure = d
      enddo
      end function column_departure

      subroutine chk(name, measured, reference, tol, nf)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, reference, tol
      integer, intent(inout) :: nf
      character(len=4) :: verdict
      verdict = 'PASS'
      if (.not. (abs(measured - reference)                                 &
                 .le. tol*max(abs(reference),1.0d0))) then
         verdict = 'FAIL';  nf = nf + 1
      endif
      write(*,'(a,1x,a,a,es13.6,a,es13.6,a,es9.2)') trim(verdict),         &
         trim(name), ' measured=', measured, ' reference=', reference,     &
         ' tol=', tol
      end subroutine chk

      subroutine chk_exact(name, measured, nf)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured
      integer, intent(inout) :: nf
      character(len=4) :: verdict
      verdict = 'PASS'
      if (measured .ne. 0.0d0) then
         verdict = 'FAIL';  nf = nf + 1
      endif
      write(*,'(a,1x,a,a,es13.6,a)') trim(verdict), trim(name),            &
         ' measured=', measured, ' reference=0 tol=0'
      end subroutine chk_exact

      subroutine chk_gt(name, measured, bound, nf)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, bound
      integer, intent(inout) :: nf
      character(len=4) :: verdict
      verdict = 'PASS'
      if (.not. (measured .gt. bound)) then
         verdict = 'FAIL';  nf = nf + 1
      endif
      write(*,'(a,1x,a,a,es13.6,a,es13.6,a)') trim(verdict), trim(name),   &
         ' measured=', measured, ' reference=>', bound, ' tol=0'
      end subroutine chk_gt

      subroutine chk_lt(name, measured, bound, nf)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, bound
      integer, intent(inout) :: nf
      character(len=4) :: verdict
      verdict = 'PASS'
      if (.not. (measured .lt. bound)) then
         verdict = 'FAIL';  nf = nf + 1
      endif
      write(*,'(a,1x,a,a,es13.6,a,es13.6,a)') trim(verdict), trim(name),   &
         ' measured=', measured, ' reference=<', bound, ' tol=0'
      end subroutine chk_lt

      end program enthalpy_flux_of_mass_divergence
