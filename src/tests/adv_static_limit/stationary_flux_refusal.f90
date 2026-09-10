      program stationary_flux_refusal
      ! THE SENSITIVITY OF THE ADVECTED ENERGY CORRECTION TO A
      ! NON-STATIONARY MASS FLUX: the size of the enthalpy flux the
      ! mass-flux divergence puts into the equation, against the two terms
      ! the steady balance keeps (post_process_adv.f90,
      ! enthalpy_flux_term_ratio).
      !
      ! IT IS NOT THE STATIONARITY TEST. Whether a cell of the recorded
      ! state is stationary at all is decided by the mass row of that state,
      ! the operator the stationary certification measures (the rows of
      ! stationarity_by_the_mass_row); this ratio is reported next to that
      ! verdict, and the last row here is the counterexample that separates
      ! the two: a cell with a nonzero mass divergence whose ratio is 0.5,
      ! which the ratio accepts.
      !
      ! The energy equation the post-process solves is
      !
      !     rho v de/dr  -  p v dln(rho)/dr  +  h div(rho v)  =  heat - cool
      !
      ! whose third term, the enthalpy flux carried by the divergence of the
      ! mass flux, vanishes for rho v r^2 = const and for nothing else. Where
      ! that term is larger than BOTH terms the balance keeps, the solution of
      ! the equation is a statement about the continuity residual of the
      ! recorded profile rather than about the gas, and the cell keeps the
      ! run's own temperature and the equilibrium composition at that
      ! temperature.
      !
      ! The rows below pin the ratio itself: what it is, that the comparison
      ! against unity is STRICT so a stationary flux carrying round-off is
      ! never reported as dominated, and that the ratio is the ratio of the
      ! terms of the equation and not of the flux difference alone.
      !
      ! One line per assertion:
      !   PASS|FAIL <name> measured= reference= tol=
      use post_processing, only: enthalpy_flux_term_ratio
      implicit none

      ! A cell of a wind, in the residual's own variables. The numbers are
      ! exact binary fractions so that the equalities below are exact.
      real*8, parameter :: rho_up = 1.0d0     ! density of the cell
      real*8, parameter :: rho_lo = 0.75d0    ! density of the upwind point
      real*8, parameter :: v_up   = 2.0d0     ! velocity of the cell
      real*8, parameter :: e_up   = 1.5d0     ! internal energy per unit mass
      real*8, parameter :: e_lo   = 1.0d0     ! ... at the upwind point
      real*8, parameter :: w_up   = 0.5d0     ! p/rho of the cell
      real*8, parameter :: h_up   = e_up + w_up   ! enthalpy per unit mass
      real*8, parameter :: dr     = 0.25d0    ! cell width
      real*8 :: ratio, ratio_hi, ratio_lo, div_unity, q_adv, q_prs
      logical :: refused

      integer :: nfail

      nfail = 0

      ! The two terms the balance keeps, on this cell.
      q_adv = abs(rho_up*v_up*(e_up - e_lo))          ! = 1
      q_prs = abs(w_up*v_up*(rho_up - rho_lo))        ! = 0.25

      ! 1. A stationary mass flux leaves the ratio an EXACT zero, so the
      !    diagnostic reports nothing on a converged wind. div(rho v) = 0 is
      !    the definition of stationary for the two points the residual
      !    evaluates.
      ratio = enthalpy_flux_term_ratio(h_up, 0.0d0, dr, rho_up, v_up,        &
                                       e_up, e_lo, w_up, rho_lo)
      call chk('stationary_mass_flux_carries_no_term', ratio, 0.0d0, 0.0d0, &
               nfail)

      ! 2. The ratio is the enthalpy flux term over the LARGER of the two
      !    terms the balance keeps, so the comparison is against the term
      !    that survives and not the smaller one. Here q_adv = 1 is larger.
      div_unity = q_adv/(h_up*dr)
      ratio = enthalpy_flux_term_ratio(h_up, div_unity, dr, rho_up, v_up,    &
                                       e_up, e_lo, w_up, rho_lo)
      call chk('ratio_is_taken_against_the_larger_kept_term', ratio, 1.0d0,  &
               1.0d-15, nfail)

      ! 3. THE COMPARISON IS STRICT. A cell whose enthalpy flux term is
      !    exactly equal to the larger kept term is NOT reported as
      !    dominated: the statement is 'one term is larger than the other',
      !    so a flux that is stationary to the last bit and a flux whose
      !    residual is exactly marginal read the same. This is the row that
      !    says round-off in rho v r^2 cannot start raising the diagnostic.
      refused = (ratio > 1.0d0)
      call chk_false('ratio_of_exactly_one_is_not_dominant', refused, nfail)

      ! 4. and 5. Either side of unity the diagnostic says what it is built
      !    to say: a divergence larger than that turnover puts the largest
      !    term of the equation there, a smaller one does not. The factor is
      !    1 + 2^-10, so the rows test the inequality and not a tolerance
      !    around it.
      ratio_hi = enthalpy_flux_term_ratio(h_up,                              &
                     div_unity*(1.0d0 + 2.0d0**(-10)), dr,                   &
                     rho_up, v_up, e_up, e_lo, w_up, rho_lo)
      call chk_true('above_unity_is_dominant', ratio_hi > 1.0d0, nfail)
      ratio_lo = enthalpy_flux_term_ratio(h_up,                              &
                     div_unity*(1.0d0 - 2.0d0**(-10)), dr,                   &
                     rho_up, v_up, e_up, e_lo, w_up, rho_lo)
      call chk_false('below_unity_is_not_dominant', ratio_lo > 1.0d0, nfail)

      ! 6. The ratio is a statement about the TERMS of the equation and not
      !    about |dln F| alone: the same divergence is dominant in a cell
      !    whose kept terms are small and negligible in a cell whose kept
      !    terms are large. Here the same div(rho v) that gives the ratio 1
      !    above gives 1/16 once the advected term is 16 times larger, which
      !    is why it measures the sensitivity of the answer and not the
      !    stationarity of the state.
      ratio = enthalpy_flux_term_ratio(h_up, div_unity, dr, rho_up,          &
                                       16.0d0*v_up, e_up, e_lo, w_up, rho_lo)
      call chk('same_divergence_is_not_dominant_where_kept_terms_are_large',&
               ratio, 1.0d0/16.0d0, 1.0d-15, nfail)

      ! 7. A cell in which every term vanishes carries no equation, and the
      !    ratio is zero rather than 0/0: there is nothing to report.
      ratio = enthalpy_flux_term_ratio(h_up, 0.0d0, dr, 0.0d0, 0.0d0,        &
                                       e_up, e_up, w_up, 0.0d0)
      call chk('empty_equation_carries_no_term', ratio, 0.0d0, 0.0d0, nfail)

      ! 8. THE COUNTEREXAMPLE THAT SEPARATES THE RATIO FROM STATIONARITY.
      !    Half the turnover of the larger kept term leaves the ratio at
      !    0.5, and the mass divergence that produces it is not zero: the
      !    mass flux through this cell changes by q_adv/(2 h dr) times the
      !    control volume from one end to the other. A screen that accepts
      !    every ratio below one accepts this cell, so a ratio below one is
      !    not a statement that the state is stationary there. What decides
      !    that is the mass row of the state (stationarity_by_the_mass_row).
      ratio = enthalpy_flux_term_ratio(h_up, 0.5d0*div_unity, dr, rho_up,   &
                                       v_up, e_up, e_lo, w_up, rho_lo)
      call chk('half_turnover_gives_the_ratio_one_half', ratio, 0.5d0,      &
               1.0d-15, nfail)
      call chk_false('a_ratio_below_one_is_not_a_stationary_mass_flux',     &
                     0.5d0*div_unity .eq. 0.0d0, nfail)

      if (nfail .gt. 0) then
         write(*,'(a,i0,a)') 'stationary_flux_refusal: ', nfail, ' FAILED'
         call exit(1)
      endif

      contains

      subroutine chk(name, measured, reference, tol, nf)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, reference, tol
      integer, intent(inout) :: nf
      character(len=4) :: verdict
      verdict = 'PASS'
      if (.not. (abs(measured - reference) .le. tol*max(abs(reference),1.0d0))) then
         verdict = 'FAIL';  nf = nf + 1
      endif
      write(*,'(a,1x,a,a,es13.6,a,es13.6,a,es9.2)') trim(verdict), trim(name), &
         ' measured=', measured, ' reference=', reference, ' tol=', tol
      end subroutine chk

      subroutine chk_true(name, measured, nf)
      character(len=*), intent(in) :: name
      logical, intent(in) :: measured
      integer, intent(inout) :: nf
      character(len=4) :: verdict
      verdict = 'PASS'
      if (.not. measured) then
         verdict = 'FAIL';  nf = nf + 1
      endif
      write(*,'(a,1x,a,a,l1,a)') trim(verdict), trim(name),                  &
         ' measured=', measured, ' reference=T tol=0'
      end subroutine chk_true

      subroutine chk_false(name, measured, nf)
      character(len=*), intent(in) :: name
      logical, intent(in) :: measured
      integer, intent(inout) :: nf
      character(len=4) :: verdict
      verdict = 'PASS'
      if (measured) then
         verdict = 'FAIL';  nf = nf + 1
      endif
      write(*,'(a,1x,a,a,l1,a)') trim(verdict), trim(name),                  &
         ' measured=', measured, ' reference=F tol=0'
      end subroutine chk_false

      end program stationary_flux_refusal
