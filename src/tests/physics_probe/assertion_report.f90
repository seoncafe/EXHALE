      module assertion_report
      ! One printed verdict line per assertion, and a failure count the
      ! calling program turns into its exit status.
      !
      ! Line format (the aggregation convention of the Phase 0 test
      ! directories):
      !     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
      !
      ! A comparison is written so that a NaN FAILS: every verdict is
      ! reached through .not.(deviation .le. tolerance), never through
      ! deviation .gt. tolerance.

      implicit none
      private
      public :: check_relative, check_absolute, check_within_factor,      &
                check_positive, check_at_most, check_at_least,            &
                assertion_failures

      integer, save :: assertion_failures = 0

      contains

      !--------------!

      subroutine verdict(name, ok, measured, reference, tol)
      character(len=*), intent(in) :: name, reference, tol
      logical, intent(in) :: ok
      real*8, intent(in)  :: measured
      character(len=8) :: tag
      if (ok) then
         tag = 'PASS'
      else
         tag = 'FAIL'
         assertion_failures = assertion_failures + 1
      endif
      write(*,'(a,1x,a,1x,a,es23.15,1x,a,1x,a)')                          &
           trim(tag), trim(name), 'measured=', measured,                  &
           'reference='//trim(reference), 'tol='//trim(tol)
      end subroutine verdict

      !--------------!

      character(len=24) function number(x)
      real*8, intent(in) :: x
      write(number,'(es23.15)') x
      number = adjustl(number)
      end function number

      !--------------!

      subroutine check_relative(name, measured, reference, tol)
      ! |measured - reference| <= tol * |reference|, reference nonzero.
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, reference, tol
      logical :: ok
      ok = abs(measured - reference) .le. tol*abs(reference)
      call verdict(name, ok, measured, number(reference),                 &
                   trim(number(tol))//' relative')
      end subroutine check_relative

      !--------------!

      subroutine check_absolute(name, measured, reference, tol)
      ! |measured - reference| <= tol, for a quantity that is already
      ! dimensionless (a relative jump, a closure residual).
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, reference, tol
      logical :: ok
      ok = abs(measured - reference) .le. tol
      call verdict(name, ok, measured, number(reference),                 &
                   trim(number(tol))//' absolute')
      end subroutine check_absolute

      !--------------!

      subroutine check_within_factor(name, measured, reference, factor)
      ! reference/factor <= measured <= reference*factor, both positive.
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, reference, factor
      logical :: ok
      character(len=32) :: tol
      ok = (measured .ge. reference/factor) .and.                         &
           (measured .le. reference*factor)
      write(tol,'(f0.3)') factor
      call verdict(name, ok, measured, number(reference),                 &
                   'factor '//trim(adjustl(tol)))
      end subroutine check_within_factor

      !--------------!

      subroutine check_at_most(name, measured, bound)
      ! measured <= bound, for a count or a cost with an upper bound.
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, bound
      call verdict(name, measured .le. bound, measured, number(bound),    &
                   'at most')
      end subroutine check_at_most

      !--------------!

      subroutine check_at_least(name, measured, bound)
      ! measured >= bound, for a quantity a test needs to stay large enough
      ! to keep the comparison it is part of discriminating.
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, bound
      call verdict(name, measured .ge. bound, measured, number(bound),    &
                   'at least')
      end subroutine check_at_least

      !--------------!

      subroutine check_positive(name, measured)
      ! measured > 0 and not a NaN.
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured
      call verdict(name, measured .gt. 0.0d0, measured, '0',              &
                   'strictly positive')
      end subroutine check_positive

      end module assertion_report
