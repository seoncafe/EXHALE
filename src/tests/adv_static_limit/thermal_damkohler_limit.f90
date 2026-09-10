      program thermal_damkohler_limit
      ! The validity condition of the advection-corrected ENERGY solve, as a
      ! function of the flow speed (post_process_adv.f90,
      ! thermal_damkohler_number).
      !
      ! The corrected temperature solves the steady advected energy balance of
      ! one cell, upwind. Whether that balance says anything about the
      ! temperature depends on how the time the gas spends in the cell compares
      ! with the time the local net radiative rate needs to rewrite the gas
      ! internal energy: the Damkohler number
      !
      !     Da = t_cross * |heating - cooling| / u_th ,   t_cross = dr/v
      !
      ! Da << 1 the temperature is carried by the flow, Da >> 1 it is the root
      ! of the local balance and what the gas carried in has been forgotten.
      ! The rows below are the statements the correction rests on, and in
      ! particular that Da grows without bound as v -> 0, so the correction
      ! switches itself off continuously in that limit instead of switching on
      ! the sign of v.
      !
      ! One line per assertion:
      !   PASS|FAIL <name> measured= reference= tol=
      use post_processing, only: thermal_damkohler_number
      implicit none

      real*8, parameter :: dr_cm = 1.0d7      ! cell width [cm]
      real*8, parameter :: cs    = 3.5d5      ! sound speed at 1140 K [cm/s]
      real*8, parameter :: q_rad = 1.0d-12    ! |heat - cool| [erg cm^-3 s^-1]
      real*8, parameter :: u_th  = 1.0d-4     ! internal energy [erg cm^-3]
      ! The two slowest outflowing cells at the base of the converged wind of
      ! backup/regression/wasp_full, MEASURED 2026-09-08 by post-processing
      ! its recorded state. r = 1.0327 Rp carries v = 28.3 cm/s and sits at
      ! Da = 0.95, just inside the advected regime, and it is the largest
      ! Damkohler number any tested cell of that wind reaches. The cell below
      ! it, r = 1.0323 Rp with v = 10.0 cm/s, is at Da = 2.7 and is in the
      ! locally balanced regime; in the run it never reaches this condition,
      ! because its inner face carries inflow (v = -8.1 cm/s) and the upwind
      ! branch has pinned it already. The two rows that use these numbers say
      ! the condition sits at the edge of a real wind rather than far from it,
      ! so it is a guard on the regime and not a switch that removes the
      ! correction.
      real*8, parameter :: t_cross_slow = 1.9311834d5   ! [s]  r = 1.0327 Rp
      real*8, parameter :: q_rad_slow   = 2.4415495d-6  ! [erg cm^-3 s^-1]
      real*8, parameter :: u_th_slow    = 4.9592176d-1  ! [erg cm^-3]
      real*8, parameter :: t_cross_stag = 5.4402086d5   ! [s]  r = 1.0323 Rp
      real*8, parameter :: q_rad_stag   = 2.4744558d-6  ! [erg cm^-3 s^-1]
      real*8, parameter :: u_th_stag    = 5.0081897d-1  ! [erg cm^-3]
      real*8 :: Da, Da_half, v_off

      integer :: nfail

      nfail = 0

      ! 1. The definition: one crossing time equal to one thermal time is
      !    Da = 1, the point at which the two terms of the equation are equal.
      Da = thermal_damkohler_number(u_th/q_rad, q_rad, u_th)
      call chk('damkohler_unity_at_equal_times', Da, 1.0d0, 1.0d-14, nfail)

      ! 2. Da is inversely proportional to the flow speed, so halving the speed
      !    doubles it: the correction leaves the advected regime smoothly and
      !    not on the sign of v.
      Da      = thermal_damkohler_number(dr_cm/(1.0d-2*cs), q_rad, u_th)
      Da_half = thermal_damkohler_number(dr_cm/(0.5d-2*cs), q_rad, u_th)
      call chk('damkohler_scales_as_inverse_speed', Da_half/Da, 2.0d0,        &
               1.0d-14, nfail)

      ! 3. and 4. The speed at which the condition turns over is finite and
      !    positive, v_off = dr q_rad / u_th, and the correction is off below
      !    it and on above it. Together with row 2 that is the v -> 0 limit of
      !    the energy correction: below v_off the cell keeps the temperature
      !    the run converged to, whatever v is.
      v_off = dr_cm*q_rad/u_th
      Da    = thermal_damkohler_number(dr_cm/(0.5d0*v_off), q_rad, u_th)
      call chk_gt('damkohler_off_below_turnover_speed', Da, 1.0d0, nfail)
      Da    = thermal_damkohler_number(dr_cm/(2.0d0*v_off), q_rad, u_th)
      call chk_lt('damkohler_on_above_turnover_speed', Da, 1.0d0, nfail)

      ! 5. A cell with no internal energy is the Da -> infinity end of the same
      !    statement: there is no temperature for the correction to start from.
      Da = thermal_damkohler_number(dr_cm/cs, q_rad, 0.0d0)
      call chk_gt('damkohler_off_without_internal_energy', Da, 1.0d0, nfail)

      ! 6. and 7. The base of a converged wind, on the two cells above: the
      !    correction stays on where the gas still crosses the cell faster
      !    than the radiative rate can rewrite its energy, and the next cell
      !    in, three times slower, is on the other side of the condition.
      Da = thermal_damkohler_number(t_cross_slow, q_rad_slow, u_th_slow)
      call chk_lt('damkohler_on_at_wind_base', Da, 1.0d0, nfail)
      Da = thermal_damkohler_number(t_cross_stag, q_rad_stag, u_th_stag)
      call chk_gt('damkohler_off_at_stagnating_cell', Da, 1.0d0, nfail)

      if (nfail .gt. 0) stop 1

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

      subroutine chk_gt(name, measured, bound, nf)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, bound
      integer, intent(inout) :: nf
      character(len=4) :: verdict
      verdict = 'PASS'
      if (.not. (measured .gt. bound)) then
         verdict = 'FAIL';  nf = nf + 1
      endif
      write(*,'(a,1x,a,a,es13.6,a,es13.6,a)') trim(verdict), trim(name),      &
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
      write(*,'(a,1x,a,a,es13.6,a,es13.6,a)') trim(verdict), trim(name),      &
         ' measured=', measured, ' reference=<', bound, ' tol=0'
      end subroutine chk_lt

      end program thermal_damkohler_limit
