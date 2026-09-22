      program h3p_cooling_limits
      ! The limits, the continuities and the published values of the H3+
      ! infrared cooling rate.
      !
      ! Origin: the H3+ lines of docs/audit_20260905/audit_probe.f90 and of
      ! roe_equal_pressure.f90 (which also prints the non-LTE factor at the
      ! last tabulated density).  Production routines exercised:
      ! h3p_cooling_rate, h3p_emission_lte and h3p_nonlte_factor of
      ! src/modules/lower_atmosphere/h3p_cooling.f90.
      !
      ! REFERENCES.
      !  * Collisional limit.  H3+ radiates from levels populated by
      !    collisions with H2.  Far below the density at which the
      !    populations thermalize, every collisional excitation is followed
      !    by a radiative decay, so the emission per H3+ ion is set by the
      !    collision rate: Lambda is proportional to n_H2 and vanishes at
      !    n_H2 = 0.  Miller et al. (2013) Table 6 places the transition to
      !    LTE near n_H2 = 1e10 cm^-3 at 1000 K (their factor is 0.4955
      !    there and 0.0013 at their lowest tabulated density, 1e6 cm^-3);
      !    the densities used below, 1 and 100 cm^-3, are four to six
      !    decades below that lowest tabulated density and far inside the
      !    collisional limit.  The factor 2 tolerance on the two-decade
      !    ratio is what an extrapolation that is only asymptotically
      !    linear is allowed to differ by.
      !  * Continuity of the LTE emission.  The reference is the JOIN RULE
      !    of h3p_cooling.f90, NOT an identity of Miller et al. (2013)
      !    Table 5.  Their four segments do not join: continued across its
      !    upper limit, the 300-800 K polynomial stands 43.05 per cent above
      !    the 30-300 K one at 300 K, the 800-1800 K one 0.113 per cent
      !    below the 300-800 K one at 800 K, and the 1800-5000 K one
      !    2.353 per cent below the 800-1800 K one at 1800 K (measured with
      !    the Table 5 coefficients; the paper's own Table 6 LTE value at
      !    300 K, 0.53503e-23 W molecule^-1 sr^-1, is the one the upper
      !    segment reproduces).  A cooling function with a 43 per cent step
      !    in it is not a function of state, so the code blends the two
      !    published polynomials in log_e E over the 5 per cent of the join
      !    temperature below each join, and the assertion below is that this
      !    stated rule holds: the evaluated emission is continuous.  The
      !    tolerance 1e-3 is far looser than the ramp needs and is kept from
      !    the four figures of the published coefficients.
      !  * Continuity of the non-LTE factor is a property of the physical
      !    function: s(T,n_H2) reaches unity from below, so it may not step
      !    at the last tabulated column.
      !  * Published values.  Table 4 gives E(H3+,T) computed with z35000(T)
      !    at 500 to 5000 K, and Table 6 gives the LTE emission at 300 K;
      !    the tolerances below are the fit errors the paper quotes for each
      !    range: less than 0.1 per cent over 300-800 K, "all but
      !    nonexistent" over 800-1800 K, a few tenths of a per cent over
      !    1800-5000 K.  The paper tabulates no E(H3+,T) below 300 K, so the
      !    30-300 K segment carries a transcription assertion instead: the
      !    Table 5 coefficients are re-typed here and evaluated at 150 K.
      !  * Domain records.  b1_target_system_20260906.md section 6.1
      !    decision 8: a cell below the tabulated collider range is
      !    evaluated on the linear collisional limit and carries an
      !    informational record; a cell outside 30-5000 K is clamped and
      !    recorded.  The records count EVALUATIONS, and they count only
      !    STRICT excursions: the Table 5 fit is defined at 30 and at
      !    5000 K and the Table 6 rows run from 300 to 5000 K with columns
      !    from 1e6 to 1e14 cm^-3, so an argument sitting exactly on one of
      !    those endpoints is inside the published domain.  An argument
      !    that is not an ordinary real is no statement about an interval
      !    at all and carries its own record.
      !
      ! At 35d9dd5 the density argument was clamped at the lowest tabulated
      ! density (log n_H2 = 6), so the rate was constant below it and
      ! nonzero at n_H2 = 0; the factor stepped from 0.9985 to 1 at
      ! n_H2 = 1e14 cm^-3 at 5000 K; the four published segments were used
      ! as published and their steps stood in the evaluated emission.
      use h3p_cooling, only: h3p_cooling_rate, h3p_emission_lte,          &
                             h3p_nonlte_factor, h3p_reset_domain_records, &
                             h3p_evaluations_below_collider,              &
                             h3p_evaluations_below_fit_T,                 &
                             h3p_evaluations_above_fit_T,                 &
                             h3p_evaluations_outside_nonlte_T,            &
                             h3p_evaluations_nonfinite_T,                 &
                             h3p_evaluations_nonfinite_collider
      use assertion_report
      implicit none

      real*8, parameter :: t_collisional = 1000.0d0
      real*8, parameter :: joins(3) = (/ 300.0d0, 800.0d0, 1800.0d0 /)
      real*8, parameter :: n_lte_edge = 1.0d14
      real*8, parameter :: t_edge = 5000.0d0

      ! Miller et al. (2013) Table 4, column "E(H3+,T) using (z35000(T))",
      ! in W molecule^-1 sr^-1, and the range each temperature falls in.
      integer, parameter :: n_pub = 6
      real*8, parameter :: t_pub(n_pub) =                                 &
           (/ 500.d0, 1000.d0, 1500.d0, 2000.d0, 3000.d0, 5000.d0 /)
      real*8, parameter :: e_pub(n_pub) = (/ 0.67786d-21, 0.33393d-19,    &
           0.19205d-18, 0.60456d-18, 0.20089d-17, 0.34119d-17 /)
      real*8, parameter :: tol_pub(n_pub) = (/ 1.0d-3, 1.0d-3, 1.0d-3,    &
           5.0d-3, 5.0d-3, 5.0d-3 /)
      ! Table 6, [H2] = 1e20 m^-3 row-entry at 300 K: the LTE emission.
      real*8, parameter :: e_pub_300 = 0.53503d-23

      ! Table 5, 30-300 K column, re-typed from the paper.
      real*8, parameter :: cA_paper(0:9) = (/ -81.9599d0, +0.886768d0,    &
           -0.0264611d0, +0.000462693d0, -4.70108d-6, +2.84979d-8,        &
           -1.03090d-10, +2.13794d-13, -2.26029d-16, +8.66357d-20 /)

      ! The IEEE-754 binary64 patterns of a quiet NaN and of +infinity,
      ! written out so that neither has to be produced by an arithmetic
      ! operation the compiler may fold or trap on.
      integer*8, parameter :: nan_bits = int(z'7FF8000000000000', 8)
      integer*8, parameter :: inf_bits = int(z'7FF0000000000000', 8)

      real*8 :: lam0, lam1, lam2, jump, tk, ssum, dummy
      real*8 :: below, at, qnan, pinf
      character(len=48) :: label
      integer :: i, n

      lam0 = h3p_cooling_rate(t_collisional, 1.0d0, 0.0d0)
      lam1 = h3p_cooling_rate(t_collisional, 1.0d0, 1.0d0)
      lam2 = h3p_cooling_rate(t_collisional, 1.0d0, 1.0d2)

      ! No collision partner, no collisionally excited emission.
      call check_absolute('h3p_cooling_vanishes_without_h2',              &
                          lam0/lam2, 0.0d0, 1.0d-12)

      ! Linear in the collision-partner density over two decades.
      call check_within_factor('h3p_cooling_collisional_slope',           &
                               lam2/lam1, 1.0d2, 2.0d0)

      ! Continuity of the evaluated LTE emission across the three joins.
      do i = 1, 3
         tk = joins(i)
         jump = h3p_emission_lte(tk + 1.0d-7)/h3p_emission_lte(tk-1.0d-7) &
                - 1.0d0
         write(label,'(a,i0)') 'h3p_lte_emission_join_', int(tk)
         call check_absolute(label, jump, 0.0d0, 1.0d-3)
      enddo

      ! Continuity of the non-LTE factor at the last tabulated density.
      below = h3p_nonlte_factor(t_edge, n_lte_edge*(1.0d0 - 1.0d-8))
      at    = h3p_nonlte_factor(t_edge, n_lte_edge)
      call check_absolute('h3p_nonlte_factor_high_density_join',          &
                          at/below - 1.0d0, 0.0d0, 1.0d-3)

      ! One published value inside each fit segment.
      do i = 1, n_pub
         write(label,'(a,i0)') 'h3p_lte_emission_table4_', int(t_pub(i))
         call check_relative(label, h3p_emission_lte(t_pub(i)),           &
                             e_pub(i), tol_pub(i))
      enddo
      call check_relative('h3p_lte_emission_table6_lte_300',              &
                          h3p_emission_lte(300.0d0), e_pub_300, 1.0d-3)

      ! 30-300 K segment: coefficient transcription, the paper tabulating
      ! no E(H3+,T) there.
      ssum = 0.0d0
      do n = 9, 0, -1
         ssum = ssum*150.0d0 + cA_paper(n)
      enddo
      call check_relative('h3p_lte_emission_table5_coefficients_150',     &
                          h3p_emission_lte(150.0d0), exp(ssum), 1.0d-12)

      ! Domain records: informational, and raised where they should be.
      call h3p_reset_domain_records()
      dummy = h3p_cooling_rate(1000.0d0, 1.0d0, 1.0d3)
      call check_absolute('h3p_record_below_collider',                    &
                          dble(h3p_evaluations_below_collider),           &
                          1.0d0, 0.5d0)
      call h3p_reset_domain_records()
      dummy = h3p_emission_lte(20.0d0)
      dummy = h3p_emission_lte(6000.0d0)
      call check_absolute('h3p_record_outside_fit_temperature',           &
                          dble(h3p_evaluations_below_fit_T                &
                               + h3p_evaluations_above_fit_T),            &
                          2.0d0, 0.5d0)
      call h3p_reset_domain_records()
      dummy = h3p_nonlte_factor(6000.0d0, 1.0d12)
      call check_absolute('h3p_record_outside_nonlte_temperature',        &
                          dble(h3p_evaluations_outside_nonlte_T),         &
                          1.0d0, 0.5d0)

      ! THE ENDPOINTS ARE INSIDE, AND ONE UNIT IN THE LAST PLACE OUTSIDE
      ! THEM IS OUT.  nearest() moves by the smallest representable step,
      ! so each pair brackets the endpoint as tightly as the arithmetic
      ! allows and no tolerance on the temperature is involved.
      call h3p_reset_domain_records()
      dummy = h3p_emission_lte(30.0d0)
      call check_absolute('h3p_record_fit_temperature_30_is_inside',      &
                          dble(h3p_evaluations_below_fit_T), 0.0d0, 0.5d0)
      call h3p_reset_domain_records()
      dummy = h3p_emission_lte(nearest(30.0d0, -1.0d0))
      call check_absolute('h3p_record_fit_temperature_below_30',          &
                          dble(h3p_evaluations_below_fit_T), 1.0d0, 0.5d0)
      call h3p_reset_domain_records()
      dummy = h3p_emission_lte(5000.0d0)
      call check_absolute('h3p_record_fit_temperature_5000_is_inside',    &
                          dble(h3p_evaluations_above_fit_T), 0.0d0, 0.5d0)
      call h3p_reset_domain_records()
      dummy = h3p_emission_lte(nearest(5000.0d0, 1.0d0))
      call check_absolute('h3p_record_fit_temperature_above_5000',        &
                          dble(h3p_evaluations_above_fit_T), 1.0d0, 0.5d0)
      call h3p_reset_domain_records()
      dummy = h3p_nonlte_factor(300.0d0, 1.0d12)
      dummy = h3p_nonlte_factor(5000.0d0, 1.0d12)
      call check_absolute('h3p_record_nonlte_row_endpoints_are_inside',   &
                          dble(h3p_evaluations_outside_nonlte_T),         &
                          0.0d0, 0.5d0)
      call h3p_reset_domain_records()
      dummy = h3p_nonlte_factor(nearest(300.0d0, -1.0d0), 1.0d12)
      call check_absolute('h3p_record_nonlte_row_below_300',              &
                          dble(h3p_evaluations_outside_nonlte_T),         &
                          1.0d0, 0.5d0)
      call h3p_reset_domain_records()
      dummy = h3p_nonlte_factor(1000.0d0, 1.0d6)
      call check_absolute('h3p_record_collider_1e6_is_inside',            &
                          dble(h3p_evaluations_below_collider),           &
                          0.0d0, 0.5d0)
      call h3p_reset_domain_records()
      dummy = h3p_nonlte_factor(1000.0d0, nearest(1.0d6, -1.0d0))
      call check_absolute('h3p_record_collider_below_1e6',                &
                          dble(h3p_evaluations_below_collider),           &
                          1.0d0, 0.5d0)

      ! AN ARGUMENT THAT IS NOT AN ORDINARY REAL IS ITS OWN CATEGORY.  The
      ! evaluation still lands on an end of the fit and on a zero collider
      ! density, which is what the clamps are for, but calling that an
      ! excursion below the published range would be a statement about an
      ! interval that a NaN or an infinity never made.
      qnan = transfer(nan_bits, 1.0d0)
      pinf = transfer(inf_bits, 1.0d0)
      call h3p_reset_domain_records()
      dummy = h3p_emission_lte(qnan)
      call check_absolute('h3p_record_nonfinite_temperature',             &
                          dble(h3p_evaluations_nonfinite_T), 1.0d0, 0.5d0)
      call check_absolute('h3p_record_nonfinite_temperature_not_below',   &
                          dble(h3p_evaluations_below_fit_T), 0.0d0, 0.5d0)
      call h3p_reset_domain_records()
      dummy = h3p_emission_lte(pinf)
      call check_absolute('h3p_record_infinite_temperature',              &
                          dble(h3p_evaluations_nonfinite_T), 1.0d0, 0.5d0)
      call check_absolute('h3p_record_infinite_temperature_not_above',    &
                          dble(h3p_evaluations_above_fit_T), 0.0d0, 0.5d0)
      call h3p_reset_domain_records()
      dummy = h3p_nonlte_factor(1000.0d0, qnan)
      call check_absolute('h3p_record_nonfinite_collider',                &
                          dble(h3p_evaluations_nonfinite_collider),       &
                          1.0d0, 0.5d0)
      call check_absolute('h3p_record_nonfinite_collider_not_below',      &
                          dble(h3p_evaluations_below_collider),           &
                          0.0d0, 0.5d0)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'h3p_cooling_limits: ', assertion_failures,  &
              ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'h3p_cooling_limits: all assertions passed'

      end program h3p_cooling_limits
