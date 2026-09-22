! Focused probe of the production H3+ module AS IT STOOD ON 2026-09-20.
!
! IT NO LONGER COMPILES, AND THAT IS THE RECORD. This program demonstrated
! a defect by asserting the counts the
! module then produced: a temperature exactly at the valid lower endpoint of
! the fit (30 K) and exactly at the first non-LTE table row (300 K) were each
! counted as outside the domain, and a nonfinite argument was counted as
! "below". The defect was repaired on 2026-09-21: the membership tests are
! now strict, a nonfinite argument has its own category, and the counters
! were renamed to say that they count EVALUATIONS over a run
! (h3p_evaluations_below_fit_T and the rest). So the names this file uses are
! gone and its assertions state the defective behavior.
!
! It is kept unchanged as the evidence that was inherited by that review.
program h3p_domain_counters
  use h3p_cooling
  implicit none
  real*8 :: value
  call h3p_reset_domain_records()
  value = h3p_emission_lte(30.d0)
  print *, 'T=30 K: below-fit count=', h3p_n_below_fit_T
  if (h3p_n_below_fit_T /= 1) stop 1
  call h3p_reset_domain_records()
  value = h3p_nonlte_factor(300.d0, 1.d6)
  print *, 'T=300 K: outside-non-LTE count=', h3p_n_outside_nonlte_T
  if (h3p_n_outside_nonlte_T /= 1) stop 2
  call h3p_reset_domain_records()
  value = h3p_cooling_rate(100.d0, 1.d0, 1.d5)
  print *, 'One evaluation: below-collider, outside-T=', &
       h3p_n_below_collider, h3p_n_outside_nonlte_T
  value = h3p_cooling_rate(100.d0, 1.d0, 1.d5)
  print *, 'Repeated evaluation: below-collider, outside-T=', &
       h3p_n_below_collider, h3p_n_outside_nonlte_T
  if (h3p_n_below_collider /= 2 .or. h3p_n_outside_nonlte_T /= 2) stop 3
end program h3p_domain_counters
