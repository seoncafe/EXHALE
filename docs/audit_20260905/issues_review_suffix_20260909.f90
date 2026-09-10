end module numerical_review

program issues_review
  use numerical_review
  use, intrinsic :: ieee_arithmetic
  implicit none
  real*8 :: su(1),sn(1),au,an,delta,r,g,ag,tau_old,tau_line
  real*8 :: y(1),f0(1),fsp(1,1),d(1),dr(1),abf(1,1),b(1),x(1),ax(1),rr
  real*8 :: ratio,heh_dev
  integer :: piv(1),iterations
  logical :: boundary,truncated

  su=-0.5d0; sn=-1d0; delta=0d0
  call dogleg_step(su,sn,delta,au,an,boundary)
  write(*,'(A,ES15.7)') 'zero_radius_step_norm=',abs(au*su(1)+an*sn(1))
  if (abs(au*su(1)+an*sn(1)) /= 0d0) error stop 1

  ! True scaled Jacobian A=1; approximate band transpose gives g=10.
  r=1d0; g=10d0; ag=10d0
  tau_old=g*g/(ag*ag)
  tau_line=r*ag/(ag*ag)
  write(*,'(A,ES15.7)') 'initial_merit=',0.5d0*r*r
  write(*,'(A,ES15.7)') 'current_cauchy_merit=',0.5d0*(r-tau_old*ag)**2
  write(*,'(A,ES15.7)') 'line_minimizing_merit=',0.5d0*(r-tau_line*ag)**2
  if (0.5d0*(r-tau_old*ag)**2 <= 0.5d0*r*r) error stop 2

  y=0d0; f0=0d0; fsp=0d0; d=1d0; dr=1d0; abf=1d0; b=1d0; piv=1
  operator_value=1d0
  call pgmres(y,f0,fsp,d,dr,abf,piv,0d0,b,x,1,0.1d0,iterations, &
              truncated,ax,rr)
  write(*,'(A,ES15.7)') 'identity_gmres_solution=',x(1)
  if (abs(x(1)-1d0)>1d-14) error stop 3

  operator_value=0d0
  call pgmres(y,f0,fsp,d,dr,abf,piv,0d0,b,x,1,0.1d0,iterations, &
              truncated,ax,rr)
  write(*,'(A,L1)') 'zero_operator_gmres_solution_finite=',ieee_is_finite(x(1))
  write(*,'(A,ES15.7)') 'zero_operator_reported_relative_residual=',rr
  if (ieee_is_finite(x(1))) error stop 4

  ! The exact production enthalpy diagnostic is inserted by the driver.
  ratio=enthalpy_flux_term_ratio(2d0,0.25d0,1d0,1d0,1d0,1d0,0d0,1d0,1d0)
  write(*,'(A,ES15.7)') 'nonzero_mass_divergence_term_ratio=',ratio
  if (ratio>=1d0) error stop 5
  heh_dev=abs(7.9307d-2-0.0793d0)/0.0793d0
  write(*,'(A,ES15.7)') 'quoted_restart_values_relative_difference=',heh_dev
  if (heh_dev<=1d-6) error stop 6
  ! Independent carrier ceilings do not enforce a shared element budget.
  write(*,'(A,F6.3)') 'hydrogen_inventory_at_individual_ceilings=',2d0*0.5d0+1d0
  if (2d0*0.5d0+1d0<=1d0) error stop 7
  ! A shell outside the stellar projected radius still crosses an inner ray.
  write(*,'(A,ES15.7)') 'outer_shell_chord_coordinate=',sqrt(12d0**2-2d0**2)
  if (.not. (12d0>10d0 .and. 2d0<10d0 .and. 12d0>=2d0)) error stop 8
  write(*,'(A)') 'All review assertions passed; defects reproduced, not repaired.'
end program issues_review
