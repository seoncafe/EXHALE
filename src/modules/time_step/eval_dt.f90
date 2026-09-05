   module eval_time_step
    
   use global_parameters
   use caloric_eos, only: caloric_mixture_active,                     &
                         adiabatic_index_from_state

   implicit none

   contains
        
   subroutine eval_dt(W,dt,dt_loc)
   ! Evaluate the time step according to the CFL condition.
   ! dt     = global CFL step (minimum over the grid; original behavior).
   ! dt_loc = cell-by-cell pseudo-time step used by the RK and source updates.
   !          With "Time stepping: Local" each cell gets its own CFL step
   !          dt_j = CFL*dr_j/(|v_j|+cs_j) (steady-state acceleration: the
   !          fixed point dF = S, heating = cooling is dt-independent);
   !          otherwise dt_loc is uniformly the global dt, which makes the
   !          update arithmetic bit-identical to the original scalar form.
   real*8, dimension(3,1-Ng:N+Ng), intent(in) :: W
   real*8, dimension(1-Ng:N+Ng) :: rho,v,p,cs
   integer :: j
   real*8, intent(out) :: dt
   real*8, dimension(1-Ng:N+Ng), intent(out) :: dt_loc

   ! Extract physical variables
   rho = W(1,:)
   v   = W(2,:)
   p   = W(3,:)

   ! Evaluate sound speed
   if (caloric_mixture_active) then
      do j = 1-Ng, N+Ng
         cs(j) = sqrt(adiabatic_index_from_state(j,rho(j),p(j))         &
                      *p(j)/rho(j))
      enddo
   else
   cs = sqrt(gamma_ad*p/rho)
   endif

   ! Evaluate time step according to CFL condition
   dt = CFL*minval(dr_j/(abs(v) + cs))

   if (use_local_dt) then
      dt_loc = CFL*(dr_j/(abs(v) + cs))
   else
      dt_loc = dt
   endif

   ! End of subroutine
   end subroutine eval_dt

   ! End of module
   end module eval_time_step