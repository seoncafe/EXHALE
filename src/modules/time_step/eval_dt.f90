   module eval_time_step

   use global_parameters
   use grid_construction, only: spherical_cell_volume
   use caloric_eos, only: caloric_mixture_active,                     &
                         adiabatic_index_from_state

   implicit none

   contains

   subroutine eval_dt(W,dt,dt_loc)
   ! Explicit stability interval for evolved cells 1..N. dt is their
   ! minimum; dt_loc uses each cell interval under "Time stepping: Local"
   ! and otherwise equals dt. The stationary fixed point is dt-independent.
   !
   ! RK_rhs uses spherical finite-volume divergence:
   !   dF_j=(A_+ F_+ - A_- F_-)/V_j,
   !   A=r^2, V_j=(r_+^3-r_-^3)/3.
   ! SSP-RK3 has forward-Euler SSP coefficient 1. Face Riemann signals
   ! satisfy S_f=max(|v|+c) over the bounding states, giving
   !   dt_j=CFL*V_j/max(A_- S_-,A_+ S_+), CFL<=1.
   ! In the Cartesian uniform limit this is CFL*dr/(|v|+c), with input
   ! "CFL:" default 0.6. Replacing max by the sum prevents opposing fans
   ! from meeting and corresponds to CFL/2 on a uniform grid.
   !
   ! Spherical geometry has V_j/A_+ < dr_j < V_j/A_-; with equal face
   ! speeds the outer face gives the tighter limit. In general the face
   ! with larger A_f S_f binds. Ghost cells supply boundary Riemann states
   ! at r_edg(0) and r_edg(N). S_f uses cell averages rather than
   ! reconstructed states; the latter can have larger |v|+c, and the
   ! selected CFL below one supplies the assumed margin.
   real*8, dimension(3,1-Ng:N+Ng), intent(in) :: W
   real*8, dimension(1-Ng:N+Ng) :: rho,v,p,cs
   real*8, dimension(0:N) :: s_face
   real*8 :: dt_j
   integer :: j
   real*8, intent(out) :: dt
   real*8, dimension(1-Ng:N+Ng), intent(out) :: dt_loc

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

   ! The fastest signal of each face of an evolved cell, the two boundary
   ! faces included.
   do j = 0, N
      s_face(j) = max(abs(v(j))   + cs(j),                              &
                      abs(v(j+1)) + cs(j+1))
   enddo

   ! The interval of each evolved cell, and the minimum over them.
   dt = huge(1.0d0)
   do j = 1, N
      dt_j = CFL*spherical_cell_volume(j)                               &
             /max(r_edg(j)*r_edg(j)*s_face(j),                          &
                  r_edg(j-1)*r_edg(j-1)*s_face(j-1))
      dt_loc(j) = dt_j
      dt = min(dt, dt_j)
   enddo

   if (use_local_dt) then
      ! A ghost row is not evolved: it is rebuilt by Apply_BC after every
      ! stage. It is given the interval of the evolved cell it borders so
      ! that a stage formed over the padded range is arithmetic on a
      ! physical interval and not on an undefined one.
      dt_loc(1-Ng:0)   = dt_loc(1)
      dt_loc(N+1:N+Ng) = dt_loc(N)
   else
      dt_loc = dt
   endif

   end subroutine eval_dt

   end module eval_time_step
