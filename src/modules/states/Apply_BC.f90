   module BC_Apply
   ! Implementation of boundary conditions

   use global_parameters
   use Conversion

   implicit none
   
   contains
   
   subroutine Apply_BC(u_in,u_out)
   ! Boundary conditions for conservative variables

   real*8, intent(in) :: u_in(1-Ng:N+Ng,3)
   real*8 :: W(1-Ng:N+Ng,3)
   real*8, intent(out) :: u_out(1-Ng:N+Ng,3)
   
   ! Convert to primitive variables
   call U_to_W(u_in,W)
   
   ! Apply bc to W's
   call Apply_BC_W(W,W)     
   
   ! Return to conservaive variables
   call W_to_U(W,u_out)     
   
   ! End of subroutine
   end subroutine Apply_BC
   
   !------------------------------------------!
   
   subroutine BC_component_constrho(W_in,index)
   ! Component-wise BC in case of zero-velocity gradient
   !     at the lower boundary

   real*8, intent(inout)  :: W_in(1-Ng:N+Ng,3)
   integer, intent(in)    :: index

   ! Density: always the base anchor rho_bc (the mass reservoir).
   W_in(index,1) = rho_bc
   ! Velocity. base_v_massflux: CETIMB-style v0 = F_c/(rho_bc r^2) with F_c the
   ! mass-flux constant from the [j_min:N] constant-momentum region (NOT the
   ! base, where rho*v*r^2 is not yet flat). Else the legacy one-way valve.
   if (base_v_massflux .and. base_flux_const .gt. 0.0d0) then
      W_in(index,2) = base_flux_const/(rho_bc*r(index)**2)
   else if (valve_eps .gt. 0.0d0) then
      ! Smooth one-way valve 0.5*(v + sqrt(v^2 + eps^2)): differentiable at
      ! v=0, -> 0 as v -> -inf, -> v for v >> eps (bias +eps/2 only near
      ! v ~ 0). Needed by the steady-state Newton solver, whose line search
      ! cannot cross the hard-max kink the breathing base sits on.
      W_in(index,2) = 0.5d0*(W_in(1,2)                                   &
                      + sqrt(W_in(1,2)**2 + valve_eps**2))
   else
      W_in(index,2) = max(W_in(1,2),0.0)
   endif
   ! Pressure. Legacy: fixed ntot_bc + dp_bc (-> T = T0 isothermal base).
   ! hydrostatic_base (momentum-consistent): extrapolate the interior pressure
   ! gradient (cells 1,2) into the ghost so dp/dr is CONTINUOUS at the base
   ! rather than flattened to 0; the base-face pressure gradient can then
   ! balance gravity (the source of the breathing momentum residual). With
   ! rho pinned to rho_bc, T_base floats slightly off T0.
   if (hydrostatic_base) then
      W_in(index,3) = W_in(1,3) + (W_in(2,3) - W_in(1,3))                 &
                      /(r(2) - r(1))*(r(index) - r(1))
   else
      W_in(index,3) = ntot_bc + dp_bc
   endif

   ! End of subroutine
   end subroutine BC_component_constrho

   !------------------------------------------!

   subroutine Apply_BC_W(W_in,W_out)
   ! Boundary conditions for primitive variables
   
   real*8, intent(in)  :: W_in(1-Ng:N+Ng,3)
   integer :: k
   real*8, intent(out) :: W_out(1-Ng:N+Ng,3)
   
   ! Copy input vector
   W_out = W_in
   
   ! BC with constant rho at lower boundary
   do k = 1,Ng
      call BC_component_constrho(W_out,1-k)
   enddo
         
   do k = 1,Ng
      ! Upper boundary
      W_out(N+k,:) = W_out(N,:)
      if (use_weno3) W_out(N+k,:) = 2.0*W_out(N+k-1,:) - W_out(N+k-2,:)
   enddo
   
   ! End of subroutine
   end subroutine Apply_BC_W
   
   !------------------------------------------!
   
   subroutine Rec_BC(WL_in,WR_in,WL_out,WR_out)
   ! Boundary conditions for reconstructed variables

   real*8, dimension(1-Ng:N+Ng,3), intent(in) :: WL_in, WR_in
   integer :: k
   real*8, dimension(1-Ng:N+Ng,3), intent(out) :: WL_out, WR_out
   
   WL_out = WL_in
   WR_out = WR_in
   
   ! Lower boundary
   call BC_component_constrho(WR_out,1-Ng)
   
   ! BC with constant density at lower boundary
   do k = 1,Ng
         call BC_component_constrho(WL_out,1-k)
   enddo
         
   ! Upper boundary
   do k = 1,Ng
      WR_out(N+k,:) = WR_out(N,:)
      if (use_weno3) WR_out(N+k,:) = 2.0*WR_out(N+k-1,:) - WR_out(N+k-2,:)
   enddo
   
   WL_out(N+2,:) = WL_out(N+1,:)
   if (use_weno3) WL_out(N+2,:) = 2.0*WL_out(N+1,:) - WL_out(N,:)

   ! End of subroutine
   end subroutine Rec_BC

   !------------------------------------------!

   subroutine shapiro_filter(u)
   ! Periodic 1-2-1 Shapiro (1970) low-pass filter on the conservative
   ! variables, to damp the gravity-unbalanced sound waves (base breathing)
   ! the way CETIMB (Koskinen et al. 2013a) does. No-op if shapiro_eps <= 0.
   ! Caller must Apply_BC before (ghosts j=0, N+1 are used) and after (to reset
   ! the ghosts to their BC values).
   real*8, intent(inout) :: u(1-Ng:N+Ng,3)
   real*8 :: f(1-Ng:N+Ng,3)
   integer :: j, k
   if (shapiro_eps .le. 0.0d0) return
   f = u
   do k = 1,3
      do j = 1,N
         u(j,k) = f(j,k) + 0.25d0*shapiro_eps                            &
                  *(f(j-1,k) - 2.0d0*f(j,k) + f(j+1,k))
      enddo
   enddo
   end subroutine shapiro_filter

   !------------------------------------------!

   subroutine viscous_accel(vel, Tcell, Fv)
   ! Leading CETIMB viscous momentum acceleration (Koskinen 2022 B5, first term):
   ! F_mu = (4/3)(1/r^2) d/dr(r^2 mu dvel/dr), with mu = visc_mu0 * Tcell^visc_s
   ! in code units. experimental (un-validated): no-op if visc_mu0 <= 0;
   ! the (dmu/dr)(dvel/dr) and -(16/3)mu vel/r^2 corrections, the viscous
   ! dissipation q_mu + heat conduction, and a stable semi-implicit time
   ! integration are deferred to a later revision (calibration + validation vs Koskinen).
   real*8, intent(in)  :: vel(1-Ng:N+Ng), Tcell(1-Ng:N+Ng)
   real*8, intent(out) :: Fv(1-Ng:N+Ng)
   real*8 :: mu(1-Ng:N+Ng), rp, rm, mup, mum, fluxp, fluxm
   integer :: j
   Fv = 0.0d0
   if (visc_mu0 .le. 0.0d0) return
   do j = 1-Ng, N+Ng
      mu(j) = visc_mu0*Tcell(j)**visc_s
   enddo
   do j = 1, N
      rp = 0.5d0*(r(j) + r(j+1));  mup = 0.5d0*(mu(j) + mu(j+1))
      rm = 0.5d0*(r(j) + r(j-1));  mum = 0.5d0*(mu(j) + mu(j-1))
      fluxp = rp*rp*mup*(vel(j+1) - vel(j))/(r(j+1) - r(j))
      fluxm = rm*rm*mum*(vel(j) - vel(j-1))/(r(j) - r(j-1))
      Fv(j) = (4.0d0/3.0d0)/(r(j)*r(j))                                  &
              *(fluxp - fluxm)/(0.5d0*(r(j+1) - r(j-1)))
   enddo
   end subroutine viscous_accel

   ! End of module
   end module BC_Apply
