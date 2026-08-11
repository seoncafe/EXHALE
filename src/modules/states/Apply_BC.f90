   module BC_Apply
   ! Implementation of boundary conditions

   use global_parameters
   use Conversion

   implicit none
   
   contains
   
   subroutine Apply_BC(u)
   ! Boundary conditions for conservative variables (in place)

   real*8, intent(inout) :: u(3,1-Ng:N+Ng)
   real*8 :: W(3,1-Ng:N+Ng)

   ! Convert to primitive variables
   call U_to_W(u,W)

   ! Apply bc to W's
   call Apply_BC_W(W)

   ! Return to conservaive variables
   call W_to_U(W,u)

   ! End of subroutine
   end subroutine Apply_BC
   
   !------------------------------------------!
   
   subroutine BC_component_constrho(W_in,index)
   ! Component-wise BC in case of zero-velocity gradient
   !     at the lower boundary

   real*8, intent(inout)  :: W_in(3,1-Ng:N+Ng)
   integer, intent(in)    :: index

   ! Density: always the base anchor rho_bc (the mass reservoir).
   W_in(1,index) = rho_bc
   ! Velocity. base_v_massflux: CETIMB-style v0 = F_c/(rho_bc r^2) with F_c the
   ! mass-flux constant from the [j_min:N] constant-momentum region (NOT the
   ! base, where rho*v*r^2 is not yet flat). Else the legacy one-way valve.
   if (base_v_massflux .and. base_flux_const .gt. 0.0d0) then
      W_in(2,index) = base_flux_const/(rho_bc*r(index)**2)
   else if (valve_eps .gt. 0.0d0) then
      ! Smooth one-way valve 0.5*(v + sqrt(v^2 + eps^2)): differentiable at
      ! v=0, -> 0 as v -> -inf, -> v for v >> eps (bias +eps/2 only near
      ! v ~ 0). Needed by the steady-state Newton solver, whose line search
      ! cannot cross the hard-max kink the breathing base sits on.
      W_in(2,index) = 0.5d0*(W_in(2,1)                                   &
                      + sqrt(W_in(2,1)**2 + valve_eps**2))
   else
      W_in(2,index) = max(W_in(2,1),0.0)
   endif
   ! Pressure. Legacy: fixed ntot_bc + dp_bc (-> T = T0 isothermal base).
   ! hydrostatic_base (momentum-consistent): extrapolate the interior pressure
   ! gradient (cells 1,2) into the ghost so dp/dr is CONTINUOUS at the base
   ! rather than flattened to 0; the base-face pressure gradient can then
   ! balance gravity (the source of the breathing momentum residual). With
   ! rho pinned to rho_bc, T_base floats slightly off T0.
   ! base_ghost_T_continuous: dT/dr = 0 instead of T = T0. The ghost keeps the
   ! base composition (ntot_bc nuclei + dp_bc electrons at rho_bc, the same
   ! particle count the isothermal pin uses) but carries the temperature of
   ! the first interior cell, T(1) = W(3,1)/n_part_cell1 in units of T0:
   !    p_ghost = (ntot_bc + dp_bc)*T(1).
   ! n_part_cell1 = n_tot(1) + n_e(1) comes from the composition solve, so the
   ! ghost pressure stays a differentiable function of the interior pressure
   ! (what the JFNK line search needs) while the ionization state it divides by
   ! is lagged exactly like every other composition quantity in a hydro step.
   if (hydrostatic_base) then
      W_in(3,index) = W_in(3,1) + (W_in(3,2) - W_in(3,1))                 &
                      /(r(2) - r(1))*(r(index) - r(1))
   else if (base_ghost_T_continuous) then
      W_in(3,index) = (ntot_bc + dp_bc)*W_in(3,1)/n_part_cell1
   else
      W_in(3,index) = ntot_bc + dp_bc
   endif

   ! End of subroutine
   end subroutine BC_component_constrho

   !------------------------------------------!

   subroutine Apply_BC_W(W)
   ! Boundary conditions for primitive variables (in place)

   real*8, intent(inout) :: W(3,1-Ng:N+Ng)
   integer :: k

   ! BC with constant rho at lower boundary
   do k = 1,Ng
      call BC_component_constrho(W,1-k)
   enddo

   do k = 1,Ng
      ! Upper boundary
      W(:,N+k) = W(:,N)
      if (use_weno3) W(:,N+k) = 2.0*W(:,N+k-1) - W(:,N+k-2)
   enddo

   ! End of subroutine
   end subroutine Apply_BC_W
   
   !------------------------------------------!
   
   subroutine Rec_BC(WL_in,WR_in,WL_out,WR_out)
   ! Boundary conditions for reconstructed variables

   real*8, dimension(3,1-Ng:N+Ng), intent(in) :: WL_in, WR_in
   integer :: k
   real*8, dimension(3,1-Ng:N+Ng), intent(out) :: WL_out, WR_out
   
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
      WR_out(:,N+k) = WR_out(:,N)
      if (use_weno3) WR_out(:,N+k) = 2.0*WR_out(:,N+k-1) - WR_out(:,N+k-2)
   enddo
   
   WL_out(:,N+2) = WL_out(:,N+1)
   if (use_weno3) WL_out(:,N+2) = 2.0*WL_out(:,N+1) - WL_out(:,N)

   ! End of subroutine
   end subroutine Rec_BC

   !------------------------------------------!

   subroutine shapiro_filter(u)
   ! Periodic 1-2-1 Shapiro (1970) low-pass filter on the conservative
   ! variables, to damp the gravity-unbalanced sound waves (base breathing)
   ! the way CETIMB (Koskinen et al. 2013a) does. No-op if shapiro_eps <= 0.
   ! Caller must Apply_BC before (ghosts j=0, N+1 are used) and after (to reset
   ! the ghosts to their BC values).
   real*8, intent(inout) :: u(3,1-Ng:N+Ng)
   real*8 :: f(3,1-Ng:N+Ng)
   integer :: j, k
   if (shapiro_eps .le. 0.0d0) return
   f = u
   do k = 1,3
      do j = 1,N
         u(k,j) = f(k,j) + 0.25d0*shapiro_eps                            &
                  *(f(k,j-1) - 2.0d0*f(k,j) + f(k,j+1))
      enddo
   enddo
   end subroutine shapiro_filter

   ! End of module
   end module BC_Apply
