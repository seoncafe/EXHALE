   module BC_Apply
   ! Implementation of boundary conditions

   use global_parameters
   use Conversion

   implicit none

   ! Number of ghost-cell states whose linear extrapolation left the
   ! admissible thermodynamic state (rho > 0, p > 0) and was dropped to a
   ! zero-gradient copy of the last admissible cell, summed over the whole
   ! run. Reported at the end of a run; zero means every extrapolated ghost
   ! was admissible and the run is the one an unguarded build would have
   ! produced.
   integer :: n_ghost_cells_positivity_limited = 0

   ! Mach number of the gas the lower boundary admits, and whether it ever
   ! became supersonic.
   !
   ! This is a validity check on the boundary condition itself, independent of
   ! the chemistry. A subsonic inflow carries one outgoing characteristic, so
   ! one piece of interior information belongs in the ghost and the other two
   ! are the reservoir's to state. A SUPERSONIC inflow carries none: all three
   ! characteristics point into the domain, every ghost variable is the
   ! reservoir's, and a boundary condition that copies an interior velocity
   ! into it is over-specified -- it is feeding the domain a state the domain
   ! itself produced. The solution downstream of such a boundary is not a
   ! solution of the stated problem, whatever the chemistry does.
   !
   ! Reported, not enforced: the run continues and says so, because the
   ! crossing is evidence about the configuration rather than a numerical
   ! fault to be repaired in place. `Base velocity: massflux` removes the
   ! copied-velocity feedback and is the usual response, but it does not make
   ! this check redundant -- a mass-flux base can still go supersonic if the
   ! reservoir state and the flux constant disagree.
   integer :: n_base_supersonic_inflow_steps  = 0
   integer :: base_supersonic_inflow_first_step = -1
   real*8  :: base_inflow_mach_max = 0.0d0

   contains

   !------------------------------------------!

   real*8 function base_inflow_mach_number(W)
   ! |v|/c of the base ghost the lower boundary condition has just written,
   ! with the code's own adiabatic sound speed c = sqrt(gamma p/rho) -- the
   ! same expression eval_dt uses for the CFL condition, so the two agree on
   ! what "sonic" means.
   !
   ! Measured on ghost 0, the cell the first interior face sees. Returns zero
   ! for a state with no sound speed (a non-positive or NaN rho or p), so a
   ! broken ghost is reported by the positivity counter above rather than
   ! twice.
   real*8, intent(in) :: W(3,1-Ng:N+Ng)
   real*8 :: cs2
   cs2 = g*W(3,0)/W(1,0)
   if (.not. (cs2 .gt. 0.0d0)) then
      base_inflow_mach_number = 0.0d0
   else
      base_inflow_mach_number = abs(W(2,0))/sqrt(cs2)
   endif
   end function base_inflow_mach_number

   !------------------------------------------!

   subroutine check_base_inflow_is_subsonic(W,step)
   ! Record the base Mach number of this step and warn on the FIRST step at
   ! which the boundary admits supersonic inflow. Inflow is v > 0 here: the
   ! radial coordinate increases outward, so gas entering the domain at the
   ! base moves outward. An outflowing (v < 0) base is a different condition
   ! and is not what this check is about, so it is measured but not warned on.
   real*8, intent(in)  :: W(3,1-Ng:N+Ng)
   integer, intent(in) :: step
   real*8 :: mach
   mach = base_inflow_mach_number(W)
   if (mach .gt. base_inflow_mach_max) base_inflow_mach_max = mach
   if (W(2,0) .gt. 0.0d0 .and. mach .gt. 1.0d0) then
      n_base_supersonic_inflow_steps = n_base_supersonic_inflow_steps + 1
      if (base_supersonic_inflow_first_step .lt. 0) then
         base_supersonic_inflow_first_step = step
         write(*,'(A,I0,A,ES10.3)')                                       &
            '     base boundary WARNING: inflow became supersonic at step ',&
            step, ', Mach = ', mach
         write(*,*) '       All three characteristics now enter the'//    &
                    ' domain, so the ghost is over-'
         write(*,*) '       specified and the interior is not a solution'//&
                    ' of the stated problem.'
      endif
   endif
   end subroutine check_base_inflow_is_subsonic

   subroutine Apply_BC(u)
   ! Boundary conditions for conservative variables (in place)

   real*8, intent(inout) :: u(3,1-Ng:N+Ng)
   real*8 :: W(3,1-Ng:N+Ng)

   ! Interior cells only: the ghosts are outputs of Apply_BC_W below, never
   ! inputs, and on entry they can still be zero (first residual evaluation
   ! of the steady solver).
   call U_to_W_interior(u,W)

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
      ! A pressure gradient steep enough to extrapolate through zero puts the
      ! ghost outside the admissible state, where the gradient continuity the
      ! extrapolation buys is worthless: the ghost is read as a real gas state
      ! by the reconstruction stencil and by eval_dt. Fall back on the
      ! zero-gradient base pressure, which is the first-order limit of the
      ! same boundary condition. Density stays pinned at rho_bc above.
      if (.not. (W_in(3,index) .gt. 0.0d0)) then
         W_in(3,index) = W_in(3,1)
         n_ghost_cells_positivity_limited =                               &
            n_ghost_cells_positivity_limited + 1
      endif
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

   ! Upper boundary: free outflow. The zero-gradient copy is the boundary
   ! condition itself; the linear extrapolation the WENO3 stencil asks for is
   ! an accuracy device on top of it, valid only where the state it
   ! extrapolates is smooth.
   !
   ! Where it is not -- a cold-start transient, a front crossing the outer
   ! cells -- the extrapolation can carry rho or p through zero, and the
   ! result is not a poor approximation but no gas state at all. These ghosts
   ! are read as real cells further downstream: calc_column_dens starts its
   ! column integral at N+Ng, the ionization sweep runs over 1-Ng..N+Ng, and
   ! eval_dt takes a minimum over the same range, so a negative ghost density
   ! becomes negative nuclei densities and negative stage populations inside
   ! the photoionization heating. Falling back on the zero-gradient state --
   ! the first-order limit of the same free-outflow condition -- keeps the
   ! boundary physical.
   !
   ! Written as the negation of "strictly positive", so a NaN ghost -- which
   ! compares false against everything -- is caught too.
   do k = 1,Ng
      W(:,N+k) = W(:,N)
      if (use_weno3) then
         W(:,N+k) = 2.0*W(:,N+k-1) - W(:,N+k-2)
         if (.not. (W(1,N+k) .gt. 0.0d0 .and. W(3,N+k) .gt. 0.0d0)) then
            W(:,N+k) = W(:,N+k-1)
            n_ghost_cells_positivity_limited =                            &
               n_ghost_cells_positivity_limited + 1
         endif
      endif
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
         
   ! Upper boundary. The extrapolations below are left unguarded: these are
   ! face states, and positivity_limited_faces already drops any face that
   ! leaves rho > 0, p > 0 back to the cell averages.
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
