   module BC_Apply
   ! Implementation of boundary conditions

   use global_parameters
   use Conversion
   use caloric_eos, only: adiabatic_index_from_state
   use base_boundary, only: base_boundary_states

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
   ! fault to be repaired in place. The (p, s) reservoir of base_boundary
   ! states two conditions and this model has no third to state, so a
   ! supersonic base is a configuration to change, not a branch to add.
   integer :: n_base_supersonic_inflow_steps  = 0
   integer :: base_supersonic_inflow_first_step = -1
   real*8  :: base_inflow_mach_max = 0.0d0

   ! What the lower boundary states, as base_boundary last derived it: the
   ! primitive state AT THE FACE r_edg(0), the ghost CELL AVERAGES, and the
   ! point state at the next face down. Apply_BC_W fills them from the
   ! interior cell averages and Rec_BC reads them, so both paths use one
   ! answer from one function.
   real*8 :: base_face_W(3)        = 0.0d0
   real*8 :: base_ghost_W(3,1-Ng:0) = 0.0d0
   real*8 :: base_face_lower_W(3)  = 0.0d0

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
   cs2 = adiabatic_index_from_state(0,W(1,0),W(3,0))*W(3,0)/W(1,0)
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
   ! Boundary conditions for conservative variables (in place).
   !
   ! WRITES THE GHOSTS AND NOTHING ELSE. A boundary condition states the cells
   ! outside the domain; the interior is the caller's, and this routine hands
   ! it back unchanged, bit for bit.
   !
   ! It used to round-trip the interior as well -- U_to_W_interior, then
   ! W_to_U over the whole array -- and that round trip is not the identity in
   ! floating point: it recomputes the kinetic energy as 0.5*W(1)*W(2)**2
   ! where the input held 0.5*u(2)*u(2)/u(1), and those are not the same
   ! double. Two consequences, both real. (i) Apply_BC was not idempotent, so
   ! calling it twice moved every interior cell at the last bit; that is what
   ! blocked the repair of the ghost composition lag of section 141.6, since
   ! the repair is a second call. (ii) The interior of a state was changing
   ! under a routine whose subject is the boundary, which is not a property
   ! anything should have to reason about.
   !
   ! Writing only the ghosts removes both. It is not byte-identical to what it
   ! replaced -- it is the last bit of every interior cell of every call -- and
   ! the goldens are refreshed for it.
   real*8, intent(inout) :: u(3,1-Ng:N+Ng)
   real*8 :: W(3,1-Ng:N+Ng)
   integer :: k

   ! Interior cells only: the ghosts are outputs of Apply_BC_W below, never
   ! inputs, and on entry they can still be zero (first residual evaluation
   ! of the steady solver).
   call U_to_W_interior(u,W)

   ! Apply bc to W's
   call Apply_BC_W(W)

   ! Return to conservative variables -- the ghosts only. W_to_U_comp is the
   ! same arithmetic W_to_U performs cell by cell, and for a gas with no
   ! molecules in it energy_density_from_pressure returns
   ! p/(gamma_ad - 1) verbatim, so an atomic ghost is the double it always was.
   do k = 1,Ng
      call W_to_U_comp(W(:,1-k), u(:,1-k), 1-k)
      call W_to_U_comp(W(:,N+k), u(:,N+k), N+k)
   enddo

   ! End of subroutine
   end subroutine Apply_BC
   

   !------------------------------------------!

   subroutine Apply_BC_W(W)
   ! Boundary conditions for primitive variables (in place)

   real*8, intent(inout) :: W(3,1-Ng:N+Ng)
   integer :: k
   real*8  :: Wzg(3), Wex(3)

   ! Lower boundary: the characteristic face condition. base_boundary_states
   ! derives the face state, the ghost cell averages and the state at the next
   ! face down from the first interior CELL AVERAGE, all at once, and stores
   ! them for Rec_BC. The ghosts written here are the volume averages of the
   ! hydrostatic isentrope through the face state, NOT copies of it: the two
   ! differ by half a cell of stratification, and writing the face value into
   ! a cell average is the defect this replaces.
   call base_boundary_states(W, base_face_W, base_ghost_W, base_face_lower_W)
   do k = 1,Ng
      W(:,1-k) = base_ghost_W(:,1-k)
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
   ! While the PLM -> WENO3 continuation is running, the two free-outflow
   ! ghosts -- the zero-gradient copy the PLM stage uses and the linear
   ! extrapolation the WENO3 stencil asks for -- are combined under the same
   ! lambda the right-hand side is combined under, so that the boundary
   ! operator travels the same homotopy the interior does. The positivity
   ! guard is applied to the extrapolated state BEFORE the blend, so a
   ! repaired extrapolation is what enters the combination and the counter
   ! means what it meant. At lambda = 1 this is the WENO3 ghost to the bit;
   ! at lambda = 0 the branch is not taken at all.
   do k = 1,Ng
      if (recon_lambda_on .and. recon_lambda .gt. 0.0d0) then
         Wzg = W(:,N)
         Wex = 2.0*W(:,N+k-1) - W(:,N+k-2)
         if (.not. (Wex(1) .gt. 0.0d0 .and. Wex(3) .gt. 0.0d0)) then
            Wex = W(:,N+k-1)
            n_ghost_cells_positivity_limited =                            &
               n_ghost_cells_positivity_limited + 1
         endif
         W(:,N+k) = (1.0d0 - recon_lambda)*Wzg + recon_lambda*Wex
      else
         W(:,N+k) = W(:,N)
         if (use_weno3) then
            W(:,N+k) = 2.0*W(:,N+k-1) - W(:,N+k-2)
            if (.not. (W(1,N+k) .gt. 0.0d0 .and. W(3,N+k) .gt. 0.0d0)) then
               W(:,N+k) = W(:,N+k-1)
               n_ghost_cells_positivity_limited =                         &
                  n_ghost_cells_positivity_limited + 1
            endif
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
   
   ! Lower boundary. The left state of the base face IS the boundary
   ! condition, and it is read here from the same call that wrote the ghost
   ! cell averages -- Apply_BC_W always runs on this state first, in the
   ! marching loop, in the residual and in init. The right state of that face,
   ! WR_out(:,0), is the interior's reconstruction of cell 1 down to it and is
   ! deliberately left alone: the Riemann solver is what mediates the two.
   !
   ! The face below the ghost carries the same continuation on both sides, so
   ! it is jump free; its flux lies outside the cells the update writes.
   !
   ! With Ng = 2 the lower boundary has exactly two faces to state, r_edg(0)
   ! and r_edg(-1), and base_boundary returns one state for each. A larger Ng
   ! would add faces between them, and base_ghost_averages would have to return
   ! the continuation at each; the loop below would then be over them rather
   ! than over the single lower face.
   WL_out(:,0)    = base_face_W
   do k = 2,Ng
      WL_out(:,1-k) = base_face_lower_W
      WR_out(:,1-k) = base_face_lower_W
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
