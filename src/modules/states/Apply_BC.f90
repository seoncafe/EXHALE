   module BC_Apply
   ! Implementation of boundary conditions

   use global_parameters
   use Conversion
   use caloric_eos, only: adiabatic_index_from_state, caloric_cell_mixture
   use base_boundary, only: base_boundary_states, base_reservoir_p,       &
                            base_reservoir_T, base_reservoir_nhat,        &
                            r_base_level,                                 &
                            base_reservoir_prescription_version

   implicit none

   ! Number of outer ghost states whose free-outflow continuation left the
   ! admissible thermodynamic state (rho > 0, p > 0, both finite) and was
   ! dropped to a zero-gradient copy of the last physical cell, summed over
   ! the whole run. Reported at the end of a run. The continuation of
   ! free_outflow_ghost is positive whenever cell N is, so a nonzero count
   ! means cell N itself was not an admissible gas state, or the potential
   ! difference across the ghost overflowed the exponential; it is a
   ! statement about the interior, not about the boundary.
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
   !
   ! THE INVARIANT THAT MAKES THAT SAFE: they are a function of the interior
   ! state alone, so whoever replaces the interior state owes the next
   ! Reconstruct a call to Apply_BC. Rec_BC puts base_face_W straight into
   ! the left slot of the base face, so a stale value is a base-face flux
   ! built from a state the domain no longer holds. The one path that
   ! replaces the interior without a flux update of its own is the restore
   ! of the attempted-step checkpoint, and rebuild_state_from_checkpoint
   ! calls Apply_BC for exactly this reason.
   real*8 :: base_face_W(3)        = 0.0d0
   real*8 :: base_ghost_W(3,1-Ng:0) = 0.0d0
   real*8 :: base_face_lower_W(3)  = 0.0d0

   ! ---- WHICH STATE THE THREE ARRAYS ABOVE WERE BUILT FROM ----
   !
   ! The boundary state is a function of a stated input set, and these hold
   ! that input set as it stood at the derivation, so that a consumer can ask
   ! whether the cache belongs to the state it is about to use it on instead
   ! of assuming it (docs/lhs1140b_stationary_D5a_20260918.md sections 5 and
   ! 9, where the cache of one composition read by the residual of another is
   ! measured at 6.92 to 28.15 rounding floors of the base continuity row).
   !
   ! THE INPUT SET, and nothing else: the interior conserved state, the
   ! composition-derived state the ghost continuation and the caloric map
   ! read (the two ghost cells and the first interior cell of
   ! caloric_cell_mixture, and the cell-1 particle count), and the prescribed
   ! reservoir with its version. The ghost rows of a restart file are not in
   ! it: a restart rebuilds the boundary from the physical column.
   !
   ! bc_W_hold is the primitive array the derivation ran on, kept so that a
   ! read which finds the composition moved can recompute the face states
   ! with their own inputs rather than insert a state of another gas.
   logical :: base_boundary_installed = .false.
   real*8, allocatable :: bc_u_hold(:,:)   ! interior conserved, 1..N
   real*8, allocatable :: bc_W_hold(:,:)   ! the primitive array derived from
   real*8  :: bc_nk_hold(1-Ng:1)  = 0.0d0
   real*8  :: bc_xh2_hold(1-Ng:1) = 0.0d0
   logical :: bc_mol_hold(1-Ng:1) = .false.
   real*8  :: bc_npart1_hold      = 0.0d0
   real*8  :: bc_res_hold(4)      = 0.0d0
   integer :: bc_res_version_hold = -1

   ! Reads at which the cached face state did not belong to the composition
   ! installed at the read, and was recomputed from its own inputs before it
   ! was used. A nonzero count means a caller refreshed the composition and
   ! did not derive the boundary again; it is a statement about the call
   ! order, reported at the end of a run.
   integer :: n_base_face_state_recomputed = 0
   logical :: base_face_recompute_announced = .false.

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
   ! THE BOUNDARY-STATE OPERATION. One call derives everything the two
   ! boundaries state, from one stated input set, and tags the result with
   ! the state it was built from.
   !
   !   INPUT   the physical conserved state u(:,1:N) and the composition
   !           installed for it (through the composition-derived state
   !           caloric_cell_mixture and n_part_cell1 carry), the PRESCRIBED
   !           reservoir of base_boundary with its version, the model
   !           options. The ghost rows of a restart file are NOT an input:
   !           a restart rebuilds the boundary from the physical column
   !           (load_IC states the same contract at the read).
   !   OUTPUT  the ghost conserved state u(:,1-Ng:0) and u(:,N+1:N+Ng), and
   !           the cached face states base_face_W, base_ghost_W and
   !           base_face_lower_W, all of ONE composition, tagged with the
   !           input set they were built from so that a stale cache is
   !           detectable and not merely improbable.
   !
   ! The ghost COMPOSITION is not written here: it is solved, with the
   ! ghost's own ionization balance, by the composition sweep
   ! (ionization_equilibrium, the base handoff block), and this routine turns
   ! the composition installed for the ghost into its conserved state through
   ! the caloric map. The two together are the boundary state.
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

   ! The input set this boundary state belongs to.
   call hold_boundary_state_inputs(u, W)

   ! End of subroutine
   end subroutine Apply_BC

   !------------------------------------------!

   subroutine hold_boundary_state_inputs(u, W)
   ! Record the input set the boundary was just derived from: the interior
   ! conserved state, the primitive array the derivation ran on, the
   ! composition-derived state of the cells the boundary reads, and the
   ! prescribed reservoir with its version.
   real*8, intent(in) :: u(3,1-Ng:N+Ng), W(3,1-Ng:N+Ng)
   if (.not. allocated(bc_u_hold)) allocate(bc_u_hold(3,1:N))
   if (.not. allocated(bc_W_hold)) allocate(bc_W_hold(3,1-Ng:N+Ng))
   bc_u_hold = u(:,1:N)
   bc_W_hold = W
   call hold_boundary_closure_inputs()
   base_boundary_installed = .true.
   end subroutine hold_boundary_state_inputs

   !------------------------------------------!

   subroutine hold_boundary_closure_inputs()
   ! The half of the input set that is reachable without an argument: the
   ! composition-derived state of the two lower ghosts and the first interior
   ! cell, the cell-1 particle count, and the prescribed reservoir.
   integer :: j
   do j = 1-Ng, 1
      call caloric_cell_mixture(j, bc_nk_hold(j), bc_xh2_hold(j),         &
                                bc_mol_hold(j))
   enddo
   bc_npart1_hold      = n_part_cell1
   bc_res_hold(1)      = base_reservoir_p
   bc_res_hold(2)      = base_reservoir_T
   bc_res_hold(3)      = base_reservoir_nhat
   bc_res_hold(4)      = r_base_level
   bc_res_version_hold = base_reservoir_prescription_version
   end subroutine hold_boundary_closure_inputs

   !------------------------------------------!

   logical function base_boundary_closure_is_current() result(ok)
   ! Is the composition (and reservoir) the cached boundary was built from
   ! the one installed now? Compared value by value, so equality is the
   ! equality of the numbers and not of a digest of them.
   integer :: j
   real*8  :: nk, x2
   logical :: ismol
   ok = base_boundary_installed
   if (.not. ok) return
   do j = 1-Ng, 1
      call caloric_cell_mixture(j, nk, x2, ismol)
      ok = ok .and. (nk .eq. bc_nk_hold(j))                               &
              .and. (x2 .eq. bc_xh2_hold(j))                              &
              .and. (ismol .eqv. bc_mol_hold(j))
   enddo
   ok = ok .and. (n_part_cell1        .eq. bc_npart1_hold)                &
           .and. (base_reservoir_p    .eq. bc_res_hold(1))                &
           .and. (base_reservoir_T    .eq. bc_res_hold(2))                &
           .and. (base_reservoir_nhat .eq. bc_res_hold(3))                &
           .and. (r_base_level        .eq. bc_res_hold(4))                &
           .and. (base_reservoir_prescription_version                     &
                  .eq. bc_res_version_hold)
   end function base_boundary_closure_is_current

   !------------------------------------------!

   logical function base_boundary_cache_is_current(u) result(ok)
   ! The whole input set: the composition and reservoir above AND the
   ! interior conserved state. A consumer that holds the state can ask this;
   ! Rec_BC, which is handed reconstructed face states and not cell averages,
   ! can ask only the first half, and recomputes on it.
   real*8, intent(in) :: u(3,1-Ng:N+Ng)
   ok = base_boundary_closure_is_current()
   if (.not. ok) return
   if (.not. allocated(bc_u_hold)) then
      ok = .false.
      return
   endif
   ok = all(u(:,1:N) .eq. bc_u_hold)
   end function base_boundary_cache_is_current

   !------------------------------------------!

   subroutine report_base_face_cache_reads(tag)
   ! How many reads found the cache built from another composition and
   ! recomputed it. Zero is the statement that every flux assembly of this
   ! run stood on the boundary of the state it was assembling.
   character(len=*), intent(in) :: tag
   write(*,'(A,I0)') ' [base boundary] '//trim(tag)//': face states'//    &
        ' recomputed at a read because the composition had moved since'// &
        ' the boundary was derived: ', n_base_face_state_recomputed
   end subroutine report_base_face_cache_reads
   

   !------------------------------------------!

   subroutine Apply_BC_W(W)
   ! Boundary conditions for primitive variables (in place)

   real*8, intent(inout) :: W(3,1-Ng:N+Ng)
   integer :: k
   real*8  :: Wex(3)

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

   ! Upper boundary: free outflow, one rule for every reconstruction.
   ! free_outflow_ghost carries the whole statement and the reason for it.
   do k = 1,Ng
      call free_outflow_ghost(W,k,Wex)
      W(:,N+k) = Wex
   enddo

   ! End of subroutine
   end subroutine Apply_BC_W
   
   !------------------------------------------!

   subroutine free_outflow_ghost(W,k,Wg)
   ! One outer ghost CELL AVERAGE of the free-outflow boundary.
   !
   ! No characteristic enters the domain at an outflowing outer face, so the
   ! boundary imposes no condition and the ghost is a continuation of the
   ! interior. It is stated ONCE, for every reconstruction: which
   ! characteristics cross a face is a property of the flow and not of the
   ! scheme that reads the ghost.
   !
   ! THE CONTINUATION is the isothermal hydrostatic one,
   !
   !    rho_g = rho_N exp[-(phi_g - phi_N)/(p_N/rho_N)],  p_g = p_N x (same),
   !
   ! the state of rest the last physical cell carries outward under the code's
   ! own potential, with p/rho -- the isothermal sound speed squared of the
   ! mixture -- held at cell N's value. Since rho and p are scaled by the same
   ! factor, p_g/rho_g = p_N/rho_N: the ghost is at the SAME TEMPERATURE as
   ! cell N and differs from it only in density, so everything downstream that
   ! reads a ghost (the ionization sweep, the column integral of
   ! calc_column_dens, the sound speed of eval_dt) sees the temperature it saw
   ! before, at the density the stratification gives.
   !
   ! WHY A GRADIENT AND NOT A COPY. A zero-gradient ghost asserts dp/dr = 0 at
   ! the outer face, which is a hydrostatic imbalance of the order of the
   ! local weight rho g that no resolution removes. Under PLM it is worse than
   ! that: the forward difference entering the MC limiter of PLM_rec is then
   ! EXACTLY zero, a zero enters the limiter's argument list, the limited slope
   ! of cell N is exactly zero and that cell carries no pressure gradient at
   ! all. MEASURED on the analytic hydrostatic column by
   ! src/tests/grid_and_gates/hydrostatic_residual, the momentum residual of
   ! cell N in units of its own weight, at Grid cells = 250 / 500 / 1000 /
   ! 2000 and its order in 1/N:
   !
   !    zero-gradient copy (PLM)     0.759 0.753 0.751 0.750   order 0.006
   !    linear extrapolation         0.059 0.021 0.008 0.0030  order 1.42
   !    the continuation above       0.0044 0.0011 0.00032 0.00010  order 1.80
   !
   ! and with the continuation the largest residual of the domain is no longer
   ! at the boundary but at the base/stretched grid junction, where the
   ! interior truncation error is. The linear extrapolation taken in r instead
   ! of in the cell index was measured too and is no better than in the index
   ! (0.062 to 0.0032).
   !
   ! THE VELOCITY is a zero-gradient copy, v_g = v_N, and stays one. The
   ! outgoing characteristics carry the velocity out of the domain and the
   ! boundary has nothing of its own to say about it; a copy states the least
   ! and can neither reverse the flow (which would make the boundary an
   ! INFLOW, carrying characteristics inward that only a reservoir may state,
   ! and this closure has none) nor manufacture a speed. It leaves the limited
   ! slope of the velocity in cell N at zero under PLM, one variable at a time
   ! -- the MC limiter zeroes the slope of any variable whose forward
   ! difference vanishes -- which is a first-order velocity in one cell and
   ! not an unbalanced weight; the pressure is what the momentum residual
   ! measures and the pressure carries its gradient.
   !
   ! TWO GRADIENT-CARRYING VELOCITY RULES WERE MEASURED AND BOTH REJECTED, on
   ! wasp_full, whose outer edge is subsonic at Mach 0.96 so that the ghost
   ! reaches the interior:
   !   - a linear extrapolation in r carries about 1 percent less mass flux
   !     than cell N, and the run then sits in the base-breathing limit cycle
   !     with `du` between 4e-3 and 1.6e-2, never reaching the 1e-3 marching
   !     stop within 43000 steps, where the copy stops at 17135 and the
   !     zero-gradient ghost this replaced at 16208;
   !   - v_g = v_N (rho_N/rho_g)(r_N/r_g)^2, the steady continuity equation,
   !     stops the same case sooner (14442 steps) but divides by a ghost
   !     density that the stratification can make arbitrarily small: where the
   !     ghost cell is thicker than a scale height it manufactures an
   !     arbitrarily large speed at the boundary. It is also the only rule of
   !     the three under which wasp_full_newton reached a segmentation fault
   !     (step 4627 of the marching stage, reproducible in the -O3 build and
   !     absent at -O2 and under -fcheck, so a latent defect its trajectory
   !     reaches rather than one it contains).
   !
   ! DOMAIN. A supersonic outflow carries no information inward THROUGH THE
   ! RIEMANN SOLUTION: with every wave speed of the face positive, HLLC
   ! returns Phys_flux(WL) and never forms an expression in the right state
   ! (Num_Fluxes.f90 clamps SL = min(0, ...), so the branch is entered with
   ! SL exactly 0). It does still reach cell N through that cell's own
   ! reconstruction stencil, which is a property of the scheme's stencil and
   ! not of the characteristics, so a supersonic run moves under a change of
   ! this rule -- by the size of cell N's limited slope, and no further.
   ! MEASURED on mol_base_handoff at 12000 steps, whose outer edge sits at
   ! Mach 1.92: the outermost cell moves by 3.8e-4 in rho and it is the only
   ! cell of 500 that moves by more than 1e-4; the mass flux above r_flux is
   ! unchanged to 1.9e-8 and `du` to five digits.
   !
   ! The continuation therefore matters where the outer edge is SUBSONIC, and
   ! that is the state it is built for: there gravity is the leading term of
   ! dp/dr, the ram contribution rho v dv/dr being smaller than rho dphi/dr by
   ! the square of the Mach number.
   !
   ! FALLBACK. The exponential of a finite argument with p_N/rho_N > 0 is
   ! positive, so an admissible cell N gives an admissible ghost. A cell N
   ! that is not itself an admissible gas state, or a potential difference
   ! large enough to overflow the exponential, leaves the ghost with the
   ! zero-gradient copy of cell N and is counted in
   ! n_ghost_cells_positivity_limited. The tests are written as the negation
   ! of "admissible", so a NaN -- which compares false against everything --
   ! is caught too.
   real*8, intent(in)  :: W(3,1-Ng:N+Ng)
   integer, intent(in) :: k
   real*8, intent(out) :: Wg(3)
   real*8 :: csq, fac

   csq = W(3,N)/W(1,N)
   fac = 0.0d0
   if (csq .gt. 0.0d0) fac = exp(-(Gphi_c(N+k) - Gphi_c(N))/csq)
   Wg(1) = W(1,N)*fac
   Wg(3) = W(3,N)*fac

   if (.not. (Wg(1) .gt. 0.0d0 .and. Wg(3) .gt. 0.0d0 .and.              &
              Wg(1) .lt. huge(1.0d0) .and. Wg(3) .lt. huge(1.0d0))) then
      Wg = W(:,N)
      n_ghost_cells_positivity_limited =                                 &
         n_ghost_cells_positivity_limited + 1
      return
   endif

   Wg(2) = W(2,N)

   ! End of subroutine
   end subroutine free_outflow_ghost

   !------------------------------------------!
   
   subroutine Rec_BC(WL_in,WR_in,WL_out,WR_out)
   ! Boundary conditions for reconstructed variables

   real*8, dimension(3,1-Ng:N+Ng), intent(in) :: WL_in, WR_in
   integer :: k
   real*8, dimension(3,1-Ng:N+Ng), intent(out) :: WL_out, WR_out
   
   WL_out = WL_in
   WR_out = WR_in

   ! THE FACE STATE INSERTED BELOW IS THE ONE OF THE COMPOSITION INSTALLED
   ! NOW. The cache carries the input set it was derived from; where the
   ! composition-derived state has moved since, the three face states are
   ! derived again here, from the interior primitive array the boundary was
   ! derived from and the composition installed at this read, which is the
   ! same arithmetic Apply_BC performs. Nothing of another gas reaches the
   ! Riemann solve.
   !
   ! A recompute means the caller refreshed the composition and did not
   ! derive the boundary again, so the GHOST CELL AVERAGES it holds are still
   ! the old composition's; only the face states are repaired here, and the
   ! count says how often a call order left that repair to be done.
   if (.not. base_boundary_closure_is_current()) call recompute_base_face_states

   ! Lower boundary. The left state of the base face IS the boundary
   ! condition, and it is read here from the same call that wrote the ghost
   ! cell averages -- Apply_BC_W always runs on this state first, in the
   ! marching loop, in the residual, in init, and on a state restored from
   ! the attempted-step checkpoint. The right state of that face,
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

   subroutine recompute_base_face_states
   ! Derive the three cached face states again, from the interior primitive
   ! array the boundary was derived from and the composition installed now.
   ! Called from the read when the two no longer belong together; it does
   ! nothing when no boundary has been derived yet, since there is then no
   ! interior state to derive one from and the first Apply_BC of the run is
   ! still ahead.
   if (.not. base_boundary_installed) return
   if (.not. allocated(bc_W_hold)) return
   call base_boundary_states(bc_W_hold, base_face_W, base_ghost_W,        &
                             base_face_lower_W)
   call hold_boundary_closure_inputs()
   n_base_face_state_recomputed = n_base_face_state_recomputed + 1
   if (.not. base_face_recompute_announced) then
      base_face_recompute_announced = .true.
      write(*,'(A)') ' [base boundary] the cached base face state'//      &
           ' belonged to another composition at a read and was'//         &
           ' recomputed'
      write(*,'(A)') '   (a caller refreshed the composition without'//   &
           ' deriving the boundary again; counted for the run)'
   endif
   end subroutine recompute_base_face_states

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
