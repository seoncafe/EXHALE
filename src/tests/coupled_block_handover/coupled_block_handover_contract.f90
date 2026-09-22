      program coupled_block_handover_contract
      ! THE TRANSITION FROM "THE ALTERNATION STOPPED AT THE MOVEMENT BOUND"
      ! TO "THE COUPLED CARRIER BLOCK IS ENTERED", AND THE ONE STATE IT
      ! CARRIES.
      !
      ! WHAT WAS WRONG. The outer iteration of
      ! steady_wind_with_element_diffusion decides the progress of a pass
      ! AHEAD of that pass's composition update, so the pass that declares
      ! the stall runs no relaxation at all: its carrier outcome token is
      ! the reset value, "nothing to advance". The handover asked for that
      ! token to read "the composition movement bound" on the same pass on
      ! which the iteration had already declared no progress, and a pass
      ! ending in no progress returns before any later pass. The two
      ! conjuncts cannot hold together, so the branch was unreachable.
      !
      ! WHAT IS ASSERTED HERE, on a sequence of passes built to be the one
      ! the defect is about (three relaxations ending on the bound, then a
      ! pass that declares no progress and runs no relaxation):
      !
      !   0. THE DEFECT ITSELF, kept as a row: on the stalling pass the
      !      LIVE outcome token is not the movement bound, so a condition
      !      that reads it there can never be met.
      !   a. that sequence reaches the transition and the block is entered;
      !   b. the transition is taken exactly once: a second attempt is
      !      refused and the count of entries stays at one;
      !   c. the residual the coupled block starts from is a fresh residual
      !      OF THE ENTRY SNAPSHOT: the restored state is the recorded
      !      state bit for bit, so the same operator returns the recorded
      !      number, and a re-evaluated number that has moved, or one
      !      evaluated at a generation the snapshot has left, is refused
      !      rather than accepted;
      !   d. with "Coupled carrier solve: On stall" not set nothing fires:
      !      the transition stays at alternation_running with the option
      !      named as the reason, and the refusal the alternation already
      !      had is what stands;
      !   e. a pass sequence whose relaxations did NOT end on the movement
      !      bound triggers nothing, however flat the scalar norm: a joint
      !      distance that stopped falling is one condition of three.
      !
      ! WHAT STANDS IN FOR THE CARRIER BALANCE. The transition contract is
      ! about WHICH state a number belongs to, not about the operator that
      ! forms it, so the balance below is a deterministic reduction of the
      ! state arrays. Its only property the assertions use is that it is a
      ! function of those arrays and of nothing else, which is the property
      ! the production carrier balance has as well.
      !
      ! Each assertion prints
      !     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
      ! and the exit status is nonzero if any of them fails.
      use coupled_block_handover
      implicit none

      integer, parameter :: ncell = 8
      integer, parameter :: nsp   = 4
      ! The number of consecutive bound endings the transition asks for,
      ! the same number the outer iteration asks of its progress rule
      ! (outer_no_fall_max in EXHALE_main.f90).
      integer, parameter :: bound_endings_required = 3

      ! The outcome tokens of the carrier relaxation, in the values
      ! diffusive_photochemistry gives them; only the movement bound and
      ! the reset value are needed here.
      integer, parameter :: relax_nothing_to_advance = 0
      integer, parameter :: relax_movement_bound     = 3

      type(coupled_entry_snapshot)     :: snap_stall, snap_off, snap_flat
      type(alternation_stall_evidence) :: ev
      real*8  :: u(3,ncell), f_sp(ncell,nsp)
      real*8  :: rho(ncell), v(ncell), p(ncell), T(ncell), Frho(ncell)
      real*8  :: u_work(3,ncell), f_work(ncell,nsp)
      real*8  :: rho_w(ncell), v_w(ncell), p_w(ncell), T_w(ncell)
      real*8  :: Frho_w(ncell)
      real*8  :: res_recorded, res_fresh, res_restored
      real*8  :: cgap, egap
      integer :: state_now, why, n_fail, gen_at_evaluation
      integer :: live_token_on_the_stalling_pass
      logical :: written, restored, fits

      n_fail = 0
      u_work = 0.0d0;  f_work = 0.0d0
      rho_w  = 0.0d0;  v_w    = 0.0d0;  p_w = 0.0d0;  T_w = 0.0d0
      Frho_w = 0.0d0

      ! ---- 0. THE DEFECT, AS A ROW ----
      ! The pass that declares the stall takes no composition update, so
      ! its outcome token is the reset value however the passes before it
      ! ended.
      live_token_on_the_stalling_pass = relax_nothing_to_advance
      call check_int('live_token_on_the_stalling_pass_is_not_the_bound',  &
           live_token_on_the_stalling_pass, relax_nothing_to_advance)
      call check_log('old_conjunction_cannot_hold',                       &
           (live_token_on_the_stalling_pass .eq. relax_movement_bound),   &
           .false.)

      ! ---- d. THE CONTROL: THE OPTION IS NOT SET ----
      ! The same pass sequence, with "Coupled carrier solve: On stall"
      ! off. Nothing is recorded by the outer iteration in that case; the
      ! snapshot is written here anyway, so that the refusal below is the
      ! option's and not the absence of a record.
      call three_bound_passes_then_a_stall(snap_off)
      call coupled_entry_stall_transition(snap_off, .false., .false.,     &
           .true., .true., bound_endings_required, state_now, why)
      call check_int('option_off_stays_at_alternation_running',           &
           state_now, alternation_running)
      call check_int('option_off_names_the_option',                       &
           why, entry_refused_option_off)
      call check_int('option_off_enters_no_block',                        &
           n_coupled_block_entries(), 0)

      ! ---- e. NO MOVEMENT-BOUND ENDING ----
      ! The joint distance stopped falling and the relaxations ended
      ! elsewhere. A flat scalar norm on its own is not the stall this
      ! transition is about.
      call three_passes_off_the_bound_then_a_stall(snap_flat)
      call coupled_entry_stall_transition(snap_flat, .true., .false.,     &
           .true., .true., bound_endings_required, state_now, why)
      call check_int('a_flat_norm_alone_stays_at_alternation_running',    &
           state_now, alternation_running)
      call check_int('a_flat_norm_alone_names_the_missing_bound_ending',  &
           why, entry_refused_not_bound_ending)
      call check_int('a_flat_norm_alone_enters_no_block',                 &
           n_coupled_block_entries(), 0)
      ! And the same record with the joint distance still falling: the
      ! stall condition itself is read and not assumed.
      call coupled_entry_stall_transition(snap_flat, .true., .false.,     &
           .true., .false., bound_endings_required, state_now, why)
      call check_int('a_falling_joint_distance_names_no_stall',           &
           why, entry_refused_no_stall)

      ! ---- a. THE CASE THAT REACHES BOTH CONDITIONS ----
      call three_bound_passes_then_a_stall(snap_stall)
      res_recorded = snap_stall%evidence%carrier_residual
      call coupled_entry_stall_transition(snap_stall, .true., .false.,    &
           .true., .true., bound_endings_required, state_now, why)
      call check_int('the_stalled_sequence_reaches_the_transition',       &
           state_now, alternation_stalled_at_bound)
      call check_int('the_stalled_sequence_is_refused_by_nothing',        &
           why, entry_refused_none)
      ! Too few consecutive bound endings is its own refusal.
      call coupled_entry_stall_transition(snap_stall, .true., .false.,    &
           .true., .true., bound_endings_required + 1, state_now, why)
      call check_int('one_bound_ending_short_is_refused',                 &
           why, entry_refused_few_bound_ends)

      ! ---- c. THE FIRST RESIDUAL OF THE BLOCK IS A FRESH RESIDUAL OF THE
      !         ENTRY SNAPSHOT ----
      ! Evaluated on the snapshot's OWN arrays, at the snapshot's own
      ! generation.
      gen_at_evaluation = snap_stall%generation
      res_fresh = carrier_balance(snap_stall%rho, snap_stall%v,           &
                                  snap_stall%f_sp)
      call check_real('the_fresh_entry_residual_is_the_recorded_one',     &
           res_fresh, res_recorded, 0.0d0)

      ! A number evaluated at a generation the snapshot has left cannot
      ! authorize it: the snapshot is rewritten by one more pass and the
      ! stale generation is refused.
      call one_more_bound_pass(snap_stall)
      fits = coupled_entry_state_fits(snap_stall, u_work, f_work, rho_w,  &
                                      v_w, p_w, T_w, Frho_w)
      call check_log('the_recorded_state_fits_the_arrays_of_this_run',    &
           fits, .true.)
      call accept_coupled_block_entry(snap_stall, gen_at_evaluation,      &
           fits, res_fresh, .true., 0.0d0, .false., state_now, why,       &
           cgap, egap)
      call check_int('a_residual_of_a_generation_that_moved_is_refused',  &
           why, entry_refused_snapshot_moved)
      call check_int('a_stale_generation_enters_no_block',                &
           n_coupled_block_entries(), 0)

      ! The same acceptance at the CURRENT generation, with a residual
      ! that has moved beyond the stated bound: refused on the number.
      gen_at_evaluation = snap_stall%generation
      res_recorded = snap_stall%evidence%carrier_residual
      call accept_coupled_block_entry(snap_stall, gen_at_evaluation,      &
           fits, res_recorded*1.01d0, .true., 0.0d0, .false.,             &
           state_now, why, cgap, egap)
      call check_int('a_residual_that_moved_is_refused',                  &
           why, entry_refused_residual_moved)
      call check_real('the_refusal_reports_its_own_distance',             &
           cgap, 1.0d-2/1.01d0, 1.0d-12)
      call check_int('a_moved_residual_enters_no_block',                  &
           n_coupled_block_entries(), 0)

      ! And the acceptance the sequence is entitled to.
      res_fresh = carrier_balance(snap_stall%rho, snap_stall%v,           &
                                  snap_stall%f_sp)
      call accept_coupled_block_entry(snap_stall, gen_at_evaluation,      &
           fits, res_fresh, .true., 0.0d0, .false., state_now, why,       &
           cgap, egap)
      call check_int('the_block_is_entered', state_now,                   &
           coupled_block_entry)
      call check_int('the_entry_is_refused_by_nothing', why,              &
           entry_refused_none)

      ! The state the block consumes is the recorded one, whole: restored
      ! into the working arrays it returns the same number the recorded
      ! pass measured, bit for bit, because it is the same operator at the
      ! same bits.
      u_work = 0.0d0;  f_work = 0.0d0
      rho_w  = 0.0d0;  v_w    = 0.0d0;  p_w = 0.0d0;  T_w = 0.0d0
      Frho_w = 0.0d0
      call restore_coupled_entry_state(snap_stall, u_work, f_work, rho_w, &
           v_w, p_w, T_w, Frho_w, restored)
      call check_log('the_entry_state_is_handed_back', restored, .true.)
      res_restored = carrier_balance(rho_w, v_w, f_work)
      call check_real('the_first_residual_of_the_block_is_the_entry_one', &
           res_restored, snap_stall%evidence%carrier_residual, 0.0d0)
      call check_real('the_restored_conserved_state_is_the_recorded_one', &
           maxval(abs(u_work - snap_stall%u)), 0.0d0, 0.0d0)
      call check_real('the_restored_composition_is_the_recorded_one',     &
           maxval(abs(f_work - snap_stall%f_sp)), 0.0d0, 0.0d0)
      call check_real('the_restored_face_mass_flux_is_the_recorded_one',  &
           maxval(abs(Frho_w - snap_stall%Frho_elem)), 0.0d0, 0.0d0)

      ! ---- b. EXACTLY ONCE ----
      call check_int('the_block_was_entered_once',                        &
           n_coupled_block_entries(), 1)
      ! A second transition, with the caller's own token set as the outer
      ! iteration sets it.
      call coupled_entry_stall_transition(snap_stall, .true., .true.,     &
           .true., .true., bound_endings_required, state_now, why)
      call check_int('a_second_transition_stays_at_running',              &
           state_now, alternation_running)
      call check_int('a_second_transition_names_the_entry_already_taken', &
           why, entry_refused_already_entered)
      ! And with the caller's token NOT set: the seal alone refuses it.
      gen_at_evaluation = snap_stall%generation
      call accept_coupled_block_entry(snap_stall, gen_at_evaluation,      &
           fits, res_fresh, .true., 0.0d0, .false., state_now, why,       &
           cgap, egap)
      call check_int('a_second_acceptance_is_refused_by_the_seal',        &
           why, entry_refused_already_entered)
      call check_int('the_block_was_still_entered_once',                  &
           n_coupled_block_entries(), 1)
      ! The sealed snapshot is immutable: a later pass cannot rewrite it.
      call one_more_bound_pass_written(snap_stall, written)
      call check_log('a_sealed_snapshot_is_not_rewritten', written,       &
           .false.)
      call check_real('a_sealed_snapshot_keeps_its_residual',             &
           snap_stall%evidence%carrier_residual, res_fresh, 0.0d0)

      if (n_fail .gt. 0) then
         write(*,'(A,I0,A)') 'coupled_block_handover_contract: ', n_fail, &
              ' assertion(s) failed'
         stop 1
      endif
      write(*,'(A)') 'coupled_block_handover_contract: every assertion'// &
           ' passed'

      contains

      ! --------------------------------------------------------------- !

      real*8 function carrier_balance(rho_in, v_in, f_in) result(b)
      ! A DETERMINISTIC REDUCTION OF THE STATE, standing for the carrier
      ! balance of the production code. What the assertions use is that it
      ! is a function of these arrays and of nothing else.
      real*8, intent(in) :: rho_in(:), v_in(:), f_in(:,:)
      integer :: j
      b = 0.0d0
      do j = 1, size(rho_in)
         b = max(b, abs(rho_in(j)*f_in(j,1) - v_in(j))/                   &
                    (1.0d0 + abs(rho_in(j))))
      enddo
      end function carrier_balance

      ! --------------------------------------------------------------- !

      subroutine state_of_pass(ipass)
      ! One state of the sequence, distinct from pass to pass so that a
      ! number of the wrong pass cannot pass for a number of this one.
      integer, intent(in) :: ipass
      integer :: j, is
      do j = 1, ncell
         rho(j)  = 1.0d0 + 0.1d0*j + 0.01d0*ipass
         v(j)    = 0.5d0*j - 0.05d0*ipass
         p(j)    = 2.0d0 + 0.2d0*j
         T(j)    = 1.0d3 + 10.0d0*j + ipass
         Frho(j) = 3.0d0 + 0.3d0*j
         u(1,j)  = rho(j)
         u(2,j)  = rho(j)*v(j)
         u(3,j)  = p(j)/0.6666666666666666d0 + 0.5d0*rho(j)*v(j)**2
         do is = 1, nsp
            f_sp(j,is) = 0.25d0 + 0.01d0*is - 0.001d0*j + 0.002d0*ipass
         enddo
      enddo
      end subroutine state_of_pass

      ! --------------------------------------------------------------- !

      subroutine write_the_pass(snap, ipass, on_the_bound, n_bound,       &
                                n_no_fall, written_out)
      ! The outer iteration's own write: one state, one set of evidence,
      ! one call.
      type(coupled_entry_snapshot), intent(inout) :: snap
      integer, intent(in)  :: ipass, n_bound, n_no_fall
      logical, intent(in)  :: on_the_bound
      logical, intent(out) :: written_out
      call state_of_pass(ipass)
      ev%pass = ipass
      if (on_the_bound) then
         ev%carrier_outcome      = relax_movement_bound
         ev%carrier_outcome_text = 'the composition movement bound'
         ev%ended_on_movement_bound = .true.
      else
         ev%carrier_outcome      = relax_nothing_to_advance
         ev%carrier_outcome_text = 'the fixed point of the carrier block'
         ev%ended_on_movement_bound = .false.
      endif
      ev%movement_bound    = 1.0d-2
      ev%bound_cell        = 3
      ev%bound_carrier     = 1
      ev%bound_measure     = 1.0d-2
      ev%bound_entry       = 0.5d0
      ev%bound_is_fraction = .true.
      ev%carrier_residual_available = .true.
      ev%carrier_residual  = carrier_balance(rho, v, f_sp)
      ev%carrier_res_abs   = 1.0d-3
      ev%carrier_scale_abs = 1.0d0
      ev%carrier_worst_cell    = 3
      ev%carrier_worst_carrier = 1
      ev%element_residual_available = .false.
      ev%progress_reference = 4.2d0
      ev%n_no_fall          = n_no_fall
      ev%n_bound_endings    = n_bound
      ev%boundary_rebuild_suppressed = .false.
      ev%reconstruction_operator = 'PLM'
      ev%reconstruction_is_plm   = .true.
      call record_coupled_entry_state(snap, ev, u, f_sp, rho, v, p, T,    &
                                      Frho, written_out)
      end subroutine write_the_pass

      ! --------------------------------------------------------------- !

      subroutine three_bound_passes_then_a_stall(snap)
      ! Passes 1 to 3 relax and end on the movement bound; pass 4 declares
      ! no progress and runs no relaxation, so it writes nothing and the
      ! record is pass 3's.
      type(coupled_entry_snapshot), intent(inout) :: snap
      integer :: ipass
      logical :: w
      do ipass = 1, 3
         call write_the_pass(snap, ipass, .true., ipass, ipass - 1, w)
      enddo
      end subroutine three_bound_passes_then_a_stall

      ! --------------------------------------------------------------- !

      subroutine three_passes_off_the_bound_then_a_stall(snap)
      ! The same four passes with relaxations that ended elsewhere: the
      ! joint distance is just as flat and the bound was never met.
      type(coupled_entry_snapshot), intent(inout) :: snap
      integer :: ipass
      logical :: w
      do ipass = 1, 3
         call write_the_pass(snap, ipass, .false., 0, ipass - 1, w)
      enddo
      end subroutine three_passes_off_the_bound_then_a_stall

      ! --------------------------------------------------------------- !

      subroutine one_more_bound_pass(snap)
      type(coupled_entry_snapshot), intent(inout) :: snap
      logical :: w
      call write_the_pass(snap, 4, .true., 4, 3, w)
      end subroutine one_more_bound_pass

      ! --------------------------------------------------------------- !

      subroutine one_more_bound_pass_written(snap, w)
      type(coupled_entry_snapshot), intent(inout) :: snap
      logical, intent(out) :: w
      call write_the_pass(snap, 5, .true., 5, 4, w)
      end subroutine one_more_bound_pass_written

      ! --------------------------------------------------------------- !

      subroutine check_int(name, measured, reference)
      character(*), intent(in) :: name
      integer,      intent(in) :: measured, reference
      if (measured .eq. reference) then
         write(*,'(A,A,A,I0,A,I0,A)') 'PASS ', name, ' measured=',        &
              measured, ' reference=', reference, ' tol=0'
      else
         write(*,'(A,A,A,I0,A,I0,A)') 'FAIL ', name, ' measured=',        &
              measured, ' reference=', reference, ' tol=0'
         n_fail = n_fail + 1
      endif
      end subroutine check_int

      ! --------------------------------------------------------------- !

      subroutine check_log(name, measured, reference)
      character(*), intent(in) :: name
      logical,      intent(in) :: measured, reference
      character(5) :: sm, sr
      sm = 'false';  sr = 'false'
      if (measured)  sm = 'true '
      if (reference) sr = 'true '
      if (measured .eqv. reference) then
         write(*,'(A,A,A,A,A,A,A)') 'PASS ', name, ' measured=',          &
              trim(sm), ' reference=', trim(sr), ' tol=0'
      else
         write(*,'(A,A,A,A,A,A,A)') 'FAIL ', name, ' measured=',          &
              trim(sm), ' reference=', trim(sr), ' tol=0'
         n_fail = n_fail + 1
      endif
      end subroutine check_log

      ! --------------------------------------------------------------- !

      subroutine check_real(name, measured, reference, tol)
      ! tol = 0 asks for the two numbers to agree bit for bit, which is
      ! what one operator at one set of bits gives; a positive tol is an
      ! absolute bound.
      character(*), intent(in) :: name
      real*8,       intent(in) :: measured, reference, tol
      logical :: ok
      if (tol .le. 0.0d0) then
         ok = (measured .eq. reference)
      else
         ok = (abs(measured - reference) .le. tol)
      endif
      if (ok) then
         write(*,'(A,A,A,ES22.15,A,ES22.15,A,ES9.2)') 'PASS ', name,      &
              ' measured=', measured, ' reference=', reference, ' tol=',  &
              tol
      else
         write(*,'(A,A,A,ES22.15,A,ES22.15,A,ES9.2)') 'FAIL ', name,      &
              ' measured=', measured, ' reference=', reference, ' tol=',  &
              tol
         n_fail = n_fail + 1
      endif
      end subroutine check_real

      end program coupled_block_handover_contract
