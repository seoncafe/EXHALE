      module coupled_block_handover
      ! THE STATE AT WHICH THE WIND-COMPOSITION ALTERNATION STOPS, AND THE
      ! ONE STATEMENT THAT LETS THE COUPLED CARRIER BLOCK BE ENTERED FROM
      ! IT.
      !
      ! WHAT THE ALTERNATION IS AND WHERE IT STOPS. The stationary outer
      ! iteration solves the wind at a held composition and then relaxes
      ! the composition at that held wind. The relaxation is bounded: one
      ! pass may move the composition by at most the movement bound, and a
      ! relaxation that meets the bound has moved as far as a held wind
      ! admits and no further. Where the joint distance of the state stops
      ! falling AND the relaxation keeps ending on that bound, the
      ! composition the carrier rows need lies outside the range in which
      ! holding the wind is admissible, and no shorter step reaches it.
      ! That is the state this module names, and the remedy is to stop
      ! holding the wind: the wind and the transported balances become one
      ! map of one vector.
      !
      ! WHY A STATE TRANSITION AND NOT A TEST ON THE LAST TOKEN. The pass
      ! that declares the stall takes no composition update, so on that
      ! pass no relaxation runs and the relaxation's outcome token holds
      ! "nothing to advance". A condition that reads the live token on that
      ! pass can never be met. The transition below is taken from a
      ! RECORDED state instead: the entry snapshot, written by the last
      ! pass in which a relaxation actually ran.
      !
      ! THE ONE COHERENT STATE THE COUPLED SOLVE CONSUMES is the pair
      ! (u, f_sp) the snapshot carries, with the primitive arrays rho, v,
      ! p, T and the face mass flux Frho_elem that belong to that pair.
      ! Both fixed-wind relaxations hold the CONSERVED state: the element
      ! relaxation is given rho and v and does not write them, and the
      ! carrier relaxation is given u and does not write it. So the
      ! conserved state the relaxation was entered at and the composition
      ! it returned are one state and not two thermodynamic stages; p and T
      ! are those of that pair, recomputed by the caller from u and f_sp
      ! after the relaxation. The residual evidence beside it was measured
      ! on that same pair. Nothing here is combined with a field of any
      ! other pass.
      !
      ! WHAT THE SNAPSHOT IS FOR, BESIDES THE STATE. The bound value, the
      ! exact ending of the relaxation, the cell the bound was refused at,
      ! the two residual norms with their own row terms and their worst
      ! cells, the progress reference of that pass and the two stall counts
      ! are recorded with it, so the transition is judged on the pass that
      ! produced the evidence and the log line names that pass.
      !
      ! THE SNAPSHOT IS SEALED AT THE TRANSITION. A generation counter is
      ! raised at every write. The acceptance below is given the generation
      ! the fresh residual was evaluated at and refuses a generation that
      ! has moved, so a residual measured on one snapshot can never
      ! authorize another. Once accepted, the snapshot is immutable: a
      ! later write is refused, and a second acceptance is refused, so the
      ! transition is taken exactly once per run.
      !
      ! WHAT THIS MODULE DOES NOT DO. It evaluates no residual. The caller
      ! owns the operators and evaluates them on the snapshot's own arrays;
      ! this module holds the state, states the transition and compares the
      ! fresh numbers with the recorded ones.
      implicit none
      private

      ! ---- THE STATES OF THE TRANSITION ----
      ! The alternation is running; it has stopped at the movement bound
      ! with the joint distance no longer falling; the coupled block has
      ! been entered from the recorded state.
      integer, parameter, public :: alternation_running          = 0
      integer, parameter, public :: alternation_stalled_at_bound = 1
      integer, parameter, public :: coupled_block_entry          = 2

      ! ---- WHY A TRANSITION WAS NOT TAKEN ----
      integer, parameter, public :: entry_refused_none             = 0
      integer, parameter, public :: entry_refused_option_off       = 1
      integer, parameter, public :: entry_refused_already_entered  = 2
      integer, parameter, public :: entry_refused_no_snapshot      = 3
      integer, parameter, public :: entry_refused_no_stall         = 4
      integer, parameter, public :: entry_refused_not_bound_ending = 5
      integer, parameter, public :: entry_refused_few_bound_ends   = 6
      integer, parameter, public :: entry_refused_no_carrier_rows  = 7
      integer, parameter, public :: entry_refused_snapshot_moved   = 8
      integer, parameter, public :: entry_refused_residual_moved   = 9
      integer, parameter, public :: entry_refused_residual_absent  = 10
      integer, parameter, public :: entry_refused_shape_mismatch   = 11

      ! HOW FAR THE RE-EVALUATED ENTRY RESIDUAL MAY STAND FROM THE
      ! RECORDED ONE. Both are the same operator at the same arrays, so
      ! they agree to rounding; the bound admits a different reduction
      ! order in a threaded sum and nothing more. It is a RELATIVE bound
      ! and never an equality of two floats: the same quantity formed by
      ! two orders of the same terms differs in the last bits.
      real*8, parameter, public :: coupled_entry_residual_rel_tol = 1.0d-6

      ! ---- THE EVIDENCE OF THE PASS THAT WROTE THE SNAPSHOT ----
      ! Everything measured on the recorded state, filled at the call site
      ! so that one write carries one pass.
      type, public :: alternation_stall_evidence
         ! The outer pass whose relaxation produced this state.
         integer :: pass = 0
         ! THE EXACT ENDING OF THE RELAXATION: the outcome token of the
         ! carrier relaxation and its text, and the statement the
         ! transition reads, which is whether that ending was the
         ! composition movement bound.
         integer            :: carrier_outcome = -1
         character(len=76)  :: carrier_outcome_text = 'not recorded'
         logical            :: ended_on_movement_bound = .false.
         ! The bound the pass ran under, and the cell it was refused at.
         real*8  :: movement_bound = 0.0d0
         integer :: bound_cell     = 0
         integer :: bound_carrier  = 0
         real*8  :: bound_measure  = 0.0d0
         real*8  :: bound_entry    = 0.0d0
         logical :: bound_is_fraction = .false.
         ! THE CARRIER BALANCE OF THE RECORDED COMPOSITION: the norm, the
         ! two absolute numbers behind it, and the cell and carrier that
         ! carry it. Unavailable is an outcome and is recorded as one.
         logical :: carrier_residual_available = .false.
         real*8  :: carrier_residual  = 0.0d0
         real*8  :: carrier_res_abs   = 0.0d0
         real*8  :: carrier_scale_abs = 0.0d0
         integer :: carrier_worst_cell    = 0
         integer :: carrier_worst_carrier = 0
         ! THE ELEMENTAL TRANSPORT BALANCE of the same composition, in the
         ! same form.
         logical :: element_residual_available = .false.
         real*8  :: element_residual  = 0.0d0
         real*8  :: element_res_abs   = 0.0d0
         real*8  :: element_scale_abs = 0.0d0
         integer :: element_worst_cell    = 0
         integer :: element_worst_element = 0
         ! THE OUTER PROGRESS REFERENCE of that pass (its joint distance,
         ! the largest certification entry over its own tolerance) and the
         ! two stall counts: consecutive passes without a fall, and
         ! consecutive relaxations ending on the bound.
         real*8  :: progress_reference = 0.0d0
         integer :: n_no_fall       = 0
         integer :: n_bound_endings = 0
         ! THE BOUNDARY AND RECONSTRUCTION STATE THE COUPLED SOLVE WILL
         ! CONSUME: whether the boundary is rebuilt from the composition
         ! at each assembly, and which discrete operator the faces are
         ! taken with, by the name of the operator the run is solving with
         ! and by the endpoint the assembly selects. All three belong to
         ! the recorded state because all three change the residual the
         ! block starts from.
         logical           :: boundary_rebuild_suppressed = .false.
         character(len=32) :: reconstruction_operator = 'not recorded'
         logical           :: reconstruction_is_plm = .false.
      end type alternation_stall_evidence

      ! ---- THE ENTRY SNAPSHOT ----
      type, public :: coupled_entry_snapshot
         logical :: recorded = .false.
         ! Sealed at the transition: nothing may write it afterwards.
         logical :: sealed = .false.
         ! Raised at every write. A fresh residual carries the generation
         ! it was evaluated at, and an acceptance at another generation is
         ! refused.
         integer :: generation = 0
         type(alternation_stall_evidence) :: evidence
         ! THE COMPLETE STATE, ghost cells included: the conserved state
         ! and the composition are the pair; the primitive arrays and the
         ! face mass flux are what the operators of that pair read.
         real*8, allocatable :: u(:,:)
         real*8, allocatable :: f_sp(:,:)
         real*8, allocatable :: rho(:)
         real*8, allocatable :: v(:)
         real*8, allocatable :: p(:)
         real*8, allocatable :: T(:)
         real*8, allocatable :: Frho_elem(:)
      end type coupled_entry_snapshot

      ! How many times a run has entered the coupled block from a recorded
      ! state. The transition is taken once; this is what says so.
      integer, save :: n_block_entries = 0

      public :: record_coupled_entry_state
      public :: coupled_entry_stall_transition
      public :: accept_coupled_block_entry
      public :: restore_coupled_entry_state
      public :: coupled_entry_state_fits
      public :: coupled_entry_refusal_text
      public :: n_coupled_block_entries
      public :: reset_coupled_block_entry_count

      contains

      ! --------------------------------------------------------------- !

      subroutine record_coupled_entry_state(snap, evidence, u, f_sp,      &
                                            rho, v, p, T, Frho_elem,      &
                                            written)
      ! WRITE THE ENTRY SNAPSHOT FROM THE PASS THAT JUST RAN A RELAXATION.
      ! One write carries one pass: the whole state and the whole evidence
      ! are replaced together, so no field of an earlier pass survives
      ! beside a field of this one. A sealed snapshot is not written and
      ! says so.
      type(coupled_entry_snapshot),     intent(inout) :: snap
      type(alternation_stall_evidence), intent(in)    :: evidence
      real*8, dimension(:,:), intent(in) :: u, f_sp
      real*8, dimension(:),   intent(in) :: rho, v, p, T, Frho_elem
      logical, intent(out) :: written

      written = .false.
      if (snap%sealed) return

      snap%evidence = evidence
      call copy_2d(snap%u,    u)
      call copy_2d(snap%f_sp, f_sp)
      call copy_1d(snap%rho,  rho)
      call copy_1d(snap%v,    v)
      call copy_1d(snap%p,    p)
      call copy_1d(snap%T,    T)
      call copy_1d(snap%Frho_elem, Frho_elem)
      snap%recorded   = .true.
      snap%generation = snap%generation + 1
      written         = .true.
      end subroutine record_coupled_entry_state

      ! --------------------------------------------------------------- !

      subroutine coupled_entry_stall_transition(snap, on_stall_enabled,   &
                       already_entered, carrier_rows_transported,         &
                       no_progress_declared, bound_endings_required,      &
                       state_now, why)
      ! IS THE ALTERNATION AT THE STATE THE COUPLED BLOCK IS ENTERED FROM?
      !
      ! The three conditions are the ones that document the stall, and
      ! every one of them is read from the RECORDED pass and not from the
      ! pass that asks:
      !   the joint distance did not fall on outer_no_fall_max consecutive
      !     passes, which is what the caller's no-progress ending says;
      !   the relaxation of the recorded pass ended on the composition
      !     movement bound;
      !   it had ended on that bound for bound_endings_required
      !     consecutive passes.
      ! A merely flat scalar norm is the first condition alone and does
      ! not reach the transition.
      type(coupled_entry_snapshot), intent(in)  :: snap
      logical, intent(in)  :: on_stall_enabled
      logical, intent(in)  :: already_entered
      logical, intent(in)  :: carrier_rows_transported
      logical, intent(in)  :: no_progress_declared
      integer, intent(in)  :: bound_endings_required
      integer, intent(out) :: state_now
      integer, intent(out) :: why

      state_now = alternation_running
      why       = entry_refused_none

      if (.not. on_stall_enabled) then
         why = entry_refused_option_off;        return
      endif
      if (already_entered) then
         why = entry_refused_already_entered;   return
      endif
      if (.not. carrier_rows_transported) then
         why = entry_refused_no_carrier_rows;   return
      endif
      if (.not. no_progress_declared) then
         why = entry_refused_no_stall;          return
      endif
      if (.not. snap%recorded) then
         why = entry_refused_no_snapshot;       return
      endif
      if (.not. snap%evidence%ended_on_movement_bound) then
         why = entry_refused_not_bound_ending;  return
      endif
      if (snap%evidence%n_bound_endings .lt. bound_endings_required) then
         why = entry_refused_few_bound_ends;    return
      endif

      state_now = alternation_stalled_at_bound
      end subroutine coupled_entry_stall_transition

      ! --------------------------------------------------------------- !

      subroutine accept_coupled_block_entry(snap, generation_evaluated,   &
                       state_fits_here,                                   &
                       carrier_residual_fresh, carrier_residual_ok,       &
                       element_residual_fresh, element_residual_ok,       &
                       state_now, why, carrier_gap, element_gap)
      ! TAKE THE TRANSITION, OR REFUSE IT ON THE RE-EVALUATED RESIDUAL.
      !
      ! The caller has evaluated the carrier balance, and the elemental
      ! one where the recorded pass had it, on the snapshot's OWN arrays.
      ! Those numbers are compared here with the ones the recorded pass
      ! measured. They are the same operator at the same state, so they
      ! agree to rounding; a disagreement says the state the transition
      ! was about is not the state that would be handed over, and the
      ! handover is refused rather than taken on a stale number.
      !
      ! generation_evaluated is the snapshot generation the fresh numbers
      ! were evaluated at. A generation that has moved means the snapshot
      ! was rewritten between the evaluation and here, so the two sets of
      ! fields belong to different passes and the acceptance is refused.
      type(coupled_entry_snapshot), intent(inout) :: snap
      integer, intent(in)  :: generation_evaluated
      ! Whether the recorded arrays fit the arrays of this run, asked
      ! BEFORE the transition is taken: a snapshot that cannot be handed
      ! back is not a state the block can be entered at, and a sealed
      ! snapshot that was never restored would say a handover happened
      ! where none did.
      logical, intent(in)  :: state_fits_here
      real*8,  intent(in)  :: carrier_residual_fresh
      logical, intent(in)  :: carrier_residual_ok
      real*8,  intent(in)  :: element_residual_fresh
      logical, intent(in)  :: element_residual_ok
      integer, intent(out) :: state_now
      integer, intent(out) :: why
      ! The relative distance between each re-evaluated norm and the
      ! recorded one, so a refusal is reported with its number.
      real*8,  intent(out) :: carrier_gap, element_gap

      state_now   = alternation_stalled_at_bound
      why         = entry_refused_none
      carrier_gap = 0.0d0
      element_gap = 0.0d0

      if (snap%sealed) then
         state_now = alternation_running
         why = entry_refused_already_entered;  return
      endif
      if (.not. snap%recorded) then
         state_now = alternation_running
         why = entry_refused_no_snapshot;      return
      endif
      if (.not. state_fits_here) then
         state_now = alternation_running
         why = entry_refused_shape_mismatch;   return
      endif
      if (generation_evaluated .ne. snap%generation) then
         why = entry_refused_snapshot_moved;   return
      endif

      ! The carrier balance is the row the stall is about, so it must be
      ! measurable on the entry state. Where the recorded pass could not
      ! measure it, nothing here can be verified.
      if (.not. snap%evidence%carrier_residual_available .or.             &
          .not. carrier_residual_ok) then
         why = entry_refused_residual_absent;  return
      endif
      carrier_gap = relative_gap(carrier_residual_fresh,                  &
                                 snap%evidence%carrier_residual)
      if (carrier_gap .gt. coupled_entry_residual_rel_tol) then
         why = entry_refused_residual_moved;   return
      endif

      if (snap%evidence%element_residual_available) then
         if (.not. element_residual_ok) then
            why = entry_refused_residual_absent;  return
         endif
         element_gap = relative_gap(element_residual_fresh,               &
                                    snap%evidence%element_residual)
         if (element_gap .gt. coupled_entry_residual_rel_tol) then
            why = entry_refused_residual_moved;   return
         endif
      endif

      snap%sealed     = .true.
      n_block_entries = n_block_entries + 1
      state_now       = coupled_block_entry
      end subroutine accept_coupled_block_entry

      ! --------------------------------------------------------------- !

      logical function coupled_entry_state_fits(snap, u, f_sp, rho, v, p, &
                                                T, Frho_elem) result(fits)
      ! DOES THE RECORDED STATE FIT THE ARRAYS OF THIS RUN? Asked before
      ! the transition is taken, so that the seal below is never set on a
      ! state that could not be handed back.
      type(coupled_entry_snapshot), intent(in) :: snap
      real*8, dimension(:,:), intent(in) :: u, f_sp
      real*8, dimension(:),   intent(in) :: rho, v, p, T, Frho_elem
      fits = .false.
      if (.not. snap%recorded)        return
      if (.not. allocated(snap%u))    return
      if (.not. allocated(snap%f_sp)) return
      if (size(snap%u,1)    .ne. size(u,1)    .or.                        &
          size(snap%u,2)    .ne. size(u,2))    return
      if (size(snap%f_sp,1) .ne. size(f_sp,1) .or.                        &
          size(snap%f_sp,2) .ne. size(f_sp,2)) return
      if (size(snap%rho) .ne. size(rho) .or. size(snap%v) .ne. size(v)    &
          .or. size(snap%p) .ne. size(p) .or. size(snap%T) .ne. size(T)   &
          .or. size(snap%Frho_elem) .ne. size(Frho_elem)) return
      fits = .true.
      end function coupled_entry_state_fits

      ! --------------------------------------------------------------- !

      subroutine restore_coupled_entry_state(snap, u, f_sp, rho, v, p, T, &
                                             Frho_elem, restored)
      ! HAND THE RECORDED STATE BACK, WHOLE. The coupled solve consumes
      ! (u, f_sp); the primitive arrays and the face mass flux are handed
      ! back with them so that every array the next assembly reads is of
      ! that one pair. A shape that does not match is not restored, and
      ! nothing is written.
      type(coupled_entry_snapshot), intent(in) :: snap
      real*8, dimension(:,:), intent(inout) :: u, f_sp
      real*8, dimension(:),   intent(inout) :: rho, v, p, T, Frho_elem
      logical, intent(out) :: restored

      restored = .false.
      if (.not. coupled_entry_state_fits(snap, u, f_sp, rho, v, p, T,     &
                                         Frho_elem)) return

      u         = snap%u
      f_sp      = snap%f_sp
      rho       = snap%rho
      v         = snap%v
      p         = snap%p
      T         = snap%T
      Frho_elem = snap%Frho_elem
      restored  = .true.
      end subroutine restore_coupled_entry_state

      ! --------------------------------------------------------------- !

      integer function n_coupled_block_entries() result(n)
      ! How many times this run entered the coupled block from a recorded
      ! state.
      n = n_block_entries
      end function n_coupled_block_entries

      ! --------------------------------------------------------------- !

      subroutine reset_coupled_block_entry_count()
      ! Only a test resets it: within a run the count is the run's.
      n_block_entries = 0
      end subroutine reset_coupled_block_entry_count

      ! --------------------------------------------------------------- !

      function coupled_entry_refusal_text(why) result(text)
      integer, intent(in) :: why
      character(len=96) :: text
      select case (why)
      case (entry_refused_none)
         text = 'nothing refused it'
      case (entry_refused_option_off)
         text = '"Coupled carrier solve: On stall" is not set'
      case (entry_refused_already_entered)
         text = 'the coupled block was already entered'
      case (entry_refused_no_snapshot)
         text = 'no pass of this run ran a relaxation, so there is no'//  &
                ' entry state'
      case (entry_refused_no_stall)
         text = 'the joint distance of the state is still falling'
      case (entry_refused_not_bound_ending)
         text = 'the last relaxation did not end on the composition'//    &
                ' movement bound'
      case (entry_refused_few_bound_ends)
         text = 'too few consecutive relaxations ended on the movement'// &
                ' bound'
      case (entry_refused_no_carrier_rows)
         text = 'this configuration transports no carrier'
      case (entry_refused_snapshot_moved)
         text = 'the entry snapshot was rewritten after the residual'//   &
                ' was evaluated on it'
      case (entry_refused_residual_moved)
         text = 'the residual re-evaluated on the entry state is not'//   &
                ' the one the stall was declared on'
      case (entry_refused_residual_absent)
         text = 'the carrier balance has no measurement on the entry'//   &
                ' state'
      case (entry_refused_shape_mismatch)
         text = 'the recorded state does not fit the arrays of this run'
      case default
         text = 'unnamed'
      end select
      end function coupled_entry_refusal_text

      ! --------------------------------------------------------------- !

      real*8 function relative_gap(a, b) result(g)
      ! |a - b| over the larger of the two magnitudes: the distance
      ! between two measurements of one quantity, which is what a
      ! tolerance on a residual norm is written against. Two zeros are at
      ! no distance.
      real*8, intent(in) :: a, b
      real*8 :: s
      s = max(abs(a), abs(b))
      if (s .le. 0.0d0) then
         g = 0.0d0
      else
         g = abs(a - b)/s
      endif
      end function relative_gap

      ! --------------------------------------------------------------- !

      subroutine copy_1d(dst, src)
      real*8, allocatable, intent(inout) :: dst(:)
      real*8,              intent(in)    :: src(:)
      if (allocated(dst)) then
         if (size(dst) .ne. size(src)) deallocate(dst)
      endif
      if (.not. allocated(dst)) allocate(dst(size(src)))
      dst = src
      end subroutine copy_1d

      ! --------------------------------------------------------------- !

      subroutine copy_2d(dst, src)
      real*8, allocatable, intent(inout) :: dst(:,:)
      real*8,              intent(in)    :: src(:,:)
      if (allocated(dst)) then
         if (size(dst,1) .ne. size(src,1) .or.                            &
             size(dst,2) .ne. size(src,2)) deallocate(dst)
      endif
      if (.not. allocated(dst)) allocate(dst(size(src,1),size(src,2)))
      dst = src
      end subroutine copy_2d

      end module coupled_block_handover
