      module boundary_state_trace_mod
      ! THE LOWER BOUNDARY AN EVALUATION STANDS ON, MEASURED. Every routine
      ! here reads state and writes a number; none decides anything, and all
      ! of it is off unless EXHALE_BOUNDARY_TRACE is set.
      !
      ! THE INVARIANT IT MEASURES: one composition, one boundary, one
      ! residual. The base face state, the ghost cell averages and the ghost
      ! conserved state are functions of the interior cell averages and of
      ! the installed composition, through the caloric map the ghost
      ! continuation uses (base_boundary.f90, continue_hydrostatic_isentrope
      ! and adiabatic_index_at_T at the ghost index). Apply_BC derives them
      ! and CACHES them in BC_Apply (base_face_W, base_ghost_W,
      ! base_face_lower_W), and Rec_BC puts base_face_W straight into the
      ! left slot of the base face. A residual whose sources, pressure map
      ! and sound speeds belong to one composition and whose cached face
      ! state belongs to another is therefore a residual of no single gas.
      !
      ! base_face_state_of_installed_state returns the distance between the
      ! cached boundary and the boundary the installed state would give now:
      ! zero exactly where the invariant holds.

      use global_parameters
      use Conversion, only: U_to_W_interior
      use BC_Apply, only: base_face_W, base_ghost_W, base_face_lower_W
      use base_boundary, only: base_boundary_states,                      &
                               ghost_state_ntot, ghost_state_ne,          &
                               ghost_state_x_h2, ghost_state_x_h2p,       &
                               ghost_state_x_h3p, ghost_state_x_hehp,     &
                               ghost_state_x_hii, ghost_state_x_heii,     &
                               ghost_state_x_heiii,                       &
                               ghost_state_charge_gap,                    &
                               ghost_state_element_gap,                   &
                               ghost_state_reaction_res,                  &
                               ghost_state_partition_res,                 &
                               ghost_state_closure_passes,                &
                               n_base_interior_inadmissible,              &
                               n_base_face_mach_limited,                  &
                               n_base_reversal_evals,                     &
                               n_base_supersonic_outflow,                 &
                               base_face_mach_last, base_face_blend_last, &
                               base_face_Mi_last, base_face_Mwind_last,   &
                               base_face_swind_last

      implicit none
      private
      public :: boundary_trace_armed, boundary_rebuild_suppressed,        &
                base_face_state_of_installed_state,                       &
                report_base_face_state_consistency, base_mass_row_trace,  &
                ghost_record_armed, write_ghost_record

      integer, parameter :: trace_unit = 78
      logical, save      :: trace_opened = .false.
      integer, parameter :: ghost_record_unit = 79
      logical, save      :: ghost_record_opened = .false.

      contains

      ! ------------------------------------------------------!

      logical function boundary_trace_armed() result(on)
      ! EXHALE_BOUNDARY_TRACE, default off.
      character(len=32) :: env
      call get_environment_variable('EXHALE_BOUNDARY_TRACE', env)
      on = (len_trim(env) .gt. 0 .and. trim(env) .ne. '0')
      end function boundary_trace_armed

      ! ------------------------------------------------------!

      logical function boundary_rebuild_suppressed() result(off)
      ! EXHALE_TRACE_EXPERIMENT=boundary_from_the_entry_composition, default
      ! off: the evaluation routes leave the boundary of the composition
      ! they were ENTERED with standing, instead of deriving it again from
      ! the composition the sweep returned. It exists so that one binary can
      ! measure both orders of the same evaluation, which is what the
      ! boundary_state suite asserts on; with the key unset the boundary is
      ! always derived from the composition the residual is assembled with.
      character(len=64) :: env
      call get_environment_variable('EXHALE_TRACE_EXPERIMENT', env)
      off = (index(env, 'boundary_from_the_entry_composition') .gt. 0)
      end function boundary_rebuild_suppressed

      ! ------------------------------------------------------!

      subroutine open_trace()
      if (trace_opened) then
         open(unit = trace_unit, file = './output/boundary_trace.txt',    &
              status = 'old', position = 'append')
      else
         open(unit = trace_unit, file = './output/boundary_trace.txt')
         write(trace_unit,'(A)') '# the lower boundary at named points'// &
              ' of one evaluation'
         trace_opened = .true.
      endif
      end subroutine open_trace

      ! ------------------------------------------------------!

      real*8 function base_face_state_of_installed_state(u_in)            &
             result(dmax)
      ! THE CACHED BOUNDARY AGAINST THE BOUNDARY OF THE STATE INSTALLED NOW,
      ! as the largest relative difference over the base face state, the
      ! ghost cell averages and the state at the next face down.
      !
      ! The derivation is the one Apply_BC performs: the interior cell
      ! averages in primitive form, then base_boundary_states. It is
      ! evaluated into local arrays, so the cache the residual reads is not
      ! touched; the branch diagnostics base_boundary_states keeps of its
      ! last evaluation are saved and put back, so a measurement does not
      ! enter the counts a run reports.
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u_in
      real*8, dimension(3,1-Ng:N+Ng) :: W_now
      real*8 :: Wf(3), Wg(3,1-Ng:0), Wl(3)
      integer :: n_inad, n_mach, n_rev, n_sup
      real*8  :: m_last, b_last, mi_last, mw_last, sw_last
      integer :: k, j

      n_inad  = n_base_interior_inadmissible
      n_mach  = n_base_face_mach_limited
      n_rev   = n_base_reversal_evals
      n_sup   = n_base_supersonic_outflow
      m_last  = base_face_mach_last
      b_last  = base_face_blend_last
      mi_last = base_face_Mi_last
      mw_last = base_face_Mwind_last
      sw_last = base_face_swind_last

      call U_to_W_interior(u_in, W_now)
      call base_boundary_states(W_now, Wf, Wg, Wl)

      n_base_interior_inadmissible = n_inad
      n_base_face_mach_limited     = n_mach
      n_base_reversal_evals        = n_rev
      n_base_supersonic_outflow    = n_sup
      base_face_mach_last          = m_last
      base_face_blend_last         = b_last
      base_face_Mi_last            = mi_last
      base_face_Mwind_last         = mw_last
      base_face_swind_last         = sw_last

      dmax = 0.0d0
      do k = 1, 3
         dmax = max(dmax, relative_gap(base_face_W(k), Wf(k)))
         dmax = max(dmax, relative_gap(base_face_lower_W(k), Wl(k)))
         do j = 1-Ng, 0
            dmax = max(dmax, relative_gap(base_ghost_W(k,j), Wg(k,j)))
         enddo
      enddo
      end function base_face_state_of_installed_state

      ! ------------------------------------------------------!

      real*8 function relative_gap(a, b) result(d)
      real*8, intent(in) :: a, b
      real*8 :: s
      s = max(abs(a), abs(b))
      d = 0.0d0
      if (s .gt. 0.0d0) d = abs(a - b)/s
      end function relative_gap

      ! ------------------------------------------------------!

      subroutine report_base_face_state_consistency(label, u_in)
      ! One line naming the point of the evaluation and the gap above.
      character(len=*),               intent(in) :: label
      real*8, dimension(3,1-Ng:N+Ng), intent(in) :: u_in
      real*8 :: d
      if (.not. boundary_trace_armed()) return
      d = base_face_state_of_installed_state(u_in)
      write(*,'(A,A,A,ES12.5)') ' [boundary trace] cached base face'//    &
           ' state against the boundary of the installed composition (',  &
           trim(label), '): ', d
      call open_trace()
      write(trace_unit,'(A,A,1X,ES26.16E3)') 'face_state_gap ',           &
           trim(label), d
      close(trace_unit)
      end subroutine report_base_face_state_consistency

      ! ------------------------------------------------------!

      subroutine base_mass_row_trace(label, jlast, R1, Ffac, sc, fl)
      ! The signed continuity row of the base cells, the two face mass
      ! fluxes it differences, its scale and its rounding floor. The caller
      ! owns the definitions; this routine writes them so that two
      ! evaluations can be compared cell by cell as dimensional quantities.
      character(len=*), intent(in) :: label
      integer,          intent(in) :: jlast
      real*8,           intent(in) :: R1(1:jlast), Ffac(0:jlast)
      real*8,           intent(in) :: sc(1:jlast), fl(1:jlast)
      integer :: j
      if (.not. boundary_trace_armed()) return
      call open_trace()
      write(trace_unit,'(A,A)') '# massrow ', trim(label)
      write(trace_unit,'(A)') '# j R1_signed F_lower F_upper scale floor'
      do j = 1, jlast
         write(trace_unit,'(A,I4,5ES26.16E3)') 'row ', j, R1(j),          &
              Ffac(j-1), Ffac(j), sc(j), fl(j)
      enddo
      close(trace_unit)
      end subroutine base_mass_row_trace

      ! ------------------------------------------------------!

      logical function ghost_record_armed() result(on)
      ! EXHALE_GHOST_RECORD, default off. When it is set every evaluation
      ! appends one record to ./output/ghost_record.txt.
      character(len=32) :: env
      call get_environment_variable('EXHALE_GHOST_RECORD', env)
      on = (len_trim(env) .gt. 0 .and. trim(env) .ne. '0')
      end function ghost_record_armed

      ! ------------------------------------------------------!

      subroutine write_ghost_record(label, p_g, rho_g, T_g, base_row,     &
                                    F_lower, F_upper, row_scale, row_floor)
      ! ONE RECORD PER EVALUATION: everything a controlled attribution has to
      ! compare between two evaluations of one state, as dimensional
      ! quantities over a stated scale.
      !
      ! Per lower ghost cell: the thermodynamic row the boundary left there
      ! (p, rho, T), the heavy-particle and electron counts of the
      ! composition the sweep returned, the molecular and ionization
      ! fractions of that composition, the charge and elemental budgets, the
      ! reaction residual the cell was ACCEPTED at, the residual of its
      ! imposed molecular partition at the returned state, and the pass the
      ! closure ended on.
      !
      ! And, in the same record, the base cell's continuity row with the two
      ! face mass fluxes it is built from. base_row is the signed row of cell
      ! 1 as the residual holds it; F_lower is the mass flux through the base
      ! face r_edg(0) and F_upper the flux through the face above it, in the
      ! same units and the same definition base_mass_row_trace writes, and
      ! their signed difference is the cancellation the row rests on. The
      ! row is judged as |base_row|/row_scale, and row_floor is the rounding
      ! floor in those same units (steady_residual.f90,
      ! mass_row_rounding_floor), so a row below its floor is the arithmetic
      ! of its terms and not a balance. The row is NOT F_upper - F_lower:
      ! the face areas and the cell volume of the spherical geometry stand
      ! between them, and the record carries both so that a move in either
      ! can be seen.
      !
      ! Off unless EXHALE_GHOST_RECORD is set; it reads state and decides
      ! nothing.
      character(len=*), intent(in) :: label
      real*8, intent(in) :: p_g(1-Ng:0), rho_g(1-Ng:0), T_g(1-Ng:0)
      real*8, intent(in) :: base_row, F_lower, F_upper, row_scale, row_floor
      integer :: j
      if (.not. ghost_record_armed()) return
      if (ghost_record_opened) then
         open(unit = ghost_record_unit, file = './output/ghost_record.txt',&
              status = 'old', position = 'append')
      else
         open(unit = ghost_record_unit, file = './output/ghost_record.txt')
         write(ghost_record_unit,'(A)') '# the lower ghost of one'//      &
              ' evaluation, and the base continuity row it produces'
         write(ghost_record_unit,'(A)') '# ghost <label> <cell> p rho'//  &
              ' T n_tot n_e x_H2 x_H2p x_H3p x_HeHp x_HII x_HeII'//       &
              ' x_HeIII charge_budget element_budget reaction_res'//      &
              ' partition_res closure_passes'
         write(ghost_record_unit,'(A)') '# base  <label> row_cell1'//    &
              ' F_lower F_upper F_upper_minus_F_lower row_scale row_floor'
         ghost_record_opened = .true.
      endif
      do j = 1-Ng, 0
         write(ghost_record_unit,'(A,1X,A,1X,I4,16ES26.16E3,1X,I6)')      &
              'ghost', trim(label), j, p_g(j), rho_g(j), T_g(j),          &
              ghost_state_ntot(j), ghost_state_ne(j),                     &
              ghost_state_x_h2(j), ghost_state_x_h2p(j),                  &
              ghost_state_x_h3p(j), ghost_state_x_hehp(j),                &
              ghost_state_x_hii(j), ghost_state_x_heii(j),                &
              ghost_state_x_heiii(j), ghost_state_charge_gap(j),          &
              ghost_state_element_gap(j), ghost_state_reaction_res(j),    &
              ghost_state_partition_res(j), ghost_state_closure_passes(j)
      enddo
      write(ghost_record_unit,'(A,1X,A,6ES26.16E3)') 'base ', trim(label),&
           base_row, F_lower, F_upper, F_upper - F_lower, row_scale,      &
           row_floor
      close(ghost_record_unit)
      end subroutine write_ghost_record

      end module boundary_state_trace_mod
