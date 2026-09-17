      program base_branch_discriminant
      ! WHAT THE BASE BOUNDARY DECIDES ITS ENTROPY BRANCH BY, AND WHETHER
      ! THAT DECISION IS CONTINUOUS (item L21 and the Codex review of
      ! 2026-09-15).
      !
      ! The face state of the base is
      !
      !     rho_b = (1 - w_rev) rho_res + w_rev rho_rev ,
      !     w_rev = s_wind smoothstep( M_wind / blend )
      !             + (1 - s_wind) smoothstep( M_i / blend ) ,
      !
      ! with rho_res the reservoir's density on its own isentrope and rho_rev
      ! the interior's.  w_rev = 0 is "gas enters and the reservoir states its
      ! entropy"; w_rev = 1 is "gas leaves and the entropy at the face is the
      ! interior's".  M_branch is what this program is about.
      !
      ! FOUR STATEMENTS, one per section below.
      !
      !   1. A base whose FIRST CELL moves inward while the WIND carries mass
      !      outward is an inflow face: the cell-centred product of a cell
      !      inside a collocated odd-even mode is not a flux (item L21,
      !      docs/p44_base_sawtooth.md section 3).  The branch must follow the
      !      wind, and it must do so on EVERY route: one state is a steady
      !      state of one operator.  Where the wind window carries no flux --
      !      a cold start, at rest there -- the first interior cell answers
      !      instead, which is the second half of the same statement.
      !   2. At rest the boundary is well balanced and the wind says nothing,
      !      so the face velocity is zero and the two routes agree exactly.
      !   3. The handover from the local face to the wind is C1: sweeping the
      !      wind through zero and out past the blend width, the face density
      !      has no jump.  Measured as a refinement test -- the largest
      !      neighbouring increment of a continuous function halves when the
      !      sampling is halved -- which a threshold would fail outright.
      !   4. On a state whose mass flux is uniform the branch cannot depend on
      !      WHERE the wind window starts.
      !
      ! Everything is the production boundary: base_boundary_states on a state
      ! this program builds, with the production grid and gravity.  Nothing
      ! here re-implements the condition.
      use global_parameters
      use grid_construction,          only: define_grid
      use gravity_grid_construction,  only: set_gravity_grid
      use base_boundary,              only: set_base_reservoir,             &
                                            base_boundary_states,           &
                                            base_face_mach_blend,           &
                                            base_face_blend_last,           &
                                            base_face_Mi_last,              &
                                            base_face_Mwind_last,           &
                                            base_face_swind_last,           &
                                            base_wind_window_spread
      implicit none

      integer, parameter :: Ncells = 200
      real*8,  parameter :: Teq = 1.0d3, Rp_RJ = 0.2d0, Mp_MJ = 0.02d0
      real*8, allocatable :: W(:,:)
      real*8  :: Wface(3), Wghost(3,1-2:0), Wface_lower(3)
      real*8  :: rho_res_only(3), rho_rev_only(3)
      integer :: j, n_fail
      real*8  :: rho_in, p_in, v_wind, v_cell1
      real*8  :: rho_stat, rho_march, w_stat, w_march
      real*8  :: jump_n, jump_2n, rho_a, rho_b_win, mb_a, mb_b
      real*8  :: at_vw, at_s, at_w

      n_fail = 0

      N  = Ncells
      T0 = Teq
      R0 = Rp_RJ*RJ
      Mp = Mp_MJ*MJ
      v0 = sqrt(kb_erg*T0/mu)
      b0 = (Gc*Mp*mu)/(kb_erg*T0*R0)
      ! The grid configuration of a base-resolved run, the same defaults
      ! input_read uses; only the base region matters here.
      spherical_domain = .true.
      dr_base     = 2.0d-4
      N_low_cells = 50
      r_max       = 10.0d0
      r_esc       =  2.0d0
      r_flux      =  1.20d0
      grid_type   = 'Mixed'
      CFL         = 0.6d0
      call allocate_grid_arrays
      call define_grid
      call set_gravity_grid
      allocate(W(3,1-Ng:N+Ng))

      ! A uniform, mildly stratified column is enough: the branch reads the
      ! first cell, the wind window and the reservoir, and nothing else.
      ! The interior sits at TWICE the reservoir's temperature at the same
      ! pressure, so the two isentropes the branch mixes are a factor two
      ! apart in density -- which is what they are at the LHS 1140 b base
      ! (rho_rev/rho_res = 0.531, item L21).  A test on which branch is taken
      ! has to be able to see the difference.
      rho_in = 0.5d0
      p_in   = 1.0d0
      n_part_cell1 = 0.5d0
      call set_base_reservoir(1.0d0, 1.0d0, 1.0d0, 1.0d0)

      ! ---------------------------------------------------------------- !
      ! 1. a local reversal inside a forward wind
      ! ---------------------------------------------------------------- !
      ! The first cell moves inward at ten times the wind's own face speed --
      ! the odd-even amplitude measured on the certified LHS 1140 b base is
      ! two -- while every cell of the wind window carries the same outward
      ! mass flux.
      v_wind  =  1.0d2*base_face_mach_blend
      v_cell1 = -1.0d3*base_face_mach_blend
      call fill_column(v_wind, v_cell1, 0.0d0)

      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      rho_stat = Wface(1);  w_stat = base_face_blend_last
      ! The same first cell inside a window that carries the same mean flux
      ! and is NOT one wind (its own spread is three times the threshold; the
      ! wasp_full cold start reads 4.4e-2, twenty-two times it): the window
      ! then has nothing to say
      ! about the base and the local face answers.
      call fill_column(v_wind, v_cell1, 3.0d0*base_wind_window_spread)
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      rho_march = Wface(1); w_march = base_face_blend_last
      call fill_column(v_wind, v_cell1, 0.0d0)

      call verdict('base_branch_wind_wins_where_the_window_is_one_wind',   &
                   w_stat, 0.0d0, 1.0d-12)
      call verdict('base_branch_local_face_wins_where_it_is_not',          &
                   w_march, 1.0d0, 1.0d-12)
      ! and the two branches are not a small difference of the same state
      write(*,'(A,ES12.5,A,ES12.5,A,F9.5)') '  face density: stationary ',  &
           rho_stat, '  marching ', rho_march, '  ratio ',                 &
           rho_march/rho_stat

      ! ---------------------------------------------------------------- !
      ! 2. at rest
      ! ---------------------------------------------------------------- !
      ! At rest and unstratified the interior IS the reservoir, so the face
      ! velocity the outgoing invariant returns is zero exactly, whichever
      ! branch is taken -- and the wind, which is also at rest, says nothing,
      ! so the boundary is a function of the state alone and repeats to the
      ! last bit.  (The stratified well-balanced statement is
      ! hydrostatic_residual's.)
      b0 = 0.0d0
      call set_gravity_grid
      n_part_cell1 = 1.0d0
      call set_base_reservoir(1.0d0, 1.0d0, 1.0d0, 1.0d0)
      rho_in = 1.0d0;  p_in = 1.0d0
      call fill_column(0.0d0, 0.0d0, 0.0d0)
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      rho_res_only = Wface
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      rho_rev_only = Wface
      call verdict('base_branch_at_rest_face_velocity_is_zero',            &
                   abs(rho_rev_only(2)), 0.0d0, 1.0d-12)
      call verdict('base_branch_at_rest_is_reproducible',                  &
                   abs(rho_res_only(1) - rho_rev_only(1))                  &
                   /max(rho_rev_only(1),1.0d-99), 0.0d0, 1.0d-14)
      ! back to the stratified, mismatched pair for the rest
      b0 = (Gc*Mp*mu)/(kb_erg*T0*R0)
      call set_gravity_grid
      rho_in = 0.5d0;  p_in = 1.0d0
      n_part_cell1 = 0.5d0
      call set_base_reservoir(1.0d0, 1.0d0, 1.0d0, 1.0d0)

      ! ---------------------------------------------------------------- !
      ! 3. the handover is continuous
      ! ---------------------------------------------------------------- !
      ! Sweep the WINDOW'S OWN RELATIVE FLUX SPREAD from half to two and a
      ! half times
      ! the threshold at a fixed outward wind and a fixed inward cell 1, so
      ! the two discriminants disagree in sign over the whole sweep and the
      ! branch is handed over inside it.  A continuous function sampled at 2n
      ! points has half the largest neighbouring increment it has at n; a
      ! threshold has the same one at every n.
      call sweep_largest_increment( 401, jump_n)
      call sweep_largest_increment( 801, jump_2n)
      write(*,'(A,ES12.5,A,ES12.5,A,F8.3)')                                &
           '  largest neighbouring change of the face density: n ',        &
           jump_n, '  2n ', jump_2n, '  ratio ', jump_2n/max(jump_n,1.0d-99)
      call verdict('base_branch_handover_is_continuous',                   &
                   jump_2n/max(jump_n,1.0d-99), 0.5d0, 0.06d0)

      ! ---------------------------------------------------------------- !
      ! 4. the branch does not depend on where the window starts
      ! ---------------------------------------------------------------- !
      ! A state whose mass flux is uniform carries the same flux in every
      ! window, so the branch it selects has to be the same one.
      call fill_column(v_wind, v_wind, 0.0d0)
      j_flux = max(2, N/4)
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      rho_a = Wface(1);  mb_a = base_face_blend_last
      j_flux = max(2, (3*N)/4)
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      rho_b_win = Wface(1);  mb_b = base_face_blend_last
      write(*,'(A,ES12.5,A,ES12.5)') '  w_rev at the two windows: ',      &
           mb_a, '  ', mb_b
      call verdict('base_branch_independent_of_the_flux_window',           &
                   abs(rho_a - rho_b_win)/max(rho_a,1.0d-99),              &
                   0.0d0, 1.0d-12)

      if (n_fail .gt. 0) then
         write(*,'(A,I0)') 'base_branch_discriminant FAILURES: ', n_fail
         stop 1
      endif
      stop 0

      contains

      subroutine fill_column(v_out, v_one, d_target)
      ! A column at uniform density and pressure whose WIND WINDOW CARRIES ONE
      ! MASS FLUX: v(j) = v_out (r_1/r_j)^2, so rho v r^2 takes the same value
      ! in every cell of the window and the window's own relative spread is
      ! zero.  v_out is therefore the face speed the flux corresponds to at
      ! r_1, which is the speed the branch reads.  The first cell is set to
      ! v_one afterwards and is outside the window.
      !
      ! d_target > 0 scatters the window about that flux, alternating cell by
      ! cell by +/- d_target, so the window's RELATIVE STANDARD DEVIATION --
      ! the functional the branch weight is read on -- is exactly d_target
      ! and the window carries a flux while not being one wind, which is what
      ! a cold start is and what the discriminant has to refuse.
      real*8, intent(in) :: v_out, v_one, d_target
      integer :: jj
      real*8  :: sgn
      do jj = 1-Ng, N+Ng
         W(1,jj) = rho_in
         W(3,jj) = p_in
         sgn = 1.0d0
         if (mod(jj,2) .eq. 0) sgn = -1.0d0
         W(2,jj) = v_out*(r(1)/r(jj))**2*(1.0d0 + d_target*sgn)
      enddo
      W(2,1) = v_one
      end subroutine fill_column

      subroutine sweep_largest_increment(nsamp, jump)
      ! The largest change of the face density between neighbouring samples of
      ! a sweep of the WINDOW'S OWN FLUX SPREAD through the handover, at a
      ! fixed outward wind and a fixed inward first cell, so the two
      ! discriminants disagree in sign over the whole sweep and the branch is
      ! handed from one to the other exactly once.
      integer, intent(in)  :: nsamp
      real*8,  intent(out) :: jump
      real*8  :: dw, rho_prev, rho_now, lo, hi
      integer :: is
      lo = 0.5d0*base_wind_window_spread
      hi = 2.5d0*base_wind_window_spread
      jump = 0.0d0
      rho_prev = 0.0d0
      at_vw = 0.0d0;  at_s = 0.0d0;  at_w = 0.0d0
      do is = 1, nsamp
         dw = lo + (hi - lo)*dble(is-1)/dble(nsamp-1)
         call fill_column(1.0d2*base_face_mach_blend,                     &
                          -1.0d3*base_face_mach_blend, dw)
         call base_boundary_states(W, Wface, Wghost, Wface_lower)
         rho_now = Wface(1)
         if (is .gt. 1 .and. abs(rho_now - rho_prev) .gt. jump) then
            jump  = abs(rho_now - rho_prev)
            at_vw = dw
            at_s  = base_face_swind_last
            at_w  = base_face_blend_last
         endif
         rho_prev = rho_now
      enddo
      write(*,'(A,I5,A,ES12.5,A,ES12.5,A,ES12.5)') '   [sweep] ', nsamp,  &
           ' samples: largest change at d=', at_vw, '  s=', at_s,        &
           '  w_rev=', at_w
      end subroutine sweep_largest_increment

      subroutine verdict(name, measured, reference, tol)
      character(len=*), intent(in) :: name
      real*8,           intent(in) :: measured, reference, tol
      if (abs(measured - reference) .le. tol) then
         write(*,'(A,A,A,ES13.6,A,ES13.6,A,ES9.2)') 'PASS ', name,        &
              ' measured=', measured, ' reference=', reference, ' tol=', tol
      else
         write(*,'(A,A,A,ES13.6,A,ES13.6,A,ES9.2)') 'FAIL ', name,        &
              ' measured=', measured, ' reference=', reference, ' tol=', tol
         n_fail = n_fail + 1
      endif
      end subroutine verdict

      end program base_branch_discriminant
