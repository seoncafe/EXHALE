      program base_branch_discriminant
      ! WHAT THE BASE BOUNDARY DECIDES THE DIRECTION OF ITS CONTACT BY, AND
      ! WHAT THE TWO DIRECTIONS GIVE.
      !
      ! The face state of the base is
      !
      !     rho_b = (1 - w_rev) rho_res + w_rev rho_rev ,     w_rev in {0,1},
      !
      ! with rho_res the reservoir's density on its own isentrope and rho_rev
      ! the interior's.  w_rev = 0 is "the reservoir states the entropy of the
      ! level"; w_rev = 1 is "the interior's trace leaves through the face".
      ! The weight takes no other value: the contact is upwinded on the
      ! direction of the mass flux the level carries, which is a discrete
      ! fact about a state.  What that direction is READ FROM is what this
      ! program is about.
      !
      ! FIVE STATEMENTS, one per section below.
      !
      !   1. A base whose FIRST CELL moves inward while the WIND carries mass
      !      outward is an inflow face: the cell-centred product of a cell
      !      inside a collocated odd-even mode is not a flux (item L21,
      !      docs/p44_base_sawtooth.md section 3).  The direction follows the
      !      wind wherever the wind window carries one flux.
      !   2. Where the window carries no single flux -- a cold start, a
      !      column at rest, a state that is not a wind -- the direction is
      !      the sign of the MATCHED FACE VELOCITY and not of the first
      !      interior cell's.  The two are told apart here: the cell moves
      !      outward while the pressure mismatch across the level drives the
      !      face inward, and the face is what the boundary follows.
      !   3. At a pressure-balanced stationary contact the face velocity is
      !      zero and the RESERVOIR owns the level: w_rev = 0 and the face
      !      carries the reservoir's density, not an average of the two
      !      isentropes.
      !   4. Over a sweep of the window's own flux spread the branch takes
      !      only its two values and the face density only the two
      !      isentropes: there is no state between them.
      !   5. On a state whose mass flux is uniform the branch cannot depend on
      !      WHERE the wind window starts.
      !
      ! THE COLUMN every section is built on is a HYDROSTATIC ISENTROPE
      ! THROUGH THE RESERVOIR'S PRESSURE AT THE FACE, at twice the
      ! reservoir's temperature there.  Two properties make it the right
      ! fixture: the two isentropes the branch chooses between are a factor
      ! two apart in density, so a test can see which was taken; and the
      ! interior's pressure at the face is the reservoir's, so the matched
      ! face velocity equals the interior's own velocity there and a
      ! deliberate pressure mismatch is the only thing that separates them.
      !
      ! Everything is the production boundary: base_boundary_states on a state
      ! this program builds, with the production grid and gravity.  Nothing
      ! here re-implements the condition.
      use global_parameters
      use grid_construction,          only: define_grid
      use gravity_grid_construction,  only: set_gravity_grid
      use base_boundary,              only: set_base_reservoir,             &
                                            base_boundary_states,           &
                                            continue_hydrostatic_isentrope, &
                                            base_face_blend_last,           &
                                            base_face_Mi_last,              &
                                            base_face_rho_res_last,         &
                                            base_face_rho_rev_last,         &
                                            base_wind_window_spread
      implicit none

      integer, parameter :: Ncells = 200
      real*8,  parameter :: Teq = 1.0d3, Rp_RJ = 0.2d0, Mp_MJ = 0.02d0
      real*8,  parameter :: p_res_level = 1.0d0, T_res_level = 1.0d0
      real*8,  parameter :: nhat_res    = 1.0d0
      ! The interior stands at this multiple of the level's temperature at
      ! the same pressure, so the two isentropes are a factor two apart.
      real*8,  parameter :: T_interior_over_level = 2.0d0
      ! The pressure mismatch of section 2, large enough that the acoustic
      ! response it drives at the face outweighs the cell velocity beside it
      ! and small enough to stay a linear perturbation.
      real*8,  parameter :: dp_over_p = 1.0d-5

      real*8, allocatable :: W(:,:)
      real*8  :: Wface(3), Wghost(3,1-2:0), Wface_lower(3)
      real*8  :: rho_res_face, T_res_face, p_res_face, rho_face, T_face
      integer :: n_fail
      real*8  :: v_wind, v_cell1
      real*8  :: w_window, w_local, rho_window, rho_local
      real*8  :: v_b_rest, w_rest, rho_rest
      real*8  :: n_values, rho_a, rho_b_win, mb_a, mb_b
      real*8  :: mi_local

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

      call set_base_reservoir(p_res_level, T_res_level, nhat_res, 1.0d0)
      call continue_hydrostatic_isentrope(0, nhat_res, 1.0d0,              &
               p_res_level/(nhat_res*T_res_level), T_res_level, r_edg(0),  &
               rho_res_face, T_res_face)
      p_res_face = nhat_res*rho_res_face*T_res_face
      T_face     = T_interior_over_level*T_res_face
      rho_face   = p_res_face/(nhat_res*T_face)

      ! The face speed the sections below drive the column at, in units of
      ! the reservoir's own sound speed at the face: small enough to be a
      ! weak flow and far above the rounding floor of the flux.
      v_wind  =  1.0d-6*sqrt(p_res_face/rho_res_face)
      v_cell1 = -1.0d1*v_wind

      ! ---------------------------------------------------------------- !
      ! 1. a local reversal inside a forward wind
      ! ---------------------------------------------------------------- !
      ! The first cell moves inward at ten times the wind's own face speed
      ! -- the odd-even amplitude measured on the certified LHS 1140 b base
      ! is two -- while every cell of the wind window carries the same
      ! outward mass flux.  The window has standing, so the face is an
      ! inflow and the reservoir states the entropy.
      call fill_column(v_wind, v_cell1, 0.0d0, 1.0d0)
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      w_window = base_face_blend_last;  rho_window = Wface(1)
      call verdict('base_branch_wind_wins_where_the_window_is_one_wind',   &
                   w_window, 0.0d0, 0.0d0)
      call verdict('base_branch_inflow_face_is_the_reservoir',             &
                   abs(rho_window - base_face_rho_res_last)                &
                   /base_face_rho_res_last, 0.0d0, 1.0d-13)

      ! ---------------------------------------------------------------- !
      ! 2. where the window says nothing, the matched face velocity answers
      ! ---------------------------------------------------------------- !
      ! The same wind speed, scattered cell by cell so the window's own
      ! relative spread is three times the threshold and it carries no
      ! single flux (the wasp_full cold start reads 4.4e-2, twenty-two times
      ! it).  The first cell now moves OUTWARD at the wind's own speed while
      ! the interior pressure stands one part in 1e5 above the reservoir's,
      ! and that mismatch drives the face inward: a condition keyed on the
      ! cell would read an inflow, the matched face velocity reads a reverse
      ! flow, and the interior's trace is what leaves.
      call fill_column(v_wind, v_wind, 3.0d0*base_wind_window_spread,      &
                       1.0d0 + dp_over_p)
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      w_local = base_face_blend_last;  rho_local = Wface(1)
      mi_local = base_face_Mi_last
      write(*,'(A,ES12.5,A,ES12.5)') '  cell 1 at the face: M_i =',        &
           mi_local, '   matched face velocity v_b =', Wface(2)
      call verdict('base_branch_cell_one_is_an_inflow_here',               &
                   sign(1.0d0, mi_local), 1.0d0, 0.0d0)
      call verdict('base_branch_matched_face_velocity_is_a_reverse_flow',  &
                   sign(1.0d0, Wface(2)), -1.0d0, 0.0d0)
      call verdict('base_branch_local_face_wins_where_the_window_is_not'// &
                   '_one_wind', w_local, 1.0d0, 0.0d0)
      call verdict('base_branch_reverse_face_is_the_interior_trace',       &
                   abs(rho_local - base_face_rho_rev_last)                 &
                   /base_face_rho_rev_last, 0.0d0, 1.0d-13)
      ! and the two branches are not a small difference of the same state
      write(*,'(A,ES12.5,A,ES12.5,A,F9.5)') '  face density: inflow ',     &
           rho_window, '  reverse ', rho_local, '  ratio ',                &
           rho_local/rho_window

      ! ---------------------------------------------------------------- !
      ! 3. the pressure-balanced stationary contact
      ! ---------------------------------------------------------------- !
      ! The column at rest, its own isentrope through the reservoir's
      ! pressure at the face: the contact is stationary and the two sides
      ! carry two entropies.  The reservoir owns the level.
      call fill_column(0.0d0, 0.0d0, 0.0d0, 1.0d0)
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      v_b_rest = Wface(2);  w_rest = base_face_blend_last
      rho_rest = Wface(1)
      call verdict('base_branch_at_rest_face_velocity_is_zero',            &
                   abs(v_b_rest), 0.0d0, 1.0d-14)
      call verdict('base_branch_at_rest_the_reservoir_owns_the_level',     &
                   w_rest, 0.0d0, 0.0d0)
      call verdict('base_branch_at_rest_face_is_the_reservoir',            &
                   abs(rho_rest - base_face_rho_res_last)                  &
                   /base_face_rho_res_last, 0.0d0, 1.0d-13)

      ! ---------------------------------------------------------------- !
      ! 4. the branch has two values and no state between them
      ! ---------------------------------------------------------------- !
      ! Sweep the window's own relative flux spread from half to two and a
      ! half times the threshold, at a fixed outward wind and a face the
      ! pressure mismatch drives inward, so the two discriminants disagree
      ! in sign over the whole sweep and the branch is handed from one to
      ! the other exactly once.  The face density must take the two
      ! isentropes and nothing else.
      call sweep_distinct_face_densities(801, n_values)
      call verdict('base_branch_sweep_takes_only_the_two_isentropes',      &
                   n_values, 2.0d0, 0.0d0)

      ! ---------------------------------------------------------------- !
      ! 5. the branch does not depend on where the window starts
      ! ---------------------------------------------------------------- !
      ! A state whose mass flux is uniform carries the same flux in every
      ! window, so the branch it selects has to be the same one.
      call fill_column(v_wind, v_wind, 0.0d0, 1.0d0)
      j_flux = max(2, N/4)
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      rho_a = Wface(1);  mb_a = base_face_blend_last
      j_flux = max(2, (3*N)/4)
      call base_boundary_states(W, Wface, Wghost, Wface_lower)
      rho_b_win = Wface(1);  mb_b = base_face_blend_last
      write(*,'(A,ES12.5,A,ES12.5)') '  w_rev at the two windows: ',       &
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

      subroutine fill_column(v_out, v_one, d_target, p_fac)
      ! The hydrostatic isentrope through the reservoir's pressure at the
      ! face at T_interior_over_level times its temperature there, carrying a
      ! UNIFORM MASS FLUX: v(j) = v_out rho_1 r_1^2/(rho_j r_j^2), so
      ! rho v r^2 takes the same value in every cell and the window's own
      ! relative spread is zero.  v_out is therefore the face speed the flux
      ! corresponds to at cell 1.  The first cell is set to v_one afterwards.
      !
      ! d_target > 0 scatters the window about that flux, alternating cell by
      ! cell by +/- d_target, so the window's RELATIVE STANDARD DEVIATION --
      ! the functional the window's standing is read on -- is exactly
      ! d_target and the window carries a flux while not being one wind,
      ! which is what a cold start is.
      !
      ! p_fac scales the pressure of every cell, which moves the interior's
      ! pressure at the face off the reservoir's and is the only thing that
      ! separates the matched face velocity from the first cell's own.
      real*8, intent(in) :: v_out, v_one, d_target, p_fac
      integer :: jj
      real*8  :: sgn, rho_q, T_q, flux1
      do jj = 1-Ng, N+Ng
         call continue_hydrostatic_isentrope(0, nhat_res, r_edg(0),        &
                  rho_face, T_face, r(jj), rho_q, T_q)
         W(1,jj) = rho_q
         W(3,jj) = p_fac*nhat_res*rho_q*T_q
      enddo
      flux1 = W(1,1)*v_out*r(1)*r(1)
      do jj = 1-Ng, N+Ng
         sgn = 1.0d0
         if (mod(jj,2) .eq. 0) sgn = -1.0d0
         W(2,jj) = flux1*(1.0d0 + d_target*sgn)/(W(1,jj)*r(jj)*r(jj))
      enddo
      W(2,1) = v_one
      n_part_cell1 = nhat_res*W(1,1)
      end subroutine fill_column

      subroutine sweep_distinct_face_densities(nsamp, ndist)
      ! How many distinct face densities the sweep of section 4 produces.
      ! A branch that takes two values gives two; a handover of finite width
      ! gives as many as there are samples inside it.
      integer, intent(in)  :: nsamp
      real*8,  intent(out) :: ndist
      real*8  :: dw, lo, hi, seen(8)
      integer :: is, nseen, k
      logical :: known
      lo = 0.5d0*base_wind_window_spread
      hi = 2.5d0*base_wind_window_spread
      nseen = 0;  seen = 0.0d0
      do is = 1, nsamp
         dw = lo + (hi - lo)*dble(is-1)/dble(nsamp-1)
         call fill_column(v_wind, v_wind, dw, 1.0d0 + dp_over_p)
         call base_boundary_states(W, Wface, Wghost, Wface_lower)
         known = .false.
         do k = 1, min(nseen,8)
            if (Wface(1) .eq. seen(k)) known = .true.
         enddo
         if (.not. known) then
            nseen = nseen + 1
            if (nseen .le. 8) seen(nseen) = Wface(1)
         endif
      enddo
      write(*,'(A,I5,A,I4,A,2ES14.6)') '   [sweep] ', nsamp,               &
           ' samples gave ', nseen, ' distinct face densities: ',          &
           seen(1), seen(2)
      ndist = dble(nseen)
      end subroutine sweep_distinct_face_densities

      subroutine verdict(name, measured, reference, tol)
      character(len=*), intent(in) :: name
      real*8,           intent(in) :: measured, reference, tol
      if (abs(measured - reference) .le. tol) then
         write(*,'(A,A,A,ES13.6,A,ES13.6,A,ES9.2)') 'PASS ', name,         &
              ' measured=', measured, ' reference=', reference, ' tol=', tol
      else
         write(*,'(A,A,A,ES13.6,A,ES13.6,A,ES9.2)') 'FAIL ', name,         &
              ' measured=', measured, ' reference=', reference, ' tol=', tol
         n_fail = n_fail + 1
      endif
      end subroutine verdict

      end program base_branch_discriminant
