      program base_closure_candidates_probe
      ! TWO DISCRIMINANTS FOR THE SAME CHARACTERISTIC BASE CONDITION, READ ON
      ! ONE STATE UNDER THE PRODUCTION STATIONARY OPERATOR (item L35 step 2 of
      ! docs/PLAN_20260917.md section 9).
      !
      ! THE CHARACTERISTIC ANALYSIS, with the sign conventions written out.
      ! The radial coordinate increases outward, so at the lower face
      ! r_edg(0) the three eigenvalues of the one-dimensional Euler system are
      !
      !     lambda_1 = v - c ,   lambda_2 = v ,   lambda_3 = v + c ,
      !
      ! and a wave ENTERS the domain through that face when its eigenvalue is
      ! positive.  Gas entering at the base therefore has v > 0.  The four
      ! regimes and the number of conditions the reservoir owes are
      !
      !     0 < v < c    subsonic inflow      2 reservoir, 1 interior
      !    -c < v < 0    flow reversal        1 reservoir, 2 interior
      !         v > c    supersonic inflow    3 reservoir, 0 interior
      !         v < -c   supersonic outflow   0 reservoir, 3 interior
      !
      ! For the SUBSONIC INFLOW branch, which is the branch every certified
      ! LHS 1140 b state carries, the single outgoing wave is lambda_1 < 0 and
      ! its compatibility condition is the linearized left-running acoustic
      ! relation
      !
      !     p_b - rho_i c_i v_b  =  p_i - rho_i c_i v_i ,                 (C-)
      !
      ! the LODI form of Thompson (1987, J. Comput. Phys. 68, 1) and Poinsot
      ! and Lele (1992, J. Comput. Phys. 101, 104), in the statement of
      ! Carlson (2011, NASA/TM-2011-217181, section 2).  The two entering
      ! waves carry the reservoir pressure and the reservoir specific
      ! entropy, continued from the stated base level to the face along the
      ! hydrostatic isentrope.  Those are the sources base_boundary.f90 uses
      ! and this program calls; it re-derives none of them.
      !
      ! WHAT IS ACTUALLY IN QUESTION.  The three branches above are told
      ! apart by the SIGN of the face velocity, and at a base Mach number of
      ! 1e-6 that sign is not resolved by the cell-centered state.  So the
      ! production closure does not read it there: it reads a mass flux where
      ! the cell-centered product IS a flux, over the remote wind window
      ! r >= r_flux, and maps it onto the face through the interior density,
      !
      !     M_wind = F_wind / (rho_i r_b^2) / c_i ,
      !
      ! with F_wind the window mean of rho v r^2.  That is a stated NONLOCAL
      ! closure.  The local alternative is the same relation read at the face
      ! itself, from the mass flux the PRODUCTION Riemann solve puts through
      ! it under the stationary operator,
      !
      !     M_face = [ r_edg(0)^2 (rho v)_0 ] / (rho_i r_b^2) / c_i ,
      !
      ! which is what stationary_face_mass_flux returns as budget%flux_base.
      ! The evidence that excluded a local closure (L26 memo section R3) was
      ! taken from PLM evaluations of WENO3 states and from cell-centered
      ! extrapolations, and both were withdrawn
      ! (docs/session_handoff_20260917_rev1.md section 4.3).
      !
      ! HOW THE TWO ARE COMPARED, so that only the discriminant differs.
      ! characteristic_base_face_state takes the flux, its relative spread
      ! and an availability flag as ARGUMENTS.  This program calls that one
      ! production routine twice on the same installed state:
      !
      !   remote : (F_wind, d_window, have_F) of wind_window_mass_flux,
      !            which is what base_boundary_states installs;
      !   local  : (budget%flux_base, 0, .true.), the production face flux,
      !            with zero spread because one face carries one flux.
      !
      ! Everything else -- the reservoir isentrope, (C-), the isentropic
      ! density at p_b, the supersonic-outflow blend, the face Mach cap and
      ! the eight-point ghost quadrature -- is the production code in both
      ! rows.  The difference printed is therefore the difference between the
      ! two discriminants and nothing else.
      !
      ! THE LOCAL ROW IS NOT A CLOSED CONDITION UNTIL IT REPRODUCES THE FLUX
      ! IT READ, so the program also takes Picard iterations on the local
      ! row: install the ghost the local row implies, re-evaluate the
      ! production face flux on that state, close again.  Whether that map
      ! has a fixed point, and where, is a measurement and not an assumption.
      !
      ! WHAT IS SWEPT.  For the retained remote closure, the window location
      ! r_flux, since the closure's answer is a function of where the window
      ! starts; a domain extension is a separate run directory and is not a
      ! sweep of this program.
      !
      ! NO BRANCH IS SELECTED HERE.  The deliverable is the comparison.
      !
      ! The state is whatever the run directory this program is started in
      ! holds: it calls input_read and init exactly as EXHALE_main does, so
      ! the grid, the reservoir, the composition of cell 1 and every boundary
      ! option are the ones that case runs with.  A second argument selects a
      ! perturbation of the loaded state, applied to the conserved array
      ! before anything is measured:
      !
      !   none      the state as it was loaded
      !   weak      the velocity of physical cells 1..6 scaled by 1e-2,
      !             which leaves a base carrying a hundredth of its own mass
      !             flux and does not touch the wind
      !   reversed  the velocity of physical cells 1..6 multiplied by -1
      !   window    the mass density of every cell at r >= r_flux scaled by
      !             1 + 1e-3, the momentum and the total energy carried with
      !             it at fixed velocity and fixed temperature, so that only
      !             the remote window moves and the base cells do not
      !   rest      the momentum of every cell set to zero at fixed mass
      !             density and fixed thermal energy, which is the state at
      !             rest the well-balanced property is stated on
      use global_parameters
      use Read_input,     only: input_read
      use Initialization, only: init
      use lya_rt,         only: lya_rt_allocate_arrays
      use excited_hydrogen, only: excited_H_allocate_arrays
      use ionization_equilibrium, only: ioniz_eq_allocate_arrays
      use mol_rates,      only: h2_thermochemistry_init
      use Conversion,     only: U_to_W, W_to_U
      use composition,    only: get_species_densities
      use species_table,  only: n_mion
      use BC_Apply,       only: Apply_BC, base_face_W, base_ghost_W,       &
                                base_face_lower_W
      use stationary_operator, only: stationary_face_mass_flux,            &
                                face_mass_flux_budget,                     &
                                select_stationary_reconstruction,          &
                                restore_reconstruction_selection,          &
                                reconstruction_selection
      use Reconstruction_step, only: Reconstruct
      use RK_integration, only: RK_rhs, face_flux
      use base_boundary,  only: read_base_branch_options,                  &
                                wind_window_mass_flux,                     &
                                characteristic_base_face_state,            &
                                base_ghost_averages,                       &
                                base_face_mach_blend,                      &
                                base_wind_window_spread,                   &
                                base_face_Mi_last,                         &
                                base_face_Mwind_last,                      &
                                base_face_swind_last,                      &
                                base_face_blend_last,                      &
                                base_face_rho_res_last,                    &
                                base_face_rho_rev_last
      implicit none

      real*8, allocatable :: W0(:,:), u(:,:), f_sp(:,:), Wq(:,:)
      real*8, allocatable :: u_state(:,:)
      real*8, allocatable :: rho_s(:), nhi_s(:), nhii_s(:), nhei_s(:)
      real*8, allocatable :: nheii_s(:), nheiii_s(:), nheiTR_s(:)
      real*8, allocatable :: ne_s(:), ntot_s(:), nm_s(:,:)
      real*8, allocatable :: r2rv(:)
      real*8, allocatable :: Wghost_r(:,:), Wghost_l(:,:)
      real*8, allocatable :: WLp(:,:), WRp(:,:), dFp(:,:), Sp(:,:)
      type(face_mass_flux_budget) :: budget_p
      real*8  :: Wface_r(3), Wface_l(3), Wi_r(3), Wi_l(3)
      real*8  :: Wlow_r(3), Wlow_l(3)
      real*8  :: F_win, d_win, F_face
      real*8  :: nhat1, T1
      real*8  :: Mi_r, Mw_r, sw_r, wrev_r, rres_r, rrev_r
      real*8  :: Mi_l, Mf_l, sw_l, wrev_l, rres_l, rrev_l
      real*8  :: r_flux_entry, rf
      logical :: have_F
      integer :: j, it, k
      character(len=64) :: label, perturb
      ! The window locations swept below, in R_p: two inside the region the
      ! window deliberately excludes, the production value 1.20, and four
      ! above it.
      real*8, parameter :: rf_sweep(8) = (/ 1.01d0, 1.05d0, 1.10d0,       &
           1.20d0, 1.50d0, 2.00d0, 5.00d0, 10.0d0 /)

      label   = 'state'
      perturb = 'none'
      if (command_argument_count() .ge. 1) call get_command_argument(1,label)
      if (command_argument_count() .ge. 2) call get_command_argument(2,perturb)

      call input_read
      open(unit = outfile, file = 'EXHALE_setup.out')
      call lya_rt_allocate_arrays
      call excited_H_allocate_arrays
      call ioniz_eq_allocate_arrays
      call h2_thermochemistry_init

      allocate(W0(3,1-Ng:N+Ng), u(3,1-Ng:N+Ng), Wq(3,1-Ng:N+Ng))
      allocate(u_state(3,1-Ng:N+Ng))
      allocate(f_sp(1-Ng:N+Ng,n_species))
      allocate(rho_s(1-Ng:N+Ng), nhi_s(1-Ng:N+Ng), nhii_s(1-Ng:N+Ng))
      allocate(nhei_s(1-Ng:N+Ng), nheii_s(1-Ng:N+Ng), nheiii_s(1-Ng:N+Ng))
      allocate(nheiTR_s(1-Ng:N+Ng), ne_s(1-Ng:N+Ng), ntot_s(1-Ng:N+Ng))
      allocate(nm_s(1-Ng:N+Ng,n_mion))
      allocate(r2rv(1-Ng:N+Ng))
      allocate(Wghost_r(3,1-Ng:0), Wghost_l(3,1-Ng:0))
      allocate(WLp(3,1-Ng:N+Ng), WRp(3,1-Ng:N+Ng),                        &
               dFp(3,1-Ng:N+Ng), Sp(3,1-Ng:N+Ng))

      call init(W0,u,f_sp)
      u_state = u
      call read_base_branch_options()
      r_flux_entry = r_flux

      call perturb_state(u_state, f_sp, perturb)

      write(*,'(A)') ''
      write(*,'(A)') '=========== base closure candidates probe ==========='
      write(*,'(A,A)') ' state label        : ', trim(label)
      write(*,'(A,A)') ' perturbation       : ', trim(perturb)
      write(*,'(A,I6,A,I6)') ' grid cells N =', N, '   ghosts Ng =', Ng
      write(*,'(A,ES16.9,A,ES16.9)') ' r_edg(0) =', r_edg(0),             &
           '   r(1) =', r(1)
      write(*,'(A,ES12.5,A,ES12.5,A,ES12.5)') ' n0 [cm^-3] =', n0,        &
           '  v0 [cm/s] =', v0, '  T0 [K] =', T0
      write(*,'(A,ES12.5,A,ES12.5)') ' base_face_mach_blend =',           &
           base_face_mach_blend, '  base_wind_window_spread =',           &
           base_wind_window_spread
      write(*,'(A)') ''

      call close_both(u_state, f_sp, .true.)

      ! ---------------------------------------------------------------- !
      ! the Picard map of the local row: does it reproduce the flux it read
      ! ---------------------------------------------------------------- !
      write(*,'(A)') ''
      write(*,'(A)') '---- the local row as a map: the ghost it implies '// &
                     'installed, the face flux re-evaluated on that '//     &
                     'state ----'
      write(*,'(A)') '  The first row reads the face flux of the state as'//&
                     ' the production boundary leaves it; every later row'
      write(*,'(A)') '  reads it after the LOCAL ghost has replaced that'// &
                     ' boundary, so a fixed point of this map is a'
      write(*,'(A)') '  closure that reproduces the flux it was closed on.'
      write(*,'(A)') '  it   F_face [code]      F_face/F_window     '//     &
                     'rho_b [mH/cm3]     v_b [cm/s]        w_rev'
      u  = u_state
      call stationary_face_mass_flux(u, f_sp, budget_p, r2rv)
      F_face = budget_p%flux_base
      do it = 1, 6
         call cell1_thermo(u, f_sp, nhat1, T1, Wq)
         call characteristic_base_face_state(Wq(:,1), nhat1, T1, F_face,   &
                       0.0d0, .true., Wface_l, Wi_l)
         wrev_l = base_face_blend_last
         call base_ghost_averages(Wface_l, Wghost_l, Wlow_l)
         if (F_win .ne. 0.0d0) then
            write(*,'(A,I3,5ES19.10)') '  ROW-PICARD ', it, F_face,        &
                 F_face/F_win, Wface_l(1)*n0, Wface_l(2)*v0, wrev_l
         else
            write(*,'(A,I3,ES19.10,A,3ES19.10)') '  ROW-PICARD ', it,      &
                 F_face, '      (no window)', Wface_l(1)*n0,               &
                 Wface_l(2)*v0, wrev_l
         endif
         F_face = base_face_flux_of_installed(u, f_sp, Wface_l, Wghost_l,  &
                                              Wlow_l)
      enddo

      ! ---------------------------------------------------------------- !
      ! the retained remote closure against the window location
      ! ---------------------------------------------------------------- !
      write(*,'(A)') ''
      write(*,'(A)') '---- the remote closure against the window inner '//  &
                     'radius r_flux ----'
      write(*,'(A)') '   r_flux   j_flux  cells   F_window        '//       &
                     'd_window     M_wind       s_wind    w_rev    '//      &
                     '  rho_b            v_b [cm/s]'
      do k = 1, 8
         rf = rf_sweep(k)
         call set_flux_window(rf)
         call remote_row(u_state, f_sp, rf)
      enddo
      call set_flux_window(r_flux_entry)

      write(*,'(A)') ''
      write(*,'(A)') '=========== end of probe ==========='

      contains

      ! ------------------------------------------------------------------ !

      subroutine set_flux_window(rf)
      ! The inner edge of the wind window, the way define_grid sets it: the
      ! first PHYSICAL cell whose center is at or above rf.
      real*8, intent(in) :: rf
      integer :: jj
      r_flux = rf
      j_flux = N
      do jj = 1, N
         if (r(jj) .ge. rf) then
            j_flux = jj
            exit
         endif
      enddo
      j_flux = max(j_flux, 1)
      end subroutine set_flux_window

      ! ------------------------------------------------------------------ !

      subroutine perturb_state(uu, fsp, what)
      ! The perturbations named in the header, applied to the conserved
      ! array.  The base ones move the momentum of physical cells 1 to 6 at
      ! fixed mass density and fixed THERMAL energy, so only the velocity and
      ! the kinetic energy move; the window one moves the mass density of the
      ! wind at fixed velocity and fixed temperature, so the base cells are
      ! untouched by construction.
      real*8, intent(inout) :: uu(3,1-Ng:N+Ng)
      real*8, intent(in)    :: fsp(1-Ng:N+Ng,n_species)
      character(len=*), intent(in) :: what
      real*8 :: vv, ek, s
      integer :: jj
      if (trim(what) .eq. 'none') return
      if (trim(what) .eq. 'weak' .or. trim(what) .eq. 'reversed') then
         s = 1.0d-2
         if (trim(what) .eq. 'reversed') s = -1.0d0
         do jj = 1, 6
            vv = uu(2,jj)/uu(1,jj)
            ek = 0.5d0*uu(1,jj)*vv*vv
            uu(3,jj) = uu(3,jj) - ek
            uu(2,jj) = uu(1,jj)*(s*vv)
            uu(3,jj) = uu(3,jj) + 0.5d0*uu(1,jj)*(s*vv)*(s*vv)
         enddo
      else if (trim(what) .eq. 'rest') then
         do jj = 1-Ng, N+Ng
            vv = uu(2,jj)/uu(1,jj)
            uu(3,jj) = uu(3,jj) - 0.5d0*uu(1,jj)*vv*vv
            uu(2,jj) = 0.0d0
         enddo
      else if (trim(what) .eq. 'window') then
         s = 1.0d0 + 1.0d-3
         do jj = max(j_flux,1), N
            uu(1,jj) = s*uu(1,jj)
            uu(2,jj) = s*uu(2,jj)
            uu(3,jj) = s*uu(3,jj)
         enddo
      else
         write(*,'(A,A)') ' (base_closure_candidates_probe) unknown '//    &
              'perturbation: ', trim(what)
         error stop 1
      endif
      end subroutine perturb_state

      ! ------------------------------------------------------------------ !

      subroutine cell1_thermo(uu, fsp, nhat, T1v, Wout)
      ! The primitive state and the two cell-1 quantities the base condition
      ! reads, derived the way base_boundary_states derives them: the
      ! composition first (it sets n_part_cell1), then p = n_part T.
      real*8, intent(in)  :: uu(3,1-Ng:N+Ng)
      real*8, intent(in)  :: fsp(1-Ng:N+Ng,n_species)
      real*8, intent(out) :: nhat, T1v
      real*8, intent(out) :: Wout(3,1-Ng:N+Ng)
      nhei_s = 0.0d0;  nheii_s = 0.0d0;  nheiii_s = 0.0d0;  nheiTR_s = 0.0d0
      rho_s = uu(1,:)
      call get_species_densities(rho_s, fsp, nhi_s, nhii_s, nhei_s,        &
                                 nheii_s, nheiii_s, nheiTR_s, nm_s,        &
                                 ne_s, ntot_s)
      call U_to_W(uu, Wout)
      nhat = n_part_cell1/Wout(1,1)
      T1v  = Wout(3,1)/n_part_cell1
      end subroutine cell1_thermo

      ! ------------------------------------------------------------------ !

      real*8 function base_face_flux_of_installed(uu, fsp, Wf, Wg, Wl)     &
                      result(Fb)
      ! r_edg(0)^2 (rho v)_0 of the state whose lower boundary is the GIVEN
      ! face state, ghost averages and lower face state, under the
      ! stationary operator.  Apply_BC is deliberately not called: it would
      ! install the production boundary in place of the one being tested.
      ! The outer ghosts stay as the caller's state carries them; they do
      ! not enter face 0.
      real*8, intent(in) :: uu(3,1-Ng:N+Ng)
      real*8, intent(in) :: fsp(1-Ng:N+Ng,n_species)
      real*8, intent(in) :: Wf(3), Wg(3,1-Ng:0), Wl(3)
      real*8 :: nh, Tc
      type(reconstruction_selection) :: entry_selection
      integer :: jj
      call cell1_thermo(uu, fsp, nh, Tc, Wq)
      do jj = 1-Ng, 0
         Wq(:,jj) = Wg(:,jj)
      enddo
      base_face_W       = Wf
      base_ghost_W      = Wg
      base_face_lower_W = Wl
      call W_to_U(Wq, u)
      call select_stationary_reconstruction(entry_selection)
      call Reconstruct(u, WLp, WRp)
      call RK_rhs(u, WLp, WRp, dFp, Sp)
      call restore_reconstruction_selection(entry_selection)
      Fb = face_flux(1,0)*r_edg(0)*r_edg(0)
      end function base_face_flux_of_installed

      ! ------------------------------------------------------------------ !

      subroutine close_both(uu, fsp, verbose)
      ! The two closures of the same installed state, side by side.
      real*8, intent(in) :: uu(3,1-Ng:N+Ng)
      real*8, intent(in) :: fsp(1-Ng:N+Ng,n_species)
      logical, intent(in):: verbose
      type(face_mass_flux_budget) :: budget
      integer :: jj

      call stationary_face_mass_flux(uu, fsp, budget, r2rv)
      F_face = budget%flux_base
      call cell1_thermo(uu, fsp, nhat1, T1, Wq)
      call wind_window_mass_flux(Wq, F_win, d_win, have_F)

      if (verbose) then
         write(*,'(A)') '---- the state, under the production stationary '//&
                        'operator ----'
         write(*,'(A,A,A,A,A,L1)') ' operator = ',                         &
              trim(budget%operator_name), '   flux = ',                    &
              trim(budget%flux_name), '   well balanced = ',               &
              budget%well_balanced
         write(*,'(A,ES19.10)') ' window mean  F_window   [code] =',       &
              budget%window_mean
         write(*,'(A,ES19.10)') ' base face    F_face     [code] =',       &
              budget%flux_base
         write(*,'(A,ES19.10,A,I6)') ' face minimum           [code] =',   &
              budget%flux_min, '   at face ', budget%j_min
         write(*,'(A,ES19.10,A,I6)') ' face maximum           [code] =',   &
              budget%flux_max, '   at face ', budget%j_max
         if (budget%window_mean .ne. 0.0d0) then
            write(*,'(A,ES19.10)') ' F_face / F_window             =',     &
                 budget%flux_base/budget%window_mean
            write(*,'(A,ES19.10)') ' (F_max-F_min)/|F_window|      =',     &
                 (budget%flux_max - budget%flux_min)                       &
                 /abs(budget%window_mean)
         endif
         write(*,'(A,ES19.10,A,L1)') ' window spread d_window        =',   &
              d_win, '   have_F = ', have_F
         write(*,'(A)') ''
      endif

      ! -- the remote-window closure, the one base_boundary_states installs
      call characteristic_base_face_state(Wq(:,1), nhat1, T1, F_win,       &
                    d_win, have_F, Wface_r, Wi_r)
      Mi_r = base_face_Mi_last;  Mw_r = base_face_Mwind_last
      sw_r = base_face_swind_last;  wrev_r = base_face_blend_last
      rres_r = base_face_rho_res_last;  rrev_r = base_face_rho_rev_last
      call base_ghost_averages(Wface_r, Wghost_r, Wlow_r)

      ! -- the local characteristic closure, on the production face flux
      call characteristic_base_face_state(Wq(:,1), nhat1, T1, F_face,      &
                    0.0d0, .true., Wface_l, Wi_l)
      Mi_l = base_face_Mi_last;  Mf_l = base_face_Mwind_last
      sw_l = base_face_swind_last;  wrev_l = base_face_blend_last
      rres_l = base_face_rho_res_last;  rrev_l = base_face_rho_rev_last
      call base_ghost_averages(Wface_l, Wghost_l, Wlow_l)

      write(*,'(A)') '---- the two closures of the same state ----'
      write(*,'(A)') ' quantity                 remote window        '//   &
                     'local face flux      (local-remote)/remote'
      call cmp(' face  rho  [mH/cm3] ', Wface_r(1)*n0,   Wface_l(1)*n0)
      call cmp(' face  v    [cm/s]   ', Wface_r(2)*v0,   Wface_l(2)*v0)
      call cmp(' face  p    [code]   ', Wface_r(3),      Wface_l(3))
      call cmp(' ghost0 rho [mH/cm3] ', Wghost_r(1,0)*n0, Wghost_l(1,0)*n0)
      call cmp(' ghost0 v   [cm/s]   ', Wghost_r(2,0)*v0, Wghost_l(2,0)*v0)
      call cmp(' ghost0 p   [code]   ', Wghost_r(3,0),   Wghost_l(3,0))
      call cmp(' ghost1 rho [mH/cm3] ', Wghost_r(1,-1)*n0, Wghost_l(1,-1)*n0)
      call cmp(' ghost1 v   [cm/s]   ', Wghost_r(2,-1)*v0, Wghost_l(2,-1)*v0)
      call cmp(' ghost1 p   [code]   ', Wghost_r(3,-1),  Wghost_l(3,-1))
      write(*,'(A)') ''
      write(*,'(A)') ' discriminants and weights'
      call cmp(' M_i                 ', Mi_r,   Mi_l)
      call cmp(' M of the branch     ', Mw_r,   Mf_l)
      call cmp(' s (window say)      ', sw_r,   sw_l)
      call cmp(' w_rev               ', wrev_r, wrev_l)
      call cmp(' rho_res [mH/cm3]    ', rres_r*n0, rres_l*n0)
      call cmp(' rho_rev [mH/cm3]    ', rrev_r*n0, rrev_l*n0)
      write(*,'(A,ES19.10,A,ES19.10)') ' interior at the face: rho =',     &
           Wi_r(1)*n0, '  v [cm/s] =', Wi_r(2)*v0
      ! the ghosts the boundary actually installed on this state, as a check
      ! that the remote row IS the installed one
      call U_to_W(uu, Wq)
      call W_to_U(Wq, u)
      call Apply_BC(u)
      call U_to_W(u, Wq)
      write(*,'(A,3ES19.10)') ' installed ghost0 (rho,v,p) [code] =',      &
           Wq(1,0), Wq(2,0), Wq(3,0)
      write(*,'(A,3ES19.10)') ' remote row ghost0          [code] =',      &
           Wghost_r(1,0), Wghost_r(2,0), Wghost_r(3,0)
      end subroutine close_both

      ! ------------------------------------------------------------------ !

      subroutine remote_row(uu, fsp, rf)
      ! The remote closure at one window location.
      real*8, intent(in) :: uu(3,1-Ng:N+Ng)
      real*8, intent(in) :: fsp(1-Ng:N+Ng,n_species)
      real*8, intent(in) :: rf
      real*8 :: Fw, dw, Wf(3), Wi(3)
      logical :: hF
      call cell1_thermo(uu, fsp, nhat1, T1, Wq)
      call wind_window_mass_flux(Wq, Fw, dw, hF)
      call characteristic_base_face_state(Wq(:,1), nhat1, T1, Fw, dw, hF,  &
                                          Wf, Wi)
      write(*,'(A,F8.2,I7,I7,7ES14.6)') '  ROW-RFLUX ', rf, j_flux,        &
           N - max(j_flux,1) + 1, Fw, dw, base_face_Mwind_last,            &
           base_face_swind_last, base_face_blend_last, Wf(1)*n0,           &
           Wf(2)*v0
      end subroutine remote_row

      ! ------------------------------------------------------------------ !

      subroutine cmp(name, a, b)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: a, b
      real*8 :: d
      d = 0.0d0
      if (a .ne. 0.0d0) d = (b - a)/a
      write(*,'(A,2ES21.12,ES15.5)') name, a, b, d
      end subroutine cmp

      end program base_closure_candidates_probe
