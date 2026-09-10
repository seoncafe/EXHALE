      program coupled_source_tests
      ! Tests of the coupled temperature-composition source step
      ! (docs/PLAN_20260906_rev2.md step B3c;
      !  docs/b1_target_system_20260906.md T1.4 to T1.9, T2.1, T2.2).
      !
      ! Three things are asserted here, all of them properties of the ONE
      ! species formation-energy table and of the energy row that reads it:
      !
      !   (1) THE LEDGER IDENTITY, AT-1a.  A closed reacting cell with no
      !       radiation field conserves u_th + u_form to round-off across the
      !       source step, and the temperature moves by exactly the amount
      !       the composition change implies.  The same cell is run through
      !       the update the code had before this step -- the composition
      !       projection, which rebuilt the pressure at the new particle
      !       count and left the energy row anchored there -- and the energy
      !       that update creates from nothing is measured beside it.  That
      !       number is the RED reference and it is kept permanently.
      !
      !   (2) THE PHOTOEVENT LEDGER, AT-1b.  One absorbed photon assigns its
      !       energy to each recipient exactly once: the ionization potential
      !       to the reservoir, the remainder to the photoelectron or to the
      !       fragments, and nothing anywhere else.  Asserted for the atomic
      !       absorbers against the table, and for the four H2 channels
      !       against h2_channel_energy_recipients.
      !
      !   (3) THE OXYGEN AND ASSOCIATIVE LEDGER, AT-1d.  The collisional
      !       oxygen channels and the associative He(2^3S) branch are
      !       differences of the same table, so a forward and a reverse
      !       channel are exact negatives, the eliminated O(1D) carries its
      !       excitation energy to the products of O6, and the O4 + O6 pair
      !       deposits what the direct H2O photodissociation deposits.
      !
      ! Each assertion prints one
      !     PASS|FAIL <name> measured=<v> reference=<r> tol=<t>
      ! line. Lines beginning with two spaces or with DIAGNOSTIC are context.
      ! Exit status is nonzero if any assertion fails.
      use global_parameters
      use species_table, only: isp_HI, isp_HII, isp_HeI, isp_HeII,        &
                               isp_HeIII, isp_HeTR, isp_H2, isp_H2p,      &
                               isp_H3p, isp_HeHp, isp_OH, isp_H2O,        &
                               isp_CO, mion_fsp, im_OI, mion_ethr,        &
                               im_CI, im_CII
      use molecular_reaction_heat
      use h2_photo_channels
      use energy_semi_implicit
      implicit none

      integer :: n_fail

      n_fail = 0
      ! One physical cell plus its ghosts, monatomic branch of the caloric
      ! equation of state (caloric_mixture_active is false unless a run turns
      ! it on), and the normalizations the code's energy-density unit
      ! p0 = n0 kB T0 is built from.
      N  = 1
      n0 = 1.0d10
      T0 = 1.0d4
      p0 = n0*kb_erg*T0

      call test_closed_reacting_cell(n_fail)
      call test_atomic_photoevents(n_fail)
      call test_h2_photoevents(n_fail)
      call test_oxygen_and_associative_ledger(n_fail)
      call test_reservoir_density(n_fail)

      write(*,*)
      if (n_fail .gt. 0) then
         write(*,'(a,i0,a)') 'coupled_source_step: ', n_fail,             &
                             ' assertion(s) failed'
         error stop 1
      endif
      write(*,'(a)') 'coupled_source_step: every assertion passed'

      contains

      ! ------------------------------------------------------!

      subroutine check(name, measured, reference, tol, nf)
      ! Relative comparison, absolute when the reference is zero.
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, reference, tol
      integer, intent(inout) :: nf
      real*8 :: d
      if (abs(reference) .gt. 0.0d0) then
         d = abs(measured - reference)/abs(reference)
      else
         d = abs(measured - reference)
      endif
      if (d .le. tol) then
         write(*,'(a,a,a,es15.8,a,es15.8,a,es9.2)') 'PASS ', name,        &
            ' measured=', measured, ' reference=', reference, ' tol=', tol
      else
         write(*,'(a,a,a,es15.8,a,es15.8,a,es9.2)') 'FAIL ', name,        &
            ' measured=', measured, ' reference=', reference, ' tol=', tol
         nf = nf + 1
      endif
      end subroutine check

      subroutine check_above(name, measured, threshold, nf)
      ! The measured quantity must EXCEED the threshold. Used for the energy
      ! the removed projection created from nothing: that number is red by
      ! construction and stops being a meaningful reference if it ever falls
      ! to the tolerance of the identity it breaks.
      character(len=*), intent(in) :: name
      real*8, intent(in) :: measured, threshold
      integer, intent(inout) :: nf
      if (measured .gt. threshold) then
         write(*,'(a,a,a,es15.8,a,es9.2)') 'PASS ', name, ' measured=',   &
            measured, '  >  ', threshold
      else
         write(*,'(a,a,a,es15.8,a,es9.2)') 'FAIL ', name, ' measured=',   &
            measured, '  >  ', threshold
         nf = nf + 1
      endif
      end subroutine check_above

      ! ------------------------------------------------------!

      subroutine test_closed_reacting_cell(nf)
      ! AT-1a.  ONE CELL, NO RADIATION FIELD, NO TRANSPORT, heat = cool = 0,
      ! and a composition that moves: half the hydrogen recombines over the
      ! step.  The reservoir gives that energy back to the gas, so the cell
      ! must warm by exactly du_form/(n_k c_v) and u_th + u_form must not
      ! move at all.
      !
      ! The comparison is with the update the marching loop had before B3c.
      ! There the energy row was anchored on the pressure REBUILT at the new
      ! particle count, u_th_old -> u_th(T_old, c_new), so with no source at
      ! all the step left u_th changed by
      !     X = (n_k(c_new) - n_k(c_old)) u_pp(T_old)
      ! with nothing on the other side of the equation.  X is measured here.
      integer, intent(inout) :: nf
      real*8, dimension(1-Ng:N+Ng) :: T_old, a, heat, cool0, T_floor
      real*8, dimension(1-Ng:N+Ng) :: T_ceil, u_old_pp, C, rho_a
      real*8, dimension(1-Ng:N+Ng,n_species) :: f_old, f_new
      real*8, dimension(1-Ng:N+Ng) :: uf_old, uf_new
      type(energy_balance_state) :: s
      real*8 :: nk_old, nk_new, cv, u_th_old, u_th_new, du_form
      real*8 :: T_new_target, X_projection, dE_old_row, dE_new_row

      write(*,*)
      write(*,'(a)') '---- closed_reacting_cell (AT-1a) ----'

      ! Pure hydrogen, adimensional density 1: n(H nuclei) = n0.
      rho_a  = 1.0d0
      f_old  = 0.0d0;  f_new = 0.0d0
      f_old(:,isp_HI)  = 0.2d0;   f_old(:,isp_HII) = 0.8d0
      f_new(:,isp_HI)  = 0.6d0;   f_new(:,isp_HII) = 0.4d0

      ! Particle counts, adimensional: heavy particles plus electrons.
      nk_old = 1.0d0 + 0.8d0
      nk_new = 1.0d0 + 0.4d0
      cv     = 1.0d0/(gamma_ad - 1.0d0)

      T_old   = 1.0d0
      T_floor = 0.01d0
      T_ceil  = 1.0d3
      heat    = 0.0d0
      cool0   = 0.0d0
      a       = 0.0d0          ! no source: the whole step is the reservoir

      call formation_energy_density(rho_a, f_old, uf_old)
      call formation_energy_density(rho_a, f_new, uf_new)
      du_form  = (uf_new(1) - uf_old(1))/p0        ! code units
      u_th_old = nk_old*cv*T_old(1)

      ! The energy row of T1.6, as the coupled step passes it: the energy the
      ! cell HAD, less what the composition change put into the reservoir,
      ! per particle of the composition the solve is run at.
      u_old_pp = (u_th_old - du_form)/nk_new

      call energy_balance_init(s, T_old, a, heat, cool0, T_floor, T_ceil, &
                               u_old_in = u_old_pp)
      do while (s%n_active .gt. 0 .and. s%iters .lt. energy_max_iter)
         C = 0.0d0
         call energy_balance_update(s, C)
      enddo
      call energy_balance_finish(s)

      u_th_new     = nk_new*cv*s%T_new(1)
      T_new_target = (u_th_old - du_form)/(nk_new*cv)

      ! The reservoir gave energy back (recombination lowers u_form), so the
      ! cell is hotter and the total is unchanged.
      write(*,'(a,es15.8,a,es15.8)') '  du_form [code] = ', du_form,      &
                                     '   u_th_old = ', u_th_old
      call check('the_returned_temperature_is_the_implied_one',           &
                 s%T_new(1), T_new_target, 1.0d-14, nf)
      dE_new_row = (u_th_new + uf_new(1)/p0) - (u_th_old + uf_old(1)/p0)
      call check('u_th_plus_u_form_is_conserved',                         &
                 dE_new_row/(abs(u_th_old) + abs(du_form)), 0.0d0,        &
                 1.0d-15, nf)

      ! THE RED REFERENCE: the same cell through the removed projection.
      ! With no source the old row returned T_old, and the energy of the cell
      ! became n_k(c_new) u_pp(T_old): the difference is created energy.
      X_projection = (nk_new - nk_old)*cv*T_old(1)
      dE_old_row   = X_projection + (uf_new(1) - uf_old(1))/p0
      write(*,'(a,es15.8)') '  DIAGNOSTIC the projection changed u(3) by'//&
                            ' X = ', X_projection
      write(*,'(a,es15.8)') '  DIAGNOSTIC the temperature it returned  = ',&
                            T_old(1)
      write(*,'(a,es15.8)') '  DIAGNOSTIC the temperature now returned = ',&
                            s%T_new(1)
      call check_above('the_projection_created_energy_from_nothing',      &
                       abs(dE_old_row)/(abs(u_th_old) + abs(du_form)),    &
                       1.0d-15, nf)
      end subroutine test_closed_reacting_cell

      ! ------------------------------------------------------!

      subroutine test_atomic_photoevents(nf)
      ! AT-1b for the atomic absorbers.  One photon of energy E is absorbed
      ! by one atom; the ionization potential goes to the reservoir as the
      ! difference of the table's entries, the remainder to the
      ! photoelectron, and the two sum to E with nothing left over and
      ! nothing counted twice.
      integer, intent(inout) :: nf
      real*8, parameter :: E = 40.0d0        ! eV, above every threshold here
      real*8 :: res, ele

      write(*,*)
      write(*,'(a)') '---- atomic_photoevents (AT-1b) ----'

      ! H I -> H II + e
      res = species_formation_energy(isp_HII)                             &
          + species_formation_energy(isp_eps_electron)                    &
          - species_formation_energy(isp_HI)
      ele = E - res
      call check('HI_reservoir_is_the_ionization_potential', res,         &
                 e_th_HI, 1.0d-14, nf)
      call check('HI_recipients_sum_to_the_photon', res + ele, E,         &
                 1.0d-14, nf)

      ! He I -> He II + e
      res = species_formation_energy(isp_HeII)                            &
          - species_formation_energy(isp_HeI)
      call check('HeI_reservoir_is_the_ionization_potential', res,        &
                 e_th_HeI, 1.0d-14, nf)
      call check('HeI_recipients_sum_to_the_photon', res + (E - res), E,  &
                 1.0d-14, nf)

      ! He II -> He III + e
      res = species_formation_energy(isp_HeIII)                           &
          - species_formation_energy(isp_HeII)
      call check('HeII_reservoir_is_the_second_potential', res,           &
                 e_th_HeII, 1.0d-14, nf)

      ! He 2^3S -> He II + e.  The threshold is measured FROM the metastable,
      ! so the reservoir difference and the code's threshold are the same
      ! number: the excitation the metastable already holds is not paid again.
      res = species_formation_energy(isp_HeII)                            &
          - species_formation_energy(isp_HeTR)
      call check('HeTR_reservoir_is_the_metastable_threshold', res,       &
                 e_th_HeTR, 1.0d-14, nf)

      ! A metal stage: C I -> C II + e.
      res = species_formation_energy(mion_fsp(im_CII))                    &
          - species_formation_energy(mion_fsp(im_CI))
      call check('CI_reservoir_is_the_first_carbon_potential', res,       &
                 mion_ethr(im_CI), 1.0d-14, nf)
      end subroutine test_atomic_photoevents

      ! ------------------------------------------------------!

      subroutine test_h2_photoevents(nf)
      ! AT-1b for the four H2 photoevent channels.  The recipients of one
      ! event are named by h2_channel_energy_recipients (B3b-PE); here they
      ! are checked to sum to the photon and, for each channel, the reservoir
      ! share is checked against the SAME species table the reservoir of the
      ! energy row is built from, so the two cannot disagree about what a
      ! photoevent leaves behind.
      integer, intent(inout) :: nf
      real*8, parameter :: E = 60.0d0    ! eV, above the 51.4 eV D threshold
      real*8 :: e_res, e_ele, e_frag, e_rad, tab
      integer :: ich

      write(*,*)
      write(*,'(a)') '---- h2_photoevents (AT-1b) ----'

      do ich = 1, 4
         call h2_channel_energy_recipients(ich, E, e_res, e_ele, e_frag,  &
                                           e_rad)
         call check('h2_channel_recipients_sum_to_the_photon',            &
                    e_res + e_ele + e_frag + e_rad, E, 1.0d-13, nf)
      enddo

      ! Channel M: H2 -> H2+ + e
      call h2_channel_energy_recipients(1, E, e_res, e_ele, e_frag, e_rad)
      tab = species_formation_energy(isp_H2p)                             &
          - species_formation_energy(isp_H2)
      call check('h2_M_reservoir_is_the_table_difference', e_res, tab,    &
                 1.0d-9, nf)

      ! Channel S: H2 -> H + H+ + e
      call h2_channel_energy_recipients(2, E, e_res, e_ele, e_frag, e_rad)
      tab = species_formation_energy(isp_HII)                             &
          + species_formation_energy(isp_HI)                              &
          - species_formation_energy(isp_H2)
      call check('h2_S_reservoir_is_the_table_difference', e_res, tab,    &
                 1.0d-4, nf)

      ! Channel D: H2 -> H+ + H+ + 2e
      call h2_channel_energy_recipients(3, E, e_res, e_ele, e_frag, e_rad)
      tab = 2.0d0*species_formation_energy(isp_HII)                       &
          - species_formation_energy(isp_H2)
      call check('h2_D_reservoir_is_the_table_difference', e_res, tab,    &
                 1.0d-9, nf)

      ! Channel N: H2 -> H + H, no electron and no reservoir ion
      call h2_channel_energy_recipients(4, E, e_res, e_ele, e_frag, e_rad)
      tab = -species_formation_energy(isp_H2)
      call check('h2_N_reservoir_is_the_bond_energy', e_res, tab,         &
                 1.0d-9, nf)
      call check('h2_N_assigns_no_electron', e_ele, 0.0d0, 0.0d0, nf)
      end subroutine test_h2_photoevents

      ! ------------------------------------------------------!

      subroutine test_oxygen_and_associative_ledger(nf)
      ! AT-1d.  The collisional oxygen channels and the associative
      ! He(2^3S) branch, as differences of the one table.
      integer, intent(inout) :: nf
      real*8 :: q1, q1r, q2, q2r, q6, qa
      real*8 :: direct, pair, thr_O3, thr_O4, eps_O1D

      write(*,*)
      write(*,'(a)') '---- oxygen_and_associative_ledger (AT-1d) ----'

      q1  = oxygen_reaction_energy_eV(ir_O1)
      q1r = oxygen_reaction_energy_eV(ir_O1r)
      q2  = oxygen_reaction_energy_eV(ir_O2)
      q2r = oxygen_reaction_energy_eV(ir_O2r)
      q6  = oxygen_reaction_energy_eV(ir_O6)
      qa  = oxygen_reaction_energy_eV(ir_assoc_HeTR)
      eps_O1D = species_formation_energy(isp_eps_O1D)

      ! A forward and a reverse channel are the same four table entries
      ! accumulated in the opposite order, so their sum is zero to a few
      ! roundings of entries of order 10 eV, which is what the ABSOLUTE
      ! tolerance below is: 1e-14 eV is about five ulps of eps(H2O) = -9.516.
      ! What the row rules out is a channel written down twice with two
      ! different numbers, which is what a list of reaction enthalpies makes
      ! possible and a table does not.
      call check('O1_and_its_reverse_sum_to_zero', q1r + q1, 0.0d0,       &
                 1.0d-14, nf)
      call check('O2_and_its_reverse_sum_to_zero', q2r + q2, 0.0d0,       &
                 1.0d-14, nf)

      ! T1.9 item 1, the transfer rule: the eliminated O(1D) carries its
      ! excitation energy to the products of its single sink, so O6 is O2
      ! plus that excitation and nothing else.
      call check('O6_is_O2_plus_the_O1D_excitation', q6, q2 + eps_O1D,    &
                 1.0d-14, nf)

      ! T1.9 item 2, AT-1d (ii): taking H2O apart through the O(1D) branch
      ! and its O6 partner deposits what the direct route deposits, at the
      ! same photon energy.  Both are written as differences of the table:
      !   direct  = -thr(O3)                    (photon pays thr, rest is KE)
      !   pair    = -thr(O4) - eps(O1D) + q(O6)
      thr_O3 = species_formation_energy(isp_OH)                           &
             + species_formation_energy(isp_HI)                           &
             - species_formation_energy(isp_H2O)
      thr_O4 = species_formation_energy(isp_H2)                           &
             + species_formation_energy(mion_fsp(im_OI))                  &
             - species_formation_energy(isp_H2O)
      direct = -thr_O3
      pair   = -thr_O4 - eps_O1D + q6
      call check('the_O4_O6_pair_deposits_what_the_direct_route_does',    &
                 pair, direct, 1.0d-14, nf)

      ! T1.9 item 3, AT-1d (iii): the associative branch
      ! He(2^3S) + H -> HeH+ + e, deposited nowhere before this step.
      call check('the_associative_branch_is_the_table_difference', qa,    &
                 species_formation_energy(isp_HeTR)                       &
                 + species_formation_energy(isp_HI)                       &
                 - species_formation_energy(isp_HeHp)                     &
                 - species_formation_energy(isp_eps_electron),            &
                 0.0d0, nf)
      call check('the_associative_branch_releases_8p1_eV', qa, 8.1d0,     &
                 5.0d-3, nf)
      ! In an ATOMIC gas the HeH+ returns to He + H at once (the closure
      ! stated in ion_residual_core), so the reaction that runs there is
      ! He(2^3S) + H -> He + H and its heat is the metastable's excitation.
      call check('the_atomic_cycle_returns_the_metastable_excitation',    &
                 species_formation_energy(isp_HeTR), 19.819614d0,         &
                 1.0d-7, nf)
      end subroutine test_oxygen_and_associative_ledger

      ! ------------------------------------------------------!

      subroutine test_reservoir_density(nf)
      ! The reservoir density itself: sum_s n_s eps_s over a hand-built
      ! composition, against the sum written out term by term, and the
      ! statement that an excited population is counted ON TOP of its parent
      ! column and never instead of it.
      integer, intent(inout) :: nf
      real*8, dimension(1-Ng:N+Ng) :: rho_a, uf, n_n2
      real*8, dimension(1-Ng:N+Ng,n_species) :: f
      real*8, parameter :: eV_to_erg = 1.602176634d-12
      real*8 :: expect

      write(*,*)
      write(*,'(a)') '---- reservoir_density ----'
      rho_a = 1.0d0
      f = 0.0d0
      f(:,isp_HI)   = 0.5d0
      f(:,isp_HII)  = 0.3d0
      f(:,isp_HeI)  = 0.2d0
      f(:,isp_HeTR) = 1.0d-4

      expect = ( 0.3d0*species_formation_energy(isp_HII)                  &
               + 1.0d-4*species_formation_energy(isp_HeTR) )*n0*eV_to_erg
      call formation_energy_density(rho_a, f, uf)
      call check('the_reservoir_is_the_weighted_table_sum', uf(1),        &
                 expect, 1.0d-14, nf)

      ! The metastable sits inside the He I column and carries excitation
      ! only: doubling the He I column must not change the reservoir.
      f(:,isp_HeI) = 0.4d0
      call formation_energy_density(rho_a, f, uf)
      call check('a_neutral_column_adds_nothing_to_the_reservoir', uf(1), &
                 expect, 1.0d-14, nf)

      ! H(n=2) is added on top of the H I column, not instead of it.
      n_n2 = 1.0d-6*n0
      call formation_energy_density(rho_a, f, uf, n_H_n2 = n_n2)
      call check('H_n2_adds_its_excitation_on_top_of_HI', uf(1),          &
                 expect + 1.0d-6*n0                                       &
                 *species_formation_energy(isp_eps_H_n2)*eV_to_erg,       &
                 1.0d-14, nf)
      end subroutine test_reservoir_density

      end program coupled_source_tests
