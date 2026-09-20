      program closure_returned_state
      ! WHAT THE THERMOCHEMICAL CLOSURE HANDS BACK, MEASURED ON THE COLUMN
      ! THE carrier_retry ROW READS (PLAN_20260918_rev2 item D4a).
      !
      ! Sections 1, 2, 3, 5 and 6 MEASURE, and assert only that their own
      ! diagnostics put the state back. Section 4 is a gate: a trial the
      ! relaxation refuses has to leave the state it was entered with, item
      ! by item. No tolerance of the production code is touched.
      !
      ! THE STATE. The twelve-cell molecular hydrogen and helium column of
      ! src/tests/carrier_retry/carrier_retry.f90, with the molecular
      ! chemistry and the transported proton on and no metal abundance, the
      ! same checkpoint (one transport interval taken from the entry
      ! composition), the same wind 20 r and the same relaxation pass at
      ! trust 1e-2. The reproduction is checked against the number the
      ! unmodified suite prints for its section-11 row.
      !
      ! WHAT IS MEASURED, in the order of the item:
      !  (1) the composition, temperature, electron density, imposed
      !      fractions, radiation field and rate state the production
      !      closure handed back, written to a file;
      !  (2) the chemical residual of that composition, row by row, both
      !      against the rate state the closure left (the temperature its
      !      last sweep ran at) and against a rate state rebuilt at the
      !      RETURNED temperature with the returned electron density;
      !  (3) three complete fixed-conserved-energy closure increments
      !      (chemistry AND the equation-of-state temperature update),
      !      taken by the production routine itself on an isolated copy;
      !  (4) what a refused trial puts back: composition, rates, caloric
      !      state, boundary context, carrier checkpoint;
      !  (5) the conditioning, as the constrained Jacobian of the
      !      normalized residual rows in the state vector and as a
      !      controlled perturbation of the transported parent;
      !  (6) the stopping rule against the residuals at the exit.
      !
      ! THE ISOLATION. Every diagnostic that advances the state does so on
      ! a copy, with the module state the next real operation reads saved
      ! and put back; the restore is PROVEN by re-evaluating the residual
      ! of (2) and requiring the same bits.

      use global_parameters
      use species_table
      use ion_cell_state, only: ion_rates
      use ionization_equilibrium, only: bg_cell, bg_ready,                &
           ioniz_eq_allocate_arrays, ioniz_eq, ieq_res_tol,               &
           ieq_rate_cell, ieq_ne_cell, ieq_TK_cell, ieq_ntot_cell,        &
           ieq_met_coef, ieq_rates_ready, ieq_neq_stored,                 &
           ieq_mbase_stored, ieq_iox_stored, ieq_nonroot_streak,          &
           ieq_state_vector, ionization_closure_residual_cell,            &
           ionization_closure_residual_profile
      use diffusive_photochemistry, only:                                 &
           carrier_set_init, carrier_transport_interval,                  &
           carrier_checkpoint, carrier_checkpoint_take,                   &
           carrier_checkpoint_restore, carrier_checkpoint_matches,        &
           carrier_verdict, n_carrier_max, ic_H2, ic_Hp,                  &
           carrier_transport_stop_on_failure,                             &
           relax_photochemical_composition,                               &
           carrier_relax_outcome_text, carrier_relax_chemistry_refused,   &
           equilibrate_chemistry_at_fixed_conserved_state,                &
           chem_cycles_cap_for_test, chem_cycles_max, chem_cycle_tol,     &
           n_chem_last_reason, chem_closure_reason_text,                  &
           chem_last_increment, chem_last_offsimplex
      use test_columns,  only: column_carrying_its_own_density
      use Conversion,     only: W_to_U
      use caloric_eos,    only: pressure_from_energy_density
      use element_census, only: n_element
      use gravity_grid_construction, only: set_gravity_grid
      use energy_vectors_construct,  only: set_energy_vectors
      use charge_exchange,           only: cx_init
      use composition,               only: get_species_densities,         &
                                           comp_p_from_T, comp_T_from_p
      use base_boundary,             only: set_base_reservoir
      use BC_Apply,                  only: Apply_BC
      use steady_residual_mod,       only: assemble_residual
      use assertion_report
      implicit none

      integer, parameter :: ncell = 12
      integer, parameter :: nxmax = 10 + 2*12

      real*8,  allocatable :: rho(:), v(:), dt_code(:)
      real*8,  allocatable :: f_sp(:,:), f_sp0(:,:), f_ret(:,:)
      real*8,  allocatable :: f_work(:,:), f_probe(:,:), f_hold(:,:)
      real*8,  allocatable :: fc_probe(:,:)
      real*8,  allocatable :: p_col(:), T_col(:)
      real*8,  allocatable :: heat_col(:), cool_col(:), eta_col(:)
      real*8,  allocatable :: u_col(:,:), W_col(:,:)
      real*8,  allocatable :: ntot0_c(:), T_ent(:)
      real*8,  allocatable :: v_relax(:)
      real*8,  allocatable :: nhi_r(:), nhii_r(:), nhei_r(:), nheii_r(:)
      real*8,  allocatable :: nheiii_r(:), nheiTR_r(:), ne_r(:), ntot_r(:)
      real*8,  allocatable :: nm_r(:,:), T_ret(:)
      real*8,  allocatable :: p_work(:), T_work(:)
      real*8,  allocatable :: heat_work(:), cool_work(:), eta_work(:)
      ! The module state the next real operation reads, held across a
      ! diagnostic and put back.
      type(ion_rates), allocatable :: bg_hold(:), rate_hold(:)
      real*8,  allocatable :: ne_hold(:), TK_hold(:), ntot_hold(:)
      real*8,  allocatable :: met_hold(:,:,:)
      integer, allocatable :: streak_hold(:)
      real*8,  allocatable :: p_cal0(:), p_cal1(:)
      real*8,  allocatable :: p_hold(:), T_hold(:)
      real*8,  allocatable :: heat_hold(:), cool_hold(:), eta_hold(:)
      logical :: ready_hold
      integer :: neq_hold, mbase_hold, iox_hold
      real*8  :: dpbc_hold, npc1_hold
      type(carrier_checkpoint) :: chk0, chk_iso

      real*8,  dimension(1:ncell) :: res_solved, res_trip_solved
      real*8,  dimension(1:ncell) :: res_retT,   res_trip_retT
      real*8,  dimension(1:ncell) :: res_back,   res_trip_back
      character(len=160) :: why

      real*8  :: xv(nxmax), rowres(nxmax), xp(nxmax), rowp(nxmax)
      real*8  :: rowm(nxmax), rescell
      integer :: nx, nrow

      real*8  :: move_rel, fref, mv, dT_inc
      real*8  :: dne(3), dTrel(3), dperp(3,3), eps3(3), eps5(5)
      real*8  :: rowsens(7,5), rowmaxs(7,5), slope(7), adm(7)
      real*8  :: chain(3), ne_unpert
      real*8  :: dTshift(4) = (/ 0.0d0, 1.0d-7, 1.0d-6, 1.0d-5 /)
      integer :: jmove, ispmove, jbind, j, k, isp, m, ii, kk
      integer :: n_cycles, reason, ns, oc, nsub
      real*8  :: drift, frac_done
      logical :: completed, okc
      type(carrier_verdict) :: verdict
      real*8  :: big
      integer :: jbig, ispbig
      real*8  :: hstep, x0
      character(len=8) :: rowname(7)

      call setup_globals()
      call build_molecular_hydrogen_column()

      ! ---- the checkpoint the section-11 pass is entered with --------- !
      allocate(fc_probe(1-Ng:N+Ng,n_carrier_max))
      f_sp = f_sp0
      call carrier_transport_interval(rho, v, f_sp, dt_code, completed,   &
                                      frac_done, nsub, verdict)
      call carrier_checkpoint_take(chk0)
      run_mode = run_mode_phys
      carrier_transport_stop_on_failure = .true.

      allocate(v_relax(1-Ng:N+Ng))
      v_relax = 20.0d0*r
      allocate(nhi_r(1-Ng:N+Ng), nhii_r(1-Ng:N+Ng))
      allocate(nhei_r(1-Ng:N+Ng), nheii_r(1-Ng:N+Ng))
      allocate(nheiii_r(1-Ng:N+Ng), nheiTR_r(1-Ng:N+Ng))
      allocate(ne_r(1-Ng:N+Ng), ntot_r(1-Ng:N+Ng), T_ret(1-Ng:N+Ng))
      allocate(nm_r(1-Ng:N+Ng,n_mion))
      allocate(f_ret(1-Ng:N+Ng,n_species), f_work(1-Ng:N+Ng,n_species))
      allocate(f_probe(1-Ng:N+Ng,n_species), f_hold(1-Ng:N+Ng,n_species))
      allocate(p_work(1-Ng:N+Ng), T_work(1-Ng:N+Ng))
      allocate(heat_work(1-Ng:N+Ng), cool_work(1-Ng:N+Ng))
      allocate(eta_work(1-Ng:N+Ng))
      allocate(p_cal0(1-Ng:N+Ng), p_cal1(1-Ng:N+Ng))
      allocate(p_hold(1-Ng:N+Ng), T_hold(1-Ng:N+Ng))
      allocate(heat_hold(1-Ng:N+Ng), cool_hold(1-Ng:N+Ng))
      allocate(eta_hold(1-Ng:N+Ng))

      ! ---- (1) the production return ---------------------------------- !
      call seed_mass_row_of_the_column(v_relax)
      call seed_background_at_the_entry_composition()
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call conserved_of(v_relax)
      call relax_photochemical_composition(u_col, v_relax, f_sp, p_col,   &
                                           T_col, heat_col, cool_col,     &
                                           eta_col, 1.0d-2, drift,        &
                                           ns, oc)
      f_ret = f_sp
      write(*,'(A,A,A,I0,A,A,A,ES12.5)') ' (D4a) the pass ends on ',      &
           trim(carrier_relax_outcome_text(oc)), ', kept ', ns,           &
           ' step(s); the closure of its last trial reports ',            &
           trim(chem_closure_reason_text(n_chem_last_reason)),            &
           ', last temperature increment ', chem_last_increment
      write(*,'(A,ES10.3,A,ES10.3,A,I0)') ' (D4a) tolerances READ from'// &
           ' source: ieq_res_tol ', ieq_res_tol, ', chem_cycle_tol ',     &
           chem_cycle_tol, ', chem_cycles_max ', chem_cycles_max

      ! The densities and the temperature OF the returned composition at
      ! the fixed conserved state, formed exactly as
      ! pressure_and_temperature_at_fixed_conserved_state forms them.
      call densities_of(f_ret)
      do j = 1-Ng, N+Ng
         p_cal0(j) = pressure_from_energy_density(j, u_col(1,j),          &
                        u_col(3,j) - 0.5d0*u_col(2,j)**2/u_col(1,j))
      enddo
      call comp_T_from_p(p_col, ntot_r, ne_r, T_ret)

      ! One further sweep on the returned composition, the quantity the
      ! carrier_retry row reads, so that this program is anchored to it.
      call save_module_state()
      f_work = f_ret
      call ioniz_eq(T_col, rho, f_work, heat_work, cool_work, eta_work)
      move_rel = 0.0d0;  jmove = 0;  ispmove = 0
      do j = 1, N
         do isp = 1, n_species
            fref = f_ret(j,isp)
            if (fref .le. 1.0d-20) cycle
            mv = abs(f_work(j,isp) - fref)/fref
            if (mv .gt. move_rel) then
               move_rel = mv;  jmove = j;  ispmove = isp
            endif
         enddo
      enddo
      call restore_module_state(f_ret)
      write(*,'(A,ES22.15,A,I0,A,I0)') ' (D4a) one further sweep moves'// &
           ' the returned composition by ', move_rel, ', f_sp column ',   &
           ispmove, ' at cell ', jmove
      jbind = max(jmove, 1)

      ! ---- (2) the chemical residual of the returned composition ------ !
      !
      ! (2a) against the rate state the closure LEFT: the rates its last
      ! sweep was run at, which is the temperature of the cycle before the
      ! one that ended it.
      call ionization_closure_residual_profile(rho, f_ret, res_solved,    &
                                               res_trip_solved, okc, why)
      if (.not. okc) write(*,'(A,A)') ' (D4a) residual (as solved)'//     &
           ' unavailable: ', trim(why)

      ! (2b) against a rate state rebuilt at the RETURNED temperature and
      ! the returned electron density: one sweep at T_ret on a COPY fills
      ! ieq_rate_cell and ieq_ne_cell from (T_ret, f_ret); the residual is
      ! then evaluated on f_ret itself.
      call save_module_state()
      f_probe = f_ret
      call ioniz_eq(T_ret, rho, f_probe, heat_work, cool_work, eta_work)
      call ionization_closure_residual_profile(rho, f_ret, res_retT,      &
                                               res_trip_retT, okc, why)
      if (.not. okc) write(*,'(A,A)') ' (D4a) residual (at T returned)'// &
           ' unavailable: ', trim(why)
      ! Row by row at the binding cell, at the returned temperature.
      call ieq_state_vector(jbind, rho, f_ret, xv, nx)
      call ionization_closure_residual_cell(jbind, xv(1:nx), nx, rescell, &
                                            rowres(1:nx))
      call restore_module_state(f_ret)

      ! THE PROOF OF THE RESTORE: the residual of (2a), re-evaluated after
      ! the two isolated sweeps above, has to come back bit for bit.
      call ionization_closure_residual_profile(rho, f_ret, res_back,      &
                                               res_trip_back, okc, why)
      call check_absolute('the_as_solved_residual_profile_is'//           &
           '_reinstated_bit_for_bit',                                     &
           maxval(abs(res_back - res_solved)), 0.0d0, 0.0d0)

      rowname(1) = 'H_II'
      rowname(2) = 'He_II'
      rowname(3) = 'He_III'
      rowname(4) = 'H2'
      rowname(5) = 'H2p'
      rowname(6) = 'H3p'
      rowname(7) = 'HeHp'

      write(*,'(A)') ' (D4a) TABLE 2: the chemical residual of the'//     &
           ' returned composition'
      write(*,'(A)') '   cell  res(stale rates)  res(at T returned)'//&
                     '  T_ieq_rate[K]   T_bg_cell[K]  T_returned[K]'//    &
                     '  |dT|/T vs ieq_rate'
      do j = 1, N
         write(*,'(I6,2ES19.9,3F15.6,ES14.4)') j, res_solved(j),          &
              res_retT(j), ieq_TK_cell(j), bg_cell(j)%T_K, T_ret(j)*T0,   &
              abs(T_ret(j)*T0 - ieq_TK_cell(j))/max(T_ret(j)*T0,1.0d-300)
      enddo
      write(*,'(A,I0,A,ES12.4,A,ES12.4)') ' (D4a) binding cell ', jbind,  &
           ': res(as solved) ', res_solved(jbind), ', res(at T'//         &
           ' returned) ', res_retT(jbind)
      write(*,'(A)') ' (D4a) TABLE 2b: the rows of the binding cell at'// &
           ' the returned temperature'
      do k = 1, min(nx,7)
         write(*,'(A,I3,1X,A8,2ES18.8)') '   row ', k, rowname(k),        &
              xv(k), rowres(k)
      enddo
      if (nx .gt. 7) write(*,'(A,ES14.6)') '   the largest metal row'//   &
           ' residual of that cell: ', maxval(rowres(8:nx))

      ! ---- (1, written out) ------------------------------------------- !
      call write_returned_state()

      ! ---- (3) three complete closure increments, isolated ------------ !
      write(*,'(A)') ' (D4a) TABLE 3: complete fixed-conserved-energy'//  &
           ' closure increments on an isolated copy'
      call save_module_state()
      f_work    = f_ret
      heat_work = heat_col
      cool_work = cool_col
      eta_work  = eta_col
      f_hold    = f_ret
      do m = 1, 3
         chem_cycles_cap_for_test = 1
         call equilibrate_chemistry_at_fixed_conserved_state(u_col,       &
                   f_work, p_work, T_work, heat_work, cool_work,          &
                   eta_work, okc, n_cycles, reason, dT_inc)
         chem_cycles_cap_for_test = -1
         dTrel(m) = dT_inc
         ! What the increment moved, species by species, relative, over
         ! the species present by the fixture's own 1e-20 rule.
         big = 0.0d0;  jbig = 0;  ispbig = 0
         do j = 1, N
            do isp = 1, n_species
               fref = f_hold(j,isp)
               if (fref .le. 1.0d-20) cycle
               mv = abs(f_work(j,isp) - fref)/fref
               if (mv .gt. big) then
                  big = mv;  jbig = j;  ispbig = isp
               endif
            enddo
         enddo
         write(*,'(A,I0,A,ES12.4,A,ES12.4,A,I0,A,I0)') '   increment ',   &
              m, ': thermal |dT|/T ', dT_inc, ', largest chemical'//      &
              ' movement ', big, ' in f_sp column ', ispbig, ' at cell ', &
              jbig
         if (m .eq. 1) then
            chain(1) = relmove(f_ret, f_work, jbind, isp_H2p)
            chain(2) = relmove(f_ret, f_work, jbind, isp_H3p)
            chain(3) = relmove(f_ret, f_work, jbind, isp_HeHp)
         endif
         write(*,'(A,4ES13.4)') '     at the binding cell, relative'//    &
              ' movement of H2+, H3+, HeH+ and H I: ',                    &
              relmove(f_hold, f_work, jbind, isp_H2p),                    &
              relmove(f_hold, f_work, jbind, isp_H3p),                    &
              relmove(f_hold, f_work, jbind, isp_HeHp),                   &
              relmove(f_hold, f_work, jbind, 1)
         f_hold = f_work
      enddo
      call restore_module_state(f_ret)
      ! The same proof again after the three increments.
      call ionization_closure_residual_profile(rho, f_ret, res_back,      &
                                               res_trip_back, okc, why)
      call check_absolute('the_as_solved_residual_profile_survives'//     &
           '_three_isolated_increments',                                  &
           maxval(abs(res_back - res_solved)), 0.0d0, 0.0d0)
      call check_absolute('and_the_returned_composition_is_untouched',    &
           maxval(abs(f_ret - f_sp)), 0.0d0, 0.0d0)

      ! ---- (5) conditioning ------------------------------------------- !
      !
      ! (5a) WHAT A ROW RESIDUAL OF ieq_res_tol ADMITS IN ABUNDANCE. The
      ! returned composition is a root of the network at the returned
      ! temperature, so the normalized residual of a row grows from zero
      ! when the unknown that row balances is displaced: |R_k| is read at
      ! five relative displacements of x_k, |R_k|/eps is the sensitivity of
      ! the row to its own unknown, and ieq_res_tol divided by that
      ! sensitivity is the RELATIVE abundance error the residual tolerance
      ! admits for that species. The code reports the ABSOLUTE value of the
      ! row residual, so the displacement is one-sided and no derivative is
      ! taken through the kink at the root.
      call save_module_state()
      f_probe = f_ret
      call ioniz_eq(T_ret, rho, f_probe, heat_work, cool_work, eta_work)
      call ieq_state_vector(jbind, rho, f_ret, xv, nx)
      nrow = min(nx, 7)
      eps5(1) = 1.0d-8
      eps5(2) = 1.0d-7
      eps5(3) = 1.0d-6
      eps5(4) = 1.0d-4
      eps5(5) = 1.0d-2
      write(*,'(A)') ' (D4a) TABLE 5a: the normalized residual of a row'//&
           ' when its own unknown is displaced, at the binding cell'
      write(*,'(A)') '   row         x         |R| at eps 1e-8    1e-7'// &
           '       1e-6       1e-4       1e-2      |R|/eps   '//          &
           ' relative abundance admitted at ieq_res_tol'
      do kk = 1, nrow
         do m = 1, 5
            xp(1:nx) = xv(1:nx)
            xp(kk) = xv(kk)*(1.0d0 + eps5(m))
            call ionization_closure_residual_cell(jbind, xp(1:nx), nx,    &
                                                  rescell, rowp(1:nx))
            rowsens(kk,m) = rowp(kk)
            rowmaxs(kk,m) = rescell
         enddo
         slope(kk) = rowsens(kk,5)/eps5(5)
         adm(kk)   = 0.0d0
         if (slope(kk) .gt. 0.0d0) adm(kk) = ieq_res_tol/slope(kk)
         write(*,'(3X,A8,ES11.3,5ES11.3,ES11.3,ES13.4)') rowname(kk),     &
              xv(kk), (rowsens(kk,m), m=1,5), slope(kk), adm(kk)
      enddo
      write(*,'(A)') '   the same displacements, largest residual of'//   &
           ' ANY row of that cell:'
      do kk = 1, nrow
         write(*,'(3X,A8,5ES11.3)') rowname(kk), (rowmaxs(kk,m), m=1,5)
      enddo
      call restore_module_state(f_ret)

      ! (5b) a controlled perturbation of the transported parent: x(H2) of
      ! the binding cell is what the sweep holds fixed there, so a relative
      ! change of it is a change the sweep propagates rather than erases.
      eps3(1) = 1.0d-8
      eps3(2) = 1.0d-7
      eps3(3) = 1.0d-6
      call save_module_state()
      f_probe = f_ret
      call ioniz_eq(T_col, rho, f_probe, heat_work, cool_work, eta_work)
      f_hold = f_probe
      call densities_of(f_probe)
      ne_unpert = ne_r(jbind)
      call restore_module_state(f_ret)
      write(*,'(A)') ' (D4a) TABLE 5c: one sweep after a relative'//      &
           ' change of x(H2) at the binding cell'
      write(*,'(A)') '     eps        d ln H2+     d ln H3+    '//        &
                     ' d ln HeH+    d ln n_e'
      do m = 1, 3
         call save_module_state()
         f_probe = f_ret
         f_probe(jbind,isp_H2) = f_ret(jbind,isp_H2)*(1.0d0 + eps3(m))
         f_probe(jbind,1)      = f_ret(jbind,1)                           &
                               - 2.0d0*f_ret(jbind,isp_H2)*eps3(m)
         call ioniz_eq(T_col, rho, f_probe, heat_work, cool_work,         &
                       eta_work)
         dperp(m,1) = relmove(f_hold, f_probe, jbind, isp_H2p)/eps3(m)
         dperp(m,2) = relmove(f_hold, f_probe, jbind, isp_H3p)/eps3(m)
         dperp(m,3) = relmove(f_hold, f_probe, jbind, isp_HeHp)/eps3(m)
         call densities_of(f_probe)
         dne(m) = abs(ne_r(jbind) - ne_unpert)/ne_unpert/eps3(m)
         call restore_module_state(f_ret)
         write(*,'(3X,ES10.2,4ES13.4)') eps3(m), dperp(m,1), dperp(m,2),  &
              dperp(m,3), dne(m)
      enddo


      ! ---- (6) the stopping rule against the residual at the exit ----- !
      !
      ! The closure ends on dT/T < chem_cycle_tol alone (the exit test of
      ! equilibrate_chemistry_at_fixed_conserved_state). Started again from
      ! the state it handed back, it is asked what it now reports, at what
      ! cost, and what its FIRST cycle still moves.
      call save_module_state()
      f_work    = f_ret
      heat_work = heat_col
      cool_work = cool_col
      eta_work  = eta_col
      call equilibrate_chemistry_at_fixed_conserved_state(u_col, f_work,  &
                p_work, T_work, heat_work, cool_work, eta_work, okc,      &
                n_cycles, reason, dT_inc)
      big = 0.0d0;  jbig = 0;  ispbig = 0
      do j = 1, N
         do isp = 1, n_species
            fref = f_ret(j,isp)
            if (fref .le. 1.0d-20) cycle
            mv = abs(f_work(j,isp) - fref)/fref
            if (mv .gt. big) then
               big = mv;  jbig = j;  ispbig = isp
            endif
         enddo
      enddo
      call restore_module_state(f_ret)
      write(*,'(A)') ' (D4a) TABLE 6: the stopping rule at the state it'//&
           ' handed back'
      write(*,'(A,L1,A,A,A,I0)') '   the closure restarted from the'//    &
           ' returned state reports ok = ', okc, ', reason ',             &
           trim(chem_closure_reason_text(reason)), ', cycles ', n_cycles
      write(*,'(A,ES12.4,A,ES10.3)') '   its temperature increment ',     &
           dT_inc, ' against chem_cycle_tol ', chem_cycle_tol
      write(*,'(A,ES12.4,A,I0,A,I0)') '   the composition it then'//      &
           ' stands at differs from the returned one by ', big,           &
           ' in f_sp column ', ispbig, ' at cell ', jbig
      write(*,'(A,ES12.4,A,ES10.3)') '   the HeH+ row residual of the'//  &
           ' returned composition at the returned temperature is ',       &
           rowres(7), ' against ieq_res_tol ', ieq_res_tol
      write(*,'(A,ES12.4)') '   and the largest row residual of any'//    &
           ' cell at the returned temperature is ', maxval(res_retT)

      ! THE CHAIN, per link, from the first complete increment.
      write(*,'(A)') ' (D4a) the chain of the first complete increment'// &
           ' at the binding cell, link by link:'
      write(*,'(A,ES11.3,A,ES11.3,A,F7.2)') '   H2+ ', chain(1),          &
           ' -> H3+ ', chain(2), ', factor ',                             &
           chain(2)/max(chain(1),1.0d-300)
      write(*,'(A,ES11.3,A,ES11.3,A,F7.2)') '   H3+ ', chain(2),          &
           ' -> HeH+ ', chain(3), ', factor ',                            &
           chain(3)/max(chain(2),1.0d-300)


      ! (6b) THE CLAIM IN THE ROUTINE'S OWN HEADER, tested at the
      ! tolerance the code uses. The header states that a temperature
      ! change below the sweep's reaction-residual tolerance cannot move a
      ! composition the sweep converged to that tolerance; the exit test is
      ! written on chem_cycle_tol, which is ten times ieq_res_tol. The
      ! returned composition is therefore swept again at temperatures
      ! displaced by 1e-7, 1e-6 and 1e-5 relative, and the movement is
      ! compared with the movement at the returned temperature itself.
      write(*,'(A)') ' (D4a) TABLE 6b: one sweep at a displaced'//        &
           ' temperature, movement of the returned composition'
      write(*,'(A)') '     dT/T        largest relative movement'//       &
                     '   species  cell    HeH+ at the binding cell'
      do m = 1, 4
         call save_module_state()
         f_probe = f_ret
         T_work = T_col*(1.0d0 + dTshift(m))
         call ioniz_eq(T_work, rho, f_probe, heat_work, cool_work,        &
                       eta_work)
         big = 0.0d0;  jbig = 0;  ispbig = 0
         do j = 1, N
            do isp = 1, n_species
               fref = f_ret(j,isp)
               if (fref .le. 1.0d-20) cycle
               mv = abs(f_probe(j,isp) - fref)/fref
               if (mv .gt. big) then
                  big = mv;  jbig = j;  ispbig = isp
               endif
            enddo
         enddo
         write(*,'(3X,ES10.2,ES22.4,I9,I7,ES20.4)') dTshift(m), big,      &
              ispbig, jbig, relmove(f_ret, f_probe, jbind, isp_HeHp)
         call restore_module_state(f_ret)
      enddo

      ! ---- (4) a refused trial and what it puts back ------------------ !
      !
      ! The refusal is forced through the closure's own cycle cap: with one
      ! cycle allowed the closure spends its budget on every trial, returns
      ! ok = .false. with the reason "cycle budget spent", and the pass
      ! refuses each trial on the chemistry. That is the production
      ! refusal path; a nonfinite or inadmissible state cannot be forced
      ! from outside, since no key manufactures one (reported in the memo).
      call seed_mass_row_of_the_column(v_relax)
      call seed_background_at_the_entry_composition()
      f_sp = f_sp0
      call carrier_checkpoint_restore(chk0)
      call conserved_of(v_relax)
      ! THE PASS IS RUN TWICE.  The first is entered with the pressure the
      ! column was seeded with through comp_p_from_T, while the pass's own
      ! first act is to take p and T from the caloric equation of state at
      ! the unchanged thermal energy, so p and T come back one round trip
      ! (1e-14) from the seeded pair and that difference is not a restore.
      ! The second is entered with exactly the pair a refused pass returns,
      ! so there every item of the list, p and T included, has to come back
      ! bit for bit.
      do m = 1, 2
         ! The state the next trial reads, recorded at the entry of the
         ! pass.
         call hold_next_operation_state()
         call carrier_checkpoint_take(chk_iso)
         chem_cycles_cap_for_test = 1
         call relax_photochemical_composition(u_col, v_relax, f_sp,       &
                                              p_col, T_col, heat_col,     &
                                              cool_col, eta_col, 1.0d-2,  &
                                              drift, ns, oc)
         chem_cycles_cap_for_test = -1
         write(*,'(A,I0,A,A,A,I0,A,A)') ' (D4a) forced-refusal pass ', m, &
              ' ends on ', trim(carrier_relax_outcome_text(oc)),          &
              ', kept ', ns, ' step(s); the closure reports ',            &
              trim(chem_closure_reason_text(n_chem_last_reason))
         if (m .eq. 1) then
            call check_absolute('a_forced_closure_refusal_ends_the'//     &
                 '_pass_on_the_chemistry', dble(oc),                      &
                 dble(carrier_relax_chemistry_refused), 0.0d0)
            call check_absolute('and_keeps_no_step', dble(ns), 0.0d0,     &
                 0.0d0)
         endif
         call report_restored_after_refusal(m .eq. 2)
      enddo


      ! ---- the reproduction, and how tightly it is pinned ------------- !
      !
      ! HOW THIS PROGRAM'S STATE DIFFERS FROM THE SUITE'S. The section-11
      ! sequence of carrier_retry.f90 re-seeds the mass row, the frozen
      ! background, the composition and the carrier checkpoint, but NOT the
      ! pressure: conserved_of builds the conserved state from whatever
      ! p_col the preceding passes of that program left, while here it is
      ! built from the pressure the column was seeded with. The sequence is
      ! therefore repeated, each repetition entering it with the pressure
      ! the previous one returned, and the spread of the row over the
      ! repetitions is the size of that dependence.
      write(*,'(A)') ' (D4a) the section-11 row over repetitions of the'//&
           ' sequence, the pressure entering it being the one the'//      &
           ' previous repetition returned:'
      do m = 1, 4
         call seed_mass_row_of_the_column(v_relax)
         call seed_background_at_the_entry_composition()
         f_sp = f_sp0
         call carrier_checkpoint_restore(chk0)
         call conserved_of(v_relax)
         call relax_photochemical_composition(u_col, v_relax, f_sp,       &
                   p_col, T_col, heat_col, cool_col, eta_col, 1.0d-2,     &
                   drift, ns, oc)
         f_work = f_sp
         call save_module_state()
         call ioniz_eq(T_col, rho, f_work, heat_work, cool_work, eta_work)
         big = 0.0d0;  jbig = 0;  ispbig = 0
         do j = 1, N
            do isp = 1, n_species
               fref = f_sp(j,isp)
               if (fref .le. 1.0d-20) cycle
               mv = abs(f_work(j,isp) - fref)/fref
               if (mv .gt. big) then
                  big = mv;  jbig = j;  ispbig = isp
               endif
            enddo
         enddo
         call restore_module_state(f_sp)
         write(*,'(A,I0,A,ES22.15,A,I0,A,I0,A,ES23.15E3)') '   pass ', m, &
              ': one further sweep moves ', big, ', f_sp column ',        &
              ispbig, ' at cell ', jbig, ', p(1) entering ', p_col(1)
      enddo

      if (assertion_failures .gt. 0) then
         write(*,'(A,I0)') ' closure_returned_state: FAILED, assertions', &
              assertion_failures
         stop 1
      endif
      write(*,'(A)') ' closure_returned_state: PASSED'

      contains

      !--------------!

      double precision function relmove(a, b, j, isp) result(d)
      ! The relative movement of one species of one cell between two
      ! compositions, zero where the reference is not present by the
      ! fixture's own 1e-20 rule.
      real*8, intent(in) :: a(1-Ng:,:), b(1-Ng:,:)
      integer, intent(in) :: j, isp
      d = 0.0d0
      if (a(j,isp) .gt. 1.0d-20) d = abs(b(j,isp) - a(j,isp))/a(j,isp)
      end function relmove

      !--------------!

      subroutine densities_of(fs)
      real*8, intent(in) :: fs(1-Ng:,:)
      nhei_r = 0.0d0;  nheii_r = 0.0d0
      nheiii_r = 0.0d0;  nheiTR_r = 0.0d0
      call get_species_densities(rho, fs, nhi_r, nhii_r, nhei_r,          &
                                 nheii_r, nheiii_r, nheiTR_r, nm_r,       &
                                 ne_r, ntot_r)
      end subroutine densities_of

      !--------------!

      subroutine hold_next_operation_state()
      ! EVERY PIECE OF MODULE STATE THE NEXT REAL OPERATION READS, held so
      ! that a refused trial can be measured against it. The list is read
      ! from the source:
      !   bg_cell                 ionization_equilibrium, the frozen cell
      !                           state the transport operator and the next
      !                           sweep read (rates, T_K, n_tot, the imposed
      !                           fractions and their flags)
      !   ieq_rate_cell, ieq_ne_cell, ieq_TK_cell, ieq_ntot_cell,
      !   ieq_met_coef, ieq_neq_stored, ieq_mbase_stored, ieq_iox_stored,
      !   ieq_rates_ready         the rate state the closure evaluator and
      !                           the certification read
      !   ieq_nonroot_streak      the non-root persistence counter
      !   nk_per_mass, x_h2, molecular_cell (caloric_eos, private): read
      !                           through pressure_from_energy_density, and
      !                           refreshed by get_species_densities
      !   dp_bc, n_part_cell1     global_parameters, the base ghost the
      !                           boundary condition reads
      !   the carrier checkpoint  diffusive_photochemistry
      call hold_arrays()
      p_hold = p_col;  T_hold = T_col
      heat_hold = heat_col;  cool_hold = cool_col;  eta_hold = eta_col
      f_hold = f_sp
      do j = 1-Ng, N+Ng
         p_cal1(j) = pressure_from_energy_density(j, u_col(1,j),          &
                        u_col(3,j) - 0.5d0*u_col(2,j)**2/u_col(1,j))
      enddo
      end subroutine hold_next_operation_state

      !--------------!

      subroutine report_restored_after_refusal(gate_primitive_pair)
      ! WHAT A REFUSED TRIAL PUTS BACK, item by item, each one read from
      ! the source rather than from the production enumeration, so that an
      ! item missing from that enumeration fails here.  Every difference is
      ! a gate: a trial that is undone leaves the state it was entered with
      ! and nothing of its own.  gate_primitive_pair is false on the pass
      ! entered with a seeded pressure, where p and T legitimately come back
      ! one caloric round trip away (see the call site).
      logical, intent(in) :: gate_primitive_pair
      real*8 :: dbg, drt, dne_, dTK, dnt, dmet, dcal, dsp_
      real*8 :: dstreak, dready, dstored
      integer :: jj
      dbg = 0.0d0;  drt = 0.0d0;  dne_ = 0.0d0;  dTK = 0.0d0
      dnt = 0.0d0;  dmet = 0.0d0;  dcal = 0.0d0
      dsp_ = maxval(abs(f_sp - f_hold))
      dstreak = dble(maxval(abs(ieq_nonroot_streak - streak_hold)))
      dready  = 0.0d0
      if (ieq_rates_ready .neqv. ready_hold) dready = 1.0d0
      dstored = dble(max(abs(ieq_neq_stored - neq_hold),                  &
                         abs(ieq_mbase_stored - mbase_hold),              &
                         abs(ieq_iox_stored - iox_hold)))
      do jj = 1-Ng, N+Ng
         dbg = max(dbg, bg_difference(bg_cell(jj), bg_hold(jj)))
         drt = max(drt, bg_difference(ieq_rate_cell(jj), rate_hold(jj)))
         dne_ = max(dne_, abs(ieq_ne_cell(jj) - ne_hold(jj)))
         dTK = max(dTK, abs(ieq_TK_cell(jj) - TK_hold(jj)))
         dnt = max(dnt, abs(ieq_ntot_cell(jj) - ntot_hold(jj)))
         dmet = max(dmet, maxval(abs(ieq_met_coef(:,:,jj)                 &
                                     - met_hold(:,:,jj))))
         dcal = max(dcal, abs(pressure_from_energy_density(jj,            &
                     u_col(1,jj), u_col(3,jj)                             &
                     - 0.5d0*u_col(2,jj)**2/u_col(1,jj)) - p_cal1(jj)))
      enddo
      write(*,'(A)') ' (D4a) TABLE 4: what the refused trials put back'// &
           ' (0 = bit for bit)'
      write(*,'(A,ES12.4)') '   composition f_sp                  ', dsp_
      write(*,'(A,ES12.4)') '   frozen cell state bg_cell         ', dbg
      write(*,'(A,ES12.4)') '   rate state ieq_rate_cell          ', drt
      write(*,'(A,ES12.4)') '   ieq_ne_cell                       ', dne_
      write(*,'(A,ES12.4)') '   ieq_TK_cell                       ', dTK
      write(*,'(A,ES12.4)') '   ieq_ntot_cell                     ', dnt
      write(*,'(A,ES12.4)') '   ieq_met_coef                      ', dmet
      write(*,'(A,ES12.4)') '   caloric state (p from rho e)      ', dcal
      write(*,'(A,ES12.4)') '   pressure p                        ',      &
           maxval(abs(p_col - p_hold))
      write(*,'(A,ES12.4)') '   temperature T                     ',      &
           maxval(abs(T_col - T_hold))
      write(*,'(A,ES12.4)') '   heating, cooling, eta             ',      &
           max(maxval(abs(heat_col - heat_hold)),                         &
               maxval(abs(cool_col - cool_hold)),                         &
               maxval(abs(eta_col - eta_hold)))
      write(*,'(A,ES12.4)') '   base ghost electron term dp_bc    ',      &
           abs(dp_bc - dpbc_hold)
      write(*,'(A,ES12.4)') '   base ghost particle count         ',      &
           abs(n_part_cell1 - npc1_hold)
      write(*,'(A,ES12.4)') '   non-root streak                   ',      &
           dstreak
      write(*,'(A,L1)')     '   carrier checkpoint matches        ',      &
           carrier_checkpoint_matches(chk_iso)
      call check_absolute('a_refused_trial_returns_the_entry'//           &
           '_composition', dsp_, 0.0d0, 0.0d0)
      call check_absolute('and_the_frozen_cell_state', dbg, 0.0d0, 0.0d0)
      call check_absolute('and_the_rate_state_the_closure_evaluator'//    &
           '_reads', drt, 0.0d0, 0.0d0)
      call check_absolute('and_the_electron_density_of_that_rate'//       &
           '_state', dne_, 0.0d0, 0.0d0)
      call check_absolute('and_its_temperature', dTK, 0.0d0, 0.0d0)
      call check_absolute('and_its_particle_count', dnt, 0.0d0, 0.0d0)
      call check_absolute('and_its_metal_coefficients', dmet, 0.0d0,      &
           0.0d0)
      call check_absolute('and_the_layout_the_rate_state_was_stored'//    &
           '_with', dready + dstored, 0.0d0, 0.0d0)
      call check_absolute('and_the_non_root_persistence_counter',         &
           dstreak, 0.0d0, 0.0d0)
      call check_absolute('and_the_caloric_map_from_thermal_energy'//     &
           '_to_pressure', dcal, 0.0d0, 0.0d0)
      call check_absolute('and_the_base_ghost_electron_term',             &
           abs(dp_bc - dpbc_hold), 0.0d0, 0.0d0)
      call check_absolute('and_the_base_ghost_particle_count',            &
           abs(n_part_cell1 - npc1_hold), 0.0d0, 0.0d0)
      call check_absolute('and_the_carrier_subsystem',                    &
           merge(0.0d0, 1.0d0, carrier_checkpoint_matches(chk_iso)),      &
           0.0d0, 0.0d0)
      if (gate_primitive_pair) then
         call check_absolute('and_the_pressure', maxval(abs(p_col         &
              - p_hold)), 0.0d0, 0.0d0)
         call check_absolute('and_the_temperature', maxval(abs(T_col      &
              - T_hold)), 0.0d0, 0.0d0)
         call check_absolute('and_the_heating_cooling_and_efficiency',    &
              max(maxval(abs(heat_col - heat_hold)),                      &
                  maxval(abs(cool_col - cool_hold)),                      &
                  maxval(abs(eta_col - eta_hold))), 0.0d0, 0.0d0)
      endif
      end subroutine report_restored_after_refusal

      !--------------!

      double precision function bg_difference(a, b) result(d)
      ! The largest absolute difference between two cell rate states, over
      ! every field the next operation reads.
      type(ion_rates), intent(in) :: a, b
      d = 0.0d0
      d = max(d, abs(a%P_HI - b%P_HI), abs(a%P_HeI - b%P_HeI))
      d = max(d, abs(a%P_HeII - b%P_HeII), abs(a%P_HeITR - b%P_HeITR))
      d = max(d, abs(a%P_H2 - b%P_H2), abs(a%k_LW - b%k_LW))
      d = max(d, abs(a%P_H2_di - b%P_H2_di), abs(a%P_H2_dd - b%P_H2_dd))
      d = max(d, abs(a%P_H2_nd - b%P_H2_nd))
      d = max(d, abs(a%rchiiB - b%rchiiB), abs(a%rcheiiB - b%rcheiiB))
      d = max(d, abs(a%rcheiiiB - b%rcheiiiB))
      d = max(d, abs(a%a_ion_HI - b%a_ion_HI))
      d = max(d, abs(a%a_ion_HeI - b%a_ion_HeI))
      d = max(d, abs(a%a_ion_HeII - b%a_ion_HeII))
      d = max(d, abs(a%kcx_He0_Hp - b%kcx_He0_Hp))
      d = max(d, abs(a%kcx_Hep_H0 - b%kcx_Hep_H0))
      d = max(d, abs(a%nh - b%nh), abs(a%nhe - b%nhe))
      d = max(d, abs(a%T_K - b%T_K), abs(a%ntot - b%ntot))
      d = max(d, abs(a%x_h2_fix - b%x_h2_fix))
      d = max(d, abs(a%x_hp_fix - b%x_hp_fix))
      d = max(d, abs(a%x_oh_fix - b%x_oh_fix))
      d = max(d, abs(a%x_h2o_fix - b%x_h2o_fix))
      end function bg_difference

      !--------------!

      subroutine hold_arrays()
      if (.not. allocated(bg_hold)) then
         allocate(bg_hold(1-Ng:N+Ng), rate_hold(1-Ng:N+Ng))
         allocate(ne_hold(1-Ng:N+Ng), TK_hold(1-Ng:N+Ng))
         allocate(ntot_hold(1-Ng:N+Ng))
         allocate(met_hold(n_melem,size(ieq_met_coef,2),1-Ng:N+Ng))
         allocate(streak_hold(1-Ng:N+Ng))
      endif
      bg_hold     = bg_cell
      rate_hold   = ieq_rate_cell
      ne_hold     = ieq_ne_cell
      TK_hold     = ieq_TK_cell
      ntot_hold   = ieq_ntot_cell
      met_hold    = ieq_met_coef
      streak_hold = ieq_nonroot_streak
      ready_hold  = ieq_rates_ready
      neq_hold    = ieq_neq_stored
      mbase_hold  = ieq_mbase_stored
      iox_hold    = ieq_iox_stored
      dpbc_hold   = dp_bc
      npc1_hold   = n_part_cell1
      end subroutine hold_arrays

      !--------------!

      subroutine save_module_state()
      call hold_arrays()
      end subroutine save_module_state

      !--------------!

      subroutine restore_module_state(fs)
      ! Put back every item of the list in hold_next_operation_state. The
      ! caloric state is private to caloric_eos and is restored the way the
      ! production code restores it, by refreshing it from the composition
      ! that is being reinstated (get_species_densities).
      real*8, intent(in) :: fs(1-Ng:,:)
      bg_cell            = bg_hold
      ieq_rate_cell      = rate_hold
      ieq_ne_cell        = ne_hold
      ieq_TK_cell        = TK_hold
      ieq_ntot_cell      = ntot_hold
      ieq_met_coef       = met_hold
      ieq_nonroot_streak = streak_hold
      ieq_rates_ready    = ready_hold
      ieq_neq_stored     = neq_hold
      ieq_mbase_stored   = mbase_hold
      ieq_iox_stored     = iox_hold
      dp_bc              = dpbc_hold
      n_part_cell1       = npc1_hold
      call densities_of(fs)
      end subroutine restore_module_state

      !--------------!

      subroutine write_returned_state()
      ! THE UNMODIFIED PRODUCTION RETURN, written so that a later reading
      ! stands on the same numbers: the composition of every physical cell,
      ! the temperature and pressure it was handed back with, the electron
      ! and total particle densities it implies, the fractions the sweep
      ! held imposed with their flags, and the radiation field and the rate
      ! state of the cell.  The file lands in the working directory, which
      ! run_closure_returned_state.sh sets to the suite's object directory.
      integer :: uu, jj, is
      open(newunit=uu, file='D4a_returned_state.txt', status='replace',   &
           action='write')
      write(uu,'(A)') '# The state equilibrate_chemistry_at_fixed'//      &
           '_conserved_state handed back on the carrier_retry column'
      write(uu,'(A,ES23.15E3)') '# ieq_res_tol   ', ieq_res_tol
      write(uu,'(A,ES23.15E3)') '# chem_cycle_tol', chem_cycle_tol
      write(uu,'(A,I0)')        '# chem_cycles_max ', chem_cycles_max
      write(uu,'(A,A)')  '# closure reason of the last trial: ',          &
           trim(chem_closure_reason_text(n_chem_last_reason))
      write(uu,'(A,ES23.15E3)') '# last temperature increment ',          &
           chem_last_increment
      write(uu,'(A)') '# columns: cell r p[code] T[code] T[K]'//          &
           ' n_e[cm^-3] n_tot[cm^-3] T_sweep[K] ne_sweep[cm^-3]'//        &
           ' x_h2_fixed x_h2_fix x_hp_fixed x_hp_fix P_HI P_HeI'//        &
           ' P_H2 k_LW rchiiB'
      do jj = 1, N
         write(uu,'(I5,17(1X,ES23.15E3))') jj, r(jj), p_col(jj),          &
              T_ret(jj), T_ret(jj)*T0, ne_r(jj)*n0, ntot_r(jj)*n0,        &
              ieq_TK_cell(jj), ieq_ne_cell(jj),                           &
              merge(1.0d0, 0.0d0, bg_cell(jj)%x_h2_fixed),                &
              bg_cell(jj)%x_h2_fix,                                       &
              merge(1.0d0, 0.0d0, bg_cell(jj)%x_hp_fixed),                &
              bg_cell(jj)%x_hp_fix, bg_cell(jj)%P_HI,                     &
              bg_cell(jj)%P_HeI, bg_cell(jj)%P_H2, bg_cell(jj)%k_LW,      &
              bg_cell(jj)%rchiiB
      enddo
      write(uu,'(A)') '# the composition, f_sp(cell, species)'
      do jj = 1, N
         do is = 1, n_species
            if (f_ret(jj,is) .eq. 0.0d0) cycle
            write(uu,'(2I6,1X,ES23.15E3)') jj, is, f_ret(jj,is)
         enddo
      enddo
      close(uu)
      write(*,'(A)') ' (D4a) the production return is written to'//       &
           ' D4a_returned_state.txt'
      end subroutine write_returned_state

      !--------------!
      ! The fixture below is the column of src/tests/carrier_retry/
      ! carrier_retry.f90, reproduced field for field so that this
      ! program measures the state that suite's section-11 row reads.
      !--------------!

      subroutine setup_globals()
      integer :: j
      N   = ncell
      n0  = 1.0d10
      R0  = 1.0d10
      v0  = 1.0d6
      T0  = 1000.0d0
      HeH = 0.0793d0
      b0  = 10.0d0
      p_base_bar         = 1.0d-6
      spherical_domain   = .true.
      thereis_He         = .true.
      thereis_HeITR      = .false.
      thereis_metals     = .true.
      thereis_mol        = .true.
      thereis_oxychem    = .false.
      eos_include_metals = .false.
      he_diffusion       = .false.
      carrier_transport  = .true.
      ionization_transport = .true.
      weno_mode          = 0
      j_min              = 1
      call allocate_grid_arrays
      do j = 1-Ng, N+Ng
         r(j)     = 1.0d0 + 0.02d0*dble(j)
         r_edg(j) = 1.0d0 + 0.02d0*(dble(j) + 0.5d0)
         dr_j(j)  = 0.02d0
      enddo
      kzz_cell = 0.0d0
      call set_gravity_grid
      grid_type  = 'Uniform'
      rec_method = 'WENO3'
      use_plm    = .false.
      use_weno3  = .true.
      flux       = 'ROE'
      CFL        = 0.6d0
      r_max      = r_edg(N)
      r_esc      = r(N)
      r_flux     = r(1)
      Mp         = 1.0d30
      q0         = 1.0d0
      allocate(melem_ab(n_melem))
      melem_ab = 0.0d0
      is_PL_sed    = .true.
      thereis_Xray = .true.
      e_low  = 13.60d0
      e_mid  = 123.98d0
      e_top  = 1.24d3
      PLind  = -1.0d0
      LX     = 27.20d0
      LEUV   = 27.93d0
      a_orb  = 0.0480d0*AU
      call set_energy_vectors
      N_eq = 7 + 2*n_melem
      lwa  = (N_eq*(3*N_eq + 13))/2
      allocate(sys_sol(N_eq), sys_x(N_eq), wa(lwa))
      call cx_init
      call carrier_set_init()
      call ioniz_eq_allocate_arrays()
      bg_ready = .true.
      end subroutine setup_globals

      !--------------!

      subroutine build_molecular_hydrogen_column()
      integer :: j
      allocate(rho(1-Ng:N+Ng), v(1-Ng:N+Ng), dt_code(1-Ng:N+Ng))
      allocate(f_sp(1-Ng:N+Ng,n_species), f_sp0(1-Ng:N+Ng,n_species))
      do j = 1-Ng, N+Ng
         rho(j) = exp(-4.0d0*(r(j) - 1.0d0))
         v(j)   = 0.02d0*r(j)
      enddo
      dt_code = 1.0d-3
      call column_carrying_its_own_density(f_sp, 0.8d0, .false., 0.0d0)
      f_sp(:,isp_HII)  = 1.0d-6
      f_sp(:,isp_HI)   = f_sp(:,isp_HI) - 1.0d-6
      f_sp(:,isp_H2p)  = 1.0d-10
      f_sp(:,isp_H3p)  = 1.0d-10
      f_sp(:,isp_H2)   = f_sp(:,isp_H2) - 2.5d-10
      f_sp(:,isp_HeII) = 1.0d-8
      f_sp(:,isp_HeI)  = f_sp(:,isp_HeI) - 1.0d-8
      do j = 1-Ng, N+Ng
         call set_frozen_cell_rates(j)
      enddo
      f_sp0 = f_sp
      allocate(p_col(1-Ng:N+Ng), T_col(1-Ng:N+Ng))
      allocate(heat_col(1-Ng:N+Ng), cool_col(1-Ng:N+Ng))
      allocate(eta_col(1-Ng:N+Ng))
      allocate(u_col(3,1-Ng:N+Ng), W_col(3,1-Ng:N+Ng))
      call seed_pressure_of_the_column()
      end subroutine build_molecular_hydrogen_column

      !--------------!

      subroutine seed_pressure_of_the_column()
      real*8, allocatable :: nhi_s(:), nhii_s(:), nhei_s(:), nheii_s(:)
      real*8, allocatable :: nheiii_s(:), nheiTR_s(:), ne_s(:), ntot_s(:)
      real*8, allocatable :: nm_s(:,:)
      integer :: j
      allocate(nhi_s(1-Ng:N+Ng), nhii_s(1-Ng:N+Ng))
      allocate(nhei_s(1-Ng:N+Ng), nheii_s(1-Ng:N+Ng))
      allocate(nheiii_s(1-Ng:N+Ng), nheiTR_s(1-Ng:N+Ng))
      allocate(ne_s(1-Ng:N+Ng), ntot_s(1-Ng:N+Ng))
      allocate(nm_s(1-Ng:N+Ng,n_mion))
      nhei_s   = 0.0d0
      nheii_s  = 0.0d0
      nheiii_s = 0.0d0
      nheiTR_s = 0.0d0
      call get_species_densities(rho, f_sp0, nhi_s, nhii_s, nhei_s,       &
                                 nheii_s, nheiii_s, nheiTR_s, nm_s,       &
                                 ne_s, ntot_s)
      do j = 1-Ng, N+Ng
         T_col(j) = bg_cell(j)%T_K/T0
      enddo
      call comp_p_from_T(T_col, ntot_s, ne_s, p_col)
      heat_col = 0.0d0
      cool_col = 0.0d0
      eta_col  = 0.0d0
      end subroutine seed_pressure_of_the_column

      !--------------!

      subroutine conserved_of(vv)
      real*8, dimension(1-Ng:N+Ng), intent(in) :: vv
      real*8, allocatable :: a1(:), a2(:), a3(:), a4(:), a5(:), a6(:)
      real*8, allocatable :: a7(:), a8(:), am(:,:)
      allocate(a1(1-Ng:N+Ng), a2(1-Ng:N+Ng), a3(1-Ng:N+Ng))
      allocate(a4(1-Ng:N+Ng), a5(1-Ng:N+Ng), a6(1-Ng:N+Ng))
      allocate(a7(1-Ng:N+Ng), a8(1-Ng:N+Ng), am(1-Ng:N+Ng,n_mion))
      a3 = 0.0d0;  a4 = 0.0d0;  a5 = 0.0d0;  a6 = 0.0d0
      call get_species_densities(rho, f_sp0, a1, a2, a3, a4, a5, a6, am,  &
                                 a7, a8)
      W_col(1,:) = rho
      W_col(2,:) = vv
      W_col(3,:) = p_col
      call W_to_U(W_col, u_col)
      end subroutine conserved_of

      !--------------!

      subroutine seed_background_at_the_entry_composition()
      real*8, allocatable :: nhi_e(:), nhii_e(:), nhei_e(:), nheii_e(:)
      real*8, allocatable :: nheiii_e(:), nheiTR_e(:), ne_e(:)
      real*8, allocatable :: nm_e(:,:)
      integer :: j
      if (.not. allocated(ntot0_c)) then
         allocate(ntot0_c(1-Ng:N+Ng), T_ent(1-Ng:N+Ng))
      endif
      allocate(nhi_e(1-Ng:N+Ng), nhii_e(1-Ng:N+Ng))
      allocate(nhei_e(1-Ng:N+Ng), nheii_e(1-Ng:N+Ng))
      allocate(nheiii_e(1-Ng:N+Ng), nheiTR_e(1-Ng:N+Ng))
      allocate(ne_e(1-Ng:N+Ng), nm_e(1-Ng:N+Ng,n_mion))
      nhei_e   = 0.0d0
      nheii_e  = 0.0d0
      nheiii_e = 0.0d0
      nheiTR_e = 0.0d0
      call get_species_densities(rho, f_sp0, nhi_e, nhii_e, nhei_e,       &
                                 nheii_e, nheiii_e, nheiTR_e, nm_e,       &
                                 ne_e, ntot0_c)
      call comp_T_from_p(p_col, ntot0_c, ne_e, T_ent)
      do j = 1-Ng, N+Ng
         call set_frozen_cell_rates(j)
         bg_cell(j)%ntot = ntot0_c(j)*n0
         bg_cell(j)%T_K  = T_ent(j)*T0
         T_col(j)        = T_ent(j)
      enddo
      heat_col = 0.0d0
      cool_col = 0.0d0
      eta_col  = 0.0d0
      end subroutine seed_background_at_the_entry_composition

      !--------------!

      subroutine seed_mass_row_of_the_column(v_wind)
      real*8, intent(in)  :: v_wind(1-Ng:N+Ng)
      real*8, allocatable :: u(:,:), Res(:,:)
      real*8, allocatable :: n_part(:), heat(:), cool(:), p_s(:), T_s(:)
      integer :: j
      allocate(u(3,1-Ng:N+Ng), Res(3,1-Ng:N+Ng))
      allocate(n_part(1-Ng:N+Ng), heat(1-Ng:N+Ng), cool(1-Ng:N+Ng))
      allocate(p_s(1-Ng:N+Ng), T_s(1-Ng:N+Ng))
      heat = 0.0d0
      cool = 0.0d0
      do j = 1-Ng, N+Ng
         T_s(j)    = bg_cell(j)%T_K/T0
         n_part(j) = rho(j)
         p_s(j)    = n_part(j)*T_s(j)
         u(1,j)    = rho(j)
         u(2,j)    = rho(j)*v_wind(j)
         u(3,j)    = 0.5d0*rho(j)*v_wind(j)**2                            &
                     + p_s(j)/(gamma_ad - 1.0d0)
      enddo
      n_part_cell1 = n_part(1)
      call set_base_reservoir(p_s(1), T_s(1), 1.0d0, 1.0d0)
      call Apply_BC(u)
      call assemble_residual(u, n_part, heat, cool, Res)
      end subroutine seed_mass_row_of_the_column

      !--------------!

      subroutine set_frozen_cell_rates(j)
      integer, intent(in) :: j
      bg_cell(j)%P_HI        = 1.0d-6
      bg_cell(j)%P_HeI       = 1.0d-7
      bg_cell(j)%P_HeII      = 0.0d0
      bg_cell(j)%P_HeITR     = 0.0d0
      bg_cell(j)%P_H2        = 1.0d-8
      bg_cell(j)%P_H2_di     = 0.0d0
      bg_cell(j)%P_H2_dd     = 0.0d0
      bg_cell(j)%P_H2_nd     = 0.0d0
      bg_cell(j)%k_LW        = 0.0d0
      bg_cell(j)%rchiiB      = 2.6d-13
      bg_cell(j)%rcheiiB     = 4.3d-13
      bg_cell(j)%rcheiiiB    = 2.2d-12
      bg_cell(j)%rcheiTR     = 0.0d0
      bg_cell(j)%a_ion_HI    = 0.0d0
      bg_cell(j)%a_ion_HeI   = 0.0d0
      bg_cell(j)%a_ion_HeII  = 0.0d0
      bg_cell(j)%a_ion_HeITR = 0.0d0
      bg_cell(j)%q13         = 0.0d0
      bg_cell(j)%q31a        = 0.0d0
      bg_cell(j)%q31b        = 0.0d0
      bg_cell(j)%Q31         = 0.0d0
      bg_cell(j)%A31         = 0.0d0
      bg_cell(j)%kcx_He0_Hp  = 0.0d0
      bg_cell(j)%kcx_Hep_H0  = 0.0d0
      bg_cell(j)%nh          = 0.0d0
      bg_cell(j)%nhe         = 0.0d0
      bg_cell(j)%n_ofam      = 0.0d0
      bg_cell(j)%n_co        = 0.0d0
      bg_cell(j)%x_h2_fixed  = .false.
      bg_cell(j)%x_ox_fixed  = .false.
      bg_cell(j)%x_hp_fixed  = .false.
      bg_cell(j)%x_h2_fix    = 0.0d0
      bg_cell(j)%x_oh_fix    = 0.0d0
      bg_cell(j)%x_h2o_fix   = 0.0d0
      bg_cell(j)%x_hp_fix    = 0.0d0
      bg_cell(j)%T_K         = 1500.0d0
      bg_cell(j)%ntot        = 0.5d0*rho(j)*n0
      end subroutine set_frozen_cell_rates

      end program closure_returned_state
