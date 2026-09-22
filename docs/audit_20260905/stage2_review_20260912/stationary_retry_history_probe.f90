! Diagnostic written on 2026-09-12.
! Production routines are unchanged. Existing fixture statements are retained;
! original full-line comments are omitted. Added AUDIT rows test returned-state contracts.
      program stationary_retry_history_probe

      use global_parameters
      use species_table
      use ionization_equilibrium, only: bg_cell, bg_ready,                &
                                        ioniz_eq_allocate_arrays, ioniz_eq
      use diffusive_photochemistry, only:                                 &
           carrier_set_init, carrier_transport_interval,                  &
           carrier_checkpoint, carrier_checkpoint_take,                   &
           carrier_checkpoint_restore, carrier_checkpoint_matches,        &
           carrier_perturb_checkpointed_state_for_test,                   &
           carrier_reject_leading_attempts_for_test,                      &
           carrier_initial_substep_for_test,                              &
           carrier_attempt_record, carrier_co_domain_record,              &
           carrier_co_domain_perturb_for_test, carrier_retry_max,         &
           carrier_verdict, n_carrier_max, ic_H2, ic_Hp,                  &
           carrier_row_scales, carrier_returned_state_verdict,            &
           carrier_solve_converged, carrier_stall_at_floor,             &
           carrier_transport_diagnostics,                               &
           carrier_last_interval_halvings, carrier_row_roundoff,        &
           carrier_reject_row_residual, carrier_refused_rows,           &
           carrier_refused_row, carrier_roundoff_limited_record,         &
           carrier_record_refused_rows, carrier_name,                     &
           photochemical_transport_step, carrier_interval_covered,        &
           carrier_interval_exhausted, carrier_last_interval_status,      &
           carrier_exhausted_record, carrier_history_certifiable,         &
           carrier_transport_stop_on_failure,                             &
           carrier_transport_stops_suppressed,                            &
           carrier_step_verdict, carrier_last_verdict_of_run,             &
           carrier_no_interval,                                           &
           relax_photochemical_composition,                               &
           carrier_relax_movement_bound, carrier_relax_fixed_point,       &
           carrier_relax_step_budget, carrier_relax_interval_refused,     &
           carrier_relax_nothing_to_advance, carrier_relax_outcome_text
      use element_census, only: element_nuclei_and_charge, n_element
      use gravity_grid_construction, only: set_gravity_grid
      use energy_vectors_construct,  only: set_energy_vectors
      use charge_exchange,           only: cx_init
      use composition,               only: get_species_densities,         &
                                           comp_p_from_T, comp_T_from_p
      use base_boundary,             only: set_base_reservoir
      use BC_Apply,                  only: Apply_BC
      use steady_residual_mod,       only: assemble_residual
      use caloric_eos, only: energy_density_from_pressure, pressure_from_energy_density
      use certification, only: mass_row_cell_verdict
      use assertion_report
      implicit none

      integer, parameter :: ncell = 12
      real*8,  allocatable :: rho(:), v(:), dt_code(:)
      real*8,  allocatable :: f_sp(:,:), f_sp0(:,:)
      real*8,  allocatable :: f_one(:,:), f_half(:,:), f_quarter(:,:)
      real*8,  allocatable :: f_retried(:,:)
      real*8,  allocatable :: fc_probe(:,:)
      type(carrier_checkpoint) :: chk0, chk_fresh, chk_after
      type(carrier_verdict)    :: verdict
      logical :: completed, ok
      real*8  :: frac_done, e1, e2
      integer :: nsub
      integer :: nattempt0, nreject0, naccept0, nretried0
      integer :: nattempt,  nreject,  naccept,  nretried
      integer :: nout_co0, nhot_co0, nhep_co0
      integer :: nout_co,  nhot_co,  nhep_co
      real*8,  allocatable :: tfull_dt(:,:), tphys_dt(:,:)
      real*8,  allocatable :: tfull_short(:,:), tphys_short(:,:)
      real*8,  allocatable :: res_dt(:,:), res_short(:,:)
      real*8,  allocatable :: dt_short(:)
      logical, allocatable :: unconstrained(:)
      type(carrier_verdict) :: vd_full_old, vd_full_new
      type(carrier_verdict) :: vd_short_old, vd_short_new
      real*8  :: ratio_co0, r_co0, form_co0
      real*8  :: ratio_co,  r_co,  form_co
      real*8,  allocatable :: dt_1024(:)
      real*8  :: rw_dt, rw_16, rw_1024
      real*8  :: rs_dt, rs_16, rs_1024, worst_lim
      integer :: nit_dt, nit_16, nit_1024, nlim
      integer :: nsub_dt, nsub_16, nsub_1024
      integer :: nhalve_1024, nrl_1024, nrl_last
      integer :: nrl_total, nrl_substeps
      integer :: nrefused, jref, icref
      real*8  :: res_ref, phys_ref, full_ref
      real*8,  allocatable :: tf_round(:,:)
      type(carrier_verdict) :: vd_round, vd_mixed
      logical :: cov_dt, cov_16, cov_1024
      integer :: st_phys, st_init, nstop0, nstop_phys, nstop_init
      integer :: nex0, nex_phys, nex_init
      integer :: ex_first, ex_last, ex_j, ex_ic
      real*8  :: ex_ratio, ex_phys, ex_full
      integer :: dout_tot, dhot_tot, dhep_tot
      integer :: dout_f1,  dhot_f1,  dhep_f1
      integer :: dout_f2,  dhot_f2,  dhep_f2
      real*8  :: drat_tot, dr_tot, dform_tot
      real*8  :: drat_f1,  dr_f1,  dform_f1
      real*8  :: drat_f2,  dr_f2,  dform_f2
      integer :: st_nothing, reason_before
      type(carrier_verdict) :: vd_step, vd_run
      real*8  :: dr_wide, dr_tight, dr_none, dr_zero, dr_ref
      integer :: ns_wide, ns_tight, ns_none, ns_zero, ns_ref
      integer :: oc_wide, oc_tight, oc_none, oc_zero, oc_ref
      real*8, allocatable :: nnuc0(:,:), nnuc1(:,:)
      real*8, allocatable :: nchg0(:), nchg1(:), rcomp0(:), rcomp1(:)
      real*8, allocatable :: v_relax(:), v_fast(:)
      real*8, allocatable :: p_col(:), T_col(:)
      real*8, allocatable :: heat_col(:), cool_col(:), eta_col(:)
      real*8, allocatable :: f_swept(:,:)
      real*8, allocatable :: nhi_r(:), nhii_r(:), nhei_r(:), nheii_r(:)
      real*8, allocatable :: nheiii_r(:), nheiTR_r(:), ne_r(:), ntot_r(:)
      real*8, allocatable :: nm_r(:,:), T_ret(:)
      real*8, allocatable :: ntot0_c(:), T_ent(:)
      real*8  :: move_pass, move_sweep, dev_ntot, dev_TK
      integer :: jj
      real*8, allocatable :: audit_energy(:)
      real*8 :: audit_tol, audit_dist, audit_p_error, audit_energy_error
      logical :: audit_within, audit_anchored

      call setup_globals()
      call build_molecular_hydrogen_column()
      allocate(nnuc0(1-Ng:N+Ng,n_element),nchg0(1-Ng:N+Ng),rcomp0(1-Ng:N+Ng))
      call element_nuclei_and_charge(rho,f_sp0,nnuc0,nchg0,rcomp0)
      do jj=1-Ng,N+Ng
         f_sp0(jj,:)=f_sp0(jj,:)*(rho(jj)*n0/rcomp0(jj))
      enddo
      f_sp=f_sp0
      call seed_pressure_of_the_column()
      call element_nuclei_and_charge(rho,f_sp0,nnuc0,nchg0,rcomp0)
      write(*,'(A,ES24.16)') 'AUDIT entry_mass_closure=', &
           maxval(abs(rcomp0(1:N)/(rho(1:N)*n0)-1d0))
      allocate(v_fast(1-Ng:N+Ng))
      v_fast=2000d0*r
      call seed_mass_row_of_the_column(v_fast)
      run_mode=run_mode_init
      f_sp=f_sp0
      carrier_reject_leading_attempts_for_test=1000000
      write(*,'(A,L1)') 'AUDIT history_before=',carrier_history_certifiable()
      call relax_photochemical_composition(rho,v_fast,p_col,T_col,f_sp, &
           heat_col,cool_col,eta_col,1d0,dr_ref,ns_ref,oc_ref)
      write(*,'(A,L1)') 'AUDIT history_after=',carrier_history_certifiable()
      write(*,'(A,ES24.16,A,I0,A,I0)') 'AUDIT composition_change=',maxval(abs(f_sp-f_sp0)), &
           ' kept_steps=',ns_ref,' outcome=',oc_ref
      carrier_reject_leading_attempts_for_test=0

      contains


      function roundoff_of(res, tfull) result(flag)
      real*8, intent(in) :: res(:,:), tfull(:,:)
      logical :: flag(size(res,1),size(res,2))
      flag = (abs(res) .le. carrier_row_roundoff*tfull)
      end function roundoff_of


      double precision function carrier_row_difference(a, b) result(d)
      real*8, intent(in) :: a(1-Ng:,:), b(1-Ng:,:)
      integer :: k, isp(6)
      real*8  :: sc
      isp = (/ isp_HI, isp_HII, isp_H2, isp_OH, isp_H2O, isp_CO /)
      d = 0.0d0
      do k = 1, 6
         sc = maxval(abs(f_sp0(1:N,isp(k))))
         if (sc .le. 0.0d0) cycle
         d = max(d, maxval(abs(a(1:N,isp(k)) - b(1:N,isp(k))))/sc)
      enddo
      end function carrier_row_difference


      double precision function element_deviation(a, b) result(d)
      real*8, intent(in) :: a(1-Ng:,:), b(1-Ng:,:)
      integer :: j, ie
      d = 0.0d0
      do ie = 1, n_element
         do j = 1, N
            if (a(j,ie) .le. 0.0d0) cycle
            d = max(d, abs(b(j,ie) - a(j,ie))/a(j,ie))
         enddo
      enddo
      end function element_deviation


      double precision function logical_as_double(l) result(x)
      logical, intent(in) :: l
      x = 0.0d0
      if (l) x = 1.0d0
      end function logical_as_double


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


      subroutine build_molecular_hydrogen_column()
      integer :: j
      allocate(rho(1-Ng:N+Ng), v(1-Ng:N+Ng), dt_code(1-Ng:N+Ng))
      allocate(f_sp(1-Ng:N+Ng,n_species), f_sp0(1-Ng:N+Ng,n_species))
      allocate(f_one(1-Ng:N+Ng,n_species), f_half(1-Ng:N+Ng,n_species))
      allocate(f_quarter(1-Ng:N+Ng,n_species))
      allocate(f_retried(1-Ng:N+Ng,n_species))
      do j = 1-Ng, N+Ng
         rho(j) = exp(-4.0d0*(r(j) - 1.0d0))
         v(j)   = 0.02d0*r(j)
      enddo
      dt_code = 1.0d-3
      f_sp = 0.0d0
      f_sp(:,isp_HI)   = 0.20d0
      f_sp(:,isp_HII)  = 1.0d-6
      f_sp(:,isp_H2)   = 0.40d0
      f_sp(:,isp_H2p)  = 1.0d-10
      f_sp(:,isp_H3p)  = 1.0d-10
      f_sp(:,isp_HeI)  = 0.0793d0
      f_sp(:,isp_HeII) = 1.0d-8
      do j = 1-Ng, N+Ng
         call set_frozen_cell_rates(j)
      enddo
      f_sp0 = f_sp
      allocate(p_col(1-Ng:N+Ng), T_col(1-Ng:N+Ng))
      allocate(heat_col(1-Ng:N+Ng), cool_col(1-Ng:N+Ng))
      allocate(eta_col(1-Ng:N+Ng))
      call seed_pressure_of_the_column()
      end subroutine build_molecular_hydrogen_column


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


      subroutine restore_frozen_background()
      integer :: j
      do j = 1-Ng, N+Ng
         call set_frozen_cell_rates(j)
         T_col(j) = bg_cell(j)%T_K/T0
      enddo
      heat_col = 0.0d0
      cool_col = 0.0d0
      eta_col  = 0.0d0
      end subroutine restore_frozen_background


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


      subroutine set_frozen_cell_rates(j)
      integer, intent(in) :: j
      bg_cell(j)%P_HI        = 1.0d-6      ! [1/s]
      bg_cell(j)%P_HeI       = 1.0d-7
      bg_cell(j)%P_HeII      = 0.0d0
      bg_cell(j)%P_HeITR     = 0.0d0
      bg_cell(j)%P_H2        = 1.0d-8
      bg_cell(j)%P_H2_di     = 0.0d0
      bg_cell(j)%P_H2_dd     = 0.0d0
      bg_cell(j)%P_H2_nd     = 0.0d0
      bg_cell(j)%k_LW        = 0.0d0
      bg_cell(j)%rchiiB      = 2.6d-13     ! [cm^3/s] at ~1e4 K
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

      end program stationary_retry_history_probe

