      program diffusion_tests
      ! Acceptance tests for the binary H/He element-diffusion operator
      ! (src/modules/functions/binary_element_diffusion.f90).  Test ids are
      ! those of docs/binary_diffusion_design.md section 6:
      !
      !   T1a  diffusive equilibrium in a closed isothermal column
      !   T1b  Dirichlet reservoir base: boundary flux accounts for the mass
      !   T3   hydrogen-trace limit (He/H = 1000)
      !   T4   elemental conservation, closed column, both eos_metals settings
      !   T5   two-component mass closure of the write-back
      !   T6   uniform mixture preserved under an arbitrary velocity field
      !   T9   ambipolar limits: relative settling mass 3, 2.5, 5/3
      !   T11  a cell whose face velocities straddle zero stays advectively
      !        coupled to its donor (the base-dip regression)
      !   T10  grid and time-step convergence of T1a
      !   T12  stage-resolved friction: the three limits of the effective
      !        binary coefficient (neutral, fully ionized, half ionized)
      !   T7   molecular closure and homopause: the Blanc carrier mixture in
      !        its two limits, and the homopause radius against K_zz
      !   T2a  recorded only: the kernel it compared against is gone
      !
      ! T0 (`make check` byte-identical with diffusion off), T2b and T8 are
      ! run outside this driver.
      !
      ! The columns are synthetic: the driver sets the global_parameters
      ! scalars and the radial grid itself, so no input.inp, no ionization
      ! solve and no hydro are involved.

      use global_parameters
      use species_table, only: isp_HI, isp_HII, isp_HeI, isp_HeII,        &
                               isp_HeIII, n_bsp, bsp_fsp, bsp_nH, bsp_nHe,&
                               bsp_is_excited_level,                       &
                               n_melem, melem_i0, mion_fsp, melem_A,     &
                               bsp_charge, isp_H2
      use composition,   only: mass_per_H_nucleus_without_He
      use lower_atmosphere_profile, only: eddy_diffusion_on_grid
      use binary_element_diffusion, only: element_diffusion_step,         &
                                          relative_settling_mass,         &
                                          helium_hydrogen_diffusion,      &
                                          alpha_HI, alpha_HeI

      implicit none

      real*8, allocatable :: rho_a(:), v_a(:), T_a(:), dt_a(:), f_a(:,:)
      real*8, allocatable :: Jf_a(:)
      integer :: n_pass, n_fail

      n_pass = 0
      n_fail = 0

      write(*,'(A)') '===== binary_element_diffusion acceptance tests ====='

      call test_T1a_and_T10()
      call test_T1b()
      call test_T3()
      call test_T4()
      call test_T5()
      call test_T6()
      call test_T9()
      call test_T11()
      call test_T12()
      call test_T7()
      call test_T2a()

      write(*,'(A)') '===================================================='
      write(*,'(A,I0,A,I0,A)') ' TOTAL: ', n_pass, ' passed, ',           &
                               n_fail, ' failed'
      if (n_fail .gt. 0) stop 1

      contains

      ! ================================================================= !
      !  harness
      ! ================================================================= !

      subroutine verdict(name, ok, val, tol, what)
      character(len=*), intent(in) :: name, what
      logical,          intent(in) :: ok
      real*8,           intent(in) :: val, tol
      if (ok) then
         n_pass = n_pass + 1
         write(*,'(A,A,A,A,ES10.3,A,ES10.3,A)') ' PASS  ', name, '  ',    &
              what, val, '  (tol ', tol, ')'
      else
         n_fail = n_fail + 1
         write(*,'(A,A,A,A,ES10.3,A,ES10.3,A)') ' FAIL  ', name, '  ',    &
              what, val, '  (tol ', tol, ')'
      endif
      end subroutine verdict

      ! ----------------------------------------------------------------- !

      subroutine setup_column(Ncell, r_top, T_base, R_pl, b_jeans, n_base)
      ! Build a synthetic column: N cells uniform in r from 1 to r_top, and
      ! the global scalars the operator reads.  b_jeans is b0 = G M mu/(k T R)
      ! set directly, so the gravity is g(r) = b0 v0^2/(R0 r^2) without a
      ! planet mass having to be invented.
      integer, intent(in) :: Ncell
      real*8,  intent(in) :: r_top, T_base, R_pl, b_jeans, n_base
      real*8  :: dr_u
      integer :: j

      N  = Ncell
      T0 = T_base
      R0 = R_pl
      n0 = n_base
      v0 = sqrt(kb_erg*T0/mu)
      t_s = R0/v0
      p0 = n0*mu*v0*v0
      q0 = n0*mu*v0*v0*v0/R0
      b0 = b_jeans
      spherical_domain = .true.
      thereis_He    = .true.
      thereis_HeITR = .false.
      thereis_mol   = .false.
      thereis_metals = .false.
      eos_include_metals = .true.
      he_diffusion  = .true.
      he_kzz        = 0.0d0
      call eddy_diffusion_on_grid
      he_ambipolar  = .false.
      he_alphaT     = 0.0d0
      he_metal_diffusion = .false.

      if (allocated(r))     deallocate(r)
      if (allocated(r_edg)) deallocate(r_edg)
      if (allocated(dr_j))  deallocate(dr_j)
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      if (.not. allocated(melem_ab)) allocate(melem_ab(n_melem))
      melem_ab = 0.0d0

      dr_u = (r_top - 1.0d0)/dble(N-1)
      do j = 1-Ng, N+Ng
         r(j) = 1.0d0 + dble(j-1)*dr_u
      enddo
      r_edg(1-Ng:N+Ng-1) = 0.5d0*(r(1-Ng:N+Ng-1) + r(2-Ng:N+Ng))
      r_edg(N+Ng) = 2.0d0*r_edg(N+Ng-1) - r_edg(N+Ng-2)
      dr_j = dr_u

      call alloc_state()

      end subroutine setup_column

      ! ----------------------------------------------------------------- !

      subroutine alloc_state()
      if (allocated(rho_a)) deallocate(rho_a, v_a, T_a, dt_a, f_a)
      if (allocated(Jf_a))  deallocate(Jf_a)
      allocate(rho_a(1-Ng:N+Ng), v_a(1-Ng:N+Ng), T_a(1-Ng:N+Ng),          &
               dt_a(1-Ng:N+Ng), f_a(1-Ng:N+Ng,n_species), Jf_a(0:N))
      rho_a = 1.0d0
      v_a   = 0.0d0
      T_a   = 1.0d0
      dt_a  = 1.0d0
      f_a   = 0.0d0
      end subroutine alloc_state

      ! ----------------------------------------------------------------- !

      subroutine set_composition(f_sp, heh_cell, xi_H, he_stage)
      ! Fill f_sp from a He/H element ratio given cell by cell, a hydrogen
      ! ionization fraction and a helium charge stage (0 = HeI, 1 = HeII,
      ! 2 = HeIII), normalized as the code does: sum_s m_s f_s = 1 over the species
      ! calc_rho counts.  The metals, when present, are seeded at the
      ! reservoir metal/H in their neutral stage.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(out) :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: heh_cell
      real*8,                                 intent(in)  :: xi_H
      integer,                                intent(in)  :: he_stage
      real*8  :: mpH, m_1
      integer :: j, ie

      m_1 = mass_per_H_nucleus_without_He()
      f_sp = 0.0d0
      do j = 1-Ng, N+Ng
         mpH = m_1 + 4.0d0*heh_cell(j)
         f_sp(j,isp_HI)  = (1.0d0 - xi_H)/mpH
         f_sp(j,isp_HII) = xi_H/mpH
         if (he_stage .eq. 0) f_sp(j,isp_HeI)   = heh_cell(j)/mpH
         if (he_stage .eq. 1) f_sp(j,isp_HeII)  = heh_cell(j)/mpH
         if (he_stage .eq. 2) f_sp(j,isp_HeIII) = heh_cell(j)/mpH
         if (thereis_metals) then
            do ie = 1, n_melem
               f_sp(j,mion_fsp(melem_i0(ie))) = melem_ab(ie)/mpH
            enddo
         endif
      enddo

      end subroutine set_composition

      ! ----------------------------------------------------------------- !

      function nucleus_fraction_He(f_sp) result(xnum)
      ! x = n_He/(n_H + n_He), the helium fraction of the collision partners.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8, dimension(1-Ng:N+Ng) :: xnum, nH_l, nHe_l
      call nucleus_counts(f_sp, nH_l, nHe_l)
      xnum = nHe_l/max(nH_l + nHe_l, 1.0d-300)
      end function nucleus_fraction_He

      ! ----------------------------------------------------------------- !

      subroutine nucleus_counts(f_sp, nH_l, nHe_l)
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: nH_l, nHe_l
      integer :: ib
      nH_l  = 0.0d0
      nHe_l = 0.0d0
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         if (bsp_nH(ib)  .gt. 0)                                          &
            nH_l  = nH_l  + dble(bsp_nH(ib)) *f_sp(:,bsp_fsp(ib))
         if (bsp_nHe(ib) .gt. 0)                                          &
            nHe_l = nHe_l + dble(bsp_nHe(ib))*f_sp(:,bsp_fsp(ib))
      enddo
      end subroutine nucleus_counts

      ! ----------------------------------------------------------------- !

      function helium_mass_in_column(f_sp, jlo) result(mtot)
      ! sum rho X r^2 dr over cells jlo..N, the discrete form of the integral
      ! the conservation tests quote (same volume weight r^2 (r_edg(j) -
      ! r_edg(j-1)) the operator's finite volumes use).
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      integer,                                intent(in) :: jlo
      real*8  :: mtot
      real*8, dimension(1-Ng:N+Ng) :: nH_l, nHe_l
      integer :: j
      call nucleus_counts(f_sp, nH_l, nHe_l)
      mtot = 0.0d0
      do j = jlo, N
         mtot = mtot + rho_a(j)*n0*mu*4.0d0*nHe_l(j)                      &
                      *(r(j)*R0)**2*(r_edg(j)-r_edg(j-1))*R0
      enddo
      end function helium_mass_in_column

      ! ================================================================= !
      !  T1a  diffusive equilibrium + T10 convergence order
      ! ================================================================= !

      subroutine test_T1a_and_T10()
      real*8  :: err(3), dr_h(3), ord21, ord32, mass0, mass1, dmass
      real*8  :: dts(4)
      integer :: ic, ncl(3)

      write(*,'(A)') ' --- T1a: diffusive equilibrium, closed isothermal '//&
                     'column (v=0, K=0)'
      ncl = [100, 200, 400]
      do ic = 1, 3
         call equilibrium_column(ncl(ic), 1.0d6, 4000, err(ic), dr_h(ic), &
                                 mass0, mass1)
         write(*,'(A,I4,A,ES10.3,A,ES10.3)') '       N = ', ncl(ic),      &
              '   max|logit(x) - analytic| = ', err(ic),                  &
              '   He mass drift = ', abs(mass1-mass0)/mass0
      enddo
      call verdict('T1a ', err(3) .lt. 2.0d-3, err(3), 2.0d-3,            &
                   'logit(x) error at N=400: ')

      ! Helium mass in the closed column.  The operator conserves it exactly;
      ! what leaks is the round-off of the tridiagonal solve, amplified by the
      ! conditioning of the implicit step, which grows with dt/t_diff(cell).
      ! The relaxation runs above deliberately use dt = 1e6 (about 1e5 cell
      ! diffusion times) to reach equilibrium quickly, so their drift is that
      ! amplified round-off, not a leak; the table below shows the scaling and
      ! the gate is applied at a resolved step.
      write(*,'(A)') '       He mass drift after 1000 steps vs step size:'
      dts = [1.0d1, 1.0d3, 1.0d5, 1.0d7]
      do ic = 1, 4
         call equilibrium_column(200, dts(ic), 1000, err(1), dr_h(1),     &
                                 mass0, mass1)
         write(*,'(A,ES9.2,A,ES10.3)') '         dt_code = ', dts(ic),    &
              '   drift = ', abs(mass1-mass0)/mass0
         if (ic .eq. 1) dmass = abs(mass1-mass0)/mass0
      enddo
      call verdict('T1a ', dmass .lt. 1.0d-13, dmass, 1.0d-13,            &
                   'closed-column He mass drift at dt_code = 10: ')

      write(*,'(A)') ' --- T10: grid convergence of T1a'
      ord21 = log(err(1)/err(2))/log(dr_h(1)/dr_h(2))
      ord32 = log(err(2)/err(3))/log(dr_h(2)/dr_h(3))
      do ic = 1, 3
         call equilibrium_column(ncl(ic), 1.0d6, 4000, err(ic), dr_h(ic), &
                                 mass0, mass1)
      enddo
      ord21 = log(err(1)/err(2))/log(dr_h(1)/dr_h(2))
      ord32 = log(err(2)/err(3))/log(dr_h(2)/dr_h(3))
      write(*,'(A,F6.3,A,F6.3)') '       order(100->200) = ', ord21,      &
                                 '   order(200->400) = ', ord32
      call verdict('T10 ', ord32 .gt. 1.7d0, ord32, 1.7d0,                &
                   'observed order in dr: ')
      ! time-step convergence: the steady state of T1a is dt-independent
      call equilibrium_column(200, 1.0d5, 40000, err(1), dr_h(1),         &
                              mass0, mass1)
      call equilibrium_column(200, 1.0d7, 400, err(2), dr_h(2),           &
                              mass0, mass1)
      write(*,'(A,ES10.3,A,ES10.3)') '       error at dt = 1e5: ',        &
           err(1), '   at dt = 1e7: ', err(2)
      call verdict('T10 ', abs(err(1)-err(2))/err(1) .lt. 1.0d-3,         &
                   abs(err(1)-err(2))/err(1), 1.0d-3,                     &
                   'dt-independence of the steady state: ')

      end subroutine test_T1a_and_T10

      ! ----------------------------------------------------------------- !

      subroutine equilibrium_column(Ncell, dt_use, nstep, err, dr_h,      &
                                    mass0, mass1)
      ! Relax a closed, isothermal, motionless neutral column to diffusive
      ! equilibrium and compare with the analytic barometric separation
      !     dx/dr = -x(1-x) (m_He - m_1) g/(kT)
      ! whose isothermal integral with g = b0 v0^2/(R0 r^2), kT = mu v0^2 and
      ! m_He - m_1 = 3 m_H is
      !     logit x(r) - logit x(r_ref) = -3 b0 (1/r_ref - 1/r).
      integer, intent(in)  :: Ncell, nstep
      real*8,  intent(in)  :: dt_use
      real*8,  intent(out) :: err, dr_h, mass0, mass1

      real*8, allocatable :: heh_l(:), xnum(:)
      real*8  :: lg_num, lg_ana, lg_ref, r_ref
      integer :: j, it, jref

      call setup_column(Ncell, 2.0d0, 1.0d3, 1.0d10, 1.33d0, 1.0d10)
      allocate(heh_l(1-Ng:N+Ng), xnum(1-Ng:N+Ng))
      heh_l = HeH_default()
      HeH   = HeH_default()
      call set_composition(f_a, heh_l, 0.0d0, 0)
      rho_a = 1.0d0
      v_a   = 0.0d0
      T_a   = 1.0d0
      dt_a  = dt_use
      dr_h  = r(2) - r(1)
      mass0 = helium_mass_in_column(f_a, 1)

      do it = 1, nstep
         call element_diffusion_step(rho_a, v_a, T_a, f_a, dt_a,          &
                                     closed_base = .true.)
      enddo
      mass1 = helium_mass_in_column(f_a, 1)

      xnum = nucleus_fraction_He(f_a)
      jref = N/2
      r_ref = r(jref)
      lg_ref = log(xnum(jref)/(1.0d0-xnum(jref)))
      err = 0.0d0
      do j = 1, N
         lg_num = log(xnum(j)/(1.0d0-xnum(j))) - lg_ref
         lg_ana = -3.0d0*b0*(1.0d0/r_ref - 1.0d0/r(j))
         err = max(err, abs(lg_num - lg_ana))
      enddo
      deallocate(heh_l, xnum)

      end subroutine equilibrium_column

      ! ----------------------------------------------------------------- !

      real*8 function HeH_default()
      HeH_default = 0.0833333333333333d0
      end function HeH_default

      ! ================================================================= !
      !  T1b  Dirichlet base: the boundary flux accounts for the mass change
      ! ================================================================= !

      subroutine test_T1b()
      real*8, allocatable :: heh_l(:)
      real*8  :: mass0, mass1, budget, dt_use, rel
      integer :: it

      write(*,'(A)') ' --- T1b: Dirichlet reservoir base, integrated '//   &
                     'boundary flux vs mass change'
      call setup_column(200, 2.0d0, 1.0d3, 1.0d10, 1.33d0, 1.0d10)
      allocate(heh_l(1-Ng:N+Ng))
      heh_l = HeH_default()
      HeH   = HeH_default()
      call set_composition(f_a, heh_l, 0.0d0, 0)
      rho_a = 1.0d0
      v_a   = 0.0d0
      T_a   = 1.0d0
      dt_use = 1.0d4
      dt_a  = dt_use

      mass0 = helium_mass_in_column(f_a, 2)
      budget = 0.0d0
      do it = 1, 500
         call element_diffusion_step(rho_a, v_a, T_a, f_a, dt_a,          &
                                     Jface_out = Jf_a)
         budget = budget + dt_use*t_s                                     &
                  *((r_edg(1)*R0)**2*Jf_a(1) - (r_edg(N)*R0)**2*Jf_a(N))
      enddo
      mass1 = helium_mass_in_column(f_a, 2)
      rel = abs((mass1 - mass0) - budget)/max(abs(mass1-mass0), 1.0d-300)
      write(*,'(A,ES16.9,A,ES16.9)') '       mass change = ',             &
           mass1-mass0, '   integrated base flux = ', budget
      call verdict('T1b ', rel .lt. 1.0d-10, rel, 1.0d-10,                &
                   'relative budget mismatch: ')
      deallocate(heh_l)

      end subroutine test_T1b

      ! ================================================================= !
      !  T3  hydrogen-trace limit
      ! ================================================================= !

      subroutine test_T3()
      real*8, allocatable :: heh_l(:), xnum(:), nH_l(:), nHe_l(:)
      real*8  :: Xmin, Xmax, hmax, X_base, m_1, Xtop
      integer :: it, j
      logical :: finite_ok

      write(*,'(A)') ' --- T3: hydrogen-trace limit, He/H = 1000 with an '//&
                     'outflow'
      call setup_column(200, 3.0d0, 5.0d3, 1.0d10, 5.0d0, 1.0d9)
      allocate(heh_l(1-Ng:N+Ng), xnum(1-Ng:N+Ng))
      allocate(nH_l(1-Ng:N+Ng), nHe_l(1-Ng:N+Ng))
      heh_l = 1.0d3
      HeH   = 1.0d3
      call set_composition(f_a, heh_l, 0.0d0, 0)
      do j = 1-Ng, N+Ng
         rho_a(j) = exp(-2.0d0*(r(j)-1.0d0))
         v_a(j)   = 0.05d0*r(j)                    ! outflowing wind
      enddo
      T_a  = 1.0d0
      dt_a = 1.0d3

      do it = 1, 3000
         call element_diffusion_step(rho_a, v_a, T_a, f_a, dt_a)
      enddo

      call nucleus_counts(f_a, nH_l, nHe_l)
      m_1 = mass_per_H_nucleus_without_He()
      X_base = 4.0d0*HeH/(m_1 + 4.0d0*HeH)
      Xmin = 1.0d300
      Xmax = -1.0d300
      hmax = 0.0d0
      finite_ok = .true.
      do j = 1, N
         Xtop = 4.0d0*nHe_l(j)/(m_1*nH_l(j) + 4.0d0*nHe_l(j))
         Xmin = min(Xmin, Xtop)
         Xmax = max(Xmax, Xtop)
         hmax = max(hmax, nHe_l(j)/max(nH_l(j),1.0d-300))
         if (.not. (Xtop .eq. Xtop)) finite_ok = .false.
         if (abs(Xtop) .gt. 1.0d300)  finite_ok = .false.
      enddo
      Xtop = 4.0d0*nHe_l(N)/(m_1*nH_l(N) + 4.0d0*nHe_l(N))
      write(*,'(A,ES12.5,A,ES12.5,A,ES12.5)') '       X range [',         &
           Xmin, ',', Xmax, ']   X_base = ', X_base
      write(*,'(A,ES12.5,A,ES12.5)') '       max He/H = ', hmax,          &
           '   (He/H)/HeH at the top = ',                                 &
           nHe_l(N)/max(nH_l(N),1.0d-300)/HeH
      call verdict('T3  ', finite_ok .and. Xmin .ge. 0.0d0 .and.          &
                   Xmax .le. 1.0d0, Xmax, 1.0d0,                          &
                   'finite and 0 <= X <= 1, max X: ')
      call verdict('T3  ', Xtop .lt. X_base, X_base - Xtop, 0.0d0,        &
                   'hydrogen rises, X_base - X_top: ')

      deallocate(heh_l, xnum, nH_l, nHe_l)

      end subroutine test_T3

      ! ================================================================= !
      !  T4  elemental conservation in a closed column
      ! ================================================================= !

      subroutine test_T4()
      write(*,'(A)') ' --- T4: closed-column elemental conservation, '//   &
                     '1e4 steps, metals present'
      call conserving_column(.true.)
      call conserving_column(.false.)
      end subroutine test_T4

      subroutine conserving_column(eos_met)
      logical, intent(in) :: eos_met
      real*8, allocatable :: heh_l(:)
      real*8  :: mass0, mass1, rel
      integer :: it, j
      character(len=3) :: tag

      call setup_column(200, 2.0d0, 1.0d3, 1.0d10, 1.33d0, 1.0d10)
      thereis_metals     = .true.
      eos_include_metals = eos_met
      melem_ab = 1.0d-4
      allocate(heh_l(1-Ng:N+Ng))
      HeH = HeH_default()
      do j = 1-Ng, N+Ng                       ! arbitrary initial profile
         heh_l(j) = HeH*(1.0d0 + 0.5d0*sin(6.0d0*(r(j)-1.0d0)))
      enddo
      call set_composition(f_a, heh_l, 0.3d0, 1)
      do j = 1-Ng, N+Ng
         rho_a(j) = exp(-1.5d0*(r(j)-1.0d0))
      enddo
      v_a  = 0.0d0
      T_a  = 1.0d0
      ! dt at about the cell diffusion time: the conditioning of the implicit
      ! step (and with it the round-off of the solve) grows with dt/t_diff,
      ! see the table under T1a.
      dt_a = 1.0d1

      mass0 = helium_mass_in_column(f_a, 1)
      do it = 1, 10000
         call element_diffusion_step(rho_a, v_a, T_a, f_a, dt_a,          &
                                     closed_base = .true.)
      enddo
      mass1 = helium_mass_in_column(f_a, 1)
      rel = abs(mass1-mass0)/mass0
      tag = 'on '
      if (.not. eos_met) tag = 'off'
      call verdict('T4  ', rel .lt. 1.0d-12, rel, 1.0d-12,                &
                   'He mass drift, eos_metals '//tag//': ')
      deallocate(heh_l)

      end subroutine conserving_column

      ! ================================================================= !
      !  T5  two-component mass closure of the write-back
      ! ================================================================= !

      subroutine test_T5()
      real*8, allocatable :: heh_l(:), nH_l(:), nHe_l(:), msum0(:)
      real*8  :: resid, m_1
      integer :: it, j

      write(*,'(A)') ' --- T5: single diffusive flux -- the two-component'//&
                     ' mass closure of the write-back'
      write(*,'(A)') '       (a binary mixture has ONE independent '//     &
                     'diffusive flux: the operator carries J_He and the'
      write(*,'(A)') '        hydrogen component -J_He, so J_He + J_1 = '//&
                     '0 identically; what the discretization'
      write(*,'(A)') '        can break is m_1 n_H + m_He n_He = rho, '//  &
                     'which is what is measured here)'
      call setup_column(200, 2.0d0, 1.0d3, 1.0d10, 1.33d0, 1.0d10)
      thereis_metals     = .true.
      eos_include_metals = .true.
      melem_ab = 1.0d-4
      allocate(heh_l(1-Ng:N+Ng), nH_l(1-Ng:N+Ng), nHe_l(1-Ng:N+Ng))
      allocate(msum0(1-Ng:N+Ng))
      HeH = HeH_default()
      do j = 1-Ng, N+Ng
         heh_l(j) = HeH*(1.0d0 + 0.5d0*sin(6.0d0*(r(j)-1.0d0)))
         rho_a(j) = exp(-1.5d0*(r(j)-1.0d0))
         v_a(j)   = 0.02d0*r(j)
      enddo
      call set_composition(f_a, heh_l, 0.3d0, 1)
      T_a  = 1.0d0
      dt_a = 1.0d3
      m_1  = mass_per_H_nucleus_without_He()
      call nucleus_counts(f_a, nH_l, nHe_l)
      msum0 = m_1*nH_l + 4.0d0*nHe_l

      resid = 0.0d0
      do it = 1, 200
         call element_diffusion_step(rho_a, v_a, T_a, f_a, dt_a)
         call nucleus_counts(f_a, nH_l, nHe_l)
         do j = 1, N
            resid = max(resid, abs(m_1*nH_l(j) + 4.0d0*nHe_l(j)           &
                                   - msum0(j))/msum0(j))
         enddo
      enddo
      call verdict('T5  ', resid .lt. 1.0d-13, resid, 1.0d-13,            &
                   'max relative mass-closure error: ')
      deallocate(heh_l, nH_l, nHe_l, msum0)

      end subroutine test_T5

      ! ================================================================= !
      !  T6  uniform mixture preserved under an arbitrary velocity field
      ! ================================================================= !

      subroutine test_T6()
      real*8, allocatable :: heh_l(:), nH_l(:), nHe_l(:)
      real*8  :: X_uni, dev, m_1, Xj
      integer :: it, j

      write(*,'(A)') ' --- T6: uniform X preserved for g = 0 and an '//    &
                     'arbitrary v(r)'
      call setup_column(200, 2.0d0, 1.0d3, 1.0d10, 0.0d0, 1.0d10)
      allocate(heh_l(1-Ng:N+Ng), nH_l(1-Ng:N+Ng), nHe_l(1-Ng:N+Ng))
      HeH   = HeH_default()
      heh_l = HeH
      call set_composition(f_a, heh_l, 0.0d0, 0)
      do j = 1-Ng, N+Ng
         rho_a(j) = exp(-1.5d0*(r(j)-1.0d0))
         ! compression and expansion in the same column
         v_a(j)   = 0.3d0*sin(5.0d0*(r(j)-1.0d0))
      enddo
      T_a  = 1.0d0 + 0.5d0*(r - 1.0d0)          ! non-uniform temperature
      dt_a = 1.0d2

      m_1 = mass_per_H_nucleus_without_He()
      X_uni = 4.0d0*HeH/(m_1 + 4.0d0*HeH)
      do it = 1, 500
         call element_diffusion_step(rho_a, v_a, T_a, f_a, dt_a)
      enddo
      call nucleus_counts(f_a, nH_l, nHe_l)
      dev = 0.0d0
      do j = 1, N
         Xj = 4.0d0*nHe_l(j)/(m_1*nH_l(j) + 4.0d0*nHe_l(j))
         dev = max(dev, abs(Xj - X_uni)/X_uni)
      enddo
      call verdict('T6  ', dev .lt. 1.0d-13, dev, 1.0d-13,                &
                   'max relative departure from uniform: ')
      deallocate(heh_l, nH_l, nHe_l)

      end subroutine test_T6

      ! ================================================================= !
      !  T11  divergence-point cell stays coupled to its donor
      ! ================================================================= !

      subroutine test_T11()
      ! Regression for the base helium hole (docs/binary_diffusion_design.md
      ! section 9, post-review findings).  A breathing base alternates the sign
      ! of v from cell to cell, so every face average (v_j + v_{j+1})/2 is
      ! ~0 while the cell velocities are not.  An advection built from face
      ! velocities then decouples such a cell from BOTH neighbours and any
      ! composition it holds is frozen at the molecular diffusion time; built
      ! from the cell velocity it stays coupled to its donor.  The column below
      ! is seeded with a 500-fold helium hole in the first free cell -- the
      ! measured HD 209458 b defect -- and must refill it.
      real*8, allocatable :: heh_l(:), nH_l(:), nHe_l(:)
      real*8  :: m_1, X_base, X2, vamp, dr_phys, dev
      integer :: it, j

      write(*,'(A)') ' --- T11: cell-to-cell alternating v -- the first free'//&
                     ' cell refills'
      ! b0 = 0: no settling, so the steady state of the column is exactly the
      ! reservoir composition and the test measures the transport alone.
      call setup_column(200, 2.0d0, 1.0d3, 1.0d10, 0.0d0, 1.0d10)
      allocate(heh_l(1-Ng:N+Ng), nH_l(1-Ng:N+Ng), nHe_l(1-Ng:N+Ng))
      HeH   = HeH_default()
      heh_l = HeH
      heh_l(2) = HeH/5.0d2                        ! the measured base hole
      call set_composition(f_a, heh_l, 0.0d0, 0)
      rho_a = 1.0d0
      T_a   = 1.0d0
      ! v alternating in sign cell by cell: every face average vanishes, the
      ! cell velocities do not.  v(2) > 0, so the donor of the first free cell
      ! is the Dirichlet base.
      vamp = 1.0d-2
      do j = 1-Ng, N+Ng
         v_a(j) = vamp*dble((-1)**j)
      enddo
      dr_phys = (r_edg(2) - r_edg(1))*R0
      dt_a    = 0.5d0*dr_phys/(vamp*v0)/t_s        ! half an advective crossing

      m_1    = mass_per_H_nucleus_without_He()
      X_base = 4.0d0*HeH/(m_1 + 4.0d0*HeH)
      do it = 1, 200
         call element_diffusion_step(rho_a, v_a, T_a, f_a, dt_a)
      enddo
      call nucleus_counts(f_a, nH_l, nHe_l)
      X2  = 4.0d0*nHe_l(2)/(m_1*nH_l(2) + 4.0d0*nHe_l(2))
      dev = abs(X2 - X_base)/X_base
      write(*,'(A,ES12.5,A,ES12.5)') '       X(first free cell) = ', X2,   &
           '   X_base = ', X_base
      call verdict('T11 ', dev .lt. 1.0d-3, dev, 1.0d-3,                   &
                   'relative departure from the reservoir: ')
      deallocate(heh_l, nH_l, nHe_l)

      end subroutine test_T11

      ! ================================================================= !
      !  T9  ambipolar limits
      ! ================================================================= !

      subroutine test_T9()
      real*8 :: dm

      write(*,'(A)') ' --- T9: ambipolar limits of the relative settling '//&
                     'mass'
      call settling_limit(0, 1.0d0,        3.0d0,          dm, 'neutral')
      call settling_limit(1, 0.5d0,        2.5d0,          dm, 'H+ plasma')
      call settling_limit(2, 4.0d0/3.0d0,  5.0d0/3.0d0,    dm, 'He++ plasma')
      call settling_front()

      end subroutine test_T9

      ! ----------------------------------------------------------------- !

      subroutine settling_limit(mode, mbar, target_dm, dmout, label)
      ! Build a DISCRETELY hydrostatic isothermal column: ln(rho) is chosen so
      ! that the operator's own centred log-derivative of n_e T reproduces
      ! -mbar m_H g/(kT) exactly, i.e. eE = mbar m_H g to round-off.  The
      ! relative settling mass must then be 3 - (Zbar_He - Zbar_1) mbar.
      ! mode 0 neutral (Zbar = 0 for both), 1 H+ plasma with trace He++,
      ! 2 He++ plasma with trace H+.
      integer,          intent(in)  :: mode
      real*8,           intent(in)  :: mbar, target_dm
      real*8,           intent(out) :: dmout
      character(len=*), intent(in)  :: label

      real*8, allocatable :: heh_l(:), dm_l(:)
      real*8  :: dr_u, err
      integer :: j

      call setup_column(200, 2.0d0, 1.0d3, 1.0d10, 1.33d0, 1.0d10)
      he_ambipolar = .true.
      allocate(heh_l(1-Ng:N+Ng), dm_l(1-Ng:N+Ng))
      if (mode .eq. 2) then
         heh_l = 1.0d3                          ! helium-dominated
      else
         heh_l = HeH_default()
      endif
      HeH = heh_l(1)
      if (mode .eq. 0) call set_composition(f_a, heh_l, 0.0d0, 0)
      if (mode .eq. 1) call set_composition(f_a, heh_l, 1.0d0, 2)
      if (mode .eq. 2) call set_composition(f_a, heh_l, 1.0d0, 2)

      dr_u = r(2) - r(1)
      rho_a(1-Ng) = 1.0d0
      rho_a(2-Ng) = 1.0d0
      do j = 2-Ng, N+Ng-1
         rho_a(j+1) = rho_a(j-1)                                          &
                      *exp(-2.0d0*dr_u*mbar*b0/r(j)**2)
      enddo
      T_a = 1.0d0

      dm_l = relative_settling_mass(rho_a, T_a, f_a)
      err = 0.0d0
      do j = 1, N
         err = max(err, abs(dm_l(j) - target_dm))
      enddo
      dmout = dm_l(N/2)
      write(*,'(A,A,A,F10.7,A,F10.7)') '       ', label,                  &
           ': dm_eff = ', dmout, '   target ', target_dm
      call verdict('T9  ', err .lt. 1.0d-12, err, 1.0d-12,                &
                   'max |dm_eff - target| ('//label//'): ')
      deallocate(heh_l, dm_l)

      end subroutine settling_limit

      ! ----------------------------------------------------------------- !

      subroutine settling_front()
      ! Partially ionized front: the ionization fraction sweeps from 0 to 1
      ! across the column, over a width comparable with the pressure scale
      ! length (a front much thinner than that carries a correspondingly large
      ! electron-pressure gradient, and with it a large ambipolar field -- that
      ! is physics, not a defect, so the test asks for a resolved front).
      ! dm_eff must stay finite and vary smoothly.
      real*8, allocatable :: heh_l(:), dm_l(:)
      real*8  :: xi, mpH, jump, m_1, dr_u
      integer :: j
      logical :: finite_ok

      call setup_column(200, 2.0d0, 1.0d3, 1.0d10, 1.33d0, 1.0d10)
      he_ambipolar = .true.
      allocate(heh_l(1-Ng:N+Ng), dm_l(1-Ng:N+Ng))
      HeH   = HeH_default()
      heh_l = HeH
      m_1 = mass_per_H_nucleus_without_He()
      f_a = 0.0d0
      do j = 1-Ng, N+Ng
         xi  = 0.5d0*(1.0d0 + tanh((r(j)-1.5d0)/0.3d0))
         mpH = m_1 + 4.0d0*HeH
         f_a(j,isp_HI)    = (1.0d0-xi)/mpH
         f_a(j,isp_HII)   = xi/mpH
         f_a(j,isp_HeI)   = HeH*(1.0d0-xi)/mpH
         f_a(j,isp_HeIII) = HeH*xi/mpH
      enddo
      dr_u = r(2) - r(1)
      rho_a(1-Ng) = 1.0d0
      rho_a(2-Ng) = 1.0d0
      do j = 2-Ng, N+Ng-1
         rho_a(j+1) = rho_a(j-1)*exp(-2.0d0*dr_u*0.7d0*b0/r(j)**2)
      enddo
      T_a = 1.0d0

      dm_l = relative_settling_mass(rho_a, T_a, f_a)
      finite_ok = .true.
      jump = 0.0d0
      do j = 2, N-1
         if (.not. (dm_l(j) .eq. dm_l(j))) finite_ok = .false.
         if (abs(dm_l(j)) .gt. 1.0d3)      finite_ok = .false.
         jump = max(jump, abs(dm_l(j+1) - dm_l(j)))
      enddo
      write(*,'(A,F10.7,A,F10.7,A,ES10.3)') '       front: dm_eff '//     &
           'neutral side ', dm_l(2), '  ionized side ', dm_l(N-1),        &
           '  max cell-to-cell jump ', jump
      call verdict('T9  ', finite_ok .and. jump .lt. 1.0d-1, jump,        &
                   1.0d-1, 'finite and continuous across the front: ')
      deallocate(heh_l, dm_l)

      end subroutine settling_front

      ! ================================================================= !
      !  T12  stage-resolved friction (memo 2.6)
      ! ================================================================= !

      subroutine test_T12()
      ! The effective binary coefficient D_eff(He, component 1) in three
      ! prescribed ionization states of one isothermal column, against
      ! closed-form values written out HERE from the equations of
      ! docs/binary_diffusion_design.md section 2.6 -- the module is not
      ! asked for the pair coefficients, only for the mixture, so the check
      ! is not circular.
      !
      !   (a) all neutral      -> the Banks & Kockarts hard sphere, eq. (8)
      !   (b) fully ionized    -> the He++/H+ Coulomb coefficient, eq. (10)
      !   (c) half ionized     -> the four-term Blanc average, eq. (11), of
      !                           one hard-sphere pair, one Coulomb pair and
      !                           two ion-neutral pairs, each ion-neutral
      !                           coefficient being the polarization and
      !                           rigid-core channels added as frictions,
      !                           1/D_in = 1/D_pol + 1/D_hs
      real*8, dimension(:), allocatable :: heh_c, Deff, Dneut
      real*8  :: TK, ntot, ne_l, kT, e2, ee, lamD, mHHe, err
      real*8  :: Dhs, Dcoul, Dph, Dphe, Dmix, lnL11, lnL12, lnL22
      real*8  :: Dinh, Dinhe, Dlo, Dpl, Dhl, Dcmb
      real*8  :: Dc11, Dc12, mpH
      integer :: j

      write(*,'(A)') ' --- T12: stage-resolved friction, three limits'

      call setup_column(64, 2.0d0, 1.0d4, 1.0d10, 5.0d0, 1.0d8)
      allocate(heh_c(1-Ng:N+Ng), Deff(1-Ng:N+Ng), Dneut(1-Ng:N+Ng))
      heh_c = 0.0833333333333333d0

      ee   = 4.803204713d-10
      e2   = ee*ee
      TK   = T0
      kT   = kb_erg*TK
      mHHe = 4.0d0/5.0d0*mu                       ! reduced mass of He-H [g]

      ! ---------------- (a) all neutral: HI + HeI ---------------------- !
      call set_composition(f_a, heh_c, 0.0d0, 0)
      call helium_hydrogen_diffusion(rho_a, T_a, f_a, Deff, Dneut)
      call column_totals(f_a, ntot, ne_l)
      Dhs = 1.52d18*sqrt(1.0d0/4.0d0 + 1.0d0/1.0d0)*sqrt(TK)/ntot
      err = abs(Deff(N/2) - Dhs)/Dhs
      call verdict('T12a', err .lt. 1.0d-12, err, 1.0d-12,                &
           'all neutral: D_eff vs Banks & Kockarts, rel. err = ')

      ! ---------------- (b) fully ionized: HII + HeIII ----------------- !
      call set_composition(f_a, heh_c, 1.0d0, 2)
      call helium_hydrogen_diffusion(rho_a, T_a, f_a, Deff, Dneut)
      call column_totals(f_a, ntot, ne_l)
      lamD  = sqrt(kT/(4.0d0*pi*ne_l*e2))
      lnL12 = max(log(3.0d0*kT*lamD/(2.0d0*e2)), 1.0d0)
      Dcoul = 3.0d0*kT**2.5d0                                             &
              /(4.0d0*sqrt(2.0d0*pi*mHHe)*ntot*(2.0d0*e2)**2*lnL12)
      err = abs(Deff(N/2) - Dcoul)/Dcoul
      call verdict('T12b', err .lt. 1.0d-12, err, 1.0d-12,                &
           'fully ionized: D_eff vs Coulomb He++/H+, rel. err = ')
      write(*,'(A,ES12.5,A,ES12.5,A,ES10.3)')                             &
           '       D_neutral = ', Dneut(N/2), '  D_eff = ', Deff(N/2),    &
           '  ratio = ', Deff(N/2)/Dneut(N/2)

      ! ---------------- (c) half ionized: 50/50 in each element -------- !
      f_a = 0.0d0
      do j = 1-Ng, N+Ng
         mpH = mass_per_H_nucleus_without_He() + 4.0d0*heh_c(j)
         f_a(j,isp_HI)    = 0.5d0/mpH
         f_a(j,isp_HII)   = 0.5d0/mpH
         f_a(j,isp_HeI)   = 0.5d0*heh_c(j)/mpH
         f_a(j,isp_HeII)  = 0.5d0*heh_c(j)/mpH
      enddo
      call helium_hydrogen_diffusion(rho_a, T_a, f_a, Deff, Dneut)
      call column_totals(f_a, ntot, ne_l)
      lamD  = sqrt(kT/(4.0d0*pi*ne_l*e2))
      lnL11 = max(log(3.0d0*kT*lamD/e2), 1.0d0)
      Dhs   = 1.52d18*sqrt(1.0d0/4.0d0 + 1.0d0/1.0d0)*sqrt(TK)/ntot
      Dc11  = 3.0d0*kT**2.5d0                                             &
              /(4.0d0*sqrt(2.0d0*pi*mHHe)*ntot*e2**2*lnL11)
      ! He+ against neutral H: polarization on alpha(H); HeI against H+:
      ! polarization on alpha(He).  Same reduced mass in both.  Each
      ! ion-neutral pair adds the rigid-core channel to the polarization one,
      ! 1/D_in = 1/D_pol + 1/D_hs, so its friction is the sum of the two.
      Dph   = kT/(2.21d0*pi*ee*ntot*sqrt(alpha_HI *mHHe))
      Dphe  = kT/(2.21d0*pi*ee*ntot*sqrt(alpha_HeI*mHHe))
      Dinh  = 1.0d0/(1.0d0/Dph  + 1.0d0/Dhs)
      Dinhe = 1.0d0/(1.0d0/Dphe + 1.0d0/Dhs)
      Dmix  = 1.0d0/(0.25d0*(1.0d0/Dhs + 1.0d0/Dinhe                      &
                           + 1.0d0/Dinh + 1.0d0/Dc11))
      err = abs(Deff(N/2) - Dmix)/Dmix
      call verdict('T12c', err .lt. 1.0d-12, err, 1.0d-12,                &
           'half ionized: D_eff vs 4-term Blanc average, rel. err = ')
      write(*,'(A,4ES12.4)') '       pair D (hs, in-He, in-H, Coul):  ',  &
           Dhs, Dinhe, Dinh, Dc11
      write(*,'(A,2ES12.4,A,2ES12.4)')                                    &
           '       polarization channel alone (He, H): ', Dphe, Dph,      &
           '   combined/polarization: ', Dinhe/Dphe, Dinh/Dph

      ! ---------------- (d) the two limits of the ion-neutral form ------ !
      ! 1/D_in = 1/D_pol + 1/D_hs must reduce to D_pol where the
      ! polarization friction dominates and to D_hs where the rigid core
      ! does, and must never exceed either channel at any temperature.  The
      ! approach to each limit is slow -- D_pol/D_hs ~ T^(1/2), crossing
      ! unity near 1.5e3 K for this pair -- so the limits are taken far from
      ! the crossover and the tolerance is the expected residual 1/(1+ratio).
      do j = 1, 2
         if (j .eq. 1) Dlo = 1.0d0                  ! K, polarization limit
         if (j .eq. 2) Dlo = 1.0d9                  ! K, rigid-core limit
         Dpl   = kb_erg*Dlo/(2.21d0*pi*ee*ntot*sqrt(alpha_HI*mHHe))
         Dhl   = 1.52d18*sqrt(1.25d0)*sqrt(Dlo)/ntot
         Dcmb  = 1.0d0/(1.0d0/Dpl + 1.0d0/Dhl)
         if (j .eq. 1) then
            err = abs(Dcmb - Dpl)/Dpl
            call verdict('T12d', err .lt. 3.0d-2 .and. Dcmb .le. Dpl      &
                         .and. Dcmb .le. Dhl, err, 3.0d-2,                &
                 'T = 1 K: combined -> polarization, rel. departure = ')
         else
            err = abs(Dcmb - Dhl)/Dhl
            call verdict('T12d', err .lt. 3.0d-3 .and. Dcmb .le. Dpl      &
                         .and. Dcmb .le. Dhl, err, 3.0d-3,                &
                 'T = 1e9 K: combined -> hard sphere, rel. departure = ')
         endif
      enddo
      write(*,'(A,F8.1,A)') '       polarization/rigid-core crossover at ',&
           (2.21d0*pi*ee*sqrt(alpha_HI*mHHe)                              &
            *1.52d18*sqrt(1.25d0)/kb_erg)**2, ' K'

      deallocate(heh_c, Deff, Dneut)

      end subroutine test_T12

      ! ----------------------------------------------------------------- !

      subroutine column_totals(f_sp, ntot, ne_l)
      ! Total nucleus density and free-electron density of the (uniform)
      ! column, in cm^-3, under the same conventions the operator uses:
      ! nuclei for n, every charge for n_e.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8,                                 intent(out) :: ntot, ne_l
      real*8, dimension(1-Ng:N+Ng) :: nH_l, nHe_l
      integer :: ib
      call nucleus_counts(f_sp, nH_l, nHe_l)
      ntot = (nH_l(N/2) + nHe_l(N/2))*rho_a(N/2)*n0
      ne_l = 0.0d0
      do ib = 1, n_bsp
         if (bsp_is_excited_level(ib)) cycle
         ne_l = ne_l + dble(bsp_charge(ib))*f_sp(N/2,bsp_fsp(ib))
      enddo
      ne_l = max(ne_l*rho_a(N/2)*n0, 1.0d0)
      end subroutine column_totals

      ! ================================================================= !
      !  T7  molecular closure and homopause
      ! ================================================================= !

      subroutine test_T7()
      ! The closure of docs/binary_diffusion_design.md section 5 (milestone
      ! M4), in three parts:
      !
      !  (a) forced ATOMIC composition   -> D_eff = D(He,H)
      !  (b) forced MOLECULAR composition (every H nucleus in H2)
      !                                  -> D_eff = D(He,H2)
      !  (c) a prescribed molecular column with a constant K_zz: the helium
      !      profile is flat where K_zz >> D_eff and separates above it, and
      !      the homopause radius follows D_eff(r_h) = K_zz.
      !
      ! Every comparison value is typed out here from the equations of the
      ! memo, not read back from the module: the hard-sphere pair coefficients
      ! (8) with the CARRIER mass and the CARRIER density, the Blanc mixture
      ! (11) over the two hydrogen carriers, and the settling coefficient with
      ! the mean carrier mass of section 5.
      real*8, dimension(:), allocatable :: qmol, Deff, Dneut, Dt, Gt, Sf
      real*8  :: ntot, Dana, err, mpH, TK, m_1, kk(3), rt(3)
      real*8  :: rh_num, rh_ana, dr_cell, flat, sep
      integer :: j, ik

      write(*,'(A)') ' --- T7: molecular closure (Blanc carriers) and'//    &
                     ' the homopause'

      ! ---------------- (a) forced atomic ------------------------------ !
      call setup_column(64, 2.0d0, 1.0d3, 1.0d10, 1.0d0, 1.0d12)
      thereis_mol = .true.
      HeH = HeH_default()
      allocate(qmol(1-Ng:N+Ng), Deff(1-Ng:N+Ng), Dneut(1-Ng:N+Ng))
      TK  = T0
      m_1 = mass_per_H_nucleus_without_He()

      qmol = 0.0d0
      call set_molecular_composition(f_a, HeH, qmol)
      call helium_hydrogen_diffusion(rho_a, T_a, f_a, Deff, Dneut)
      ! carriers = nuclei in the atomic limit
      mpH  = m_1 + 4.0d0*HeH
      ntot = (1.0d0 + HeH)/mpH*rho_a(N/2)*n0
      Dana = 1.52d18*sqrt(1.0d0/4.0d0 + 1.0d0/1.0d0)*sqrt(TK)/ntot
      err  = abs(Deff(N/2) - Dana)/Dana
      call verdict('T7a ', err .lt. 1.0d-12, err, 1.0d-12,                 &
           'forced atomic: D_eff vs D(He,H), rel. err = ')

      ! ---------------- (b) forced fully molecular --------------------- !
      qmol = 1.0d0
      call set_molecular_composition(f_a, HeH, qmol)
      call helium_hydrogen_diffusion(rho_a, T_a, f_a, Deff, Dneut)
      ! one carrier per two H nuclei, so the carrier density is (1/2 + He/H)
      ntot = (0.5d0 + HeH)/mpH*rho_a(N/2)*n0
      Dana = 1.52d18*sqrt(1.0d0/4.0d0 + 1.0d0/2.0d0)*sqrt(TK)/ntot
      err  = abs(Deff(N/2) - Dana)/Dana
      call verdict('T7b ', err .lt. 1.0d-12, err, 1.0d-12,                 &
           'forced molecular: D_eff vs D(He,H2), rel. err = ')
      write(*,'(A,ES12.5,A,ES12.5)') '       D(He,H) = ',                  &
           1.52d18*sqrt(1.25d0)*sqrt(TK)/((1.0d0+HeH)/mpH*rho_a(N/2)*n0),  &
           '   D(He,H2) = ', Dana
      deallocate(qmol, Deff, Dneut)

      ! ---------------- (c) homopause against K_zz --------------------- !
      ! A motionless isothermal column, neutral throughout, molecular below
      ! r = 1.8 and atomic above, with a density falling steeply enough that
      ! D_eff spans several decades.  With v = 0, a Dirichlet reservoir base and
      ! a zero-flux top the steady state has J = 0 everywhere, i.e.
      !
      !     d logit(X)/dr = - [ D/(D+K) ] G ,
      !
      ! so the measured slope divided by the typed-out G is exactly the
      ! suppression factor D/(D+K), and the radius where it passes 1/2 is the
      ! radius where D = K.  That is the homopause, and the test is whether it
      ! sits where the typed-out D_eff(r) equals the K_zz that was imposed.
      call setup_column(300, 3.0d0, 1.0d3, 1.0d10, 2.0d0, 1.0d14)
      thereis_mol = .true.
      HeH = HeH_default()
      allocate(qmol(1-Ng:N+Ng), Dt(1-Ng:N+Ng), Gt(1-Ng:N+Ng),              &
               Sf(1-Ng:N+Ng))
      do j = 1-Ng, N+Ng
         qmol(j)  = 0.5d0*(1.0d0 - tanh((r(j) - 1.8d0)/0.2d0))
         rho_a(j) = exp(-25.0d0*(1.0d0 - 1.0d0/max(r(j), 0.5d0)))
      enddo
      v_a  = 0.0d0
      T_a  = 1.0d0
      dt_a = 1.0d9

      ! K_zz values placed so that the three homopauses land at 1.5, 2.0, 2.5
      call set_molecular_composition(f_a, HeH, qmol)
      call molecular_column_coefficients(f_a, Dt, Gt)
      rt = [1.5d0, 2.0d0, 2.5d0]
      do ik = 1, 3
         j = 1 + nint((rt(ik) - 1.0d0)/(r(2) - r(1)))
         kk(ik) = Dt(j)
      enddo
      dr_cell = (r(2) - r(1))

      write(*,'(A)') '       K_zz [cm2/s]   r_h(D_eff=K_zz)   r_h(profile)'//&
                     '   (He/H)/HeH at 1.05   at 2.9'
      do ik = 1, 3
         call set_molecular_composition(f_a, HeH, qmol)
         he_kzz = kk(ik)
         call eddy_diffusion_on_grid
         do j = 1, 600
            call element_diffusion_step(rho_a, v_a, T_a, f_a, dt_a)
         enddo
         call molecular_column_coefficients(f_a, Dt, Gt)
         call suppression_profile(f_a, Gt, Sf)
         rh_num = crossing_radius(Sf, 0.5d0)
         rh_ana = crossing_radius(Dt/kk(ik), 1.0d0)
         flat   = helium_ratio_at(f_a, 1.05d0)/HeH
         sep    = helium_ratio_at(f_a, 2.9d0)/HeH
         write(*,'(4X,ES12.4,4F16.5)') kk(ik), rh_ana, rh_num, flat, sep
         err = abs(rh_num - rh_ana)
         call verdict('T7c ', err .lt. 2.0d0*dr_cell, err, 2.0d0*dr_cell,  &
              'homopause radius, |profile - D_eff=K_zz| [R_p] = ')
      enddo
      ! The largest K_zz mixes the column right through the molecular front:
      ! below its homopause the ratio must sit on the reservoir value, while
      ! the top of the domain, above it, must still be separated.
      call set_molecular_composition(f_a, HeH, qmol)
      he_kzz = kk(3)
      call eddy_diffusion_on_grid
      do j = 1, 600
         call element_diffusion_step(rho_a, v_a, T_a, f_a, dt_a)
      enddo
      flat = helium_ratio_at(f_a, 1.05d0)/HeH
      sep  = helium_ratio_at(f_a, 2.9d0)/HeH
      call verdict('T7c ', abs(flat - 1.0d0) .lt. 1.0d-2 .and.             &
                   sep .lt. 0.9d0, abs(flat - 1.0d0), 1.0d-2,              &
           'largest K_zz: mixed below the homopause, |He/H / HeH - 1| = ')
      he_kzz = 0.0d0
      call eddy_diffusion_on_grid
      deallocate(qmol, Dt, Gt, Sf)

      ! ---------------- (d) the chemistry driver on its own ------------ !
      ! With no gravity, no field and no eddy, the only driver left in G is
      ! -dln(psi)/dr, and the steady state it produces is the statement that
      ! the closure rests on: diffusion drives the MOLE fraction of helium
      ! among the collision partners uniform, NOT its mass fraction.  A column
      ! that is molecular below and atomic above therefore ends up with a mass
      ! fraction that varies by the factor psi across the front and a mole
      ! fraction x = n_He/(n_1c + n_He) that is flat.  Nothing in this check
      ! comes from the module's algebra: it is the invariant of the
      ! Chapman-Cowling mole-fraction force.
      call setup_column(300, 3.0d0, 1.0d3, 1.0d10, 0.0d0, 1.0d14)
      thereis_mol = .true.
      HeH = HeH_default()
      allocate(qmol(1-Ng:N+Ng), Sf(1-Ng:N+Ng))
      do j = 1-Ng, N+Ng
         qmol(j)  = 0.5d0*(1.0d0 - tanh((r(j) - 1.8d0)/0.2d0))
         rho_a(j) = exp(-25.0d0*(1.0d0 - 1.0d0/max(r(j), 0.5d0)))
      enddo
      v_a  = 0.0d0
      T_a  = 1.0d0
      dt_a = 1.0d9
      he_kzz = 0.0d0
      call eddy_diffusion_on_grid
      call set_molecular_composition(f_a, HeH, qmol)
      do j = 1, 600
         call element_diffusion_step(rho_a, v_a, T_a, f_a, dt_a)
      enddo
      call carrier_fraction_He(f_a, Sf)
      flat = 0.0d0
      sep  = 0.0d0
      do j = 1, N
         flat = max(flat, abs(Sf(j)/Sf(1) - 1.0d0))
         sep  = max(sep, abs(helium_ratio_at(f_a, r(j))/HeH - 1.0d0))
      enddo
      write(*,'(A,ES10.3,A,ES10.3)')                                       &
           '       g = 0: max|x/x_base - 1| = ', flat,                     &
           '   while max|(He/H)/HeH - 1| = ', sep
      call verdict('T7d ', flat .lt. 1.0d-3 .and. sep .gt. 0.1d0, flat,    &
                   1.0d-3,                                                 &
           'g = 0: the MOLE fraction is what levels out, max departure = ')
      deallocate(qmol, Sf)

      end subroutine test_T7

      ! ----------------------------------------------------------------- !

      subroutine carrier_fraction_He(f_sp, xc)
      ! x = n_He/(n_1c + n_He), helium's fraction of the COLLISION PARTNERS,
      ! for the neutral H/H2/He column of T7 (an H2 is one partner).
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: xc
      xc = f_sp(:,isp_HeI)/max(f_sp(:,isp_HI) + f_sp(:,isp_H2)             &
                               + f_sp(:,isp_HeI), 1.0d-300)
      end subroutine carrier_fraction_He

      ! ----------------------------------------------------------------- !

      subroutine set_molecular_composition(f_sp, heh, qm)
      ! Neutral column with a prescribed molecular fraction qm(r) = the share
      ! of the hydrogen NUCLEI bound into H2, normalized as the code counts
      ! f_sp (sum_s m_s f_s = 1).  One H2 holds two nuclei, hence the 1/2.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(out) :: f_sp
      real*8,                                 intent(in)  :: heh
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: qm
      real*8  :: mpH
      integer :: j
      mpH = mass_per_H_nucleus_without_He() + 4.0d0*heh
      f_sp = 0.0d0
      do j = 1-Ng, N+Ng
         f_sp(j,isp_HI)  = (1.0d0 - qm(j))/mpH
         f_sp(j,isp_H2)  = 0.5d0*qm(j)/mpH
         f_sp(j,isp_HeI) = heh/mpH
      enddo
      end subroutine set_molecular_composition

      ! ----------------------------------------------------------------- !

      subroutine molecular_column_coefficients(f_sp, Dt, Gt)
      ! The two coefficients of the steady state, written out here from the
      ! memo rather than taken from the module:
      !
      !   1/D_eff = y_H/D(He,H) + y_H2/D(He,H2),   D(He,c) the hard sphere (8)
      !             with the CARRIER mass and the CARRIER density
      !   G       = (m_He - m_c1) m_H g/(kT) - dln(psi)/dr
      !
      ! for the neutral column of T7c (no ambipolar field, no ionization).
      ! m_c1 = m_1 mcar is the mean carrier mass and psi = 1/mcar the
      ! collision partners per hydrogen nucleus, mcar = (n_H + 2 n_H2)/
      ! (n_H + n_H2).  With g = b0 v0^2/(R0 r^2) and kT = mu v0^2 the force
      ! part is (4 - m_c1) b0/(R0 r^2) in cm^-1, and the chemistry part is the
      ! same central difference the operator takes.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: Dt, Gt
      real*8, dimension(1-Ng:N+Ng) :: mcar
      real*8  :: nHI, nH2, nHe, ncar, yH1, yH2, DHeH, DHeH2, TKl
      integer :: j
      do j = 1-Ng, N+Ng
         nHI  = f_sp(j,isp_HI) *rho_a(j)*n0
         nH2  = f_sp(j,isp_H2) *rho_a(j)*n0
         nHe  = f_sp(j,isp_HeI)*rho_a(j)*n0
         ncar = max(nHI + nH2 + nHe, 1.0d0)
         yH1  = nHI/max(nHI + nH2, 1.0d-30)
         yH2  = nH2/max(nHI + nH2, 1.0d-30)
         mcar(j) = (nHI + 2.0d0*nH2)/max(nHI + nH2, 1.0d-30)
         TKl  = T_a(j)*T0
         DHeH  = 1.52d18*sqrt(1.0d0/4.0d0 + 1.0d0/1.0d0)*sqrt(TKl)/ncar
         DHeH2 = 1.52d18*sqrt(1.0d0/4.0d0 + 1.0d0/2.0d0)*sqrt(TKl)/ncar
         Dt(j) = 1.0d0/(yH1/DHeH + yH2/DHeH2)
         Gt(j) = (4.0d0 - mcar(j)*mass_per_H_nucleus_without_He())         &
                 *b0/(R0*r(j)*r(j))
      enddo
      do j = 2-Ng, N+Ng-1
         Gt(j) = Gt(j) + (log(mcar(j+1)) - log(mcar(j-1)))                 &
                 /max((r(j+1) - r(j-1))*R0, 1.0d0)
      enddo
      end subroutine molecular_column_coefficients

      ! ----------------------------------------------------------------- !

      subroutine suppression_profile(f_sp, Gt, Sf)
      ! S(r) = -[d logit(X)/dr]/G at each face, which at the zero-flux steady
      ! state of T7c is exactly D_eff/(D_eff + K_zz):
      !
      !     rho (D + K) dX/dr = -rho X(1-X) D G
      !
      ! with X the helium MASS fraction, which is what the operator
      ! transports, and G the settling coefficient of
      ! molecular_column_coefficients (forces plus the chemistry term).
      ! Faces where the profile has no gradient left to measure get S = 0.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in)  :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(in)  :: Gt
      real*8, dimension(1-Ng:N+Ng),           intent(out) :: Sf
      real*8, dimension(1-Ng:N+Ng) :: nH_l, nHe_l, Xl
      real*8  :: slope, Gf
      integer :: j
      call nucleus_counts(f_sp, nH_l, nHe_l)
      Xl = 4.0d0*nHe_l/max(mass_per_H_nucleus_without_He()*nH_l            &
                           + 4.0d0*nHe_l, 1.0d-300)
      Sf = 0.0d0
      do j = 1, N-1
         if (Xl(j) .le. 0.0d0 .or. Xl(j+1) .le. 0.0d0) cycle
         if (Xl(j) .ge. 1.0d0 .or. Xl(j+1) .ge. 1.0d0) cycle
         slope = (log(Xl(j+1)/(1.0d0-Xl(j+1)))                             &
                - log(Xl(j)  /(1.0d0-Xl(j))))/((r(j+1)-r(j))*R0)
         Gf = 0.5d0*(Gt(j) + Gt(j+1))
         if (abs(Gf) .gt. 0.0d0) Sf(j) = -slope/Gf
      enddo
      end subroutine suppression_profile

      ! ----------------------------------------------------------------- !

      real*8 function crossing_radius(a, target_v)
      ! The radius at which the monotonically decreasing (S) or increasing
      ! (D/K) profile a(j) passes target_v, by linear interpolation on the
      ! first crossing found from the base outward.  Returns -1 if there is
      ! none.
      real*8, dimension(1-Ng:N+Ng), intent(in) :: a
      real*8,                       intent(in) :: target_v
      integer :: j
      crossing_radius = -1.0d0
      do j = 2, N-2
         if ((a(j) - target_v)*(a(j+1) - target_v) .le. 0.0d0 .and.        &
             a(j) .ne. a(j+1)) then
            crossing_radius = r(j) + (target_v - a(j))/(a(j+1) - a(j))     &
                              *(r(j+1) - r(j))
            return
         endif
      enddo
      end function crossing_radius

      ! ----------------------------------------------------------------- !

      real*8 function helium_ratio_at(f_sp, rq)
      ! He/H (nuclei) at the cell nearest radius rq.
      real*8, dimension(1-Ng:N+Ng,n_species), intent(in) :: f_sp
      real*8,                                 intent(in) :: rq
      real*8, dimension(1-Ng:N+Ng) :: nH_l, nHe_l
      integer :: j, jq
      call nucleus_counts(f_sp, nH_l, nHe_l)
      jq = 1
      do j = 1, N
         if (abs(r(j) - rq) .lt. abs(r(jq) - rq)) jq = j
      enddo
      helium_ratio_at = nHe_l(jq)/max(nH_l(jq), 1.0d-300)
      end function helium_ratio_at

      ! ================================================================= !
      !  T2a  convergence to the trace kernel (recorded, no longer run)
      ! ================================================================= !

      subroutine test_T2a()
      ! T2a compared the binary operator against the trace-helium kernel it
      ! replaces, as helium is made trace.  That kernel was removed with
      ! milestone M3, so the comparison can no longer be run and its result
      ! stands as the measurement made while both existed
      ! (docs/Update_EXHALE.md section 68): the relative difference of the
      ! one-step increments on the LHS 1140 b wind was
      !
      !     He/H = 8.33e-2 : 2.27625e-03
      !     He/H = 1.00e-2 : 3.00403e-04
      !     He/H = 1.00e-3 : 3.04126e-05
      !     He/H = 1.00e-4 : 3.04503e-06
      !
      ! i.e. first order in the helium fraction and 3.05e-6 at He/H = 1e-4,
      ! against the 1e-3 criterion of docs/binary_diffusion_design.md.
      ! No verdict is counted: nothing is measured here.
      write(*,'(A)') ' --- T2a: convergence to the trace kernel'//        &
                     ' -- SKIPPED (comparison kernel removed at M3)'
      write(*,'(A)') '       recorded M2 result: 3.045E-06 at He/H ='//   &
                     ' 1e-4, first order in He/H'
      write(*,'(A)') '       (docs/Update_EXHALE.md section 68)'

      end subroutine test_T2a

      end program diffusion_tests
