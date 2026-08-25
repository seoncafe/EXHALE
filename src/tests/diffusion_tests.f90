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
      !   T2a  recorded only: the kernel it compared against is gone
      !
      ! T0 (`make check` byte-identical with diffusion off), T2b, T7 and T8
      ! are run outside this driver.
      !
      ! The columns are synthetic: the driver sets the global_parameters
      ! scalars and the radial grid itself, so no input.inp, no ionization
      ! solve and no hydro are involved.

      use global_parameters
      use species_table, only: isp_HI, isp_HII, isp_HeI, isp_HeII,        &
                               isp_HeIII, n_bsp, bsp_fsp, bsp_nH, bsp_nHe,&
                               n_melem, melem_i0, mion_fsp, melem_A
      use composition,   only: mass_per_H_nucleus_without_He
      use binary_element_diffusion, only: element_diffusion_step,         &
                                          relative_settling_mass

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
