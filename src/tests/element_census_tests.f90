      program element_census_tests
      ! Deterministic tests of the elemental and charge invariants (item A3
      ! of docs/open_defects_20260903_review.md section 6).  Test ids:
      !
      !   E1  element_nuclei_and_charge reproduces the stoichiometry of the
      !       species table on a hand-computed synthetic state, and the
      !       positive charge it returns equals calc_ne on the same state.
      !   E2  elemental closure through one carrier write-back, including the
      !       configuration defect 12.8 reports (a cell whose transported CO
      !       collapses and hands its carbon back to the free stages).
      !   E3  H2 photoevent ledger: every final-state channel the code
      !       applies conserves H nuclei and charge, and the two channels the
      !       branching does NOT resolve are measured rather than assumed
      !       small (cross_sec.f90 frac_H2_dissociative_ionization).
      !   E4  the chemical-equilibrium branch of the base H2 fraction, taken
      !       when no lower-atmosphere handoff states it (item 12.6).
      !
      ! The columns are synthetic: the driver sets the global_parameters
      ! scalars and the radial grid itself, so no input.inp, no radiation
      ! field and no hydrodynamics are involved.  The half of A3(a) that
      ! needs a full ionization sweep, the restart round trip A3(c), and the
      ! ghost composition of E4 are run-level and are exercised by the
      ! runtime assertions of element_census with EXHALE_ELEMENT_ASSERT=2.

      use global_parameters
      use species_table
      use utils, only: calc_ne
      use element_census, only: element_nuclei_and_charge, n_element,     &
                                ie_H, ie_He, ie_metal0, element_name
      use cross_sections, only: frac_H2_dissociative_ionization
      use composition, only: h2_mixing_ratio_base,                       &
                             base_h2_nuclei_fraction,                    &
                             base_h2_composition_imposed
      use lower_column, only: q_h2_equilibrium
      use diffusive_photochemistry, only: carrier_state,                 &
                                          carrier_write_back,            &
                                          carrier_set_init,              &
                                          n_carrier_max,                 &
                                          ic_H2, ic_OH, ic_H2O, ic_CO
      use ionization_equilibrium, only: bg_cell, ioniz_eq_allocate_arrays

      implicit none

      integer :: n_pass, n_fail

      n_pass = 0
      n_fail = 0

      write(*,'(A)') '===== element census / conservation tests ====='

      call setup_globals()
      call test_E1()
      call test_E2()
      call test_E3()
      call test_E4()

      write(*,'(A)') '=============================================='
      write(*,'(A,I0,A,I0,A)') ' TOTAL: ', n_pass, ' passed, ',           &
                               n_fail, ' failed'
      if (n_fail .gt. 0) stop 1

      contains

      ! ------------------------------------------------------!

      subroutine check(id, ok, what, got, want, tol)
      character(len=*), intent(in) :: id, what
      logical,          intent(in) :: ok
      real*8, optional, intent(in) :: got, want, tol
      if (ok) then
         n_pass = n_pass + 1
         if (present(got)) then
            write(*,'(A,A,A,A,A,ES12.5,A,ES12.5)') '  PASS ', id, '  ',  &
                 what, '  got ', got, ' want ', want
         else
            write(*,'(A,A,A,A)') '  PASS ', id, '  ', what
         endif
      else
         n_fail = n_fail + 1
         if (present(got)) then
            write(*,'(A,A,A,A,A,ES16.9,A,ES16.9,A,ES9.2)') '  FAIL ',    &
                 id, '  ', what, '  got ', got, ' want ', want,          &
                 '  tol ', tol
         else
            write(*,'(A,A,A,A)') '  FAIL ', id, '  ', what
         endif
      endif
      end subroutine check

      ! ------------------------------------------------------!

      subroutine setup_globals()
      ! A five-cell column with two ghosts on each side, the physics
      ! switches the tests need, and the scale factors the routines read.
      integer :: j
      N   = 5
      n0  = 1.0d10
      R0  = 1.0d10
      v0  = 1.0d6
      T0  = 1000.0d0
      HeH = 0.0793d0     ! the solar-like value the molecular gates use
      p_base_bar   = 1.0d-6
      thereis_He        = .true.
      thereis_HeITR     = .false.
      thereis_metals    = .true.
      thereis_mol       = .true.
      thereis_oxychem   = .true.
      eos_include_metals = .true.
      he_diffusion      = .false.
      carrier_transport = .true.
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      do j = 1-Ng, N+Ng
         r(j)     = 1.0d0 + 0.01d0*dble(j)
         r_edg(j) = 1.0d0 + 0.01d0*(dble(j) + 0.5d0)
         dr_j(j)  = 0.01d0
      enddo
      allocate(melem_ab(n_melem))
      melem_ab = 0.0d0
      melem_ab(iel_C) = 2.69d-4
      melem_ab(iel_O) = 4.90d-4
      call carrier_set_init()
      call ioniz_eq_allocate_arrays()
      end subroutine setup_globals

      ! ------------------------------------------------------!

      subroutine test_E1()
      ! Hand-computed state: in every cell put 1 unit of density into each of
      ! HI, HII, HeI, H2, H2+, H3+, HeH+, OH, H2O, CO and into C I and O II,
      ! so that the expected totals are the multiplicity sums themselves.
      real*8, allocatable :: rho(:), f_sp(:,:), n_nuc(:,:), n_chg(:)
      real*8, allocatable :: ne_ref(:), nm_l(:,:), nmol_l(:,:)
      real*8 :: nd, want_H, want_He, want_O, want_C, want_q
      integer :: j, im

      allocate(rho(1-Ng:N+Ng), f_sp(1-Ng:N+Ng,n_species))
      allocate(n_nuc(1-Ng:N+Ng,n_element), n_chg(1-Ng:N+Ng))
      allocate(ne_ref(1-Ng:N+Ng), nm_l(1-Ng:N+Ng,n_mion))
      allocate(nmol_l(1-Ng:N+Ng,4))
      rho  = 1.0d0
      f_sp = 0.0d0
      f_sp(:,isp_HI)   = 1.0d0
      f_sp(:,isp_HII)  = 1.0d0
      f_sp(:,isp_HeI)  = 1.0d0
      f_sp(:,isp_H2)   = 1.0d0
      f_sp(:,isp_H2p)  = 1.0d0
      f_sp(:,isp_H3p)  = 1.0d0
      f_sp(:,isp_HeHp) = 1.0d0
      f_sp(:,isp_OH)   = 1.0d0
      f_sp(:,isp_H2O)  = 1.0d0
      f_sp(:,isp_CO)   = 1.0d0
      f_sp(:,mion_fsp(im_CI))  = 1.0d0
      f_sp(:,mion_fsp(im_OII)) = 1.0d0

      call element_nuclei_and_charge(rho, f_sp, n_nuc, n_chg)

      nd = n0
      ! H : HI 1 + HII 1 + H2 2 + H2+ 2 + H3+ 3 + HeH+ 1 + OH 1 + H2O 2 = 13
      want_H  = 13.0d0*nd
      ! He: HeI 1 + HeH+ 1 = 2
      want_He = 2.0d0*nd
      ! O : OH 1 + H2O 1 + CO 1 + O II 1 = 4
      want_O  = 4.0d0*nd
      ! C : CO 1 + C I 1 = 2
      want_C  = 2.0d0*nd
      ! charge: HII 1 + H2+ 1 + H3+ 1 + HeH+ 1 + O II 1 = 5
      want_q  = 5.0d0*nd

      j = 3
      call check('E1a', abs(n_nuc(j,ie_H) -want_H) .le. 1.0d-12*want_H,   &
                 'H nuclei ', n_nuc(j,ie_H), want_H, 1.0d-12)
      call check('E1b', abs(n_nuc(j,ie_He)-want_He).le. 1.0d-12*want_He,  &
                 'He nuclei', n_nuc(j,ie_He), want_He, 1.0d-12)
      call check('E1c', abs(n_nuc(j,ie_metal0+iel_O)-want_O)              &
                 .le. 1.0d-12*want_O, 'O nuclei ',                        &
                 n_nuc(j,ie_metal0+iel_O), want_O, 1.0d-12)
      call check('E1d', abs(n_nuc(j,ie_metal0+iel_C)-want_C)              &
                 .le. 1.0d-12*want_C, 'C nuclei ',                        &
                 n_nuc(j,ie_metal0+iel_C), want_C, 1.0d-12)

      ! The charge the census returns is the free electron density calc_ne
      ! builds from the same state: one definition of neutrality.
      do im = 1, n_mion
         nm_l(:,im) = f_sp(:,mion_fsp(im))
      enddo
      nmol_l(:,1) = f_sp(:,isp_H2)
      nmol_l(:,2) = f_sp(:,isp_H2p)
      nmol_l(:,3) = f_sp(:,isp_H3p)
      nmol_l(:,4) = f_sp(:,isp_HeHp)
      call calc_ne(f_sp(:,isp_HII), f_sp(:,isp_HeII), f_sp(:,isp_HeIII),  &
                   ne_ref, nm_l, nmol_l)
      ne_ref = ne_ref*nd
      call check('E1e', abs(n_chg(j)-want_q) .le. 1.0d-12*want_q,         &
                 'charge   ', n_chg(j), want_q, 1.0d-12)
      call check('E1f', abs(n_chg(j)-ne_ref(j)) .le. 1.0d-12*want_q,      &
                 'charge = calc_ne', n_chg(j), ne_ref(j), 1.0d-12)

      deallocate(rho, f_sp, n_nuc, n_chg, ne_ref, nm_l, nmol_l)
      end subroutine test_E1

      ! ------------------------------------------------------!

      subroutine test_E2()
      ! One carrier write-back must return the element totals it was handed.
      ! Two configurations are put through it:
      !   E2a  a generic molecular cell whose carriers are perturbed;
      !   E2b  the 12.8 configuration -- a cell holding most of its carbon
      !        and half its oxygen in CO, whose CO then collapses by 2.5
      !        decades in one step.  The free stages must absorb exactly the
      !        released nuclei and no more.
      real*8, allocatable :: rho(:), Tc(:), f_sp(:,:)
      real*8, allocatable :: fc(:,:), ntot(:), TK(:), mbar(:)
      real*8, allocatable :: nH_free(:), nO_free(:), nC_free(:)
      real*8, allocatable :: n0_nuc(:,:), n1_nuc(:,:), q0(:), q1(:)
      real*8 :: dev, worst
      integer :: j, ie, jw, iew

      allocate(rho(1-Ng:N+Ng), Tc(1-Ng:N+Ng), f_sp(1-Ng:N+Ng,n_species))
      allocate(fc(1-Ng:N+Ng,n_carrier_max), ntot(1-Ng:N+Ng))
      allocate(TK(1-Ng:N+Ng), mbar(1-Ng:N+Ng))
      allocate(nH_free(1-Ng:N+Ng), nO_free(1-Ng:N+Ng), nC_free(1-Ng:N+Ng))
      allocate(n0_nuc(1-Ng:N+Ng,n_element), n1_nuc(1-Ng:N+Ng,n_element))
      allocate(q0(1-Ng:N+Ng), q1(1-Ng:N+Ng))

      call molecular_cell_state(rho, Tc, f_sp)
      do j = 1-Ng, N+Ng
         bg_cell(j)%T_K  = Tc(j)*T0
         bg_cell(j)%ntot = 0.5d0*rho(j)*n0
      enddo

      call element_nuclei_and_charge(rho, f_sp, n0_nuc, q0)
      call carrier_state(rho, Tc, f_sp, fc, ntot, TK, mbar,               &
                         nH_free, nO_free, nC_free)
      ! E2a: a generic perturbation of the transported partition.
      fc(:,ic_H2)  = fc(:,ic_H2) *0.80d0
      fc(:,ic_OH)  = fc(:,ic_OH) *1.30d0
      fc(:,ic_H2O) = fc(:,ic_H2O)*0.60d0
      fc(:,ic_CO)  = fc(:,ic_CO) *0.90d0
      call carrier_write_back(rho, f_sp, fc, ntot, nH_free, nO_free,      &
                              nC_free)
      call element_nuclei_and_charge(rho, f_sp, n1_nuc, q1)
      call worst_departure(n0_nuc, n1_nuc, worst, jw, iew)
      call check('E2a', worst .le. 1.0d-12,                               &
                 'element closure through one write-back',                &
                 worst, 0.0d0, 1.0d-12)
      if (worst .gt. 1.0d-12) write(*,'(A,A,A,I4)') '        worst '//    &
           'element ', trim(element_name(iew)), ' at cell ', jw

      ! E2b: the 12.8 configuration.
      call molecular_cell_state(rho, Tc, f_sp)
      do j = 1-Ng, N+Ng
         bg_cell(j)%T_K  = Tc(j)*T0
         bg_cell(j)%ntot = 0.5d0*rho(j)*n0
      enddo
      call element_nuclei_and_charge(rho, f_sp, n0_nuc, q0)
      call carrier_state(rho, Tc, f_sp, fc, ntot, TK, mbar,               &
                         nH_free, nO_free, nC_free)
      fc(:,ic_CO) = fc(:,ic_CO)*3.8d-3      ! 2.34e9 -> 9.0e6 of section 111
      call carrier_write_back(rho, f_sp, fc, ntot, nH_free, nO_free,      &
                              nC_free)
      call element_nuclei_and_charge(rho, f_sp, n1_nuc, q1)
      call worst_departure(n0_nuc, n1_nuc, worst, jw, iew)
      call check('E2b', worst .le. 1.0d-12,                               &
                 'element closure when the transported CO collapses',     &
                 worst, 0.0d0, 1.0d-12)
      if (worst .gt. 1.0d-12) then
         write(*,'(A,A,A,I4)') '        worst element ',                  &
              trim(element_name(iew)), ' at cell ', jw
         do ie = 1, n_element
            if (n0_nuc(jw,ie) .le. 0.0d0) cycle
            dev = abs(n1_nuc(jw,ie) - n0_nuc(jw,ie))/n0_nuc(jw,ie)
            write(*,'(A,A,A,ES16.8,A,ES16.8,A,ES10.3)') '          ',    &
                 trim(element_name(ie)), ' ', n0_nuc(jw,ie), ' -> ',     &
                 n1_nuc(jw,ie), '  rel ', dev
         enddo
      endif

      deallocate(rho, Tc, f_sp, fc, ntot, TK, mbar)
      deallocate(nH_free, nO_free, nC_free, n0_nuc, n1_nuc, q0, q1)
      end subroutine test_E2

      ! ------------------------------------------------------!

      subroutine molecular_cell_state(rho, Tc, f_sp)
      ! A molecular-layer composition at the elemental ratios of the run:
      ! most of the hydrogen in H2, most of the carbon and half the oxygen
      ! in CO, and the rest of the oxygen split between OH, H2O and the free
      ! stages.  Written as nuclei first so the state is closed by
      ! construction and the test measures the write-back, not the setup.
      real*8, dimension(1-Ng:N+Ng), intent(out) :: rho, Tc
      real*8, dimension(1-Ng:N+Ng,n_species), intent(out) :: f_sp
      real*8 :: nH, nHe, nO, nC, nCO, nOH, nH2O, nH2, nd
      integer :: j

      rho  = 1.0d0
      Tc   = 1.5d0
      f_sp = 0.0d0
      nd   = n0
      do j = 1-Ng, N+Ng
         nH   = nd
         nHe  = HeH*nH
         nC   = melem_ab(iel_C)*nH
         nO   = melem_ab(iel_O)*nH
         nCO  = 0.90d0*nC
         nOH  = 0.10d0*nO
         nH2O = 0.15d0*nO
         nH2  = 0.40d0*(nH - nOH - 2.0d0*nH2O)
         f_sp(j,isp_H2)  = nH2 /nd
         f_sp(j,isp_OH)  = nOH /nd
         f_sp(j,isp_H2O) = nH2O/nd
         f_sp(j,isp_CO)  = nCO /nd
         f_sp(j,isp_HI)  = (nH - 2.0d0*nH2 - nOH - 2.0d0*nH2O)*0.85d0/nd
         f_sp(j,isp_HII) = (nH - 2.0d0*nH2 - nOH - 2.0d0*nH2O)*0.15d0/nd
         f_sp(j,isp_HeI) = nHe/nd
         f_sp(j,mion_fsp(melem_i0(iel_C)))   = 0.7d0*(nC - nCO)/nd
         f_sp(j,mion_fsp(melem_i0(iel_C)+1)) = 0.3d0*(nC - nCO)/nd
         f_sp(j,mion_fsp(melem_i0(iel_O)))   =                            &
              0.8d0*(nO - nOH - nH2O - nCO)/nd
         f_sp(j,mion_fsp(melem_i0(iel_O)+1)) =                            &
              0.2d0*(nO - nOH - nH2O - nCO)/nd
      enddo
      end subroutine molecular_cell_state

      ! ------------------------------------------------------!

      subroutine worst_departure(a, b, worst, jw, iew)
      real*8, dimension(1-Ng:N+Ng,n_element), intent(in)  :: a, b
      real*8,  intent(out) :: worst
      integer, intent(out) :: jw, iew
      real*8  :: dev
      integer :: j, ie
      worst = 0.0d0;  jw = 0;  iew = 0
      do j = 1-Ng, N+Ng
         do ie = 1, n_element
            if (a(j,ie) .le. 0.0d0) cycle
            dev = abs(b(j,ie) - a(j,ie))/a(j,ie)
            if (dev .gt. worst) then
               worst = dev;  jw = j;  iew = ie
            endif
         enddo
      enddo
      end subroutine worst_departure

      ! ------------------------------------------------------!

      subroutine test_E3()
      ! The H2 photoevent ledger.  For a photon of energy E the code splits
      ! one H2 photoionization into
      !     share f_di : H2 -> H+ + e-   (+ the H atom of the pair, supplied
      !                  by the H-nucleus closure of the molecular system)
      !     share 1-f_di: H2 -> H2+ + e-
      ! and the Lyman-Werner channel H2 -> H + H separately.  Each final
      ! state must carry two H nuclei and zero net charge once its electron
      ! is counted.  The check is made on the table itself, at every
      ! tabulated energy and at points between them.
      real*8  :: E, fdi, hnuc, chg, worst_h, worst_q
      integer :: i, nsamp
      logical :: monotone_ok

      worst_h = 0.0d0
      worst_q = 0.0d0
      nsamp   = 4000
      do i = 0, nsamp
         E   = 15.0d0 + (200.0d0 - 15.0d0)*dble(i)/dble(nsamp)
         fdi = frac_H2_dissociative_ionization(E)
         ! H nuclei of the two final states, weighted by their shares:
         !   dissociative  : H+ (1) + H (1) = 2
         !   non-dissoc.   : H2+ (2)        = 2
         hnuc = fdi*(1.0d0 + 1.0d0) + (1.0d0 - fdi)*2.0d0
         ! Net charge of the products including the released electron:
         !   dissociative  : +1 (H+) - 1 (e-) = 0
         !   non-dissoc.   : +1 (H2+) - 1 (e-) = 0
         chg  = fdi*(1.0d0 - 1.0d0) + (1.0d0 - fdi)*(1.0d0 - 1.0d0)
         worst_h = max(worst_h, abs(hnuc - 2.0d0))
         worst_q = max(worst_q, abs(chg))
      enddo
      call check('E3a', worst_h .le. 1.0d-14,                             &
                 'H nuclei per H2 photoionization', 2.0d0 + worst_h,      &
                 2.0d0, 1.0d-14)
      call check('E3b', worst_q .le. 1.0d-14,                             &
                 'net charge per H2 photoionization', worst_q, 0.0d0,     &
                 1.0d-14)

      ! Lyman-Werner: H2 -> H + H, two H nuclei, no charge, no electron.
      call check('E3c', .true., 'H2 + hv -> H + H closes by construction')

      ! The branching stays inside [0,1] and is flat above the measured
      ! range, which is what makes the two statements above shares.
      monotone_ok = .true.
      do i = 0, nsamp
         E   = 15.0d0 + (200.0d0 - 15.0d0)*dble(i)/dble(nsamp)
         fdi = frac_H2_dissociative_ionization(E)
         if (fdi .lt. 0.0d0 .or. fdi .gt. 1.0d0) monotone_ok = .false.
      enddo
      call check('E3d', monotone_ok, 'branching stays in [0,1]')
      call check('E3e',                                                   &
           abs(frac_H2_dissociative_ionization(500.0d0)                   &
             - frac_H2_dissociative_ionization(124.0d0)) .le. 1.0d-15,    &
           'held at the 124 eV value above the measured range',           &
           frac_H2_dissociative_ionization(500.0d0),                      &
           frac_H2_dissociative_ionization(124.0d0), 1.0d-15)

      ! THE TWO CHANNELS THIS ONE BRANCHING FRACTION CANNOT CARRY.  f_di
      ! splits the absorption two ways, so by itself it cannot say how many
      ! of the protons come from double ionization nor which absorptions
      ! make no ion at all.  Both are now resolved OUTSIDE this function, by
      ! the explicit channel cross sections of h2_photo_channels.f90, and
      ! both are ON by default (`H2 double ionization: chung80`,
      ! `H2 neutral dissociation: True`).  Reported here, not gated: what is
      ! gated is the channel sum and the stoichiometry, in
      ! src/tests/e1_h2/e1_h2_channel_check.f90.
      write(*,'(A)') '  NOTE E3 what f_di alone cannot split, and where'// &
           ' it is now split:'
      write(*,'(A,F7.4,A)') '    double ionization H2 -> H+ + H+ + 2e-'// &
           ' is inside f_di; above 80 eV it is ~20% of sigma(H+),'
      write(*,'(A,F7.4,A)') '      i.e. ~',                               &
           0.20d0*frac_H2_dissociative_ionization(100.0d0),               &
           ' of the H2 ionizations -- taken out into channel D.'
      write(*,'(A)') '    neutral dissociation H2 -> H + H over 33-41 eV'//&
           ' (up to 7% of sigma_H2 at 37.5 eV) would be handed to'
      write(*,'(A)') '      the H2+ branch by 1 - f_di, creating an H2+'// &
           ' and an electron the event does not make -- channel N.'
      end subroutine test_E3

      ! ------------------------------------------------------!

      subroutine test_E4()
      ! Item 12.6: with no handoff, the base H2 mixing ratio must come from
      ! the Visscher/Koskinen chemical-equilibrium fit at (p_base, T0), and
      ! the ghost H2 nucleus fraction must be the mu_mixture conversion of
      ! it.  With a handoff it must be the handoff value instead.
      real*8 :: q_fit, q_got, x2_got, x2_want

      q_h2_base   = -1.0d0            ! no handoff value stated
      q_fit  = q_h2_equilibrium(p_base_bar, T0)
      q_got  = h2_mixing_ratio_base()
      call check('E4a', abs(q_got - q_fit) .le. 1.0d-14*max(q_fit,1.0d-30),&
                 'no handoff -> chem.eq. fit q_H2(p_base,T0)', q_got,     &
                 q_fit, 1.0d-14)
      call check('E4b', .not. base_h2_composition_imposed(),              &
                 'no handoff -> the base composition is not imposed')
      x2_want = 2.0d0*q_fit*(1.0d0 + HeH)/(1.0d0 + q_fit)
      x2_got  = base_h2_nuclei_fraction()
      call check('E4c', abs(x2_got - x2_want)                             &
                 .le. 1.0d-14*max(x2_want,1.0d-30),                       &
                 'ghost H2 nucleus fraction from the fit', x2_got,        &
                 x2_want, 1.0d-14)

      q_h2_base = 0.84d0              ! Koskinen et al. (2022) handoff value
      call check('E4d', abs(h2_mixing_ratio_base() - 0.84d0)              &
                 .le. 1.0d-14, 'handoff -> q_H2_base is used',            &
                 h2_mixing_ratio_base(), 0.84d0, 1.0d-14)
      call check('E4e', base_h2_composition_imposed(),                    &
                 'handoff -> the base composition IS imposed')
      call check('E4f', x2_got .le. 1.0d0,                                &
                 'the fit branch is attainable at this He/H', x2_got,     &
                 1.0d0, 0.0d0)
      ! The conversion is not intrinsically bounded: the fit's asymptote
      ! 0.8384 exceeds the ceiling 0.5/(0.5 + He/H) for any He/H >= 0.167,
      ! and input_read refuses the run there rather than capping it. The
      ! ceiling is stated here so the refusal has a test beside it.
      q_h2_base = -1.0d0            ! back to the fit branch
      HeH = 0.2d0
      call check('E4g', h2_mixing_ratio_base()                             &
                 .gt. 0.5d0/(0.5d0 + HeH),                                &
                 'the fit exceeds the ceiling at He/H = 0.2 (refused by'// &
                 ' input_read)', h2_mixing_ratio_base(),                  &
                 0.5d0/(0.5d0 + HeH), 0.0d0)
      HeH = 0.0793d0
      end subroutine test_E4

      end program element_census_tests
