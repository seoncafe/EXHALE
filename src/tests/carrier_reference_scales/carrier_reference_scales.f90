      program carrier_reference_scales
      ! The element a transported carrier is made of, and whether the
      ! chemistry rows of that carrier can be differentiated at all.
      !
      ! WHAT IS BEING TESTED, AND WHY IT IS PHYSICS AND NOT BOOKKEEPING.
      ! The carrier transport-chemistry Newton
      ! (src/modules/lower_atmosphere/diffusive_photochemistry.f90) needs two
      ! densities of every carrier that are not the carrier's own value:
      !
      !   * the perturbation with which its source rows are differentiated
      !     (carrier_source_derivative_step), and
      !   * the absolute floor below which a row is declared converged
      !     (the 1e-20 nfl term of carrier_residual).
      !
      ! Both are fractions of the free density of the element the carrier is
      ! MADE OF, because that is the only density in the cell of the size of
      ! the terms of that carrier's rows.  Stoichiometry fixes the pairing:
      ! H2 and H+ carry hydrogen nuclei, OH and H2O are bounded by their
      ! oxygen, CO needs one nucleus each of oxygen and carbon.
      !
      ! In a hydrogen and helium atmosphere the free oxygen density is
      ! identically zero.  A proton charged to the oxygen reservoir
      ! therefore gets a differentiation step of 1e-300 wherever its own
      ! density is small, and a forward difference over 1e-300 of a quantity
      ! built from densities of 1e10 returns the round-off of the source and
      ! not its derivative.  That is the state this suite measures.
      !
      ! ASSERTIONS.
      !  (1) the element table: which free density each carrier is scaled
      !      against.  An exact identity, tested to round-off.
      !  (2) the diagonal derivative d(src_H+)/d(n_H+) obtained with the
      !      production step, at three proton densities spanning the range
      !      the transported ionization front passes through, against a
      !      CENTRAL difference taken with an independently chosen step.
      !      The reference step is 1e-7 of the free hydrogen density, chosen
      !      here and nowhere read from the code under test.
      !  (3) the production step itself: 1e-12 of the free hydrogen density
      !      wherever the proton's own density is smaller than that, and
      !      1e-6 of the proton density where it is larger.
      !
      ! The gas is metal-free by construction (the O and C columns of f_sp
      ! are zero), which is the configuration of the hp_* regression cases
      ! and the one in which the proton row has no oxygen to be scaled by.
      !
      ! TOLERANCE OF (2), AND WHERE IT COMES FROM.  A forward difference
      ! over a step dn carries a truncation error of order dn |f''| / 2 and
      ! a cancellation error of order eps |f| / (dn |f'|).  At the two small
      ! proton densities the second dominates: the proton row is about
      ! 2e3 cm^-3 s^-1 here, its derivative about 3e-6, and the step is
      ! 1e-2 cm^-3, so the cancellation floor is eps*2e3/(1e-2*3e-6), about
      ! 4e-5 relative.  MEASURED: 3.7e-5 and 4.5e-5.  The tolerance is
      ! therefore 1e-3 relative -- one and a half decades above the
      ! cancellation floor of a correctly scaled step, and far below what a
      ! step of 1e-300 or of 1e-8 can reach: with the proton scaled to the
      ! oxygen reservoir the same two points return a derivative of exactly
      ! zero and one dominated by round-off, both at a relative error of
      ! order unity.

      use global_parameters
      use species_table
      use ion_cell_state, only: ion_rates
      use ionization_equilibrium, only: bg_cell, ioniz_eq_allocate_arrays
      use diffusive_photochemistry, only: carrier_set_init, carrier_state, &
                                          carrier_source,                  &
                                          carrier_element_reference_density,&
                                          carrier_source_derivative_step,   &
                                          n_carrier_max,                    &
                                          ic_H2, ic_OH, ic_H2O, ic_CO, ic_Hp
      use assertion_report
      implicit none

      ! Free element densities of the table test [cm^-3].  Three distinct
      ! values so that a routine returning the wrong one cannot pass.
      real*8, parameter :: nH_tab = 3.0d10
      real*8, parameter :: nO_tab = 7.0d6
      real*8, parameter :: nC_tab = 2.0d6
      ! The cell the derivative is taken in, and the reference step of the
      ! central difference as a fraction of the free hydrogen density.
      integer, parameter :: jcell = 3
      real*8, parameter  :: h_frac = 1.0d-7

      real*8, allocatable :: rho(:), f_sp(:,:), fc(:,:)
      real*8, allocatable :: ntot(:), TK(:), mbar(:), nrho(:), wfac(:)
      real*8, allocatable :: nH_free(:), nO_free(:), nC_free(:)
      real*8 :: nc(n_carrier_max), sp(n_carrier_max), sm(n_carrier_max)
      real*8 :: s0(n_carrier_max)
      real*8 :: dn, hh, dfwd, dcen, nhp, nHc
      character(len=64) :: label
      integer :: k

      call setup_globals()

      ! ---- (1) the element table ------------------------------------- !
      call check_relative('carrier_reference_H2_is_hydrogen',             &
           carrier_element_reference_density(ic_H2, nH_tab, nO_tab,       &
                                             nC_tab), nH_tab, 1.0d-12)
      call check_relative('carrier_reference_Hp_is_hydrogen',             &
           carrier_element_reference_density(ic_Hp, nH_tab, nO_tab,       &
                                             nC_tab), nH_tab, 1.0d-12)
      call check_relative('carrier_reference_OH_is_oxygen',               &
           carrier_element_reference_density(ic_OH, nH_tab, nO_tab,       &
                                             nC_tab), nO_tab, 1.0d-12)
      call check_relative('carrier_reference_H2O_is_oxygen',              &
           carrier_element_reference_density(ic_H2O, nH_tab, nO_tab,      &
                                             nC_tab), nO_tab, 1.0d-12)
      ! CO needs one nucleus of each, so the scarcer element bounds it.
      call check_relative('carrier_reference_CO_is_scarcer_of_O_and_C',   &
           carrier_element_reference_density(ic_CO, nH_tab, nO_tab,       &
                                             nC_tab), nC_tab, 1.0d-12)

      ! ---- the column the derivatives are taken on -------------------- !
      call build_molecular_hydrogen_column()

      ! The gas must really be metal-free, or assertion (2) tests nothing.
      call check_absolute('free_oxygen_density_is_zero',                  &
                          nO_free(jcell), 0.0d0, 0.0d0)
      call check_absolute('free_carbon_density_is_zero',                  &
                          nC_free(jcell), 0.0d0, 0.0d0)
      call check_positive('free_hydrogen_density_is_positive',            &
                          nH_free(jcell))

      ! ---- (2) and (3) the proton row --------------------------------- !
      nHc = nH_free(jcell)
      do k = 1, 3
         if (k .eq. 1) then
            nhp   = 0.0d0
            label = 'proton_zero'
         else if (k .eq. 2) then
            nhp   = 1.0d-12*nHc
            label = 'proton_1e-12_of_nH'
         else
            nhp   = 1.0d-3*nHc
            label = 'proton_1e-3_of_nH'
         endif

         dn = carrier_source_derivative_step(ic_Hp, nhp, nH_free(jcell),  &
                                             nO_free(jcell),             &
                                             nC_free(jcell))
         ! (3) the step itself.  Below 1e-6 of the free hydrogen density the
         ! proton's own value cannot set a resolvable step and the element
         ! floor has to; above it the step follows the carrier.
         if (1.0d-6*nhp .gt. 1.0d-12*nHc) then
            call check_relative('carrier_step_follows_carrier_'//         &
                                trim(label), dn, 1.0d-6*nhp, 1.0d-15)
         else
            call check_relative('carrier_step_on_hydrogen_floor_'//       &
                                trim(label), dn, 1.0d-12*nHc, 1.0d-15)
         endif

         ! Forward difference with the production step, exactly as the
         ! Jacobian block of solve_carriers builds its diagonal entry.
         nc = 0.0d0
         nc(ic_H2) = fc(jcell,ic_H2)*nrho(jcell)
         nc(ic_Hp) = nhp
         call carrier_source(jcell, nc, nH_free(jcell), nO_free(jcell), s0)
         nc(ic_Hp) = nhp + dn
         call carrier_source(jcell, nc, nH_free(jcell), nO_free(jcell), sp)
         dfwd = (sp(ic_Hp) - s0(ic_Hp))/dn

         ! Central difference with a step chosen here.  The rows are smooth
         ! algebraic functions of n(H+) around these points, so a symmetric
         ! difference is second-order accurate and independent of the code
         ! under test.
         hh = h_frac*nHc
         nc(ic_Hp) = nhp + hh
         call carrier_source(jcell, nc, nH_free(jcell), nO_free(jcell), sp)
         nc(ic_Hp) = nhp - hh
         call carrier_source(jcell, nc, nH_free(jcell), nO_free(jcell), sm)
         dcen = (sp(ic_Hp) - sm(ic_Hp))/(2.0d0*hh)

         call check_relative('proton_row_diagonal_derivative_'//          &
                             trim(label), dfwd, dcen, 1.0d-3)
      enddo

      if (assertion_failures .gt. 0) stop 1

      contains

      !--------------!

      subroutine setup_globals()
      ! A five-cell hydrogen and helium column with the molecular chemistry
      ! and the transported proton on, and no metals in the gas.
      integer :: j
      N   = 5
      n0  = 1.0d10
      R0  = 1.0d10
      v0  = 1.0d6
      T0  = 1000.0d0
      HeH = 0.0793d0
      p_base_bar         = 1.0d-6
      thereis_He         = .true.
      thereis_HeITR      = .false.
      thereis_metals     = .true.
      thereis_mol        = .true.
      thereis_oxychem    = .false.
      eos_include_metals = .false.
      he_diffusion       = .false.
      carrier_transport  = .true.
      ionization_transport = .true.
      allocate(r(1-Ng:N+Ng), r_edg(1-Ng:N+Ng), dr_j(1-Ng:N+Ng))
      do j = 1-Ng, N+Ng
         r(j)     = 1.0d0 + 0.01d0*dble(j)
         r_edg(j) = 1.0d0 + 0.01d0*(dble(j) + 0.5d0)
         dr_j(j)  = 0.01d0
      enddo
      allocate(melem_ab(n_melem))
      melem_ab = 0.0d0
      call carrier_set_init()
      call ioniz_eq_allocate_arrays()
      end subroutine setup_globals

      !--------------!

      subroutine build_molecular_hydrogen_column()
      ! Molecular base gas: 80 per cent of the hydrogen nuclei in H2, the
      ! rest atomic, helium at the solar-like ratio, a trace of the
      ! molecular ions, and no oxygen or carbon anywhere.
      integer :: j
      allocate(rho(1-Ng:N+Ng), f_sp(1-Ng:N+Ng,n_species))
      allocate(fc(1-Ng:N+Ng,n_carrier_max))
      allocate(ntot(1-Ng:N+Ng), TK(1-Ng:N+Ng), mbar(1-Ng:N+Ng))
      allocate(nrho(1-Ng:N+Ng), wfac(1-Ng:N+Ng))
      allocate(nH_free(1-Ng:N+Ng), nO_free(1-Ng:N+Ng), nC_free(1-Ng:N+Ng))
      rho  = 1.0d0
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
      call carrier_state(rho, f_sp, fc, ntot, nrho, wfac, TK, mbar,       &
                         nH_free, nO_free, nC_free)
      end subroutine build_molecular_hydrogen_column

      !--------------!

      subroutine set_frozen_cell_rates(j)
      ! Every field of the frozen cell state, set here so that none is read
      ! undefined.  The photoionization rate and the case-B recombination
      ! coefficient are the two that make the proton row depend on n(H+) at
      ! all: production from H I and the recombination sink alpha n_e n_H+,
      ! with n_e itself carrying the trial proton.
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

      end program carrier_reference_scales
