      program opacity_model_p_ledger
      ! QUANTITY UNDER TEST
      !   the radiative ledger of a cell whose photoionization cross
      !   sections are broadened by the factor f(p) of opacity model 'P'
      !   (opacity_models.f90, opacity_pT_factor; the run stores it in
      !   opa_pf).  Two conservation statements, both of the SAME cell:
      !
      !     (1) ENERGY.  The energy the cell absorbs,
      !           q_abs(j)*dr(j) = INT F(tau_j) sigma_eff n dr dE ,
      !         equals the energy the beam loses across it,
      !           INT F(tau_{j+1}) [1 - exp(-dtau_j)] dE .
      !
      !     (2) PHOTON NUMBER, H I channel.  The H I photoionizations,
      !           P_HI(j) n_HI(j) dr(j) ,
      !         equal the photons the beam loses across the cell times the
      !         H I share dtau_HI/dtau of its optical depth.
      !
      !   Both hold for ANY value of f only if f multiplies the local
      !   absorption exactly as it multiplies the optical depth.  The
      !   columns (calc_column_dens*) carry f; the assertions here are on
      !   the integrands of PH_heat_HHe, which is where the photons removed
      !   from the beam become ionizations and heat.
      !
      ! HOW THE CODE FORMS THEM.  PH_heat_HHe (util_ion_eq.f90) is called
      ! with the composition below; nothing is re-implemented.  q_abs of the
      ! cell is recovered from the two returned quantities that define it,
      ! q = heat/q_abs, as q_abs = heat/q.
      !
      ! REFERENCE.  The beam itself, on the same energy grid, with the same
      ! production cross sections (s_hi, s_hei, s_heii, s_heiTR, sigma_tab)
      ! and the same columns the code builds: dtau_j(E) = sigma(E) n dr
      ! opa_pf(j), tau_j = sum over the cells outside j inclusive of j, and
      ! the flux F(tau) = F_XUV exp(-tau) of a_tau = 0.  The reference is
      ! therefore the definition of absorption, not a second model of it.
      !
      ! TOLERANCE 1e-5 relative.  The identity holds at ANY cell optical
      ! depth: the code illuminates a cell with the MEAN of the field over
      ! it, exp(-tau_out) (1 - exp(-dtau))/dtau (utils,
      ! cell_mean_attenuation), which is exactly the beam loss divided by
      ! dtau.  What is left inside the tolerance is arithmetic: the
      ! configuration is thin, the largest optical depth of the broadened
      ! cell being 6.3e-8 (printed as a diagnostic), so the cancellation in
      ! 1 - exp(-dtau) leaves about 2e-9.  The defect this test targets is a
      ! FACTOR f = 10, five orders of magnitude outside it.
      !
      ! CONFIGURATION.  A twelve-cell uniform slab (N = 8 plus two ghost
      ! cells at each end) of H I, He I (ground singlet), He II, He 2^3S,
      ! C I and Mg II, illuminated by the power-law spectrum of
      ! backup/regression/wasp_full/input.inp.  opa_pf = 1 everywhere except
      ! the broadened cell j = 4, where it is 10.  opa_pf is set directly
      ! rather than through opacity.inp: it is the state model 'P' leaves
      ! (ionization_equilibrium.f90 fills it from opacity_pT_factor), and
      ! the ledger asserted here is a statement about that state, not about
      ! the pressure parameterization that produces it.  Both cells are
      ! asserted: the broadened one, and an unbroadened neighbour that shows
      ! the default path (f = 1, every model but 'P') is unaffected.
      !
      ! Secondary ionization is off, so P_HI is the primary rate.

      use global_parameters
      use species_table,          only: n_melem, n_mion, n_mphot,          &
                                        mion_isphot, mion_iphot,           &
                                        im_CI, im_MgII
      use energy_vectors_construct, only: set_energy_vectors
      use utils_ion_eq,           only: PH_heat_HHe
      use assertion_report

      implicit none

      ! Slab geometry.  dr_geo is the geometric width of every cell [cm];
      ! the densities and the width are chosen together so that the optical
      ! depth of one cell stays below 1e-6 at every photon energy (see the
      ! tolerance note above).
      real*8, parameter :: dr_geo   = 1.0d4
      real*8, parameter :: n_HI     = 1.0d5
      real*8, parameter :: n_HeIS   = 1.0d4
      real*8, parameter :: n_HeII   = 1.0d4
      real*8, parameter :: n_HeTR   = 1.0d1
      real*8, parameter :: n_CI     = 1.0d1
      real*8, parameter :: n_MgII   = 1.0d1
      real*8, parameter :: f_broad  = 10.0d0
      integer, parameter :: j_broad = 4
      integer, parameter :: j_plain = 6

      real*8, dimension(:), allocatable :: nhi,nhei,nheii,nheiTR,xion
      real*8, dimension(:), allocatable :: P_HI,P_HeI,P_HeII,P_HeITR
      real*8, dimension(:), allocatable :: heat,q
      real*8, dimension(:,:), allocatable :: nm, P_m
      ! Optical depth of one cell at unit opa_pf, and its H I part.
      real*8, dimension(:), allocatable :: dtau1, dtau1_HI
      real*8 :: dtau_max

      ! --- the run-wide input, as input_read resolves it ------------------
      is_PL_sed    = .true.
      do_read_sed  = .false.
      is_monochr   = .false.
      PLind        = -1.0d0
      thereis_Xray = .true.
      e_low        = 13.60d0
      e_mid        = 123.98d0
      e_top        = 1.24d3
      LX           = 29.46d0
      LEUV         = 30.42d0
      a_orb        = 0.02544d0*AU
      appx_mth     = 'Rate/2 + Mdot/2'
      ! a_tau = 0 makes the field a pure exponential, so the beam the
      ! reference differences is the beam the code integrates.
      a_tau        = 0.0d0
      thereis_He   = .true.
      thereis_HeITR       = .true.
      thereis_lowIP_metal = .false.
      use_sec_ion    = .false.
      sec_ion_active = .false.
      allocate(melem_ab(n_melem))
      melem_ab     = 0.0d0

      call set_energy_vectors

      ! --- the grid -------------------------------------------------------
      N  = 8
      R0 = 1.0d10
      call allocate_grid_arrays
      dr_j   = dr_geo/R0
      opa_pf = 1.0d0
      opa_pf(j_broad) = f_broad

      allocate(nhi(1-Ng:N+Ng), nhei(1-Ng:N+Ng), nheii(1-Ng:N+Ng),          &
               nheiTR(1-Ng:N+Ng), xion(1-Ng:N+Ng))
      allocate(P_HI(1-Ng:N+Ng), P_HeI(1-Ng:N+Ng), P_HeII(1-Ng:N+Ng),       &
               P_HeITR(1-Ng:N+Ng), heat(1-Ng:N+Ng), q(1-Ng:N+Ng))
      allocate(nm(1-Ng:N+Ng,n_mion), P_m(1-Ng:N+Ng,n_mion))
      allocate(dtau1(Nl), dtau1_HI(Nl))

      nhi    = n_HI
      nheiTR = n_HeTR
      ! PH_heat_HHe takes the SUMMED neutral helium and subtracts the
      ! metastable itself (he_ground_singlet_density), so the singlet the
      ! reference uses is n_HeIS by construction.
      nhei   = n_HeIS + n_HeTR
      nheii  = n_HeII
      xion   = 0.5d0
      nm     = 0.0d0
      nm(:,im_CI)   = n_CI
      nm(:,im_MgII) = n_MgII

      call PH_heat_HHe(nhi,nhei,nheii,nheiTR, nm, xion,                    &
                       P_HI,P_HeI,P_HeII,P_HeITR, P_m, heat,q)

      call cell_optical_depth(dtau1, dtau1_HI)
      dtau_max = maxval(dtau1)*f_broad
      write(*,'(a,es10.3,a,i0,a,es10.3)')                                  &
           '  DIAGNOSTIC max optical depth of the broadened cell = ',      &
           dtau_max, '  over Nl=', Nl, ' points; grid floor [eV] = ', e_v(1)

      call ledger_of_cell(j_broad, 'opa_pf_10')
      call ledger_of_cell(j_plain, 'opa_pf_1')

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'opacity_model_p_ledger: ',                   &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'opacity_model_p_ledger: all assertions passed'

      contains

      !--------------!

      subroutine cell_optical_depth(dtau_unit, dtau_unit_HI)
      ! Optical depth of one cell of the slab at opa_pf = 1, and its H I
      ! part, from the production cross sections.  The cross-section arrays
      ! carry the code's internal unit of 1e-18 cm^2.
      real*8, dimension(Nl), intent(out) :: dtau_unit, dtau_unit_HI
      integer :: i,k
      dtau_unit_HI = s_hi*n_HI*dr_geo*1.0d-18
      dtau_unit    = dtau_unit_HI                                          &
                   + (s_hei*n_HeIS + s_heii*n_HeII + s_heiTR*n_HeTR)       &
                     *dr_geo*1.0d-18
      do i = 1,n_mion
         if (.not. mion_isphot(i)) cycle
         k = mion_iphot(i)
         dtau_unit = dtau_unit                                             &
                   + sigma_tab(:,k)*nm(1,i)*dr_geo*1.0d-18
      enddo
      end subroutine cell_optical_depth

      !--------------!

      subroutine ledger_of_cell(jc, label)
      ! The two conservation statements for cell jc.
      integer,          intent(in) :: jc
      character(len=*), intent(in) :: label
      integer :: k
      real*8, dimension(Nl) :: tau_out, dtau_c, loss, share_HI
      real*8 :: e_absorbed, e_beam_loss, ph_beam_loss_HI, ionizations
      real*8 :: q_abs

      ! Column of everything OUTSIDE the cell, the code's tau at cell jc+1.
      tau_out = 0.0d0
      do k = jc+1, N+Ng
         tau_out = tau_out + dtau1*opa_pf(k)
      enddo
      dtau_c = dtau1*opa_pf(jc)

      ! Energy and photons the beam loses across the cell.  Written as
      ! F(tau_out)*(1 - exp(-dtau)) rather than as a difference of two
      ! nearly equal exponentials: dtau is formed exactly, so the only
      ! cancellation left is in 1 - exp(-dtau), of order 1e-16/dtau.
      loss = F_XUV*exp(-tau_out)*(1.0d0 - exp(-dtau_c))

      ! Share of the removed photons that H I takes.  A photon energy at
      ! which nothing absorbs removes no photons either, so the share there
      ! multiplies a zero and is set to zero rather than formed as 0/0.
      share_HI = 0.0d0
      where (dtau_c .gt. 0.0d0) share_HI = dtau1_HI*opa_pf(jc)/dtau_c

      e_beam_loss     = sum(loss*de_v)
      ph_beam_loss_HI = sum(loss*erg2eV/e_v*share_HI*de_v)

      ! What the code says the cell absorbs.  q is defined as heat/q_abs by
      ! PH_heat_HHe, so this recovers its own absorbed-energy integral.
      q_abs       = heat(jc)/q(jc)
      e_absorbed  = q_abs*dr_geo
      ionizations = P_HI(jc)*nhi(jc)*dr_geo

      call check_relative('energy_absorbed_equals_beam_loss['//label//']', &
                          e_absorbed, e_beam_loss, 1.0d-5)
      call check_relative('HI_ionizations_equal_beam_photon_loss['         &
                          //label//']', ionizations, ph_beam_loss_HI,      &
                          1.0d-5)
      end subroutine ledger_of_cell

      end program opacity_model_p_ledger
