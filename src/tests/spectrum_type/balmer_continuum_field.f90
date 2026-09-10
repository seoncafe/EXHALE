      program balmer_continuum_field
      ! QUANTITY UNDER TEST
      !   the two scalar Balmer-continuum rates of H(n=2),
      !     gamma_2 = INT_{3.40}^{13.6 eV} F_E(E)/(E) sigma_2(E) dE   [s^-1]
      !     h_2     = the same integrand weighted by (E - 3.40 eV) [erg s^-1]
      !   with sigma_2(E) = 1.4e-17 (3.40 eV/E)^3 cm^2, as the production
      !   functions gamma_n2_balmer and heat_n2_balmer form them, and the
      !   field F_E they read.
      !
      ! WHAT DECISION 13 REQUIRES (docs/development_plan_20260905_rev3.md
      ! section 10.5): ONE spectrum type builds every band, the Balmer
      ! continuum included.  Until this change the two integrals were a
      ! photospheric blackbody in EVERY run, including the power-law runs of
      ! the whole regression matrix, so one run described the 3.4-13.6 eV
      ! band two ways.  They must now read the run's own type: the power law
      ! for "Power-law", the blackbody for "Planck", the loaded table for
      ! "Load", and a band no type states must stop the run.
      !
      ! REFERENCES OF THE ASSERTIONS
      !   1,2  power law   the two integrals against the same integral of the
      !                    production J_inc on 40001 uniform points, written
      !                    out here.  Tolerance 1e-3 (a quadrature against a
      !                    finer one; the production rule uses 400 points).
      !   3,4  Planck      the two integrals against the FREQUENCY writing of
      !                    the same band integral, INT F_nu/(h nu) sigma_2
      !                    dnu with nu = E*eV2Hz, on 400001 uniform points of
      !                    the band.  F_E dE = F_nu dnu, so the two writings
      !                    are the same number and what is left is the
      !                    quadrature error of the reference: tolerance 1e-6
      !                    relative, against a MEASURED 3.7e-10.
      !   5    Load        the field returned inside a synthetic two-decade
      !                    power-law table is that power law, tolerance 1e-12
      !                    relative (log-log interpolation of a power law is
      !                    exact), and the dayside dilution is divided out
      !                    exactly once.
      !   6-9  coverage    spectrum_covers_eV at the n=2 edge: true for the
      !                    power law and for Planck, false for a monochromatic
      !                    run and for a loaded table that stops at 13.6 eV.
      !                    Tolerance 0.  These are the configurations in which
      !                    the production functions stop the run.
      !   10,11 floor      a loaded table read down to the grid floor, whose
      !                    lowest selected row sits a fraction of a row above
      !                    it: the floor is covered (tolerance 0) and the
      !                    field there is the table's own, continued from the
      !                    first segment (1e-12 relative).
      !
      ! CONFIGURATION: the stellar and spectral input of
      ! backup/regression/wasp_full/input.inp -- T_eff 6459 K,
      ! R_star 1.458 Rsun, a 0.02544 AU, power law index -1 on
      ! [13.60, 123.98, 1.24e3] eV, LX 29.46, LEUV 30.42, "2D approximate
      ! method: Rate/2 + Mdot/2" (xi = 1/2).
      !
      ! EXPECTED BEFORE THE CHANGE: RED (the program does not build:
      ! stellar_flux_eV and spectrum_covers_eV do not exist and the two
      ! production functions take a temperature and a dilution factor).
      ! After: GREEN.

      use global_parameters
      use species_table,     only: n_melem
      use J_incident,        only: J_inc, planck_stellar_flux_nu,           &
                                   planck_stellar_flux_eV, stellar_flux_eV, &
                                   spectrum_covers_eV, eV2Hz, e_th_HI_n2
      use excited_hydrogen,  only: gamma_n2_balmer, heat_n2_balmer

      implicit none

      ! The n=2 edge and the threshold cross section, as the code states
      ! them (Osterbrock & Ferland 2006 hydrogenic n=2 continuum): the edge
      ! from J_incident, which is where the one definition lives, so this
      ! reference cannot drift away from the production band.
      real*8, parameter :: E2 = e_th_HI_n2
      real*8, parameter :: s2 = 1.4d-17

      integer :: nfail
      real*8  :: g_code, h_code, g_ref, h_ref
      real*8  :: g_pl, g_bb

      nfail = 0

      ! --- the run-wide input, as input_read resolves it for wasp_full ----
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
      a_tau        = 0.0d0
      T_star_eff   = 6459.0d0
      R_star       = 1.458d0*Rsun
      J_XUV        = (10.0d0**LX + 10.0d0**LEUV)/(4.0d0*pi*a_orb**2.0)
      allocate(melem_ab(n_melem))
      melem_ab     = 0.0d0

      ! ---- 1,2  the power-law type -------------------------------------
      sp_type   = 'Power-law'
      is_PL_sed = .true.
      g_code = gamma_n2_balmer()
      h_code = heat_n2_balmer()
      call power_law_reference(40001, g_ref, h_ref)
      call relative_verdict('balmer_rate_power_law', g_code, g_ref,         &
                            1.0d-3, nfail)
      call relative_verdict('balmer_heat_power_law', h_code, h_ref,         &
                            1.0d-3, nfail)
      g_pl = g_code

      ! ---- 3,4  the Planck type, against the frequency writing ----------
      sp_type   = 'Planck'
      is_PL_sed = .false.
      g_code = gamma_n2_balmer()
      h_code = heat_n2_balmer()
      call planck_frequency_reference(400001, g_ref, h_ref)
      call relative_verdict('balmer_rate_planck_frequency_form', g_code,    &
                            g_ref, 1.0d-6, nfail)
      call relative_verdict('balmer_heat_planck_frequency_form', h_code,    &
                            h_ref, 1.0d-6, nfail)
      g_bb = g_code

      write(*,'(A,ES11.4,A,ES11.4,A,ES10.3)')                              &
         'DIAGNOSTIC n=2 photoionization rate [s^-1] before dilution:'//    &
         ' power law ', g_pl, '  Planck ', g_bb, '  ratio ', g_pl/g_bb

      ! ---- 5  the loaded-table branch ------------------------------------
      call loaded_table_field(nfail)

      ! ---- 6-9  which types state a field at the n=2 edge -----------------
      sp_type   = 'Power-law'
      is_PL_sed = .true.
      do_read_sed = .false.
      is_monochr  = .false.
      call logical_verdict('balmer_band_covered_power_law',                &
                           spectrum_covers_eV(E2), .true., nfail)
      sp_type   = 'Planck'
      is_PL_sed = .false.
      call logical_verdict('balmer_band_covered_planck',                   &
                           spectrum_covers_eV(E2), .true., nfail)
      sp_type    = 'Monochromatic'
      is_monochr = .true.
      call logical_verdict('balmer_band_uncovered_monochromatic',          &
                           spectrum_covers_eV(E2), .false., nfail)
      is_monochr  = .false.
      sp_type     = 'Load'
      do_read_sed = .true.
      call short_table(13.6d0, 1.24d3)
      call logical_verdict('balmer_band_uncovered_short_table',            &
                           spectrum_covers_eV(E2), .false., nfail)

      ! ---- 10,11  a table that reaches the floor, sampled at the floor ---
      call table_reaching_the_grid_floor(nfail)

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'balmer_continuum_field: ', nfail,            &
                             ' assertion(s) FAILED'
         call exit(1)
      endif
      write(*,'(A)') 'balmer_continuum_field: all assertions PASSED'

      contains

      subroutine power_law_reference(npt, gam, hea)
      ! The two integrals of the power-law field, written out here from the
      ! production J_inc on npt uniform points of the band.
      integer, intent(in)  :: npt
      real*8,  intent(out) :: gam, hea
      integer :: i
      real*8  :: E, dE, f, sig, wg, wh

      dE  = (e_th_HI - E2)/dble(npt-1)
      gam = 0.0d0
      hea = 0.0d0
      do i = 1, npt
         E   = E2 + dble(i-1)*dE
         f   = J_inc(E)
         sig = s2*(E2/E)**3.0
         wg  = f/(hp_erg*E*eV2Hz)*sig
         wh  = wg*(hp_erg*(E - E2)*eV2Hz)
         if (i .eq. 1 .or. i .eq. npt) then
            wg = 0.5d0*wg
            wh = 0.5d0*wh
         endif
         gam = gam + wg
         hea = hea + wh
      enddo
      gam = gam*dE
      hea = hea*dE
      end subroutine power_law_reference

      subroutine planck_frequency_reference(n_nu, gam, hea)
      ! The same band integral written in FREQUENCY, of pi B_nu (R/a)^2, on
      ! n_nu uniform points of the band.  The production functions integrate
      ! the field on its own nodes, so the reference has to be fine enough
      ! that its own truncation error is what this comparison measures.
      integer, intent(in) :: n_nu
      real*8, intent(out) :: gam, hea
      integer :: i
      real*8  :: nu, dnu, Fnu, sig, nu_a, nu_b, wg, wh, nu2

      nu2  = E2*eV2Hz
      nu_a = E2*eV2Hz
      nu_b = e_th_HI*eV2Hz
      dnu  = (nu_b - nu_a)/dble(n_nu-1)
      gam  = 0.0d0
      hea  = 0.0d0
      do i = 1, n_nu
         nu  = nu_a + dble(i-1)*dnu
         Fnu = planck_stellar_flux_nu(T_star_eff,                          &
                                      R_star/max(a_orb,1.0d-30), nu)
         sig = s2*(nu2/nu)**3.0
         wg  = Fnu/(hp_erg*nu)*sig
         wh  = wg*hp_erg*(nu - nu2)
         if (i .eq. 1 .or. i .eq. n_nu) then
            wg = 0.5d0*wg
            wh = 0.5d0*wh
         endif
         gam = gam + wg
         hea = hea + wh
      enddo
      gam = gam*dnu
      hea = hea*dnu
      end subroutine planck_frequency_reference

      subroutine loaded_table_field(nfail)
      ! A synthetic table F(E) = A E^-2 on 41 logarithmic points from 3.0 to
      ! 300 eV, as read_sed leaves the rows of a file: e_sed_node ascending
      ! in energy, F_sed_node the undiluted flux per unit photon energy.
      ! The field the code returns between two rows must be that power law.
      integer, intent(inout) :: nfail
      integer, parameter :: npt = 41
      integer :: i
      real*8  :: A, E_test, ref, measured

      sp_type     = 'Load'
      is_PL_sed   = .false.
      is_monochr  = .false.
      do_read_sed = .true.
      A           = 1.0d3

      if (allocated(e_sed_node)) deallocate(e_sed_node)
      if (allocated(F_sed_node)) deallocate(F_sed_node)
      allocate(e_sed_node(npt), F_sed_node(npt))
      do i = 1, npt
         e_sed_node(i) = 3.0d0*(3.0d2/3.0d0)**(dble(i-1)/dble(npt-1))
         F_sed_node(i) = A*e_sed_node(i)**(-2.0d0)
      enddo

      E_test   = sqrt(e_sed_node(7)*e_sed_node(8))  ! inside one segment
      ref      = A*E_test**(-2.0d0)
      measured = stellar_flux_eV(E_test)
      call relative_verdict('loaded_table_field_interpolated', measured,    &
                            ref, 1.0d-12, nfail)
      end subroutine loaded_table_field

      subroutine table_reaching_the_grid_floor(nfail)
      ! A loaded table read down to the photon-grid floor e_low, as read_sed
      ! leaves one: its rows straddle the floor, so the LOWEST SELECTED row
      ! sits a fraction of a row ABOVE it (the selection stops at the first
      ! row below the floor). The band between the two is stated by the
      ! file, and the code must
      !   (a) call the floor covered, so that the Balmer integral is not
      !       refused for a file that does reach 3647 A, and
      !   (b) return the table's own field there, by continuing the first
      !       segment.
      ! The table is the same A E^-2 power law as assertion 5, on rows one
      ! part in 3600 apart, the spacing of the 1 A rows of the eps Eri
      ! spectrum at 3647 A; log-log continuation of a power law is exact, so
      ! the field is an identity, tolerance 1e-12.
      integer, intent(inout) :: nfail
      integer, parameter :: npt = 41
      integer :: i
      real*8  :: A, ref, measured, e_bot

      sp_type     = 'Load'
      is_PL_sed   = .false.
      is_monochr  = .false.
      do_read_sed = .true.
      A           = 1.0d3
      e_low       = e_th_HI_n2
      e_bot       = e_th_HI_n2*(1.0d0 + 1.0d0/3.6d3)

      if (allocated(e_sed_node)) deallocate(e_sed_node)
      if (allocated(F_sed_node)) deallocate(F_sed_node)
      allocate(e_sed_node(npt), F_sed_node(npt))
      do i = 1, npt
         e_sed_node(i) = e_bot*(1.24d3/e_bot)**(dble(i-1)/dble(npt-1))
         F_sed_node(i) = A*e_sed_node(i)**(-2.0d0)
      enddo

      call logical_verdict('balmer_band_covered_table_reaching_floor',     &
                           spectrum_covers_eV(e_th_HI_n2), .true., nfail)
      ref      = A*e_th_HI_n2**(-2.0d0)
      measured = stellar_flux_eV(e_th_HI_n2)
      call relative_verdict('loaded_table_field_at_the_grid_floor',        &
                            measured, ref, 1.0d-12, nfail)
      end subroutine table_reaching_the_grid_floor

      subroutine short_table(e_bot, e_up)
      ! A loaded table whose rows start at e_bot, i.e. above the n=2 edge.
      real*8, intent(in) :: e_bot, e_up
      integer, parameter :: npt = 11
      integer :: i
      if (allocated(e_sed_node)) deallocate(e_sed_node)
      if (allocated(F_sed_node)) deallocate(F_sed_node)
      allocate(e_sed_node(npt), F_sed_node(npt))
      do i = 1, npt
         e_sed_node(i) = e_bot*(e_up/e_bot)**(dble(i-1)/dble(npt-1))
         F_sed_node(i) = 1.0d0
      enddo
      end subroutine short_table

      subroutine relative_verdict(name, measured, reference, tol, nfail)
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: measured, reference, tol
      integer,          intent(inout) :: nfail
      real*8 :: dev
      dev = abs(measured - reference)/max(abs(reference), 1.0d-300)
      if (dev .le. tol) then
         write(*,'(A,ES22.15,A,ES22.15,A,ES9.2,A,ES10.3)') 'PASS '//name//  &
            ' measured=', measured, ' reference=', reference, ' tol=', tol, &
            ' relative deviation=', dev
      else
         write(*,'(A,ES22.15,A,ES22.15,A,ES9.2,A,ES10.3)') 'FAIL '//name//  &
            ' measured=', measured, ' reference=', reference, ' tol=', tol, &
            ' relative deviation=', dev
         nfail = nfail + 1
      endif
      end subroutine relative_verdict

      subroutine logical_verdict(name, measured, reference, nfail)
      character(len=*), intent(in)    :: name
      logical,          intent(in)    :: measured, reference
      integer,          intent(inout) :: nfail
      if (measured .eqv. reference) then
         write(*,'(A,L2,A,L2,A)') 'PASS '//name//' measured=', measured,    &
            ' reference=', reference, ' tol=0'
      else
         write(*,'(A,L2,A,L2,A)') 'FAIL '//name//' measured=', measured,    &
            ' reference=', reference, ' tol=0'
         nfail = nfail + 1
      endif
      end subroutine logical_verdict

      end program balmer_continuum_field
