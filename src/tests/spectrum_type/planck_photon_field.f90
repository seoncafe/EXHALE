      program planck_photon_field
      ! QUANTITY UNDER TEST
      !   the incident flux F_XUV(E) that set_energy_vectors puts on the
      !   photon grid for "Spectrum type: Planck", i.e. the photospheric
      !   field of a star of effective temperature T_eff and radius R_star
      !   seen at the orbital distance a,
      !
      !     F_E(E) = xi * pi B_nu(T_eff) (R_star/a)^2 * (dnu/dE) ,
      !     B_nu   = 2 h nu^3/c^2 / (exp(h nu / k T_eff) - 1) ,
      !     nu     = E * eV2Hz ,  dnu/dE = eV2Hz ,
      !
      !   in erg cm^-2 s^-1 eV^-1, xi the run's dayside dilution.  This is
      !   decision 13 of docs/development_plan_20260905_rev3.md section 10.5:
      !   ONE spectrum type builds every band of the grid, and the Planck
      !   type is the photospheric alternative to the extrapolated power law
      !   below 13.6 eV, where the He 2^3S metastable and the low-IP metals
      !   absorb.
      !
      ! REFERENCES OF THE FOUR ASSERTIONS
      !   1 grid       the Planck type builds the SAME photon grid as the
      !                power law (same floor rule, same points): e_v and de_v
      !                are compared point by point, tolerance 0 (identity).
      !   2 field      F_XUV(i) against the expression above, re-evaluated
      !                here from the production constants, tolerance 1e-12
      !                relative (an identity between two writings of one
      !                formula, so round-off is the only allowance).
      !   3 total      the frequency integral of the same field, INT pi B_nu
      !                (R/a)^2 dnu over 0 < x < 60, against the Stefan-
      !                Boltzmann law sigma T_eff^4 (R_star/a)^2 with
      !                sigma = 5.670374419e-5 erg cm^-2 s^-1 K^-4 (CODATA
      !                2018), tolerance 1e-6 relative.  This is the
      !                independent check of the normalization and of the
      !                units of the shared Planck function; the field of
      !                assertion 2 is that function times dnu/dE.
      !   4 admissible every F_XUV(i) finite and >= 0 on the whole grid,
      !                tolerance 0.  The Wien tail of a 6459 K star at the
      !                1.24 keV top of the grid has h nu / k T = 2228, so a
      !                writing that forms exp(+x) before dividing returns
      !                Inf/NaN there; this assertion is that gate.
      !
      ! CONFIGURATION: the stellar and spectral input of
      ! backup/regression/wasp_full/input.inp -- T_eff 6459 K,
      ! R_star 1.458 Rsun, a 0.02544 AU, power law index -1 on
      ! [13.60, 123.98, 1.24e3] eV, LX 29.46, LEUV 30.42, "2D approximate
      ! method: Rate/2 + Mdot/2" (xi = 1/2), He 2^3S on, so that the grid
      ! carries the sub-Lyman points down to 4.80 eV.
      !
      ! EXPECTED BEFORE THE CHANGE: RED (the program does not build; there is
      ! no Planck branch and no shared Planck function).  After: GREEN.

      use global_parameters
      use species_table,            only: n_melem
      use energy_vectors_construct, only: set_energy_vectors
      use J_incident,               only: J_inc, planck_stellar_flux_nu,    &
                                          eV2Hz

      implicit none

      ! Stefan-Boltzmann constant, CODATA 2018 (2 pi^5 k^4 / (15 c^2 h^3)),
      ! erg cm^-2 s^-1 K^-4.  The published reference of assertion 3.
      real*8, parameter :: sigma_SB = 5.670374419d-5

      integer :: nfail, i
      real*8, dimension(:), allocatable :: e_pl, de_pl
      real*8  :: xi, R_over_a, ref, dev, dev_max, e_dev
      real*8  :: f_int, f_sb, x_lo, x_hi, nu_a, nu_b, lg, nu, w
      integer :: n_nu, i_nu
      logical :: ok

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
      thereis_HeITR       = .true.
      thereis_lowIP_metal = .false.
      allocate(melem_ab(n_melem))
      melem_ab     = 0.0d0

      xi       = dayside_dilution()
      R_over_a = R_star/a_orb

      ! ---- the power-law grid, for the grid identity of assertion 1 ------
      sp_type   = 'Power-law'
      is_PL_sed = .true.
      call clear_spectral_arrays
      call set_energy_vectors
      allocate(e_pl(Nl), de_pl(Nl))
      e_pl  = e_v
      de_pl = de_v

      ! ---- the Planck grid and field -------------------------------------
      sp_type   = 'Planck'
      is_PL_sed = .false.
      call clear_spectral_arrays
      call set_energy_vectors

      ! --- 1 the grid is the power law's, point for point -----------------
      if (Nl .ne. size(e_pl)) then
         write(*,'(A,I0,A,I0,A)') 'FAIL planck_grid_points measured=', Nl,  &
            ' reference=', size(e_pl), ' tol=0'
         nfail = nfail + 1
      else
         dev_max = 0.0d0
         do i = 1,Nl
            dev_max = max(dev_max, abs(e_v(i) - e_pl(i)),                   &
                                   abs(de_v(i) - de_pl(i)))
         enddo
         call verdict('planck_grid_identical_to_power_law', dev_max,        &
                      0.0d0, 0.0d0, nfail)
      endif

      ! --- 2 the field is pi B_nu (R/a)^2 dnu/dE, diluted ------------------
      dev_max = 0.0d0
      e_dev   = 0.0d0
      do i = 1,Nl
         ref = xi*planck_flux_reference(T_star_eff, R_over_a, e_v(i))
         if (ref .gt. 0.0d0) then
            dev = abs(F_XUV(i) - ref)/ref
         else
            dev = abs(F_XUV(i) - ref)
         endif
         if (dev .gt. dev_max) then
            dev_max = dev
            e_dev   = e_v(i)
         endif
      enddo
      write(*,'(A,ES12.5,A)') '  largest deviation at E = ', e_dev, ' eV'
      call verdict('planck_field_on_grid', dev_max, 0.0d0, 1.0d-12, nfail)

      ! --- 3 the frequency integral is sigma T^4 (R/a)^2 -------------------
      ! Trapezoid on a logarithmic frequency grid over 1e-6 <= h nu / k T
      ! <= 60; the tail beyond x = 60 is 1e-24 of the total.
      x_lo  = 1.0d-6
      x_hi  = 6.0d1
      nu_a  = x_lo*kb_erg*T_star_eff/hp_erg
      nu_b  = x_hi*kb_erg*T_star_eff/hp_erg
      n_nu  = 400001
      lg    = log(nu_b/nu_a)/dble(n_nu-1)
      f_int = 0.0d0
      do i_nu = 1,n_nu
         nu = nu_a*exp(dble(i_nu-1)*lg)
         w  = nu*lg                              ! dnu = nu dln(nu)
         if (i_nu .eq. 1 .or. i_nu .eq. n_nu) w = 0.5d0*w
         f_int = f_int + planck_stellar_flux_nu(T_star_eff, R_over_a, nu)*w
      enddo
      f_sb = sigma_SB*T_star_eff**4*R_over_a**2
      write(*,'(A,ES14.7,A,ES14.7)') '  integral = ', f_int,               &
         ' erg cm^-2 s^-1   Stefan-Boltzmann = ', f_sb
      call verdict('planck_flux_stefan_boltzmann', abs(f_int-f_sb)/f_sb,    &
                   0.0d0, 1.0d-6, nfail)

      ! --- 4 finite and non-negative over the whole grid --------------------
      ok = .true.
      do i = 1,Nl
         if (F_XUV(i) .lt. 0.0d0) ok = .false.
         if (.not. (F_XUV(i) .eq. F_XUV(i))) ok = .false.      ! NaN
         if (F_XUV(i) .gt. huge(1.0d0)) ok = .false.           ! Inf
      enddo
      if (ok) then
         write(*,'(A,I0,A,I0,A)') 'PASS planck_field_admissible measured=', &
            Nl, ' reference=', Nl, ' tol=0'
         write(*,'(A)') '   (grid points finite and >= 0)'
      else
         write(*,'(A)') 'FAIL planck_field_admissible measured=not_finite'//&
            ' reference=finite_nonnegative tol=0'
         nfail = nfail + 1
      endif

      ! --- context: what the two types mean for the metastable ------------
      ! The unattenuated He 2^3S photoionization rate the code forms,
      !   P = 1e-18 erg2eV sum(F sigma_HeITR / E dE)   (util_ion_eq),
      ! under the two spectrum types on the same grid and cross section.
      f_int = 0.0d0
      f_sb  = 0.0d0
      do i = 1,Nl
         f_int = f_int + F_XUV(i)*s_heiTR(i)/e_v(i)*de_v(i)
         f_sb  = f_sb  + xi*J_inc(e_v(i))*s_heiTR(i)/e_v(i)*de_v(i)
      enddo
      f_int = f_int*1.0d-18*erg2eV
      f_sb  = f_sb *1.0d-18*erg2eV
      write(*,'(A,ES11.4,A,ES11.4,A,ES10.3)')                              &
         'DIAGNOSTIC He 2^3S photoionization rate [s^-1]: Planck ', f_int,  &
         '  power law ', f_sb, '  ratio ', f_int/max(f_sb,1.0d-300)

      write(*,'(A)') 'DIAGNOSTIC Planck / power law at four energies:'
      do i = 1,Nl
         if (e_v(i) .lt. 4.81d0 .or.                                        &
             (e_v(i) .gt. 5.99d0 .and. e_v(i) .lt. 6.30d0) .or.             &
             (e_v(i) .gt. 9.90d0 .and. e_v(i) .lt. 10.4d0) .or.             &
             (e_v(i) .gt. 13.0d0 .and. e_v(i) .lt. 13.7d0)) then
            write(*,'(A,F8.3,A,ES11.4,A,ES11.4,A,ES10.3)')                  &
               '   E = ', e_v(i), ' eV  Planck ', F_XUV(i),                 &
               '  power law ', xi*J_inc(e_v(i)),                            &
               '  ratio ', F_XUV(i)/max(xi*J_inc(e_v(i)),1.0d-300)
         endif
      enddo

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'planck_photon_field: ', nfail,                &
                             ' assertion(s) FAILED'
         call exit(1)
      endif
      write(*,'(A)') 'planck_photon_field: all assertions PASSED'

      contains

      subroutine clear_spectral_arrays
      ! set_energy_vectors allocates these unconditionally.
      if (allocated(e_v))       deallocate(e_v)
      if (allocated(de_v))      deallocate(de_v)
      if (allocated(F_XUV))     deallocate(F_XUV)
      if (allocated(sigma_tab)) deallocate(sigma_tab)
      if (allocated(s_h2_di))   deallocate(s_h2_di)
      if (allocated(s_h2_dd))   deallocate(s_h2_dd)
      if (allocated(s_h2_nd))   deallocate(s_h2_nd)
      end subroutine clear_spectral_arrays

      double precision function planck_flux_reference(Tstar, R_a, E)
      ! The field of decision 13, written out here independently of the
      ! production function: F_E = pi B_nu (R/a)^2 dnu/dE with nu = E*eV2Hz.
      real*8, intent(in) :: Tstar, R_a, E
      real*8 :: nu, x, Bnu
      nu  = E*eV2Hz
      x   = hp_erg*nu/(kb_erg*Tstar)
      if (x .gt. 7.0d2) then
         Bnu = (2.0d0*hp_erg*nu**3.0/c_light**2.0)*exp(-x)
      else
         Bnu = (2.0d0*hp_erg*nu**3.0/c_light**2.0)/(exp(x) - 1.0d0)
      endif
      planck_flux_reference = pi*Bnu*R_a**2.0*eV2Hz
      end function planck_flux_reference

      subroutine verdict(name, measured, reference, tol, nfail)
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: measured, reference, tol
      integer,          intent(inout) :: nfail
      if (abs(measured - reference) .le. tol) then
         write(*,'(A,ES12.5,A,ES12.5,A,ES9.2)') 'PASS '//name//             &
            ' measured=', measured, ' reference=', reference, ' tol=', tol
      else
         write(*,'(A,ES12.5,A,ES12.5,A,ES9.2)') 'FAIL '//name//             &
            ' measured=', measured, ' reference=', reference, ' tol=', tol
         nfail = nfail + 1
      endif
      end subroutine verdict

      end program planck_photon_field
