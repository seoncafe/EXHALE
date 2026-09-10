      program photon_grid_quadrature
      ! QUANTITY UNDER TEST
      !   the photoionization rates of one unattenuated cell,
      !     P_X = 1e-18 * erg2eV * INT F(E) sigma_X(E) / E  dE
      !   for X = H I, He I, He II, and the H I photoheating rate of one atom,
      !     h1_HI = 1e-18 * INT F(E) (1 - E_th/E) sigma_HI(E) dE ,
      !   as the production code forms them, against a fine reference
      !   integral of the SAME integrand.
      !
      ! HOW THE CODE FORMS THEM (util_ion_eq.f90, ionization_equilibrium
      ! branch). At zero optical depth int_f = F_XUV, and
      !   int_1  = int_f*s_hi/e_v          (770-772)
      !   PIR_1  = sum(int_1*de_v)         (1015-1017)
      !   P_HI   = PIR_1*1.0e-18*erg2eV    (1055)
      !   acc_HI = photoelectron_share(e_th_HI,e_v)*fhv*s_hi   (826, fhv = 1
      !            with the secondary-ionization coupling off)
      !   h1_HI  = sum(int_f*acc_HI*de_v)*1.0e-18              (1043)
      ! This program evaluates those same expressions on the arrays that the
      ! PRODUCTION routine set_energy_vectors filled (e_v, de_v, F_XUV, s_hi,
      ! s_hei, s_heii) and with the PRODUCTION photoelectron_share.
      !
      ! REFERENCE. The code's sum is a quadrature of the same integrand
      ! over the photon band of the run, which is [the grid floor, e_top]:
      ! the floor is the lowest ionization threshold of an active absorber
      ! (4.80 eV for the He 2^3S metastable, the lowest neutral-metal edge
      ! with a low-IP metal, 13.6 eV with neither) and e_top is the X-ray
      ! top J_inc normalizes over. The reference is that same integrand,
      ! from the same production functions J_inc and photoion_sigma,
      ! integrated by the trapezoid rule over that band on 400001
      ! logarithmically spaced points. Nothing is re-implemented: the
      ! reference differs from the code only in the number of quadrature
      ! points and in carrying no bin that straddles a threshold. A second
      ! reference on 200001 points is reported so that the fine integral can
      ! be seen to be converged.
      !
      ! TOLERANCE 1e-3 relative, the tolerance section 10.1 item 2 of
      ! docs/development_plan_20260905_rev3.md states for this gate.
      !
      ! Whether the code's own weights span that whole band is a separate
      ! statement about the grid and is reported as a DIAGNOSTIC (sum(de_v)
      ! against e_top minus the floor), not folded into the assertions on
      ! the quadrature.
      !
      ! CONFIGURATIONS (the grid floor of set_energy_vectors.f90 50-73):
      !   1  Include He23S? True             -> floor e_th_HeTR = 4.80 eV
      !   2  no He 2^3S, Na I active         -> floor mion_ethr(NaI) = 5.139 eV
      !   3  neither                         -> floor e_th_HI = 13.6 eV
      ! The spectrum is the one of backup/regression/wasp_full/input.inp:
      ! Power-law, index -1, [13.60, 123.98, 1.24e3] eV, LX 29.46,
      ! LEUV 30.42, a_orb 0.02544 AU, "2D approximate method: Rate/2 + Mdot/2".

      use global_parameters
      use species_table,          only: n_melem, melem_i0, mion_ethr
      use energy_vectors_construct, only: set_energy_vectors
      use J_incident,             only: J_inc
      use opacity_models,         only: photoion_sigma
      use utils_ion_eq,           only: photoelectron_share

      implicit none

      integer :: nfail
      ! The photon band of the configuration under test: the grid floor (the
      ! lowest ionization threshold of an active absorber) and e_top.  Host
      ! variables so that the reference integral and the coverage diagnostic
      ! read the same band the code is asked to cover.
      real*8  :: e_grid_lo, e_grid_hi
      nfail = 0

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
      a_tau        = 0.0d0
      allocate(melem_ab(n_melem))
      melem_ab     = 0.0d0

      call one_configuration(1, 'He23S_on_floor_4.80eV',  nfail)
      call one_configuration(2, 'NaI_metal_floor_5.139eV', nfail)
      call one_configuration(3, 'no_sub_Lyman_floor_13.6eV', nfail)

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'photon_grid_quadrature: ', nfail,              &
                             ' assertion(s) FAILED'
         call exit(1)
      endif
      write(*,'(A)') 'photon_grid_quadrature: all assertions PASSED'

      contains

      subroutine one_configuration(icfg, label, nfail)
      integer,          intent(in)    :: icfg
      character(len=*), intent(in)    :: label
      integer,          intent(inout) :: nfail

      real*8 :: p_hi_code, p_hei_code, p_heii_code, h1_hi_code
      real*8 :: p_hi_ref,  p_hei_ref,  p_heii_ref,  h1_hi_ref
      real*8 :: p_hi_ref2, p_hei_ref2, p_heii_ref2, h1_hi_ref2

      ! Free what a previous configuration allocated; set_energy_vectors
      ! allocates these unconditionally.
      if (allocated(e_v))       deallocate(e_v)
      if (allocated(de_v))      deallocate(de_v)
      if (allocated(F_XUV))     deallocate(F_XUV)
      if (allocated(sigma_tab)) deallocate(sigma_tab)
      if (allocated(s_h2_di))   deallocate(s_h2_di)
      if (allocated(s_h2_dd))   deallocate(s_h2_dd)
      if (allocated(s_h2_nd))   deallocate(s_h2_nd)

      melem_ab            = 0.0d0
      thereis_HeITR       = .false.
      thereis_lowIP_metal = .false.
      e_grid_lo           = e_th_HI
      select case (icfg)
         case (1)
            thereis_HeITR = .true.
            e_grid_lo     = e_th_HeTR
         case (2)
            ! Na is element 7 of the species table; melem_i0(7) is the
            ! neutral stage and mion_ethr of it is 5.139 eV, below the
            ! 13.6 eV edge, which is what arms thereis_lowIP_metal in
            ! input_read.f90 1248-1252.
            melem_ab(7)         = 1.0d-6
            thereis_lowIP_metal = .true.
            e_grid_lo           = mion_ethr(melem_i0(7))
      end select
      e_grid_hi = e_top

      call set_energy_vectors

      write(*,'(A,A,A,I0,A,F9.4,A,F9.4,A)') '  configuration ', label,       &
           ': Nl=', Nl, '  grid [', e_v(1), ',', e_v(Nl), '] eV'
      write(*,'(A,I0,A,F9.4,A)') '     NlTR=', NlTR,                         &
           '  lowest metal threshold used=', mion_ethr(melem_i0(7)), ' eV'

      ! --- the code's own quadrature -------------------------------------
      p_hi_code   = sum(F_XUV*s_hi  /e_v*de_v)*1.0d-18*erg2eV
      p_hei_code  = sum(F_XUV*s_hei /e_v*de_v)*1.0d-18*erg2eV
      p_heii_code = sum(F_XUV*s_heii/e_v*de_v)*1.0d-18*erg2eV
      h1_hi_code  = sum(F_XUV*photoelectron_share(e_th_HI,e_v)*s_hi*de_v)    &
                    *1.0d-18

      ! --- the fine reference of the same integrand ----------------------
      call fine_reference(400001, p_hi_ref,  p_hei_ref,  p_heii_ref,         &
                          h1_hi_ref)
      call fine_reference(200001, p_hi_ref2, p_hei_ref2, p_heii_ref2,        &
                          h1_hi_ref2)
      write(*,'(A,4(1X,ES10.3))')                                            &
           '     reference self-convergence (400001 vs 200001, relative):',  &
           abs(p_hi_ref-p_hi_ref2)/p_hi_ref,                                 &
           abs(p_hei_ref-p_hei_ref2)/p_hei_ref,                              &
           abs(p_heii_ref-p_heii_ref2)/p_heii_ref,                           &
           abs(h1_hi_ref-h1_hi_ref2)/h1_hi_ref

      call grid_band_coverage

      call ratio_verdict('P_HI['  //label//']', p_hi_code,   p_hi_ref,  nfail)
      call ratio_verdict('P_HeI[' //label//']', p_hei_code,  p_hei_ref, nfail)
      call ratio_verdict('P_HeII['//label//']', p_heii_code, p_heii_ref,nfail)
      call ratio_verdict('heat_HI['//label//']',h1_hi_code,  h1_hi_ref, nfail)

      end subroutine one_configuration

      subroutine fine_reference(npt, p_hi, p_hei, p_heii, h1_hi)
      ! Trapezoid rule on npt logarithmic points over the photon band of the
      ! run, [e_grid_lo, e_grid_hi], of the SAME integrand, from the same
      ! production functions.
      integer, intent(in)  :: npt
      real*8,  intent(out) :: p_hi, p_hei, p_heii, h1_hi
      integer :: i
      real*8  :: e, w, dlog, xi_day, flux, s1, s15, s2

      dlog   = log(e_grid_hi/e_grid_lo)/dble(npt-1)
      xi_day = dayside_dilution()
      p_hi   = 0.0d0
      p_hei  = 0.0d0
      p_heii = 0.0d0
      h1_hi  = 0.0d0
      do i = 1, npt
         e = e_grid_lo*exp(dlog*dble(i-1))
         ! trapezoid weight in E of a logarithmic abscissa
         w = e*dlog
         if (i .eq. 1 .or. i .eq. npt) w = 0.5d0*w
         flux = xi_day*J_inc(e)
         s1   = photoion_sigma('HI',   e)
         s15  = photoion_sigma('HeI',  e)
         s2   = photoion_sigma('HeII', e)
         p_hi   = p_hi   + flux*s1 /e*w
         p_hei  = p_hei  + flux*s15/e*w
         p_heii = p_heii + flux*s2 /e*w
         h1_hi  = h1_hi  + flux*photoelectron_share(e_th_HI,e)*s1*w
      enddo
      p_hi   = p_hi  *1.0d-18*erg2eV
      p_hei  = p_hei *1.0d-18*erg2eV
      p_heii = p_heii*1.0d-18*erg2eV
      h1_hi  = h1_hi *1.0d-18
      end subroutine fine_reference

      subroutine grid_band_coverage
      ! DIAGNOSTIC, not an assertion: whether the code's own weights span
      ! the photon band [e_grid_lo, e_grid_hi].  With contiguous bins
      ! sum(de_v) is the band width exactly; with the central-difference
      ! weights of a grid whose points sit on the band boundaries it
      ! telescopes to e_v(Nl) - e_v(1), which is short of the band at both
      ! ends.
      real*8 :: cover, band
      cover = sum(de_v)
      band  = e_grid_hi - e_grid_lo
      write(*,'(A,F9.3,A,F9.3,A,ES10.3)')                                    &
           '     DIAGNOSTIC sum(de_v) = ', cover, ' eV, band width ', band,  &
           ' eV, relative shortfall ', (band - cover)/band
      end subroutine grid_band_coverage

      subroutine ratio_verdict(name, code_value, ref_value, nfail)
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: code_value, ref_value
      integer,          intent(inout) :: nfail
      real*8, parameter :: tol = 1.0d-3
      real*8 :: ratio
      ratio = code_value/ref_value
      if (abs(ratio - 1.0d0) .le. tol) then
         write(*,'(A,A,A,ES22.15,A,ES22.15,A,ES9.2,A,F10.6)') 'PASS ', name, &
              ' measured=', code_value, ' reference=', ref_value,            &
              ' tol=', tol, ' code/exact=', ratio
      else
         write(*,'(A,A,A,ES22.15,A,ES22.15,A,ES9.2,A,F10.6)') 'FAIL ', name, &
              ' measured=', code_value, ' reference=', ref_value,            &
              ' tol=', tol, ' code/exact=', ratio
         nfail = nfail + 1
      endif
      end subroutine ratio_verdict

      end program photon_grid_quadrature
