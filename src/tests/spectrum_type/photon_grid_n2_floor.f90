      program photon_grid_n2_floor
      ! QUANTITY UNDER TEST
      !   the photon grid set_energy_vectors builds when the
      !   excited-hydrogen coupling is armed, against the grid of the same
      !   run with it off.  H(n=2) is an absorber like the metastable helium
      !   and the low-IP metals, so its threshold e_th_HI_n2 = 3.400 eV
      !   (3647 A) floors the grid; the band it opens, [3.400, 4.768] eV
      !   here, carries NO cross section of the grid, because every
      !   photoionization cross section of the code turns on at or above the
      !   lowest bulk-gas absorber edge e_abs_low.
      !
      ! WHY IT MATTERS
      !   The H(n=2) photoionization rate is an integral of stellar_flux_eV
      !   over its own band (excited_hydrogen), not a sum over these bins.
      !   The bins therefore only STATE the band, and lowering the floor must
      !   leave every rate the grid does compute exactly where it was.  This
      !   is the statement that the floor rule costs the analytic spectrum
      !   types nothing (development plan rev 3, section 10.5 decisions 13
      !   and 17).
      !
      ! ASSERTIONS (references and tolerances)
      !   1 sum(de_v) against e_top - e_th_HI_n2, 1e-12 relative: the grid
      !     is contiguous and it reaches the n=2 edge.
      !   2 the number of bins the floor adds, against num_n2 = 20 of
      !     set_energy_vectors, tolerance 0.
      !   3 e_v, de_v and F_XUV of every bin above e_abs_low, against the
      !     same bins of the armed-off grid, tolerance 0 (exact equality):
      !     the partition above the split is untouched.
      !   4 s_hi, s_hei, s_heii, s_heiTR, s_h2, the three H2 channels and
      !     every metal cross section on the added bins, against 0,
      !     tolerance 0.
      !   5 the four H and He photoionization rate integrals
      !     sum(F sigma/E de_v) of the two grids, tolerance 0 (exact
      !     equality).
      !   6 the same integral for every metal ion of an element metals.inp
      !     carries, tolerance 0.
      !
      ! CONFIGURATION: the stellar and spectral input of
      ! backup/regression/wasp_full/input.inp -- power law of index -1 on
      ! [13.60, 123.98, 1.24e3] eV, LX 29.46, LEUV 30.42, a 0.02544 AU,
      ! "2D approximate method: Rate/2 + Mdot/2" -- with the He 2^3S
      ! metastable and solar Na, so that e_abs_low is the metastable edge
      ! 4.768 eV and the band the n=2 floor opens is [3.400, 4.768] eV.
      !
      ! EXPECTED BEFORE THE 2c-N2FLOOR CHANGE: RED on assertions 1 and 2
      ! (the grid stopped at e_abs_low, so it added no bin and its width was
      ! e_top - e_th_HeTR).

      use global_parameters
      use species_table,            only: n_melem, n_mion, n_mphot,          &
                                          melem_i0, melem_top, mion_ethr,    &
                                          mion_iphot
      use energy_vectors_construct, only: set_energy_vectors
      use J_incident,               only: e_th_HI_n2

      implicit none

      integer :: nfail, nl_off, nl_on, n_added, i, k
      real*8, allocatable :: e_off(:), de_off(:), f_off(:)
      real*8 :: p_hi_off, p_hei_off, p_heii_off, p_hetr_off
      real*8 :: metal_rate_off(n_mphot)
      real*8 :: p_hi_on,  p_hei_on,  p_heii_on,  p_hetr_on
      real*8 :: dev, worst, width

      nfail = 0

      ! --- the run-wide input, as input_read resolves it for wasp_full ----
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
      ! Na is element 7 of the species table; its neutral threshold 5.139 eV
      ! is below the H I edge, which is what arms thereis_lowIP_metal in
      ! input_read.  The metastable edge 4.768 eV is lower still, so
      ! e_abs_low is the metastable's.
      melem_ab(7)         = 1.74d-6
      thereis_lowIP_metal = .true.
      thereis_HeITR       = .true.

      ! ---- the grid without the excited-hydrogen coupling ---------------
      use_excited_H = .false.
      call fresh_grid
      nl_off = Nl
      allocate(e_off(nl_off), de_off(nl_off), f_off(nl_off))
      e_off  = e_v
      de_off = de_v
      f_off  = F_XUV
      call photoionization_rates(p_hi_off, p_hei_off, p_heii_off, p_hetr_off)
      do k = 1, n_mphot
         metal_rate_off(k) = sum(F_XUV*sigma_tab(:,k)/e_v*de_v)
      enddo
      write(*,'(A,I0,A,F9.4,A,F9.4,A)') '  armed off: Nl=', nl_off,          &
           '  grid [', e_v(1), ',', e_v(Nl), '] eV'

      ! ---- the same run with it armed -----------------------------------
      use_excited_H = .true.
      call fresh_grid
      nl_on = Nl
      call photoionization_rates(p_hi_on, p_hei_on, p_heii_on, p_hetr_on)
      write(*,'(A,I0,A,F9.4,A,F9.4,A)') '  armed on : Nl=', nl_on,           &
           '  grid [', e_v(1), ',', e_v(Nl), '] eV'
      n_added = nl_on - nl_off

      ! ---- 1  the grid reaches the n=2 edge -----------------------------
      width = sum(de_v)
      call relative_verdict('n2_floor_grid_reaches_the_n2_edge', width,      &
                            e_top - e_th_HI_n2, 1.0d-12, nfail)

      ! ---- 2  how many bins it added ------------------------------------
      call integer_verdict('n2_floor_added_bin_count', n_added, 20, nfail)

      ! ---- 3  the bins above e_abs_low are untouched --------------------
      worst = 0.0d0
      if (n_added .ge. 0 .and. nl_on - n_added .eq. nl_off) then
         do i = 1, nl_off
            worst = max(worst, abs(e_v(n_added+i)  - e_off(i)))
            worst = max(worst, abs(de_v(n_added+i) - de_off(i)))
            worst = max(worst, abs(F_XUV(n_added+i)/max(f_off(i),1.0d-300)   &
                                   - 1.0d0))
         enddo
      else
         worst = huge(1.0d0)
      endif
      call absolute_verdict('n2_floor_bins_above_the_absorber_edge_identical',&
                            worst, 0.0d0, 0.0d0, nfail)

      ! ---- 4  nothing ACTIVE absorbs on the added bins ------------------
      ! The metal cross-section table is filled for all n_mphot ions
      ! whatever metals.inp carries, so only the columns of an element with
      ! a non-zero abundance are absorbers of this run. K I (4.341 eV) sits
      ! in the table below the metastable edge and is non-zero here, but it
      ! is inert unless K is present -- and if it were present it would be
      ! the lowest bulk-gas threshold itself, so the split would move down
      ! to it and these bins would again carry nothing.
      worst = 0.0d0
      do i = 1, n_added
         worst = max(worst, abs(s_hi(i)),   abs(s_hei(i)))
         worst = max(worst, abs(s_heii(i)), abs(s_heiTR(i)))
         worst = max(worst, abs(s_h2(i)),   abs(s_h2_di(i)))
         worst = max(worst, abs(s_h2_dd(i)), abs(s_h2_nd(i)))
         do k = 1, n_mion
            if (mion_iphot(k) .gt. 0 .and. ion_is_active(k))                 &
               worst = max(worst, abs(sigma_tab(i,mion_iphot(k))))
         enddo
      enddo
      call absolute_verdict(                                                 &
           'n2_floor_cross_sections_vanish_below_the_absorber_edge',         &
           worst, 0.0d0, 0.0d0, nfail)

      ! ---- 5  no rate the grid computes moved ---------------------------
      worst = 0.0d0
      worst = max(worst, reldev(p_hi_on,   p_hi_off))
      worst = max(worst, reldev(p_hei_on,  p_hei_off))
      worst = max(worst, reldev(p_heii_on, p_heii_off))
      worst = max(worst, reldev(p_hetr_on, p_hetr_off))
      write(*,'(A,4(1X,ES14.7))') '     P_HI, P_HeI, P_HeII, P_HeTR [s^-1]'//&
           ' armed on:', p_hi_on, p_hei_on, p_heii_on, p_hetr_on
      call absolute_verdict('n2_floor_photoionization_rates_unchanged',      &
                            worst, 0.0d0, 0.0d0, nfail)

      ! ---- 6  nor did any metal rate of an active element ---------------
      worst = 0.0d0
      do k = 1, n_mion
         if (mion_iphot(k) .le. 0) cycle
         if (.not. ion_is_active(k)) cycle
         worst = max(worst, reldev(                                          &
            sum(F_XUV*sigma_tab(:,mion_iphot(k))/e_v*de_v),                  &
            metal_rate_off(mion_iphot(k))))
      enddo
      call absolute_verdict('n2_floor_metal_rates_unchanged',                &
                            worst, 0.0d0, 0.0d0, nfail)

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'photon_grid_n2_floor: ', nfail,                &
                             ' assertion(s) FAILED'
         call exit(1)
      endif
      write(*,'(A)') 'photon_grid_n2_floor: all assertions PASSED'

      contains

      subroutine fresh_grid
      ! set_energy_vectors allocates the grid and the cross-section vectors
      ! unconditionally, so a second call needs a clean slate.
      if (allocated(e_v))       deallocate(e_v)
      if (allocated(de_v))      deallocate(de_v)
      if (allocated(F_XUV))     deallocate(F_XUV)
      if (allocated(sigma_tab)) deallocate(sigma_tab)
      if (allocated(s_h2_di))   deallocate(s_h2_di)
      if (allocated(s_h2_dd))   deallocate(s_h2_dd)
      if (allocated(s_h2_nd))   deallocate(s_h2_nd)
      call set_energy_vectors
      end subroutine fresh_grid

      subroutine photoionization_rates(p_hi, p_hei, p_heii, p_hetr)
      ! The four rate integrals of the grid as util_ion_eq forms them at zero
      ! optical depth, sum(F sigma/E de_v) [s^-1] (the 1e-18 of the cross
      ! sections and the eV-to-erg factor are the same in both grids and are
      ! kept, so the numbers are the rates themselves).
      real*8, intent(out) :: p_hi, p_hei, p_heii, p_hetr
      p_hi   = sum(F_XUV*s_hi   /e_v*de_v)*1.0d-18*erg2eV
      p_hei  = sum(F_XUV*s_hei  /e_v*de_v)*1.0d-18*erg2eV
      p_heii = sum(F_XUV*s_heii /e_v*de_v)*1.0d-18*erg2eV
      p_hetr = sum(F_XUV*s_heiTR/e_v*de_v)*1.0d-18*erg2eV
      end subroutine photoionization_rates

      logical function ion_is_active(iion)
      ! Whether ion iion belongs to an element metals.inp gave a non-zero
      ! abundance, i.e. whether it is an absorber of this run at all.
      integer, intent(in) :: iion
      integer :: je
      ion_is_active = .false.
      do je = 1, n_melem
         if (melem_ab(je) .le. 0.0d0) cycle
         if (iion .ge. melem_i0(je) .and. iion .le. melem_top(je))           &
            ion_is_active = .true.
      enddo
      end function ion_is_active

      double precision function reldev(a, b)
      real*8, intent(in) :: a, b
      reldev = abs(a - b)/max(abs(b), 1.0d-300)
      end function reldev

      subroutine relative_verdict(name, measured, reference, tol, nfail)
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: measured, reference, tol
      integer,          intent(inout) :: nfail
      dev = abs(measured - reference)/max(abs(reference), 1.0d-300)
      if (dev .le. tol) then
         write(*,'(A,ES22.15,A,ES22.15,A,ES9.2,A,ES10.3)') 'PASS '//name//   &
            ' measured=', measured, ' reference=', reference, ' tol=', tol,  &
            ' relative deviation=', dev
      else
         write(*,'(A,ES22.15,A,ES22.15,A,ES9.2,A,ES10.3)') 'FAIL '//name//   &
            ' measured=', measured, ' reference=', reference, ' tol=', tol,  &
            ' relative deviation=', dev
         nfail = nfail + 1
      endif
      end subroutine relative_verdict

      subroutine absolute_verdict(name, measured, reference, tol, nfail)
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: measured, reference, tol
      integer,          intent(inout) :: nfail
      if (abs(measured - reference) .le. tol) then
         write(*,'(A,ES22.15,A,ES22.15,A,ES9.2)') 'PASS '//name//            &
            ' measured=', measured, ' reference=', reference, ' tol=', tol
      else
         write(*,'(A,ES22.15,A,ES22.15,A,ES9.2)') 'FAIL '//name//            &
            ' measured=', measured, ' reference=', reference, ' tol=', tol
         nfail = nfail + 1
      endif
      end subroutine absolute_verdict

      subroutine integer_verdict(name, measured, reference, nfail)
      character(len=*), intent(in)    :: name
      integer,          intent(in)    :: measured, reference
      integer,          intent(inout) :: nfail
      if (measured .eq. reference) then
         write(*,'(A,I0,A,I0,A)') 'PASS '//name//' measured=', measured,     &
            ' reference=', reference, ' tol=0'
      else
         write(*,'(A,I0,A,I0,A)') 'FAIL '//name//' measured=', measured,     &
            ' reference=', reference, ' tol=0'
         nfail = nfail + 1
      endif
      end subroutine integer_verdict

      end program photon_grid_n2_floor
