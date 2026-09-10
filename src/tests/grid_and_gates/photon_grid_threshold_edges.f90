      program photon_grid_threshold_edges
      ! QUANTITY UNDER TEST
      !   (a) the position of every ionization threshold of an ACTIVE
      !       absorber relative to the bin partition of the photon grid, and
      !   (b) the metal photoionization rates of one unattenuated cell,
      !         P_m = 1e-18 * erg2eV * INT F(E) sigma_m(E) / E  dE ,
      !       as the production code forms them, against a fine reference
      !       integral of the SAME integrand.
      !
      ! WHY (a) IS AN EXACT IDENTITY AND NOT A TOLERANCE. A photoionization
      ! rate is a rectangle rule sum(F sigma/E de_v) over the bins. A bin
      ! that CONTAINS a threshold charges that absorber's cross section, read
      ! at the bin centre, over the whole bin -- including the part of the bin
      ! that lies BELOW the threshold, where the absorber cannot absorb at
      ! all. The only partition free of that error is the one in which every
      ! threshold of an active absorber is a bin EDGE, so that no bin
      ! straddles a turn-on of any cross section the sum contracts with.
      ! Item 2c-QUAD imposed this for H I, He I and He II; the thresholds
      ! tested here are the remaining ones -- the He 2^3S metastable when a
      ! lower metal floors the grid, H2, and each active metal ion.
      !
      ! HOW THE CODE FORMS THE RATES (util_ion_eq.f90 1150-1157): at zero
      ! optical depth int_f = F_XUV and
      !   int_m     = int_f*sigma_tab(:,k)/e_v
      !   Pm_loc(i) = sum(int_m*de_v)*1.0e-18*erg2eV
      ! with k = mion_iphot(i). This program evaluates that same expression on
      ! the arrays the PRODUCTION routine set_energy_vectors filled.
      !
      ! REFERENCE for (b). The same integrand from the same production
      ! functions J_inc and metal_photoion_sigma, integrated by the trapezoid
      ! rule on 400001 logarithmically spaced points over [threshold, e_top],
      ! which is the support of the integrand inside the run's photon band.
      ! Nothing is re-implemented: the reference differs from the code only in
      ! the number of quadrature points and in carrying no bin that straddles
      ! the threshold. A 200001-point integral is reported so that the fine
      ! one can be seen to be converged.
      !
      ! TOLERANCES. (a) 1e-12 relative: a bin edge either IS the threshold or
      ! it is not, and the grid is built from the constant itself, so the only
      ! admissible difference is the round-off of reconstructing the edges
      ! from the widths.
      !
      ! (b) 1e-3 relative, the tolerance section 10.1 item 2 of
      ! docs/development_plan_20260905_rev3.md states for this gate, for every
      ! ion. The four NEUTRAL metals whose thresholds lie in the sub-Lyman
      ! band [e_abs_low, 13.6 eV] are the hardest: their cross sections fall
      ! by up to a factor 2.5 over 0.4 eV above threshold (Fe I: the Verner
      ! et al. 1996 fit has E_0 = 0.0546 eV), so what is left once every
      ! threshold is an edge is the rectangle-rule error of that band's
      ! partition, bounded by its bin budget num_TR. MEASURED with this
      ! program, code/exact against num_TR:
      !          num_TR    20        30        40        50        80        89
      !   Mg I          0.997527  0.998711  0.999150  0.999550  0.999824  0.999853
      !   Ca I          0.997368  0.998772  0.999369  0.999612  0.999840  0.999870
      !   Na I          0.997523  0.999122  0.999582  0.999666  0.999869  0.999902
      !   Fe I          0.994172  0.996849  0.998048  0.998825  0.999525  0.999615
      ! num_TR = 89 gives the band the bin ratio of the pure-H I band and puts
      ! all four inside 5e-4, so the 1e-3 gate applies to them as well.
      !
      ! THE EDGES ARE RECONSTRUCTED, not read: set_energy_vectors deallocates
      ! its edge array. The bins are contiguous and de_v is the exact width,
      ! so edge(1) = photon_grid_floor_eV() and edge(j+1) = edge(j) + de_v(j)
      ! rebuild the partition to round-off.
      !
      ! CONFIGURATIONS, both from backup/regression:
      !   wasp_full   Power-law, index -1, [13.60, 123.98, 1.24e3] eV,
      !               LX 29.46, LEUV 30.42, a_orb 0.02544 AU,
      !               "Include He23S? True", stellar Teff and radius set (so
      !               use_excited_H is armed and the grid floor is the H(n=2)
      !               edge), metals.inp = solar C, N, O, Mg, Ca, Na, Fe.
      !   mol_metals  the same metals and He 2^3S on the hot Uranus,
      !               LX 27.20, LEUV 27.93, a_orb 0.0480 AU, no stellar
      !               lines, "Molecular chemistry: True" -- so H2 is an
      !               active absorber and its 15.4259 eV threshold must be an
      !               edge too.

      use global_parameters
      use species_table,            only: n_melem, n_mion, mion_ethr,        &
                                          mion_iphot, mion_isphot,           &
                                          mion_elem, mion_name,              &
                                          iel_C, iel_N, iel_O, iel_Mg,       &
                                          iel_Ca, iel_Na, iel_Fe
      use sed_reader,               only: photon_grid_floor_eV
      use energy_vectors_construct, only: set_energy_vectors
      use J_incident,               only: J_inc
      use Cross_sections,           only: metal_photoion_sigma

      implicit none

      integer :: nfail
      nfail = 0

      ! --- the run-wide input both configurations share -------------------
      is_PL_sed    = .true.
      do_read_sed  = .false.
      is_monochr   = .false.
      PLind        = -1.0d0
      thereis_Xray = .true.
      e_low        = 13.60d0
      e_mid        = 123.98d0
      e_top        = 1.24d3
      appx_mth     = 'Rate/2 + Mdot/2'
      a_tau        = 0.0d0
      allocate(melem_ab(n_melem))
      melem_ab     = 0.0d0

      call one_configuration('wasp_full',  nfail)
      call one_configuration('mol_metals', nfail)

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'photon_grid_threshold_edges: ', nfail,         &
                             ' assertion(s) FAILED'
         call exit(1)
      endif
      write(*,'(A)') 'photon_grid_threshold_edges: all assertions PASSED'

      contains

      subroutine one_configuration(label, nfail)
      character(len=*), intent(in)    :: label
      integer,          intent(inout) :: nfail

      ! The reconstructed bin edges of the grid this configuration builds.
      real*8, allocatable :: e_edge(:)
      real*8  :: e_floor, worst, band
      integer :: i, j
      character(len=8) :: name_worst

      ! Free what a previous configuration allocated; set_energy_vectors
      ! allocates these unconditionally.
      if (allocated(e_v))       deallocate(e_v)
      if (allocated(de_v))      deallocate(de_v)
      if (allocated(F_XUV))     deallocate(F_XUV)
      if (allocated(sigma_tab)) deallocate(sigma_tab)
      if (allocated(s_h2_di))   deallocate(s_h2_di)
      if (allocated(s_h2_dd))   deallocate(s_h2_dd)
      if (allocated(s_h2_nd))   deallocate(s_h2_nd)

      ! metals.inp of both cases: solar C, N, O, Mg, Ca, Na, Fe
      ! (Asplund et al. 2009 number ratios to H).
      melem_ab            = 0.0d0
      melem_ab(iel_C)     = 2.69d-4
      melem_ab(iel_N)     = 6.76d-5
      melem_ab(iel_O)     = 4.90d-4
      melem_ab(iel_Mg)    = 3.98d-5
      melem_ab(iel_Ca)    = 2.19d-6
      melem_ab(iel_Na)    = 1.74d-6
      melem_ab(iel_Fe)    = 3.16d-5
      thereis_metals      = .true.
      thereis_lowIP_metal = .true.
      thereis_HeITR       = .true.

      select case (label)
         case ('wasp_full')
            LX            = 29.46d0
            LEUV          = 30.42d0
            a_orb         = 0.02544d0*AU
            use_excited_H = .true.
            thereis_mol   = .false.
         case ('mol_metals')
            LX            = 27.20d0
            LEUV          = 27.93d0
            a_orb         = 0.0480d0*AU
            use_excited_H = .false.
            thereis_mol   = .true.
      end select

      call set_energy_vectors

      e_floor = photon_grid_floor_eV()
      allocate(e_edge(Nl+1))
      e_edge(1) = e_floor
      do j = 1,Nl
         e_edge(j+1) = e_edge(j) + de_v(j)
      enddo

      write(*,'(A,A,A,I0,A,I0,A,F9.4,A,F9.4,A)') '  configuration ', label,  &
           ': Nl=', Nl, '  NlTR=', NlTR, '  grid [', e_floor, ',',           &
           e_edge(Nl+1), '] eV'

      ! --- (a) every active threshold is a bin edge ----------------------
      worst      = 0.0d0
      name_worst = '        '
      call threshold_distance('H I  ', e_th_HI,   e_edge, worst, name_worst)
      call threshold_distance('He I ', e_th_HeI,  e_edge, worst, name_worst)
      call threshold_distance('He II', e_th_HeII, e_edge, worst, name_worst)
      if (thereis_HeITR)                                                     &
         call threshold_distance('He2^3S', e_th_HeTR, e_edge, worst,         &
                                 name_worst)
      if (thereis_mol)                                                       &
         call threshold_distance('H2   ', e_th_H2, e_edge, worst, name_worst)
      do i = 1,n_mion
         if (.not. mion_isphot(i))                    cycle
         if (mion_ethr(i) .le. 0.0d0)                 cycle
         if (melem_ab(mion_elem(i)) .le. 0.0d0)       cycle
         call threshold_distance(mion_name(i), mion_ethr(i), e_edge, worst,  &
                                 name_worst)
      enddo
      call verdict('threshold_on_bin_edge['//label//']', worst, 0.0d0,       &
                   1.0d-12, 'farthest '//trim(name_worst), nfail)

      ! --- the two grid invariants the partition carries -----------------
      band = e_top - e_floor
      call verdict('grid_covers_the_band['//label//']',                      &
                   abs(sum(de_v) - band)/band, 0.0d0, 1.0d-12, ' ', nfail)
      if (Nl .eq. Nl_fix + NlTR) then
         write(*,'(A,A,A,I0,A,I0,A)') 'PASS bin_count_unchanged[', label,    &
              '] measured=', Nl, ' reference=', Nl_fix + NlTR, ' tol=0'
      else
         write(*,'(A,A,A,I0,A,I0,A)') 'FAIL bin_count_unchanged[', label,    &
              '] measured=', Nl, ' reference=', Nl_fix + NlTR, ' tol=0'
         nfail = nfail + 1
      endif

      ! --- (b) the metal photoionization rates ---------------------------
      ! Only in the atomic configuration: the quadrature statement is about
      ! the grid, and the two configurations differ in it only through the
      ! H2 edge, which assertion (a) already covers.
      if (label .eq. 'wasp_full') then
         call rate_verdict('CI   ', 1,  1.0d-3, nfail)
         call rate_verdict('CII  ', 2,  1.0d-3, nfail)
         call rate_verdict('OI   ', 3,  1.0d-3, nfail)
         call rate_verdict('MgII ', 8,  1.0d-3, nfail)
         call rate_verdict('MgI  ', 7,  1.0d-3, nfail)
         call rate_verdict('CaI  ', 11, 1.0d-3, nfail)
         call rate_verdict('NaI  ', 13, 1.0d-3, nfail)
         call rate_verdict('FeI  ', 16, 1.0d-3, nfail)
      endif

      deallocate(e_edge)

      end subroutine one_configuration

      subroutine threshold_distance(name, e_thr, e_edge, worst, name_worst)
      ! Relative distance from one ionization threshold to the NEAREST bin
      ! edge, printed as context and accumulated into the configuration's
      ! worst case.  A threshold outside the grid band is skipped: it is not
      ! a partition point of a grid that does not reach it.
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: e_thr
      real*8,           intent(in)    :: e_edge(:)
      real*8,           intent(inout) :: worst
      character(len=8), intent(inout) :: name_worst
      integer :: j
      real*8  :: dist, d

      if (e_thr .lt. e_edge(1) .or. e_thr .gt. e_edge(Nl+1)) then
         write(*,'(A,A,A,F10.4,A)') '    threshold ', name, ' = ', e_thr,    &
              ' eV is outside the grid band, not a partition point'
         return
      endif
      dist = huge(1.0d0)
      do j = 1,Nl+1
         d = abs(e_edge(j) - e_thr)/e_thr
         if (d .lt. dist) dist = d
      enddo
      write(*,'(A,A,A,F10.4,A,ES10.3)') '    threshold ', name, ' = ',       &
           e_thr, ' eV, relative distance to the nearest bin edge ', dist
      if (dist .gt. worst) then
         worst      = dist
         name_worst = name
      endif

      end subroutine threshold_distance

      subroutine rate_verdict(name, k, tol, nfail)
      ! The code's P_m against the fine reference of the same integrand.
      character(len=*), intent(in)    :: name
      integer,          intent(in)    :: k
      real*8,           intent(in)    :: tol
      integer,          intent(inout) :: nfail
      integer :: i
      real*8  :: p_code, p_ref, p_ref2, e_thr, ratio

      p_code = sum(F_XUV*sigma_tab(:,k)/e_v*de_v)*1.0d-18*erg2eV
      e_thr  = 0.0d0
      do i = 1,n_mion
         if (mion_iphot(i) .eq. k) e_thr = mion_ethr(i)
      enddo
      call fine_metal_rate(k, e_thr, 400001, p_ref)
      call fine_metal_rate(k, e_thr, 200001, p_ref2)
      ratio = p_code/p_ref
      write(*,'(A,ES10.3)') '    reference self-convergence (400001 vs '//   &
           '200001, relative): ', abs(p_ref-p_ref2)/p_ref
      if (abs(ratio - 1.0d0) .le. tol) then
         write(*,'(A,A,A,ES22.15,A,ES22.15,A,ES9.2,A,F10.6)')                &
              'PASS P_m[', trim(name), '] measured=', p_code,                &
              ' reference=', p_ref, ' tol=', tol, ' code/exact=', ratio
      else
         write(*,'(A,A,A,ES22.15,A,ES22.15,A,ES9.2,A,F10.6)')                &
              'FAIL P_m[', trim(name), '] measured=', p_code,                &
              ' reference=', p_ref, ' tol=', tol, ' code/exact=', ratio
         nfail = nfail + 1
      endif

      end subroutine rate_verdict

      subroutine fine_metal_rate(k, e_thr, npt, p_m)
      ! Trapezoid rule on npt logarithmic points over [e_thr, e_top], the
      ! support of F sigma_m/E inside the run's photon band, from the same
      ! production functions the grid was filled with.
      integer, intent(in)  :: k, npt
      real*8,  intent(in)  :: e_thr
      real*8,  intent(out) :: p_m
      integer :: i
      real*8  :: e, w, dlog, xi_day, e_bot

      e_bot  = max(e_thr, photon_grid_floor_eV())
      dlog   = log(e_top/e_bot)/dble(npt-1)
      xi_day = dayside_dilution()
      p_m    = 0.0d0
      do i = 1, npt
         e = e_bot*exp(dlog*dble(i-1))
         ! trapezoid weight in E of a logarithmic abscissa
         w = e*dlog
         if (i .eq. 1 .or. i .eq. npt) w = 0.5d0*w
         p_m = p_m + xi_day*J_inc(e)*metal_photoion_sigma(k,e)/e*w
      enddo
      p_m = p_m*1.0d-18*erg2eV

      end subroutine fine_metal_rate

      subroutine verdict(name, measured, reference, tol, note, nfail)
      character(len=*), intent(in)    :: name, note
      real*8,           intent(in)    :: measured, reference, tol
      integer,          intent(inout) :: nfail
      if (abs(measured - reference) .le. tol) then
         write(*,'(A,A,A,ES22.15,A,ES22.15,A,ES9.2,1X,A)') 'PASS ', name,    &
              ' measured=', measured, ' reference=', reference,              &
              ' tol=', tol, trim(note)
      else
         write(*,'(A,A,A,ES22.15,A,ES22.15,A,ES9.2,1X,A)') 'FAIL ', name,    &
              ' measured=', measured, ' reference=', reference,              &
              ' tol=', tol, trim(note)
         nfail = nfail + 1
      endif
      end subroutine verdict

      end program photon_grid_threshold_edges
