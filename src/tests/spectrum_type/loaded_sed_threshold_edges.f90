      program loaded_sed_threshold_edges
      ! QUANTITY UNDER TEST
      !   the photon grid read_sed builds from a LOADED SED, and the
      !   photoionization rates the code forms on it.
      !
      ! WHY IT MATTERS
      !   A photoionization rate is a rectangle rule sum(F sigma/E de_v) over
      !   the bins, and a cross section turns ON at its threshold.  A bin
      !   that CONTAINS a threshold charges the absorber's cross section,
      !   read at the bin point, over the part of the bin lying BELOW the
      !   threshold, where the absorber cannot absorb at all.  The analytic
      !   spectrum types cut their grid at every threshold for that reason
      !   (set_energy_vectors; development plan rev 3 section 10.1 item 2).
      !   Every planet calculation of record runs on a literature SED
      !   (decision 16), so the loaded path is the one that matters most.
      !
      ! THE GRID UNDER TEST
      !   rows E_1 < ... < E_N of the file; interior edges sqrt(E_k E_k+1),
      !   end edges E_1 and E_N, each row's flux over its whole bin, and
      !   every active ionization threshold strictly inside the span
      !   inserted as a further edge (splitting one bin into two halves that
      !   both carry that row's flux).
      !
      !   QUADRATURE IDENTITY.  Splitting a bin [b,b'] at t and giving both
      !   halves the flux F leaves F(t-b) + F(b'-t) = F(b'-b): the split is
      !   exact for the integrated flux.  Hence sum(de_v) = E_N - E_1 and
      !   sum(F_XUV de_v) = sum_k F_k (b_k - b_k-1) whatever is inserted.
      !
      ! ASSERTIONS (references and tolerances)
      !   1 every active ionization threshold inside the table span, against
      !     the nearest bin edge, 1e-12 relative.  The edges are
      !     RECONSTRUCTED from the grid floor and the cumulative de_v, as
      !     photon_grid_threshold_edges does.
      !   2 sum(de_v) against E_N - E_1, 1e-12 relative.
      !   3 sum(F_XUV de_v) against the histogram integral of the table over
      !     its own span, times the dayside dilution, 1e-12 relative.
      !   4 P_HI, P_HeI, P_HeII, P_HeTR and every active metal rate, formed
      !     exactly as util_ion_eq forms them, sum(F sigma/e_v de_v)*1e-18
      !     *erg2eV, against the reference integral of the SAME table over
      !     [threshold, E_N]: piecewise-constant flux per table bin, the
      !     cross section on 1000 midpoints of each clipped bin.  1e-3
      !     relative.
      !
      !   The PREVIOUS grid -- the table rows as bin points with
      !   central-difference widths, no threshold an edge -- is reproduced in
      !   this driver and its rates are printed against the same reference as
      !   a DIAGNOSTIC, so the straddling error each rate carried is on the
      !   record next to the value that replaces it.
      !
      ! CONFIGURATIONS
      !   wasp52        benchmarks/wasp52/input.inp: the eps Eri table
      !                 inputdata/sed/wasp52b_epseri_yan2022.txt, a 0.0272 AU,
      !                 [13.60, 123.98, 1.24e3] eV, X-rays on, He 2^3S on,
      !                 stellar Teff and radius set (so the excited-hydrogen
      !                 coupling is armed and the grid floor is the H(n=2)
      !                 edge 3.400 eV), no metals.inp.
      !   wasp52_metals the same, with solar C, N, O, Mg, Ca, Na and Fe, so
      !                 that the metal thresholds are exercised as well.
      !
      ! EXPECTED BEFORE THE 2c-SEDGRID CHANGE: RED -- the driver does not
      ! build (read_sed took no argument and the table rows were not kept),
      ! and the grid it tests did not exist.

      use global_parameters
      use species_table,            only: n_melem, n_mion, mion_ethr,        &
                                          mion_iphot, mion_elem, mion_name,  &
                                          iel_C, iel_N, iel_O, iel_Mg,       &
                                          iel_Ca, iel_Na, iel_Fe
      use energy_vectors_construct, only: set_energy_vectors,                &
                                          ionization_thresholds_active,     &
                                          n_thr_max
      use Cross_sections,           only: metal_photoion_sigma
      use opacity_models,           only: photoion_sigma

      implicit none

      integer, parameter :: n_sub = 1000     ! midpoints per table bin
      integer :: nfail
      character(len=512) :: root

      nfail = 0

      ! The repository root, seen from the directory the driver runs in
      ! (build/tests/spectrum_type/<name>_run), so that the SED file of
      ! inputdata/ is found without copying it.
      call get_environment_variable('EXHALE_TEST_ROOT', root)
      if (len_trim(root) .eq. 0) root = '../../../..'

      ! --- the run-wide input, as input_read resolves benchmarks/wasp52 ---
      sp_type      = 'Load'
      is_PL_sed    = .false.
      do_read_sed  = .true.
      is_monochr   = .false.
      thereis_Xray = .true.
      e_low        = 13.60d0
      e_mid        = 123.98d0
      e_top        = 1.24d3
      LX           = 28.0d0
      LEUV         = 29.0d0
      a_orb        = 0.0272d0*AU
      appx_mth     = 'Rate/2 + Mdot/2'
      a_tau        = 0.0d0
      sed_file     = trim(root)//                                           &
                     '/inputdata/sed/wasp52b_epseri_yan2022.txt'
      allocate(melem_ab(n_melem))
      melem_ab     = 0.0d0

      call one_configuration('wasp52',        nfail)
      call one_configuration('wasp52_metals', nfail)

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'loaded_sed_threshold_edges: ', nfail,          &
                             ' assertion(s) FAILED'
         call exit(1)
      endif
      write(*,'(A)') 'loaded_sed_threshold_edges: all assertions PASSED'

      contains

      subroutine one_configuration(label, nfail)
      character(len=*), intent(in)    :: label
      integer,          intent(inout) :: nfail

      integer :: i, j, k, n_node, n_thr
      ! Sized by the module's own bound: every active threshold, the H2
      ! channel edges and, since METALS-INNER, every metal subshell turn-on.
      real*8  :: e_thr(n_thr_max)
      real*8, allocatable :: e_edge(:), b_lo(:), b_hi(:)
      real*8  :: worst, e_worst, dev_i, xi, span, f_int_ref, f_int_code
      real*8  :: p_code(4), p_ref(4), p_old(4)
      character(len=6), parameter :: p_name(4) =                            &
         (/ 'P_HI  ', 'P_HeI ', 'P_HeII', 'P_HeTR' /)

      ! Free what a previous configuration allocated; set_energy_vectors and
      ! read_sed allocate these unconditionally.
      if (allocated(e_v))        deallocate(e_v)
      if (allocated(de_v))       deallocate(de_v)
      if (allocated(F_XUV))      deallocate(F_XUV)
      if (allocated(sigma_tab))  deallocate(sigma_tab)
      if (allocated(s_h2_di))    deallocate(s_h2_di)
      if (allocated(s_h2_dd))    deallocate(s_h2_dd)
      if (allocated(s_h2_nd))    deallocate(s_h2_nd)
      if (allocated(e_sed_node)) deallocate(e_sed_node)
      if (allocated(F_sed_node)) deallocate(F_sed_node)

      melem_ab            = 0.0d0
      thereis_metals      = .false.
      thereis_lowIP_metal = .false.
      thereis_HeITR       = .true.
      thereis_mol         = .false.
      use_excited_H       = .true.
      e_low               = 13.60d0
      if (label .eq. 'wasp52_metals') then
         ! Asplund et al. (2009) number ratios to hydrogen.
         melem_ab(iel_C)     = 2.69d-4
         melem_ab(iel_N)     = 6.76d-5
         melem_ab(iel_O)     = 4.90d-4
         melem_ab(iel_Mg)    = 3.98d-5
         melem_ab(iel_Ca)    = 2.19d-6
         melem_ab(iel_Na)    = 1.74d-6
         melem_ab(iel_Fe)    = 3.16d-5
         thereis_metals      = .true.
         thereis_lowIP_metal = .true.
      endif

      call set_energy_vectors

      n_node = size(e_sed_node)
      xi     = dayside_dilution()
      span   = e_sed_node(n_node) - e_sed_node(1)
      write(*,'(A,A,A,I0,A,I0,A,F8.4,A,F9.3,A)') '  ', label,               &
         ': table rows=', n_node, '  grid bins Nl=', Nl, '  span [',        &
         e_sed_node(1), ',', e_sed_node(n_node), '] eV'

      ! The table's own bins, from which every reference below is built.
      allocate(b_lo(n_node), b_hi(n_node))
      do k = 1,n_node
         if (k .eq. 1) then
            b_lo(k) = e_sed_node(1)
         else
            b_lo(k) = sqrt(e_sed_node(k-1)*e_sed_node(k))
         endif
         if (k .eq. n_node) then
            b_hi(k) = e_sed_node(n_node)
         else
            b_hi(k) = sqrt(e_sed_node(k)*e_sed_node(k+1))
         endif
      enddo

      ! ---- 1  every active threshold inside the span is a bin edge -------
      ! The edges are reconstructed: the bins are contiguous and de_v is the
      ! exact width, so edge(1) = E_1 and edge(j+1) = edge(j) + de_v(j).
      allocate(e_edge(Nl+1))
      e_edge(1) = e_sed_node(1)
      do j = 1,Nl
         e_edge(j+1) = e_edge(j) + de_v(j)
      enddo
      call ionization_thresholds_active(e_thr, n_thr)
      worst   = 0.0d0
      e_worst = 0.0d0
      do i = 1,n_thr
         if (e_thr(i) .le. e_sed_node(1)) cycle
         if (e_thr(i) .ge. e_sed_node(n_node)) cycle
         dev_i = nearest_edge_deviation(e_thr(i), e_edge)
         if (e_worst .le. 0.0d0 .or. dev_i .gt. worst) then
            worst   = dev_i
            e_worst = e_thr(i)
         endif
      enddo
      write(*,'(A,I0,A,F10.5,A,ES10.3)') '     thresholds inside the span: ',&
         n_inside(e_thr, n_thr), ', worst at ', e_worst,                    &
         ' eV, relative distance to its edge ', worst
      call absolute_verdict(label//'_thresholds_are_bin_edges', worst,      &
                            0.0d0, 1.0d-12, nfail)

      ! ---- 2  the grid spans the table exactly ---------------------------
      call relative_verdict(label//'_grid_spans_the_table', sum(de_v),      &
                            span, 1.0d-12, nfail)

      ! ---- 3  the split is exact for the integrated flux -----------------
      f_int_ref = 0.0d0
      do k = 1,n_node
         f_int_ref = f_int_ref + F_sed_node(k)*(b_hi(k) - b_lo(k))
      enddo
      f_int_ref  = xi*f_int_ref
      f_int_code = sum(F_XUV*de_v)
      call relative_verdict(label//'_integrated_flux_preserved',            &
                            f_int_code, f_int_ref, 1.0d-12, nfail)

      ! ---- 4  the rates against the reference integral -------------------
      p_code(1) = sum(F_XUV*s_hi   /e_v*de_v)*1.0d-18*erg2eV
      p_code(2) = sum(F_XUV*s_hei  /e_v*de_v)*1.0d-18*erg2eV
      p_code(3) = sum(F_XUV*s_heii /e_v*de_v)*1.0d-18*erg2eV
      p_code(4) = sum(F_XUV*s_heiTR/e_v*de_v)*1.0d-18*erg2eV
      p_ref(1)  = reference_rate('HI   ', 0, e_th_HI,   b_lo, b_hi, xi)
      p_ref(2)  = reference_rate('HeI  ', 0, e_th_HeI,  b_lo, b_hi, xi)
      p_ref(3)  = reference_rate('HeII ', 0, e_th_HeII, b_lo, b_hi, xi)
      p_ref(4)  = reference_rate('HeITR', 0, e_th_HeTR, b_lo, b_hi, xi)
      p_old(1) = previous_grid_rate('HI   ', 0, xi)
      p_old(2) = previous_grid_rate('HeI  ', 0, xi)
      p_old(3) = previous_grid_rate('HeII ', 0, xi)
      p_old(4) = previous_grid_rate('HeITR', 0, xi)
      do i = 1,4
         write(*,'(A,A,A,ES14.7,A,ES10.3,A,ES10.3)')                        &
            '     DIAGNOSTIC ', p_name(i), ' [s^-1] reference=', p_ref(i),  &
            '  previous grid deviation=', reldev(p_old(i), p_ref(i)),       &
            '  this grid deviation=', reldev(p_code(i), p_ref(i))
         call relative_verdict(label//'_'//trim(p_name(i)), p_code(i),      &
                               p_ref(i), 1.0d-3, nfail)
      enddo

      call metal_rates(label, b_lo, b_hi, xi, nfail)

      deallocate(e_edge, b_lo, b_hi)
      end subroutine one_configuration

      ! ----------------------------------------------------------------- !

      integer function n_inside(e_thr, n_thr)
      real*8,  intent(in) :: e_thr(:)
      integer, intent(in) :: n_thr
      integer :: i, n_node
      n_node   = size(e_sed_node)
      n_inside = 0
      do i = 1,n_thr
         if (e_thr(i) .gt. e_sed_node(1) .and.                              &
             e_thr(i) .lt. e_sed_node(n_node)) n_inside = n_inside + 1
      enddo
      end function n_inside

      double precision function nearest_edge_deviation(e, e_edge)
      ! Relative distance from e to the nearest reconstructed bin edge.
      real*8, intent(in) :: e, e_edge(:)
      integer :: j
      nearest_edge_deviation = huge(1.0d0)
      do j = 1,size(e_edge)
         nearest_edge_deviation = min(nearest_edge_deviation,               &
                                      abs(e_edge(j) - e)/e)
      enddo
      end function nearest_edge_deviation

      double precision function reference_rate(species, kmetal, e_th,       &
                                               b_lo, b_hi, xi) result(p)
      ! The exact photoionization rate of the SAME table [s^-1]: the flux is
      ! the table's histogram, one value per row, and the cross section is
      ! integrated over each bin clipped to [e_th, E_N] on n_sub midpoints.
      ! kmetal > 0 selects a column of the metal cross-section table;
      ! kmetal = 0 selects the H/He absorber named by `species`.
      character(len=*), intent(in) :: species
      integer,          intent(in) :: kmetal
      real*8,           intent(in) :: e_th, b_lo(:), b_hi(:), xi
      integer :: k, i, n_node
      real*8  :: lo, hi, dE, E, acc
      n_node = size(b_lo)
      p      = 0.0d0
      do k = 1,n_node
         lo = max(b_lo(k), e_th)
         hi = b_hi(k)
         if (hi .le. lo) cycle
         dE  = (hi - lo)/dble(n_sub)
         acc = 0.0d0
         do i = 1,n_sub
            E = lo + (dble(i) - 0.5d0)*dE
            if (kmetal .gt. 0) then
               acc = acc + metal_photoion_sigma(kmetal, E)/E
            else
               acc = acc + photoion_sigma(trim(species), E)/E
            endif
         enddo
         p = p + F_sed_node(k)*acc*dE
      enddo
      p = xi*p*1.0d-18*erg2eV
      end function reference_rate

      double precision function previous_grid_rate(species, kmetal, xi)   &
                               result(p)
      ! One photoionization rate [s^-1] on the grid the loaded path built
      ! BEFORE the threshold-edge layout: the table rows themselves as the
      ! bin points, with central-difference widths and half-widths at the
      ! two ends, and no threshold a bin edge.  A transcription of that
      ! arithmetic, kept here so the straddling error it carried is measured
      ! against the same reference as the grid that replaces it.
      character(len=*), intent(in) :: species
      integer,          intent(in) :: kmetal
      real*8,           intent(in) :: xi
      integer :: k, n_node
      real*8, allocatable :: e_old(:), d_old(:)
      real*8  :: sig
      n_node = size(e_sed_node)
      allocate(e_old(n_node), d_old(n_node))
      e_old = e_sed_node
      d_old(1)          = 0.5d0*(e_old(2) - e_old(1))
      d_old(2:n_node-1) = 0.5d0*(e_old(3:n_node) - e_old(1:n_node-2))
      d_old(n_node)     = 0.5d0*(e_old(n_node) - e_old(n_node-1))
      p = 0.0d0
      do k = 1,n_node
         if (kmetal .gt. 0) then
            sig = metal_photoion_sigma(kmetal, e_old(k))
         else
            sig = photoion_sigma(trim(species), e_old(k))
         endif
         p = p + F_sed_node(k)*sig/e_old(k)*d_old(k)
      enddo
      p = xi*p*1.0d-18*erg2eV
      deallocate(e_old, d_old)
      end function previous_grid_rate

      subroutine metal_rates(label, b_lo, b_hi, xi, nfail)
      ! Every photo-ionizable ion of an element the run carries, as
      ! util_ion_eq forms its rate off sigma_tab, against the same reference
      ! integral.
      character(len=*), intent(in)    :: label
      real*8,           intent(in)    :: b_lo(:), b_hi(:), xi
      integer,          intent(inout) :: nfail
      integer :: k, kp, n_active
      real*8  :: pc, pr, po, worst, worst_old
      character(len=8) :: name_worst, name_worst_old
      worst          = 0.0d0
      worst_old      = 0.0d0
      name_worst     = 'none'
      name_worst_old = 'none'
      n_active       = 0
      do k = 1,n_mion
         kp = mion_iphot(k)
         if (kp .le. 0) cycle
         if (melem_ab(mion_elem(k)) .le. 0.0d0) cycle
         n_active = n_active + 1
         pc = sum(F_XUV*sigma_tab(:,kp)/e_v*de_v)*1.0d-18*erg2eV
         pr = reference_rate('     ', kp, mion_ethr(k), b_lo, b_hi, xi)
         po = previous_grid_rate('     ', kp, xi)
         if (reldev(pc, pr) .gt. worst) then
            worst      = reldev(pc, pr)
            name_worst = mion_name(k)
         endif
         if (reldev(po, pr) .gt. worst_old) then
            worst_old      = reldev(po, pr)
            name_worst_old = mion_name(k)
         endif
      enddo
      if (n_active .eq. 0) then
         write(*,'(A)') '     no metal is active in this configuration'
         return
      endif
      write(*,'(A,I0,A,A,A,ES10.3)') '     DIAGNOSTIC metal ions checked: ',&
         n_active, ', worst on the previous grid is ',                      &
         trim(name_worst_old), ', deviation ', worst_old
      write(*,'(A,A)') '     worst on this grid is ', trim(name_worst)
      call absolute_verdict(label//'_metal_rates', worst, 0.0d0, 1.0d-3,    &
                            nfail)
      end subroutine metal_rates

      double precision function reldev(a, b)
      real*8, intent(in) :: a, b
      reldev = abs(a - b)/max(abs(b), 1.0d-300)
      end function reldev

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

      subroutine absolute_verdict(name, measured, reference, tol, nfail)
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: measured, reference, tol
      integer,          intent(inout) :: nfail
      if (abs(measured - reference) .le. tol) then
         write(*,'(A,ES22.15,A,ES22.15,A,ES9.2)') 'PASS '//name//           &
            ' measured=', measured, ' reference=', reference, ' tol=', tol
      else
         write(*,'(A,ES22.15,A,ES22.15,A,ES9.2)') 'FAIL '//name//           &
            ' measured=', measured, ' reference=', reference, ' tol=', tol
         nfail = nfail + 1
      endif
      end subroutine absolute_verdict

      end program loaded_sed_threshold_edges
