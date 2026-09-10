      program wasp121b_sed_sensitivity
      ! QUANTITY UNDER TEST
      !   the production photoionization, Balmer-continuum and FUV band
      !   quantities of the THREE WASP-121 b spectrum candidates, formed on
      !   ONE prescribed unattenuated field at the planet's orbit:
      !     A  inputdata/sed/wasp121b_solar_huang2023.txt
      !        solar (WHI 2008) XUV below 1700 A, 6459 K blackbody above
      !     B  inputdata/sed/wasp121b_wasp17_huang2023.txt
      !        WASP-17 (F6V, MUSCLES v24) throughout
      !     C  inputdata/sed/wasp121b_composite_huang2023.txt
      !        A's solar XUV below 1700 A, B's WASP-17 photosphere above
      !
      ! WHAT THE COMPARISON IS, AND WHAT IT IS NOT
      !   The three differ only in the input spectrum: the same photon grid
      !   is built by the same set_energy_vectors, from the same input keys,
      !   with no atmosphere in front of it (no column, no optical depth) and
      !   only the run's own dayside dilution applied.  The differences
      !   between the three rates are therefore a SENSITIVITY of the
      !   production rates to the choice of candidate input.  They are NOT a
      !   bias of any candidate: a bias would need a reference field for
      !   WASP-121, and the workspace has none -- the synthetic spectrum of
      !   Huang et al. (2023) and an LLmodels run for T_eff = 6459 K are both
      !   unavailable (inputdata/sed/README.md, WASP-121 b section).
      !
      ! WHY IT MATTERS
      !   The He 2^3S metastable photoionizes out of the 1700-2600.5 A band,
      !   which is the band in which candidate A carries a blackbody and
      !   candidates B and C carry a measured F-star photosphere; the H, He
      !   and He+ ground-state rates are set below 912 A, where candidates A
      !   and C carry the same solar shape and candidate B does not.  The
      !   composite exists so that each band is taken from the candidate
      !   whose construction is defensible there, and the two assertions
      !   below are the statement that it is exactly that and nothing else.
      !
      ! ASSERTIONS (references and tolerances)
      !   1 P_HI, P_HeI and P_HeII of C against A, 1e-6 relative.  Their
      !     thresholds are 13.598, 24.587 and 54.418 eV, so the integrands
      !     live below 912 A, where C carries A's rows unchanged.
      !   2 the metastable rate, in three statements.  P_HeTR is NOT a
      !     photosphere-only quantity: its threshold is 4.767775 eV, so the
      !     integrand runs from 2600.5 A all the way down to the 10 A top of
      !     the grid and takes part of the rate from the XUV.  So
      !     2a P_HeTR of C against the sum of B's part below the 1700 A join
      !        and A's part above it, 1e-3 relative: the composite rate is
      !        exactly its two components and nothing else;
      !     2b the part of P_HeTR formed below the join (bins with
      !        e_v < hc/1700 A) of C against B's, 1e-3 relative;
      !     2c the part formed above the join of C against A's, 1e-2
      !        relative.
      !     Neither part is exact, because the one bin that straddles 1700 A
      !     has different edges in each file: the geometric-mean edge is
      !     built from the last row below the join and the first row above
      !     it, and those pairs are (1699.5, 1700.5) in A, (1699.0, 1700.0)
      !     in B and (1699.5, 1700.0) in C.  That single bin is the whole of
      !     the departure, and the tolerances are its size, not a fit.
      !     The full P_HeTR of C is 8 per cent below B's, and the table
      !     prints it: B's XUV is WASP-17's and C's is solar.
      !   3 F(1-912 A) of C, integrated from the file by the production row
      !     reader, against the 1.6e6 erg cm^-2 s^-1 its header states,
      !     1e-3 relative.  This is Huang et al. (2023) section 2.2 with the
      !     band definition of Salz et al. (2019), A&A 623, A57.
      !   4 F(10-912 A) of C against A, 1e-12 relative: the whole XUV
      !     component is A's.
      !   5 F(1700-2600.5 A) of C against B, 1e-12 relative: the whole
      !     metastable band is B's.
      !   6 the driver's band integral against the production
      !     lyman_werner_band_flux_from_sed(), 1e-12 relative on all three
      !     candidates, so the band fluxes it prints are the production
      !     prescription and not a second one.  The LW band is 912-1201 A
      !     on both sides since 2026-09-06.
      !
      !   Everything else is printed as DIAGNOSTIC: the three rates side by
      !   side with their C/A and C/B ratios, the Balmer-continuum rate and
      !   heating of excited_hydrogen, and the five FUV band fluxes.
      !
      ! CONFIGURATION: WASP-121b/input.inp and its metals.inp --
      !   a 0.02544 AU, [13.60, 123.98, 1.24e3] eV, X-rays on,
      !   "2D approximate method: Rate/2 + Mdot/2", "Include He23S? True",
      !   Stellar Teff 6459 K and Stellar radius 1.458 R_sun (both positive,
      !   so the excited-hydrogen coupling is armed and the photon-grid floor
      !   is the H(n=2) edge 3.400 eV), "Stellar Lya flux 1.0e5",
      !   solar C, N, O, Mg, Ca, Na and Fe -- with "Spectrum type: Load" in
      !   place of the power law the input carries today, which is the switch
      !   this comparison is the prerequisite for.
      !
      ! EXPECTED BEFORE THE C2 CHANGE: RED -- the composite file does not
      ! exist, so set_energy_vectors stops on the missing table.

      use global_parameters
      use species_table,            only: n_melem, iel_C, iel_N, iel_O,      &
                                          iel_Mg, iel_Ca, iel_Na, iel_Fe
      use energy_vectors_construct, only: set_energy_vectors
      use sed_reader,               only: sed_next_row,                      &
                                          lyman_werner_band_flux_from_sed
      use excited_hydrogen,         only: gamma_n2_balmer, heat_n2_balmer

      implicit none

      integer, parameter :: n_cand = 3
      ! The FUV band edges [A] of water_photolysis, in its band order
      ! (LW, B2 = Ly-alpha, B3, B4); parameters.f90 states B3 and B4.  The
      ! Ly-alpha row is integrated over 1206-1226 A rather than the band's
      ! nominal 1202-1230 A, so that the line and not its continuum wings is
      ! what the row reports.
      integer, parameter :: n_band = 4
      real*8,  parameter :: w_band_lo(n_band) =                              &
         (/  912.0d0, 1206.0d0, 1231.0d0, 1451.0d0 /)
      real*8,  parameter :: w_band_hi(n_band) =                              &
         (/ 1201.0d0, 1226.0d0, 1450.0d0, 2304.0d0 /)
      character(len=2), parameter :: band_name(n_band) =                     &
         (/ 'LW', 'B2', 'B3', 'B4' /)
      ! The join of the two components, and the long edge of the metastable
      ! band (hc/e_th_HeTR).
      real*8,  parameter :: w_join = 1700.0d0

      character(len=1),  parameter :: cand_tag(n_cand)  = (/ 'A','B','C' /)
      character(len=40)            :: cand_file(n_cand)

      real*8  :: p(4,n_cand)          ! P_HI, P_HeI, P_HeII, P_HeTR [s^-1]
      real*8  :: p_tr_lo(n_cand)      ! the part of P_HeTR below the join
      real*8  :: p_tr_hi(n_cand)      ! the part of P_HeTR above the join
      real*8  :: g_bal(n_cand), h_bal(n_cand)
      real*8  :: f_band(n_band,n_cand), f_lw_prod(n_cand)
      real*8  :: f_xuv1(n_cand), f_xuv10(n_cand), f_tr_band(n_cand)
      real*8  :: w_tr

      integer :: nfail, i, k
      character(len=512) :: root
      character(len=6), parameter :: rate_label(4) =                             &
         (/ 'P_HI  ', 'P_HeI ', 'P_HeII', 'P_HeTR' /)

      nfail = 0
      w_tr  = hp_eV*c_light*1.0d8/e_th_HeTR

      cand_file(1) = 'wasp121b_solar_huang2023.txt'
      cand_file(2) = 'wasp121b_wasp17_huang2023.txt'
      cand_file(3) = 'wasp121b_composite_huang2023.txt'

      call get_environment_variable('EXHALE_TEST_ROOT', root)
      if (len_trim(root) .eq. 0) root = '../../../..'

      ! --- the run-wide input, as input_read resolves WASP-121b -----------
      sp_type      = 'Load'
      is_PL_sed    = .false.
      do_read_sed  = .true.
      is_monochr   = .false.
      thereis_Xray = .true.
      e_low        = 13.60d0
      e_mid        = 123.98d0
      e_top        = 1.24d3
      PLind        = -1.0d0
      LX           = 29.46d0
      LEUV         = 30.42d0
      a_orb        = 0.02544d0*AU
      T_star_eff   = 6459.0d0
      R_star       = 1.458d0*Rsun
      F_Lya_star   = 1.0d5
      appx_mth     = 'Rate/2 + Mdot/2'
      a_tau        = 0.0d0

      ! metals.inp of WASP-121b: Asplund et al. (2009) number ratios to H.
      allocate(melem_ab(n_melem))
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
      thereis_mol         = .false.
      use_excited_H       = .true.

      write(*,'(A,F8.5,A,ES10.3)') '  one prescribed unattenuated field:'// &
         ' a = ', a_orb/AU, ' AU, dayside dilution = ', dayside_dilution()
      write(*,'(A,F9.4,A,F9.4,A)') '  join at ', w_join,                    &
         ' A; metastable band edge hc/e_th_HeTR = ', w_tr, ' A'

      do i = 1,n_cand
         call one_candidate(i)
      enddo

      call sensitivity_table

      ! ---- 1  the XUV rates of C are A's --------------------------------
      do k = 1,3
         call relative_verdict('C_equals_A_'//trim(rate_label(k)),              &
                               p(k,3), p(k,1), 1.0d-6, nfail)
      enddo

      ! ---- 2  the metastable rate of C is B's below the join and A's
      !         above it, and the sum of the two is all of it --------------
      call relative_verdict('C_P_HeTR_is_B_below_join_plus_A_above',        &
                            p(4,3), p_tr_lo(2) + p_tr_hi(1), 1.0d-3, nfail)
      call relative_verdict('C_equals_B_P_HeTR_below_the_join',             &
                            p_tr_lo(3), p_tr_lo(2), 1.0d-3, nfail)
      call relative_verdict('C_equals_A_P_HeTR_above_the_join',             &
                            p_tr_hi(3), p_tr_hi(1), 1.0d-2, nfail)

      ! ---- 3  the composite integrates to its stated F_XUV --------------
      call relative_verdict('C_F_1_912_is_the_stated_F_XUV', f_xuv1(3),     &
                            1.6d6, 1.0d-3, nfail)

      ! ---- 4, 5  the two components are carried unchanged ---------------
      call relative_verdict('C_F_10_912_equals_A', f_xuv10(3), f_xuv10(1),  &
                            1.0d-12, nfail)
      call relative_verdict('C_F_1700_2600_equals_B', f_tr_band(3),         &
                            f_tr_band(2), 1.0d-12, nfail)

      ! ---- 6  the driver's band integral is the production one ----------
      do i = 1,n_cand
         call relative_verdict('band_integral_is_production_LW_'//          &
                               cand_tag(i), f_band(1,i), f_lw_prod(i),      &
                               1.0d-12, nfail)
      enddo

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'wasp121b_sed_sensitivity: ', nfail,           &
                             ' assertion(s) FAILED'
         call exit(1)
      endif
      write(*,'(A)') 'wasp121b_sed_sensitivity: all assertions PASSED'

      contains

      ! ----------------------------------------------------------------- !

      subroutine one_candidate(i)
      ! Build the production photon grid from candidate i and measure every
      ! quantity of the comparison on it.
      integer, intent(in) :: i
      integer :: j, ib
      real*8  :: e_join, xi

      if (allocated(e_v))        deallocate(e_v)
      if (allocated(de_v))       deallocate(de_v)
      if (allocated(F_XUV))      deallocate(F_XUV)
      if (allocated(sigma_tab))  deallocate(sigma_tab)
      if (allocated(s_h2_di))    deallocate(s_h2_di)
      if (allocated(s_h2_dd))    deallocate(s_h2_dd)
      if (allocated(s_h2_nd))    deallocate(s_h2_nd)
      if (allocated(e_sed_node)) deallocate(e_sed_node)
      if (allocated(F_sed_node)) deallocate(F_sed_node)

      sed_file = trim(root)//'/inputdata/sed/'//trim(cand_file(i))
      ! read_sed lowers e_low to the photon-grid floor, so the input value
      ! is restored before every call.
      e_low    = 13.60d0
      call set_energy_vectors

      xi = dayside_dilution()
      write(*,'(A,A,A,A,A,I0,A,I0)') '  candidate ', cand_tag(i), '  ',     &
         trim(cand_file(i)), ': table rows kept=', size(e_sed_node),        &
         '  grid bins Nl=', Nl

      ! The four rates, formed exactly as util_ion_eq forms them:
      ! sum(F sigma/e_v de_v)*1e-18*erg2eV, the cross sections in Mb and
      ! F_XUV already diluted by set_energy_vectors.
      p(1,i) = sum(F_XUV*s_hi   /e_v*de_v)*1.0d-18*erg2eV
      p(2,i) = sum(F_XUV*s_hei  /e_v*de_v)*1.0d-18*erg2eV
      p(3,i) = sum(F_XUV*s_heii /e_v*de_v)*1.0d-18*erg2eV
      p(4,i) = sum(F_XUV*s_heiTR/e_v*de_v)*1.0d-18*erg2eV

      ! The same metastable rate split at the join, so that the part formed
      ! from the photosphere is separated from the part formed from the XUV.
      e_join     = hp_eV*c_light*1.0d8/w_join
      p_tr_lo(i) = 0.0d0
      p_tr_hi(i) = 0.0d0
      do j = 1,Nl
         if (e_v(j) .lt. e_join) then
            p_tr_lo(i) = p_tr_lo(i) + F_XUV(j)*s_heiTR(j)/e_v(j)*de_v(j)
         else
            p_tr_hi(i) = p_tr_hi(i) + F_XUV(j)*s_heiTR(j)/e_v(j)*de_v(j)
         endif
      enddo
      p_tr_lo(i) = p_tr_lo(i)*1.0d-18*erg2eV
      p_tr_hi(i) = p_tr_hi(i)*1.0d-18*erg2eV

      ! The Balmer continuum of excited_hydrogen, on this run's own field.
      ! Both are undiluted at the source (excited_H_update applies the
      ! dilution), so it is applied here to keep the table one convention.
      g_bal(i) = xi*gamma_n2_balmer()
      h_bal(i) = xi*heat_n2_balmer()

      ! Band fluxes of the file itself, by the prescription
      ! lyman_werner_band_flux_from_sed uses: trapezoid on the rows, clipped
      ! to the band, at the planet, undiluted (the file is already at the
      ! planet's orbit).  Assertion 6 checks the transcription against the
      ! production function on the band that function owns.
      do ib = 1,n_band
         f_band(ib,i) = file_band_flux(w_band_lo(ib), w_band_hi(ib))
      enddo
      f_lw_prod(i) = lyman_werner_band_flux_from_sed()
      f_xuv1(i)    = file_band_flux(1.0d0,  911.6d0)
      f_xuv10(i)   = file_band_flux(10.0d0, 911.6d0)
      f_tr_band(i) = file_band_flux(w_join, w_tr)

      end subroutine one_candidate

      ! ----------------------------------------------------------------- !

      double precision function file_band_flux(w_lo, w_hi) result(F_band)
      ! Band-integrated flux of sed_file over [w_lo, w_hi] at the planet
      ! [erg cm^-2 s^-1].  A transcription of lyman_werner_band_flux_from_sed
      ! with the band an argument: the production row reader, a trapezoid on
      ! consecutive rows clipped to the band, and no dilution, because an
      ! EXHALE SED file is already the flux at the planet's orbit.
      real*8, intent(in) :: w_lo, w_hi
      real*8  :: w, f, w_prev, f_prev, wa, wb
      integer :: io, nin
      logical :: have_prev
      F_band    = 0.0d0
      nin       = 0
      have_prev = .false.
      w_prev    = 0.0d0
      f_prev    = 0.0d0
      open(unit = 73, file = sed_file, status = 'old', iostat = io)
      if (io .ne. 0) return
      do
         if (.not. sed_next_row(73, w, f, io)) exit
         if (io .ne. 0) exit
         if (have_prev .and. w .gt. w_prev) then
            wa = max(w_prev, w_lo)
            wb = min(w,      w_hi)
            if (wb .gt. wa) then
               F_band = F_band + 0.5d0*(f_prev + f)*(wb - wa)
               nin    = nin + 1
            endif
         endif
         w_prev    = w
         f_prev    = f
         have_prev = .true.
         if (w .gt. w_hi) exit
      enddo
      close(73)
      if (nin .lt. 2) F_band = 0.0d0
      end function file_band_flux

      ! ----------------------------------------------------------------- !

      subroutine sensitivity_table
      ! The three candidates side by side.  These are DIFFERENCES BETWEEN
      ! CANDIDATE INPUTS.  Calling any of them a bias would need a reference
      ! field for WASP-121, which the workspace does not have.
      integer :: k, ib
      write(*,'(A)') ''
      write(*,'(A)') '  DIAGNOSTIC SENSITIVITY of the production rates to'//&
                     ' the candidate input.  These are'
      write(*,'(A)') '  DIAGNOSTIC differences between candidate inputs,'// &
                     ' not a bias of any candidate:'
      write(*,'(A)') '  DIAGNOSTIC a bias needs a reference field for'//    &
                     ' WASP-121, which the workspace has none of.'
      write(*,'(A)') '  DIAGNOSTIC quantity            A'//                 &
                     '             B             C          C/A      C/B'
      do k = 1,4
         call sens_row(trim(rate_label(k))//' [s^-1]', p(k,1), p(k,2), p(k,3))
      enddo
      call sens_row('P_HeTR below join', p_tr_lo(1), p_tr_lo(2), p_tr_lo(3))
      call sens_row('P_HeTR above join', p_tr_hi(1), p_tr_hi(2), p_tr_hi(3))
      call sens_row('Balmer rate [s^-1]', g_bal(1), g_bal(2), g_bal(3))
      call sens_row('Balmer heat [erg/s]', h_bal(1), h_bal(2), h_bal(3))
      call sens_row('F(1-912 A)', f_xuv1(1), f_xuv1(2), f_xuv1(3))
      call sens_row('F(10-912 A)', f_xuv10(1), f_xuv10(2), f_xuv10(3))
      call sens_row('F(1700-2600.5 A)', f_tr_band(1), f_tr_band(2),         &
                    f_tr_band(3))
      do ib = 1,n_band
         call sens_row('FUV '//band_name(ib)//' band', f_band(ib,1),        &
                       f_band(ib,2), f_band(ib,3))
      enddo
      write(*,'(A)') '  DIAGNOSTIC the FUV band fluxes above are the'//     &
                     ' spectrum integrated over each band.  What'
      write(*,'(A)') '  DIAGNOSTIC the equilibrium solve uses is'//         &
                     ' fuv_band_flux(ib), which returns the input'
      write(*,'(A)') '  DIAGNOSTIC keys (Stellar LW / Lya / B3'//          &
                     ' / B4 flux), and only the LW band is'
      write(*,'(A)') '  DIAGNOSTIC taken from the spectrum file, and only'//&
                     ' in a molecular run.  B3 and B4 are'
      write(*,'(A)') '  DIAGNOSTIC therefore identical for the three'//     &
                     ' candidates as the input stands.'
      write(*,'(A)') ''
      end subroutine sensitivity_table

      subroutine sens_row(name, va, vb, vc)
      character(len=*), intent(in) :: name
      real*8,           intent(in) :: va, vb, vc
      character(len=20) :: nm
      nm = name
      write(*,'(A,A,3ES14.6,2F9.4)') '  DIAGNOSTIC ', nm, va, vb, vc,       &
         ratio(vc, va), ratio(vc, vb)
      end subroutine sens_row

      double precision function ratio(x, y)
      real*8, intent(in) :: x, y
      if (y .eq. 0.0d0) then
         ratio = 0.0d0
      else
         ratio = x/y
      endif
      end function ratio

      ! ----------------------------------------------------------------- !

      subroutine relative_verdict(name, measured, reference, tol, nfail)
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: measured, reference, tol
      integer,          intent(inout) :: nfail
      real*8 :: dev
      dev = reldev(measured, reference)
      if (dev .le. tol) then
         write(*,'(A,A,A,ES16.9,A,ES16.9,A,ES9.2,A,ES9.2)')                 &
            'PASS ', name, ' measured=', measured, ' reference=',           &
            reference, ' tol=', tol, ' deviation=', dev
      else
         write(*,'(A,A,A,ES16.9,A,ES16.9,A,ES9.2,A,ES9.2)')                 &
            'FAIL ', name, ' measured=', measured, ' reference=',           &
            reference, ' tol=', tol, ' deviation=', dev
         nfail = nfail + 1
      endif
      end subroutine relative_verdict

      double precision function reldev(a, b)
      real*8, intent(in) :: a, b
      if (b .eq. 0.0d0) then
         reldev = abs(a)
      else
         reldev = abs(a - b)/abs(b)
      endif
      end function reldev

      end program wasp121b_sed_sensitivity
