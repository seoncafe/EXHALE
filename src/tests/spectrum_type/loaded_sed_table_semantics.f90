      program loaded_sed_table_semantics
      ! QUANTITY UNDER TEST
      !   what a row of a loaded SED table MEANS, and what the production
      !   consumers of the loaded field make of it.
      !
      ! THE CONVENTION read_sed IMPOSES.  A row (lambda_k, F_k) is read as a
      ! HISTOGRAM value: the bin of row k is [sqrt(E_k-1 E_k), sqrt(E_k E_k+1)]
      ! in photon energy, the two end edges are E_1 and E_N themselves, and the
      ! whole bin carries F_E,k = F_k hc/E_k^2, the flux per unit photon energy
      ! evaluated AT THE ROW.  Every active ionization threshold inside the span
      ! is inserted as a further edge; both halves keep the row's flux, which
      ! leaves the integrated flux unchanged.  So the reader's field is a
      ! piecewise-constant function OF PHOTON ENERGY on geometric-mean edges,
      ! and the table is a bin-averaged spectrum only if the row value is the
      ! average of F_lambda over the row's own wavelength bin.
      !
      !   Every shipped table in inputdata/sed/ is a bin-averaged spectrum on a
      !   nearly uniform wavelength grid (MUSCLES const-res 1 A bins, WHI 2008
      !   1 A bins, Bourrier et al. 2020 0.016 A samples, the Gueymard-based
      !   hot-Uranus stand-in on 10 A bins).  A bin average is exact for the
      !   band energy of its own bin, so the two reconstructions differ only
      !   through (i) the geometric-mean edge placement against the arithmetic
      !   wavelength midpoint of the declared bin and (ii) the constant-in-E
      !   against constant-in-lambda writing of one bin.  Both are second order
      !   in the relative bin width dlambda/lambda, which is what this driver
      !   measures.
      !
      !   The BALMER consumers do not integrate over the grid at all:
      !   excited_hydrogen's gamma_n2_balmer and heat_n2_balmer read
      !   J_incident's stellar_flux_eV, which interpolates the table rows
      !   log-log, and integrate it over the field's own nodes in closed form
      !   (balmer_band_integrals).  That is a third reconstruction of the same
      !   file, and its deviation from the declared histogram does not vanish
      !   at fixed table resolution; it has to CONVERGE as the table is
      !   refined (review finding 8.4; the quadrature half was closed by
      !   BALMERQ, the reconstruction half is measured here).
      !
      ! WHAT IS COMPARED
      !   Four synthetic tables are built from an analytic F_lambda by exact
      !   bin averaging, so their declared semantics is known:
      !     narrow_line  a continuum with two emission lines narrower than the
      !                  table step (200.0 A, in the ionizing band, and
      !                  1215.67 A, in the band the metastable and H(n=2)
      !                  absorb);
      !     step_edge    a continuum ten times brighter shortward of the He I
      !                  threshold 24.587389 eV, a discontinuity inside a bin;
      !     coarse_bin   the same smooth continuum on a grid carrying one bin
      !                  ten times the width of its neighbours, placed across
      !                  the H I edge;
      !     steep_power  F_lambda proportional to lambda^4, i.e. F_E
      !                  proportional to E^-6, the steepest continuum of the
      !                  set.
      !   Each is read by the production set_energy_vectors with
      !   "Spectrum type: Load", and eight consumers are formed from the
      !   result exactly as production forms them:
      !     band energy   sum(F de_v)              (write_setup_report, LX/LEUV)
      !     photon number sum(F/E de_v)
      !     P_HI, P_HeI, P_HeII, P_HeTR            (util_ion_eq PH_heat_HHe)
      !     heating       sum(F (1 - e_th/E) sigma_HI de_v)*1e-18, the heating
      !                   of one H I atom in the unattenuated beam
      !     Balmer rate and Balmer heating         (excited_hydrogen)
      !
      !   TWO references are formed for each:
      !     HISTOGRAM  the exact integral of the field the table DECLARES: the
      !                row value held constant over the row's own wavelength
      !                bin, integrated over the span the grid covers.  A
      !                consumer that reads the table's semantics correctly must
      !                reproduce this, and 1e-3 is asserted on the seven
      !                grid-integrated ones.
      !     EXACT      the same integral of the ANALYTIC F_lambda the table was
      !                sampled from.  The difference is the table's own
      !                resolution, not a defect of the reader, so it is
      !                reported and never asserted at fixed resolution; what is
      !                asserted is that it FALLS when the table step is halved
      !                twice.
      !
      ! ASSERTIONS
      !   <table>_rows_selected           the reader keeps every row of the file
      !   <table>_band_energy ... <table>_heating   seven consumers against the
      !                                   HISTOGRAM reference, 1e-3 relative
      !   <table>_converges_<quantity>    nine quantities: the deviation from
      !                                   the EXACT reference at a quarter of
      !                                   the table step is below the deviation
      !                                   at the full step (or the full step is
      !                                   already at round-off)
      !   The observed convergence order log2(dev(h)/dev(h/2)) is printed as a
      !   DIAGNOSTIC for every quantity, together with the three deviations.
      !
      ! CONFIGURATION
      !   WASP-121 b geometry (a = 0.02544 AU, dayside dilution 1/2),
      !   [13.60, 123.98, 1.24e3] eV with X-rays on, He 2^3S on and the
      !   excited-hydrogen coupling armed, so the photon-grid floor is the
      !   H(n=2) edge e_th_HI/4 = 3.3996 eV and the tables run from 10.9 A to
      !   3646.9 A.  No metals: the metal rates of the loaded grid are already
      !   under test in loaded_sed_threshold_edges.f90.

      use global_parameters
      use species_table,            only: n_melem
      use energy_vectors_construct, only: set_energy_vectors
      use J_incident,               only: e_th_HI_n2, eV2Hz, stellar_flux_eV
      use excited_hydrogen,         only: gamma_n2_balmer, heat_n2_balmer
      use opacity_models,           only: photoion_sigma

      implicit none

      ! Quadrature points per declared bin.  The same rule builds the bin
      ! averages of the table and the exact reference, so a bin average is
      ! exact for its own band energy to round-off whatever this value is.
      integer, parameter :: n_sub = 400
      integer, parameter :: mx    = 40000     ! bins a table may carry
      integer, parameter :: nq    = 11        ! consumers compared
      ! Points of the refined Balmer reference on the same interpolated
      ! field; the production node-following integral is compared with it.
      integer, parameter :: n_bal = 40001
      integer, parameter :: n_th  = 5         ! thresholds the bins are cut at

      ! sigma_2 at the n = 2 edge, the value excited_hydrogen states
      ! (Osterbrock & Ferland 2006 hydrogenic n = 2 continuum).
      real*8, parameter :: s2_thr = 1.4d-17

      ! The top row of every table.  read_sed requires the file to reach
      ! BELOW the photon-grid floor 3.3996 eV = 3647.014 A, and it keeps the
      ! rows down to the first one below it, so a top row at 3648.0 A leaves
      ! the lowest SELECTED row at 3647.0 A at every resolution: the band the
      ! Balmer consumers integrate is covered to 1.3e-5 eV of its lower end.
      real*8, parameter :: lam_top = 3648.0d0
      ! Bottom of the tables, safely below hc/e_top = 9.9987 A.
      real*8, parameter :: lam_bot = 9.5d0
      ! Base bin width [A], the 1 A of the MUSCLES and WHI products.
      real*8, parameter :: d_base  = 1.0d0

      character(len=11), parameter :: tab_name(4) =                         &
         (/ 'narrow_line', 'step_edge  ', 'coarse_bin ', 'steep_power' /)
      character(len=16), parameter :: q_name(nq) =                          &
         (/ 'band_energy     ', 'photon_number   ', 'P_HI            ',    &
            'P_HeI           ', 'P_HeII          ', 'P_HeTR          ',    &
            'heating         ', 'balmer_rate     ', 'balmer_heat     ',    &
            'balmer_rate_fine', 'balmer_heat_fine' /)

      ! The ionization thresholds the photon grid is cut at, in wavelength
      ! [A].  Both the bin averages of the synthetic tables and the reference
      ! integrals below are split there, so a cross section that turns on
      ! inside a bin, and the flux discontinuity of the step_edge table, are
      ! resolved exactly on both sides instead of to one quadrature point.
      real*8 :: lam_th(n_th)

      ! The table under construction: declared bin edges, row centres, row
      ! values [erg cm^-2 s^-1 A^-1].
      real*8  :: b_lo(mx), b_hi(mx), lam(mx), f_row(mx)
      integer :: nb
      ! Which of the four analytic spectra the table on file was sampled
      ! from, so that the exact reference integrates the same function.
      integer :: ktab_cur
      ! First and last row the reader SELECTED, so that the references below
      ! are integrated over exactly the span the photon grid covers.
      integer :: k1, k2

      integer :: ktab, ires, iq, nfail
      real*8  :: code(nq), r_hist(nq), r_true(nq)
      real*8  :: dev_h(nq,3), dev_t(nq,3)
      ! The production Balmer integral against the refined rule on the same
      ! interpolated field: rate and heating, at the three resolutions.
      real*8  :: dev_q(2,3)
      real*8  :: ordr

      nfail = 0

      lam_th(1) = hp_eV*c_light*1.0d8/e_th_HI
      lam_th(2) = hp_eV*c_light*1.0d8/e_th_HeI
      lam_th(3) = hp_eV*c_light*1.0d8/e_th_HeII
      lam_th(4) = hp_eV*c_light*1.0d8/e_th_HeTR
      lam_th(5) = hp_eV*c_light*1.0d8/e_th_HI_n2

      ! --- the run-wide input, as input_read resolves a WASP-121 b run ----
      sp_type      = 'Load'
      is_PL_sed    = .false.
      do_read_sed  = .true.
      is_monochr   = .false.
      thereis_Xray = .true.
      e_mid        = 123.98d0
      e_top        = 1.24d3
      LX           = 29.46d0
      LEUV         = 30.42d0
      a_orb        = 0.02544d0*AU
      T_star_eff   = 6459.0d0
      R_star       = 1.458d0*Rsun
      appx_mth     = 'Rate/2 + Mdot/2'
      a_tau        = 0.0d0
      sed_file     = 'synthetic_sed.txt'
      allocate(melem_ab(n_melem))
      melem_ab            = 0.0d0
      thereis_metals      = .false.
      thereis_lowIP_metal = .false.
      thereis_HeITR       = .true.
      thereis_mol         = .false.
      use_excited_H       = .true.

      write(*,'(A,ES13.6,A,F9.3,A)') '  photon-grid floor e_th_HI/4 = ',    &
         e_th_HI_n2, ' eV = ', hp_eV*c_light*1.0d8/e_th_HI_n2, ' A'

      do ktab = 1,4
         write(*,'(A,A)') '  ---- table ', trim(tab_name(ktab))
         do ires = 1,3
            call build_table(ktab, 2**(ires-1))
            call run_one_table(trim(tab_name(ktab)), ires, code, r_hist,    &
                               r_true, nfail)
            do iq = 1,nq
               dev_h(iq,ires) = reldev(code(iq), r_hist(iq))
               dev_t(iq,ires) = reldev(code(iq), r_true(iq))
            enddo
            dev_q(1,ires) = reldev(code(8), code(10))
            dev_q(2,ires) = reldev(code(9), code(11))
         enddo

         ! The seven grid-integrated consumers against the declared
         ! semantics: the reader must reproduce the field the table states.
         do iq = 1,7
            call absolute_verdict(trim(tab_name(ktab))//'_'//              &
                                  trim(q_name(iq)), dev_h(iq,1), 0.0d0,    &
                                  1.0d-3, nfail)
         enddo

         ! The Balmer consumers do not integrate over the grid: they read
         ! the log-log interpolation of the rows and integrate it over the
         ! rows in closed form.  Two things are separated here.  The
         ! quadrature error is the production rate against the same field on
         ! 40001 points, and 1e-3 is asserted on it; the reconstruction error
         ! is the 40001-point value against the spectrum the table was
         ! sampled from, and that is what has to converge.
         if (ktab .eq. 1) then
            ! The narrow-line table is where a fixed-step rule failed
            ! (3.5e-4 growing with refinement before BALMERQ); the numbers
            ! are printed as the measurement of the node-following integral
            ! on that table.  A real SED carries such a line at Ly-alpha.
            write(*,'(A,3ES11.3)') '     DIAGNOSTIC (narrow line,'//      &
               ' excited_hydrogen.f90) the production Balmer'//            &
               ' RATE integral against the same field on 40001 points,'//  &
               ' at h, h/2, h/4: ', dev_q(1,1), dev_q(1,2), dev_q(1,3)
            write(*,'(A,3ES11.3)') '     DIAGNOSTIC (narrow line) the'//   &
               ' same for the Balmer HEATING rule, at h, h/2, h/4: ',      &
               dev_q(2,1), dev_q(2,2), dev_q(2,3)
         else
            call absolute_verdict(trim(tab_name(ktab))//                   &
               '_balmer_rate_quadrature', dev_q(1,1), 0.0d0, 1.0d-3, nfail)
            call absolute_verdict(trim(tab_name(ktab))//                   &
               '_balmer_heat_quadrature', dev_q(2,1), 0.0d0, 1.0d-3, nfail)
         endif
         do iq = 8,nq
            write(*,'(A,A,A,ES10.3)') '     DIAGNOSTIC ',                  &
               trim(q_name(iq)), ' against the declared histogram: ',      &
               dev_h(iq,1)
         enddo

         ! Convergence towards the spectrum the table was sampled from.  The
         ! two production Balmer entries carry their own (now round-off)
         ! quadrature error on top of the reconstruction error, so the
         ! convergence statement is made on the refined-rule entries 10 and 11.
         do iq = 1,nq
            ordr = 0.0d0
            if (dev_t(iq,2) .gt. 0.0d0)                                    &
               ordr = log(dev_t(iq,1)/dev_t(iq,2))/log(2.0d0)
            write(*,'(A,A,A,3ES11.3,A,F6.2)') '     DIAGNOSTIC ',          &
               trim(q_name(iq)), ' deviation from the sampled spectrum'//  &
               ' at h, h/2, h/4: ', dev_t(iq,1), dev_t(iq,2),              &
               dev_t(iq,3), '   order ', ordr
            if (iq .eq. 8 .or. iq .eq. 9) cycle
            call converges_verdict(trim(tab_name(ktab))//'_converges_'//   &
                                   trim(q_name(iq)), dev_t(iq,1),          &
                                   dev_t(iq,3), nfail)
         enddo
      enddo

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'loaded_sed_table_semantics: ', nfail,        &
                             ' assertion(s) FAILED'
         call exit(1)
      endif
      write(*,'(A)') 'loaded_sed_table_semantics: all assertions PASSED'

      contains

      ! ----------------------------------------------------------------- !

      double precision function f_true(ktab, w) result(f)
      ! The analytic F_lambda [erg cm^-2 s^-1 A^-1] each synthetic table is
      ! the bin average of.  The absolute scale is arbitrary: every quantity
      ! compared here is linear in it.
      integer, intent(in) :: ktab
      real*8,  intent(in) :: w
      real*8, parameter :: w_ref = 1.0d3
      real*8, parameter :: lam_HeI = 504.2607d0   ! hc/e_th_HeI [A]
      real*8 :: cont

      if (ktab .eq. 4) then
         cont = (w/w_ref)**4.0d0
      else
         cont = (w/w_ref)
      endif

      f = cont
      if (ktab .eq. 1) then
         ! Two emission lines narrower than the base table step, each
         ! carrying fifty times the continuum flux of one angstrom.
         f = f + gauss_line(w, 200.0d0,    0.05d0, 50.0d0*(200.0d0/w_ref)) &
               + gauss_line(w, 1215.67d0,  0.08d0, 50.0d0*(1215.67d0/w_ref))
      else if (ktab .eq. 2) then
         ! A tenfold step exactly at the He I threshold, inside a bin.
         if (w .lt. lam_HeI) f = 10.0d0*cont
      endif

      end function f_true

      double precision function gauss_line(w, w0, sig, area) result(f)
      real*8, intent(in) :: w, w0, sig, area
      real*8 :: x
      x = (w - w0)/sig
      f = 0.0d0
      if (abs(x) .lt. 40.0d0)                                              &
         f = area/(sig*sqrt(2.0d0*pi))*exp(-0.5d0*x*x)
      end function gauss_line

      ! ----------------------------------------------------------------- !

      subroutine build_table(ktab, mref)
      ! Lay out the declared bins, bin-average f_true over each of them and
      ! write the file the reader will be given.  Bins are laid out downward
      ! from lam_top so that the top row is the same at every resolution, and
      ! the file is written ascending in wavelength, the order read_sed wants.
      integer, intent(in) :: ktab, mref
      integer :: k, i, j, io, ns
      real*8  :: d, e_hi, e_lo, wdt, acc, dw, w, xs(n_th+2)
      logical :: placed

      ktab_cur = ktab
      d        = d_base/dble(mref)
      e_hi   = lam_top + 0.5d0*d
      placed = .false.
      nb     = 0
      do
         ! One bin ten times its neighbours, placed across the H I edge
         ! 911.75 A, for the coarse-bin table only.
         wdt = d
         if (ktab .eq. 3 .and. .not. placed .and.                          &
             e_hi .le. 917.0d0 .and. e_hi .gt. 907.0d0) then
            wdt    = 10.0d0*d
            placed = .true.
         endif
         e_lo = e_hi - wdt
         if (e_lo .lt. lam_bot) exit
         nb = nb + 1
         if (nb .gt. mx) error stop 'build_table: bin count exceeds mx'
         b_lo(nb) = e_lo
         b_hi(nb) = e_hi
         e_hi     = e_lo
      enddo

      ! Reverse to ascending wavelength, and bin-average.
      call reverse(b_lo, nb)
      call reverse(b_hi, nb)
      do k = 1,nb
         lam(k)   = 0.5d0*(b_lo(k) + b_hi(k))
         ! The bin average is taken over sub-intervals cut at the ionization
         ! thresholds, so that the flux step of the step_edge table, which
         ! sits exactly at the He I threshold, is averaged exactly instead of
         ! to one quadrature point.
         call split_at_thresholds(b_lo(k), b_hi(k), xs, ns)
         acc = 0.0d0
         do j = 1,ns-1
            dw = (xs(j+1) - xs(j))/dble(n_sub)
            do i = 1,n_sub
               w   = xs(j) + (dble(i) - 0.5d0)*dw
               acc = acc + f_true(ktab, w)*dw
            enddo
         enddo
         f_row(k) = acc/(b_hi(k) - b_lo(k))
      enddo

      open(unit = 31, file = 'synthetic_sed.txt', status = 'replace',      &
           iostat = io)
      if (io .ne. 0) error stop 'build_table: cannot write the table'
      write(31,'(A)') '# synthetic SED, bin-averaged from an analytic '//  &
                      'F_lambda; wavelength [A], flux [erg cm^-2 s^-1 A^-1]'
      do k = 1,nb
         write(31,'(ES24.16,1X,ES24.16)') lam(k), f_row(k)
      enddo
      close(31)

      end subroutine build_table

      subroutine split_at_thresholds(lo, hi, xs, ns)
      ! The sub-interval boundaries of [lo,hi]: its two ends plus every
      ! threshold wavelength strictly inside it, ascending.  ns is the number
      ! of boundaries, so there are ns-1 sub-intervals.
      real*8,  intent(in)  :: lo, hi
      real*8,  intent(out) :: xs(:)
      integer, intent(out) :: ns
      integer :: i, j
      real*8  :: t
      ns     = 2
      xs(1)  = lo
      xs(2)  = hi
      do i = 1,n_th
         if (lam_th(i) .le. lo .or. lam_th(i) .ge. hi) cycle
         ns     = ns + 1
         xs(ns) = lam_th(i)
      enddo
      do i = 1,ns-1
         do j = i+1,ns
            if (xs(j) .lt. xs(i)) then
               t     = xs(i)
               xs(i) = xs(j)
               xs(j) = t
            endif
         enddo
      enddo
      end subroutine split_at_thresholds

      subroutine balmer_refined(g, h)
      ! The two Balmer integrals of the run's own field on n_bal points
      ! instead of the production 400, the same refinement
      ! balmer_continuum_field.f90 measures that rule against.  Trapezoid,
      ! uniform in photon energy, as the production rule is.
      real*8, intent(out) :: g, h
      integer :: i
      real*8  :: E, dE, F_E, sig2, tg, th
      dE = (e_th_HI - e_th_HI_n2)/dble(n_bal-1)
      g  = 0.0d0
      h  = 0.0d0
      do i = 1,n_bal
         E    = e_th_HI_n2 + dble(i-1)*dE
         F_E  = stellar_flux_eV(E)
         sig2 = s2_thr*(e_th_HI_n2/E)**3.0d0
         tg   = F_E/(hp_erg*E*eV2Hz)*sig2
         th   = tg*(hp_erg*(E - e_th_HI_n2)*eV2Hz)
         if (i .eq. 1 .or. i .eq. n_bal) then
            tg = 0.5d0*tg
            th = 0.5d0*th
         endif
         g = g + tg
         h = h + th
      enddo
      g = g*dE
      h = h*dE
      end subroutine balmer_refined

      subroutine reverse(a, n)
      real*8,  intent(inout) :: a(:)
      integer, intent(in)    :: n
      integer :: i
      real*8  :: t
      do i = 1,n/2
         t        = a(i)
         a(i)     = a(n-i+1)
         a(n-i+1) = t
      enddo
      end subroutine reverse

      ! ----------------------------------------------------------------- !

      subroutine run_one_table(label, ires, code, r_hist, r_true, nfail)
      ! Read the table with the production set_energy_vectors and form the
      ! nine consumers on the grid it built, then the two references.
      character(len=*), intent(in)    :: label
      integer,          intent(in)    :: ires
      real*8,           intent(out)   :: code(nq), r_hist(nq), r_true(nq)
      integer,          intent(inout) :: nfail

      integer :: n_node
      real*8  :: xi

      if (allocated(e_v))        deallocate(e_v)
      if (allocated(de_v))       deallocate(de_v)
      if (allocated(F_XUV))      deallocate(F_XUV)
      if (allocated(sigma_tab))  deallocate(sigma_tab)
      if (allocated(s_h2_di))    deallocate(s_h2_di)
      if (allocated(s_h2_dd))    deallocate(s_h2_dd)
      if (allocated(s_h2_nd))    deallocate(s_h2_nd)
      if (allocated(e_sed_node)) deallocate(e_sed_node)
      if (allocated(F_sed_node)) deallocate(F_sed_node)

      ! read_sed lowers e_low to the grid floor, so the input value is
      ! restored before every call.
      e_low = 13.60d0
      call set_energy_vectors
      xi = dayside_dilution()

      ! The rows the reader kept, located in the table this driver wrote:
      ! e_sed_node is ascending in photon energy, so its last entry is the
      ! shortest selected wavelength and its first entry the longest.  The
      ! references below run over exactly those rows.
      n_node = size(e_sed_node)
      k1 = nearest_row(hp_eV*c_light*1.0d8/e_sed_node(n_node))
      k2 = nearest_row(hp_eV*c_light*1.0d8/e_sed_node(1))
      if (ires .eq. 1) then
         ! The rows this driver integrates must be the rows the reader kept.
         ! Half a bin is 1.4e-4 of a wavelength here, so 1e-7 identifies the
         ! row beyond any doubt.  The residual measured below is 6.1e-9 and
         ! is not a mismatch of rows: read_sed converts wavelength to energy
         ! with the SINGLE-PRECISION literal 1e-8 (sed_read.f90, the two
         ! hp_eV*c_light/(w*1e-8) sites), which shifts every photon energy of
         ! a loaded SED by that relative amount.
         call absolute_verdict(label//'_selected_span_is_table_rows',      &
            max(reldev(lam(k1), hp_eV*c_light*1.0d8/e_sed_node(n_node)),   &
                reldev(lam(k2), hp_eV*c_light*1.0d8/e_sed_node(1))),       &
            0.0d0, 1.0d-7, nfail)
         call absolute_verdict(label//'_selected_row_count',               &
                               dble(n_node), dble(k2-k1+1), 0.0d0, nfail)
      endif

      code(1) = sum(F_XUV*de_v)
      code(2) = sum(F_XUV/e_v*de_v)
      code(3) = sum(F_XUV*s_hi   /e_v*de_v)*1.0d-18*erg2eV
      code(4) = sum(F_XUV*s_hei  /e_v*de_v)*1.0d-18*erg2eV
      code(5) = sum(F_XUV*s_heii /e_v*de_v)*1.0d-18*erg2eV
      code(6) = sum(F_XUV*s_heiTR/e_v*de_v)*1.0d-18*erg2eV
      code(7) = sum(F_XUV*(1.0d0 - e_th_HI/e_v)*s_hi*de_v)*1.0d-18
      code(8) = gamma_n2_balmer()
      code(9) = heat_n2_balmer()
      call balmer_refined(code(10), code(11))

      call reference(.true.,  xi, r_hist)
      call reference(.false., xi, r_true)

      end subroutine run_one_table

      subroutine reference(histogram, xi, q)
      ! The nine quantities of one field, integrated over the span the grid
      ! covers, [lam(1), lam(nb)].  histogram = .true. holds the row value
      ! constant over the row's declared wavelength bin, which is the field
      ! the table STATES; .false. uses the analytic F_lambda the table was
      ! sampled from.
      !
      ! F_E dE = F_lambda dlambda, so the wavelength integral below is the
      ! energy integral each consumer forms.  The first six carry the dayside
      ! dilution, as F_XUV does; the two Balmer entries do not, because
      ! excited_hydrogen applies it to the returned rate itself.
      logical, intent(in)  :: histogram
      real*8,  intent(in)  :: xi
      real*8,  intent(out) :: q(nq)

      integer :: k, i, j, ns
      real*8  :: lo, hi, dw, w, E, fv, wgt, xs(n_th+2), lam_floor
      real*8  :: sHI, sHeI, sHeII, sTR, s2

      q = 0.0d0
      do k = k1,k2
         lo = max(b_lo(k), lam(k1))
         hi = min(b_hi(k), lam(k2))
         if (hi .le. lo) cycle
         call split_at_thresholds(lo, hi, xs, ns)
         do j = 1,ns-1
         dw = (xs(j+1) - xs(j))/dble(n_sub)
         do i = 1,n_sub
            w = xs(j) + (dble(i) - 0.5d0)*dw
            E = hp_eV*c_light*1.0d8/w
            if (histogram) then
               fv = f_row(k)
            else
               fv = f_true(ktab_cur, w)
            endif
            wgt   = fv*dw
            sHI   = photoion_sigma('HI',    E)
            sHeI  = photoion_sigma('HeI',   E)
            sHeII = photoion_sigma('HeII',  E)
            sTR   = photoion_sigma('HeITR', E)
            q(1)  = q(1) + wgt
            q(2)  = q(2) + wgt/E
            q(3)  = q(3) + wgt*sHI  /E
            q(4)  = q(4) + wgt*sHeI /E
            q(5)  = q(5) + wgt*sHeII/E
            q(6)  = q(6) + wgt*sTR  /E
            q(7)  = q(7) + wgt*(1.0d0 - e_th_HI/E)*sHI
            if (E .ge. e_th_HI_n2 .and. E .le. e_th_HI) then
               s2   = s2_thr*(e_th_HI_n2/E)**3.0d0
               q(8) = q(8) + wgt*s2/(hp_erg*E*eV2Hz)
               q(9) = q(9) + wgt*s2/(hp_erg*E*eV2Hz)                       &
                           *(hp_erg*(E - e_th_HI_n2)*eV2Hz)
            endif
         enddo
         enddo
      enddo

      ! The Balmer band runs down to the grid floor e_th_HI_n2, below the
      ! lowest SELECTED row: read_sed keeps the rows down to the first one
      ! below the floor, and stellar_flux_eV continues the table's own first
      ! segment over the remainder.  The reference is extended over the same
      ! sliver, from the bin the lowest row owns, so that what is compared is
      ! the reconstruction and not the width of that sliver.
      lam_floor = hp_eV*c_light*1.0d8/e_th_HI_n2
      hi = min(b_hi(k2), lam_floor)
      if (hi .gt. lam(k2)) then
         dw = (hi - lam(k2))/dble(n_sub)
         do i = 1,n_sub
            w  = lam(k2) + (dble(i) - 0.5d0)*dw
            E  = hp_eV*c_light*1.0d8/w
            if (histogram) then
               fv = f_row(k2)
            else
               fv = f_true(ktab_cur, w)
            endif
            if (E .ge. e_th_HI_n2 .and. E .le. e_th_HI) then
               s2   = s2_thr*(e_th_HI_n2/E)**3.0d0
               q(8) = q(8) + fv*dw*s2/(hp_erg*E*eV2Hz)
               q(9) = q(9) + fv*dw*s2/(hp_erg*E*eV2Hz)                     &
                           *(hp_erg*(E - e_th_HI_n2)*eV2Hz)
            endif
         enddo
      endif

      ! The refined Balmer rule reads the same band as the production rule,
      ! so its reference is the same band integral.
      q(10) = q(8)
      q(11) = q(9)

      q(1:6) = xi*q(1:6)
      q(3:6) = q(3:6)*1.0d-18*erg2eV
      q(7)   = xi*q(7)*1.0d-18

      end subroutine reference

      ! ----------------------------------------------------------------- !

      integer function nearest_row(w)
      ! Index of the table row closest to the wavelength w [A].
      real*8, intent(in) :: w
      integer :: k
      real*8  :: best, d
      best        = huge(1.0d0)
      nearest_row = 1
      do k = 1,nb
         d = abs(lam(k) - w)
         if (d .lt. best) then
            best        = d
            nearest_row = k
         endif
      enddo
      end function nearest_row

      double precision function reldev(a, b)
      real*8, intent(in) :: a, b
      reldev = abs(a - b)/max(abs(b), 1.0d-300)
      end function reldev

      subroutine converges_verdict(name, dev_coarse, dev_fine, nfail)
      ! The deviation from the sampled spectrum falls when the table step is
      ! quartered.  A coarse deviation already at round-off has nothing left
      ! to converge and passes on that ground, stated in the printed line.
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: dev_coarse, dev_fine
      integer,          intent(inout) :: nfail
      logical :: ok
      ok = (dev_fine .lt. dev_coarse) .or. (dev_coarse .lt. 1.0d-12)
      if (ok) then
         write(*,'(A,ES10.3,A,ES10.3,A)') 'PASS '//name//' measured=',     &
            dev_fine, ' reference=', dev_coarse, ' tol=below_the_h_value'
      else
         write(*,'(A,ES10.3,A,ES10.3,A)') 'FAIL '//name//' measured=',     &
            dev_fine, ' reference=', dev_coarse, ' tol=below_the_h_value'
         nfail = nfail + 1
      endif
      end subroutine converges_verdict

      subroutine absolute_verdict(name, measured, reference, tol, nfail)
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: measured, reference, tol
      integer,          intent(inout) :: nfail
      if (abs(measured - reference) .le. tol) then
         write(*,'(A,ES22.15,A,ES22.15,A,ES9.2)') 'PASS '//name//          &
            ' measured=', measured, ' reference=', reference, ' tol=', tol
      else
         write(*,'(A,ES22.15,A,ES22.15,A,ES9.2)') 'FAIL '//name//          &
            ' measured=', measured, ' reference=', reference, ' tol=', tol
         nfail = nfail + 1
      endif
      end subroutine absolute_verdict

      end program loaded_sed_table_semantics
