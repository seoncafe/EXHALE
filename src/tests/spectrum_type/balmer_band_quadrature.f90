      program balmer_band_quadrature
      ! QUANTITY UNDER TEST
      !   the QUADRATURE of the Balmer-continuum band integrals of
      !   excited_hydrogen, gamma_n2_balmer and heat_n2_balmer, against the
      !   exact integral of the very field they read, stellar_flux_eV.
      !
      !   This is not the reconstruction question.  What field a loaded table
      !   states, and how well the log-log interpolation of its rows
      !   reproduces the spectrum the table was sampled from, is measured in
      !   loaded_sed_table_semantics.f90.  Here the field is taken as given
      !   and only the rule that integrates it is under test.
      !
      ! WHY IT MATTERS
      !   A quadrature has to resolve the field it integrates.  A rule whose
      !   nodes are fixed independently of the spectrum cannot: a measured
      !   stellar table carries emission lines narrower than any fixed step,
      !   Ly-alpha among them at 10.20 eV, inside this band, and refining the
      !   table narrows the interpolated line further while the step stays
      !   put, so the error GROWS under refinement instead of falling.
      !
      ! THE REFERENCE
      !   Every deviation below is against a reference integral of the SAME
      !   field, formed by a rule that shares no expression with the code
      !   under test: 5-point Gauss-Legendre on n_gsub geometric
      !   sub-intervals of each interval between consecutive NODES of the
      !   field (the table's own rows inside the band, plus the two band
      !   edges; a geometric subdivision of the band for the analytic
      !   types).  The field is smooth inside one such interval whatever it
      !   does across the band, so this converges to the exact integral for
      !   any line width, and doubling n_gsub is asserted to leave it
      !   unchanged.
      !
      !   The 400001-point uniform-in-energy trapezoid the review named is
      !   formed as well and reported next to it, so that the accuracy of
      !   that reference is on the record too.
      !
      ! WHAT IS COMPARED, and against what
      !   production   gamma_n2_balmer / heat_n2_balmer
      !   uniform400   the FIXED 400-point uniform-in-energy trapezoid the
      !                production functions used before this driver was
      !                written, reproduced here so that the rule replaced is
      !                measured by the same driver on the same fields
      !   trap400k     the 400001-point uniform-in-energy trapezoid
      !   reference    the node-following Gauss-Legendre rule above
      !
      ! CONFIGURATIONS
      !   narrow_line  the synthetic bin-averaged table of
      !                loaded_sed_table_semantics.f90, continuum plus two
      !                emission lines narrower than the table step (200.0 A
      !                and 1215.67 A, the second one inside this band), built
      !                at the 1 A step of the shipped products and at a half
      !                and a quarter of it.  Same construction, so the
      !                numbers of the two drivers are comparable.
      !   wasp52       inputdata/sed/wasp52b_epseri_yan2022.txt, the eps Eri
      !                table of benchmarks/wasp52, which carries a measured
      !                Ly-alpha line inside the band.
      !   power_law    "Spectrum type: Power-law", index -1, the field of
      !                backup/regression/wasp_full.
      !   planck       "Spectrum type: Planck", the 5000 K photosphere of
      !                benchmarks/wasp52 at 0.0272 AU.
      !
      ! ASSERTIONS
      !   <config>_rate_quadrature, <config>_heat_quadrature
      !                production against the reference, 1e-6 relative
      !   <config>_reference_converged
      !                the reference against itself at twice n_gsub, 1e-12
      !   narrow_line_quadrature_converges
      !                the production deviation at a quarter of the table
      !                step is not above the deviation at the full step
      !   The deviation of uniform400 from the reference is printed as a
      !   DIAGNOSTIC for every configuration, at every resolution.
      !
      ! EXPECTED AGAINST THE FIXED 400-POINT RULE: RED.  The production
      ! functions were that rule, so <config>_rate_quadrature and
      ! <config>_heat_quadrature fail on every configuration and
      ! narrow_line_quadrature_converges fails as well.

      use global_parameters
      use species_table,            only: n_melem
      use energy_vectors_construct, only: set_energy_vectors
      use J_incident,               only: e_th_HI_n2, eV2Hz, stellar_flux_eV
      use excited_hydrogen,         only: gamma_n2_balmer, heat_n2_balmer

      implicit none

      ! Quadrature points per declared bin of the synthetic table, and the
      ! largest bin count it may carry.  Both are the values
      ! loaded_sed_table_semantics.f90 builds the same table with.
      integer, parameter :: n_sub = 400
      integer, parameter :: mx    = 40000
      integer, parameter :: n_th  = 5
      ! Points of the uniform trapezoid rules: the rule replaced, and the
      ! refined uniform rule the review named.
      integer, parameter :: n_u400 = 400
      integer, parameter :: n_trap = 400001
      ! Geometric sub-intervals per node interval of the reference rule, and
      ! the 5-point Gauss-Legendre nodes and weights on [-1,1]
      ! (Abramowitz & Stegun Table 25.4).
      integer, parameter :: n_gsub = 16
      integer, parameter :: n_gl   = 5
      real*8, parameter :: x_gl(n_gl) = (/ -9.0617984593866400d-1,          &
                                           -5.3846931010568309d-1,          &
                                            0.0d0,                          &
                                            5.3846931010568309d-1,          &
                                            9.0617984593866400d-1 /)
      real*8, parameter :: w_gl(n_gl) = (/  2.3692688505618909d-1,          &
                                            4.7862867049936647d-1,          &
                                            5.6888888888888889d-1,          &
                                            4.7862867049936647d-1,          &
                                            2.3692688505618909d-1 /)

      ! sigma_2 at the n = 2 edge, the value excited_hydrogen states
      ! (Osterbrock & Ferland 2006 hydrogenic n = 2 continuum).
      real*8, parameter :: s2_thr = 1.4d-17

      ! The synthetic table, as loaded_sed_table_semantics.f90 lays it out.
      real*8, parameter :: lam_top = 3648.0d0
      real*8, parameter :: lam_bot = 9.5d0
      real*8, parameter :: d_base  = 1.0d0

      real*8  :: b_lo(mx), b_hi(mx), lam(mx), f_row(mx)
      integer :: nb
      real*8  :: lam_th(n_th)

      integer :: ires, nfail
      real*8  :: dev_g(3), dev_h(3), dev_g1, dev_h1
      character(len=512) :: root

      nfail = 0

      lam_th(1) = hp_eV*c_light*1.0d8/e_th_HI
      lam_th(2) = hp_eV*c_light*1.0d8/e_th_HeI
      lam_th(3) = hp_eV*c_light*1.0d8/e_th_HeII
      lam_th(4) = hp_eV*c_light*1.0d8/e_th_HeTR
      lam_th(5) = hp_eV*c_light*1.0d8/e_th_HI_n2

      call get_environment_variable('EXHALE_TEST_ROOT', root)
      if (len_trim(root) .eq. 0) root = '../../../..'

      ! --- the run-wide input, as input_read resolves benchmarks/wasp52 ---
      sp_type      = 'Load'
      is_PL_sed    = .false.
      do_read_sed  = .true.
      is_monochr   = .false.
      thereis_Xray = .true.
      e_mid        = 123.98d0
      e_top        = 1.24d3
      PLind        = -1.0d0
      LX           = 28.0d0
      LEUV         = 29.0d0
      a_orb        = 0.0272d0*AU
      T_star_eff   = 5000.0d0
      R_star       = 0.79d0*Rsun
      appx_mth     = 'Rate/2 + Mdot/2'
      a_tau        = 0.0d0
      allocate(melem_ab(n_melem))
      melem_ab            = 0.0d0
      thereis_metals      = .false.
      thereis_lowIP_metal = .false.
      thereis_HeITR       = .true.
      thereis_mol         = .false.
      use_excited_H       = .true.

      write(*,'(A,ES13.6,A,F9.3,A)') '  band ', e_th_HI_n2, ' eV to ',      &
         e_th_HI, ' eV'

      ! ---- the synthetic narrow-line table, at three resolutions --------
      sed_file = 'synthetic_sed.txt'
      do ires = 1,3
         call build_narrow_line_table(2**(ires-1))
         call load_and_measure('narrow_line', ires, dev_g(ires),            &
                               dev_h(ires), nfail)
      enddo
      ! The rule under test must not lose accuracy as the table is refined:
      ! that is the signature the fixed-step rule showed.
      call converges_verdict('narrow_line_quadrature_converges_rate',       &
                             dev_g(1), dev_g(3), nfail)
      call converges_verdict('narrow_line_quadrature_converges_heat',       &
                             dev_h(1), dev_h(3), nfail)

      ! ---- the measured eps Eri table of benchmarks/wasp52 --------------
      sed_file = trim(root)//'/inputdata/sed/wasp52b_epseri_yan2022.txt'
      call load_and_measure('wasp52', 1, dev_g1, dev_h1, nfail)

      ! ---- the two analytic types ---------------------------------------
      sp_type     = 'Power-law'
      is_PL_sed   = .true.
      do_read_sed = .false.
      e_low       = 13.60d0
      call measure_analytic('power_law', nfail)

      sp_type     = 'Planck'
      is_PL_sed   = .false.
      call measure_analytic('planck', nfail)

      if (nfail .gt. 0) then
         write(*,'(A,I0,A)') 'balmer_band_quadrature: ', nfail,             &
                             ' assertion(s) FAILED'
         call exit(1)
      endif
      write(*,'(A)') 'balmer_band_quadrature: all assertions PASSED'

      contains

      ! ----------------------------------------------------------------- !

      subroutine load_and_measure(label, ires, dg, dh, nfail)
      ! Read the table with the production set_energy_vectors, then measure
      ! the production Balmer integrals of the field it built.
      character(len=*), intent(in)    :: label
      integer,          intent(in)    :: ires
      real*8,           intent(out)   :: dg, dh
      integer,          intent(inout) :: nfail

      if (allocated(e_v))        deallocate(e_v)
      if (allocated(de_v))       deallocate(de_v)
      if (allocated(F_XUV))      deallocate(F_XUV)
      if (allocated(sigma_tab))  deallocate(sigma_tab)
      if (allocated(s_h2_di))    deallocate(s_h2_di)
      if (allocated(s_h2_dd))    deallocate(s_h2_dd)
      if (allocated(s_h2_nd))    deallocate(s_h2_nd)
      if (allocated(e_sed_node)) deallocate(e_sed_node)
      if (allocated(F_sed_node)) deallocate(F_sed_node)

      ! read_sed lowers e_low to the photon-grid floor, so the input value
      ! is restored before every call.
      e_low = 13.60d0
      call set_energy_vectors

      call measure(label, ires, dg, dh, nfail)

      end subroutine load_and_measure

      subroutine measure_analytic(label, nfail)
      ! The analytic types state a field at every energy and build no table;
      ! the photon grid is rebuilt so that the run is complete, and the same
      ! measurement is made.
      character(len=*), intent(in)    :: label
      integer,          intent(inout) :: nfail
      real*8 :: dg, dh

      if (allocated(e_v))        deallocate(e_v)
      if (allocated(de_v))       deallocate(de_v)
      if (allocated(F_XUV))      deallocate(F_XUV)
      if (allocated(sigma_tab))  deallocate(sigma_tab)
      if (allocated(s_h2_di))    deallocate(s_h2_di)
      if (allocated(s_h2_dd))    deallocate(s_h2_dd)
      if (allocated(s_h2_nd))    deallocate(s_h2_nd)
      if (allocated(e_sed_node)) deallocate(e_sed_node)
      if (allocated(F_sed_node)) deallocate(F_sed_node)

      e_low = 13.60d0
      call set_energy_vectors
      call measure(label, 1, dg, dh, nfail)

      end subroutine measure_analytic

      subroutine measure(label, ires, dg, dh, nfail)
      ! The four rules on the field the run now has, and the verdicts.
      character(len=*), intent(in)    :: label
      integer,          intent(in)    :: ires
      real*8,           intent(out)   :: dg, dh
      integer,          intent(inout) :: nfail

      real*8 :: cg, ch, rg, rh, rg2, rh2, tg, th, ug, uh
      character(len=32) :: tag

      cg = gamma_n2_balmer()
      ch = heat_n2_balmer()
      call reference_integrals(n_gsub,   rg,  rh)
      call reference_integrals(2*n_gsub, rg2, rh2)
      call uniform_trapezoid(n_trap, tg, th)
      call uniform_trapezoid(n_u400, ug, uh)

      write(tag,'(A,A,I0)') label, '_h', 2**(ires-1)
      if (ires .eq. 1) tag = label

      write(*,'(A,A,A,2ES14.6)') '  ', trim(tag),                           &
         ': reference rate, heating = ', rg, rh
      write(*,'(A,ES11.3,A,ES11.3)') '     DIAGNOSTIC the 400-point'//      &
         ' uniform rule against the reference: rate ', reldev(ug, rg),      &
         '  heating ', reldev(uh, rh)
      write(*,'(A,ES11.3,A,ES11.3)') '     DIAGNOSTIC the 400001-point'//   &
         ' uniform rule against the reference: rate ', reldev(tg, rg),      &
         '  heating ', reldev(th, rh)

      call absolute_verdict(trim(tag)//'_reference_converged',              &
                            max(reldev(rg, rg2), reldev(rh, rh2)),          &
                            0.0d0, 1.0d-12, nfail)
      dg = reldev(cg, rg)
      dh = reldev(ch, rh)
      call absolute_verdict(trim(tag)//'_rate_quadrature', dg, 0.0d0,       &
                            1.0d-6, nfail)
      call absolute_verdict(trim(tag)//'_heat_quadrature', dh, 0.0d0,       &
                            1.0d-6, nfail)

      end subroutine measure

      ! ----------------------------------------------------------------- !

      subroutine reference_integrals(nsub, g, h)
      ! The two band integrals of the run's field by 5-point Gauss-Legendre
      ! on nsub geometric sub-intervals of each interval between consecutive
      ! NODES of the field: the table rows inside the band plus the two band
      ! edges, or the band itself when the run has no table.  The field is
      ! smooth inside one such interval however narrow a line it carries, so
      ! this converges to the exact integral of the field the code reads.
      integer, intent(in)  :: nsub
      real*8,  intent(out) :: g, h
      integer :: k, n_node, i, m
      real*8  :: Ea, Eb, a, b, q, x1, x2

      g  = 0.0d0
      h  = 0.0d0
      Ea = e_th_HI_n2
      Eb = e_th_HI

      if (do_read_sed .and. allocated(e_sed_node)) then
         n_node = size(e_sed_node)
         do k = 1, n_node-1
            a = e_sed_node(k)
            b = e_sed_node(k+1)
            if (k .eq. 1) a = min(a, Ea)
            a = max(a, Ea)
            b = min(b, Eb)
            if (b .le. a) cycle
            call gauss_on(a, b, nsub, g, h)
         enddo
      else
         ! No table: one interval, subdivided the same way.  e_mid is a node
         ! of the power law and is inserted if it falls inside the band.
         a = Ea
         b = Eb
         if (e_mid .gt. Ea .and. e_mid .lt. Eb) then
            call gauss_on(Ea, e_mid, nsub*8, g, h)
            call gauss_on(e_mid, Eb,  nsub*8, g, h)
         else
            call gauss_on(Ea, Eb, nsub*8, g, h)
         endif
      endif

      end subroutine reference_integrals

      subroutine gauss_on(a, b, nsub, g, h)
      ! Accumulate the two integrands over [a,b] with 5-point
      ! Gauss-Legendre on nsub geometric sub-intervals.
      real*8,  intent(in)    :: a, b
      integer, intent(in)    :: nsub
      real*8,  intent(inout) :: g, h
      integer :: i, k
      real*8  :: q, x1, x2, xm, xh, xx, F_E, sig2, w

      q  = (b/a)**(1.0d0/dble(nsub))
      x1 = a
      do i = 1, nsub
         x2 = x1*q
         if (i .eq. nsub) x2 = b
         xm = 0.5d0*(x2 + x1)
         xh = 0.5d0*(x2 - x1)
         do k = 1, n_gl
            xx   = xm + xh*x_gl(k)
            F_E  = stellar_flux_eV(xx)
            sig2 = s2_thr*(e_th_HI_n2/xx)**3.0d0
            w    = w_gl(k)*xh*F_E*sig2
            g    = g + w/(hp_erg*xx*eV2Hz)
            h    = h + w*(1.0d0 - e_th_HI_n2/xx)
         enddo
         x1 = x2
      enddo

      end subroutine gauss_on

      subroutine uniform_trapezoid(npt, g, h)
      ! The two band integrals on npt points uniform in photon energy: with
      ! npt = 400 this is the rule the production functions used before, and
      ! with npt = 400001 the refined uniform rule.
      integer, intent(in)  :: npt
      real*8,  intent(out) :: g, h
      integer :: i
      real*8  :: E, dE, F_E, sig2, tg, th

      dE = (e_th_HI - e_th_HI_n2)/dble(npt-1)
      g  = 0.0d0
      h  = 0.0d0
      do i = 1, npt
         E    = e_th_HI_n2 + dble(i-1)*dE
         F_E  = stellar_flux_eV(E)
         sig2 = s2_thr*(e_th_HI_n2/E)**3.0d0
         tg   = F_E/(hp_erg*E*eV2Hz)*sig2
         th   = tg*(hp_erg*(E - e_th_HI_n2)*eV2Hz)
         if (i .eq. 1 .or. i .eq. npt) then
            tg = 0.5d0*tg
            th = 0.5d0*th
         endif
         g = g + tg
         h = h + th
      enddo
      g = g*dE
      h = h*dE

      end subroutine uniform_trapezoid

      ! ----------------------------------------------------------------- !

      subroutine build_narrow_line_table(mref)
      ! The narrow_line table of loaded_sed_table_semantics.f90, built by the
      ! same construction so that the two drivers measure the same file: a
      ! continuum F_lambda proportional to lambda with two emission lines
      ! narrower than the table step, bin-averaged exactly over bins laid out
      ! downward from lam_top.
      integer, intent(in) :: mref
      integer :: k, i, j, io, ns
      real*8  :: d, e_hi, e_lo, acc, dw, w, xs(n_th+2)

      d    = d_base/dble(mref)
      e_hi = lam_top + 0.5d0*d
      nb   = 0
      do
         e_lo = e_hi - d
         if (e_lo .lt. lam_bot) exit
         nb = nb + 1
         if (nb .gt. mx) error stop 'build_narrow_line_table: nb > mx'
         b_lo(nb) = e_lo
         b_hi(nb) = e_hi
         e_hi     = e_lo
      enddo

      call reverse(b_lo, nb)
      call reverse(b_hi, nb)
      do k = 1,nb
         lam(k) = 0.5d0*(b_lo(k) + b_hi(k))
         call split_at_thresholds(b_lo(k), b_hi(k), xs, ns)
         acc = 0.0d0
         do j = 1,ns-1
            dw = (xs(j+1) - xs(j))/dble(n_sub)
            do i = 1,n_sub
               w   = xs(j) + (dble(i) - 0.5d0)*dw
               acc = acc + f_line(w)*dw
            enddo
         enddo
         f_row(k) = acc/(b_hi(k) - b_lo(k))
      enddo

      open(unit = 31, file = 'synthetic_sed.txt', status = 'replace',       &
           iostat = io)
      if (io .ne. 0) error stop 'build_narrow_line_table: cannot write'
      write(31,'(A)') '# synthetic SED, bin-averaged from an analytic '//   &
                      'F_lambda; wavelength [A], flux [erg cm^-2 s^-1 A^-1]'
      do k = 1,nb
         write(31,'(ES24.16,1X,ES24.16)') lam(k), f_row(k)
      enddo
      close(31)

      end subroutine build_narrow_line_table

      double precision function f_line(w) result(f)
      ! The analytic F_lambda [erg cm^-2 s^-1 A^-1] the table is the bin
      ! average of: a continuum proportional to lambda with two emission
      ! lines narrower than the base step, each carrying fifty times the
      ! continuum flux of one angstrom.  The absolute scale is arbitrary.
      real*8, intent(in) :: w
      real*8, parameter  :: w_ref = 1.0d3

      f = (w/w_ref)                                                         &
        + gauss_line(w, 200.0d0,   0.05d0, 50.0d0*(200.0d0/w_ref))          &
        + gauss_line(w, 1215.67d0, 0.08d0, 50.0d0*(1215.67d0/w_ref))

      end function f_line

      double precision function gauss_line(w, w0, sig, area) result(f)
      real*8, intent(in) :: w, w0, sig, area
      real*8 :: x
      x = (w - w0)/sig
      f = 0.0d0
      if (abs(x) .lt. 40.0d0)                                               &
         f = area/(sig*sqrt(2.0d0*pi))*exp(-0.5d0*x*x)
      end function gauss_line

      subroutine split_at_thresholds(lo, hi, xs, ns)
      ! The sub-interval boundaries of [lo,hi]: its two ends plus every
      ! threshold wavelength strictly inside it, ascending, so that a bin
      ! average is taken on both sides of a threshold.
      real*8,  intent(in)  :: lo, hi
      real*8,  intent(out) :: xs(:)
      integer, intent(out) :: ns
      integer :: i, j
      real*8  :: t
      ns    = 2
      xs(1) = lo
      xs(2) = hi
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

      double precision function reldev(a, b) result(d)
      real*8, intent(in) :: a, b
      d = 0.0d0
      if (b .ne. 0.0d0) then
         d = abs(a - b)/abs(b)
      else if (a .ne. 0.0d0) then
         d = 1.0d0
      endif
      end function reldev

      subroutine absolute_verdict(name, measured, ref, tol, nfail)
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: measured, ref, tol
      integer,          intent(inout) :: nfail
      if (abs(measured - ref) .le. tol) then
         write(*,'(A,A,A,ES11.3,A,ES11.3,A,ES11.3)') 'PASS ', name,         &
            ' measured=', measured, ' reference=', ref, ' tol=', tol
      else
         write(*,'(A,A,A,ES11.3,A,ES11.3,A,ES11.3)') 'FAIL ', name,         &
            ' measured=', measured, ' reference=', ref, ' tol=', tol
         nfail = nfail + 1
      endif
      end subroutine absolute_verdict

      subroutine converges_verdict(name, d_coarse, d_fine, nfail)
      ! The deviation at a quarter of the table step must not exceed the
      ! deviation at the full step.  Round-off is not a failure: a pair
      ! already at 1e-14 carries no information about the rule.
      character(len=*), intent(in)    :: name
      real*8,           intent(in)    :: d_coarse, d_fine
      integer,          intent(inout) :: nfail
      if (d_fine .le. max(d_coarse, 1.0d-13)) then
         write(*,'(A,A,A,ES11.3,A,ES11.3,A)') 'PASS ', name,                &
            ' measured=', d_fine, ' reference=', d_coarse, ' tol=0'
      else
         write(*,'(A,A,A,ES11.3,A,ES11.3,A)') 'FAIL ', name,                &
            ' measured=', d_fine, ' reference=', d_coarse, ' tol=0'
         nfail = nfail + 1
      endif
      end subroutine converges_verdict

      end program balmer_band_quadrature
