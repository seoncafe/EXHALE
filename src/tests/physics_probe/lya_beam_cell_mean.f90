      program lya_beam_cell_mean
      ! The stellar Ly-alpha field a cell is pumped by is the MEAN of the
      ! penetrating stellar beam over that cell's own line-centre optical
      ! depth, not the beam's value at one of the cell's faces.
      !
      ! Production routines exercised:
      ! src/modules/radiation/lya_rt.f90 --
      !   lya_line_center_optical_depth (the cell-by-cell line-centre depths),
      !   lya_stellar_beam_transmission (the point transmission),
      !   lya_stellar_beam_transmission_cell_mean (its mean over a cell),
      !   lya_photosphere_attenuation_cell_mean (the mean of the
      !   parameterized 1/(1+tau) attenuation the jlya_mode = 0 field uses),
      !   jlya_escape_prob (the field the H(n=2) pumping reads).
      !
      ! PHYSICS AND REFERENCES.
      !  * Ly-alpha is resonantly SCATTERED, so the stellar beam is not
      !    attenuated by exp(-tau): the broad stellar line penetrates
      !    through its own wings, and the fraction of it that reaches
      !    line-centre depth tau is T_s = erfc(c sqrt(tau)) with
      !    c = sqrt(a/sqrt(pi))/(sqrt(2) Xs), the one-flight penetration
      !    condition recorded at its definition site in lya_rt.f90.  The
      !    transmission is therefore NOT exponential in the column and its
      !    cell mean has no closed form.
      !  * The rate a cell undergoes is what its atoms undergo averaged over
      !    the cell.  The depth runs linearly across a cell of uniform
      !    absorber density, so that rate carries
      !      <T_s> = (1/dtau) int_{tau_out}^{tau_out+dtau} T_s(tau) dtau ,
      !    and the sum over the column of <T_s> dtau is then exactly the
      !    depth integral of the beam, which is the number of scatterings
      !    the column performs.  A face value is one-sided at every cell in
      !    the same direction and misses that budget.
      !  * The size of the correction is set by the change of the erfc
      !    ARGUMENT across the cell, Dx = c[sqrt(tau_out+dtau) -
      !    sqrt(tau_out)], not by dtau itself: a cell of line-centre depth
      !    10 sitting deep in a column of depth 1e7 changes x by 1e-6 and
      !    the beam not at all, while the base cell of an exponential layer
      !    holds a large fraction of the whole column and changes x by
      !    order unity.  Both regimes are tested below.
      !
      ! Tolerances: the depth identities and the closed-form 1/(1+tau) mean
      ! are exact and are tested to round-off, 1e-12 relative.  The beam
      ! mean is a three-point Gauss-Legendre quadrature and is tested at
      ! 1e-6 relative against an INDEPENDENT fine reference (composite
      ! Simpson, 8192 intervals, in the variable u = sqrt(tau)).  The
      ! statements about the face value are reference-only.
      use global_parameters, only: N, Ng, R0, r, dr_j, v0, pi, appx_mth,   &
                                   F_Lya_star, dv_star_lya,                &
                                   lya_star_boost, lya_bottom_absorber,    &
                                   gamma2_bal, dayside_dilution, c_light,  &
                                   kb_erg, mu
      use hydrogen_n2_rates, only: nu_lya, C_lya, A_2p1s
      use lya_rt, only: lya_rt_allocate_arrays, jstar_arr,                 &
                        lya_line_center_optical_depth,                     &
                        lya_stellar_beam_transmission,                     &
                        lya_stellar_beam_transmission_cell_mean,           &
                        lya_photosphere_attenuation_cell_mean,             &
                        jlya_escape_prob
      use assertion_report
      implicit none

      integer, parameter :: ncell   = 40
      real*8,  parameter :: dr_grid = 0.05d0     ! cell width [R0]
      real*8,  parameter :: r0_cm   = 1.0d10     ! planet radius [cm]
      real*8,  parameter :: t_layer = 8.0d3      ! layer temperature [K]
      real*8,  parameter :: h_scale = 0.05d0     ! H I scale height [R0]
      real*8,  parameter :: tau_base_target = 1.0d7   ! column at the base

      integer, parameter :: ncase = 4
      real*8, parameter :: dtau_case(ncase) = (/ 0.3d0, 1.0d0, 3.0d0,      &
                                                 10.0d0 /)

      real*8, allocatable :: nHI(:), T_K(:), nHII(:), ne(:), v_in(:)
      real*8, allocatable :: tau(:), dtau(:), tau_out(:), Jlya(:)
      real*8 :: pref, xi, Dnu_star, DnuD, kap_base, tau_col
      real*8 :: tbar, tbar2, want, want2, worst, wface
      real*8 :: budget_mean, budget_face, budget_ref, dx_worst
      real*8 :: face, x_step, t_out
      integer :: j, jlo, jhi, ic
      character(len=64) :: label

      ! ---- Run configuration: the stellar beam alone.
      N  = ncell
      R0 = r0_cm
      v0 = 1.0d5
      appx_mth   = 'None'
      dv_star_lya = 70.0d0        ! stellar Ly-alpha half-width [km/s]
      F_Lya_star  = 4.8092d3      ! [erg cm^-2 s^-1]
      ! Boost 1 removes the quadratic trapping buildup, so that J_star is the
      ! beam mean times a constant and this test reads the beam alone; the
      ! quadratic term is checked directly against its own reference below.
      lya_star_boost      = 1.0d0
      lya_bottom_absorber = .false.
      gamma2_bal          = 0.0d0

      jlo = 1 - Ng
      jhi = N + Ng
      allocate(r(jlo:jhi), dr_j(jlo:jhi))
      allocate(nHI(jlo:jhi), T_K(jlo:jhi), nHII(jlo:jhi), ne(jlo:jhi))
      allocate(v_in(jlo:jhi), tau(jlo:jhi), dtau(jlo:jhi))
      allocate(tau_out(jlo:jhi), Jlya(jlo:jhi))
      call lya_rt_allocate_arrays

      do j = jlo, jhi
         r(j)    = 1.0d0 + dble(j - 1)*dr_grid
         dr_j(j) = dr_grid
      enddo
      T_K  = t_layer
      ! No electrons and no protons: the internal (recombination and
      ! collisional) source of J_lya vanishes and J_lya is the stellar beam.
      nHII = 0.0d0
      ne   = 0.0d0
      v_in = 0.0d0

      xi       = dayside_dilution()
      Dnu_star = nu_lya*(dv_star_lya*1.0d5)/c_light
      pref     = xi*F_Lya_star/(4.0d0*pi*Dnu_star)

      ! ------------------------------------------------------------------ !
      ! (1) The line-centre depths are the cell's own.
      ! ------------------------------------------------------------------ !
      ! An exponential H I layer, the profile a molecular base carries: the
      ! bottom cells hold most of the column, which is where a cell's own
      ! depth is a large fraction of the depth above it.
      DnuD     = nu_lya*sqrt(2.0d0*kb_erg*t_layer/mu)/c_light
      kap_base = tau_base_target/(h_scale*R0)
      do j = jlo, jhi
         nHI(j) = kap_base*DnuD/C_lya*exp(-(r(j) - r(jlo))/h_scale)
      enddo

      call lya_line_center_optical_depth(T_K, nHI, tau, dtau, tau_out)

      worst = 0.0d0
      do j = jlo, jhi
         worst = max(worst, abs(dtau(j)                                    &
                 - nHI(j)*C_lya/DnuD*dr_j(j)*R0)/max(dtau(j),1.0d-99))
      enddo
      call check_absolute('lya_cell_depth_is_own_absorber', worst, 0.0d0,  &
                          1.0d-12)
      call check_absolute('lya_top_cell_face_depth_is_zero',               &
                          tau_out(jhi), 0.0d0, 1.0d-99)
      call check_relative('lya_top_cell_carries_its_own_depth',            &
                          tau(jhi), dtau(jhi), 1.0d-12)
      worst = 0.0d0
      do j = jlo, jhi - 1
         worst = max(worst, abs(tau_out(j) - tau(j+1))                     &
                            /max(tau(j+1),1.0d-99))
      enddo
      call check_absolute('lya_faces_of_neighbouring_cells_agree', worst,  &
                          0.0d0, 1.0d-12)

      ! ------------------------------------------------------------------ !
      ! (2) The field the H(n=2) pumping reads is the cell mean.
      ! ------------------------------------------------------------------ !
      call jlya_escape_prob(T_K, nHI, nHII, ne, v_in, Jlya, tau)

      worst = 0.0d0
      wface = 0.0d0
      dx_worst = 0.0d0
      do j = jlo, jhi
         want = beam_mean_reference(t_layer, tau_out(j), dtau(j), 1)
         worst = max(worst, abs(jstar_arr(j)/pref - want)                  &
                            /max(want,1.0d-99))
         face  = lya_stellar_beam_transmission(t_layer, tau_out(j))
         wface = max(wface, abs(face - want)/max(want,1.0d-99))
         x_step = erfc_argument(t_layer, tau(j))                           &
                - erfc_argument(t_layer, tau_out(j))
         dx_worst = max(dx_worst, x_step)
      enddo
      call check_absolute('lya_pumping_field_is_cell_mean', worst, 0.0d0,  &
                          1.0d-6)
      ! Reference-only: how far the star-ward face value is from that mean
      ! in this column, and the erfc-argument step that sets it.
      call check_positive('lya_face_value_error_exceeds_half',             &
                          wface - 0.5d0)
      call check_positive('lya_erfc_argument_step_of_base_cell', dx_worst)

      ! ------------------------------------------------------------------ !
      ! (3) The scattering budget of the column.
      ! ------------------------------------------------------------------ !
      ! sum_j <T_s>_j dtau_j is the depth integral of the beam, which is the
      ! number of scatterings the column performs; the face rule is not.
      budget_mean = 0.0d0
      budget_face = 0.0d0
      do j = jlo, jhi
         budget_mean = budget_mean + jstar_arr(j)/pref*dtau(j)
         budget_face = budget_face                                         &
                     + lya_stellar_beam_transmission(t_layer, tau_out(j))  &
                       *dtau(j)
      enddo
      tau_col    = tau(jlo)
      budget_ref = beam_mean_reference(t_layer, 0.0d0, tau_col, 1)*tau_col
      call check_relative('lya_column_scattering_budget_cell_mean',        &
                          budget_mean, budget_ref, 1.0d-6)
      call check_positive('lya_column_scattering_budget_face_error',       &
                          abs(budget_face - budget_ref)/budget_ref - 0.01d0)

      ! ------------------------------------------------------------------ !
      ! (4) The beam mean at the four requested cell depths.
      ! ------------------------------------------------------------------ !
      ! Each case is tested at the top of the column (tau_out = 0, where the
      ! square-root behaviour of erfc(c sqrt(tau)) is strongest) and at the
      ! depth where the beam is half absorbed, x = c sqrt(tau) = 1.
      do ic = 1, ncase
         t_out = 0.0d0
         call one_cell_case(dtau_case(ic), t_out, 'top')
         t_out = (1.0d0/erfc_argument(t_layer, 1.0d0))**2
         call one_cell_case(dtau_case(ic), t_out, 'photosphere')
      enddo

      ! ------------------------------------------------------------------ !
      ! (5) The parameterized photosphere factor of the jlya_mode = 0 field.
      ! ------------------------------------------------------------------ !
      ! Its cell mean is the closed form ln[(1+tau_in)/(1+tau_out)]/dtau, and
      ! the sum of <1/(1+tau)> dtau over the column telescopes to
      ! ln(1 + tau_total) exactly.
      worst = 0.0d0
      wface = 0.0d0
      budget_mean = 0.0d0
      do j = jlo, jhi
         want = photosphere_mean_reference(tau_out(j), dtau(j))
         tbar = lya_photosphere_attenuation_cell_mean(tau_out(j), dtau(j))
         worst = max(worst, abs(tbar - want)/want)
         wface = max(wface, abs(1.0d0/(1.0d0 + tau(j)) - want)/want)
         budget_mean = budget_mean + tbar*dtau(j)
      enddo
      call check_absolute('lya_photosphere_factor_is_cell_mean', worst,    &
                          0.0d0, 1.0d-12)
      call check_relative('lya_photosphere_factor_column_identity',        &
                          budget_mean, log(1.0d0 + tau(jlo)), 1.0d-12)
      ! Reference-only: the inner-face value's error against that mean.
      call check_positive('lya_photosphere_face_value_error', wface)
      ! The thin limit is the unattenuated field.
      call check_relative('lya_photosphere_factor_thin_limit',             &
                     lya_photosphere_attenuation_cell_mean(0.0d0,1.0d-14), &
                     1.0d0, 1.0d-12)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'lya_beam_cell_mean: ',                       &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'lya_beam_cell_mean: all assertions passed'

      contains

      ! ---------------------------------------------------------------- !

      subroutine one_cell_case(d, t_o, where_at)
      ! One cell of line-centre depth d whose star-ward face is at depth
      ! t_o: the production mean of the beam and of its square against the
      ! fine reference, and the face value's error stated beside them.
      real*8, intent(in) :: d, t_o
      character(len=*), intent(in) :: where_at

      call lya_stellar_beam_transmission_cell_mean(t_layer, t_o, d,        &
                                                   tbar, tbar2)
      want  = beam_mean_reference(t_layer, t_o, d, 1)
      want2 = beam_mean_reference(t_layer, t_o, d, 2)
      write(label,'(a,f5.2,a,a)') 'lya_beam_mean_dtau', d, '_',            &
           trim(where_at)
      call check_relative(label, tbar, want, 1.0d-6)
      write(label,'(a,f5.2,a,a)') 'lya_beam_mean_square_dtau', d, '_',     &
           trim(where_at)
      call check_relative(label, tbar2, want2, 1.0d-6)
      ! Reference-only: the relative error of the star-ward face value of
      ! the beam against that mean, for this cell.
      write(label,'(a,f5.2,a,a)') 'lya_beam_face_error_dtau', d, '_',      &
           trim(where_at)
      call check_positive(label,                                           &
           abs(lya_stellar_beam_transmission(t_layer, t_o) - want)/want)

      end subroutine one_cell_case

      ! ---------------------------------------------------------------- !

      double precision function erfc_argument(t_gas, tau_in) result(x)
      ! x = c sqrt(tau), the argument of the wing-penetration transmission
      ! erfc(x), written out here from the atomic data rather than taken
      ! from the production routine, so that the depth at which the beam is
      ! half absorbed is an independent statement of this test.
      real*8, intent(in) :: t_gas, tau_in
      real*8 :: avoigt, vth, Xs, DnuD_l
      DnuD_l = nu_lya*sqrt(2.0d0*kb_erg*max(t_gas,1.0d0)/mu)/c_light
      avoigt = (A_2p1s/(4.0d0*pi))/DnuD_l
      vth    = sqrt(2.0d0*kb_erg*max(t_gas,1.0d0)/mu)
      Xs     = (dv_star_lya*1.0d5)/vth
      x      = sqrt(avoigt*max(tau_in,0.0d0)/sqrt(pi))/(sqrt(2.0d0)*Xs)
      end function erfc_argument

      ! ---------------------------------------------------------------- !

      double precision function beam_mean_reference(t_gas, t_o, d, ipow)   &
                               result(m)
      ! Fine reference for the mean over a cell of the production
      ! transmission (ipow = 1) or of its square (ipow = 2): composite
      ! Simpson with 8192 intervals in u = sqrt(tau), where the integrand
      ! 2 u T_s(u^2)^ipow is smooth.  An independent rule from the
      ! three-point Gauss-Legendre the production routine uses.
      real*8, intent(in) :: t_gas, t_o, d
      integer, intent(in) :: ipow
      integer, parameter :: nint_ref = 8192
      real*8  :: u0, u1, h, u, s, f
      integer :: i
      if (d .le. 0.0d0) then
         m = lya_stellar_beam_transmission(t_gas, t_o)**ipow
         return
      endif
      u0 = sqrt(max(t_o,0.0d0))
      u1 = sqrt(max(t_o,0.0d0) + d)
      h  = (u1 - u0)/dble(nint_ref)
      s  = 0.0d0
      do i = 0, nint_ref
         u = u0 + dble(i)*h
         f = 2.0d0*u*lya_stellar_beam_transmission(t_gas, u*u)**ipow
         if (i .eq. 0 .or. i .eq. nint_ref) then
            s = s + f
         else if (mod(i,2) .eq. 1) then
            s = s + 4.0d0*f
         else
            s = s + 2.0d0*f
         endif
      enddo
      m = s*h/3.0d0/d
      end function beam_mean_reference

      ! ---------------------------------------------------------------- !

      double precision function photosphere_mean_reference(t_o, d)         &
                               result(m)
      ! Fine reference for the mean over a cell of the parameterized
      ! photosphere factor 1/(1+tau): composite Simpson with 8192 intervals
      ! in tau.  Independent of the closed form the production routine
      ! evaluates, and free of the cancellation that
      ! log[(1+tau_in)/(1+tau_out)] suffers when the two depths are large
      ! and close.
      real*8, intent(in) :: t_o, d
      integer, parameter :: nint_ref = 8192
      real*8  :: h, t, s, f
      integer :: i
      if (d .le. 0.0d0) then
         m = 1.0d0/(1.0d0 + max(t_o,0.0d0))
         return
      endif
      h = d/dble(nint_ref)
      s = 0.0d0
      do i = 0, nint_ref
         t = max(t_o,0.0d0) + dble(i)*h
         f = 1.0d0/(1.0d0 + t)
         if (i .eq. 0 .or. i .eq. nint_ref) then
            s = s + f
         else if (mod(i,2) .eq. 1) then
            s = s + 4.0d0*f
         else
            s = s + 2.0d0*f
         endif
      enddo
      m = s*h/3.0d0/d
      end function photosphere_mean_reference

      end program lya_beam_cell_mean
