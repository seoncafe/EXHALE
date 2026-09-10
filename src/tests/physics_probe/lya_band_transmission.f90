      program lya_band_transmission
      ! Band B2 of the FUV photolysis set is the H I Ly-alpha resonance
      ! line, so the stellar flux that drives it must reach a cell through
      ! the atomic hydrogen column above that cell -- and it must reach the
      ! cell as the MEAN of the penetrating beam over that cell's own
      ! line-centre depth, not as the beam's value at one of its faces.
      !
      ! Production routines exercised: fuv_lw_photon_field of
      ! src/modules/radiation/util_ion_eq.f90 (the shared 912-2304 A beam),
      ! lya_stellar_beam_transmission, lya_stellar_beam_transmission_cell_mean
      ! and lya_line_center_optical_depth of
      ! src/modules/radiation/lya_rt.f90, water_photolysis_rate of
      ! src/modules/lower_atmosphere/water_photolysis.f90.
      !
      ! PHYSICS AND REFERENCES.
      !  * Ly-alpha (1215.67 A) is the resonance line of the dominant
      !    absorber of these atmospheres.  It is longward of the 912 A
      !    Lyman edge, so the H I photoionization CONTINUUM does not reach
      !    it; the line does.  A star-ward neutral column of 1e19-1e21
      !    cm^-2 gives a line-centre optical depth of 1e6-1e8, so the
      !    stellar line reaching a molecular base is not the stellar line
      !    at the planet's orbit.
      !  * The beam is resonantly SCATTERED rather than destroyed, so its
      !    transmission is not exp(-tau) but the fraction of the broad
      !    stellar profile whose wings penetrate to depth tau,
      !    erfc(x1/(sqrt(2) Xs)) with x1 = sqrt(a tau/sqrt(pi)), the depth
      !    a photon of Doppler offset x crosses in one flight.  That is
      !    what lya_rt.f90 already builds its penetrating stellar beam
      !    from, and the two treatments of the one stellar line must be
      !    the same number.  What this driver tests is that identity, not
      !    the expression: x1 is not the escape frequency of the published
      !    static-slab solutions and carries the caveat recorded at the
      !    definition site in lya_rt.f90.
      !  * A photolysis RATE is what the molecules of a cell undergo,
      !    averaged over the cell, and the beam falls across a cell by that
      !    cell's own depth.  The continuum side of the same rate has been
      !    the exact cell mean since the band set existed
      !    (water_photolysis.f90); the line side of band B2 must be one
      !    too, or the two factors of one rate are taken at two different
      !    points of the same cell.  The size of that correction is set by
      !    the change of the erfc ARGUMENT across the cell,
      !    Dx = c[sqrt(tau_in) - sqrt(tau_out)], and NOT by dtau: it is
      !    nothing for a thin cell high in the column and a FACTOR at the
      !    base of a molecular layer, where one cell holds a large part of
      !    the whole atomic hydrogen column.  Both regimes are below.
      !  * The three bands with no line absorber (B1, B3, B4) and the
      !    Lyman-Werner band, whose line absorber is H2, must not move when
      !    the atomic hydrogen column changes.
      !
      ! The configurations below are isothermal H I layers on a uniform
      ! grid.  Sections 1-4 set the H2O and OH densities to zero so that the
      ! continuum optical depth of every band vanishes and each rate is the
      ! free-streaming rate times the line factor alone; that isolates the
      ! quantity under test.  Section 5 adds an H2O column and measures the
      ! photon budget of the band, which is the statement the two factors
      ! are multiplied for.
      !
      ! Tolerances: the identities (thin limit, band isolation) are exact
      ! and are tested to round-off, 1e-12 relative.  The beam mean is a
      ! three-point Gauss-Legendre quadrature and is tested at 1e-6
      ! relative against an INDEPENDENT fine reference (composite Simpson,
      ! 8192 intervals, in the variable u = sqrt(tau)).  The photon budget
      ! of section 5 is tested against the bound its own discretization
      ! implies, dtau_cont/12, and not against a fitted number.  The
      ! statements about the face value are reference-only.
      use global_parameters, only: N, Ng, R0, r, dr_j, opa_pf, appx_mth,   &
                                   thereis_mol, thereis_oxychem,           &
                                   F_Lya_star, F_LW_star,                 &
                                   F_FUV_B3, F_FUV_B4, dv_star_lya,        &
                                   pi, c_light, kb_erg, mu
      use hydrogen_n2_rates, only: nu_lya, C_lya, A_2p1s
      use utils_ion_eq, only: fuv_lw_photon_field, fuv_band_flux
      use water_photolysis, only: n_fuv_band, ib_LW, ib_B2, ib_B3,         &
                                  ib_B4, sigma_H2O_band,                   &
                                  fuv_band_photon_flux
      use lya_rt, only: lya_line_center_optical_depth,                     &
                        lya_stellar_beam_transmission
      use assertion_report
      implicit none

      integer, parameter :: ncell = 40
      real*8,  parameter :: dr_grid = 0.05d0     ! cell width [R0]
      real*8,  parameter :: r0_cm = 1.0d10       ! planet radius [cm]
      real*8,  parameter :: t_layer = 1500.0d0   ! layer temperature [K]
      real*8,  parameter :: nhi_base = 1.0d9     ! H I at the base [cm^-3]
      real*8,  parameter :: h_scale = 0.30d0     ! H I scale height [R0]

      ! The molecular base of a hot Uranus: a star-ward H I column whose
      ! line-centre depth reaches 1e7 at the base, held in an exponential
      ! layer, so that the bottom cells each hold a large fraction of the
      ! depth above them.  This is the regime in which the face value and
      ! the cell mean differ by a factor.
      real*8,  parameter :: h_base = 0.05d0      ! scale height [R0]
      real*8,  parameter :: tau_base_target = 1.0d7

      ! The graded column: cell line-centre depths running from 1 at the top
      ! of the layer to 1e3 at its base.
      real*8,  parameter :: dtau_top = 1.0d0, dtau_bot = 1.0d3

      ! H2O density of the photon-budget column [cm^-3], flat.  It puts a
      ! continuum absorber in the same cells as the line absorber while
      ! leaving the band LINE-dominated, which is the regime the two factors
      ! of the rate are multiplied in.
      real*8,  parameter :: nh2o_cont = 4.0d5

      real*8, allocatable :: nH2(:), nH2O(:), nOH(:), nHI(:), T_K(:)
      real*8, allocatable :: nH_nuc(:)
      real*8, allocatable :: NH2col(:), NH2Ocol(:), NOHcol(:)
      ! CO is the fourth absorber of this beam; this driver carries none,
      ! so its density is zero and its rate comes back zero.
      real*8, allocatable :: nCO(:), NCOcol(:), k_co(:), theta_co(:)
      real*8, allocatable :: f_shield(:), tr_lines(:), k_lw(:)
      real*8, allocatable :: p_single(:), p_absorbed(:)
      real*8, allocatable :: tau_b(:,:), j_h2o(:,:), j_oh(:,:)
      real*8, allocatable :: j_thin(:,:)
      real*8, allocatable :: tau_lya(:), dtau_lya(:), taulya_out(:)
      real*8 :: a_lines_max, col_over_overlap
      real*8 :: free_stream(5), ratio, worst, col_hi
      real*8 :: DnuD, kap_base, wface
      real*8 :: absorbed_mean, absorbed_face, absorbed_ref, cov_bound
      integer :: j, ib, jlo, jhi
      character(len=64) :: label

      ! ---- Grid and run configuration.
      N  = ncell
      R0 = r0_cm
      appx_mth = 'None'
      thereis_mol     = .false.      ! no H2: the LW line transmission is 1
      thereis_oxychem = .true.       ! the band set is solved
      dv_star_lya = 70.0d0           ! stellar Ly-alpha width [km/s], default
      F_LW_star  = 4.809d2
      F_Lya_star = 4.8092d3
      F_FUV_B3   = 6.948d2
      F_FUV_B4   = 6.588d5

      jlo = 1 - Ng
      jhi = N + Ng
      allocate(r(jlo:jhi), dr_j(jlo:jhi), opa_pf(jlo:jhi))
      allocate(nH2(jlo:jhi), nH2O(jlo:jhi), nOH(jlo:jhi), nHI(jlo:jhi))
      allocate(T_K(jlo:jhi), nH_nuc(jlo:jhi))
      allocate(NH2col(jlo:jhi), NH2Ocol(jlo:jhi), NOHcol(jlo:jhi))
      allocate(nCO(jlo:jhi), NCOcol(jlo:jhi), k_co(jlo:jhi),              &
               theta_co(jlo:jhi))
      nCO = 0.0d0
      allocate(f_shield(jlo:jhi), tr_lines(jlo:jhi), k_lw(jlo:jhi))
      allocate(p_single(jlo:jhi), p_absorbed(jlo:jhi))
      allocate(tau_b(jlo:jhi,n_fuv_band), j_h2o(jlo:jhi,n_fuv_band))
      allocate(j_oh(jlo:jhi,n_fuv_band), j_thin(jlo:jhi,n_fuv_band))
      allocate(tau_lya(jlo:jhi), dtau_lya(jlo:jhi), taulya_out(jlo:jhi))

      do j = jlo, jhi
         r(j)    = 1.0d0 + dble(j - 1)*dr_grid
         dr_j(j) = dr_grid
      enddo
      opa_pf = 1.0d0
      T_K    = t_layer
      nH2    = 0.0d0
      nH2O   = 0.0d0        ! no continuum absorber: tau_b = 0 in every band
      nOH    = 0.0d0
      nH_nuc = 0.0d0

      ! Doppler width of the layer, used to set a column to a target
      ! line-centre depth.
      DnuD = nu_lya*sqrt(2.0d0*kb_erg*t_layer/mu)/c_light

      ! Free-streaming rate of each band: what the cell would see with no
      ! absorber above it at all.
      do ib = 1, n_fuv_band
         free_stream(ib) = sigma_H2O_band(ib)                              &
                           *fuv_band_photon_flux(fuv_band_flux(ib), ib)
      enddo

      ! ---- (1) No atomic hydrogen: every band free-streams.
      nHI = 0.0d0
      call solve_field
      j_thin = j_h2o
      call check_relative('lya_band_free_streams_without_hydrogen',        &
                          j_h2o(jlo,ib_B2), free_stream(ib_B2), 1.0d-12)

      ! ---- (2) A thin exponential atomic hydrogen layer above the base.
      do j = jlo, jhi
         nHI(j) = nhi_base*exp(-(r(j) - r(jlo))/h_scale)
      enddo
      call solve_field

      ! The rate of band B2 is the free-streaming rate times the MEAN over
      ! the cell of the stellar line's transmission, cell by cell.  With no
      ! continuum absorber the continuum factor of the rate is exactly 1, so
      ! this reads the line factor alone.
      call band_against_cell_mean('lya_band_carries_beam_cell_mean',       &
                                  'thin_layer')

      ! The suppression is physical and large: the column of this layer is
      ! about 3e18 cm^-2, well inside the 1e19-1e21 range of a real
      ! molecular base, and the beam that reaches the base is a small
      ! fraction of the incident one.
      col_hi = 0.0d0
      do j = jhi-1, jlo, -1
         col_hi = col_hi + 0.5d0*(nHI(j) + nHI(j+1))*(r(j+1) - r(j))*R0
      enddo
      call check_positive('lya_band_test_column_cm2', col_hi)
      ratio = j_h2o(jlo,ib_B2)/free_stream(ib_B2)
      call check_absolute('lya_band_base_transmission_below_half',         &
                          ratio, 0.0d0, 0.5d0)

      ! ---- (3) The bands with no atomic-hydrogen line are untouched.
      do ib = 1, n_fuv_band
         if (ib .eq. ib_B2) cycle
         worst = 0.0d0
         do j = jlo, jhi
            worst = max(worst, abs(j_h2o(j,ib) - j_thin(j,ib))            &
                               /max(j_thin(j,ib), 1.0d-99))
         enddo
         write(label,'(a,a)') 'lya_band_leaves_other_band_unchanged_',    &
              trim(band_label(ib))
         call check_absolute(label, worst, 0.0d0, 1.0d-12)
      enddo

      ! ---- (4) The transmission falls inward, monotonically.
      worst = 0.0d0
      do j = jlo, jhi - 1
         worst = max(worst, j_h2o(j,ib_B2) - j_h2o(j+1,ib_B2))
      enddo
      call check_absolute('lya_band_transmission_falls_inward',            &
                          worst/free_stream(ib_B2), 0.0d0, 1.0d-12)

      ! ---- (5) The molecular base: the regime the correction matters in.
      ! An exponential H I layer of total line-centre depth 1e7, so that the
      ! bottom cells each hold most of the depth above them and the erfc
      ! argument moves by order unity across one cell.
      kap_base = tau_base_target/(h_base*R0)
      do j = jlo, jhi
         nHI(j) = kap_base*DnuD/C_lya*exp(-(r(j) - r(jlo))/h_base)
      enddo
      call solve_field
      call band_against_cell_mean('lya_band_carries_beam_cell_mean',       &
                                  'molecular_base')

      ! ---- (6) The graded column: cell line-centre depths 1 to 1e3.
      do j = jlo, jhi
         nHI(j) = cell_depth_target(j)*DnuD/(C_lya*dr_j(j)*R0)
      enddo
      call solve_field
      call band_against_cell_mean('lya_band_carries_beam_cell_mean',       &
                                  'graded_1_to_1e3')

      ! ---- (7) The photon budget of the band, with a continuum absorber
      ! in the same cells as the line absorber: the extended H I layer of
      ! section 2 plus a flat H2O density.  That layer is chosen for this
      ! measurement because its beam falls across each cell by a fraction
      ! rather than by a factor, which is the regime in which the product of
      ! the two cell means is a good form: the covariance the product drops
      ! is then small, while the face rule's first-order error is not.
      !
      ! The photons cell j takes out of the B2 beam are, as the band is
      ! discretized, the PRODUCT of the two cell means,
      !     dN_j = N_ph <T_lya>_j exp(-tau_c,out)(1 - exp(-dtau_c)) ,
      ! read here from the production rate of that cell.  The exact loss of
      ! the beam across the same cell is the integral of
      ! T_lya(tau_l) exp(-tau_c) over the cell's own continuum depth, with
      ! both depths running linearly across a cell of uniform density.  The
      ! two differ only by the COVARIANCE of the two factors across the
      ! cell, which the product of their means drops.  Both factors fall
      ! monotonically, so the covariance is positive (the product of means
      ! is a lower bound) and is at most a quarter of the product of the two
      ! relative variations across the cell -- for the exponential factor
      ! that variation is exactly dtau_c.  The assertion is that bound,
      ! measured on this column, and not a fitted number.  The face rule
      ! misses the budget by the first-order error of the LINE factor
      ! instead, which this bound does not contain.
      do j = jlo, jhi
         nHI(j) = nhi_base*exp(-(r(j) - r(jlo))/h_scale)
      enddo
      nH2O = nh2o_cont
      call solve_field
      call beam_loss_budget(absorbed_mean, absorbed_face, absorbed_ref,   &
                            cov_bound)
      write(*,'(a,4es14.6)') '# lya_band_budget mean/face/exact/bound ',  &
           absorbed_mean, absorbed_face, absorbed_ref, cov_bound
      call check_relative('lya_band_photon_budget_cell_mean',             &
                          absorbed_mean, absorbed_ref, cov_bound)
      ! Reference-only: the face rule's error against the same budget, which
      ! is first order in the erfc-argument step and is not bounded by the
      ! covariance of the two cell means.
      call check_positive('lya_band_photon_budget_face_error',            &
                          abs(absorbed_face - absorbed_ref)/absorbed_ref  &
                          - cov_bound)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'lya_band_transmission: ',                    &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'lya_band_transmission: all assertions passed'

      contains

      ! ---------------------------------------------------------------- !

      subroutine solve_field
      ! The production 912-2304 A field on the current densities.
      call fuv_lw_photon_field(nH2, nH2O, nOH, nCO, nHI, T_K, nH_nuc,      &
                               NH2col, NH2Ocol, NOHcol, NCOcol,            &
                               f_shield, tr_lines, tau_b,                  &
                               k_lw, p_single, p_absorbed,                 &
                               k_co, theta_co,                            &
                               j_h2o, j_oh, a_lines_max, col_over_overlap)
      call lya_line_center_optical_depth(T_K, nHI, tau_lya, dtau_lya,      &
                                         taulya_out)
      end subroutine solve_field

      ! ---------------------------------------------------------------- !

      subroutine band_against_cell_mean(name, where_at)
      ! The B2 rate of every cell against the fine reference for the mean of
      ! the beam over that cell, with the face value's ratio to the same
      ! mean printed cell by cell beside it.  Valid only while the continuum
      ! depth of the band is zero, which is what makes the rate the line
      ! factor alone.
      character(len=*), intent(in) :: name, where_at
      real*8  :: rw, mn, fc
      integer :: jj
      worst = 0.0d0
      wface = 0.0d0
      write(*,'(a,a)') '# lya_band face/mean per cell, column ', where_at
      do jj = jlo, jhi
         mn = beam_mean_reference(taulya_out(jj), dtau_lya(jj))
         rw = j_h2o(jj,ib_B2)/free_stream(ib_B2)
         worst = max(worst, abs(rw - mn)/max(mn,1.0d-99))
         fc = lya_stellar_beam_transmission(T_K(jj), taulya_out(jj))
         wface = max(wface, fc/max(mn,1.0d-99))
         write(*,'(a,i5,4es14.6)') '#   j,dtau,mean,face,face/mean ',      &
              jj, dtau_lya(jj), mn, fc, fc/max(mn,1.0d-99)
      enddo
      write(label,'(a,a,a)') name, '_', trim(where_at)
      call check_absolute(label, worst, 0.0d0, 1.0d-6)
      write(label,'(a,a)') 'lya_band_face_over_mean_', trim(where_at)
      call check_positive(label, wface - 1.0d0)
      end subroutine band_against_cell_mean

      ! ---------------------------------------------------------------- !

      subroutine beam_loss_budget(a_mean, a_face, a_ref, bound)
      ! Photons the continuum absorbers of the column take out of band B2,
      ! per unit area and per unit incident photon flux.
      !
      !   a_mean  what the band actually removes, READ BACK FROM THE
      !           PRODUCTION RATE: the H2O rate of a cell is
      !           sigma N_ph (line factor) <exp(-tau_c)>, so
      !           j_H2O dtau_c/(sigma N_ph) is the photons that cell takes
      !           out of the beam per incident photon, whatever the band
      !           uses for its line factor.
      !   a_face  the same sum with the star-ward FACE value of the line
      !           transmission in place of its cell mean.
      !   a_ref   the exact loss: a fine integration of
      !           T_lya(tau_l) exp(-tau_c) over each cell's own continuum
      !           depth, both depths running linearly across a cell because
      !           the densities are piecewise constant.
      !   bound   the covariance bound of the product-of-means form, a
      !           quarter of the product of the two relative variations
      !           across the cell, maximized over the column.
      real*8, intent(out) :: a_mean, a_face, a_ref, bound
      integer, parameter :: nsub = 4096
      real*8  :: t_c_out, d_c, t_l_out, d_l, s, f, hstep, tl, tc, integ
      real*8  :: mn, fc, fin, nph
      integer :: jj, i
      a_mean = 0.0d0
      a_face = 0.0d0
      a_ref  = 0.0d0
      bound  = 0.0d0
      nph    = fuv_band_photon_flux(fuv_band_flux(ib_B2), ib_B2)
      do jj = jlo, jhi
         if (jj .lt. jhi) then
            t_c_out = tau_b(jj+1,ib_B2)
         else
            t_c_out = 0.0d0
         endif
         d_c = max(tau_b(jj,ib_B2) - t_c_out, 0.0d0)
         t_l_out = taulya_out(jj)
         d_l     = dtau_lya(jj)
         mn  = beam_mean_reference(t_l_out, d_l)
         fc  = lya_stellar_beam_transmission(T_K(jj), t_l_out)
         fin = lya_stellar_beam_transmission(T_K(jj), t_l_out + d_l)
         a_mean = a_mean                                                   &
                + j_h2o(jj,ib_B2)*d_c/(sigma_H2O_band(ib_B2)*nph)
         a_face = a_face + fc*exp(-t_c_out)*(1.0d0 - exp(-d_c))
         bound  = max(bound, 0.25d0*d_c*(fc - fin)/max(mn,1.0d-99))
         ! Composite Simpson in the fractional path s across the cell, of
         ! T_lya exp(-tau_c) dtau_c.
         hstep = 1.0d0/dble(nsub)
         f     = 0.0d0
         do i = 0, nsub
            s     = dble(i)*hstep
            tl    = t_l_out + s*d_l
            tc    = t_c_out + s*d_c
            integ = lya_stellar_beam_transmission(T_K(jj), tl)             &
                    *exp(-tc)*d_c
            if (i .eq. 0 .or. i .eq. nsub) then
               f = f + integ
            else if (mod(i,2) .eq. 1) then
               f = f + 4.0d0*integ
            else
               f = f + 2.0d0*integ
            endif
         enddo
         a_ref = a_ref + f*hstep/3.0d0
      enddo
      end subroutine beam_loss_budget

      ! ---------------------------------------------------------------- !

      double precision function cell_depth_target(j_in) result(d)
      ! Line-centre depth asked of cell j_in: a geometric ramp from
      ! dtau_bot at the base of the grid to dtau_top at its top.
      integer, intent(in) :: j_in
      real*8 :: frac
      frac = dble(j_in - jlo)/dble(jhi - jlo)
      d = dtau_bot*(dtau_top/dtau_bot)**frac
      end function cell_depth_target

      ! ---------------------------------------------------------------- !

      double precision function beam_mean_reference(t_o, d) result(m)
      ! Fine reference for the mean over one cell of the production stellar
      ! beam transmission: composite Simpson with 8192 intervals in
      ! u = sqrt(tau), where the integrand 2 u T_s(u^2) is smooth.  An
      ! independent rule from the three-point Gauss-Legendre the production
      ! routine uses, and the same reference lya_beam_cell_mean applies to
      ! that routine directly.
      real*8, intent(in) :: t_o, d
      integer, parameter :: nint_ref = 8192
      real*8  :: u0, u1, h, u, s, f
      integer :: i
      if (d .le. 0.0d0) then
         m = lya_stellar_beam_transmission(t_layer, t_o)
         return
      endif
      u0 = sqrt(max(t_o,0.0d0))
      u1 = sqrt(max(t_o,0.0d0) + d)
      h  = (u1 - u0)/dble(nint_ref)
      s  = 0.0d0
      do i = 0, nint_ref
         u = u0 + dble(i)*h
         f = 2.0d0*u*lya_stellar_beam_transmission(t_layer, u*u)
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

      character(len=8) function band_label(ib_in)
      integer, intent(in) :: ib_in
      if (ib_in .eq. ib_LW) then
         band_label = 'LW'
      else if (ib_in .eq. ib_B3) then
         band_label = 'B3'
      else if (ib_in .eq. ib_B4) then
         band_label = 'B4'
      else
         band_label = 'B2'
      endif
      end function band_label

      end program lya_band_transmission
