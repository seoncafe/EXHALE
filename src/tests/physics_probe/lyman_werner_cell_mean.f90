      program lyman_werner_cell_mean
      ! The Lyman-Werner photodissociation rate a cell carries must be the
      ! MEAN of the local rate over that cell, as the H2O and OH rates of
      ! the same beam already are.
      !
      ! Production routines exercised: fuv_lw_photon_field of
      ! src/modules/radiation/util_ion_eq.f90 (the shared 912-2304 A beam),
      ! lyman_werner_dissociation_rate and
      ! lyman_werner_dissociation_rate_cell_mean of
      ! src/modules/lower_atmosphere/lyman_werner.f90,
      ! h2_lw_dissociation_cross_section of
      ! src/modules/lower_atmosphere/h2_self_shielding_table.f90.
      !
      ! PHYSICS AND REFERENCES.
      !  * calc_column_dens_one accumulates from the top down, so the
      !    column NH2col(j) already contains the whole of cell j: it is the
      !    column at that cell's INNER (planet-ward) face, and NH2col(j+1)
      !    is the column at its star-ward face.  The rate the cell carries
      !    is therefore not the rate at NH2col(j) but the mean of the local
      !    rate between the two faces.
      !  * Within a cell the H2 density is uniform, which is the rectangle
      !    rule the column integration itself uses, so the column runs
      !    linearly across the cell and the mean over the cell's radial
      !    extent is the mean over its column interval.  That is the
      !    reference integral evaluated below, with 1000 midpoints of the
      !    production local rate.
      !  * The cross section falls by four decades across the
      !    self-shielding transition (Draine & Bertoldi 1996 sec. 5.2 for
      !    the mechanism; the level-resolved table of
      !    h2_self_shielding_table for the values), monotonically with
      !    column, so the inner-face value is the SMALLEST rate anywhere in
      !    the cell and using it for the whole cell is low, in one
      !    direction, in every cell where the column grows.
      !
      ! The configuration is an isothermal exponential H2 layer on a
      ! uniform grid whose cell width is 2.5 scale heights, so that the
      ! column grows by e^2.5 across each cell of the layer and the cross
      ! section falls by more than a decade across the deepest ones: the
      ! regime the H2 front is in.  The oxygen chemistry is off, so the
      ! continuum depth of the band vanishes and the line term is what is
      ! being tested.
      !
      ! Section (6) is a second quantity of the same beam: the fraction of
      ! the band the H2 lines remove, which is the column integral of
      ! sigma_diss/p_eff.  It identifies the wavelength interval the
      ! tabulated pumping was computed over, and prints the photon budget
      ! of the band that follows from it.
      !
      ! Tolerances: 1e-3 relative on the quadrature, which is the stated
      ! accuracy of the composite Gauss rule in the production function
      ! (measured there to 2.5e-4 over the whole table); 1e-12 on the
      ! identity that a cell with no column growth carries the local rate;
      ! and, for the size of the correction, the factor by which the
      ! inner-face value falls below the mean in this configuration, which
      ! is a measured number of this test and is pinned to 1 per cent.
      use global_parameters, only: N, Ng, R0, r, dr_j, opa_pf, appx_mth,   &
                                   thereis_mol, thereis_oxychem,           &
                                   F_Lya_star, F_LW_star,                 &
                                   F_FUV_B3, F_FUV_B4
      use utils_ion_eq, only: fuv_lw_photon_field, fuv_band_flux
      use water_photolysis, only: n_fuv_band, ib_LW
      use lyman_werner_photodissociation, only:                           &
                        lyman_werner_dissociation_rate,                    &
                        lyman_werner_dissociation_rate_cell_mean,          &
                        h2_lw_dissociation_cross_section,                  &
                        h2_lw_dissociation_per_absorbed_photon,            &
                        h2_lw_band_photon_fraction_absorbed,               &
                        lyman_werner_band_absorption_rate_cell_mean,       &
                        h2_shield_max_column, e_lw_photon_erg
      use assertion_report
      implicit none

      integer, parameter :: ncell = 40
      real*8,  parameter :: dr_grid = 0.05d0     ! cell width [R0]
      real*8,  parameter :: r0_cm = 1.0d10       ! planet radius [cm]
      real*8,  parameter :: t_layer = 1300.0d0   ! layer temperature [K]
      real*8,  parameter :: nh2_base = 5.0d11    ! H2 at the base [cm^-3]
      real*8,  parameter :: h_scale = 0.020d0    ! H2 scale height [R0]
      real*8,  parameter :: n_nuc = 1.0d13       ! hydrogen nuclei [cm^-3]
      integer, parameter :: nref = 1000          ! reference midpoints
      ! The nine temperatures of the self-shielding table's own grid, so
      ! that the bound on the band share is checked where the table is
      ! defined and not only where this driver runs.
      real*8,  parameter :: t_grid_probe(9) =                             &
           (/ 700.0d0, 900.0d0, 1100.0d0, 1300.0d0, 1600.0d0,             &
              1800.0d0, 2200.0d0, 2700.0d0, 3200.0d0 /)

      real*8, allocatable :: nH2(:), nH2O(:), nOH(:), nHI(:), T_K(:)
      real*8, allocatable :: nH_nuc(:)
      real*8, allocatable :: NH2col(:), NH2Ocol(:), NOHcol(:)
      ! CO is the fourth absorber of this beam; this driver carries none,
      ! so its density is zero and its rate comes back zero.
      real*8, allocatable :: nCO(:), NCOcol(:), k_co(:), theta_co(:)
      real*8, allocatable :: f_shield(:), tr_lines(:), k_lw(:)
      real*8, allocatable :: p_single(:), p_absorbed(:)
      real*8, allocatable :: tau_b(:,:), j_h2o(:,:), j_oh(:,:)
      real*8 :: a_lines_max, col_over_overlap
      real*8 :: flux_lw, ref, worst, worst_ratio, dev, ratio
      real*8 :: col_in, col_out, sig_in, sig_out, decade_span
      real*8 :: syn_out, syn_in, syn_lo, syn_hi, syn_mean, syn_ref
      real*8 :: dr_cm, n_beam, rated_ph, beam_loss_ph, a_top_lo, a_top_hi
      integer :: j, jlo, jhi, jworst, it


      ! ---- Grid and run configuration.
      N  = ncell
      R0 = r0_cm
      appx_mth = 'None'
      thereis_mol     = .true.       ! H2 lines absorb
      thereis_oxychem = .false.      ! no continuum absorber in the band
      F_LW_star  = 4.809d2           ! the mol_lyman_werner case value
      F_Lya_star = 0.0d0
      F_FUV_B3   = 0.0d0
      F_FUV_B4   = 0.0d0

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
      allocate(j_oh(jlo:jhi,n_fuv_band))

      do j = jlo, jhi
         r(j)    = 1.0d0 + dble(j - 1)*dr_grid
         dr_j(j) = dr_grid
      enddo
      opa_pf = 1.0d0
      T_K    = t_layer
      nH2O   = 0.0d0
      nOH    = 0.0d0
      nHI    = 0.0d0
      nH_nuc = n_nuc
      do j = jlo, jhi
         nH2(j) = nh2_base*exp(-(r(j) - r(jlo))/h_scale)
      enddo

      call fuv_lw_photon_field(nH2, nH2O, nOH, nCO, nHI, T_K, nH_nuc,      &
                               NH2col, NH2Ocol, NOHcol, NCOcol,            &
                               f_shield, tr_lines, tau_b,                  &
                               k_lw, p_single, p_absorbed,                 &
                               k_co, theta_co,                            &
                               j_h2o, j_oh, a_lines_max, col_over_overlap)
      flux_lw = fuv_band_flux(ib_LW)

      ! ---- (1) The configuration is in the regime under test: the cross
      ! section falls by close to a decade across the deepest cell.
      col_in  = NH2col(jlo)
      col_out = NH2col(jlo+1)
      sig_in  = h2_lw_dissociation_cross_section(col_in,  t_layer, n_nuc)
      sig_out = h2_lw_dissociation_cross_section(col_out, t_layer, n_nuc)
      decade_span = log10(sig_out/sig_in)
      call check_positive('lw_deepest_cell_drop_above_0p8_dex_margin',     &
                          decade_span - 0.8d0)

      ! ---- (2) Every cell of the layer carries the mean of the local rate
      ! over its own column interval, to the accuracy of the quadrature.
      ! The outermost cell is excluded: its star-ward column is zero and
      ! the uniform midpoint reference used here does not resolve the
      ! cross section over the decades below the table's lowest column.
      worst  = 0.0d0
      jworst = jlo
      do j = jlo, jhi-1
         ref = reference_cell_mean(NH2col(j+1), NH2col(j))
         dev = abs(k_lw(j) - ref)/max(ref, 1.0d-99)
         if (dev .gt. worst) then
            worst  = dev
            jworst = j
         endif
      enddo
      call check_absolute('lw_rate_is_the_cell_mean', worst, 0.0d0,        &
                          1.0d-3)

      ! ---- (3) The size of the correction: the inner-face value, which is
      ! the rate at the largest column anywhere in the cell, is below the
      ! mean everywhere and by a factor of several in the deepest cells.
      worst_ratio = 0.0d0
      do j = jlo, jhi-1
         ratio = lyman_werner_dissociation_rate(flux_lw, NH2col(j),        &
                     t_layer, n_nuc, 0.0d0)/max(k_lw(j), 1.0d-99)
         worst_ratio = max(worst_ratio, ratio)
         if (j .eq. jlo) then
            ! MEASURED (2026-09-05) in the deepest cell of this
            ! configuration, whose cross section falls by 1.12 decades
            ! across the cell: the inner-face value is 3.18 times below
            ! the mean.  Before the correction this ratio was 1 by
            ! construction, the rate being the inner-face value itself.
            call check_within_factor('lw_inner_face_rate_below_the_mean',  &
                                     ratio, 3.1406d-1, 1.01d0)
         endif
      enddo
      ! One-signed: nowhere above the mean.
      call check_absolute('lw_inner_face_rate_never_above_the_mean',       &
                          max(worst_ratio - 1.0d0, 0.0d0), 0.0d0, 1.0d-12)

      ! ---- (4) A cell across which the column does not grow carries the
      ! local rate exactly: the mean of a constant is that constant.
      syn_out = 3.0d17
      call check_relative('lw_cell_mean_of_a_uniform_column',              &
           lyman_werner_dissociation_rate_cell_mean(flux_lw, syn_out,      &
                syn_out, t_layer, n_nuc, 0.0d0, 0.0d0),                    &
           lyman_werner_dissociation_rate(flux_lw, syn_out, t_layer,       &
                n_nuc, 0.0d0), 1.0d-12)

      ! ---- (5) The stated test of the item, on a synthetic cell built to
      ! span exactly one decade of cross section: bisect for the inner-face
      ! column that gives sigma(N_out)/sigma(N_in) = 10.
      syn_out = 1.0d16
      sig_out = h2_lw_dissociation_cross_section(syn_out, t_layer, n_nuc)
      syn_lo  = syn_out
      syn_hi  = 1.0d22
      do it = 1, 200
         syn_in = sqrt(syn_lo*syn_hi)
         sig_in = h2_lw_dissociation_cross_section(syn_in, t_layer, n_nuc)
         if (sig_out/sig_in .gt. 1.0d1) then
            syn_hi = syn_in
         else
            syn_lo = syn_in
         endif
      enddo
      syn_in = sqrt(syn_lo*syn_hi)
      sig_in = h2_lw_dissociation_cross_section(syn_in, t_layer, n_nuc)
      call check_relative('lw_synthetic_cell_spans_one_decade',            &
                          sig_out/sig_in, 1.0d1, 1.0d-6)
      syn_mean = lyman_werner_dissociation_rate_cell_mean(flux_lw,         &
                     syn_out, syn_in, t_layer, n_nuc, 0.0d0, 0.0d0)
      syn_ref  = reference_cell_mean(syn_out, syn_in)
      call check_relative('lw_synthetic_cell_mean_matches_reference',      &
                          syn_mean, syn_ref, 1.0d-3)

      ! ---- (6) THE BEAM LOSES WHAT THE TABLE RATES.
      !
      ! sigma_pump = sigma_diss/p_eff is the cross section for taking a
      ! band photon OUT of the beam (h2_self_shielding_table header:
      ! sigma_diss = p_eff x sigma_pump, and p_eff is the dissociations per
      ! photon removed from the beam).  Its column integral is therefore
      ! the fraction of the incident beam the H2 lines have absorbed,
      !
      !     A(N) = int_0^N sigma_pump(N') dN' ,
      !
      ! and written out over frequency it is int dnu (F_nu/h nu)
      ! (1 - e^-tau(nu,N)) divided by the band photon flux, so for a beam
      ! confined to the band it CANNOT EXCEED 1: the lines can at most take
      ! every photon there is.
      !
      ! Three statements are asserted here.
      !
      ! (i) The production closed form is that integral.
      ! h2_lw_band_photon_fraction_absorbed integrates the module's own
      ! interpolant exactly, knot interval by knot interval; the probe
      ! integrates the two public accessors numerically on 20000
      ! logarithmic steps.  The two are independent evaluations of one
      ! quantity and must agree to the quadrature error of the second.
      !
      ! (ii) A stays below 1 at the top of the column axis.  Until
      ! 2026-09-06 it did not: MEASURED then, it reached 1.5154, because
      ! the table's line list ran to 1200 A while its normalization band
      ! stopped at 1110 A, and 1.5154 is the photon content of 912-1200 A
      ! over that of 912-1110 A for the flat-F_lambda spectrum the table is
      ! built on.  The band is now 912-1201 A on both sides.
      !
      ! (iii) The photon budget of the beam closes: the line absorptions
      ! the model rates over the column equal what the beam loses.  With no
      ! continuum absorber in the band -- this driver runs with the oxygen
      ! chemistry off -- the two are the same integral written two ways,
      ! and the residual is the cell-mean quadrature of the rate against
      ! the closed-form column integral.
      a_top_lo = pump_column_integral(h2_shield_max_column(), 3200.0d0,   &
                                      1.0d12)
      a_top_hi = pump_column_integral(h2_shield_max_column(), 3200.0d0,   &
                                      1.0d14)
      ! The trapping is the only density axis of sigma_diss and p_eff
      ! carries the same factor, so their ratio has none: sigma_pump is a
      ! line-by-line quantity of (N_H2, T) alone.  The 1e-9 is the
      ! round-off of the 20001-term quadrature sum, not a physical spread.
      call check_relative('lw_pump_cross_section_has_no_density_axis',    &
                          a_top_lo, a_top_hi, 1.0d-9)
      ! (i) the production closed form against this quadrature.
      call check_relative('lw_band_fraction_matches_column_quadrature',   &
                          h2_lw_band_photon_fraction_absorbed(            &
                              h2_shield_max_column(), 3200.0d0, 1.0d14),  &
                          a_top_hi, 1.0d-4)
      call check_relative('lw_band_fraction_matches_quadrature_mid',      &
                          h2_lw_band_photon_fraction_absorbed(            &
                              2.0d20, t_layer, n_nuc),                    &
                          pump_column_integral(2.0d20, t_layer, n_nuc),   &
                          1.0d-4)
      ! (ii) the share of the band cannot pass 1, at the top of the axis
      ! and at every node of the temperature grid.
      call check_absolute('lw_band_fraction_at_axis_top_below_one',       &
                          a_top_hi, 0.5d0, 0.5d0)
      worst = 0.0d0
      do it = 1, 9
         worst = max(worst,                                               &
              h2_lw_band_photon_fraction_absorbed(                        &
                  h2_shield_max_column(), t_grid_probe(it), 1.0d13))
      enddo
      call check_absolute('lw_band_fraction_below_one_on_the_T_grid',     &
                          worst, 0.5d0, 0.5d0)

      ! (iii) the photon budget of the beam.
      dr_cm  = dr_grid*r0_cm
      n_beam = flux_lw/e_lw_photon_erg
      rated_ph = 0.0d0
      do j = jlo, jhi
         if (j .lt. jhi) then
            col_out = NH2col(j+1)
         else
            col_out = 0.0d0
         endif
         rated_ph = rated_ph + nH2(j)*dr_cm                               &
                  * lyman_werner_band_absorption_rate_cell_mean(          &
                        flux_lw, col_out, NH2col(j), t_layer, n_nuc,      &
                        0.0d0, 0.0d0)
      enddo
      beam_loss_ph = n_beam*(1.0d0 - tr_lines(jlo))
      write(*,'(a,es12.5,a,es12.5,a,f9.6)')                               &
           'note: LW line photons rated ', rated_ph,                      &
           ' against a beam loss of ', beam_loss_ph,                      &
           ', ratio ', rated_ph/beam_loss_ph
      ! The tolerance is the quadrature of the cell mean, not a fit: the
      ! rated side is the three-point Gauss-Legendre rule on geometric
      ! segments (lyman_werner.f90 measures it at 2.5e-4 against a
      ! 4000-segment reference) and the beam side is the exact column
      ! integral of the same interpolant.  With sigma_diss/p_absorbed at
      ! the cell's own column in place of this cell mean, the same sum
      ! stood 8.0 per cent above the beam loss on these cells, which span
      ! 2.5 H2 scale heights each.
      call check_relative('lw_rated_photons_equal_the_beam_loss',         &
                          rated_ph, beam_loss_ph, 1.0d-3)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'lyman_werner_cell_mean: ',                   &
              assertion_failures, ' assertion(s) failed'
         write(*,'(a,i0,a,es12.5)') 'note: worst cell j=', jworst,         &
              ' relative departure from the reference mean ', worst
         flush(6)
         stop 1
      endif
      write(*,'(a,es12.5)') 'note: worst relative departure from the'//    &
           ' 1000-point reference mean ', worst
      write(*,'(a)') 'lyman_werner_cell_mean: all assertions passed'

      contains

      double precision function reference_cell_mean(col_star_ward,         &
                                col_inner) result(kbar)
      ! Mean of the production local rate over one cell, by the midpoint
      ! rule on nref equal steps of the column between the two faces.  The
      ! column is linear in radius inside a cell, so this is the mean over
      ! the cell's radial extent.
      real*8, intent(in) :: col_star_ward, col_inner
      real*8  :: dcol, colx
      integer :: i
      dcol = col_inner - col_star_ward
      kbar = 0.0d0
      do i = 1, nref
         colx = col_star_ward + dcol*(dble(i) - 0.5d0)/dble(nref)
         kbar = kbar + lyman_werner_dissociation_rate(flux_lw, colx,       &
                           t_layer, n_nuc, 0.0d0)
      enddo
      kbar = kbar/dble(nref)
      end function reference_cell_mean

      double precision function pump_column_integral(N_H2, T, n_H)        &
                                result(A)
      ! int_0^N sigma_pump dN' with sigma_pump = sigma_diss/p_eff, by the
      ! trapezoidal rule on 20000 logarithmic steps of the column between
      ! 1 cm^-2 and N_H2, plus the column below 1 cm^-2 over which the
      ! table returns its edge value and the integrand is constant.  The
      ! integrand N sigma_pump(N) varies by four decades over the range, so
      ! the logarithmic rule is what converges: doubling the step count
      ! moves the result by under 1e-6.
      real*8, intent(in) :: N_H2, T, n_H
      real*8  :: dl, Nx, term
      integer :: i
      integer, parameter :: nlog = 20000
      dl = log10(max(N_H2, 1.0d0))/dble(nlog)
      A  = 0.0d0
      do i = 0, nlog
         Nx   = 10.0d0**(dble(i)*dl)
         term = h2_lw_dissociation_cross_section(Nx, T, n_H)              &
              / h2_lw_dissociation_per_absorbed_photon(Nx, T, n_H)        &
              * Nx*log(10.0d0)
         if (i .eq. 0 .or. i .eq. nlog) term = 0.5d0*term
         A = A + term*dl
      enddo
      A = A + h2_lw_dissociation_cross_section(1.0d-30, T, n_H)           &
            / h2_lw_dissociation_per_absorbed_photon(1.0d-30, T, n_H)
      end function pump_column_integral

      end program lyman_werner_cell_mean
