      program absorbed_fraction_switch
      ! QUANTITY UNDER TEST
      !   (1 - exp(-d))/d, the fraction of the photons entering a cell that
      !   the cell absorbs, per unit of its own optical depth d, and the
      !   statement that the XUV beam and the FUV bands take it from ONE
      !   definition, absorbed_fraction_per_unit_depth (utilities.f90).
      !
      !   Production routines exercised: absorbed_fraction_per_unit_depth
      !   and cell_mean_attenuation of src/modules/functions/utilities.f90
      !   (the XUV cell mean of util_ion_eq.f90) and water_photolysis_rate
      !   and hydroxyl_photolysis_rate of
      !   src/modules/lower_atmosphere/water_photolysis.f90 (the FUV band
      !   rates).
      !
      ! WHAT IS ASSERTED
      !   (1) ACCURACY.  The value against a QUADRUPLE-PRECISION evaluation
      !       of (1 - exp(-d))/d, over 15 decades of d.  The closed form
      !       subtracts two numbers that approach each other as d -> 0, so
      !       its relative error is about eps/d with eps = 2.2e-16: it is
      !       1e-8 just above the switch point at d = 1e-8 and falls as d
      !       grows.  Below the switch the series 1 - d/2 + d^2/6 truncates
      !       at d^3/24 in relative terms, 4e-25 there.  The bound over the
      !       whole grid is therefore about 1e-8, attained on the
      !       closed-form side of the switch, and it is asserted at 2e-8;
      !       for d >= 1e-6 the same error is below 1e-9.  A switch point
      !       moved down, or a series branch removed, breaks the first
      !       bound; a switch point moved up far enough that the series
      !       truncation shows breaks neither, which is why the identities
      !       below and not this bound are what pin the shared definition.
      !
      !   (2) ONE DEFINITION, TWO BEAMS.  The FUV band rates and the XUV
      !       cell mean must be BITWISE the same function of d.  With
      !       tau_out = 0 and no line absorber the band rate is
      !       s_b N_b (1 - exp(-d))/d and the XUV cell mean at a_tau = 0 is
      !       (1 - exp(-d))/d, so
      !
      !         water_photolysis_rate(F,ib,0,d,1)
      !             == s_b N_b cell_mean_attenuation(0,d)
      !
      !       holds to the last bit at every d, and likewise for OH.  It
      !       holds because both sides reach the same function through the
      !       same use association; a second copy of the expression, with
      !       any other switch point or series, fails it somewhere in the 15
      !       decades scanned.  This is the assertion that would have caught
      !       the duplicate this test was written with.
      !
      ! GRID.  1501 points logarithmically spaced over d = 1e-12 to 1e3,
      ! plus the switch point itself and its two neighbours in the double
      ! precision grid, so that both branches and the crossing are covered.
      ! d = 0 (a cell of no optical depth, which every empty cell of a run
      ! carries) is checked separately: the fraction is 1 there by the
      ! limit, and the branch returns it exactly.
      use global_parameters, only: a_tau
      use utils, only: absorbed_fraction_per_unit_depth,                  &
                       cell_mean_attenuation
      use water_photolysis, only: water_photolysis_rate,                  &
                                  hydroxyl_photolysis_rate,               &
                                  fuv_band_photon_flux,                   &
                                  sigma_H2O_band, sigma_OH_band, ib_B3
      use assertion_report

      implicit none

      integer, parameter :: npt = 1501
      ! Band flux [erg cm^-2 s^-1] of the point the two rates are formed
      ! at: any positive number does, the identity is in d alone.
      real*8, parameter :: f_band = 6.948d2

      integer :: i, nbad_h2o, nbad_oh
      real*8  :: d, err, err_all, err_hi
      real*8  :: nph_h2o, nph_oh, j_h2o, j_oh, ref_h2o, ref_oh, mean_xuv

      ! The XUV cell mean is the pure exponential only at a_tau = 0 (the
      ! default of the input reader, and the only value for which the field
      ! is exponential in the column); above it the 2D rate correction is
      ! quadratured instead and the identity of (2) is not the statement.
      a_tau = 0.0d0

      err_all  = 0.0d0
      err_hi   = 0.0d0
      nbad_h2o = 0
      nbad_oh  = 0

      nph_h2o = fuv_band_photon_flux(f_band, ib_B3)
      nph_oh  = nph_h2o

      do i = 0,npt+2
         if (i .le. npt) then
            d = 10.0d0**(-12.0d0 + dble(i)*15.0d0/dble(npt))
         else if (i .eq. npt+1) then
            d = 1.0d-8
         else
            d = nearest(1.0d-8, -1.0d0)
         endif

         ! (1) accuracy against the quadruple-precision value.
         err = relative_error_vs_quad(d)
         err_all = max(err_all, err)
         if (d .ge. 1.0d-6) err_hi = max(err_hi, err)

         ! (2) the two beams, bitwise.  The reference repeats the
         ! production expression of water_photolysis_rate with the XUV
         ! cell mean in place of its own factor; nothing else differs.
         mean_xuv = cell_mean_attenuation(0.0d0, d)
         j_h2o    = water_photolysis_rate(f_band, ib_B3, 0.0d0, d, 1.0d0)
         j_oh     = hydroxyl_photolysis_rate(f_band, ib_B3, 0.0d0, d, 1.0d0)
         ref_h2o  = sigma_H2O_band(ib_B3)*nph_h2o*1.0d0*exp(-0.0d0)       &
                    *mean_xuv
         ref_oh   = sigma_OH_band(ib_B3) *nph_oh *1.0d0*exp(-0.0d0)       &
                    *mean_xuv
         if (j_h2o .ne. ref_h2o) nbad_h2o = nbad_h2o + 1
         if (j_oh  .ne. ref_oh ) nbad_oh  = nbad_oh  + 1
      enddo

      call check_absolute('absorbed_fraction_error_vs_quad',              &
                          err_all, 0.0d0, 2.0d-8)
      call check_absolute('absorbed_fraction_error_above_dtau_1em6',      &
                          err_hi, 0.0d0, 1.0d-9)
      call check_absolute('water_band_rate_is_xuv_cell_mean_bitwise',     &
                          dble(nbad_h2o), 0.0d0, 0.0d0)
      call check_absolute('hydroxyl_band_rate_is_xuv_cell_mean_bitwise',  &
                          dble(nbad_oh), 0.0d0, 0.0d0)
      ! The limit at a cell of no optical depth, exactly.
      call check_absolute('absorbed_fraction_at_zero_depth',              &
                          absorbed_fraction_per_unit_depth(0.0d0),        &
                          1.0d0, 0.0d0)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'absorbed_fraction_switch: ',                &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'absorbed_fraction_switch: all assertions passed'

      contains

      !--------------!

      double precision function relative_error_vs_quad(d) result(e)
      ! |f(d) - f_quad(d)|/f_quad(d), with f_quad evaluated in quadruple
      ! precision, where the cancellation of 1 - exp(-d) costs 18 digits
      ! fewer than the value is asked to.
      real*8, intent(in) :: d
      real(kind=16) :: dq, fq
      dq = real(d, kind=16)
      fq = (real(1.0d0, kind=16) - exp(-dq))/dq
      e  = abs(absorbed_fraction_per_unit_depth(d) - real(fq, kind(1.0d0)))&
           /real(fq, kind(1.0d0))
      end function relative_error_vs_quad

      end program absorbed_fraction_switch
