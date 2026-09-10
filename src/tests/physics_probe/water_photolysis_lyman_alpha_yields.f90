      program water_photolysis_lyman_alpha_yields
      ! The H2O photodissociation branching on the band that contains
      ! Ly-alpha, against the published yields it transcribes.
      !
      ! Production data exercised: qy_H2O_OH_H, qy_H2O_H2_O1D and
      ! qy_H2O_O_H_H on the fuv_band table of
      ! src/modules/lower_atmosphere/oxygen_rates.f90.  The three-body
      ! branch is open in no other band, and it is the only photolytic
      ! source of ground-state O in the option, so a shift of the triplet
      ! moves where O and OH come from at the base.
      !
      ! REFERENCE.
      !   Slanger & Black (1982), J. Chem. Phys. 77, 2432, measured H(2S)
      !   and O(1D, 3P) yields from H2O photolysis at 1216 A.  Their
      !   Eqs. (1)-(4) are the four channels
      !
      !       H2O -> OH(X 2Pi) + H(2S)            (1)
      !           -> OH(A 2Sigma+) + H(2S)        (2)
      !           -> O(1D) + H2(X 1Sigma_g+)      (3)
      !           -> O(3P) + 2H(2S)               (4)
      !
      !   and p. 2435 reads "We may thus set yields of 78% for processes
      !   1 + 2, 10% for process 3, and 12% for process 4".  The 0.78 is
      !   therefore Phi1 + Phi2, not a single ground-state channel, and the
      !   three are not independent: the 10% is adopted from Stief et al.,
      !   the 12% is derived against it, and the 78% is what normalization
      !   leaves.
      !
      ! WHAT IS ASSERTED.
      !   The three yields of the Ly-alpha band, digit for digit, and that
      !   they sum to one, which is what makes the merged OH channel a
      !   partition of the absorptions rather than a subset of them.  The
      !   band that carries them is required to contain 1215.67 A and to be
      !   the only band with the three-body branch open, since the yields
      !   are a single-wavelength measurement and mean nothing on a band
      !   that does not contain the line.
      !
      !   Tolerance is round-off, 1e-12 relative: these are transcribed
      !   constants, not fits.

      use oxygen_rates, only: n_fuv_band, fuv_band_lo_A, fuv_band_hi_A,   &
                              qy_H2O_OH_H, qy_H2O_H2_O1D, qy_H2O_O_H_H
      use assertion_report

      implicit none

      ! Slanger & Black (1982) p. 2435, processes 1 + 2, 3 and 4.
      real*8, parameter :: phi_OH_H = 0.78d0, phi_H2_O1D = 0.10d0,        &
                           phi_O_H_H = 0.12d0
      ! Ly-alpha [A], the wavelength of the measurement.
      real*8, parameter :: lambda_lya_A = 1215.67d0

      integer :: b, band_lya, n_open
      real*8  :: total

      ! Locate the band that contains the line, and count the bands in
      ! which the three-body branch is open.
      band_lya = 0
      n_open   = 0
      do b = 1, n_fuv_band
         if (lambda_lya_A .ge. fuv_band_lo_A(b) .and.                     &
             lambda_lya_A .le. fuv_band_hi_A(b)) band_lya = b
         if (qy_H2O_O_H_H(b) .gt. 0.0d0) n_open = n_open + 1
      enddo

      call check_relative('lyman_alpha_band_count', dble(min(band_lya,1)),&
                          1.0d0, 1.0d-12)
      call check_relative('bands_with_three_body_branch_open',            &
                          dble(n_open), 1.0d0, 1.0d-12)

      if (band_lya .eq. 0) then
         write(*,'(a)') 'water_photolysis_lyman_alpha_yields: no band '// &
              'contains Ly-alpha'
         flush(6)
         stop 1
      endif

      call check_relative('h2o_photolysis_yield_oh_plus_h_lyman_alpha',   &
                          qy_H2O_OH_H(band_lya), phi_OH_H, 1.0d-12)
      call check_relative('h2o_photolysis_yield_h2_plus_o1d_lyman_alpha', &
                          qy_H2O_H2_O1D(band_lya), phi_H2_O1D, 1.0d-12)
      call check_relative('h2o_photolysis_yield_o_plus_2h_lyman_alpha',   &
                          qy_H2O_O_H_H(band_lya), phi_O_H_H, 1.0d-12)

      total = qy_H2O_OH_H(band_lya) + qy_H2O_H2_O1D(band_lya)             &
            + qy_H2O_O_H_H(band_lya)
      call check_relative('h2o_photolysis_yield_sum_lyman_alpha',         &
                          total, 1.0d0, 1.0d-12)

      write(*,'(a,i0,a,f7.1,a,f7.1,a)')                                    &
           'note: lyman_alpha_band_index = ', band_lya, ', edges = ',      &
           fuv_band_lo_A(band_lya), ' to ', fuv_band_hi_A(band_lya), ' A'

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'water_photolysis_lyman_alpha_yields: ',      &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)')                                                       &
           'water_photolysis_lyman_alpha_yields: all assertions passed'

      end program water_photolysis_lyman_alpha_yields
