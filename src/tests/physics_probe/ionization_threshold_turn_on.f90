      program ionization_threshold_turn_on
      ! ONE definition of each ionization threshold, and every photoionization
      ! cross section turning on exactly there.
      !
      ! Production data exercised: the threshold constants of module
      ! global_parameters (src/modules/init/parameters.f90), the cross
      ! sections of module Cross_sections (src/modules/functions/cross_sec.f90)
      ! reached through the production dispatcher photoion_sigma
      ! (src/modules/radiation/opacity_models.f90), and the channel thresholds
      ! of module h2_photo_channels (src/modules/functions/h2_photo_channels.f90).
      !
      ! WHY THIS IS A TEST AND NOT A STYLE CHECK.  The photon grid puts every
      ! threshold on a bin EDGE, the photoelectron energy is h nu - e_th and
      ! the collisional-ionization cooling removes e_th per event.  If a cross
      ! section turns on at its own copy of the threshold instead, the band
      ! between the two copies is charged wrongly: a cross section that turns
      ! on ABOVE the grid edge leaves that band integrated as zero (measured
      ! 2026-09-05: e_th_HeI = 24.6 against the He I fit's own E_th = 24.59
      ! lost 0.115 per cent of the He I photoionization rate of the default
      ! power law), and one that turns on BELOW it lets a photon that cannot
      ! ionize be absorbed.
      !
      ! REFERENCE for the values: the measured ionization energies of the NIST
      ! Atomic Spectra Database, and for H2 the adiabatic ionization energy of
      ! the NIST Chemistry WebBook.  The metastable is not a separate
      ! measurement: it is the He I ionization energy minus the excitation
      ! energy of the 1s2s 3S1 level, 159855.9743 cm^-1 = 19.819614 eV, so
      ! that identity is asserted to round-off rather than to a digit count.
      ! Tolerance on the four measured values is 1e-9 relative, the precision
      ! to which they are written down; it is far tighter than the 0.01 to
      ! 0.68 per cent by which the retired round numbers 13.6 / 24.6 / 54.4 /
      ! 4.80 / 15.4 miss them, which is what the test is there to catch.
      !
      ! TURN-ON.  Each cross section is evaluated 1e-6 eV on either side of
      ! its threshold: exactly zero below (tolerance 0, an identity), strictly
      ! positive above.  1e-6 eV is far narrower than any gap this test is
      ! meant to catch (the smallest was 0.003 eV) and far wider than the
      ! double-precision spacing of these energies (7e-15 eV), so the two
      ! probes bracket the threshold without touching it.
      use global_parameters, only: e_th_HI, e_th_HeI, e_th_HeII,          &
                                   e_th_HeTR, e_th_H2,                    &
                                   e_th_H2_di, e_th_H2_dd,                &
                                   e_th_HI_erg, e_th_HeI_erg,             &
                                   e_th_HeII_erg, e_th_HeTR_erg,          &
                                   erg2eV, ates_photoion_rate,            &
                                   opacity_model
      use J_incident,        only: e_th_HI_n2
      use Cross_sections,    only: sigma_H2
      use opacity_models,    only: photoion_sigma
      use h2_photo_channels, only: h2_channel_threshold,                  &
                                   ICH_M, ICH_S, ICH_D
      use assertion_report
      implicit none

      real*8, parameter :: d_probe = 1.0d-6      ! [eV] on either side
      ! NIST Atomic Spectra Database ionization energies [eV] and the NIST
      ! Chemistry WebBook adiabatic ionization energy of H2.
      real*8, parameter :: IE_HI_nist   = 13.598434599d0
      real*8, parameter :: IE_HeI_nist  = 24.587389d0
      real*8, parameter :: IE_HeII_nist = 54.417765d0
      real*8, parameter :: IE_H2_nist   = 15.425927d0
      ! Excitation energy of He I 1s2s 3S1, 159855.9743 cm^-1 (NIST ASD).
      real*8, parameter :: E_exc_He23S  = 19.819614d0

      ! The dispatcher must be on its analytic branch: the tabulated and
      ! constant opacity models are the run's choice, not this test's.
      opacity_model = 'A'
      ates_photoion_rate = .false.

      !---- the values, one definition each ----!

      call check_relative('h_i_threshold_is_the_measured_potential',      &
           e_th_HI, IE_HI_nist, 1.0d-9)
      call check_relative('he_i_threshold_is_the_measured_potential',     &
           e_th_HeI, IE_HeI_nist, 1.0d-9)
      call check_relative('he_ii_threshold_is_the_measured_potential',    &
           e_th_HeII, IE_HeII_nist, 1.0d-9)
      call check_relative('h2_threshold_is_the_measured_potential',       &
           e_th_H2, IE_H2_nist, 1.0d-9)
      call check_relative('he_metastable_threshold_is_singlet_minus_exc', &
           e_th_HeTR, IE_HeI_nist - E_exc_He23S, 1.0d-12)
      ! The H(n=2) edge of the Balmer continuum is the H I threshold of the
      ! n = 2 level, e_th_HI/n^2, not a number of its own.
      call check_relative('h_n2_edge_is_a_quarter_of_the_h_i_threshold',   &
           e_th_HI_n2, e_th_HI/4.0d0, 0.0d0)

      ! The erg forms are the same numbers, converted once.
      call check_relative('h_i_threshold_in_erg',                         &
           e_th_HI_erg*erg2eV, e_th_HI, 1.0d-15)
      call check_relative('he_i_threshold_in_erg',                        &
           e_th_HeI_erg*erg2eV, e_th_HeI, 1.0d-15)
      call check_relative('he_ii_threshold_in_erg',                       &
           e_th_HeII_erg*erg2eV, e_th_HeII, 1.0d-15)
      call check_relative('he_metastable_threshold_in_erg',               &
           e_th_HeTR_erg*erg2eV, e_th_HeTR, 1.0d-15)

      ! The H2 photoionization channels are charged the same thresholds.
      call check_relative('h2_non_dissociative_channel_threshold',        &
           h2_channel_threshold(ICH_M), e_th_H2, 0.0d0)
      call check_relative('h2_dissociative_channel_threshold',            &
           h2_channel_threshold(ICH_S), e_th_H2_di, 0.0d0)
      call check_relative('h2_double_channel_threshold',                  &
           h2_channel_threshold(ICH_D), e_th_H2_dd, 0.0d0)

      !---- the turn-on of each cross section ----!

      call check_absolute('h_i_cross_section_zero_below_threshold',       &
           photoion_sigma('HI', e_th_HI - d_probe), 0.0d0, 0.0d0)
      call check_positive('h_i_cross_section_positive_above_threshold',   &
           photoion_sigma('HI', e_th_HI + d_probe))

      call check_absolute('he_i_cross_section_zero_below_threshold',      &
           photoion_sigma('HeI', e_th_HeI - d_probe), 0.0d0, 0.0d0)
      call check_positive('he_i_cross_section_positive_above_threshold',  &
           photoion_sigma('HeI', e_th_HeI + d_probe))

      call check_absolute('he_ii_cross_section_zero_below_threshold',     &
           photoion_sigma('HeII', e_th_HeII - d_probe), 0.0d0, 0.0d0)
      call check_positive('he_ii_cross_section_positive_above_threshold', &
           photoion_sigma('HeII', e_th_HeII + d_probe))

      call check_absolute('he_23s_cross_section_zero_below_threshold',    &
           photoion_sigma('HeITR', e_th_HeTR - d_probe), 0.0d0, 0.0d0)
      call check_positive('he_23s_cross_section_positive_above_thresh',   &
           photoion_sigma('HeITR', e_th_HeTR + d_probe))

      call check_absolute('h2_cross_section_zero_below_threshold',        &
           sigma_H2(e_th_H2 - d_probe), 0.0d0, 0.0d0)
      call check_positive('h2_cross_section_positive_above_threshold',    &
           sigma_H2(e_th_H2 + d_probe))

      ! The legacy ATES He I fit is the same absorber and turns on at the
      ! same energy; it is selected by the input key `ates_photoion_rate`.
      ates_photoion_rate = .true.
      call check_absolute('he_i_ates_fit_zero_below_threshold',           &
           photoion_sigma('HeI', e_th_HeI - d_probe), 0.0d0, 0.0d0)
      call check_positive('he_i_ates_fit_positive_above_threshold',       &
           photoion_sigma('HeI', e_th_HeI + d_probe))
      ates_photoion_rate = .false.

      ! The constant and tabulated opacity models gate on the same
      ! thresholds.  With no table loaded the tabulated model returns the
      ! analytic cross section, so only the constant model is a new gate.
      opacity_model = 'C'
      call check_absolute('constant_opacity_model_zero_below_threshold',  &
           photoion_sigma('HeI', e_th_HeI - d_probe), 0.0d0, 0.0d0)
      call check_positive('constant_opacity_positive_above_threshold',    &
           photoion_sigma('HeI', e_th_HeI + d_probe))
      opacity_model = 'A'

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'ionization_threshold_turn_on: ',            &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)')                                                      &
           'ionization_threshold_turn_on: all assertions passed'

      end program ionization_threshold_turn_on
