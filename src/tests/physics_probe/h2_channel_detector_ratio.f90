      program h2_channel_detector_ratio
      ! What the measured proton-to-H2+ ratio of Chung et al. (1993) says
      ! about the four mutually exclusive H2 photoabsorption channels, the
      ! invariants every channel obeys, and the yields the consumers of the
      ! channel arrays derive from them.
      !
      ! Origin: docs/audit_20260905/h2_detector_ratio.f90, which prints the
      ! two readings of the ratio side by side.  This program asserts the
      ! one that the published measurement supports.
      !
      ! Production routines exercised: h2_channel_cross_sections and
      ! h2_channel_stoichiometry of src/modules/functions/h2_photo_channels.f90
      ! and frac_H2_dissociative_ionization of
      ! src/modules/functions/cross_sec.f90.
      !
      ! REFERENCE.  Chung et al. (1993), J. Chem. Phys. 99, 885, Table II
      ! tabulates r(E) = sigma(H+)/sigma(H2+), a ratio of two DETECTOR
      ! SIGNALS.  The two ions of a double-ionization event arrive almost
      ! simultaneously and are not resolved: the paper's p. 886 states that
      ! the detector "responds with a single pulse", so one double event
      ! contributes ONE count to the H+ signal, and the tabulated ratio is
      ! the event ratio
      !     (N_S + N_D)/N_M = r
      ! with N_M, N_S, N_D the counts of the H2+, the single dissociative
      ! and the double channel.  (The reading of the published text is the
      ! second review's, docs/development_plan_20260905_review2.md Q1,
      ! against the published version in references/.)  Chung's statement
      ! that "by 80 eV about 20% of sigma(H+) comes from double ionization"
      ! is then a statement about the same signal, i.e. the EVENT fraction
      !     N_D/(N_S + N_D) = 0.20 .
      ! Before the detector inversion of plan item A8 the module read the
      ! same two numbers as a proton-number ratio, (N_S + 2 N_D)/N_M = r
      ! and 2 N_D/(N_S + 2 N_D) = 0.20; the first two assertions below
      ! were the RED that item was gated on (measured at 35d9dd5:
      ! 0.2544530 against 0.2827255, and 0.1111111 against 0.20).
      !
      ! Both tolerances are round-off, not the 4-5 per cent the source
      ! quotes for its data: what is tested is that the channels reproduce
      ! the number the code itself was given, r from the tabulated fit and
      ! 0.20 from the 'chung80' model, not the accuracy of either.
      !
      ! The channel sum and the two closure relations are invariants of the
      ! construction, not measurements, and are tested at round-off.
      !
      ! WHAT ELSE IS ASSERTED HERE.
      !  * the 'yan_rho' model is admissible: its event fraction q_D, built
      !    from Yan et al.'s double-to-single ratio rho_D = N_D/(N_M + N_S)
      !    and the detector ratio r, stays in [0,1] over 51.4-2000 eV, and
      !    reproduces q_D = rho_D/(1 + rho_D) (1 + r)/r at 110 eV, where
      !    Yan et al. quote rho_D = 0.038.  The production routine STOPS
      !    the run on an inadmissible q_D instead of clipping it, so an
      !    inadmissible scan aborts this program rather than printing FAIL.
      !  * the channel arrays as set_energy_vectors.f90 (lines 228-246)
      !    stores them: s_h2_di = sigma_S, s_h2_dd = sigma_D,
      !    s_h2_nd = sigma_N, and the H2+ row that util_ion_eq.f90 (986)
      !    and System_HeH_mol.f90 (279-280) rebuild as
      !    s_h2 - s_h2_di - s_h2_dd - s_h2_nd.  The four must add back to
      !    the opacity, and the proton and electron yields the consumers
      !    write out by hand (System_HeH_mol.f90 245-246: P_H2_di
      !    + 2 P_H2_dd) must equal what h2_channel_stoichiometry gives.
      !    The assignment is mirrored here, not called: set_energy_vectors
      !    needs a configured run.  A change in that assignment is not
      !    caught by this program.
      use Cross_sections, only: frac_H2_dissociative_ionization, sigma_H2
      use h2_photo_channels
      use assertion_report
      implicit none

      real*8, parameter :: e_photon = 80.0d0
      real*8 :: f_di, ratio_measured, sig(n_h2_channels)
      real*8 :: dH2, dH2p, dH, dHp, dele, nuclei, charge
      real*8 :: q_min, q_max, q_scan, e_scan, rho_D, q_yan
      real*8 :: s_tot, s_h2_di, s_h2_dd, s_h2_nd, s_h2_mol
      real*8 :: y_H2, y_Hp, y_ele
      character(len=48) :: label
      integer :: ich

      f_di = frac_H2_dissociative_ionization(e_photon)
      ratio_measured = f_di/(1.0d0 - f_di)

      h2_double_ionization_model = 'chung80'
      call h2_channel_cross_sections(e_photon, 1.0d0, sig)

      ! The detector ratio the channels reproduce, as an event ratio.
      call check_relative('h2_detector_event_ratio',                      &
           (sig(ICH_S) + sig(ICH_D))/sig(ICH_M), ratio_measured, 1.0d-12)

      ! Chung's 20 per cent, as a fraction of the proton-producing events.
      call check_relative('h2_double_event_fraction',                     &
           sig(ICH_D)/(sig(ICH_S) + sig(ICH_D)), 0.20d0, 1.0d-12)

      ! Invariants.
      call check_relative('h2_channel_sum', sum(sig), 1.0d0, 1.0d-12)

      do ich = 1, n_h2_channels
         call h2_channel_stoichiometry(ich, dH2, dH2p, dH, dHp, dele)
         nuclei = 2.0d0*dH2 + 2.0d0*dH2p + dH + dHp
         charge = dH2p + dHp - dele
         write(label,'(a,i0)') 'h2_channel_nuclei_closure_', ich
         call check_absolute(label, nuclei, 0.0d0, 1.0d-12)
         write(label,'(a,i0)') 'h2_channel_charge_closure_', ich
         call check_absolute(label, charge, 0.0d0, 1.0d-12)
      enddo

      ! The double channel as a share of the absorption at 80 eV, against
      ! the value the detector reading gives from Chung's own two numbers,
      ! r(80) from Table II and q_D = 0.20 from the text of p. 888:
      ! sigma_D/sigma_abs = q_D f_di = 0.20 x 0.22041.  The tolerance is
      ! the rounding of the quoted reference, not a fit.
      call check_relative('h2_double_channel_share_80eV',                 &
           sig(ICH_D), 0.044082d0, 1.0d-5)

      ! ---------------------------------------------------------------
      ! The 'yan_rho' model: admissibility over the whole grid, and the
      ! conversion of Yan et al.'s published ratio at the energy they
      ! quote it.  q_D leaves [0,1] where rho_D exceeds r, which is a
      ! restriction on combining the two sources and not a value to clip.
      h2_double_ionization_model = 'yan_rho'
      q_min =  1.0d30
      q_max = -1.0d30
      e_scan = 51.4d0
      do while (e_scan .le. 2000.0d0)
         q_scan = frac_H2_double_of_proton_events(e_scan)
         if (q_scan .lt. q_min) q_min = q_scan
         if (q_scan .gt. q_max) q_max = q_scan
         e_scan = e_scan + 0.01d0
      enddo
      ! [0,1] written as centre 0.5, half-width 0.5, so that the printed
      ! measured value is the extremum itself.
      call check_absolute('h2_yan_rho_double_event_fraction_min',         &
           q_min, 0.5d0, 0.5d0)
      call check_absolute('h2_yan_rho_double_event_fraction_max',         &
           q_max, 0.5d0, 0.5d0)

      ! q_D = rho_D/(1 + rho_D) * (1 + r)/r at 110 eV, r from Table II.
      f_di  = frac_H2_dissociative_ionization(110.0d0)
      rho_D = 0.038d0
      ratio_measured = f_di/(1.0d0 - f_di)
      q_yan = rho_D/(1.0d0 + rho_D)*(1.0d0 + ratio_measured)/ratio_measured
      call check_relative('h2_yan_rho_double_event_fraction_110eV',       &
           frac_H2_double_of_proton_events(110.0d0), q_yan, 1.0d-12)

      ! ---------------------------------------------------------------
      ! The channel arrays as set_energy_vectors stores them, and the
      ! yields their consumers derive, at 80 eV with the default models.
      h2_double_ionization_model = 'chung80'
      s_tot = sigma_H2(e_photon)
      call h2_channel_cross_sections(e_photon, s_tot, sig)
      s_h2_di  = sig(ICH_S)
      s_h2_dd  = sig(ICH_D)
      s_h2_nd  = sig(ICH_N)
      s_h2_mol = s_tot - s_h2_di - s_h2_dd - s_h2_nd

      ! Opacity is the sum: the H2+ row the consumers rebuild by
      ! subtraction is the module's own molecular channel.
      call check_relative('h2_stored_molecular_channel',                  &
           s_h2_mol, sig(ICH_M), 1.0d-12)

      ! One H2 destroyed per event, and the proton and electron yields of
      ! the stoichiometric table against the source terms written by hand.
      y_H2 = 0.0d0; y_Hp = 0.0d0; y_ele = 0.0d0
      do ich = 1, n_h2_channels
         call h2_channel_stoichiometry(ich, dH2, dH2p, dH, dHp, dele)
         y_H2  = y_H2  - dH2 *sig(ich)
         y_Hp  = y_Hp  + dHp *sig(ich)
         y_ele = y_ele + dele*sig(ich)
      enddo
      call check_relative('h2_destroyed_per_absorption', y_H2, s_tot,     &
           1.0d-12)
      call check_relative('h2_proton_yield_of_channels', y_Hp,            &
           s_h2_di + 2.0d0*s_h2_dd, 1.0d-12)
      call check_relative('h2_electron_yield_of_channels', y_ele,         &
           s_h2_mol + s_h2_di + 2.0d0*s_h2_dd, 1.0d-12)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'h2_channel_detector_ratio: ',               &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'h2_channel_detector_ratio: all assertions passed'

      end program h2_channel_detector_ratio
