      program photoevent_energy_ledger
      ! The energy an isolated photoabsorption event has to account for.
      !
      ! Origin: docs/audit_20260905/photoevent_energy_check.py.  That script
      ! evaluates a SUPERSEDED equation as a counterexample and is therefore
      ! not an assertion about the code; it is not ported.  The equation
      ! asserted here is the accepted one of
      ! docs/development_plan_20260905_rev3.md section 4.2 item 3 and
      ! section 4.3 (A5), sharpened by docs/b1_target_system_20260906.md
      ! T1.4 and AT-1b into a statement about every recipient of the photon:
      !
      !     e_reservoir + e_electron + e_fragment + e_radiated  =  h nu ,
      !
      ! for one event in an isolated cell, with
      !
      !   e_reservoir  the formation and ionization energy stored in the
      !                products, on the T1.2 reference of neutral,
      !                ground-state, free atoms at rest,
      !   e_electron   the kinetic energy of the photoelectron(s), which
      !                electron_energy_degradation divides into heat,
      !                secondary ionizations and escaping line radiation,
      !   e_fragment   the kinetic energy of the heavy fragments of a
      !                dissociative channel, prompt heat in full,
      !   e_radiated   internal excitation of a product that no evolved
      !                population carries, which leaves as prompt line
      !                radiation.
      !
      ! Nothing else absorbs the photon: there is no radiation field in the
      ! isolated-event limit, so any remainder is energy the accounting has
      ! lost.
      !
      ! Production quantities exercised:
      !   photoelectron_share  (src/modules/radiation/util_ion_eq.f90), the
      !     share of the photon energy the code hands to the photoelectron
      !     of an ATOMIC absorber; the same function the H, He and metal
      !     heating integrals call.
      !   h2_channel_energy_recipients, h2_channel_electron_share,
      !     h2_channel_prompt_heat_share, h2_channel_radiated_share
      !     (src/modules/functions/h2_photo_channels.f90), the single
      !     definition of what each H2 channel charges each recipient.
      !   e_th_HI, e_th_HeI, e_th_HeII, e_th_HeTR, e_th_H2, e_th_H2_di,
      !     e_th_H2_dd (src/modules/init/parameters.f90) and mion_ethr
      !     (src/modules/init/species_table.f90), the thresholds those calls
      !     pass.
      !   h2_dissociation_energy_eV  (src/modules/lower_atmosphere/
      !     mol_rates.f90), the single definition of D0(H2) in this code.
      !
      ! REFERENCES.
      !  * ATOMIC absorbers (H I, He I, He II, He 2^3S, every photoionizable
      !    metal ion): the products are the next ion stage and a free
      !    electron, so the reservoir is charged the ionization potential
      !    that the same threshold names and the electron carries E - E_th.
      !    Nothing else is a recipient.  For He 2^3S the potential is
      !    measured from the METASTABLE, and the identity
      !    e_th_HeTR = e_th_HeI - E(2^3S) is asserted separately.
      !  * H2 double photoionization at 80 eV: the products are two protons
      !    and two electrons, so
      !      e_reservoir = D0(H2) + 2 I(H) = 31.675 eV
      !    (plan section 4.3), the two electrons carry E - 51.4 eV, and the
      !    remaining 19.725 eV is the kinetic energy release of the two
      !    receding protons -- 51.4 eV is a VERTICAL threshold (Yan,
      !    Sadeghpour & Dalgarno 1998, p. 1048, section 4), so at threshold
      !    the electrons are at rest and the whole vertical-to-asymptotic
      !    difference must be with the nuclei.  Before that recipient was
      !    named the three known ones summed to 60.28 of the 80 eV.
      use global_parameters, only: e_th_HI, e_th_HeI, e_th_HeII,          &
                                   e_th_HeTR, e_th_H2, e_th_H2_di,        &
                                   e_th_H2_dd
      use species_table, only: n_mion, mion_ethr, mion_isphot, mion_name
      use mol_rates, only: h2_dissociation_energy_eV
      use h2_photo_channels, only: n_h2_channels, ICH_M, ICH_S, ICH_D,    &
                                   ICH_N, h2_channel_threshold,           &
                                   h2_channel_name,                       &
                                   h2_channel_energy_recipients,          &
                                   h2_channel_electron_share,             &
                                   h2_channel_prompt_heat_share,          &
                                   h2_channel_radiated_share
      use utils_ion_eq, only: photoelectron_share
      use assertion_report
      implicit none

      real*8, parameter :: e_h  = 20.0d0
      real*8, parameter :: e_h2 = 80.0d0
      ! He 1s2s 3S1 excitation energy, 159855.9743 cm^-1 (NIST ASD), the
      ! value parameters.f90 quotes when it builds e_th_HeTR.
      real*8, parameter :: e_exc_HeTR = 19.819614d0
      ! Where each ledger is tested: just above the threshold, at twice the
      ! threshold, and at 100 eV, the top of the band that carries most of
      ! the absorbed XUV energy in these models.
      integer, parameter :: n_at = 4
      character(len=8), parameter :: at_name(n_at) =                      &
           (/ 'H I     ', 'He I    ', 'He II   ', 'He 2^3S ' /)
      real*8 :: at_eth(n_at)

      real*8 :: thermal, chemical, d0, e, tot
      real*8 :: e_res, e_ele, e_frg, e_rad
      integer :: i, k, ich
      character(len=64) :: label

      at_eth = (/ e_th_HI, e_th_HeI, e_th_HeII, e_th_HeTR /)
      d0     = h2_dissociation_energy_eV()

      ! --- H photoionization, the reference case of the plan ---------- !
      thermal  = photoelectron_share(e_th_HI, e_h)*e_h
      chemical = e_th_HI
      call check_relative('h_photoionization_event_energy_sum',           &
                          thermal + chemical, e_h, 1.0d-12)

      ! --- every atomic absorber, at three photon energies ------------- !
      ! Recipients: the ionization potential to the reservoir, E - E_th to
      ! the photoelectron, nothing else.  A failure here would mean the
      ! heating integrand charges an energy that is not the threshold the
      ! reservoir is credited with.
      do i = 1, n_at
         do k = 1, 3
            e = photon_energy(at_eth(i), k)
            tot = photoelectron_share(at_eth(i), e)*e + at_eth(i)
            write(label,'(a,a,a,i0)') 'atomic_photoevent_sum_',           &
                 trim(at_name(i)), '_E', k
            call check_relative(trim(label), tot, e, 1.0d-12)
         enddo
      enddo

      ! He 2^3S is charged the potential of the METASTABLE, not of the
      ! ground singlet: the reservoir difference eps(He II) - eps(He 2^3S)
      ! is I(He I) - E(2^3S), and the photoelectron threshold must be the
      ! same number or the two disagree by the excitation energy.
      call check_relative('he_triplet_threshold_is_reservoir_difference', &
                          e_th_HeTR, e_th_HeI - e_exc_HeTR, 1.0d-9)

      ! --- every photoionizable metal ion, at three photon energies ---- !
      do i = 1, n_mion
         if (.not. mion_isphot(i)) cycle
         do k = 1, 3
            e = photon_energy(mion_ethr(i), k)
            tot = photoelectron_share(mion_ethr(i), e)*e + mion_ethr(i)
            write(label,'(a,a,a,i0)') 'metal_photoevent_sum_',            &
                 trim(mion_name(i)), '_E', k
            call check_relative(trim(label), tot, e, 1.0d-12)
         enddo
      enddo

      ! --- every H2 channel, at three photon energies ------------------ !
      ! The recipient table is at h2_channel_energy_recipients; this is the
      ! assertion that it is complete.
      do ich = 1, n_h2_channels
         do k = 1, 3
            e = photon_energy(h2_channel_threshold(ich), k)
            call h2_channel_energy_recipients(ich, e, e_res, e_ele,       &
                                              e_frg, e_rad)
            write(label,'(a,i0,a,i0)') 'h2_channel_', ich,                &
                 '_photoevent_sum_E', k
            call check_relative(trim(label),                              &
                                e_res + e_ele + e_frg + e_rad, e, 1.0d-12)
            ! No recipient may be negative: a negative one would mean the
            ! channel is charged more than the photon carries and the
            ! balance is taken out of another recipient.
            write(label,'(a,i0,a,i0)') 'h2_channel_', ich,                &
                 '_recipients_nonnegative_E', k
            call check_absolute(trim(label),                              &
                 min(0.0d0, min(min(e_res,e_ele), min(e_frg,e_rad))),     &
                 0.0d0, 1.0d-30)
            ! The three share functions the heating integrand consumes are
            ! the same partition divided by E.
            write(label,'(a,i0,a,i0)') 'h2_channel_', ich,                &
                 '_shares_sum_E', k
            call check_relative(trim(label),                              &
                 h2_channel_electron_share(ich, e)                        &
                 + h2_channel_prompt_heat_share(ich, e)                   &
                 + h2_channel_radiated_share(ich, e)                      &
                 + e_res/e, 1.0d0, 1.0d-12)
         enddo
      enddo

      ! --- the channel thresholds against their thermochemistry -------- !
      ! (S) is an ASYMPTOTIC threshold: Chung et al.'s 18.076 eV is the
      ! energy of H + H+ + e- measured from H2, so the reservoir charge and
      ! the photoelectron threshold are one number and no fragment kinetic
      ! energy is left over.  The tolerance is the 7.1e-4 eV by which the
      ! rounded published threshold differs from D0 + I(H).
      call check_absolute('h2_single_dissociative_threshold_is_asymptotic',&
                          e_th_H2_di, d0 + e_th_HI, 1.0d-3)

      ! (D) is a VERTICAL threshold and is 19.7 eV above the asymptotic
      ! energy of the same products.  That difference is the recipient the
      ! ledger names: the Coulomb-explosion kinetic energy of the two
      ! protons.  Cross-check against the Franck-Condon picture, two bare
      ! protons at the H2 equilibrium separation R_e = 1.401 a0 repelling
      ! with e^2/R_e = 27.2114/1.401 eV; agreement to 2 per cent is what is
      ! asserted, since 51.4 eV is quoted to three figures.
      call h2_channel_energy_recipients(ICH_D, e_h2, e_res, e_ele,        &
                                        e_frg, e_rad)
      call check_relative('h2_double_coulomb_explosion_energy',           &
                          e_frg, 27.211386d0/1.401d0, 2.0d-2)

      ! --- H2 double photoionization at 80 eV -------------------------- !
      chemical = e_res
      thermal  = e_ele + e_frg

      ! The product energy of the accepted equation, against the value the
      ! plan quotes.  The 0.003 eV difference is e_th_HI = 13.6 eV in the
      ! plan against the 13.598434599 eV of parameters.f90.
      call check_relative('h2_double_chemical_product_energy',            &
                          chemical, 31.675d0, 1.0d-3)

      call check_relative('h2_double_photoevent_energy_sum',              &
                          thermal + chemical, e_h2, 1.0d-12)

      ! What the heating integrand of PH_heat_HHe charges this channel
      ! today, for comparison with the ledger above: the photoelectron
      ! share alone, with no fragment term.  The difference is the
      ! Coulomb-explosion energy.
      write(*,'(a,es23.15)') 'note: energy of one H2 double event at '//  &
           '80 eV not charged by photoelectron_share alone [eV] = ',      &
           e_h2 - photoelectron_share(e_th_H2_dd, e_h2)*e_h2 - chemical

      ! The neutral window, same comparison: the integrand charges
      ! 1 - D0/E and so hands the fragments the n = 2 excitation energy of
      ! both atoms as well as their recoil.
      call h2_channel_energy_recipients(ICH_N, 37.0d0, e_res, e_ele,      &
                                        e_frg, e_rad)
      write(*,'(a,es23.15)') 'note: energy of one H2 neutral-window '//   &
           'event at 37 eV charged to heat but radiated [eV] = ', e_rad

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'photoevent_energy_ledger: ',                &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'photoevent_energy_ledger: all assertions passed'

      contains

      real*8 function photon_energy(e_th, k) result(e)
      ! The three photon energies every ledger is asserted at: one eV above
      ! the threshold, twice the threshold, and 100 eV.
      real*8, intent(in)  :: e_th
      integer, intent(in) :: k
      select case (k)
      case (1)
         e = e_th + 1.0d0
      case (2)
         e = 2.0d0*e_th
      case default
         e = 100.0d0
      end select
      end function photon_energy

      end program photoevent_energy_ledger
