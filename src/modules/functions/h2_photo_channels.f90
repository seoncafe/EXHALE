      module h2_photo_channels
      ! H2 photoabsorption resolved into MUTUALLY EXCLUSIVE final-state
      ! channels, with one stoichiometric table that every consumer reads.
      !
      ! WHY.  sigma_H2 (cross_sec.f90) is one number and the code used to
      ! split it with a single branching fraction f_di.  Three different
      ! reactions were hiding inside that split -- neutral dissociation over
      ! 33-41 eV, double ionization above 51.4 eV, and the single
      ! dissociative channel it is named for -- so the H2 DESTRUCTION rate
      ! was right while the fragments were not.  This module makes the four
      ! channels explicit:
      !
      !   (M) H2 + hv -> H2+ + e-           threshold 15.4  eV
      !   (S) H2 + hv -> H  + H+ + e-       threshold 18.076 eV
      !   (D) H2 + hv -> H+ + H+ + 2e-      threshold 51.4  eV
      !   (N) H2 + hv -> H  + H             (no threshold of its own; it is
      !                                      the sub-unity photoionization
      !                                      yield over 33-41 eV)
      !
      ! OPACITY IS THE SUM.  sigma_M + sigma_S + sigma_D + sigma_N is
      ! identically sigma_H2(E) by construction, so no consumer of the
      ! opacity changes.  What changes is that the SOURCE terms are now
      ! derived from stoich() below instead of being written out by hand.
      !
      ! ENERGY IS ALSO A SUM.  h2_channel_energy_recipients below names every
      ! recipient of one absorbed photon in each channel -- the formation and
      ! ionization energy of the products, the photoelectron kinetic energy,
      ! the kinetic energy of the heavy fragments, and the internal
      ! excitation that leaves as prompt line radiation -- so that the four
      ! sum to h nu, once each, for every channel.  Its recipient table is
      ! the single statement of what each channel charges, and it carries the
      ! source of every entry.
      !
      ! ------------------------------------------------------------------
      ! HOW THE SPLIT IS BUILT, and what each piece is measured from.
      !
      ! Write the event counts per unit time as N_M, N_S, N_D (channels M,
      ! S, D) and N_N (channel N).  Then
      !
      !   sigma_N = f_n(E) * sigma_H2(E)                             (i)
      !   sigma_M + sigma_S + sigma_D = (1 - f_n(E)) * sigma_H2(E)   (ii)
      !
      ! f_n is Chung et al. (1993) Table 1 as a FRACTION of their own
      ! sigma(abs), carried here as a fraction so that it rides on this
      ! code's sigma_H2 normalization rather than importing theirs.
      !
      ! Inside (ii), Chung et al.'s Table II ratio column
      !
      !   r(E) = sigma(H+)/sigma(H2+)
      !
      ! is the ratio of two TIME-OF-FLIGHT DETECTOR SIGNALS, and that
      ! detector delivers one count per EVENT, not one count per proton.
      ! Chung et al. (1993), p. 886, writing about the two processes
      ! H2 + hv -> H+ + H(nl) and H2 + hv -> H+ + H+ of their Eqs. (1) and
      ! (2): "In each of the above processes our detector responds with a
      ! single pulse (in double ionization both ions arrive almost
      ! simultaneously and cannot be resolved). Thus, above the double
      ! ionization threshold (~47-48 eV) our ratio contains H+ counts from
      ! both processes given in Eqs. (1) and (2)."  The tabulated column is
      ! therefore the EVENT ratio
      !
      !   r(E) = (N_S + N_D)/N_M ,
      !
      ! and with q_D(E) = N_D/(N_S + N_D), the share of the H+-signal
      ! events that releases two protons,
      !
      !   N_M = sigma_ion/(1 + r)
      !   N_S = (1 - q_D) r sigma_ion/(1 + r)
      !   N_D =       q_D r sigma_ion/(1 + r) .
      !
      ! The code below evaluates the equivalent form in the tabulated
      ! fraction f_di = r/(1 + r) = sigma(H+)/[sigma(H+) + sigma(H2+)],
      ! which is what frac_H2_dissociative_ionization returns:
      !
      !   sigma_M = (1 - f_di) sigma_ion
      !   sigma_S = (1 - q_D) f_di sigma_ion
      !   sigma_D =       q_D f_di sigma_ion .
      !
      ! With q_D = 0 this collapses to sigma_M = (1 - f_di) sigma_ion and
      ! sigma_S = f_di sigma_ion, the two-way branching that preceded these
      ! channels.  THAT IS WHAT MODEL 'off' GIVES, which is why a run with
      ! the new channels off reproduces the previous arithmetic.
      !
      ! CHUNG ET AL.'S NORMALIZATION IS CONSISTENT WITH THAT READING.
      ! Their Eqs. (3)-(5) set sigma(H+) + sigma(H2+) = sigma_i =
      ! Y sigma(abs).  With one count per event sigma(H+) counts N_S + N_D
      ! and sigma(H2+) counts N_M, so their sum is the event count
      ! N_M + N_S + N_D and nothing is counted twice.  The two protons of a
      ! D event enter through h2_channel_stoichiometry below, where the
      ! particle multiplicity belongs; detector response belongs in the
      ! data inversion and particle multiplicity in the reaction source.
      !
      ! ------------------------------------------------------------------
      ! RANGE OF VALIDITY, stated because none of the three inputs covers
      ! the whole grid.
      !
      !  * f_n: measured 33-41 eV only (Chung et al. Table 1, twelve rows,
      !    their yields taken from Backx, Wight & van der Wiel), ramped
      !    linearly to zero at the continuity anchors 32 and 41.5 eV and
      !    zero outside.  Chung et al.: "No yield measurements were made
      !    above 28 eV because it was assumed that the ionization yield
      !    would be unity at higher energies."  Zero outside 32-41.5 eV is
      !    therefore an assumption of the source, not a measurement.
      !  * r, equivalently f_di: measured 18.076-124 eV (Table II).  Above
      !    124 eV it is HELD at its 124 eV value; the source quotes nothing
      !    higher.  See frac_H2_dissociative_ionization in cross_sec.f90,
      !    which carries the same table and the same caveat.
      !  * q_D: NOT measured as a function of energy anywhere this code
      !    has.  It is a MODEL, selected by h2_double_ionization_model,
      !    whose input key `H2 double ionization` DEFAULTS TO 'chung80';
      !    'off' exists to reproduce a pre-E1 result and is not a statement
      !    that the reaction does not happen.  See the model list at
      !    frac_H2_double_of_proton_events.
      !
      ! Uncertainty quoted by the source for the underlying data: "An
      ! estimate of the total uncertainty in the data is about +/-4%-5%".

      use global_parameters, only: e_th_H2, e_th_H2_di, e_th_H2_dd,      &
                                   e_th_HI
      use Cross_sections, only: frac_H2_dissociative_ionization
      use mol_rates, only: h2_dissociation_energy_eV

      implicit none
      private

      public :: h2_channel_cross_sections
      public :: h2_photoabsorption_cross_sections
      public :: h2_channel_stoichiometry
      public :: frac_H2_neutral_dissociation
      public :: frac_H2_double_of_proton_events
      public :: h2_double_ionization_model
      public :: n_h2_channels
      public :: ICH_M, ICH_S, ICH_D, ICH_N
      public :: h2_channel_threshold, h2_channel_name
      public :: h2_channel_energy_recipients
      public :: h2_channel_electron_share, h2_channel_prompt_heat_share
      public :: h2_channel_radiated_share
      public :: h2_double_fragment_kinetic_energy
      public :: e_rad_H2_neutral

      integer, parameter :: n_h2_channels = 4
      integer, parameter :: ICH_M = 1     ! H2+ + e-
      integer, parameter :: ICH_S = 2     ! H + H+ + e-
      integer, parameter :: ICH_D = 3     ! H+ + H+ + 2e-
      integer, parameter :: ICH_N = 4     ! H + H

      ! Threshold of each channel [eV].  The first three are read from the
      ! global constants, so a channel is charged the same threshold here,
      ! in the photoelectron energy h nu - e_th and in the turn-on of the
      ! cross section that feeds it.
      !  e_th_H2    : H2 ionization energy, the sigma_H2 turn-on.
      !  e_th_H2_di : Chung et al. (1993) Table II first row, 18.076 eV.
      !  e_th_H2_dd : Yan et al. (1998) sec. 4, "has a vertical threshold of
      !               51.4 eV".
      !  32.0       : where f_n turns on, the lower continuity anchor of
      !               frac_H2_neutral_dissociation, one eV below the first
      !               measured row (33 eV). The channel has no threshold of
      !               its own, so this is not an ionization threshold and
      !               stays here; it has to be the turn-on of f_n, or the
      !               32-33 eV events would have a cross section and no
      !               energy recipients.
      real*8, parameter :: h2_channel_threshold(n_h2_channels) =          &
           (/ e_th_H2, e_th_H2_di, e_th_H2_dd, 32.000d0 /)

      ! Excitation energy of H(n = 2) above the ground state [eV], the
      ! Bohr value (3/4) I(H) = 10.19883 eV.  It is the internal energy each
      ! fragment of the neutral window carries away (see the recipient table
      ! at h2_channel_energy_recipients) and it is the same number as
      ! e_th_HI - e_th_HI_n2, the Balmer continuum edge of J_inc.f90, which
      ! is written there as e_th_HI/4.
      real*8, parameter :: e_exc_H_n2 = 0.75d0*e_th_HI

      ! Internal excitation carried away by the two fragments of the neutral
      ! window, H2 + hv -> H(2p) + H(2s), and radiated rather than
      ! thermalized.  Source and the reading of it: the (N) row of the
      ! recipient table at h2_channel_energy_recipients.
      real*8, parameter :: e_rad_H2_neutral = 2.0d0*e_exc_H_n2

      character(len=24), parameter :: h2_channel_name(n_h2_channels) =    &
           (/ 'H2 -> H2+ + e           ',                                 &
              'H2 -> H + H+ + e        ',                                 &
              'H2 -> H+ + H+ + 2e      ',                                 &
              'H2 -> H + H             ' /)

      ! Which double-ionization model is in force.  The value below is only
      ! the state before a run is configured: set_energy_vectors copies the
      ! resolved input key `H2 double ionization` here, and that key
      ! DEFAULTS TO 'chung80'.  Selecting 'off' reproduces the pre-E1
      ! arithmetic exactly.
      !   'off'      q_D(E) = 0.  Double ionization stays inside the single
      !              dissociative channel, as it was before E1.
      !   'chung80'  q_D(E) = 0 below 51.4 eV, ramped linearly in E to 0.20
      !              at 80 eV and held at 0.20 above.  Chung et al. (1993),
      !              p. 888, discussing their Fig. 4: "the contribution to
      !              H+ from double ionization is small near the double
      !              ionization threshold.  However, by 80 eV about 20% of
      !              sigma(H+) comes from double ionization.  This
      !              percentage remains fairly constant towards higher
      !              energies."  sigma(H+) is the detector signal of the
      !              header, one count per event, so 20 per cent of it is
      !              the EVENT fraction q_D = N_D/(N_S + N_D) and not a
      !              proton share.  The ramp between 51.4 and 80 eV is NOT
      !              in the source: "small" is all it says there.
      !   'yan_rho'  q_D(E) built from the double-to-single ionization RATIO
      !              rho_D = N_D/(N_M + N_S) of Yan et al. (1998) sec. 4:
      !              "At 110 eV, the measured ratio of double to single
      !              ionization is 0.038. The asymptotic ratio is predicted
      !              to be 0.0225 (Sadeghpour & Dalgarno 1993)."
      !              rho_D is ramped from 0 at 51.4 eV to 0.038 at 110 eV
      !              and relaxed to 0.0225 by 300 eV.  The ramp and the
      !              relaxation are BOTH interpolations between two
      !              published numbers, not measurements.  The conversion
      !              to the event fraction is exact:
      !                q_D = rho_D/(1 + rho_D) * (1 + r)/r
      !              with r = (N_S + N_D)/N_M the detector ratio, which is
      !              q_D = [rho_D/(1 + rho_D)]/f_di in the tabulated
      !              fraction.  It is admissible only where rho_D <= r; an
      !              inadmissible pair of published numbers is a stated
      !              model restriction and stops the run, never a clip,
      !              because clipping would redefine the measured
      !              observable.  The two published numbers are admissible
      !              everywhere the code uses them (MEASURED over
      !              51.4-2000 eV in steps of 0.01 eV: q_D runs from 0 at
      !              the threshold to 0.17463 at 110 eV).
      !
      ! THE TWO MODELS DISAGREE and that is the point of having both: at
      ! 110 eV 'chung80' gives q_D = 0.200 against 'yan_rho' 0.1746, i.e.
      ! rho_D = 0.0438 against 0.0380, 15 per cent more double ionization.
      ! Neither is "the" answer; running both is how the uncertainty is
      ! reported.
      character(len=16) :: h2_double_ionization_model = 'off'

      contains

      !--------------!

      real*8 function frac_H2_neutral_dissociation(E)
      ! sigma_n/sigma(abs), the share of H2 photoabsorptions that leave two
      ! neutral H atoms and no ion.  Chung et al. (1993) Table 1, linearly
      ! interpolated in E, zero outside 32-41.5 eV (the measured 33-41 eV
      ! rows and the two continuity anchors below).  See the range note in
      ! the module header: zero outside the window is the source's
      ! assumption of unit ionization yield, not a measurement.
      real*8, intent(in) :: E
      integer :: k
      integer, parameter :: n_n = 14
      ! The twelve measured rows, plus a zero at each end.  THE TWO ZEROS
      ! ARE CONTINUITY ANCHORS, NOT MEASUREMENTS.  They rest on the source's
      ! own sentence about the yields it used: "This is in general agreement
      ! with the data of Backx, Wight, and van der Wiel except between 32
      ! and 41 eV, where their yields drop from 100% to 92%", i.e. the yield
      ! is unity at both ends of the window.  Without them f_n would step
      ! from 0 to 0.0078 at 33 eV and back at 41 eV, and a discontinuous
      ! cross section is differentiated by the JFNK residual.
      real*8, parameter :: e_n(n_n) = (/                                  &
          32.000d0,                                                       &
          33.000d0, 34.000d0, 35.000d0, 35.500d0, 36.000d0, 36.500d0,     &
          37.000d0, 37.500d0, 38.000d0, 39.000d0, 40.000d0, 41.000d0,     &
          41.500d0 /)
      ! sigma_n/sigma_abs from their columns sigma_n and sigma(abs).
      real*8, parameter :: f_n(n_n) = (/                                  &
          0.0000000d0,                                                    &
          0.0077720d0, 0.0218391d0, 0.0401274d0, 0.0493333d0,             &
          0.0577465d0, 0.0659259d0, 0.0713178d0, 0.0739837d0,             &
          0.0726496d0, 0.0601852d0, 0.0272727d0, 0.0076923d0,             &
          0.0000000d0 /)

      if (E .le. e_n(1) .or. E .ge. e_n(n_n)) then
         frac_H2_neutral_dissociation = 0.0d0
         return
      endif
      k = 1
      do while (k .lt. n_n-1 .and. E .gt. e_n(k+1))
         k = k + 1
      enddo
      frac_H2_neutral_dissociation = f_n(k)                               &
           + (f_n(k+1) - f_n(k))*(E - e_n(k))/(e_n(k+1) - e_n(k))

      end function frac_H2_neutral_dissociation

      !--------------!

      real*8 function frac_H2_double_of_proton_events(E)
      ! q_D(E) = N_D/(N_S + N_D): among the H2 photoionization events that
      ! release a proton, the share in which BOTH nuclei come away as
      ! protons.  It is an EVENT fraction, which is the observable Chung
      ! et al.'s time-of-flight detector delivers (see the module header),
      ! and not a share of the protons.  A MODEL, not a measurement -- see
      ! the model list in the module header.
      real*8, intent(in) :: E
      real*8 :: rho_D, f_di
      real*8, parameter :: E_thr = 51.4d0

      frac_H2_double_of_proton_events = 0.0d0
      if (E .le. E_thr) return

      select case (trim(h2_double_ionization_model))
      case ('off')
         frac_H2_double_of_proton_events = 0.0d0

      case ('chung80')
         if (E .ge. 80.0d0) then
            frac_H2_double_of_proton_events = 0.20d0
         else
            frac_H2_double_of_proton_events =                             &
                 0.20d0*(E - E_thr)/(80.0d0 - E_thr)
         endif

      case ('yan_rho')
         ! rho_D = N_D/(N_M + N_S), the published observable, converted to
         ! the event fraction at this energy's detector ratio
         ! r = (N_S + N_D)/N_M:
         !    q_D = rho_D/(1 + rho_D) * (1 + r)/r
         !        = [rho_D/(1 + rho_D)]/f_di ,   f_di = r/(1 + r).
         if (E .le. 110.0d0) then
            rho_D = 0.038d0*(E - E_thr)/(110.0d0 - E_thr)
         else if (E .ge. 300.0d0) then
            rho_D = 0.0225d0
         else
            rho_D = 0.038d0 + (0.0225d0 - 0.038d0)                        &
                            *(E - 110.0d0)/(300.0d0 - 110.0d0)
         endif
         f_di = frac_H2_dissociative_ionization(E)
         if (f_di .le. 0.0d0) return
         frac_H2_double_of_proton_events = (rho_D/(1.0d0 + rho_D))/f_di

      case default
         frac_H2_double_of_proton_events = 0.0d0
      end select

      ! q_D is a fraction of events, so it lies in [0,1].  It can leave
      ! that range only when the published numbers a model combines are
      ! mutually inadmissible -- for 'yan_rho', where rho_D exceeds the
      ! measured detector ratio r.  Clipping it would silently redefine
      ! the measured observable, so the run stops and states the
      ! restriction instead.
      if (.not. (frac_H2_double_of_proton_events .ge. 0.0d0 .and.         &
                 frac_H2_double_of_proton_events .le. 1.0d0)) then
         write(*,'(a)') ' (h2_photo_channels) ERROR: the double-'//       &
              'ionization event fraction q_D = N_D/(N_S + N_D)'
         write(*,'(a)') '   left [0,1]. The published numbers combined'// &
              ' by model '//trim(h2_double_ionization_model)
         write(*,'(a)') '   are inadmissible at this photon energy.'
         write(*,'(a,es13.5,a,es13.5)') '   E [eV] = ', E,                &
              '   q_D = ', frac_H2_double_of_proton_events
         error stop 1
      endif

      end function frac_H2_double_of_proton_events

      !--------------!

      subroutine h2_channel_cross_sections(E, sigma_tot, sig)
      ! The four mutually exclusive channel cross sections at photon energy
      ! E, given the TOTAL H2 absorption cross section sigma_tot at that
      ! energy (the caller passes sigma_H2(E) so this module does not need
      ! to know the fit).  By construction sum(sig) = sigma_tot to
      ! round-off.
      real*8, intent(in)  :: E, sigma_tot
      real*8, intent(out) :: sig(n_h2_channels)
      real*8 :: f_n, f_di, q_D, s_ion

      sig = 0.0d0
      if (sigma_tot .le. 0.0d0) return

      f_n   = frac_H2_neutral_dissociation(E)
      s_ion = (1.0d0 - f_n)*sigma_tot
      sig(ICH_N) = f_n*sigma_tot

      ! The detector inversion of the module header, written in the
      ! tabulated fraction f_di = r/(1 + r): the H2+ channel takes the
      ! (1 - f_di) share of the ionizing events, and the events that show
      ! up in the H+ signal divide between the single dissociative and the
      ! double channel in the ratio (1 - q_D) : q_D.
      f_di = frac_H2_dissociative_ionization(E)
      if (f_di .le. 0.0d0) then
         sig(ICH_M) = s_ion
         return
      endif
      q_D = frac_H2_double_of_proton_events(E)

      sig(ICH_M) = (1.0d0 - f_di)*s_ion
      sig(ICH_S) = (1.0d0 - q_D)*f_di*s_ion
      sig(ICH_D) = q_D*f_di*s_ion

      end subroutine h2_channel_cross_sections

      !--------------!

      subroutine h2_photoabsorption_cross_sections(E, sigma_tot,          &
                                                   resolve_neutral, sig)
      ! The channel cross sections a run actually uses at photon energy E,
      ! given the TOTAL H2 absorption cross section sigma_tot there: the
      ! selection of h2_channel_cross_sections by the two run switches, the
      ! neutral window (resolve_neutral) and h2_double_ionization_model.
      ! Both the stellar field (set_energy_vectors) and the recombination
      ! photons absorbed on the spot (util_ion_eq) take their channels from
      ! here, so the two fields destroy H2 into the same final states.
      ! sum(sig) = sigma_tot to round-off in every branch.
      real*8,  intent(in)  :: E, sigma_tot
      logical, intent(in)  :: resolve_neutral
      real*8,  intent(out) :: sig(n_h2_channels)
      real*8 :: sig_ion

      if (.not. resolve_neutral .and.                                    &
          trim(h2_double_ionization_model) .eq. 'off') then
         ! Neither of the two channels beyond the single dissociative one
         ! is resolved. The split is then the single branching
         ! frac_H2_dissociative_ionization, the (M)/(S) split written
         ! directly, without the neutral and double shares that are zero.
         sig        = 0.0d0
         sig(ICH_S) = sigma_tot*frac_H2_dissociative_ionization(E)
         sig(ICH_M) = sigma_tot - sig(ICH_S)
         return
      endif

      call h2_channel_cross_sections(E, sigma_tot, sig)
      if (.not. resolve_neutral) then
         ! Fold the neutral share back into the three ionizing channels in
         ! their own proportion, i.e. a unit photoionization yield. Each of
         ! the three is linear in the ionizing part of the cross section,
         ! so rescaling them is the same as evaluating them with a neutral
         ! fraction of zero.
         sig_ion = sig(ICH_M) + sig(ICH_S) + sig(ICH_D)
         if (sig_ion .gt. 0.0d0) then
            sig(ICH_M) = sig(ICH_M)*sigma_tot/sig_ion
            sig(ICH_S) = sig(ICH_S)*sigma_tot/sig_ion
            sig(ICH_D) = sig(ICH_D)*sigma_tot/sig_ion
         endif
         sig(ICH_N) = 0.0d0
      endif

      end subroutine h2_photoabsorption_cross_sections

      !--------------!

      subroutine h2_channel_stoichiometry(ich, dH2, dH2p, dH, dHp, dele)
      ! Change in the count of each species per ONE photon event in channel
      ! ich.  This is the single table every source term is derived from;
      ! nothing downstream should write these numbers out by hand.
      !
      !   channel        H2   H2+   H    H+    e-
      !   (M)            -1   +1    0    0     +1
      !   (S)            -1    0   +1   +1     +1
      !   (D)            -1    0    0   +2     +2
      !   (N)            -1    0   +2    0      0
      !
      ! Both invariants hold row by row, which is what the ledger test
      ! checks: H nuclei  2*dH2 + 2*dH2p + dH + dHp = 0
      !         charge    dH2p + dHp - dele        = 0
      integer, intent(in)  :: ich
      real*8,  intent(out) :: dH2, dH2p, dH, dHp, dele

      dH2 = -1.0d0; dH2p = 0.0d0; dH = 0.0d0; dHp = 0.0d0; dele = 0.0d0
      select case (ich)
      case (ICH_M)
         dH2p = 1.0d0; dele = 1.0d0
      case (ICH_S)
         dH   = 1.0d0; dHp  = 1.0d0; dele = 1.0d0
      case (ICH_D)
         dHp  = 2.0d0; dele = 2.0d0
      case (ICH_N)
         dH   = 2.0d0
      case default
         dH2  = 0.0d0
      end select

      end subroutine h2_channel_stoichiometry

      !--------------!

      ! ==================================================================
      ! THE PHOTON ENERGY LEDGER OF ONE H2 PHOTOEVENT
      !
      ! One absorbed photon of energy E delivers E and nothing else.  The
      ! routine below names every recipient of that energy for each of the
      ! four channels, once each, so that
      !
      !     E = e_reservoir + e_electron + e_fragment + e_radiated
      !
      ! holds to round-off, event by event and channel by channel.  The four
      ! recipients are physically distinct and go to different places in the
      ! code, which is why they are returned separately and not summed here:
      !
      !   e_reservoir  Formation and ionization energy STORED IN THE
      !                PRODUCTS, measured from the reference state in which
      !                every element is a neutral, ground-state, free atom at
      !                rest.  On that
      !                reference eps(H2) = -D0(H2), eps(H) = 0, eps(H+) =
      !                I(H), eps(H2+) = I(H) - D0(H2+); the reservoir charge
      !                of a channel is sum(eps of products) - eps(H2), and
      !                it is the same convention as the species energy
      !                table h(X) of molecular_reaction_heat.f90, whose zero
      !                is H + He + e- at rest.  It is NOT heat.
      !   e_electron   Total kinetic energy of the photoelectron(s).  This is
      !                what the code hands to the degradation partition of
      !                electron_energy_degradation.f90, which divides it into
      !                heat, secondary ionizations and escaping line
      !                radiation.
      !   e_fragment   Kinetic energy of the HEAVY fragments of a
      !                dissociative channel.  Translational from the instant
      !                it is released, so it is prompt heat in full, with no
      !                degradation partition and no quench factor.
      !   e_radiated   Internal excitation left in a product that is not
      !                carried by any evolved population, and therefore
      !                leaves the cell as prompt line radiation rather than
      !                as heat.
      !
      ! ------------------------------------------------------------------
      ! RECIPIENT TABLE, one row per channel, with the source of each entry.
      !
      !  (M) H2 + hv -> H2+ + e-
      !      e_reservoir = e_th_H2 = 15.425927 eV, the adiabatic ionization
      !        energy of H2 (NIST Chemistry WebBook, parameters.f90).  It is
      !        the same number as h(H2+) - h(H2) in the species energy
      !        table of molecular_reaction_heat.f90, whose zero is H + He + e-
      !        at rest, where h(H2+) = IP(H2) - D0(H2) and h(H2) = -D0(H2),
      !        so the two conventions charge one energy.
      !      e_electron  = E - e_th_H2.
      !      e_fragment  = 0: nothing dissociates.
      !      e_radiated  = 0: H2+ is left in its ground electronic state at
      !        the adiabatic threshold, and any vibrational excitation above
      !        it is inside the measured cross section, not resolved here.
      !
      !  (S) H2 + hv -> H + H+ + e-
      !      e_reservoir = e_th_H2_di = 18.076 eV, Chung et al. (1993)
      !        Table II first row.  THIS IS AN ASYMPTOTIC THRESHOLD, not a
      !        vertical one: the thermochemical energy of the same products,
      !        h(H) + h(H+) - h(H2) = D0(H2) + I(H) = 4.478280 + 13.598435 =
      !        18.076715 eV, agrees
      !        with it to 7.1e-4 eV, which is below the precision of either
      !        number.  Chung et al., p. 886: "The threshold energy for
      !        H + H+ production lies at 18.076 eV and is, initially,
      !        accessible only by transitions into the repulsive side of the
      !        stable H2+ 1s sigma_g state."
      !      e_electron  = E - e_th_H2_di.
      !      e_fragment  = 0, AND THIS IS AN APPROXIMATION with a stated
      !        consequence.  Dissociation on the repulsive wall of the
      !        1s sigma_g state does give the H and the H+ a real kinetic
      !        energy release, and no source this code has resolves how the
      !        excess above 18.076 eV divides between that release and the
      !        electron.  Charging all of it to the electron is what the code
      !        does.  It does not lose energy -- both recipients end as heat
      !        plus, for the electron, secondary ionizations -- so the only
      !        error is in the heat/ionization split of the excess, biased
      !        toward too many secondary ionizations.
      !      e_radiated  = 0, for the same want of a resolved source.
      !
      !  (D) H2 + hv -> H+ + H+ + 2e-
      !      e_reservoir = 2 h(H+) - h(H2) = D0(H2) + 2 I(H) = 31.675 eV,
      !        the thermochemical energy of two free protons and two free
      !        electrons measured from H2.  This is the ONLY admissible
      !        reservoir charge: the products are bare nuclei, so their
      !        stored energy is fixed by thermochemistry and cannot be the
      !        51.4 eV threshold.
      !      e_fragment  = e_th_H2_dd - [D0(H2) + 2 I(H)] = 19.725 eV, the
      !        kinetic energy release of the Coulomb explosion.  DERIVED, not
      !        fitted: 51.4 eV is a VERTICAL threshold (Yan, Sadeghpour &
      !        Dalgarno 1998, sec. 4: double photoionization "has a vertical
      !        threshold of 51.4 eV"), so at a photon of exactly 51.4 eV the
      !        two electrons come off with no kinetic energy and the whole
      !        difference between the vertical and the asymptotic threshold
      !        must be carried by the two receding protons.  Energy
      !        conservation at the threshold fixes it with no free parameter.
      !        Physical check against the Franck-Condon picture: two bare
      !        protons at the H2 equilibrium separation R_e = 1.401 a0 repel
      !        with e^2/R_e = 27.2114/1.401 = 19.42 eV, which reproduces the
      !        19.725 eV to 1.5 per cent.  It is taken INDEPENDENT OF E, the
      !        reflection approximation: the separation at which the
      !        transition occurs is set by the ground-state vibrational wave
      !        function, not by the photon, so the excess above the vertical
      !        threshold goes to the electrons.
      !      e_electron  = E - e_th_H2_dd, the total kinetic energy of the
      !        two photoelectrons.
      !      e_radiated  = 0: the products are bare nuclei with no levels.
      !
      !  (N) H2 + hv -> H + H
      !      e_reservoir = -h(H2) = D0(H2), the bond broken.
      !      e_radiated  = 2 E_21 = 20.398 eV, the n = 2 excitation of BOTH
      !        fragments.  The neutral window is not a ground-state
      !        dissociation.  Chung et al. (1993), p. 887, on the same
      !        sigma_n this channel is built from: "fluorescence cross
      !        sections have been measured by Glass-Maujean between 30 and 40
      !        eV by observing the Ly-alpha radiation produced by H(n=2).
      !        They have attributed their observations to the dissociation of
      !        the Q2 1Pi_u(1) state, which dissociates into H(2p) + H(2s)",
      !        and, of the same measurement, "The fluorescence studies of
      !        Glass-Maujean verify that the 1Pi_u(1) state dissociates into
      !        excited neutral atoms and that their curve overlaps sigma_n."
      !        The overlap of the two curves is why the share is taken as
      !        one: the window and the fluorescence are the same reaction.
      !        MAGNITUDE UNCERTAINTY, quoted by the source in the next
      !        sentences: "They obtained a peak value for this state of
      !        approximately 0.05 Mb.  Our value of 0.09 Mb for sigma_n
      !        depends on the accuracy of the ionization yields quoted by
      !        Backx, Wight, and van der Wiel", i.e. read as a magnitude
      !        ratio the excited share would be 0.05/0.09 = 0.56 rather than
      !        1, and Chung et al. attribute the difference to the yields
      !        their sigma_n came from, not to a ground-state branch.
      !        Neither H(2s) nor H(2p) is a product of
      !        h2_channel_stoichiometry, which makes two ground-state H, so
      !        this energy is in no evolved population: it leaves as prompt
      !        Ly-alpha from 2p and as the two-photon continuum from 2s.
      !        E_21 = (3/4) I(H) is the Bohr n = 1 to n = 2 energy, the same
      !        number as e_th_HI - e_th_HI_n2 in the Balmer continuum module.
      !      e_fragment  = E - D0(H2) - 2 E_21, the recoil of the two atoms.
      !        Positive wherever f_n is nonzero (32-41.5 eV), where the
      !        smallest value is 32 - 24.876 = 7.1 eV.
      !      e_electron  = 0: the channel makes no electron.
      ! ==================================================================

      subroutine h2_channel_energy_recipients(ich, E, e_reservoir,        &
                                     e_electron, e_fragment, e_radiated)
      ! The four recipients of one photon of energy E absorbed in channel
      ! ich, in eV.  They sum to E to round-off; see the table above for
      ! what each is and where its value comes from.  Below the channel's
      ! threshold the channel does not exist and every recipient is zero.
      integer, intent(in)  :: ich
      real*8,  intent(in)  :: E
      real*8,  intent(out) :: e_reservoir, e_electron, e_fragment, e_radiated
      real*8 :: D0

      e_reservoir = 0.0d0
      e_electron  = 0.0d0
      e_fragment  = 0.0d0
      e_radiated  = 0.0d0
      if (ich .lt. 1 .or. ich .gt. n_h2_channels) return
      if (E .le. h2_channel_threshold(ich)) return

      D0 = h2_dissociation_energy_eV()

      select case (ich)
      case (ICH_M)
         e_reservoir = e_th_H2
         e_electron  = E - e_th_H2

      case (ICH_S)
         e_reservoir = e_th_H2_di
         e_electron  = E - e_th_H2_di

      case (ICH_D)
         e_reservoir = D0 + 2.0d0*e_th_HI
         e_fragment  = h2_double_fragment_kinetic_energy()
         e_electron  = E - e_th_H2_dd

      case (ICH_N)
         e_reservoir = D0
         e_radiated  = e_rad_H2_neutral
         e_fragment  = E - e_reservoir - e_radiated
      end select

      end subroutine h2_channel_energy_recipients

      !--------------!

      real*8 function h2_double_fragment_kinetic_energy() result(e_ker)
      ! Kinetic energy release [eV] of the two protons of one H2 double
      ! photoionization event: the difference between the VERTICAL threshold
      ! the channel is charged and the ASYMPTOTIC energy of its products.
      ! Derived, with no free parameter, from energy conservation at the
      ! threshold itself, where the two electrons come off at rest; see the
      ! (D) row of the recipient table above for the source of the vertical
      ! threshold and for the Franck-Condon cross-check.  It is independent
      ! of the photon energy (reflection approximation), which is why the
      ! heating integrand can evaluate it once per cell like D0(H2).
      e_ker = e_th_H2_dd - (h2_dissociation_energy_eV() + 2.0d0*e_th_HI)
      end function h2_double_fragment_kinetic_energy

      !--------------!

      real*8 function h2_channel_electron_share(ich, E) result(f)
      ! Share of the absorbed photon energy that leaves as photoelectron
      ! kinetic energy in channel ich.  This is the quantity the heating
      ! integrand multiplies by the degradation heat fraction f_heat, and it
      ! replaces photoelectron_share(e_th_of_this_channel, E) with the same
      ! value, so that one routine states the threshold each channel is
      ! charged.
      integer, intent(in) :: ich
      real*8,  intent(in) :: E
      real*8 :: e_res, e_ele, e_frg, e_rad
      call h2_channel_energy_recipients(ich, E, e_res, e_ele, e_frg, e_rad)
      f = 0.0d0
      if (E .gt. 0.0d0) f = e_ele/E
      end function h2_channel_electron_share

      !--------------!

      real*8 function h2_channel_prompt_heat_share(ich, E) result(f)
      ! Share of the absorbed photon energy that is deposited as heat
      ! IMMEDIATELY, without passing through the electron degradation
      ! partition: the kinetic energy of the heavy fragments.  Zero for the
      ! two channels that leave the nuclei bound or singly charged; the
      ! Coulomb-explosion release for the double channel; the whole recoil
      ! for the neutral window.
      integer, intent(in) :: ich
      real*8,  intent(in) :: E
      real*8 :: e_res, e_ele, e_frg, e_rad
      call h2_channel_energy_recipients(ich, E, e_res, e_ele, e_frg, e_rad)
      f = 0.0d0
      if (E .gt. 0.0d0) f = e_frg/E
      end function h2_channel_prompt_heat_share

      !--------------!

      real*8 function h2_channel_radiated_share(ich, E) result(f)
      ! Share of the absorbed photon energy that leaves the cell as prompt
      ! line radiation: the n = 2 excitation of the two fragments of the
      ! neutral window, zero in every other channel.
      integer, intent(in) :: ich
      real*8,  intent(in) :: E
      real*8 :: e_res, e_ele, e_frg, e_rad
      call h2_channel_energy_recipients(ich, E, e_res, e_ele, e_frg, e_rad)
      f = 0.0d0
      if (E .gt. 0.0d0) f = e_rad/E
      end function h2_channel_radiated_share

      end module h2_photo_channels
