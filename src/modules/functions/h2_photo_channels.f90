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
      ! ------------------------------------------------------------------
      ! HOW THE SPLIT IS BUILT, and what each piece is measured from.
      !
      ! Write the event counts per unit time as N_m, N_s, N_d (channels M,
      ! S, D) and N_n (channel N).  Then
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
      !   r_pm(E) = sigma(H+)/sigma(H2+)
      !
      ! is a ratio of a PROTON count to an H2+ count, so with
      ! P = N_s + 2 N_d the measured content of that column is P/N_m.
      ! Letting d(E) = 2 N_d / P be the share of the protons that comes from
      ! the double channel,
      !
      !   N_m = sigma_ion / [ 1 + r_pm (1 - d/2) ]
      !   P   = r_pm N_m ,   N_s = (1 - d) P ,   N_d = d P / 2 .
      !
      ! With d = 0 this collapses to N_m = sigma_ion/(1 + r_pm) and
      ! N_s = r_pm sigma_ion/(1 + r_pm), i.e. exactly the old
      ! f_di = r_pm/(1 + r_pm) branching.  THAT IS THE DEFAULT, and it is
      ! why turning the new channels off reproduces the previous arithmetic
      ! bit for bit.
      !
      ! NOTE ON CHUNG ET AL.'S OWN NORMALIZATION.  Their Eqs. (3)-(5) set
      ! sigma(H+) + sigma(H2+) = sigma_i = Y sigma(abs), which equates
      ! (N_s + 2 N_d) + N_m to an event count N_m + N_s + N_d and so
      ! over-counts by N_d.  This module does not adopt that identification:
      ! it takes the RATIO column as the measurement and imposes the event
      ! count separately, which is the only way the four channels can be
      ! mutually exclusive and still sum to the absorption.
      !
      ! ------------------------------------------------------------------
      ! RANGE OF VALIDITY, stated because none of the three inputs covers
      ! the whole grid.
      !
      !  * f_n: measured 33-41 eV only (Chung et al. Table 1, twelve rows,
      !    their yields taken from Backx, Wight & van der Wiel).  Set to
      !    zero outside.  Chung et al.: "No yield measurements were made
      !    above 28 eV because it was assumed that the ionization yield
      !    would be unity at higher energies."  Zero outside 33-41 eV is
      !    therefore an assumption of the source, not a measurement.
      !  * r_pm: measured 18.076-124 eV (Table II).  Above 124 eV it is
      !    HELD at its 124 eV value; the source quotes nothing higher.  See
      !    frac_H2_dissociative_ionization in cross_sec.f90, which carries
      !    the same table and the same caveat.
      !  * d: NOT measured as a function of energy anywhere this code has.
      !    It is a MODEL, selected by h2_double_ionization_model, whose
      !    input key `H2 double ionization` DEFAULTS TO 'chung80'; 'off'
      !    exists to reproduce a pre-E1 result and is not a statement that
      !    the reaction does not happen.  See the model list at
      !    frac_H2_double_of_protons.
      !
      ! Uncertainty quoted by the source for the underlying data: "An
      ! estimate of the total uncertainty in the data is about +/-4%-5%".

      use Cross_sections, only: frac_H2_dissociative_ionization

      implicit none
      private

      public :: h2_channel_cross_sections
      public :: h2_channel_stoichiometry
      public :: frac_H2_neutral_dissociation
      public :: frac_H2_double_of_protons
      public :: h2_double_ionization_model
      public :: n_h2_channels
      public :: ICH_M, ICH_S, ICH_D, ICH_N
      public :: h2_channel_threshold, h2_channel_name

      integer, parameter :: n_h2_channels = 4
      integer, parameter :: ICH_M = 1     ! H2+ + e-
      integer, parameter :: ICH_S = 2     ! H + H+ + e-
      integer, parameter :: ICH_D = 3     ! H+ + H+ + 2e-
      integer, parameter :: ICH_N = 4     ! H + H

      ! Threshold of each channel [eV].
      !  15.4   : sigma_H2 threshold, Yan et al. (1998) Eq. 17.
      !  18.076 : Chung et al. (1993) Table II first row.
      !  51.4   : Yan et al. (1998) sec. 4, "has a vertical threshold of
      !           51.4 eV".
      !  33.0   : the low edge of the measured neutral window; the channel
      !           has no threshold of its own, this is where f_n turns on.
      real*8, parameter :: h2_channel_threshold(n_h2_channels) =          &
           (/ 15.400d0, 18.076d0, 51.400d0, 33.000d0 /)

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
      !   'off'      d(E) = 0.  Double ionization stays inside the single
      !              dissociative channel, as it was before E1.
      !   'chung80'  d(E) = 0 below 51.4 eV, ramped linearly in E to 0.20
      !              at 80 eV and held at 0.20 above.  Chung et al. (1993),
      !              discussion of their Fig. 4: "the contribution to H+
      !              from double ionization is small near the double
      !              ionization threshold.  However, by 80 eV about 20% of
      !              sigma(H+) comes from double ionization.  This
      !              percentage remains fairly constant towards higher
      !              energies."  The ramp between 51.4 and 80 eV is NOT in
      !              the source: "small" is all it says there.
      !   'yan_rho'  d(E) built from the double-to-single ionization RATIO
      !              rho = N_d/(N_m + N_s) of Yan et al. (1998) sec. 4:
      !              "At 110 eV, the measured ratio of double to single
      !              ionization is 0.038.  The asymptotic ratio is
      !              predicted to be 0.0225 (Sadeghpour & Dalgarno 1993)."
      !              rho is ramped from 0 at 51.4 eV to 0.038 at 110 eV and
      !              relaxed to 0.0225 by 300 eV.  The ramp and the
      !              relaxation are BOTH interpolations between two
      !              published numbers, not measurements.
      !
      ! THE TWO MODELS DISAGREE and that is the point of having both: at
      ! 110 eV 'chung80' gives rho = 0.023 where 'yan_rho' gives 0.038, a
      ! factor 1.6.  Neither is "the" answer; running both is how the
      ! uncertainty is reported.
      character(len=16) :: h2_double_ionization_model = 'off'

      contains

      !--------------!

      real*8 function frac_H2_neutral_dissociation(E)
      ! sigma_n/sigma(abs), the share of H2 photoabsorptions that leave two
      ! neutral H atoms and no ion.  Chung et al. (1993) Table 1, linearly
      ! interpolated in E, zero outside 33-41 eV.  See the range note in
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

      real*8 function frac_H2_double_of_protons(E)
      ! d(E) = 2 N_d / (N_s + 2 N_d): the share of the protons released by
      ! H2 photoionization that comes from the DOUBLE channel.  A MODEL,
      ! not a measurement -- see the model list in the module header.
      real*8, intent(in) :: E
      real*8 :: rho, rpm, f_di
      real*8, parameter :: E_thr = 51.4d0

      frac_H2_double_of_protons = 0.0d0
      if (E .le. E_thr) return

      select case (trim(h2_double_ionization_model))
      case ('off')
         frac_H2_double_of_protons = 0.0d0

      case ('chung80')
         if (E .ge. 80.0d0) then
            frac_H2_double_of_protons = 0.20d0
         else
            frac_H2_double_of_protons = 0.20d0*(E - E_thr)/(80.0d0 - E_thr)
         endif

      case ('yan_rho')
         ! rho = N_d/(N_m + N_s), the published observable.
         if (E .le. 110.0d0) then
            rho = 0.038d0*(E - E_thr)/(110.0d0 - E_thr)
         else if (E .ge. 300.0d0) then
            rho = 0.0225d0
         else
            rho = 0.038d0 + (0.0225d0 - 0.038d0)                          &
                            *(E - 110.0d0)/(300.0d0 - 110.0d0)
         endif
         ! Invert to d at this energy's measured proton/H2+ ratio.
         ! rho = (d rpm/2)/(1 + rpm(1 - d))  =>  d = 2 rho (1 + rpm)
         !                                          / (rpm (1 + 2 rho))
         f_di = frac_H2_dissociative_ionization(E)
         if (f_di .ge. 1.0d0 .or. f_di .le. 0.0d0) return
         rpm = f_di/(1.0d0 - f_di)
         if (rpm .le. 0.0d0) return
         frac_H2_double_of_protons =                                      &
              2.0d0*rho*(1.0d0 + rpm)/(rpm*(1.0d0 + 2.0d0*rho))

      case default
         frac_H2_double_of_protons = 0.0d0
      end select

      if (frac_H2_double_of_protons .lt. 0.0d0)                           &
          frac_H2_double_of_protons = 0.0d0
      if (frac_H2_double_of_protons .gt. 1.0d0)                           &
          frac_H2_double_of_protons = 1.0d0

      end function frac_H2_double_of_protons

      !--------------!

      subroutine h2_channel_cross_sections(E, sigma_tot, sig)
      ! The four mutually exclusive channel cross sections at photon energy
      ! E, given the TOTAL H2 absorption cross section sigma_tot at that
      ! energy (the caller passes sigma_H2(E) so this module does not need
      ! to know the fit).  By construction sum(sig) = sigma_tot exactly.
      real*8, intent(in)  :: E, sigma_tot
      real*8, intent(out) :: sig(n_h2_channels)
      real*8 :: f_n, f_di, rpm, d, s_ion, Nm, P

      sig = 0.0d0
      if (sigma_tot .le. 0.0d0) return

      f_n   = frac_H2_neutral_dissociation(E)
      s_ion = (1.0d0 - f_n)*sigma_tot
      sig(ICH_N) = f_n*sigma_tot

      f_di = frac_H2_dissociative_ionization(E)
      if (f_di .le. 0.0d0) then
         sig(ICH_M) = s_ion
         return
      endif
      if (f_di .ge. 1.0d0) then
         rpm = 1.0d30
      else
         rpm = f_di/(1.0d0 - f_di)
      endif
      d = frac_H2_double_of_protons(E)

      ! N_m (1 + rpm (1 - d/2)) = s_ion
      Nm = s_ion/(1.0d0 + rpm*(1.0d0 - 0.5d0*d))
      P  = rpm*Nm
      sig(ICH_M) = Nm
      sig(ICH_S) = (1.0d0 - d)*P
      sig(ICH_D) = 0.5d0*d*P

      end subroutine h2_channel_cross_sections

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

      end module h2_photo_channels
