      program fine_structure_escape_probability
      ! The escape probability of the ground-term fine-structure lines is the
      ! published plane-parallel Doppler-line one, and the optical depth it
      ! is fed is the one that expression takes as its argument.
      !
      ! Production routine exercised:
      ! src/modules/radiation/Cool_coeff.f90 --
      !   line_escape_probability_one_face, which fine_structure_line_transfer
      !   evaluates on the outward and inward line-centre columns of each cell
      !   and hands to the [C I], [C II], [N II] and [O I] statistical
      !   equilibrium as A_ul -> beta A_ul.
      !
      ! PHYSICS AND REFERENCES.
      !  * These are forbidden lines of a static, thermally broadened layer.
      !    The Voigt parameter is of order 1e-8, so a photon escapes out of
      !    the Doppler CORE, not out of a damping wing, and the relevant
      !    published solution is the complete-redistribution Doppler-slab one,
      !    not the damping-wing solution the Ly-alpha closure uses.
      !  * de Jong, Boland & Dalgarno (1980), A&A 91, 68, appendix B define
      !      beta(tau) = (1/2) int dx phi(x) E_2[tau phi(x)]
      !    -- the chance that a photon escapes through the NEAREST boundary of
      !    a plane-parallel layer -- and give, as their eq. (B-7),
      !      beta(tau) = [1 - exp(-2.34 tau)] / (4.68 tau)          tau < 7
      !      beta(tau) = 1 / (4 tau [ln(tau/sqrt(pi))]^(1/2))       tau >= 7
      !    "accurate to within 10% for small and intermediate tau and exact at
      !    very large tau".  This is the expression the production routine
      !    carries, and the first block of assertions is that identity, with
      !    the two constants written out here rather than taken from the
      !    module.  beta(0) = 1/2: half of the photons leave through the near
      !    face, which is why the cell closure adds the two faces.
      !  * phi is a NORMALIZED profile, int phi dx = 1, so phi(0) = 1/sqrt(pi)
      !    and the optical depth at line centre is tau/sqrt(pi).  The argument
      !    of (B-7) is therefore the frequency-integrated depth,
      !    tau_dJ = sqrt(pi) tau_centre.
      !  * The independent check of that convention is Hollenbach & McKee
      !    (1979), ApJS 41, 555, eq. (5.10),
      !      eps(tau) = 1 / (1 + tau [2 pi ln(2.13 + tau^2)]^(1/2)),
      !    the escape probability from the mid-plane of a slab of LINE-CENTRE
      !    half depth tau, counting both faces (their eq. 5.8 integrates over
      !    one hemisphere and is normalized to 1 at tau = 0), "exact at
      !    tau = 1 and with the correct asymptotic behavior at large tau".
      !    Written in one convention the two published expressions are the
      !    same function: 2 beta_dJ(sqrt(pi) tau) and eps_HM(tau) agree to
      !    better than 5% for tau >= 3 and to five digits at tau >= 1e3.
      !    That agreement is asserted here and is what fixes the convention.
      !  * OPEN DEFECT, recorded and not fixed: the production routine is fed
      !    the LINE-CENTRE column (fine_structure_line_transfer integrates
      !    line_center_opacity_lte), so its argument is short of the one
      !    (B-7) defines by sqrt(pi) and the beta it returns is too large --
      !    by a factor rising from 1 in the thin limit to 1.77 asymptotically,
      !    with a maximum of 2.11 at the branch point.  The last block of
      !    assertions measures that departure, so this probe fails the day the
      !    argument is corrected and has to be rewritten with it.  On the
      !    profiles the shipped cases reach the lines stay thin (largest
      !    line-centre depth 0.14, [O I] 63um at the base of the hot Uranus
      !    gate), where the departure is 12%; the bound asserted here is the
      !    one that holds over that whole range.
      !
      ! Tolerances: the identity with (B-7) is an exact relation between the
      ! production routine and an expression written out here and is tested at
      ! 1e-14 relative on the thin branch and 1e-12 on the thick one, where
      ! the routine and the expression group the same factors differently; the branch point is continuous only up to the
      ! exp(-2.34 tau_c) term the thick branch drops, which is 8e-8 relative;
      ! the agreement between the two published forms is a statement about
      ! two independent fits and is tested at their published accuracy.
      use Cooling_Coefficients, only: line_escape_probability_one_face
      use assertion_report
      implicit none

      real*8, parameter :: pi_l  = 4.0d0*atan(1.0d0)
      ! de Jong, Boland & Dalgarno (1980) eq. (B-7), written out here.
      real*8, parameter :: a_dJ  = 2.34d0
      real*8, parameter :: tau_c = sqrt(pi_l)*exp(0.25d0*a_dJ*a_dJ)

      integer, parameter :: nthin = 7, nthick = 7, ncmp = 8
      real*8, parameter :: tau_thin(nthin) =                              &
           (/ 1.0d-4, 1.0d-3, 1.0d-2, 1.0d-1, 5.0d-1, 2.0d0, 6.0d0 /)
      real*8, parameter :: tau_thick(nthick) =                            &
           (/ 7.0d0, 1.0d1, 3.0d1, 1.0d2, 1.0d3, 1.0d5, 1.0d7 /)
      real*8, parameter :: tau_cmp(ncmp) =                                &
           (/ 3.0d0, 7.0d0, 1.0d1, 1.0d2, 1.0d3, 1.0d4, 1.0d5, 1.0d6 /)

      integer :: i
      real*8  :: worst, dev, got, want, ratio
      character(len=64) :: label

      ! ------------------------------------------------------------------ !
      ! (1) The routine is eq. (B-7) of de Jong, Boland & Dalgarno (1980).
      ! ------------------------------------------------------------------ !
      worst = 0.0d0
      do i = 1, nthin
         got  = line_escape_probability_one_face(tau_thin(i))
         want = (1.0d0 - exp(-a_dJ*tau_thin(i)))/(2.0d0*a_dJ*tau_thin(i))
         dev  = abs(got - want)/want
         if (.not. (dev .le. worst)) worst = dev
      enddo
      call check_at_most('fs_escape_is_de_jong_b7_thin_branch', worst,    &
                         1.0d-14)

      worst = 0.0d0
      do i = 1, nthick
         got  = line_escape_probability_one_face(tau_thick(i))
         want = 1.0d0/(4.0d0*tau_thick(i)                                 &
                       *sqrt(log(tau_thick(i)/sqrt(pi_l))))
         dev  = abs(got - want)/want
         if (.not. (dev .le. worst)) worst = dev
      enddo
      call check_at_most('fs_escape_is_de_jong_b7_thick_branch', worst,   &
                         1.0d-12)

      ! Half of the photons leave through the near face of a transparent
      ! layer: beta(0) = 1/2, the value (B-7) takes at the surface.
      call check_relative('fs_escape_transparent_layer_is_one_half',      &
           line_escape_probability_one_face(0.0d0), 0.5d0, 1.0d-14)

      ! The branch point of the production routine is where the two
      ! expressions of (B-7) meet, so the switch loses only the
      ! exp(-2.34 tau_c) term that the thick branch drops.  (B-7) quotes the
      ! meeting point as tau = 7; it is sqrt(pi) exp(2.34^2/4) = 6.9676.
      call check_relative('fs_escape_branch_point',                       &
           tau_c, 6.967558963071725d0, 1.0d-12)
      call check_relative('fs_escape_branch_is_continuous',               &
           line_escape_probability_one_face(tau_c*(1.0d0 + 1.0d-12)),     &
           line_escape_probability_one_face(tau_c*(1.0d0 - 1.0d-12)),     &
           1.0d-6)

      ! ------------------------------------------------------------------ !
      ! (2) (B-7) and Hollenbach & McKee (1979) eq. (5.10) are one function.
      ! ------------------------------------------------------------------ !
      ! Evaluated at the argument (B-7) defines, tau_dJ = sqrt(pi) tau_centre,
      ! twice the one-face value is the both-face escape probability of a slab
      ! of line-centre half depth tau_centre, which is what (5.10) gives.
      worst = 0.0d0
      do i = 1, ncmp
         got  = 2.0d0*line_escape_probability_one_face                    &
                        (sqrt(pi_l)*tau_cmp(i))
         want = hollenbach_mckee_5_10(tau_cmp(i))
         dev  = abs(got - want)/want
         if (.not. (dev .le. worst)) worst = dev
      enddo
      call check_at_most('de_jong_b7_agrees_with_hollenbach_mckee_5_10',  &
                         worst, 5.0d-2)
      ! Asymptotically the two are the same expression, not two fits.
      call check_relative('de_jong_b7_is_hollenbach_mckee_at_large_tau',  &
           2.0d0*line_escape_probability_one_face(sqrt(pi_l)*1.0d5),      &
           hollenbach_mckee_5_10(1.0d5), 1.0d-4)

      ! ------------------------------------------------------------------ !
      ! (3) The departure of the production argument from the published one.
      ! ------------------------------------------------------------------ !
      ! fine_structure_line_transfer passes the LINE-CENTRE column.  These
      ! assertions record how much larger the returned beta is than the
      ! published value at the same physical depth; they are a measurement of
      ! an open defect, not a specification, and go red when it is fixed.
      do i = 1, ncmp
         ratio = line_escape_probability_one_face(tau_cmp(i))             &
               / line_escape_probability_one_face(sqrt(pi_l)*tau_cmp(i))
         write(label,'(a,i0)')                                            &
              'fs_escape_line_centre_argument_departure_case', i
         call check_at_least(trim(label), ratio, 1.75d0)
      enddo
      ! Over the depths the shipped cases reach the lines are thin and the
      ! departure is bounded well below that.
      call check_at_most('fs_escape_departure_at_shipped_depths',         &
           line_escape_probability_one_face(0.14d0)                       &
           /line_escape_probability_one_face(sqrt(pi_l)*0.14d0), 1.13d0)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'fine_structure_escape_probability: ',       &
              assertion_failures, ' assertion(s) failed'
         stop 1
      endif
      write(*,'(a)')                                                      &
           'fine_structure_escape_probability: all assertions passed'

      contains

      !--------------!

      double precision function hollenbach_mckee_5_10(tau) result(eps)
      ! Hollenbach & McKee (1979), ApJS 41, 555, eq. (5.10): the single-flight
      ! escape probability, averaged over the line profile and over angle,
      ! from the mid-plane of a Doppler slab whose LINE-CENTRE optical depth
      ! to either face is tau.
      real*8, intent(in) :: tau
      eps = 1.0d0/(1.0d0 + tau*sqrt(2.0d0*pi_l*log(2.13d0 + tau*tau)))
      end function hollenbach_mckee_5_10

      end program fine_structure_escape_probability
