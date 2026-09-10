      program voronov_carbon_ionization
      ! The electron-impact ionization rate coefficients of neutral and
      ! singly ionized carbon, against the published Voronov (1997) rows.
      !
      ! Origin: the C I finding of
      ! docs/audit_20260905/source_audit/group_B_report.md.  Production
      ! routines exercised: ion_coeff_CI_range and ion_coeff_CII_range of
      ! src/modules/radiation/Cool_coeff.f90.
      !
      ! REFERENCE.
      !   Voronov (1997, ADNDT 65, 1) fits every electron-impact
      !   ionization rate coefficient of an ion by
      !
      !       k(T) = A (1 + P sqrt(U)) U^K exp(-U) / (X + U),
      !       U    = dE / (k_B T),
      !
      !   with one row of (dE, P, A, X, K) per ion in his Table 1.  The
      !   row of neutral carbon (C 0) is (11.3 eV, 0, 6.85e-8, 0.193,
      !   0.25), so P is zero and the bracket is absent:
      !
      !       k_CI(T) = 6.85e-8 U^0.25 exp(-U) / (0.193 + U).
      !
      !   The row of C+ is (24.4 eV, 1, 1.86e-8, 0.286, 0.24) and does
      !   carry the bracket.  The fit is stated for 1e3 K to 1e9 K, which
      !   contains the three temperatures used here.
      !
      !   dE is the ionization potential.  This file writes the potentials
      !   to four figures throughout (C I 11.26 eV, C II 24.38 eV, N I
      !   14.53 eV, O I 13.62 eV) where Voronov's table prints three, so
      !   the reference below uses 11.26 and 24.38 as the code does; the
      !   note lines print what Voronov's own three-figure value would give
      !   instead, since exp(-U) turns a 0.35% change of dE into a percent
      !   level change of the rate.
      !
      !   Second witness for the C I row, read not run:
      !   p-winds/p_winds/carbon.py:160-162 (Dos Santos et al.),
      !     ionization_rate_ci = 6.85E-8 * (0.193 + energy_ratio_ci) ** (-1) \
      !         * energy_ratio_ci ** 0.25 * np.exp(-energy_ratio_ci)
      !   with energy_ratio_ci = 11.3 / (k_B T): the same P = 0 form.
      !
      ! WHAT IS ASSERTED.
      !   The two rows are transcriptions, not fits: the assertion is that
      !   the routine returns the published expression evaluated with the
      !   published constants, so the tolerance is round-off, 1e-12
      !   relative.  The reference constants are written with the same
      !   literal kinds the production rows use, so that what is compared
      !   is the algebraic form and the constants, not the rounding of a
      !   default-real literal.
      !
      !   C II is asserted alongside C I because both rows share the one
      !   functional form: it is the control that shows the form itself is
      !   right and that only the C I constants were in question.
      !
      ! At 35d9dd5 the C I row read (P, X, K) = (0.193, 0.25, 0.25): X's
      ! published value had been taken as P and K's as X, adding a
      ! (1 + 0.193 sqrt(U)) factor that does not belong to neutral carbon
      ! and moving the pole of the denominator.  The rate came out 2.56x
      ! the published one at 2000 K, 1.69x at 1e4 K and 1.29x at 5e4 K.

      use global_parameters,    only: N, Ng, kb_eV
      use Cooling_Coefficients, only: ion_coeff_CI_range, ion_coeff_CII_range
      use assertion_report

      implicit none

      ! Voronov (1997) Table 1, row C 0, with the four-figure ionization
      ! potential this file uses.  Default-real literals, as the row has.
      real*8, parameter :: dE_CI = 11.26, A_CI = 6.85e-8, X_CI = 0.193,   &
                           K_CI = 0.25
      ! Row C+ (P = 1).
      real*8, parameter :: dE_CII = 24.38, A_CII = 1.86e-8,               &
                           X_CII = 0.286,  K_CII = 0.24
      ! Voronov's own three-figure potentials, for the note lines only.
      real*8, parameter :: dE_CI_table = 11.3d0

      real*8, dimension(:), allocatable :: T_K, rate_CI, rate_CII
      real*8 :: u, reference, table_dE_ratio
      character(len=56) :: label
      integer :: j

      N = 3
      allocate(T_K(1-Ng:N+Ng), rate_CI(1-Ng:N+Ng), rate_CII(1-Ng:N+Ng))

      ! Three temperatures spanning the range in which carbon ionizes in
      ! an escaping atmosphere: the molecular base, the hydrogen-ionized
      ! wind, and the hot outer flow where the collisional term overtakes
      ! recombination.
      T_K    = 1.0d4
      T_K(1) = 2.0d3
      T_K(2) = 1.0d4
      T_K(3) = 5.0d4

      call ion_coeff_CI_range (T_K, rate_CI,  1, N)
      call ion_coeff_CII_range(T_K, rate_CII, 1, N)

      do j = 1, N
         u = dE_CI/(kb_eV*T_K(j))
         reference = A_CI*u**K_CI*exp(-u)/(X_CI + u)
         write(label,'(a,i0,a)') 'c_i_electron_impact_ionization_T',       &
              int(T_K(j)), 'K'
         call check_relative(label, rate_CI(j), reference, 1.0d-12)

         ! Diagnostic, not an assertion: what Voronov's three-figure
         ! potential 11.3 eV would give relative to the 11.26 eV used.
         table_dE_ratio = (dE_CI_table/(kb_eV*T_K(j)))**K_CI               &
              *exp(-dE_CI_table/(kb_eV*T_K(j)))                            &
              /(X_CI + dE_CI_table/(kb_eV*T_K(j)))                         &
              /(u**K_CI*exp(-u)/(X_CI + u))
         write(*,'(a,i0,a,f8.5)')                                          &
              'note: c_i_rate_ratio_with_voronov_table_dE_11.3eV_T',       &
              int(T_K(j)), 'K = ', table_dE_ratio
      enddo

      do j = 1, N
         u = dE_CII/(kb_eV*T_K(j))
         reference = A_CII*(1.0 + sqrt(u))*u**K_CII*exp(-u)/(X_CII + u)
         write(label,'(a,i0,a)') 'c_ii_electron_impact_ionization_T',      &
              int(T_K(j)), 'K'
         call check_relative(label, rate_CII(j), reference, 1.0d-12)
      enddo

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'voronov_carbon_ionization: ',                &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'voronov_carbon_ionization: all assertions passed'

      end program voronov_carbon_ionization
