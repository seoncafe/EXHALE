      program recombination_coefficient_fits
      ! Every recombination coefficient the ionization balance runs on,
      ! against the published fit evaluated here from the published
      ! parameters.
      !
      ! Production routines exercised, all of
      ! src/modules/radiation/Cool_coeff.f90: alpha_rr_metal,
      ! alpha_dr_metal and alpha_rec_metal for the metals,
      ! alpha_rec_FeI_Huang and alpha_rec_FeII_Huang for iron, and
      ! alphaB_HII_new, alphaB_HeII_new, alphaB_HeIII_new and
      ! alpha1_HeII_mao for hydrogen and helium.
      !
      ! WHICH SOURCE EACH STAGE TAKES, and the form it takes.
      !
      !  * Badnell radiative recombination, total to the ground state of
      !    the recombining ion (Badnell 2006, ApJS 167, 334; the fits used
      !    are the author's current tabulation, references/badnell_rr):
      !
      !        alpha_RR = A / [ sqrt(T/T0) (1+sqrt(T/T0))^(1-b)
      !                                    (1+sqrt(T/T1))^(1+b) ],
      !        b        = B + C exp(-T2/T).
      !
      !    Fourteen metal daughters (C I, C II, N I, N II, O I, O II, Mg I,
      !    Mg II, Si I, Si II, Ca II, Na I, K I, S I) and H I, He I, He II.
      !
      !  * Badnell dielectronic recombination (Badnell et al. 2003 and the
      !    series that follows it; references/badnell_dr):
      !
      !        alpha_DR = T^(-3/2) sum_i c_i exp(-E_i/T).
      !
      !    The same fourteen metal daughters, and He I.
      !
      !  * Mao & Kaastra (2016, A&A 587, A84) Eq. (9), the level-resolved
      !    direct-capture coefficient, at n = 1:
      !
      !        R(T) = 1e-10 cm^3/s a0 T^(-b0-c0 ln T) (1+a2 T^-b2)
      !                                              / (1+a1 T^-b1),
      !
      !    T in eV.  Subtracted from the total to make case B for H I,
      !    He I and He II, and used on its own as the He I ground-capture
      !    coefficient.
      !
      !  * Verner's rrfit.f (version 4, 1999) power law, whose Ca rows are
      !    those of Shull & Van Steenberg (1982, ApJS 48, 95):
      !
      !        alpha_RR = rrec(1,Z,N) (T/1e4 K)^(-rrec(2,Z,N)),
      !
      !    N the electron count of the RECOMBINED ion.  Ca I only: the
      !    Badnell tabulation reaches the Mg-like sequence and has no
      !    K-like row, and neither set has a dielectronic fit there.
      !
      !  * Huang et al. (2023, ApJ 951, 123) Eqs. (5) and (6), their fits
      !    to the Nahar iron data, for Fe I and Fe II, where the Badnell
      !    dielectronic project stops short of the Mn-like and Cr-like
      !    sequences.
      !
      ! WHAT IS ASSERTED.  Every row is a transcription, so the assertion
      ! is that the production routine returns the published expression
      ! evaluated with the published constants: 1e-10 relative, round-off
      ! and nothing else.  Three temperatures, 3e3, 1e4 and 3e4 K, which
      ! span the wind: the low one is where the dielectronic term of a
      ! near-neutral stage is still switching on, the high one where it
      ! dominates.
      !
      ! At the time this was written the Ca I row of alpha_rr_metal read
      ! (6.78e-13, 0.80), which is rrec(:,20,19), the rate forming Ca II.
      ! Ca I has twenty electrons and its row is rrec(:,20,20) =
      ! (1.120e-13, 0.9000).  The rate was 6.05 times the published one at
      ! 1e4 K.

      use global_parameters,    only: kb_eV
      use Cooling_Coefficients, only: alpha_rr_metal, alpha_dr_metal,     &
                                      alpha_rec_metal,                    &
                                      alpha_rec_FeI_Huang,                &
                                      alpha_rec_FeII_Huang,               &
                                      alphaB_HII_new, alphaB_HeII_new,    &
                                      alphaB_HeIII_new, alpha1_HeII_mao
      use assertion_report

      implicit none

      integer, parameter :: n_bad = 14, n_dr = 7
      ! Badnell RR rows (A, B, T0, T1, C, T2) and DR rows (c_i, E_i),
      ! keyed by the recombined daughter, in the order of the names below.
      ! Read from references/badnell_rr/clist_K and
      ! references/badnell_dr/clist_K, ground-state entries (M = 1).
      character(len=4), parameter :: dname(n_bad) = &
           [ 'CI  ','CII ','NI  ','NII ','OI  ','OII ','MgI ',            &
             'MgII','SiI ','SiII','CaII','NaI ','KI  ','SI  ' ]
      real*8, parameter :: rr_A(n_bad) = [                                &
           2.995d-9, 2.067d-9, 6.387d-10, 2.410d-9, 6.622d-11, 2.096d-9,  &
           5.452d-11, 1.345d-11, 3.262d-11, 1.964d-10, 2.248d-10,         &
           5.095d-12, 4.528d-11, 1.384d-10 ]
      real*8, parameter :: rr_B(n_bad) = [                                &
           0.7849d0, 0.8012d0, 0.7308d0, 0.7948d0, 0.6109d0, 0.7668d0,    &
           0.6845d0, 0.1074d0, 0.6270d0, 0.6287d0, 0.6605d0,             &
           0.0000d0, 0.4234d0, 0.6886d0 ]
      real*8, parameter :: rr_T0(n_bad) = [                               &
           6.670d-3, 1.643d-1, 9.467d-2, 1.231d-1, 4.136d0, 1.602d-1,     &
           5.637d0, 7.877d2, 1.590d1, 7.712d0, 6.175d0,                   &
           3.546d2, 5.931d0, 1.074d0 ]
      real*8, parameter :: rr_T1(n_bad) = [                               &
           1.943d6, 2.172d6, 2.954d6, 3.016d6, 4.214d6, 4.377d6,          &
           1.551d6, 7.925d7, 4.237d7, 2.951d7, 6.032d6,                   &
           2.310d6, 2.897d9, 7.159d5 ]
      real*8, parameter :: rr_C(n_bad) = [                                &
           0.1597d0, 0.0427d0, 0.2440d0, 0.0774d0, 0.4093d0, 0.1070d0,    &
           0.3945d0, 0.4631d0, 0.2333d0, 0.1523d0, 0.3158d0,              &
           0.9395d0, 0.3049d0, 0.1845d0 ]
      real*8, parameter :: rr_T2(n_bad) = [                               &
           4.955d4, 6.341d4, 6.739d4, 1.016d5, 8.770d4, 1.392d5,          &
           8.360d5, 5.027d5, 5.828d4, 4.804d5, 2.100d5,                   &
           4.297d5, 1.645d5, 1.858d4 ]
      integer, parameter :: dr_n(n_bad) = &
           [ 5, 6, 6, 6, 4, 6, 4, 3, 6, 5, 3, 3, 4, 7 ]
      real*8, parameter :: dr_c(n_dr,n_bad) = reshape( [                  &
       6.346d-9,9.793d-9,1.634d-6,8.369d-4,3.355d-4,0d0,0d0,              &
       3.489d-6,2.222d-7,1.954d-5,4.212d-3,2.037d-4,2.936d-4,0d0,         &
       1.658d-8,2.760d-8,2.391d-9,7.585d-7,3.012d-4,7.132d-4,0d0,         &
       7.712d-8,4.839d-8,2.218d-6,1.536d-3,3.647d-3,4.234d-5,0d0,         &
       5.629d-8,2.550d-7,6.173d-4,1.627d-4,0d0,0d0,0d0,                   &
       1.627d-7,1.262d-7,6.663d-7,3.925d-6,2.406d-3,1.146d-3,0d0,         &
       3.871d-8,4.732d-7,1.599d-3,2.628d-5,0d0,0d0,0d0,                   &
       6.269d-6,9.181d-4,3.082d-4,0d0,0d0,0d0,0d0,                        &
       3.408d-8,1.913d-7,1.679d-7,7.523d-7,8.386d-5,4.083d-3,0d0,         &
       2.930d-6,2.803d-6,9.023d-5,6.909d-3,2.582d-5,0d0,0d0,              &
       3.843d-4,8.040d-3,8.670d-3,0d0,0d0,0d0,0d0,                        &
       2.673d-6,1.918d-4,1.491d-5,0d0,0d0,0d0,0d0,                        &
       6.292d-4,2.350d-3,1.165d-2,2.470d-3,0d0,0d0,0d0,                   &
       7.300d-8,2.577d-7,4.961d-8,9.520d-7,9.586d-7,6.849d-4,6.539d-4 ],  &
       [n_dr,n_bad] )
      real*8, parameter :: dr_e(n_dr,n_bad) = reshape( [                  &
       1.217d1,7.380d1,1.523d4,1.207d5,2.144d5,0d0,0d0,                   &
       2.660d3,3.756d3,2.566d4,1.400d5,1.801d6,4.307d6,0d0,               &
       1.265d1,8.425d1,2.964d2,5.923d3,1.278d5,2.184d5,0d0,               &
       7.113d1,2.765d2,1.439d4,1.347d5,2.496d5,2.204d6,0d0,               &
       5.395d3,1.770d4,1.671d5,2.687d5,0d0,0d0,0d0,                       &
       4.535d1,2.847d2,4.166d3,2.877d4,1.953d5,3.646d5,0d0,               &
       8.415d3,1.682d4,5.000d4,2.759d5,0d0,0d0,0d0,                       &
       4.104d5,5.766d5,7.310d5,0d0,0d0,0d0,0d0,                           &
       2.431d1,1.293d2,4.272d2,3.729d3,5.514d4,1.295d5,0d0,               &
       1.162d2,5.721d3,3.477d4,1.176d5,3.505d6,0d0,0d0,                   &
       2.282d5,3.682d5,4.479d5,0d0,0d0,0d0,0d0,                           &
       3.027d5,3.732d5,4.787d5,0d0,0d0,0d0,0d0,                           &
       2.451d5,3.504d5,4.094d5,4.766d5,0d0,0d0,0d0,                       &
       5.077d2,6.007d2,2.342d3,7.269d3,2.190d4,1.483d5,1.906d5 ],         &
       [n_dr,n_bad] )

      real*8, parameter :: T_probe(3) = [ 3.0d3, 1.0d4, 3.0d4 ]

      integer :: i, k
      real*8  :: T, reference
      character(len=48) :: label

      do i = 1,n_bad
         do k = 1,3
            T = T_probe(k)
            reference = rr_reference(T, rr_A(i), rr_B(i), rr_T0(i),       &
                                     rr_T1(i), rr_C(i), rr_T2(i))
            write(label,'(a,a,a,i0,a)') 'alpha_rr_', trim(dname(i)),      &
                 '_at_', nint(T), 'K'
            call check_relative(trim(label), alpha_rr_metal(dname(i),T),  &
                 reference, 1.0d-10)

            reference = dr_reference(T, dr_n(i), dr_c(:,i), dr_e(:,i))
            write(label,'(a,a,a,i0,a)') 'alpha_dr_', trim(dname(i)),      &
                 '_at_', nint(T), 'K'
            call check_relative(trim(label), alpha_dr_metal(dname(i),T),  &
                 reference, 1.0d-10)

            reference = rr_reference(T, rr_A(i), rr_B(i), rr_T0(i),       &
                                     rr_T1(i), rr_C(i), rr_T2(i))         &
                      + dr_reference(T, dr_n(i), dr_c(:,i), dr_e(:,i))
            write(label,'(a,a,a,i0,a)') 'alpha_rec_', trim(dname(i)),     &
                 '_at_', nint(T), 'K'
            call check_relative(trim(label), alpha_rec_metal(dname(i),T), &
                 reference, 1.0d-10)
         enddo
      enddo

      ! Ca I: Verner rrfit.f rrec(:,20,20), the row of the twenty-electron
      ! (neutral) calcium the recombination forms.  No dielectronic term.
      do k = 1,3
         T = T_probe(k)
         reference = 1.120d-13*(T/1.0d4)**(-0.900d0)
         write(label,'(a,i0,a)') 'alpha_rr_CaI_at_', nint(T), 'K'
         call check_relative(trim(label), alpha_rr_metal('CaI',T),        &
              reference, 1.0d-10)
         write(label,'(a,i0,a)') 'alpha_dr_CaI_at_', nint(T), 'K'
         call check_absolute(trim(label), alpha_dr_metal('CaI',T),        &
              0.0d0, 0.0d0)
      enddo

      ! Iron, Huang et al. (2023) Eqs. (5) and (6).
      do k = 1,3
         T = T_probe(k)
         reference = 2.833d-8*T**(-1.5d0)*exp(-5.731d4/T)                 &
                     *(1.0d0 + 1.383d4*exp(-120.4d0/T))                   &
                   + 1.248d-12*(T/1.0d4)**(-0.485d0)
         write(label,'(a,i0,a)') 'alpha_rec_FeI_at_', nint(T), 'K'
         call check_relative(trim(label), alpha_rec_FeI_Huang(T),         &
              reference, 1.0d-10)
         reference = 1.094d-5*T**(-1.5d0)*exp(-1.490d4/T)                 &
                     *(1.0d0 + 36.74d0*exp(-1.153d5/T))                   &
                   + 1.728d-12*(T/1.0d4)**(-0.618d0)
         write(label,'(a,i0,a)') 'alpha_rec_FeII_at_', nint(T), 'K'
         call check_relative(trim(label), alpha_rec_FeII_Huang(T),        &
              reference, 1.0d-10)
      enddo

      ! Hydrogen and helium: Badnell total (with the He I dielectronic
      ! series) minus the Mao and Kaastra n = 1 direct capture.
      do k = 1,3
         T = T_probe(k)
         reference = rr_reference(T, 8.318d-11, 0.7472d0, 2.965d0,        &
                                  7.001d5, 0.0d0, 0.0d0)                  &
                   - mao_reference(T, 2.3390d-2, 1.2310d0, 1.0380d-2,     &
                                   1.6430d1, 7.5080d-1, 8.9970d-2,        &
                                   3.5560d-1)
         write(label,'(a,i0,a)') 'alphaB_HII_at_', nint(T), 'K'
         call check_relative(trim(label), alphaB_HII_new(T), reference,   &
              1.0d-10)

         reference = rr_reference(T, 5.235d-11, 0.6988d0, 7.301d0,        &
                                  4.475d6, 0.0829d0, 1.682d5)             &
                   + T**(-1.5d0)*( 1.417d-3*exp(-4.633d5/T)               &
                                 + 2.235d-4*exp(-5.532d5/T)               &
                                 - 2.185d-5*exp(-8.887d5/T) )             &
                   - mao_reference(T, 3.7790d-2, 1.1400d0, 3.9760d-3,     &
                                   1.2030d2, 9.9590d-1, 3.6800d0,         &
                                   4.1030d-1)
         write(label,'(a,i0,a)') 'alphaB_HeII_at_', nint(T), 'K'
         call check_relative(trim(label), alphaB_HeII_new(T), reference,  &
              1.0d-10)

         reference = mao_reference(T, 3.7790d-2, 1.1400d0, 3.9760d-3,     &
                                   1.2030d2, 9.9590d-1, 3.6800d0,         &
                                   4.1030d-1)
         write(label,'(a,i0,a)') 'alpha1_HeII_at_', nint(T), 'K'
         call check_relative(trim(label), alpha1_HeII_mao(T), reference,  &
              1.0d-10)

         reference = rr_reference(T, 1.818d-10, 0.7492d0, 1.017d1,        &
                                  2.786d6, 0.0d0, 0.0d0)                  &
                   - mao_reference(T, 1.6140d-1, 1.1180d0, 1.4260d-2,     &
                                   3.7370d1, 7.3750d-1, 4.4050d-1,        &
                                   3.6220d-1)
         write(label,'(a,i0,a)') 'alphaB_HeIII_at_', nint(T), 'K'
         call check_relative(trim(label), alphaB_HeIII_new(T), reference, &
              1.0d-10)
      enddo

      ! A daughter the tabulation does not carry is not a rate.
      call check_absolute('alpha_rec_unknown_daughter',                   &
           alpha_rec_metal('CaIII',1.0d4), 0.0d0, 0.0d0)

      if (assertion_failures .gt. 0) stop 1

      contains

      !--------------!

      double precision function rr_reference(T,A,B,T0,T1,C,T2)
      ! The Badnell radiative recombination form, written here from the
      ! key of his tabulation and not from the production module.
      real*8, intent(in) :: T,A,B,T0,T1,C,T2
      real*8 :: b_eff
      b_eff = B + C*exp(-T2/T)
      rr_reference = A / ( sqrt(T/T0)*(1.0d0+sqrt(T/T0))**(1.0d0-b_eff)   &
                                     *(1.0d0+sqrt(T/T1))**(1.0d0+b_eff) )
      end function rr_reference

      !--------------!

      double precision function dr_reference(T,n,c,e)
      ! The Badnell dielectronic recombination series.
      real*8,  intent(in) :: T, c(:), e(:)
      integer, intent(in) :: n
      integer :: m
      dr_reference = 0.0d0
      do m = 1,n
         dr_reference = dr_reference + c(m)*exp(-e(m)/T)
      enddo
      dr_reference = T**(-1.5d0)*dr_reference
      end function dr_reference

      !--------------!

      double precision function mao_reference(T,a0,b0,c0,a1,b1,a2,b2)
      ! Mao and Kaastra (2016) Eq. (9), with the electron temperature in
      ! eV as the paper states.
      real*8, intent(in) :: T,a0,b0,c0,a1,b1,a2,b2
      real*8 :: te
      te = T*kb_eV
      mao_reference = 1.0d-10*a0*te**(-b0-c0*log(te))                     &
                      *(1.0d0 + a2*te**(-b2))/(1.0d0 + a1*te**(-b1))
      end function mao_reference

      end program recombination_coefficient_fits
