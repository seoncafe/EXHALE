      module h3p_cooling
      ! Optically-thin H3+ infrared cooling (Tier 2 of the lower-atmosphere
      ! plan; docs/lower_atmosphere_coupling.*).
      !
      ! Per-molecule LTE emission E(H3+,T) from Miller, Stallard, Tennyson &
      ! Melin (2013), J. Phys. Chem. A 117, 9770 (Table 5, z(T) fits):
      !
      !   log_e E(T) = sum_n C_n T^n     [E in W molecule^-1 sr^-1],
      !
      ! piecewise over 30-300, 300-800, 800-1800, 1800-5000 K (fit errors
      ! < 0.1% over 300-5000 K; +-5% below 300 K).  The optically-thin
      ! volumetric cooling rate is
      !
      !   Lambda_H3+ = n_H3+ * 4 pi * E(T) * s_nonLTE(T, n_H2)   [W cm^-3]
      !   (multiply by 1e7 for erg s^-1 cm^-3),
      !
      ! where s_nonLTE is the departure factor of Miller et al. (2013)
      ! Table 6 (their Oka & Epp detailed-balance approach), bilinearly
      ! interpolated in (T, log10 n_H2); s -> 1 for n_H2 >~ 1e14 cm^-3
      ! (LTE) and collapses at low density where radiative depopulation
      ! wins.  Table 6 values are upper limits per the paper's caveat on
      ! the proton-hopping rate coefficient.
      !
      ! Above 5000 K the fit is held at E(5000 K): H3+ is thermally
      ! destroyed well below that temperature, so the frozen tail only
      ! guards against transients.  Same NaN-safe bracketing style as the
      ! cooling-table interpolators (code review 2026-07-02).
      !
      ! Reference PDF: references/Miller_2013_JPCA_117_9770.pdf (Tables 5, 6).

      implicit none
      private
      public :: h3p_emission_lte, h3p_nonlte_factor, h3p_cooling_rate

      real*8, parameter :: fourpi = 12.566370614359172d0

      ! --- Table 5 coefficients: log_e E = sum C_n T^n  [W/molecule/sr] ---
      ! 30-300 K (n = 0..9)
      real*8, parameter :: cA(0:9) = [ -81.9599d0,      0.886768d0,      &
           -0.0264611d0,    0.000462693d0, -4.70108d-6,   2.84979d-8,    &
           -1.03090d-10,    2.13794d-13,   -2.26029d-16,  8.66357d-20 ]
      ! 300-800 K (n = 0..6)
      real*8, parameter :: cB(0:6) = [ -92.2048d0,      0.298920d0,      &
           -0.000962580d0,  1.82712d-6,    -2.04420d-9,   1.24970d-12,   &
           -3.22212d-16 ]
      ! 800-1800 K (n = 0..6)
      real*8, parameter :: cC(0:6) = [ -62.7016d0,      0.0526104d0,     &
           -7.22431d-5,     5.93118d-8,    -2.83755d-11,  7.35415d-15,   &
           -8.01994d-19 ]
      ! 1800-5000 K, z(T) variant (n = 0..5)
      real*8, parameter :: cD(0:5) = [ -55.7672d0,      0.0162530d0,     &
           -7.68583d-6,     1.98412d-9,    -2.68044d-13,  1.47026d-17 ]

      ! --- Table 6 non-LTE scaling s(T, n_H2) --------------------------- !
      ! Grid: T rows (K) x log10(n_H2 [cm^-3]) columns 6,8,10,12,14
      ! ([H2] = 1e12..1e20 m^-3 in the paper).  s -> 1 at high density.
      integer, parameter :: nTs = 11, nNs = 5
      real*8, parameter :: sT(nTs) = [ 300.d0, 600.d0, 1000.d0, 1500.d0, &
           2000.d0, 2500.d0, 3000.d0, 3500.d0, 4000.d0, 4500.d0, 5000.d0 ]
      real*8, parameter :: sLogN(nNs) = [ 6.d0, 8.d0, 10.d0, 12.d0, 14.d0 ]
      real*8, parameter :: sTab(nTs,nNs) = reshape( [                    &
      ! log n_H2=6      (T = 300..5000)
        0.0067d0, 0.0011d0, 0.0013d0, 0.0013d0, 0.0011d0, 0.0010d0,      &
        0.0009d0, 0.0008d0, 0.0007d0, 0.0007d0, 0.0007d0,                &
      ! log n_H2=8
        0.3931d0, 0.0313d0, 0.0108d0, 0.0064d0, 0.0049d0, 0.0042d0,      &
        0.0038d0, 0.0036d0, 0.0034d0, 0.0033d0, 0.0033d0,                &
      ! log n_H2=10
        0.9848d0, 0.7449d0, 0.4955d0, 0.3640d0, 0.2894d0, 0.2469d0,      &
        0.2220d0, 0.2064d0, 0.1961d0, 0.1889d0, 0.1836d0,                &
      ! log n_H2=12
        0.9998d0, 0.9965d0, 0.9866d0, 0.9701d0, 0.9499d0, 0.9299d0,      &
        0.9136d0, 0.9014d0, 0.8923d0, 0.8854d0, 0.8802d0,                &
      ! log n_H2=14
        1.0000d0, 1.0000d0, 0.9999d0, 0.9997d0, 0.9995d0, 0.9992d0,      &
        0.9990d0, 0.9988d0, 0.9987d0, 0.9985d0, 0.9985d0 ],              &
        [nTs,nNs] )

      contains

      ! ------------------------------------------------------------------ !

      ! LTE per-molecule emission E(H3+,T) [W molecule^-1 sr^-1].
      double precision function h3p_emission_lte(T) result(E)
      real*8, intent(in) :: T
      real*8 :: tt, lnE
      tt = T
      if (.not. (tt .gt. 30.0d0))  tt = 30.0d0     ! NaN-safe lower clamp
      if (tt .gt. 5000.0d0)        tt = 5000.0d0   ! frozen high-T tail
      if (tt .le. 300.0d0) then
         lnE = poly(cA, 9, tt)
      else if (tt .le. 800.0d0) then
         lnE = poly(cB, 6, tt)
      else if (tt .le. 1800.0d0) then
         lnE = poly(cC, 6, tt)
      else
         lnE = poly(cD, 5, tt)
      endif
      E = exp(lnE)
      end function h3p_emission_lte

      ! Horner evaluation of sum_{n=0..m} c(n) T^n.
      double precision function poly(c, m, T)
      integer, intent(in) :: m
      real*8, intent(in)  :: c(0:m), T
      integer :: n
      poly = c(m)
      do n = m-1, 0, -1
         poly = poly*T + c(n)
      enddo
      end function poly

      ! ------------------------------------------------------------------ !

      ! Non-LTE departure factor s(T, n_H2) from Miller+2013 Table 6,
      ! bilinear in (T, log10 n_H2 [cm^-3]); edges clamped (s=1 at high
      ! density, table edge at low density / low-high T).
      double precision function h3p_nonlte_factor(T, nH2) result(s)
      real*8, intent(in) :: T, nH2
      real*8 :: tt, ln, ft, fn
      integer :: it, in
      tt = T
      if (.not. (tt .gt. sT(1))) tt = sT(1)
      if (tt .gt. sT(nTs))       tt = sT(nTs)
      ln = log10(max(nH2, 1.0d0))
      if (.not. (ln .gt. sLogN(1))) ln = sLogN(1)
      if (ln .ge. sLogN(nNs)) then
         s = 1.0d0                                  ! LTE at high density
         return
      endif
      ! bracket T (non-uniform rows)
      it = 1
      do while (it .lt. nTs-1 .and. sT(it+1) .lt. tt)
         it = it + 1
      enddo
      ft = (tt - sT(it))/(sT(it+1) - sT(it))
      ! bracket log n (uniform step 2)
      in = 1 + int((ln - sLogN(1))/2.0d0)
      if (in .gt. nNs-1) in = nNs-1
      fn = (ln - sLogN(in))/2.0d0
      s  = (1.d0-ft)*(1.d0-fn)*sTab(it,  in  )                            &
         +       ft *(1.d0-fn)*sTab(it+1,in  )                            &
         + (1.d0-ft)*      fn *sTab(it,  in+1)                            &
         +       ft *      fn *sTab(it+1,in+1)
      end function h3p_nonlte_factor

      ! ------------------------------------------------------------------ !

      ! Volumetric H3+ cooling rate [erg s^-1 cm^-3]:
      !   Lambda = n_H3+ * 4 pi * E_LTE(T) * s_nonLTE(T,n_H2) * 1e7.
      double precision function h3p_cooling_rate(T, nH3p, nH2) result(lam)
      real*8, intent(in) :: T, nH3p, nH2
      lam = nH3p*fourpi*h3p_emission_lte(T)*h3p_nonlte_factor(T,nH2)*1.0d7
      end function h3p_cooling_rate

      ! End of module
      end module h3p_cooling
