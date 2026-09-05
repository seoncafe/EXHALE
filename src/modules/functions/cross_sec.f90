      module Cross_sections
	   ! Energy-dependent photoionization cross sections (1e-18 cm^2)
	
      use global_parameters
	
      implicit none
	
	   contains
	
	   !----- Hydrogenic atoms -----! 
	
      double precision function sigma(E,Z)
      real*8, intent(in) :: Z,E
      real*8 :: E_0,eps,cut

      ! Ionization treshold
      E_0 = 13.6*Z*Z

      if (E.gt.E_0) then      ! For E > E_th
            
         ! Substitution
         eps = sqrt(E/E_0-1.0)               
         
         ! Cross section value
         sigma = 6.3/(Z*Z)*(E_0/E)**4.0        &
                  *exp(4.0-4.0*atan(eps)/eps)    &
                  /(1.0-exp(-2.0*pi/eps))
            
      else
         ! Cross section value
         sigma = 6.3/(Z*Z) 
      
      endif

      ! Correct if E = E_th
      cut = 0.99999*E_0

      if(E.lt.cut) sigma = 0.0

      ! End of function
      end function
      
      !----------------------------------------
      
      !----- HeI -----! 
      
      ! He I GROUND state (1^1S, singlet) photoionization.  Default: Verner+1996
      ! (VFKY96) single-shell fit.  ates_photoion_rate=.true. reverts to the legacy
      ! ATES two-term fit.  (The 2^3S metastable TRIPLET is handled separately by
      ! sigma_HeI23S; this routine is the singlet ground state only.)
      double precision function sigma_HeI(E)
      real*8,intent(in) :: E
      real*8 :: eth

      if (ates_photoion_rate) then
         ! --- Legacy ATES two-term fit ---
         eth = 24.6*0.999
         if (E.ge.eth) then
            sigma_HeI = 0.6935/((E*1.0d-2)**1.82+(E*1.0d-2)**3.23)
         else
            sigma_HeI = 0.0
         endif
      else
         ! --- Default: Verner, Ferland, Korista & Yakovlev 1996 He I 1^1S ---
         ! [E_th,E_0,sigma_0,y_a,P,y_w,y_0,y_1] (sigma_0 in 1e-18 cm^2)
         sigma_HeI = sigma_VFKY96(E, 24.59d0, 1.361d1, 9.492d2, 1.469d0, &
                                     3.188d0, 2.039d0, 4.434d-1, 2.136d0)
      endif

      ! End of function
      end function sigma_HeI
      
      !----------------------------------------
      
      ! HeI(23S) triplet photoionization cross section.
      ! Two Verner+1996 (VFKY96) single-shell forms joined by a power-law
      ! (log-linear) bridge, used piecewise:
      !   wing A = threshold -> Cooper minimum   (E <= x3 ~ 34.7 eV),
      !   bridge = log-linear rise to the bump   (x3 < E < x4 ~ 45.6 eV),
      !   wing B = bump -> high-E tail           (E >= x4).
      ! Each wing is a ~4% fit to the former Norcross-1971 broken PL / TOPbase
      ! background; the bridge joins them continuously.  Threshold (E<4.78 eV)
      ! and the high-E tail are handled by wing A/B themselves (no hard cutoff).
      ! See docs/photoion_cross_sections.tex (sec. "He 2^3S as two VFKY96
      ! forms") and docs/Update_EXHALE.  The former Norcross broken-PL +
      ! TOPbase-tail implementation is retained, commented out, below.
      !	Fit of data from Norcross (1971)
      ! Note: Oklopcic calculations include HeITR only in the
      !	range [4.8-13.6] eV
	
      double precision function sigma_HeI23S(E)
      real*8, intent(in) :: E
      real*8 :: logE, x3, x4, y3b, y4b, mb
      ! VFKY96 fit parameters [E_th,E_0,sigma_0,y_a,P,y_w,y_0,y_1] (Mb, eV):
      !   wing A -- former segments 1-2 (threshold -> Cooper min, ~4.8-34.7 eV)
      real*8, parameter :: EthA=4.78d0,  E0A=2.645d0, s0A=2.08d1,  yaA=1.0d12,   &
                           PA=3.42d0,    ywA=2.681d0, y0A=1.956d0, y1A=2.603d0
      !   wing B -- former segments 4-5 (bump -> tail, ~45.6-320 eV and beyond)
      real*8, parameter :: EthB=45.59d0, E0B=49.68d0, s0B=1.052d3, yaB=4.393d-2, &
                           PB=2.941d0,   ywB=1.717d0, y0B=5.488d-5,y1B=1.118d0

      ! Transition (bridge) node energies: Cooper minimum x3 and bump x4
      x3 = log10(hp_eV*c_light/(357.340*1e-8))   ! ~34.70 eV
      x4 = log10(hp_eV*c_light/(271.940*1e-8))    ! ~45.59 eV

      logE = log10(E)

      if (logE .le. x3) then
         ! Wing A: threshold -> Cooper minimum (VFKY96; =0 below E=4.78 eV)
         sigma_HeI23S = sigma_VFKY96(E, EthA,E0A,s0A,yaA,PA,ywA,y0A,y1A)
      else if (logE .ge. x4) then
         ! Wing B: resonance-averaged bump -> high-E tail (VFKY96)
         sigma_HeI23S = sigma_VFKY96(E, EthB,E0B,s0B,yaB,PB,ywB,y0B,y1B)
      else
         ! Bridge: log-linear (power-law) join of the two wings
         y3b = log10(sigma_VFKY96(10.0d0**x3, EthA,E0A,s0A,yaA,PA,ywA,y0A,y1A))
         y4b = log10(sigma_VFKY96(10.0d0**x4, EthB,E0B,s0B,yaB,PB,ywB,y0B,y1B))
         mb  = (y4b - y3b)/(x4 - x3)
         sigma_HeI23S = 10.0d0**(mb*(logE - x3) + y3b)
      endif

      ! ================================================================
      ! FORMER implementation (Norcross 1971 broken power law + re-aimed
      ! last segment + TOPbase tail), retained for reference:
      ! ----------------------------------------------------------------
      !   real*8 :: logsigma, a1,c1,a2,c2,a3,c3, x1,x2,xJ, y3,y4,m,
      !  &          a3n,c3n,ptail,ctail,tailJ
      !   x1 = log10(hp_eV*c_light/(2593.01*1e-8))*0.9999   ! ~4.78 eV
      !   x2 = log10(hp_eV*c_light/(1655.63*1e-8))          ! ~7.49 eV
      !   xJ = log10(70.0d0)                                ! 70 eV junction
      !   a1=-0.8134; a2=-1.772; a3=-3.039; c1=1.240
      !   c2 = c1 + x2*(a1-a2); c3 = 5.470                  ! continuity
      !   y3 = a2*x3 + c2; y4 = a3*x4 + c3; m = (y4-y3)/(x4-x3)
      !   ptail=-2.9964; ctail=5.4623; tailJ=ptail*xJ+ctail
      !   a3n=(tailJ-y4)/(xJ-x4); c3n=y4-a3n*x4
      !   if (logE .lt. x1) then; sigma_HeI23S=0.0; return; endif
      !   if (logE .le. x2)                    logsigma = a1*logE + c1
      !   if (logE .gt. x2 .and. logE .le. x3) logsigma = a2*logE + c2
      !   if (logE .gt. x3 .and. logE .lt. x4) logsigma = m*(logE-x3) + y3
      !   if (logE .ge. x4 .and. logE .le. xJ) logsigma = a3n*logE + c3n
      !   if (logE .gt. xJ)                    logsigma = ptail*logE + ctail
      !   sigma_HeI23S = 10.0**logsigma
      !   [pre-extension variant: hard cutoff at x5=log10(59.2 eV),
      !    last segment slope a3=-3.039]
      ! ================================================================

      ! End of function
      end function sigma_HeI23S
      
      !----------------------------------------

      !----- Metal photoionization (Verner et al. 1996) -----!
      !
      ! Single-shell analytic fit from Verner, Ferland, Korista &
      ! Yakovlev 1996, ApJ 465, 487, Eqs. (1)-(4). Returns sigma in
      ! units of 1e-18 cm^2 (consistent with rest of this module).
      !
      ! Two forms are provided:
      !
      !  sigma_VFKY96 -- the FULL published Eq. (1) with the (y_0,y_1)
      !  offset/asymptote parameters (photo.dat Table 1, phfit2.f ph2
      !  branch):
      !       x     = E/E_0 - y_0
      !       z     = sqrt(x^2 + y_1^2)
      !       sigma = sigma_0 * ((x-1)^2 + y_w^2) * z^(-Q)
      !                       * (1 + sqrt(z/y_a))^(-P)
      !       Q     = 5.5 - 0.5*P                     [no orbital-l term]
      !  photo.dat / ph2 columns: E_0, sigma_0, y_a, P, y_w, y_0, y_1.
      !
      !  sigma_Verner96 -- the SIMPLER Verner & Yakovlev 1995 (A&AS 109,
      !  125) inner-shell form, retained for ions whose published fit has
      !  y_0=y_1=0 (Mg I/II, Na I) and for the non-OP K I (phfit2 ph1):
      !       y     = E/E_0
      !       Q     = 5.5 + l - 0.5*P
      !       sigma = sigma_0 * ((y-1)^2 + y_w^2) * y^(-Q)
      !                       * (1 + sqrt(y/y_a))^(-P)
      !  The two agree exactly when y_0 = y_1 = 0 and l = 0.

      double precision function sigma_VFKY96(E,E_th,E_0,s_0,y_a,P,y_w, &
                                             y_0,y_1)
      real*8, intent(in) :: E,E_th,E_0,s_0,y_a,P,y_w,y_0,y_1
      real*8 :: x,z,Q,Fy

      if (E .lt. E_th) then
         sigma_VFKY96 = 0.0
         return
      endif

      x  = E/E_0 - y_0
      z  = sqrt(x*x + y_1*y_1)
      Q  = 5.5d0 - 0.5d0*P
      Fy = ((x-1.0d0)**2 + y_w**2) * z**(-Q)          &
                                   * (1.0d0 + sqrt(z/y_a))**(-P)
      sigma_VFKY96 = s_0 * Fy

      end function sigma_VFKY96

      !--------------!

      double precision function sigma_Verner96(E,E_th,E_0,s_0,y_a,P,y_w,l)
      real*8, intent(in) :: E,E_th,E_0,s_0,y_a,P,y_w
      integer, intent(in) :: l
      real*8 :: y,Q,Fy

      if (E .lt. E_th) then
         sigma_Verner96 = 0.0
         return
      endif

      y  = E/E_0
      Q  = 5.5 + dble(l) - 0.5*P
      Fy = ((y-1.0)**2 + y_w**2) * y**(-Q)            &
                                 * (1.0 + sqrt(y/y_a))**(-P)
      sigma_Verner96 = s_0 * Fy

      end function sigma_Verner96

      !--------------!

      ! H2 TOTAL PHOTOABSORPTION cross section.
      !
      ! WHAT THIS REPLACED, AND WHY (E1b).  It used to be the three-piece
      ! analytic fit of Yan, Sadeghpour & Dalgarno (1998), ApJ 496, 1044,
      ! Eqs. 17-19 throughout.  Those three pieces DO NOT MEET at the two
      ! energies where the paper joins them; evaluated on either side, in
      ! units of 1e-18 cm^2,
      !
      !     18 eV   Eq. 17 gives 5.693, Eq. 18 gives 9.787   (x 1.719)
      !     85 eV   Eq. 18 gives 0.1281, Eq. 19 gives 0.0793 (x 0.619)
      !
      ! Both formulas were re-derived by hand from the published
      ! coefficients, so this is a property of the fit as published and not
      ! a transcription error.  A 72 per cent step in an opacity is a
      ! physical defect: the JFNK residual differentiates it, and the H2
      ! destruction rate jumps across it.
      !
      ! WHAT IS USED NOW.  Over the range it was measured, the MEASUREMENT
      ! is used, taken from the paper that made it:
      !
      !   Samson, J. A. R., & Haddad, G. N. 1994, J. Opt. Soc. Am. B 11,
      !   277, "Total photoabsorption cross sections of H2 from 18 to
      !   113 eV", their TABLE 1, "Recommended Total Photoabsorption Cross
      !   Sections of H2 from 18 eV to 300 eV", 72 rows.
      !
      ! Their stated accuracy is +/-2% to +/-3% over 18-113 eV, the measured
      ! range, and +/-3% to +/-4% from 113 to 300 eV.  Yan Eq. 18 is not
      ! used at all: the table covers its whole range.
      !
      !   15.4 <= E < 18 eV : k_bx * Backx et al. (1976) Table 1
      !   18 <= E <= 300 eV : Samson & Haddad (1994) Table 1
      !   E > 300 eV        : k_hi * Yan Eq. 19     (E^-7/2 sum-rule tail)
      !
      ! all tables interpolated in log sigma against log E, and
      !   k_bx = 9.85/Backx(18)     = 1.016415
      !   k_hi = 0.00154/Eq19(300)  = 0.885805
      !
      ! NO PIECE OF THE YAN FIT SURVIVES BELOW 300 eV.
      !
      ! WHY THE TABLE IS THIS ONE AND NOT CHUNG'S.  Chung, Lee, Masuoka &
      ! Samson (1993) Table II carries a sigma(abs) column over
      ! 18.076-124 eV, and this code used to take it from there.  It is the
      ! SAME data: their text says "Values for a(abs) were taken from the
      ! data of Samson and Haddad" and their reference 31 is
      ! "J. A. R. Samson and G. N. Haddad, J. Opt. Soc. Am. (submitted)",
      ! the paper above.  Checked row by row at the 50 energies the two
      ! tables share, 49 agree exactly and one differs in its last digit
      ! (62 eV, 0.237 against 0.236, 0.42 per cent).  The original is used
      ! because it reaches 300 eV where Chung's stops at 124, so the
      ! high-energy join moves out of the band the wind actually absorbs
      ! in.  Chung's Table II is still the source of the H+/H2+ BRANCHING
      ! (frac_H2_dissociative_ionization) and of the neutral fraction
      ! (h2_photo_channels), both of which are ratios and are unaffected.
      !
      ! 15.4-18 eV: A TABLE, from the experiment the other two normalize to.
      !
      !   Backx, C., Wight, G. R., & Van der Wiel, M. J. 1976, J. Phys. B 9,
      !   315, "Oscillator strengths (10-70 eV) for absorption, ionization
      !   and dissociation in H2, HD and D2, obtained by an electron-ion
      !   coincidence method", their TABLE 1, column f^0(E).
      !
      ! UNITS AND NORMALIZATION, taken from the paper and not assumed.  Their
      ! Fig. 2 plots f^0(E) in eV^-1 with a x10^-1 axis multiplier and the
      ! table is in 10^-2 eV^-1; the absolute scale is fixed by the TRK sum
      ! rule, "normalized on an integral value of two (TRK sum rule) with a
      ! 1% correction for energies beyond 80 eV" (their Fig. 2 caption).  The
      ! conversion used here is the standard sigma[Mb] = 109.75 df/dE[1/eV].
      ! That constant is CONFIRMED against the other measurement rather than
      ! assumed: it puts Backx at 7.628 Mb at 20 eV and 2.700 at 30 eV
      ! against Samson & Haddad's 7.65 and 2.71, i.e. 0.3 and 0.4 per cent.
      ! (A factor 100 instead would put them 9 per cent low.)
      !
      ! THE JOIN.  Backx gives 9.691 Mb at 18 eV where Samson & Haddad give
      ! 9.850, so the branch is scaled by k_bx = 1.016415 to meet the table
      ! exactly.  THAT 1.6 PER CENT IS A MEASURED DIFFERENCE BETWEEN TWO
      ! INDEPENDENT EXPERIMENTS, not a fitted constant, and it is inside
      ! both quoted accuracies.
      !
      ! HOW BIG THE CORRECTION IS.  The branch this replaced was Yan Eq. 17
      ! scaled by a continuity constant.  Eq. 17 falls to 0.10 (0.17 after
      ! scaling) at 15.4 eV, where the measurement is 14.5: the analytic fit
      ! went to ZERO at the ionization threshold, which no photoionization
      ! cross section does.  The error ran from a factor 85 at threshold to
      ! 1.7 at 18 eV.
      !
      ! AN EARLIER VERSION OF THIS BRANCH was digitized from Fig. 1 of Chan,
      ! Cooper & Brion (1992), Chem. Phys. 168, 375, whose Tables 1-4 give
      ! only the discrete Lyman and Werner transitions.  That reading is kept
      ! in section 151.10 of docs/Update_EXHALE.md as a cross-check -- it
      ! agrees with the table used here to 2 per cent at 16 eV and 6 per cent
      ! at 17 eV -- and is not used in the code.  A table beats a figure.
      !
      ! WHAT IS EXCLUDED.  The table is evaluated only at and above the code
      ! threshold, 15.4 eV.  Its 15.0 eV row is carried so the interpolation
      ! into 15.4 has a left end; it is never returned.  Below the
      ! ionization potential the spectrum is the discrete Lyman and Werner
      ! band system, which this code carries separately as Lyman-Werner
      ! photodissociation and must not add to an ionization-continuum
      ! opacity twice.
      !
      ! THE SPREAD BETWEEN EXPERIMENTS IS THE REAL UNCERTAINTY.  Compared at
      ! the energies all three cover (section 151.10 has the table):
      ! Backx/Samson-Haddad runs 0.98 to 1.03 over 18-31 eV and 1.06 to 1.23
      ! over 32-68 eV; Lee, Carlson & Judge (1976), JQSRT 16, 873 runs 0.92
      ! to 0.96 over 18-31 eV and 1.0 to 1.5 above.  So sigma_H2 is good to
      ! a few per cent below about 30 eV and to 10-20 per cent above it,
      ! whatever any single paper's quoted accuracy says.
      !
      ! WHAT THE HIGH JOIN COSTS, and it is a real tension.  Eq. 19 agrees
      ! with the table to 0.33 per cent at 113 eV (0.033270 against
      ! 0.033160), so over the measured range the sum-rule tail and the
      ! measurement are the same curve.  They part company above it: at
      ! 300 eV the table reads 0.00154 and Eq. 19 gives 0.00174, and Yan et
      ! al. say why in their section 4 -- they RAISED the 300 eV point from
      ! 1.54e-21 to 1.75e-21 cm^2 precisely because the S(2) sum rules were
      ! not satisfied with the measured value.  So k_hi = 0.886 carries the
      ! measured normalization into the tail and gives up 12 per cent of
      ! their sum-rule correction.  Continuity is chosen over the sum rule
      ! here because a step in an opacity is what this whole change exists
      ! to remove, and because 300 eV is four decades below the peak of the
      ! cross section: sigma_H2(300) is 1.5e-3 against 9.85 at 18 eV.
      !
      ! Chung and Samson & Haddad tabulate megabarns, which are already
      ! 1e-18 cm^2; the Yan pieces are in barns and are converted
      ! (1 barn = 1e-6 Mb).
      double precision function sigma_H2(E)
      real*8, intent(in) :: E
      real*8 :: x, EkeV, sb
      integer :: k
      real*8, parameter :: eth = 15.4d0
      integer, parameter :: n_ab = 72
      ! Samson & Haddad (1994) Table 1, photon energy [eV] ...
      real*8, parameter :: e_ab(n_ab) = (/                                &
              18.0000d0,     19.0000d0,     20.0000d0,     21.0000d0,     22.0000d0,       &
              23.0000d0,     24.0000d0,     25.0000d0,     26.0000d0,     27.0000d0,       &
              28.0000d0,     29.0000d0,     30.0000d0,     31.0000d0,     32.0000d0,       &
              33.0000d0,     34.0000d0,     35.0000d0,     36.0000d0,     37.0000d0,       &
              38.0000d0,     39.0000d0,     40.0000d0,     42.0000d0,     44.0000d0,       &
              46.0000d0,     48.0000d0,     50.0000d0,     52.0000d0,     54.0000d0,       &
              56.0000d0,     58.0000d0,     60.0000d0,     62.0000d0,     64.0000d0,       &
              66.0000d0,     68.0000d0,     70.0000d0,     72.0000d0,     74.0000d0,       &
              76.0000d0,     78.0000d0,     80.0000d0,     85.0000d0,     90.0000d0,       &
              95.0000d0,    100.0000d0,    105.0000d0,    110.0000d0,    115.0000d0,       &
             120.0000d0,    125.0000d0,    130.0000d0,    135.0000d0,    140.0000d0,       &
             145.0000d0,    150.0000d0,    160.0000d0,    170.0000d0,    180.0000d0,       &
             190.0000d0,    200.0000d0,    210.0000d0,    220.0000d0,    230.0000d0,       &
             240.0000d0,    250.0000d0,    260.0000d0,    270.0000d0,    280.0000d0,       &
             290.0000d0,    300.0000d0 /)
      ! ... and their recommended sigma [1e-18 cm^2].
      real*8, parameter :: a_ab(n_ab) = (/                                &
              9.850000d0,     8.680000d0,     7.650000d0,     6.770000d0,     6.030000d0,       &
              5.400000d0,     4.830000d0,     4.310000d0,     3.850000d0,     3.500000d0,       &
              3.200000d0,     2.950000d0,     2.710000d0,     2.410000d0,     2.130000d0,       &
              1.930000d0,     1.740000d0,     1.570000d0,     1.420000d0,     1.290000d0,       &
              1.170000d0,     1.080000d0,     0.990000d0,     0.840000d0,     0.722000d0,       &
              0.626000d0,     0.543000d0,     0.477000d0,     0.420000d0,     0.373000d0,       &
              0.331000d0,     0.294000d0,     0.262000d0,     0.236000d0,     0.212000d0,       &
              0.191000d0,     0.173000d0,     0.157000d0,     0.143000d0,     0.130000d0,       &
              0.118000d0,     0.108000d0,     0.098000d0,     0.079200d0,     0.065800d0,       &
              0.055600d0,     0.047900d0,     0.041300d0,     0.035700d0,     0.031600d0,       &
              0.027100d0,     0.023800d0,     0.021100d0,     0.018600d0,     0.016500d0,       &
              0.014800d0,     0.013300d0,     0.010800d0,     0.008950d0,     0.007550d0,       &
              0.006400d0,     0.005450d0,     0.004700d0,     0.004030d0,     0.003500d0,       &
              0.003080d0,     0.002700d0,     0.002400d0,     0.002130d0,     0.001900d0,       &
              0.001710d0,     0.001540d0 /)
      integer, parameter :: n_lo = 7
      ! Backx et al. (1976) Table 1, photon energy [eV] ...
      real*8, parameter :: e_lo(n_lo) = (/                                &
              15.0000d0,     15.5000d0,     16.0000d0,     16.5000d0,     17.0000d0,       &
              17.5000d0,     18.0000d0 /)
      ! ... and k_bx * 109.75 * f^0(E), their absorption column [1e-18 cm^2].
      real*8, parameter :: a_lo(n_lo) = (/                                &
              14.8364d0,     14.0555d0,     13.0515d0,     12.0476d0,     11.2667d0,       &
              10.6866d0,      9.8500d0 /)

      if (E .lt. eth) then
         sigma_H2 = 0.0d0
         return
      endif

      if (E .lt. e_ab(1)) then
         ! The Backx table, interpolated in log sigma against log E.
         k = 1
         do while (k .lt. n_lo-1 .and. E .gt. e_lo(k+1))
            k = k + 1
         enddo
         sigma_H2 = exp(log(a_lo(k))                                      &
              + (log(a_lo(k+1)) - log(a_lo(k)))                           &
                *(log(E) - log(e_lo(k)))/(log(e_lo(k+1)) - log(e_lo(k))))
         return
      endif

      if (E .le. e_ab(n_ab)) then
         ! The measurement, interpolated in log sigma against log E.  Both
         ! are positive and smooth over the whole table, and a power law
         ! between adjacent rows is what the cross section is between them.
         k = 1
         do while (k .lt. n_ab-1 .and. E .gt. e_ab(k+1))
            k = k + 1
         enddo
         sigma_H2 = exp(log(a_ab(k))                                      &
              + (log(a_ab(k+1)) - log(a_ab(k)))                           &
                *(log(E) - log(e_ab(k)))/(log(e_ab(k+1)) - log(e_ab(k))))
         return
      endif

      ! Sum-rule tail, Eq. 19, rescaled to meet the table at 300 eV.
      x    = E/eth
      EkeV = E*1.0d-3
      sb = 45.57d0*(1.0d0 - 2.003d0/sqrt(x) - 4.806d0/x                  &
                    + 50.577d0/x**1.5d0 - 171.044d0/(x*x)                &
                    + 231.608d0/x**2.5d0 - 81.885d0/x**3)                &
           / EkeV**3.5d0
      if (sb .lt. 0.0d0) sb = 0.0d0
      sigma_H2 = sb*1.0d-6*sigma_H2_k_hi()

      end function sigma_H2

      !--------------!

      ! The one surviving continuity constant of sigma_H2: the ratio of the
      ! measured table endpoint at 300 eV to Eq. 19 evaluated there.

      double precision function sigma_H2_k_hi()
      real*8 :: x, EkeV, sb
      real*8, parameter :: E0 = 300.0d0, a0 = 0.00154d0
      x    = E0/15.4d0
      EkeV = E0*1.0d-3
      sb = 45.57d0*(1.0d0 - 2.003d0/sqrt(x) - 4.806d0/x                  &
                    + 50.577d0/x**1.5d0 - 171.044d0/(x*x)                &
                    + 231.608d0/x**2.5d0 - 81.885d0/x**3)                &
           / EkeV**3.5d0
      sigma_H2_k_hi = a0/(sb*1.0d-6)
      end function sigma_H2_k_hi

      !--------------!

      ! Fraction of an H2 PHOTOIONIZATION that is DISSOCIATIVE,
      !
      !     H2 + hv -> H + H+ + e-      (threshold 18.076 eV),
      !
      ! as opposed to the non-dissociative channel
      !
      !     H2 + hv -> H2+ + e-         (threshold 15.4 eV).
      !
      ! sigma_H2 above is the TOTAL photoionization cross section, so this is
      ! a BRANCHING of a cross section the code already carries and not an
      ! extra opacity: the rate at which photons destroy H2 does not change,
      ! only what the destruction leaves behind. The dissociative channel
      ! puts a proton and an H atom into the gas instead of an H2+ ion, so it
      ! bypasses the H3+ chain (H2+ + H2 -> H3+ + H) altogether.
      !
      ! SOURCE. Chung, Lee, Masuoka & Samson (1993), J. Chem. Phys. 99, 885,
      ! "Dissociative photoionization of H2 from 18 to 124 eV", their
      ! TABLE II, which tabulates sigma(H+) and sigma(H2+) at 69 energies
      ! from the threshold to 124 eV. The table below is
      !
      !     f_di(E) = sigma(H+) / [ sigma(H+) + sigma(H2+) ] ,
      !
      ! the share of the IONIZATIONS that release a proton, linearly
      ! interpolated in E. Their two partial columns reproduce their own
      ! ratio column to better than 0.5% at every row, and their sum
      ! reproduces their sigma(abs) column exactly except at the ten rows
      ! flagged with their footnote b (33-41 eV), where the photoionization
      ! yield is below unity and the difference is their Table I neutral
      ! cross section -- both checks were run on the transcription.
      !
      ! WHY NOT THE NUMBERS QUOTED BY YAN, SADEGHPOUR & DALGARNO (1998).
      ! Their section 4 says the dissociative "ratio to the total
      ! photoionization cross section increases from zero at the threshold of
      ! 18.08 eV to 0.284 at 76 eV and then decreases slowly to 0.258 at
      ! 124 eV". Those two numbers are Chung et al.'s H+/H2+ COLUMN, whose
      ! caption reads "the partial cross sections sigma(H+), sigma(H2+), and
      ! the ratio H+/H2+": they are a ratio to the H2+ partial cross section,
      ! not to the total. The ratio to the total is r/(1+r) -- 0.221 at
      ! 76 eV and 0.205 at 124 eV -- so taking the quoted numbers as a
      ! branching of sigma_H2 overstates the channel by a factor 1.26-1.28
      ! above 70 eV, and, because the true curve rises within 0.03 eV of
      ! threshold to a 0.02 plateau instead of ramping linearly, understates
      ! it by up to 20x below 20 eV. It also misses the Q1 Rydberg
      ! autoionization resonance, which takes f_di to 0.109 at 35 eV and back
      ! down to 0.091 at 37.5 eV.
      !
      ! MEASUREMENT UNCERTAINTY, as the source states it: "An estimate of the
      ! total uncertainty in the data is about +/-4%-5%"; the ionization
      ! yields used between 18 and 28 eV are unity "within an accuracy of
      ! about +/-5%". Against other determinations their H+ data lie 0%-5%
      ! below Backx, Wight & van der Wiel over 18-33 eV (26% at 70 eV), and
      ! within 5% of Kossmann et al. over 25-75 eV (11%-20% below them over
      ! 75-110 eV).
      !
      ! ABOVE 124 eV, the top of the measured range, f_di is held at its
      ! 124 eV value, 0.2048. The source quotes no higher-energy data. The
      ! curve is nearly flat over its last 50 eV (0.2219 at 74 eV to 0.2048
      ! at 124 eV), so holding it is a mild extrapolation, but it is an
      ! extrapolation.
      !
      ! TWO THINGS THIS BRANCHING CARRIES THAT IT DOES NOT NAME, both stated
      ! by the source and both confined to where sigma_H2 is already small.
      !  * DOUBLE IONIZATION. Chung et al.'s eqs. (3)-(4) note that "the
      !    branching ratio H+/(H+ + H2+) is understood to include both single
      !    and double ionization where appropriate", and their section on
      !    Fig. 4 says that "by 80 eV about 20% of sigma(H+) comes from
      !    double ionization. This percentage remains fairly constant towards
      !    higher energies." Because their sigma(H+) counts PROTONS and a
      !    double ionization yields two, using f_di as a share of the H2
      !    ionization rate gives the PROTON source rate correctly; what it
      !    mis-assigns is the co-product, since those events leave no H atom
      !    and no H2+. Above 80 eV that is about 4% of the H2 ionizations.
      !  * NEUTRAL DISSOCIATION over 33-41 eV. sigma_H2 above follows Samson
      !    & Haddad's ABSORPTION cross section there (it is 1.76 against
      !    Chung et al.'s sigma(abs) 1.74 and their ionization sum 1.699 at
      !    34 eV), so the H2 destruction rate the code applies in that band
      !    includes the neutral channel H2 + hv -> H + H, up to 7% of it at
      !    37.5 eV, and the 1 - f_di branch hands that share to H2+. This is
      !    a property of the cross-section fit, not of the branching.
      double precision function frac_H2_dissociative_ionization(E)
      real*8, intent(in) :: E
      integer :: k
      integer, parameter :: n_di = 69
      ! Chung et al. (1993) Table II: photon energy [eV] ...
      real*8, parameter :: e_di_tab(n_di) = (/                            &
          18.076d0,   18.100d0,   18.150d0,   18.200d0,   18.300d0,       &
          18.400d0,   18.500d0,   18.600d0,   18.800d0,   19.000d0,       &
          19.500d0,   20.000d0,   20.500d0,   21.000d0,   21.500d0,       &
          22.000d0,   23.000d0,   24.000d0,   25.000d0,   26.000d0,       &
          27.000d0,   28.000d0,   29.000d0,   30.000d0,   31.000d0,       &
          32.000d0,   33.000d0,   34.000d0,   35.000d0,   35.500d0,       &
          36.000d0,   36.500d0,   37.000d0,   37.500d0,   38.000d0,       &
          39.000d0,   40.000d0,   41.000d0,   42.000d0,   43.000d0,       &
          44.000d0,   45.000d0,   46.000d0,   48.000d0,   50.000d0,       &
          52.000d0,   54.000d0,   56.000d0,   58.000d0,   60.000d0,       &
          62.000d0,   64.000d0,   66.000d0,   68.000d0,   70.000d0,       &
          72.000d0,   74.000d0,   76.000d0,   78.000d0,   80.000d0,       &
          85.000d0,   90.000d0,   95.000d0,  100.000d0,  105.000d0,       &
         110.000d0,  115.000d0,  120.000d0,  124.000d0 /)
      ! ... and sigma(H+)/[sigma(H+) + sigma(H2+)] at those energies.
      real*8, parameter :: f_di_tab(n_di) = (/                            &
         0.00000d0,  0.00219d0,  0.00457d0,  0.00626d0,  0.00853d0,       &
         0.01009d0,  0.01133d0,  0.01268d0,  0.01461d0,  0.01613d0,       &
         0.01867d0,  0.01961d0,  0.02018d0,  0.02082d0,  0.02114d0,       &
         0.02156d0,  0.02222d0,  0.02298d0,  0.02365d0,  0.02454d0,       &
         0.02674d0,  0.03216d0,  0.03938d0,  0.04832d0,  0.06107d0,       &
         0.07512d0,  0.08759d0,  0.09947d0,  0.10845d0,  0.10940d0,       &
         0.10448d0,  0.09524d0,  0.09167d0,  0.09091d0,  0.09225d0,       &
         0.10049d0,  0.10892d0,  0.11641d0,  0.12381d0,  0.12949d0,       &
         0.13194d0,  0.13419d0,  0.13633d0,  0.14470d0,  0.15375d0,       &
         0.16370d0,  0.17203d0,  0.18152d0,  0.19020d0,  0.19817d0,       &
         0.20354d0,  0.20867d0,  0.21342d0,  0.21739d0,  0.22045d0,       &
         0.22160d0,  0.22188d0,  0.22119d0,  0.22037d0,  0.22041d0,       &
         0.21843d0,  0.21581d0,  0.21403d0,  0.21294d0,  0.21123d0,       &
         0.20964d0,  0.20786d0,  0.20635d0,  0.20477d0 /)

      if (E .le. e_di_tab(1)) then
         frac_H2_dissociative_ionization = 0.0d0
      else if (E .ge. e_di_tab(n_di)) then
         frac_H2_dissociative_ionization = f_di_tab(n_di)
      else
         k = 1
         do while (k .lt. n_di-1 .and. E .gt. e_di_tab(k+1))
            k = k + 1
         enddo
         frac_H2_dissociative_ionization = f_di_tab(k)                    &
              + (f_di_tab(k+1) - f_di_tab(k))                             &
                *(E - e_di_tab(k))/(e_di_tab(k+1) - e_di_tab(k))
      endif

      end function frac_H2_dissociative_ionization

      !--------------!

      ! C I photoionization (outer 2p shell). Verner+1996 Table 1 (full
      ! VFKY96 form): [E_0,sigma_0,y_a,P,y_w,y_0,y_1].
      double precision function sigma_CI(E)
      real*8, intent(in) :: E
      sigma_CI = sigma_VFKY96(E, 11.26d0, 2.144d0, 5.027d2, 6.216d1,  &
                                 5.101d0, 9.157d-2, 1.133d0, 1.607d0)
      end function sigma_CI

      !--------------!

      ! C II photoionization (outer 2p shell)
      double precision function sigma_CII(E)
      real*8, intent(in) :: E
      sigma_CII = sigma_VFKY96(E, 24.38d0, 4.058d-1, 8.709d0, 1.261d2, &
                                  8.578d0, 2.093d0, 4.929d1, 3.234d0)
      end function sigma_CII

      !--------------!

      ! O I photoionization (outer 2p shell)
      double precision function sigma_OI(E)
      real*8, intent(in) :: E
      sigma_OI = sigma_VFKY96(E, 13.62d0, 1.240d0, 1.745d3, 3.784d0,  &
                                 1.764d1, 7.589d-2, 8.698d0, 1.271d-1)
      end function sigma_OI

      !--------------!

      ! O II photoionization (outer 2p shell)
      double precision function sigma_OII(E)
      real*8, intent(in) :: E
      sigma_OII = sigma_VFKY96(E, 35.12d0, 1.386d0, 5.967d1, 3.175d1, &
                                  8.943d0, 1.934d-2, 2.131d1, 1.503d-2)
      end function sigma_OII

      !--------------!

      ! N I photoionization (outer 2p shell)
      double precision function sigma_NI(E)
      real*8, intent(in) :: E
      sigma_NI = sigma_VFKY96(E, 14.53d0, 4.034d0, 8.235d2, 8.033d1,  &
                                 3.928d0, 9.097d-2, 8.598d-1, 2.325d0)
      end function sigma_NI

      !--------------!

      ! N II photoionization (outer 2p shell)
      double precision function sigma_NII(E)
      real*8, intent(in) :: E
      sigma_NII = sigma_VFKY96(E, 29.60d0, 6.128d-2, 1.944d0, 8.163d2, &
                                  8.773d0, 1.043d1, 4.280d2, 2.030d1)
      end function sigma_NII

      !--------------!

      ! Mg I photoionization (outer 3s shell, l=0)
      ! Verner+1996 Table 1: E_th=7.646, E_0=11.97, sigma_0=1.372e8 Mb,
      ! y_a=0.2228, P=15.74, y_w=0.2805.
      double precision function sigma_MgI(E)
      real*8, intent(in) :: E
      sigma_MgI = sigma_Verner96(E, 7.646d0, 1.197d1, 1.372d8,    &
                                    2.228d-1, 1.574d1, 2.805d-1, 0)
      end function sigma_MgI

      !--------------!

      ! Mg II photoionization (outer 3s shell, l=0)
      ! Verner+1996 Table 1: E_th=15.04, E_0=8.139, sigma_0=3.278 Mb,
      ! y_a=4.341e7, P=3.610, y_w=0.0.
      double precision function sigma_MgII(E)
      real*8, intent(in) :: E
      sigma_MgII = sigma_Verner96(E, 1.504d1, 8.139d0, 3.278d0,   &
                                     4.341d7, 3.610d0, 0.0d0, 0)
      end function sigma_MgII

      !--------------!

      ! Si I photoionization (outer 3p shell). Verner+1996 Table 1
      ! (full VFKY96): E_th=8.152 eV; [E_0,sigma_0,y_a,P,y_w,y_0,y_1].
      double precision function sigma_SiI(E)
      real*8, intent(in) :: E
      sigma_SiI = sigma_VFKY96(E, 8.152d0, 2.317d1, 2.506d1, 2.057d1, &
                                  3.546d0, 2.837d-1, 1.672d-5, 4.207d-1)
      end function sigma_SiI

      !--------------!

      ! Si II photoionization (outer 3p shell)
      double precision function sigma_SiII(E)
      real*8, intent(in) :: E
      sigma_SiII = sigma_VFKY96(E, 16.35d0, 2.556d0, 4.140d0, 1.337d1, &
                                   1.191d1, 1.570d0, 6.634d0, 1.272d-1)
      end function sigma_SiII

      !--------------!

      ! Ca I photoionization (outer 4s shell). The very large sigma_0 is
      ! offset by the formula; threshold value is ~7 Mb.
      double precision function sigma_CaI(E)
      real*8, intent(in) :: E
      sigma_CaI = sigma_VFKY96(E, 6.113d0, 1.278d1, 5.370d5, 3.162d-1, &
                                  1.242d1, 4.477d-1, 1.012d-3, 1.851d-2)
      end function sigma_CaI

      !--------------!

      ! Ca II photoionization (outer 4s shell)
      double precision function sigma_CaII(E)
      real*8, intent(in) :: E
      sigma_CaII = sigma_VFKY96(E, 11.87d0, 1.553d1, 1.064d7, 7.790d-1, &
                                   2.130d1, 6.453d-1, 2.161d-3, 6.706d-2)
      end function sigma_CaII

      !--------------!

      ! Na I photoionization (outer 3s shell). Verner+1996 Table 1 has
      ! y_0=y_1=0, so the VY95 form (l=0) reproduces it exactly.
      double precision function sigma_NaI(E)
      real*8, intent(in) :: E
      sigma_NaI = sigma_Verner96(E, 5.139d0, 6.139d0, 1.601d0,    &
                                    6.148d3, 3.839d0, 0.0d0, 0)
      end function sigma_NaI

      !--------------!

      ! K I photoionization (outer 4s shell). K I is not in the VFKY96
      ! Opacity-Project table; use the Verner & Yakovlev 1995 (A&AS 109,
      ! 125) fit from phfit2.f ph1(19,19,7) [E_th,E_0,sigma_0,y_a,P,y_w],
      ! l=0.
      double precision function sigma_KI(E)
      real*8, intent(in) :: E
      sigma_KI = sigma_Verner96(E, 4.341d0, 3.824d0, 7.363d-1,    &
                                   2.410d7, 4.427d0, 2.049d-4, 0)
      end function sigma_KI

      !--------------!

      ! S I photoionization (outer 3p shell)
      double precision function sigma_SI(E)
      real*8, intent(in) :: E
      sigma_SI = sigma_VFKY96(E, 10.36d0, 1.808d1, 4.564d4, 1.000d0, &
                                 1.361d1, 6.385d-1, 9.935d-1, 2.486d-1)
      end function sigma_SI

      !--------------!

      ! Fe I photoionization. Verner+1996 Table 1 single-shell VFKY96 fit
      ! [E_0,sigma_0,y_a,P,y_w,y_0,y_1] = PH2(:,26,26) (outer 3d/4s shell);
      ! threshold = NIST IP 7.902 eV. !To Be Checked/AIOLOS tuning?  This
      ! single-shell fit underestimates inner-shell (L,M) XUV absorption and
      ! differs from Huang+2023, who adopt the Zatsarinny+2019 R-matrix Fe I
      ! cross section. Adequate for the ionization balance (Fe is mostly
      ! ionized in the thermosphere); revisit if XUV heating is sensitive.
      double precision function sigma_FeI(E)
      real*8, intent(in) :: E
      sigma_FeI = sigma_VFKY96(E, 7.902d0, 5.461d-2, 3.062d-1, 2.671d7, &
                                  7.923d0, 2.069d1, 1.382d2, 2.481d-1)
      end function sigma_FeI

      !--------------!

      ! Fe II photoionization. Verner+1996 Table 1 single-shell VFKY96 fit
      ! [E_0,sigma_0,y_a,P,y_w,y_0,y_1] = PH2(:,26,25); threshold = NIST IP
      ! 16.199 eV.
      double precision function sigma_FeII(E)
      real*8, intent(in) :: E
      sigma_FeII = sigma_VFKY96(E, 16.199d0, 1.761d-1, 4.365d3, 6.298d3, &
                                   5.204d0, 1.141d1, 9.272d1, 1.075d2)
      end function sigma_FeII

      !----------------------------------------

      ! Photoionization cross section [1e-18 cm^2] of the metal ion whose
      ! photo-table column index is k, i.e. species_table's mion_iphot(i).
      ! THE COLUMN ORDER IS DEFINED HERE, and set_energy_vectors fills
      ! sigma_tab through this function, so the ordering exists once.
      ! It is needed off the e_v grid because the He recombination channels
      ! (util_ion_eq, he_rec_coupling) emit at their own photon energies.
      ! Returns zero for k outside 1..n_mphot, and each fit returns zero
      ! below its own threshold.
      double precision function metal_photoion_sigma(k,E)
      integer, intent(in) :: k
      real*8,  intent(in) :: E

      select case (k)
         case ( 1) ; metal_photoion_sigma = sigma_CI  (E)
         case ( 2) ; metal_photoion_sigma = sigma_CII (E)
         case ( 3) ; metal_photoion_sigma = sigma_OI  (E)
         case ( 4) ; metal_photoion_sigma = sigma_OII (E)
         case ( 5) ; metal_photoion_sigma = sigma_NI  (E)
         case ( 6) ; metal_photoion_sigma = sigma_NII (E)
         case ( 7) ; metal_photoion_sigma = sigma_MgI (E)
         case ( 8) ; metal_photoion_sigma = sigma_MgII(E)
         case ( 9) ; metal_photoion_sigma = sigma_SiI (E)
         case (10) ; metal_photoion_sigma = sigma_SiII(E)
         case (11) ; metal_photoion_sigma = sigma_CaI (E)
         case (12) ; metal_photoion_sigma = sigma_CaII(E)
         case (13) ; metal_photoion_sigma = sigma_NaI (E)
         case (14) ; metal_photoion_sigma = sigma_KI  (E)
         case (15) ; metal_photoion_sigma = sigma_SI  (E)
         case (16) ; metal_photoion_sigma = sigma_FeI (E)
         case (17) ; metal_photoion_sigma = sigma_FeII(E)
         case default ; metal_photoion_sigma = 0.0d0
      end select

      end function metal_photoion_sigma

      !----------------------------------------

      ! End of module
      end module Cross_sections