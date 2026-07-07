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

      ! H2 total photoionization cross section (molecular extension).
      ! Yan, Sadeghpour & Dalgarno (1998), ApJ 496, 1044, Eqs. 17-19 (fit
      ! coefficients verified against the PDF, archived as
      ! references/Yan_1998_ApJ_496_1044.pdf).  Threshold 15.4 eV; pieces
      ! join at 18 and 85 eV; sum-rule-consistent E^-7/2 tail (their Eq. 16,
      ! sigma -> 45.6/E[keV]^3.5 barns; ~2.8x atomic H at high E).  Paper
      ! units are barns; converted here to EXHALE units of 1e-18 cm^2
      ! (1 barn = 1e-6 Mb).  Near-threshold resonance structure is smoothed
      ! (as in the source fits).
      double precision function sigma_H2(E)
      real*8, intent(in) :: E
      real*8 :: x, EkeV, sb
      real*8, parameter :: eth = 15.4d0
      real*8, parameter :: s   = 0.252d0

      if (E .lt. eth) then
         sigma_H2 = 0.0d0
         return
      endif

      x = E/15.4d0
      if (E .lt. 18.0d0) then
         ! Eq. 17 (15.4 < E < 18 eV), barns
         sb = 1.0d8*(-37.895d0 + 99.723d0*x - 87.227d0*x*x               &
                     + 25.400d0*x*x*x)
      else if (E .lt. 85.0d0) then
         ! Eq. 18 (18 < E < 85 eV), barns
         sb = 2.0d7*( 0.071d0*x**(-s)         - 0.673d0*x**(-(s+1.0d0))  &
                    + 1.977d0*x**(-(s+2.0d0)) - 0.692d0*x**(-(s+3.0d0)) )
      else
         ! Eq. 19 (E > 85 eV), barns
         EkeV = E*1.0d-3
         sb = 45.57d0*(1.0d0 - 2.003d0/sqrt(x) - 4.806d0/x               &
                       + 50.577d0/x**1.5d0 - 171.044d0/(x*x)             &
                       + 231.608d0/x**2.5d0 - 81.885d0/x**3)             &
              / EkeV**3.5d0
      endif
      if (sb .lt. 0.0d0) sb = 0.0d0
      sigma_H2 = sb*1.0d-6                 ! barns -> 1e-18 cm^2

      end function sigma_H2

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

      ! End of module
      end module Cross_sections