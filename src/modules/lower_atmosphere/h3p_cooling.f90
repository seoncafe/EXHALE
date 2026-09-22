      module h3p_cooling
      ! Optically-thin H3+ infrared cooling (the molecular level of the lower-atmosphere
      ! plan; docs/lower_atmosphere_coupling.*).
      !
      ! WHAT THE PUBLISHED FIT IS.  Miller, Stallard, Tennyson & Melin
      ! (2013), J. Phys. Chem. A 117, 9770, Table 5 gives
      !
      !   log_e E(H3+,T) = sum_n C_n T^n     [E in W molecule^-1 sr^-1],
      !
      ! the emission of ONE H3+ molecule into ONE steradian, computed for an
      ! LTE internal population at T with the partition function z(T) of
      ! their Table 1 (the "using z(t)" column of Table 5 is the one used
      ! here; the alternative zfp(T) column excludes their extrapolated
      ! levels).  The population convention is therefore LTE at the local
      ! gas temperature, and the isotropic volumetric rate is
      !
      !   Lambda_H3+ = n_H3+ * 4 pi * E(T) * s(T, n_H2)   [W cm^-3]
      !   (multiply by 1e7 for erg s^-1 cm^-3),
      !
      ! the 4 pi being the one conversion from the per-steradian fit.
      !
      ! THE NON-LTE FACTOR.  s(T, [H2]) is the departure factor of their
      ! Table 6, defined by their eq. 8 as the ratio of the vibrational-level
      ! sum with the detailed-balance (Oka & Epp) populations to the same sum
      ! with the LTE populations.  It is referenced to the SAME LTE emission
      ! E(T), so it multiplies the fit and is not a separate emission model.
      ! Its only collider is H2, through the proton-hopping reaction
      ! H3+(v=m) + H2 -> H2 + H3+(v=n) with the Oka & Epp rate coefficient
      ! 2e-15 m^3 s^-1; no other collision partner enters, and the density
      ! it depends on is the H2 number density alone.  The paper states that
      ! the Table 6 values MAY BE UPPER LIMITS, especially at low densities,
      ! because that rate coefficient is its largest uncertainty.
      !
      ! COLLIDER DOMAIN AND THE LOW-COLLIDER LIMIT.  Table 6 is tabulated at
      ! [H2] = 1e12 to 1e20 m^-3, that is 1e6 to 1e14 cm^-3.  Below 1e6 cm^-3
      ! radiative decay empties the emitting levels faster than collisions
      ! populate them, so each collisional excitation is followed by a decay
      ! and the emission per H3+ is set by the collision rate: s is linear in
      ! n(H2) and vanishes with it.  Decision 8 of
      ! docs/b1_target_system_20260906.md section 6.1 adopts that limit as
      ! the evaluation below the table,
      !
      !   s(T, n_H2 < 1e6) = s(T, 1e6) * n_H2/1e6,
      !
      ! and records the evaluation as made below the tabulated collider
      ! range (h3p_evaluations_below_collider below).  The record is
      ! informational: the coolant stays defined where the model has a known
      ! analytic limit.
      ! Above 1e14 cm^-3 the table edge is held (s = s(T,1e14), which is
      ! 0.9985 to 1.0000, so the residual departure from LTE at and above the
      ! last tabulated column is at most 0.15 per cent); no extrapolation is
      ! made above the table, and the evaluation is recorded there too
      ! (h3p_evaluations_above_collider).
      !
      ! TEMPERATURE DOMAIN AND THE JOINS.  The four Table 5 segments cover
      ! 30-300, 300-800, 800-1800 and 1800-5000 K.  Quoted fit errors: up to
      ! +-5 per cent over 30-300 K, generally < 0.1 per cent over 300-800 K,
      ! "all but nonexistent" over 800-1800 K, a few tenths of a per cent
      ! over 1800-5000 K; above 5000 K the paper does not report values.
      ! Outside 30-5000 K the evaluation is clamped to the end of the fit and
      ! recorded (h3p_evaluations_below_fit_T, h3p_evaluations_above_fit_T);
      ! H3+ and its feedstock H2 are thermally destroyed well below 5000 K,
      ! so the frozen tail only guards transients.  30 K and 5000 K are ON
      ! the fit and are inside it: the records below count STRICT excursions
      ! only, and the clamp, which is written NaN-safely and sends a
      ! nonfinite temperature to an end of the fit, is a separate statement
      ! from the record.
      !
      ! The published segments do NOT join, and the mismatch is theirs, not a
      ! transcription error.  MEASURED here with the Table 5 coefficients,
      ! continuing each segment across its upper limit:
      !
      !   at  300 K   the 300-800 K polynomial is  +43.05 per cent
      !   at  800 K   the 800-1800 K polynomial is  -0.113 per cent
      !   at 1800 K   the 1800-5000 K polynomial is -2.353 per cent
      !
      ! relative to the segment below it.  At 300 K the upper segment is the
      ! one the paper's own tables support: Table 6 gives the LTE emission at
      ! 300 K as 0.53503e-23 W molecule^-1 sr^-1, which the 300-800 K
      ! polynomial reproduces to 0.002 per cent while the 30-300 K one stands
      ! 30.1 per cent below it (MEASURED).
      !
      ! JOIN RULE OF THIS CODE.  A cooling function with a 43 per cent step
      ! in it is not a function of state, so the code interpolates across the
      ! published gap and says so here rather than hiding it.  Over
      ! [T_join*(1 - 0.05), T_join) the two published polynomials are blended
      ! linearly in log_e E, weight rising from the lower segment to the
      ! upper one; at and above T_join the upper published segment is exact,
      ! and below the ramp the lower published segment is exact.  The ramp is
      ! this code's rule, not the paper's.  It leaves every value the paper
      ! tabulates untouched: the Table 4 entries at 500-5000 K and the
      ! Table 6 LTE value at 300 K all lie at or outside the ramps, and the
      ! 500-5000 K fit reproduces Table 4 to within 0.45 per cent (MEASURED,
      ! worst point 5000 K).
      !
      ! Reference PDF: references/Miller_2013_JPCA_117_9770.pdf
      ! (Tables 4, 5, 6).  NaN-safe bracketing style of the cooling-table
      ! interpolators (code review 2026-07-02).

      implicit none
      private
      public :: h3p_emission_lte, h3p_nonlte_factor, h3p_cooling_rate,   &
                h3p_net_cooling_rate, h3p_reset_domain_records
      public :: h3p_evaluations_below_collider,                           &
                h3p_evaluations_above_collider,                           &
                h3p_evaluations_below_fit_T,                              &
                h3p_evaluations_above_fit_T,                              &
                h3p_evaluations_outside_nonlte_T,                         &
                h3p_evaluations_nonfinite_T,                              &
                h3p_evaluations_nonfinite_collider
      ! Where one argument stands relative to a published range, for a
      ! caller that maps a STATE rather than counting the evaluations of a
      ! run: the same three tests the records below are taken with, with no
      ! side effect of their own.
      public :: h3p_fit_temperature_domain,                               &
                h3p_nonlte_row_temperature_domain, h3p_collider_domain
      public :: H3P_DOMAIN_INSIDE, H3P_DOMAIN_BELOW, H3P_DOMAIN_ABOVE,    &
                H3P_DOMAIN_NONFINITE

      real*8, parameter :: fourpi = 12.566370614359172d0

      ! WHAT THESE COUNT.  Each is the number of EVALUATIONS, accumulated
      ! over the whole run, whose argument met the stated condition.  A cell
      ! visited many times contributes many times, one cell can meet several
      ! conditions, and an evaluation made inside a trial that was then
      ! rejected contributes like one inside an accepted state.  They are
      ! therefore the history of the run and are not a property of any one
      ! state; the domain map of a single state is a separate measurement
      ! (certification.f90, h3p_cooling_domain_map_of_state), taken with the
      ! three classifying functions below.  Informational either way: none
      ! of them changes a rate.
      integer, save :: h3p_evaluations_below_collider     = 0  ! n_H2 < 1e6 cm^-3 (collisional limit used)
      integer, save :: h3p_evaluations_above_collider     = 0  ! n_H2 > 1e14 cm^-3 (table edge held)
      integer, save :: h3p_evaluations_below_fit_T        = 0  ! T < 30 K   (emission clamped)
      integer, save :: h3p_evaluations_above_fit_T        = 0  ! T > 5000 K (emission clamped)
      integer, save :: h3p_evaluations_outside_nonlte_T   = 0  ! T outside 300-5000 K (Table 6 row clamped)
      integer, save :: h3p_evaluations_nonfinite_T        = 0  ! T not an ordinary real
      integer, save :: h3p_evaluations_nonfinite_collider = 0  ! n_H2 not an ordinary real

      ! Where an argument stands relative to a published range.  BELOW and
      ! ABOVE are strict excursions, so an argument sitting exactly on an
      ! endpoint is INSIDE: the fit is defined there and the table has a row
      ! there.  NONFINITE is its own answer and is never folded into BELOW.
      integer, parameter :: H3P_DOMAIN_INSIDE    = 0
      integer, parameter :: H3P_DOMAIN_BELOW     = 1
      integer, parameter :: H3P_DOMAIN_ABOVE     = 2
      integer, parameter :: H3P_DOMAIN_NONFINITE = 3

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

      ! Segment limits and the width of the join ramp, as a fraction of the
      ! join temperature (this code's rule; see JOIN RULE above).
      real*8, parameter :: T_fit_lo = 30.0d0, T_fit_hi = 5000.0d0
      real*8, parameter :: T_join(3) = [ 300.0d0, 800.0d0, 1800.0d0 ]
      real*8, parameter :: join_ramp = 0.05d0

      ! --- Table 6 non-LTE scaling s(T, n_H2) --------------------------- !
      ! Grid: T rows (K) x log10(n_H2 [cm^-3]) columns 6,8,10,12,14
      ! ([H2] = 1e12..1e20 m^-3 in the paper).
      integer, parameter :: nTs = 11, nNs = 5
      real*8, parameter :: sT(nTs) = [ 300.d0, 600.d0, 1000.d0, 1500.d0, &
           2000.d0, 2500.d0, 3000.d0, 3500.d0, 4000.d0, 4500.d0, 5000.d0 ]
      real*8, parameter :: sLogN(nNs) = [ 6.d0, 8.d0, 10.d0, 12.d0, 14.d0 ]
      real*8, parameter :: n_tab_lo = 1.0d6   ! lowest tabulated collider density [cm^-3]
      real*8, parameter :: n_tab_hi = 1.0d14  ! highest, 10**sLogN(nNs) [cm^-3]
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

      ! Zero the evaluation history above, so that what follows is counted
      ! on its own.
      subroutine h3p_reset_domain_records()
      h3p_evaluations_below_collider     = 0
      h3p_evaluations_above_collider     = 0
      h3p_evaluations_below_fit_T        = 0
      h3p_evaluations_above_fit_T        = 0
      h3p_evaluations_outside_nonlte_T   = 0
      h3p_evaluations_nonfinite_T        = 0
      h3p_evaluations_nonfinite_collider = 0
      end subroutine h3p_reset_domain_records

      ! ------------------------------------------------------------------ !

      ! Whether x is an ordinary real: no NaN and no infinity.  NaN fails
      ! every comparison including with itself and an infinity exceeds the
      ! largest representable finite value, so the two tests together cover
      ! both.  Written out rather than taken from ieee_arithmetic so that the
      ! generated module dependency graph stays over the source tree, as
      ! finite_real (ionization_equilibrium.f90) is.
      pure logical function ordinary_real(x) result(ok)
      real*8, intent(in) :: x
      ok = (x .eq. x) .and. (abs(x) .le. huge(1.0d0))
      end function ordinary_real

      ! Where T stands relative to the 30-5000 K range of the Table 5 fits.
      pure integer function h3p_fit_temperature_domain(T) result(code)
      real*8, intent(in) :: T
      if (.not. ordinary_real(T)) then
         code = H3P_DOMAIN_NONFINITE
      else if (T .lt. T_fit_lo) then
         code = H3P_DOMAIN_BELOW
      else if (T .gt. T_fit_hi) then
         code = H3P_DOMAIN_ABOVE
      else
         code = H3P_DOMAIN_INSIDE
      endif
      end function h3p_fit_temperature_domain

      ! Where T stands relative to the 300-5000 K rows of Table 6.
      pure integer function h3p_nonlte_row_temperature_domain(T)          &
                            result(code)
      real*8, intent(in) :: T
      if (.not. ordinary_real(T)) then
         code = H3P_DOMAIN_NONFINITE
      else if (T .lt. sT(1)) then
         code = H3P_DOMAIN_BELOW
      else if (T .gt. sT(nTs)) then
         code = H3P_DOMAIN_ABOVE
      else
         code = H3P_DOMAIN_INSIDE
      endif
      end function h3p_nonlte_row_temperature_domain

      ! Where the collider density stands relative to the 1e6-1e14 cm^-3
      ! columns of Table 6.  BELOW is the analytic collisional limit of the
      ! module header, which is a statement of the model and not an
      ! extrapolation outside it; ABOVE is the held table edge.
      pure integer function h3p_collider_domain(nH2) result(code)
      real*8, intent(in) :: nH2
      if (.not. ordinary_real(nH2)) then
         code = H3P_DOMAIN_NONFINITE
      else if (nH2 .lt. n_tab_lo) then
         code = H3P_DOMAIN_BELOW
      else if (nH2 .gt. n_tab_hi) then
         code = H3P_DOMAIN_ABOVE
      else
         code = H3P_DOMAIN_INSIDE
      endif
      end function h3p_collider_domain

      ! ------------------------------------------------------------------ !

      ! LTE emission per molecule E(H3+,T) [W molecule^-1 sr^-1], Table 5,
      ! with the join ramp of the module header.
      double precision function h3p_emission_lte(T) result(E)
      real*8, intent(in) :: T
      real*8 :: tt, lnE, t_ramp, f
      integer :: k
      ! THE CLAMP, which a nonfinite temperature also has to land on: the
      ! .not.(x .gt. y) form sends a NaN to the lower end of the fit.
      tt = T
      if (.not. (tt .gt. T_fit_lo)) then            ! NaN-safe lower clamp
         tt = T_fit_lo
      else if (tt .gt. T_fit_hi) then               ! frozen high-T tail
         tt = T_fit_hi
      endif
      ! THE RECORD, which is interval membership and not the clamp: T = 30 K
      ! and T = 5000 K are on the fit and count as nothing.
      select case (h3p_fit_temperature_domain(T))
      case (H3P_DOMAIN_BELOW)
!$omp atomic update
         h3p_evaluations_below_fit_T = h3p_evaluations_below_fit_T + 1
      case (H3P_DOMAIN_ABOVE)
!$omp atomic update
         h3p_evaluations_above_fit_T = h3p_evaluations_above_fit_T + 1
      case (H3P_DOMAIN_NONFINITE)
!$omp atomic update
         h3p_evaluations_nonfinite_T = h3p_evaluations_nonfinite_T + 1
      end select
      ! The segment whose published range contains tt; at a join the upper
      ! segment owns the point (see JOIN RULE).
      k = 1
      if (tt .ge. T_join(1)) k = 2
      if (tt .ge. T_join(2)) k = 3
      if (tt .ge. T_join(3)) k = 4
      lnE = ln_emission_segment(k, tt)
      if (k .le. 3) then
         t_ramp = T_join(k)*(1.0d0 - join_ramp)
         if (tt .gt. t_ramp) then
            f   = (tt - t_ramp)/(T_join(k) - t_ramp)
            lnE = (1.0d0 - f)*lnE + f*ln_emission_segment(k+1, tt)
         endif
      endif
      E = exp(lnE)
      end function h3p_emission_lte

      ! log_e E of Table 5 segment k (1: 30-300, 2: 300-800, 3: 800-1800,
      ! 4: 1800-5000 K), evaluated at T without regard to that range.
      double precision function ln_emission_segment(k, T) result(lnE)
      integer, intent(in) :: k
      real*8, intent(in)  :: T
      if (k .eq. 1) then
         lnE = poly(cA, 9, T)
      else if (k .eq. 2) then
         lnE = poly(cB, 6, T)
      else if (k .eq. 3) then
         lnE = poly(cC, 6, T)
      else
         lnE = poly(cD, 5, T)
      endif
      end function ln_emission_segment

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

      ! Non-LTE departure factor s(T, n_H2), Miller+2013 Table 6, bilinear in
      ! (T, log10 n_H2 [cm^-3]).  Below the lowest tabulated collider density
      ! the collisional limit of the module header is used, so s and the
      ! cooling vanish linearly with n_H2; above the highest tabulated one the
      ! table edge is held; outside 300-5000 K the nearest tabulated row is
      ! used.  Each of those three is recorded, and so is a nonfinite
      ! argument.
      double precision function h3p_nonlte_factor(T, nH2) result(s)
      real*8, intent(in) :: T, nH2
      real*8 :: tt, nn, ln
      ! The records first, taken on the arguments as given: 300 K and
      ! 5000 K are tabulated rows and 1e6 and 1e14 cm^-3 tabulated columns,
      ! so an argument sitting on one of them counts as nothing.
      select case (h3p_nonlte_row_temperature_domain(T))
      case (H3P_DOMAIN_BELOW, H3P_DOMAIN_ABOVE)
!$omp atomic update
         h3p_evaluations_outside_nonlte_T =                               &
              h3p_evaluations_outside_nonlte_T + 1
      case (H3P_DOMAIN_NONFINITE)
!$omp atomic update
         h3p_evaluations_nonfinite_T = h3p_evaluations_nonfinite_T + 1
      end select
      select case (h3p_collider_domain(nH2))
      case (H3P_DOMAIN_BELOW)
!$omp atomic update
         h3p_evaluations_below_collider =                                 &
              h3p_evaluations_below_collider + 1
      case (H3P_DOMAIN_ABOVE)
!$omp atomic update
         h3p_evaluations_above_collider =                                 &
              h3p_evaluations_above_collider + 1
      case (H3P_DOMAIN_NONFINITE)
!$omp atomic update
         h3p_evaluations_nonfinite_collider =                             &
              h3p_evaluations_nonfinite_collider + 1
      end select
      ! The evaluation, with the clamps as they are: a nonfinite argument
      ! lands on the low end of the rows and on a zero collider density,
      ! where the collisional limit returns zero.
      tt = T
      if (.not. (tt .gt. sT(1))) tt = sT(1)         ! NaN-safe
      if (tt .gt. sT(nTs))       tt = sT(nTs)
      nn = nH2
      if (.not. (nn .gt. 0.0d0)) nn = 0.0d0         ! NaN-safe
      if (nn .lt. n_tab_lo) then
         s = s_table(tt, sLogN(1))*(nn/n_tab_lo)    ! collisional limit
         return
      endif
      ln = log10(nn)
      if (ln .gt. sLogN(nNs)) ln = sLogN(nNs)       ! table edge held
      s = s_table(tt, ln)
      end function h3p_nonlte_factor

      ! Bilinear interpolation of Table 6 on its own grid; both arguments
      ! must already lie inside it.
      double precision function s_table(tt, ln) result(s)
      real*8, intent(in) :: tt, ln
      real*8 :: ft, fn
      integer :: it, in
      it = 1
      do while (it .lt. nTs-1 .and. sT(it+1) .lt. tt)
         it = it + 1
      enddo
      ft = (tt - sT(it))/(sT(it+1) - sT(it))
      in = 1 + int((ln - sLogN(1))/2.0d0)           ! uniform step 2 in log10
      if (in .gt. nNs-1) in = nNs-1
      fn = (ln - sLogN(in))/2.0d0
      s  = (1.d0-ft)*(1.d0-fn)*sTab(it,  in  )                            &
         +       ft *(1.d0-fn)*sTab(it+1,in  )                            &
         + (1.d0-ft)*      fn *sTab(it,  in+1)                            &
         +       ft *      fn *sTab(it+1,in+1)
      end function s_table

      ! ------------------------------------------------------------------ !

      ! Volumetric H3+ cooling rate [erg s^-1 cm^-3]:
      !   Lambda = n_H3+ * 4 pi * E_LTE(T) * s_nonLTE(T,n_H2) * 1e7.
      double precision function h3p_cooling_rate(T, nH3p, nH2) result(lam)
      real*8, intent(in) :: T, nH3p, nH2
      lam = nH3p*fourpi*h3p_emission_lte(T)*h3p_nonlte_factor(T,nH2)*1.0d7
      end function h3p_cooling_rate

      ! ------------------------------------------------------------------ !

      ! NET H3+ infrared cooling [erg s^-1 cm^-3]: the emission above minus
      ! the absorption of the thermal infrared radiation of the lower
      ! atmosphere, a blackbody at T_rad filling a fraction W_dil of the
      ! solid angle (W_dil = 0 recovers h3p_cooling_rate exactly).
      !
      ! CLOSURE.  The same net-exchange form the other infrared channels of
      ! the molecular layer use (molecular_infrared_cooling.f90, eq. 1):
      !
      !   Lambda_net = n_H3+ 4 pi s(T,n_H2) [ E(T) - W_dil E(T_rad) ]
      !
      ! The absorbed power per molecule is
      ! 4 pi W_dil Int sigma_nu B_nu(T_rad) dnu, and Kirchhoff's law at the
      ! radiating temperature identifies that integral, evaluated with the
      ! cross section AT T_rad, with the emission fit itself: E(T_rad).  So
      ! the emission fit does double duty and no line list is needed.
      !
      ! WHY NOT ONE EFFECTIVE TRANSITION.  Between 2026-08-13 and 2026-08-31
      ! this routine collapsed the emission onto a single nu2-band transition
      ! of energy hc*2521.3 cm^-1 / k = 3627.5 K and wrote the absorption as
      ! Lambda_emit*nbar*[exp(Ek/T) - 1].  That pairs the TOTAL emission fit
      ! of Miller et al. -- a sum over the nu2 fundamental, its hot bands, the
      ! overtones and the forbidden rotational lines -- with the Boltzmann
      ! factor of ONE of them.  It is self-consistent only where
      ! E(T) exp(Ek/T) is independent of T, which is what a single transition
      ! with a saturated lower level would give; measured on the Table 5 fit
      ! that product runs over forty decades, 6.7e-19 at 300 K to 1.5e22 at
      ! 30 K, because the low-temperature emission is carried by transitions
      ! far below 3627.5 K.  The bracket was therefore already -3.9e3 at
      ! 300 K and -1.2e14 at 100 K -- heating the gas at 1e14 times its own
      ! emission -- and the exp() cap at 700 only stopped it from being Inf.
      ! The form above is bounded for any T because the absorption does not
      ! reference the gas temperature at all.
      !
      ! The three limits it has to have, and does:
      !   T = T_rad, W_dil = 1  ->  exactly zero, to machine precision:
      !       gas buried in a blackbody at its own temperature neither cools
      !       nor heats.  This is the fixed point section 110 of
      !       docs/Update_EXHALE_stage1.md relies on for TO_BE_DONE item (G).
      !   T -> 0                ->  emission -> 0, absorption finite: the
      !       band heats at the rate the field supplies, not faster.
      !   W_dil = 0             ->  the emission-only rate, unchanged.
      ! The radiative equilibrium temperature is the root of
      ! E(T_eq) = W_dil E(T_rad); it has no closed form and is 975.4 K for
      ! T_rad = 1140 K, W_dil = 1/2 (the single-transition form gave 941 K).
      !
      ! WHERE IT DIFFERS FROM THE OTHER THREE CHANNELS.  H2O and CO are band
      ! absorbers with tabulated sigma_b(T), so their absorption keeps the
      ! cross section at the GAS temperature and only the Planck function at
      ! T_rad; H2 is summed transition by transition with the full two-level
      ! factor.  Neither is available here: Miller et al. (2013) publish the
      ! product of cross section and Planck function already integrated, so
      ! the cross section cannot be separated from the weighting and is
      ! frozen at T_rad.  The error that freezing makes is the temperature
      ! dependence of the band strength alone -- the fraction of molecules in
      ! the absorbing lower states, which rises toward 1 as T falls, so the
      ! absorption above is if anything a slight underestimate at low T.  It
      ! is exact at T = T_rad, which is where the fixed point sits.
      !
      ! VALIDITY.  (i) The lower atmosphere is taken to be black at 3-4 um,
      ! which an H2 atmosphere at the microbar base and below is (H2
      ! collision-induced absorption plus the H2O/CH4/CO bands), but the
      ! incident field is NOT attenuated by the intervening H3+ column, so
      ! where those lines are self-shielding this overestimates the
      ! absorption -- the same trapping would also reduce the emission term,
      ! which is not modeled either.  (ii) The non-LTE departure factor
      ! multiplies BOTH terms.  It has to: s < 1 with the absorption left at
      ! its LTE value would heat a gas that is HOTTER than the field, which
      ! no exchange between two temperatures may do.  Read this way s is the
      ! fraction of molecules radiatively coupled at the local density, and
      ! it scales the rate without moving the equilibrium temperature.
      ! Miller et al.'s Table 6 is a vacuum departure factor -- collisions
      ! against spontaneous decay, no incident field -- and a field strong
      ! enough to matter here pumps the levels back toward LTE, so the true
      ! departure is smaller than s and this underestimates the magnitude of
      ! the exchange, not its sign or its fixed point.
      double precision function h3p_net_cooling_rate                      &
                                 (T, nH3p, nH2, T_rad, W_dil) result(lam)
      real*8, intent(in) :: T, nH3p, nH2, T_rad, W_dil
      lam = h3p_cooling_rate(T, nH3p, nH2)
      if (.not. (W_dil .gt. 0.0d0)) return
      lam = lam - nH3p*fourpi*h3p_emission_lte(T_rad)                     &
                             *h3p_nonlte_factor(T,nH2)*W_dil*1.0d7
      end function h3p_net_cooling_rate

      ! End of module
      end module h3p_cooling
