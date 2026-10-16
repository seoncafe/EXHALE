!==============================================================================
! molecular_infrared_cooling
!
! Net infrared exchange of the molecular layer with the atmosphere below it:
! the H2O and CO vibration-rotation bands and the H2 quadrupole plus magnetic
! dipole line spectrum, each emitting in LTE and each absorbing the diluted
! blackbody the unresolved lower atmosphere presents.
!
! WHY THIS EXISTS.  A converged Tier-2 molecular layer radiates itself down to
! 190-400 K against an equilibrium temperature of 1100-1400 K, because the only
! infrared coolants the code carried below the H2 -> H front were H3+ and the
! ground-term fine-structure lines, treated as pure emitters into vacuum.  The
! `Base IR field` closure (Update_EXHALE_stage1.pdf section 55) gave those two families
! their incident field; it could not reach the rest of the layer, where the
! coolants that matter in a real H2 atmosphere -- water, carbon monoxide and
! molecular hydrogen itself -- were simply absent.  This module supplies them.
!
! THE PHYSICS.  For an optically thin layer in LTE the net radiative loss per
! unit volume of species X is
!
!     Lambda_net = n_X [ E_X(T) - W E_X^abs(T; T_rad) ]                    (1)
!
! with, for a band absorber,
!
!     E_X(T)            = 4 pi Sum_b sigma_b(T) Int_b B_nu(T)     dnu      (2)
!     E_X^abs(T; T_rad) = 4 pi Sum_b sigma_b(T) Int_b B_nu(T_rad) dnu      (3)
!
! and, for the H2 line sum,
!
!     E_H2(T)     = Sum_l f_u(T) A_l dE_l                                  (4)
!     E_H2^abs    = Sum_l f_u(T) A_l dE_l nbar_l [exp(dE_l/T) - 1]         (5)
!     nbar_l      = 1 / [exp(dE_l/T_rad) - 1]                              (6)
!
! where f_u(T) = g_u exp(-T_u/T) / Q(T) is the LTE population of the upper
! level and W is the geometric dilution of the incident field, the fraction
! of the sky the lower atmosphere subtends seen from radius r,
! W = [1 - sqrt(1 - (R_p/r)^2)]/2 (0.5*base_sky_fraction(r, 1), util_ion_eq),
! and zero with `Base IR field` off.  W carries no screening by the
! intervening column: the layer is treated as optically thin (VALIDITY).
! Q is the internal partition function of H2 of caloric_eos, the one
! Boltzmann sum the code forms over the rovibrational ladder, so the
! populations that radiate here are the populations whose energy the equation
! of state carries and whose free energy sets the H + H <-> H2 equilibrium of
! mol_rates.
!
! Equation (5) is the absorption MINUS the stimulated emission: the full net
! factor of a line is [1 + nbar - nbar exp(dE/T)].  Keeping the stimulated term
! is what makes (1) vanish identically at T = T_rad, W = 1 -- a layer buried in
! a blackbody at its own temperature neither cools nor heats.  Equations (2)
! and (3) have that property automatically, because the HITRAN line intensities
! the cross sections are built from already carry the stimulated-emission
! correction, so sigma_b is the coefficient that pairs with B_nu.
!
! It is that fixed point, not the magnitude of the rate, that answers item (G):
! each channel stops cooling at its own radiative equilibrium temperature
! instead of running the layer down to zero.  Because W multiplies only the
! absorption, the equilibrium temperature is set by the ratio E/E^abs and is
! therefore insensitive, to first order, to a grey escape probability applied
! to both -- which is why the optically thin limit is enough to fix the
! temperature even where the bands are marginally thick.
!
! VALIDITY.
!  * Optically thin, and the run measures whether it still is.  No escape
!    probability is applied.  The Planck-mean optical depth of the H2O and CO
!    columns from each cell to the top of the domain is written to
!    output/Cooling_breakdown.txt: on the 12000-step HD 189733 b
!    oxygen-chemistry run the maxima are 2.5e-3 and 2.8e-3, and the
!    column-integrated absorbed power reported in output/FUV_bands.txt,
!    7.5e5 erg cm^-2 s^-1, sits 73x below the incident infrared flux
!    W sigma T0^4 of 5.4e7 that bounds it.  The H2 lines are
!    thinner still: the line-centre opacity of 0-0 S(1) at 17 um gives
!    tau = 3e-4 over a scale height at n(H2) = 5e13 cm^-3 and 1000 K.  A
!    configuration with a much deeper base has not been checked; if the
!    reported depths approach unity the closure has left its range.  Note also
!    that a grey escape probability would multiply BOTH terms of eq. (1) and
!    so would barely move the equilibrium temperature, only the rate.
!  * LTE level populations.  The rotational levels of all three species
!    thermalize far below the densities of this layer.  The H2 and H2O
!    VIBRATIONAL levels are the constraint, and for H2 the margin is large:
!    the critical density of v = 1, n_crit = A_10/gamma_10 (Hollenbach &
!    McKee 1979, ApJS 41, 555, p. 576), is 2e6 to 6e6 cm^-3 over 900-1350 K
!    with H2 as the collider and 5e4 to 8e4 cm^-3 with atomic H, from their
!    A_10 = 8.3e-7 s^-1 (p. 583, J-averaged, Turner, Kirby-Docken & Dalgarno
!    1977) and their eq. (6.29) rate coefficients.  The base carries
!    1e13-1e14 cm^-3, seven to nine decades above it.  (Neither Hollenbach &
!    McKee nor Burton, Hollenbach & Tielens 1990 gives a HELIUM collision
!    rate for H2 v = 1, so a helium-dominated gas is outside both; leaving
!    helium out of the collider sum understates the de-excitation rate,
!    which errs on the side of NOT claiming LTE.)  Above the base the
!    density falls,
!    but so does the temperature, and the vibrational bands are a small share
!    of the H2 emission there: 0.7% at 400 K, 5.7% at 500 K, 18% at 600 K and
!    66% at 1000 K.  The LTE overestimate is therefore confined to the warm,
!    dense part of the layer, where LTE is also best justified.
!  * Cross-section temperature range.  The H2O and CO tables run 50-2000 K and
!    are clamped outside it.  Both molecules are largely dissociated above
!    2000 K, and the oxygen option caps its transported CO at the chemical
!    equilibrium of the local (n,T), so the clamped region carries little
!    density.  The H2 line sum has no such limit in the ladder; it is
!    tabulated on the 30 K to about 30100 K grid below (t_tab_hi) and the
!    edge value is returned above that, where no H2 survives.
!  * The lower atmosphere is taken to be black at these wavelengths and to
!    radiate B_nu(T_rad).  H2-H2 and H2-He collision-induced absorption is what
!    makes it black in the windows between the bands, and it is deliberately
!    NOT a channel here: measured on the petitRADTRANS/HITRAN CIA tables, its
!    optical depth per pressure scale height at a 1 microbar base is 2e-9
!    (HD 189733 b) to 8e-9 (hot Uranus), and, scaling as n^2, reaches unity only
!    near 0.2-0.4 bar.  It is a property of the reservoir below the base, not a
!    local coolant of this domain.
!  * The incident field is not attenuated by the intervening molecular column
!    of the domain itself, in the same approximation `Base IR field` makes for
!    H3+.
!
! SOURCES.  H2O and CO band-mean cross sections: correlated-k coefficients
! distributed with Photochem (Wogan et al. 2025, PSJ 6, 256), computed with
! HELIOS-K (Grimm & Heng 2015, ApJ 808, 182) from HITEMP 2010 (H2O; Rothman et
! al. 2010, JQSRT 111, 2139) and HITEMP 2019 (CO; Li et al. 2015, ApJS 216, 15).
! H2 line list: Roueff et al. (2019), A&A 630, A58, table 2 -- the transition
! wavenumbers and the quadrupole plus magnetic dipole Einstein A coefficients
! of that table, over the level energies of the same table's ladder, so the
! lines and the partition function that normalizes them are one level set.
! The tables themselves are in molecular_infrared_data.f90; their provenance,
! validity ranges and cross-checks are in
! cooling_data/molecular_infrared_bands.py.
!==============================================================================
      module molecular_infrared_cooling

      use molecular_infrared_data
      ! The physical constants have exactly one owner, global_parameters; do
      ! not restate them here.
      use global_parameters, only: kb_erg, hp_erg, c_light
      ! The code forms one Boltzmann sum over the H2 rovibrational ladder,
      ! in caloric_eos; the ladder itself lives in molecular_infrared_data,
      ! which this module already reads for the line list built on the same
      ! levels.
      use caloric_eos, only: h2_partition_function

      implicit none
      private

      public :: molecular_infrared_init
      public :: molecular_infrared_ready
      public :: h2o_band_emission_lte,  h2o_band_net_cooling_rate
      public :: co_band_emission_lte,   co_band_net_cooling_rate
      public :: h2_line_emission_lte,   h2_line_net_cooling_rate
      public :: molecular_planck_cross_section

      real*8, parameter :: fourpi = 12.566370614359172d0

      ! Evaluation grid of the derived emission / absorption functions.  It is
      ! log-spaced and wider than the cross-section tables so that the H2 line
      ! sum, which has no upper temperature limit, keeps its own range.  The
      ! spacing is that of 401 nodes over 30-8000 K (one node per 1.40
      ! percent in T); the grid continues at that spacing to t_tab_hi, above
      ! every temperature a layer that still holds H2 reaches, so that the
      ! clamp of tab_value (edge value outside the grid) never binds on the H2
      ! line sum.  MEASURED before the extension: the clamp at 8000 K gave
      ! 1.181e-18 against the direct sum's 1.579e-18 erg s^-1 at 10000 K.
      integer, parameter :: n_tab_ref = 401
      real*8,  parameter :: t_tab_lo  = 3.0d1
      real*8,  parameter :: t_tab_ref = 8.0d3
      integer, parameter :: n_tab     = n_tab_ref + 95
      real*8,  parameter :: t_tab_hi  = t_tab_lo*(t_tab_ref/t_tab_lo)    &
                                        **(dble(n_tab - 1)/dble(n_tab_ref - 1))

      ! Simpson nodes per band for the Planck integrals (must be odd).
      integer, parameter :: n_quad = 25

      real*8  :: t_tab(n_tab)
      real*8  :: emis_h2o(n_tab), absr_h2o(n_tab), sigp_h2o(n_tab)
      real*8  :: emis_co(n_tab),  absr_co(n_tab),  sigp_co(n_tab)
      real*8  :: emis_h2(n_tab),  absr_h2(n_tab)
      real*8  :: t_rad_built = -1.0d0
      logical :: tab_built   = .false.

      ! Floor used when a rate is carried through a logarithm.  Nothing in the
      ! layer can radiate below this and it keeps log(0) out of the tables.
      real*8, parameter :: rate_floor = 1.0d-300

      contains

!------------------------------------------------------------------------------
! Build the emission and absorption functions for one radiating temperature.
! Called once per run, before any cell is evaluated: the tables are read
! concurrently by the OpenMP cell sweep and must not be built inside it.
!------------------------------------------------------------------------------
      subroutine molecular_infrared_init(T_rad)
      implicit none
      real*8, intent(in) :: T_rad
      integer :: k, b, l
      real*8  :: tg, tr, z, fu, ex, nbar, arg, pwr
      real*8  :: bg_h2o(nb_h2o), br_h2o(nb_h2o)
      real*8  :: bg_co(nb_co),   br_co(nb_co)
      real*8  :: sg_h2o(nb_h2o), sg_co(nb_co)
      real*8  :: bsum

      tr = max(T_rad, 1.0d0)
      if (tab_built .and. tr .eq. t_rad_built) return

      do k = 1, n_tab
         t_tab(k) = t_tab_lo*(t_tab_ref/t_tab_lo)                         &
                             **(dble(k - 1)/dble(n_tab_ref - 1))
      enddo

      ! Incident field, one value per band: it does not depend on the gas
      ! temperature, so it is built once.
      do b = 1, nb_h2o
         br_h2o(b) = planck_band_integral(h2o_wl_lo(b), h2o_wl_hi(b), tr)
      enddo
      do b = 1, nb_co
         br_co(b)  = planck_band_integral(co_wl_lo(b),  co_wl_hi(b),  tr)
      enddo

      do k = 1, n_tab
         tg = t_tab(k)

         ! --- H2O and CO bands -------------------------------------------
         call band_cross_sections(tg, sg_h2o, sg_co)
         bsum = 0.0d0
         do b = 1, nb_h2o
            bg_h2o(b) = planck_band_integral(h2o_wl_lo(b), h2o_wl_hi(b), tg)
            bsum      = bsum + bg_h2o(b)
         enddo
         emis_h2o(k) = fourpi*dot_product(sg_h2o, bg_h2o)
         absr_h2o(k) = fourpi*dot_product(sg_h2o, br_h2o)
         sigp_h2o(k) = dot_product(sg_h2o, bg_h2o)/max(bsum, 1.0d-300)

         bsum = 0.0d0
         do b = 1, nb_co
            bg_co(b) = planck_band_integral(co_wl_lo(b), co_wl_hi(b), tg)
            bsum     = bsum + bg_co(b)
         enddo
         emis_co(k) = fourpi*dot_product(sg_co, bg_co)
         absr_co(k) = fourpi*dot_product(sg_co, br_co)
         sigp_co(k) = dot_product(sg_co, bg_co)/max(bsum, 1.0d-300)

         ! --- H2 lines ----------------------------------------------------
         ! The normalization of the line populations is the internal
         ! partition function of H2 held by caloric_eos: the levels these
         ! lines connect are the levels whose energy the equation of state
         ! carries and whose free energy the H + H <-> H2 equilibrium of
         ! mol_rates is built from, so f_u below is a population of the
         ! EOS ladder and not of a second one.
         z = h2_partition_function(tg)
         emis_h2(k) = 0.0d0
         absr_h2(k) = 0.0d0
         do l = 1, n_h2_line
            ! Upper level unpopulated: the line contributes to neither term.
            if (h2_line_Tu(l)/tg .gt. 7.0d2) cycle
            fu   = h2_line_gu(l)*exp(-h2_line_Tu(l)/tg)/z
            pwr  = fu*h2_line_A(l)*h2_line_dE(l)*kB_erg
            emis_h2(k) = emis_h2(k) + pwr
            ! The absorption is guarded separately: a line the incident field
            ! carries no photons at still EMITS, so this guard must not skip
            ! the term above with it.
            arg = h2_line_dE(l)/tr
            if (arg .gt. 7.0d2) cycle
            nbar = 1.0d0/(exp(arg) - 1.0d0)
            ex   = exp(min(h2_line_dE(l)/tg, 7.0d2))
            absr_h2(k) = absr_h2(k) + pwr*nbar*(ex - 1.0d0)
         enddo
      enddo

      t_rad_built = tr
      tab_built   = .true.

      end subroutine molecular_infrared_init

!------------------------------------------------------------------------------
      logical function molecular_infrared_ready()
      implicit none
      molecular_infrared_ready = tab_built
      end function molecular_infrared_ready

!------------------------------------------------------------------------------
! Band-mean cross sections at one gas temperature, interpolated in log T on the
! tabulated grid and clamped outside it (see the validity note in the header).
!------------------------------------------------------------------------------
      subroutine band_cross_sections(T, sg_h2o, sg_co)
      implicit none
      real*8, intent(in)  :: T
      real*8, intent(out) :: sg_h2o(nb_h2o), sg_co(nb_co)
      integer :: b, j
      real*8  :: w

      call locate_log(mir_T, n_mir_T, T, j, w)
      do b = 1, nb_h2o
         sg_h2o(b) = (1.0d0 - w)*h2o_sig(b,j) + w*h2o_sig(b,j+1)
      enddo
      do b = 1, nb_co
         sg_co(b)  = (1.0d0 - w)*co_sig(b,j)  + w*co_sig(b,j+1)
      enddo

      end subroutine band_cross_sections

!------------------------------------------------------------------------------
! Bracket x in an ascending grid and return the linear weight in log x, clamped
! at both ends.
!------------------------------------------------------------------------------
      subroutine locate_log(grid, n, x, j, w)
      implicit none
      integer, intent(in)  :: n
      real*8,  intent(in)  :: grid(n), x
      integer, intent(out) :: j
      real*8,  intent(out) :: w
      integer :: lo, hi, mid

      if (.not. (x .gt. grid(1))) then
         j = 1
         w = 0.0d0
         return
      endif
      if (.not. (x .lt. grid(n))) then
         j = n - 1
         w = 1.0d0
         return
      endif
      lo = 1
      hi = n
      do while (hi - lo .gt. 1)
         mid = (lo + hi)/2
         if (x .ge. grid(mid)) then
            lo = mid
         else
            hi = mid
         endif
      enddo
      j = lo
      w = log(x/grid(j))/log(grid(j+1)/grid(j))

      end subroutine locate_log

!------------------------------------------------------------------------------
! Int B_nu(T) dnu over [wl_lo, wl_hi] in micron, by Simpson in log lambda.
! [erg cm^-2 s^-1 sr^-1]
!------------------------------------------------------------------------------
      real*8 function planck_band_integral(wl_lo, wl_hi, T) result(bint)
      implicit none
      real*8, intent(in) :: wl_lo, wl_hi, T
      integer :: i
      real*8  :: rat, lam, nu, x, bnu, wgt, dlog

      bint = 0.0d0
      if (.not. (wl_hi .gt. wl_lo)) return
      dlog = log(wl_hi/wl_lo)/dble(n_quad - 1)
      do i = 1, n_quad
         if (i .eq. 1 .or. i .eq. n_quad) then
            wgt = 1.0d0
         else if (mod(i,2) .eq. 0) then
            wgt = 4.0d0
         else
            wgt = 2.0d0
         endif
         lam = wl_lo*exp(dlog*dble(i - 1))*1.0d-4          ! micron -> cm
         nu  = c_light/lam
         x   = hp_erg*nu/(kB_erg*max(T, 1.0d-3))
         if (x .lt. 6.0d2) then
            bnu = 2.0d0*hp_erg*nu**3/c_light**2/(exp(x) - 1.0d0)
         else
            bnu = 0.0d0
         endif
         ! dnu = -nu dln(lambda), so the magnitude of the integrand in
         ! ln(lambda) is B_nu * nu.
         bint = bint + wgt*bnu*nu
      enddo
      bint = bint*dlog/3.0d0

      end function planck_band_integral

!------------------------------------------------------------------------------
! Interpolate one of the derived functions at a gas temperature.  The functions
! span many decades, so the interpolation is log-log.
!------------------------------------------------------------------------------
      real*8 function tab_value(tabv, T) result(y)
      implicit none
      real*8, intent(in) :: tabv(n_tab), T
      integer :: j
      real*8  :: w, y1, y2

      call locate_log(t_tab, n_tab, T, j, w)
      y1 = max(tabv(j),   rate_floor)
      y2 = max(tabv(j+1), rate_floor)
      y  = exp((1.0d0 - w)*log(y1) + w*log(y2))
      if (y .le. rate_floor) y = 0.0d0

      end function tab_value

!------------------------------------------------------------------------------
! Optically thin LTE emission per molecule [erg s^-1].
!------------------------------------------------------------------------------
      real*8 function h2o_band_emission_lte(T) result(e)
      implicit none
      real*8, intent(in) :: T
      e = 0.0d0
      if (.not. tab_built) return
      e = tab_value(emis_h2o, T)
      end function h2o_band_emission_lte

      real*8 function co_band_emission_lte(T) result(e)
      implicit none
      real*8, intent(in) :: T
      e = 0.0d0
      if (.not. tab_built) return
      e = tab_value(emis_co, T)
      end function co_band_emission_lte

      real*8 function h2_line_emission_lte(T) result(e)
      implicit none
      real*8, intent(in) :: T
      e = 0.0d0
      if (.not. tab_built) return
      e = tab_value(emis_h2, T)
      end function h2_line_emission_lte

!------------------------------------------------------------------------------
! Planck-mean cross section [cm^2/molecule], for the optical-depth diagnostic.
! ispec = 1 water, 2 carbon monoxide.
!------------------------------------------------------------------------------
      real*8 function molecular_planck_cross_section(ispec, T) result(s)
      implicit none
      integer, intent(in) :: ispec
      real*8,  intent(in) :: T
      s = 0.0d0
      if (.not. tab_built) return
      if (ispec .eq. 1) then
         s = tab_value(sigp_h2o, T)
      else if (ispec .eq. 2) then
         s = tab_value(sigp_co, T)
      endif
      end function molecular_planck_cross_section

!------------------------------------------------------------------------------
! Net radiative loss [erg cm^-3 s^-1], equation (1).  W_dil = 0 recovers the
! optically thin emission into vacuum; W_dil = 1 with T = T_rad returns zero.
!------------------------------------------------------------------------------
      real*8 function h2o_band_net_cooling_rate(T, nH2O, W_dil) result(lam)
      implicit none
      real*8, intent(in) :: T, nH2O, W_dil
      lam = 0.0d0
      if (.not. tab_built) return
      if (.not. (nH2O .gt. 0.0d0)) return
      lam = nH2O*(tab_value(emis_h2o, T) - W_dil*tab_value(absr_h2o, T))
      end function h2o_band_net_cooling_rate

      real*8 function co_band_net_cooling_rate(T, nCO, W_dil) result(lam)
      implicit none
      real*8, intent(in) :: T, nCO, W_dil
      lam = 0.0d0
      if (.not. tab_built) return
      if (.not. (nCO .gt. 0.0d0)) return
      lam = nCO*(tab_value(emis_co, T) - W_dil*tab_value(absr_co, T))
      end function co_band_net_cooling_rate

      real*8 function h2_line_net_cooling_rate(T, nH2, W_dil) result(lam)
      implicit none
      real*8, intent(in) :: T, nH2, W_dil
      lam = 0.0d0
      if (.not. tab_built) return
      if (.not. (nH2 .gt. 0.0d0)) return
      lam = nH2*(tab_value(emis_h2, T) - W_dil*tab_value(absr_h2, T))
      end function h2_line_net_cooling_rate

      end module molecular_infrared_cooling
