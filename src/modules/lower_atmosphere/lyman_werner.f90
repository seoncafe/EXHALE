      module lyman_werner_photodissociation
      ! Photodissociation of H2 in the Lyman and Werner bands (the Solomon
      ! process):
      !
      !     H2(X) + hv(11.2-13.6 eV)  ->  H2(B,C)  ->  H2(X, continuum)
      !                               ->  H + H .
      !
      ! This is the process that destroys H2 above the level where it can
      ! shield itself, and it is the H2 sink the Koskinen et al. (2022)
      ! Table-1 network (mol_rates) leaves out -- their photo-rates start at
      ! the 15.4 eV photoionization edge, so without this term the only H2
      ! losses below 15.4 eV are thermal (R12), electron impact (R14) and
      ! ion chemistry.  Absent it, EXHALE's network keeps a base far more
      ! molecular than photochemical models of the same layer give
      ! (docs/base_composition_handoff_plan.md sec. 5).
      !
      ! ---------------------------------------------------------------
      ! 1. Band and unattenuated rate
      !
      ! Draine & Bertoldi (1996), ApJ 468, 269 (hereafter DB96) work with
      ! the 912-1110 A interval: 912 A is the H Lyman edge (shortward of it
      ! atomic H absorbs everything), and essentially all H2 pumping out of
      ! v = 0 happens longward of 1110 A only through very weak lines
      ! (their footnote 4).  We use the same interval.
      !
      ! DB96 characterize a radiation field by the photon flux in that band
      ! (their eq. [21], F = c n_phot for a beam) and tabulate F and the
      ! unshielded dissociation rate for several spectra.  For the flat-F_lam
      ! spectrum (u_nu ~ nu^-2, their eq. [24], T_color = 2.9e4 K) at chi = 1:
      !
      !     F = 1.208e7 photons cm^-2 s^-1        (their Table 1)
      !     zeta_pump = 3.09e-10 s^-1, <p_diss> = 0.135   (their Table 2)
      !     zeta_diss(0) = 4.17e-11 s^-1          (their Fig. 7 caption)
      !
      ! so the dissociation rate per band photon is an effective cross
      ! section
      !
      !     sigma_LW = 4.17e-11 / 1.208e7 = 3.452e-18 cm^2 .
      !
      ! The code takes the band-integrated ENERGY flux at the planet,
      ! F_LW_star [erg cm^-2 s^-1], and converts with the mean photon energy
      ! of a flat F_lam band, <hv> = 2hc/(912 A + 1110 A) = 12.2635 eV, so
      !
      !     k_LW,thin = sigma_LW * F_LW_star / <hv> .
      !
      ! APPROXIMATION AND ITS RANGE.  sigma_LW depends on the shape of the
      ! spectrum WITHIN the band, because the Lyman/Werner lines sample it
      ! unevenly.  Repeating the arithmetic above for the (much softer)
      ! Draine 1978 field of DB96 eq. (23) -- zeta_pump = 2.78e-10,
      ! <p_diss> = 0.119, F = 1.232e7 at chi = 1 -- gives 2.685e-18 cm^2,
      ! 22% lower.  The two spectra bracket color temperatures 1.3e4-2.9e4 K,
      ! and DB96 note (their sec. 3) that PDR properties are insensitive to
      ! the spectrum for T_color >~ 1e4 K.  The flat-F_lam value adopted here
      ! is therefore good to about +-25% for a stellar FUV band that is not
      ! strongly tilted; a spectrum concentrated at one end of the band
      ! (a single dominant emission line) is outside that range.
      !
      ! ---------------------------------------------------------------
      ! 2. Self-shielding
      !
      ! The pumping lines saturate, so the rate falls far faster than any
      ! continuum opacity would give.  DB96 eq. (37) is their fit to the
      ! full multiline calculation including line overlap:
      !
      !   f_shield(N_H2) = 0.965/(1 + x/b5)^2
      !                  + 0.035/(1+x)^0.5 * exp[-8.5e-4 (1+x)^0.5] ,
      !   x = N_H2/5e14 cm^-2 ,   b5 = b/1e5 cm s^-1 ,
      !
      ! and their eq. (40) is zeta_diss = f_shield * exp(-tau_dust) *
      ! zeta_diss(0).  It reproduces the exact self-shielding function over
      ! 1e14 < N_H2 < 3e21 cm^-2 (their Figs. 1-5), which covers the whole
      ! molecular layer of a hot-Neptune/hot-Jupiter base.
      !
      ! b is the H2 Doppler parameter, FWHM/(4 ln 2)^(1/2); we take it
      ! thermal, b = (2kT/m_H2)^(1/2), with no turbulent term.
      !
      ! WHAT IS NOT SHIELDING THE BAND HERE:
      !  - Dust: EXHALE's metals are atomic and trace, so there are no
      !    grains.  The exp(-tau_dust) of DB96 eq. (40) is identically 1.
      !  - Trace-metal continuum: the neutral low-IP metals (Mg I, Fe I,
      !    Si I, Ca I, Na I, K I) do photoionize inside the band, but at
      !    solar abundance and sigma ~ 1e-18 cm^2 their optical depth is
      !    ~4e-23 N_H, i.e. <~ 0.05 at the base column where f_shield is
      !    already below 1e-4.  Neglected.
      !  - H Lyman-series lines (Ly-beta 1025.7, Ly-gamma 972.5, ... all lie
      !    in the band): DB96 include them in the equivalent width their fit
      !    was built on and state (their sec. 4.2) that "absorption by the H
      !    Lyman lines has only a small effect on the H2 pumping rates".
      !    They are not treated separately here.  In a wind the H I / H2
      !    column ratio is much larger than in the PDRs DB96 fitted, so this
      !    is the least controlled of the three; it can only reduce the rate.
      !
      ! ---------------------------------------------------------------
      ! 3. Heating
      !
      ! The fragments of a Solomon-process dissociation leave with kinetic
      ! energy.  Black & Dalgarno (1977), ApJS 34, 405, p. 418: "Fluorescent
      ! dissociation of H2 gives rise to a pair of energetic hydrogen atoms
      ! (Milgrom, Panagia, and Salpeter 1973; Stephens and Dalgarno 1973);
      ! for a typical ultraviolet radiation field, the yield is about 0.4 eV
      ! per atom pair, but it varies slightly with depth."  We adopt the
      ! 0.4 eV, held constant with depth.
      !
      ! The 4.48 eV H-H bond energy is paid by the absorbed photon, not by
      ! the gas, so it is NOT a thermal sink of this channel.  (A thermal
      ! dissociation-energy sink for R12/R14 is a separate open item.)
      !
      ! ---------------------------------------------------------------
      ! Reference: Draine & Bertoldi (1996) ApJ 468, 269 (published
      ! version, eqs. 20, 21, 24, 37, 40, Tables 1-2, Fig. 7 caption);
      ! Black & Dalgarno (1977) ApJS 34, 405, p. 418.

      implicit none
      private
      public :: h2_doppler_parameter, h2_self_shielding_factor,          &
                lyman_werner_dissociation_rate, e_lw_fragment_erg

      ! Effective H2 dissociation cross section per band photon [cm^2],
      ! DB96 flat-F_lam (u_nu ~ nu^-2) calibration; see sec. 1 above.
      real*8, parameter :: sigma_lw = 3.4520d-18

      ! Mean photon energy of a flat-F_lam 912-1110 A band [erg]:
      ! 2hc/(912+1110 A) = 12.2635 eV.
      real*8, parameter :: e_lw_photon_erg = 1.96483d-11

      ! Unattenuated dissociation rate per unit band energy flux
      ! [s^-1 / (erg cm^-2 s^-1)] = [cm^2 erg^-1].
      real*8, parameter :: k_lw_per_flux = sigma_lw/e_lw_photon_erg

      ! Kinetic energy released to the H + H pair, 0.4 eV in erg
      ! (Black & Dalgarno 1977, p. 418).
      real*8, parameter :: e_lw_fragment_erg = 6.40871d-13

      ! Boltzmann constant and the H2 mass, in the same cgs values the rest
      ! of the code uses (global_parameters kb_erg, mu); repeated here so
      ! the module has no dependency on the parameter block.
      real*8, parameter :: kb_lw = 1.38d-16
      real*8, parameter :: m_h2  = 2.0d0*1.673d-24

      contains

      ! Thermal Doppler parameter of H2, b = (2 k T / m_H2)^(1/2) [cm s^-1].
      ! No turbulent contribution: the model has no sub-grid velocity field,
      ! and b enters f_shield only through the saturated core, where a larger
      ! b would raise the rate.
      double precision function h2_doppler_parameter(T) result(b)
      real*8, intent(in) :: T
      b = sqrt(2.0d0*kb_lw*max(T, 1.0d0)/m_h2)
      end function h2_doppler_parameter

      ! DB96 eq. (37): H2 self-shielding factor for a star-ward H2 column
      ! N_H2 [cm^-2] and Doppler parameter b [cm s^-1].  Valid (fitted) over
      ! 1e14 < N_H2 < 3e21 cm^-2; it tends to 1 as N_H2 -> 0 and falls as
      ! N_H2^-3/4 in the saturated regime, steeper than N_H2^-1/2 because of
      ! line overlap.
      double precision function h2_self_shielding_factor(N_H2, b)         &
                                result(f_shield)
      real*8, intent(in) :: N_H2, b
      real*8 :: x, b5, s
      x  = max(N_H2, 0.0d0)/5.0d14
      b5 = max(b, 1.0d0)/1.0d5
      s  = sqrt(1.0d0 + x)
      f_shield = 0.965d0/(1.0d0 + x/b5)**2                                &
               + 0.035d0/s*exp(-8.5d-4*s)
      end function h2_self_shielding_factor

      ! Photodissociation rate of H2 [s^-1] for a band-integrated stellar
      ! energy flux F_LW [erg cm^-2 s^-1] at the planet, a star-ward H2
      ! column N_H2 [cm^-2] and a gas temperature T [K].
      ! DB96 eq. (40) with no dust term.
      double precision function lyman_werner_dissociation_rate(F_LW,     &
                                N_H2, T) result(k)
      real*8, intent(in) :: F_LW, N_H2, T
      if (F_LW .le. 0.0d0) then
         k = 0.0d0
         return
      endif
      k = k_lw_per_flux*F_LW                                             &
        * h2_self_shielding_factor(N_H2, h2_doppler_parameter(T))
      end function lyman_werner_dissociation_rate

      ! End of module
      end module lyman_werner_photodissociation
