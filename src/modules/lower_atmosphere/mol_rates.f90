      module mol_rates
      ! molecular (H2 / H2+ / H3+ / HeH+) reaction-rate coefficients,
      ! transcribed VERBATIM from Koskinen et al. (2022), ApJ 929:52, Table 1
      ! (page 19; verified against the PDF in references/).  All rates in
      ! cgs (cm^3 s^-1; the two three-body rates R13/R15 in cm^6 s^-1 --
      ! multiply by the total density n as in the table's "n" factor).
      ! T = heavy-particle temperature [K]; Te = electron temperature [K]
      ! (EXHALE uses a common T).
      !
      ! Photo-rates P1-P5 (H, He, H2 photoionization; H2 dissociative and
      ! double photoionization) are "SC" in the paper -- computed from cross
      ! sections x stellar flux x column densities.  In EXHALE these follow
      ! the existing PH_heat routines (util_ion_eq) once H2 cross sections
      ! are added; they are NOT part of this module.
      !
      ! NOTE (Koskinen 2022 baseline): neutral H2 photodissociation through
      ! the Lyman-Werner bands is NOT part of their Table 1; their
      ! sensitivity test with a Backx et al. (1976) cross section and
      ! dissociation probability 0.125 changed Mdot by <= 1.4x.  EXHALE adds
      ! it separately and opt-in, in
      ! src/modules/lower_atmosphere/lyman_werner.f90 (Draine & Bertoldi
      ! 1996; key "Stellar LW flux"), so it is deliberately absent here.

      implicit none
      private
      public :: rk_R1_Hp_rec,   rk_R2_Hep_rec,  rk_R3_H_cion,            &
                rk_R4_He_cion,  rk_R5_H2p_dr,   rk_R6_H3p_dr_H2,         &
                rk_R7_H3p_dr_3H, rk_R8_H2p_H2,  rk_R9_H2p_H,             &
                rk_R10_Hp_H2v4, rk_R11_H3p_H,   rk_R12_H2_thdis,         &
                rk_R13_Hp_H2_M, rk_R14_H2_edis, rk_R15_3body_H2,         &
                rk_R16_HeHp_dr, rk_R17_Hep_H2_diss, rk_R18_HeHp_H2,      &
                rk_R19_HeHp_H,  rk_R20_Hep_H2_HeHp, rk_R21_H_Hep_cx,     &
                rk_R22_Hp_He_cx, rk_R23_H2_Hep_cx

      contains

      ! R1: H+ + e -> H + hv                    (Storey & Hummer 1995)
      double precision function rk_R1_Hp_rec(Te) result(k)
      real*8, intent(in) :: Te
      k = 4.0d-12*(300.0d0/Te)**0.64d0
      end function

      ! R2: He+ + e -> He + hv                  (Storey & Hummer 1995)
      double precision function rk_R2_Hep_rec(Te) result(k)
      real*8, intent(in) :: Te
      k = 4.6d-12*(300.0d0/Te)**0.64d0
      end function

      ! R3: H + e -> H+ + 2e                    (Voronov 1997)
      ! U = 13.6 eV / E_e;  E_e = kB Te in eV.
      double precision function rk_R3_H_cion(Te) result(k)
      real*8, intent(in) :: Te
      real*8 :: U
      U = 13.6d0/(8.6173d-5*Te)
      k = 2.91d-8*U**0.39d0*exp(-U)/(0.232d0 + U)
      end function

      ! R4: He + e -> He+ + 2e                  (Voronov 1997)
      double precision function rk_R4_He_cion(Te) result(k)
      real*8, intent(in) :: Te
      real*8 :: U
      U = 24.6d0/(8.6173d-5*Te)
      k = 1.75d-8*U**0.35d0*exp(-U)/(0.180d0 + U)
      end function

      ! R5: H2+ + e -> H + H                    (Auerbach et al. 1977)
      double precision function rk_R5_H2p_dr(Te) result(k)
      real*8, intent(in) :: Te
      k = 2.3d-8*(300.0d0/Te)**0.4d0
      end function

      ! R6: H3+ + e -> H2 + H                   (Larsson et al. 2008)
      double precision function rk_R6_H3p_dr_H2(Te) result(k)
      real*8, intent(in) :: Te
      k = 2.16d-8*(300.0d0/Te)**0.65d0
      end function

      ! R7: H3+ + e -> H + H + H                (Larsson et al. 2008)
      double precision function rk_R7_H3p_dr_3H(Te) result(k)
      real*8, intent(in) :: Te
      k = 5.04d-8*(300.0d0/Te)**0.65d0
      end function

      ! R8: H2+ + H2 -> H3+ + H                 (Theard & Huntress 1974)
      double precision function rk_R8_H2p_H2() result(k)
      k = 2.0d-9
      end function

      ! R9: H2+ + H -> H+ + H2                  (Karpas et al. 1979)
      double precision function rk_R9_H2p_H() result(k)
      k = 6.4d-10
      end function

      ! R10: H+ + H2(v>=4) -> H2+ + H           (Yelle 2004)
      double precision function rk_R10_Hp_H2v4(T) result(k)
      real*8, intent(in) :: T
      k = 1.0d-9*exp(-21900.0d0/T)
      end function

      ! R11: H3+ + H -> H2+ + H2                (Harada et al. 2010)
      double precision function rk_R11_H3p_H(T) result(k)
      real*8, intent(in) :: T
      k = 2.1d-9*exp(-20000.0d0/T)
      end function

      ! R12: H2 + M -> H + H + M                (Baulch et al. 1992)
      double precision function rk_R12_H2_thdis(T) result(k)
      real*8, intent(in) :: T
      k = 1.5d-9*exp(-48350.0d0/T)
      end function

      ! R13: H+ + H2 + M -> H3+ + M             (Miller et al. 1968)
      ! Returns the 2-body-equivalent rate: 3.2e-29 * n  [cm^3 s^-1].
      double precision function rk_R13_Hp_H2_M(n) result(k)
      real*8, intent(in) :: n
      k = 3.2d-29*n
      end function

      ! R14: H2 + e -> H + H + e                (Stibbe & Tennyson 1999)
      double precision function rk_R14_H2_edis(Te) result(k)
      real*8, intent(in) :: Te
      k = 1.33d-6*(300.0d0/Te)**0.91d0*exp(-55800.0d0/Te)
      end function

      ! R15: H + H + M -> H2 + M                (Ham et al. 1970)
      ! Returns the 2-body-equivalent rate: 8e-33 (300/T)^0.6 * n.
      double precision function rk_R15_3body_H2(T, n) result(k)
      real*8, intent(in) :: T, n
      k = 8.0d-33*(300.0d0/T)**0.6d0*n
      end function

      ! R16: HeH+ + e -> He + H                 (Yousif & Mitchell 1989)
      double precision function rk_R16_HeHp_dr(Te) result(k)
      real*8, intent(in) :: Te
      k = 1.0d-8*(300.0d0/Te)**0.6d0
      end function

      ! R17: He+ + H2 -> H+ + H + He            (Moses & Bass 2000)
      double precision function rk_R17_Hep_H2_diss(T) result(k)
      real*8, intent(in) :: T
      k = 1.0d-9*exp(-5700.0d0/T)
      end function

      ! R18: HeH+ + H2 -> H3+ + He              (Bohme et al. 1980)
      double precision function rk_R18_HeHp_H2() result(k)
      k = 1.5d-9
      end function

      ! R19: HeH+ + H -> H2+ + He               (Karpas et al. 1979)
      double precision function rk_R19_HeHp_H() result(k)
      k = 9.1d-10
      end function

      ! R20: He+ + H2 -> HeH+ + H               (Schauer et al. 1989)
      double precision function rk_R20_Hep_H2_HeHp() result(k)
      k = 4.2d-13
      end function

      ! R21: H + He+ -> H+ + He                 (Stancil et al. 1998)
      double precision function rk_R21_H_Hep_cx(T) result(k)
      real*8, intent(in) :: T
      k = 1.2d-15*(300.0d0/T)**(-0.25d0)
      end function

      ! R22: H+ + He -> H + He+                 (Glover & Jappsen 2007)
      double precision function rk_R22_Hp_He_cx(T) result(k)
      real*8, intent(in) :: T
      k = 1.75d-11*(300.0d0/T)**0.75d0*exp(-128000.0d0/T)
      end function

      ! R23: H2 + He+ -> H2+ + He               (Barlow 1984)
      double precision function rk_R23_H2_Hep_cx() result(k)
      k = 7.2d-15
      end function

      ! End of module
      end module mol_rates
