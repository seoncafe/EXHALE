      program a2_m1_rate_check
      ! Standalone check of src/modules/lower_atmosphere/oxygen_rates.f90,
      ! the milestone-M1 deliverable of docs/a2_oxygen_option_design.md.
      ! It is NOT part of the EXHALE build: it is compiled by hand, e.g.
      !
      !   gfortran -O2 -o a2_m1_rate_check.x \
      !       ../../modules/lower_atmosphere/oxygen_rates.f90 \
      !       a2_m1_rate_check.f90
      !   ./a2_m1_rate_check.x > a2_m1_rate_check.out
      !
      ! run from src/tests/a2_m1/.  Its saved output is a2_m1_rate_check.out
      ! and the audit table it feeds is docs/a2_reaction_audit.md.
      !
      ! What it tests, in order:
      !   1. Each adopted rate expression against the independent
      !      transcription in the other reference network, over that
      !      network's stated temperature range.
      !   2. The Shomate table against the NIST Chemistry WebBook values at
      !      298.15 K (the WebBook's own H and S entries) and against
      !      Burcat's NASA-9 set at 298.15, 1000 and 2000 K.
      !   3. The round-trip residual of rate_from_detailed_balance
      !      (k_fwd -> k_rev -> k_fwd), which must be at machine precision.
      !   4. The physical detailed-balance test: reversing O10 with the
      !      thermodynamic data must reproduce the SEPARATELY published
      !      rate of the opposite direction, O + H2O -> OH + OH.
      !   5. The measured size of every reaction at the HD 189733 b base,
      !      from the abundances of the P1 Photochem/NCHO run, which is what
      !      the include/exclude verdicts rest on.

      use oxygen_rates
      implicit none
      integer, parameter :: dp = kind(1.0d0)
      real(dp), parameter :: kb_erg = 1.380649d-16
      real(dp) :: T
      integer  :: i, j

      ! ---- Reference-network expressions, transcribed here so that the
      !      module and the comparison cannot drift apart silently.
      ! VULCAN NCHO_photo_network.txt:
      !   id   1  OH + H2  -> H2O + H   3.57e-16 T^1.52  exp(-1740/T)
      !   id   3  O  + H2  -> OH  + H   8.52e-20 T^2.67  exp(-3160/T)
      !   id   5  O  + H2O -> OH  + OH  8.20e-14 T^0.95  exp(-8570/T)
      !   id 615  O(1D)+H2 -> OH  + H   2.87e-10
      !   id 619  O(1D)+H2O-> OH  + OH  1.62e-10 exp(+65/T)
      !   id 657  H  + O + M -> OH + M  k0 = 1.30e-29 T^-1
      !   id 659  OH + H + M -> H2O+ M  k0 = 3.89e-25 T^-2
      ! Photochem zahnle_earth.yaml:
      !   O1D + H2 <=> OH + H           1.5e-10
      !   O1D + H2O <=> OH + OH         2.2e-10

      ! ---- Measured state of the HD 189733 b base, level 121 of the P1
      !      Photochem/NCHO run vulcan_work/pc_compare_p1/hd189_toa1e-2/
      !      pc_ncho_solution.pkl (P = 0.9386 dyn cm^-2, T = 863.911 K).
      real(dp), parameter :: T_base = 863.9113825339738d0
      real(dp), parameter :: n_base = 7.8690385368d12
      real(dp), parameter :: n_H    = 1.4741d-01*n_base
      real(dp), parameter :: n_H2   = 7.0132d-01*n_base
      real(dp), parameter :: n_O    = 1.0848d-07*n_base
      real(dp), parameter :: n_OH   = 2.4664d-07*n_base
      real(dp), parameter :: n_H2O  = 5.1164d-04*n_base
      real(dp), parameter :: n_O1D  = 6.0386d-11*n_base

      real(dp) :: kf, kr, kf_back, kv, kc, r1
      integer  :: react2(2), prod2(2), react3(3), prod3(3)

      write(*,'(a)') '# a2_m1_rate_check -- audit driver for oxygen_rates.f90'
      write(*,'(a)') '# All rates in cm^3 s^-1 unless marked; k0 in cm^6 s^-1.'
      write(*,*)

! ---------------------------------------------------------------------
      write(*,'(a)') '=== 0. Element and charge conservation, reaction by reaction ==='
      write(*,'(a)') '    Every species of the A2 set, its element and charge counts,'
      write(*,'(a)') '    and the balance of each reaction.  A non-zero residual is a'
      write(*,'(a)') '    failure of the M1 gate.'
      write(*,*)
      call conservation_report()
      write(*,*)

! ---------------------------------------------------------------------
      write(*,'(a)') '=== 1. Adopted expression against the other network''s transcription ==='
      write(*,*)
      write(*,'(a)') '--- O1  OH + H2 -> H2O + H'
      write(*,'(a)') '    adopted  Baulch et al. (2005) p. 1029   3.6e-16 T^1.52 exp(-1740/T), 250-2500 K'
      write(*,'(a)') '    check A  Baulch et al. (1992) p. 552    1.7e-16 T^1.6  exp(-1660/T), 300-2500 K'
      write(*,'(a)') '             (this is what zahnle_earth.yaml carries, as 1.740708e-16)'
      write(*,'(a)') '    check B  VULCAN id 1, Oldenborg & Loge  3.57e-16 T^1.52 exp(-1740/T), 250-2580 K'
      write(*,'(a)') '       T[K]        adopted        Baulch92        VULCAN     B92/ad   VUL/ad'
      do i = 1, 7
        T = tgrid(i)
        kf = rk_O1_OH_H2_water(T)
        kc = 1.7d-16*T**1.6d0*exp(-1660.0d0/T)
        kv = 3.57d-16*T**1.52d0*exp(-1740.0d0/T)
        write(*,'(f10.1,3es16.6,2f9.4)') T, kf, kc, kv, kc/kf, kv/kf
      end do
      write(*,*)

      write(*,'(a)') '--- O2  O + H2 -> OH + H'
      write(*,'(a)') '    adopted  Baulch et al. (2005) p. 804    6.34e-12 exp(-4000/T) + 1.46e-9 exp(-9650/T), 298-3300 K'
      write(*,'(a)') '    check    Baulch et al. (1992) p. 430    8.5e-20 T^2.67 exp(-3163/T), 300-2500 K'
      write(*,'(a)') '             (both networks carry this with exp(-3160/T))'
      write(*,'(a)') '       T[K]        adopted        Baulch92      networks    B92/ad   net/ad'
      do i = 1, 7
        T = tgrid(i)
        kf = rk_O2_O_H2_hydroxyl(T)
        kc = 8.5d-20*T**2.67d0*exp(-3163.0d0/T)
        kv = 8.5d-20*T**2.67d0*exp(-3160.0d0/T)
        write(*,'(f10.1,3es16.6,2f9.4)') T, kf, kc, kv, kc/kf, kv/kf
      end do
      write(*,*)

      write(*,'(a)') '--- O6  O(1D) + H2 -> OH + H     (temperature-independent)'
      write(*,'(a,es13.6)') '    adopted  Atkinson et al. (2004) I.A2.18, 200-350 K:            ', &
                            rk_O6_O1D_H2_hydroxyl()
      write(*,'(a,es13.6)') '    check A  VULCAN id 615, Tully (1975), 100-2100 K:              ', 2.87d-10
      write(*,'(a,es13.6)') '    check B  zahnle_earth.yaml, key Ba92 (untraceable, see module):', 1.5d-10
      write(*,'(a,2f8.4)')  '    A/adopted, B/adopted = ', 2.87d-10/rk_O6_O1D_H2_hydroxyl(),      &
                            1.5d-10/rk_O6_O1D_H2_hydroxyl()
      write(*,*)

      write(*,'(a)') '--- O8  O + H + M -> OH + M   (low-pressure limit k0, cm^6 s^-1)'
      write(*,'(a)') '    adopted  Tsang & Hampson (1986) p. 1111  1.3e-29 T^-1, no range, factor-10 estimate'
      write(*,'(a)') '    check    zahnle_earth.yaml                1.29e-29 T^-1'
      write(*,'(a)') '       T[K]        adopted        Photochem     PC/adopted'
      do i = 1, 7
        T = tgrid(i)
        kf = rk_O8_O_H_assoc(T)
        kv = 1.29d-29*T**(-1.0d0)
        write(*,'(f10.1,2es16.6,f14.5)') T, kf, kv, kv/kf
      end do
      write(*,*)

      write(*,'(a)') '--- O9  H + OH + M -> H2O + M (low-pressure limit k0, cm^6 s^-1)'
      write(*,'(a)') '    adopted  Baulch (1992) p. 498 / (2005) p. 913, M = N2   6.1e-26 T^-2.0, 300-3000 K'
      write(*,'(a)') '    same source, M = Ar                                     2.3e-26 T^-2.0'
      write(*,'(a)') '    same source, M = H2O                                    3.9e-25 T^-2.0'
      write(*,'(a)') '    check A  VULCAN id 659   = Baulch M = H2O               3.89e-25 T^-2.0'
      write(*,'(a)') '    check B  zahnle_earth.yaml = Javoy et al. (2003) sec 3.4, M = Ar, 2790-3200 K'
      write(*,'(a)') '       T[K]      adopted(N2)        Ar            H2O         VULCAN      Photochem   VUL/ad  PC/ad'
      do i = 1, 7
        T = tgrid(i)
        kf = rk_O9_H_OH_assoc(T)
        kc = 2.3d-26*T**(-2.0d0)
        kr = 3.9d-25*T**(-2.0d0)
        kv = 3.89d-25*T**(-2.0d0)
        kf_back = 1.050748d-26*T**(-2.1d0)
        write(*,'(f10.1,5es14.4,2f9.4)') T, kf, kc, kr, kv, kf_back,       &
             kv/kf, kf_back/kf
      end do
      write(*,*)

      write(*,'(a)') '--- O10 OH + OH -> H2O + O'
      write(*,'(a)') '    adopted  Baulch et al. (2005) p. 1032   5.56e-20 T^2.42 exp(+970/T), 250-2400 K'
      write(*,'(a)') '    check    Baulch et al. (1992) p. 555    2.5e-15  T^1.14 exp(-50/T),  250-2500 K'
      write(*,'(a)') '             (this is what zahnle_earth.yaml carries, as 2.549944e-15)'
      write(*,'(a)') '       T[K]        adopted        Baulch92     B92/adopted'
      do i = 1, 7
        T = tgrid(i)
        kf = rk_O10_OH_OH_water(T)
        kv = 2.5d-15*T**1.14d0*exp(-50.0d0/T)
        write(*,'(f10.1,2es16.6,f14.5)') T, kf, kv, kv/kf
      end do
      write(*,*)

      write(*,'(a)') '--- O12 O(1D) + H2O -> OH + OH'
      write(*,'(a)') '    adopted  Atkinson et al. (2004) I.A2.19  2.2e-10, 200-350 K'
      write(*,'(a)') '    check    VULCAN id 619                   1.62e-10 exp(+65/T), 235-370 K'
      write(*,'(a)') '       T[K]        adopted          other      other/adopted'
      do i = 1, 7
        T = tgrid(i)
        kf = rk_O12_O1D_H2O_hydroxyl()
        kv = 1.62d-10*exp(65.0d0/T)
        write(*,'(f10.1,2es16.6,f14.5)') T, kf, kv, kv/kf
      end do
      write(*,*)

      write(*,'(a)') '--- Cross-check of the EXISTING network: R15 of mol_rates.f90,'
      write(*,'(a)') '    H + H + M -> H2 + M, 8e-33 (300/T)^0.6 (Ham et al. 1970),'
      write(*,'(a)') '    against Baulch et al. (1992) Table 3 p. 428, M = H2:'
      write(*,'(a)') '    k0 = 2.7e-31 T^-0.6, 100-5000 K, dlog k = +-0.5.'
      write(*,'(a)') '       T[K]      mol_rates R15      Baulch 1992     Baulch/R15'
      do i = 1, 7
        T = tgrid(i)
        kf = 8.0d-33*(300.0d0/T)**0.6d0
        kv = 2.7d-31*T**(-0.6d0)
        write(*,'(f10.1,2es16.6,f14.5)') T, kf, kv, kv/kf
      end do
      write(*,*)

! ---------------------------------------------------------------------
      write(*,'(a)') '=== 2. Shomate table against the NIST Chemistry WebBook ==='
      write(*,'(a)') '    WebBook entries at 298.15 K (Chase 1998), read from the'
      write(*,'(a)') '    gas-phase thermochemistry pages:'
      write(*,'(a)') '      species   dfH(298.15) [kJ/mol]     S(298.15) [J/mol/K]'
      write(*,'(a)') '      H            217.9994                139.8711 (Shomate G)'
      write(*,'(a)') '      H2             0.0                   130.680'
      write(*,'(a)') '      O            249.18  +- 0.10         161.059 +- 0.003'
      write(*,'(a)') '      OH            38.98706               183.708 (Shomate)'
      write(*,'(a)') '      H2O         -241.8264                188.835 (Shomate)'
      write(*,'(a)') '      CO          -110.5271                197.663 (Shomate)'
      write(*,*)
      write(*,'(a)') '    The table below is also what a2_m1_thermo_xcheck.py reads back,'
      write(*,'(a)') '    so that the NASA-9 comparison cannot drift away from the module.'
      write(*,*)
      write(*,'(a)') 'SHOMATE_TABLE  T[K]  species  H[kJ/mol]  S[J/mol/K]  G[kJ/mol]'
      do j = 0, 7
        if (j == 0) then
          T = 298.15d0
        else
          T = tgrid(j)
        end if
        do i = 1, n_ox_sp
          write(*,'(a,f10.2,2x,a4,3f20.6)') 'SHOMATE_TABLE ', T,           &
               ox_species_name(i), enthalpy_shomate(i, T),                 &
               entropy_shomate(i, T), gibbs_energy_shomate(i, T)/1.0d3
        end do
      end do
      write(*,*)

! ---------------------------------------------------------------------
      write(*,'(a)') '=== 3. Round-trip residual of rate_from_detailed_balance ==='
      write(*,'(a)') '    k_fwd -> k_rev -> k_fwd must return the input exactly.'
      write(*,'(a)') '       T[K]      reaction                 |k_back/k_fwd - 1|'
      do i = 1, 7
        T = tgrid(i)
        react2 = (/ ith_OH, ith_H2 /); prod2 = (/ ith_H2O, ith_H /)
        kf = rk_O1_OH_H2_water(T)
        kr = rate_from_detailed_balance(kf, react2, prod2, T)
        kf_back = rate_from_detailed_balance(kr, prod2, react2, T)
        write(*,'(f10.1,a,es22.4)') T, '  O1  OH+H2 <-> H2O+H  ', abs(kf_back/kf - 1.0d0)
      end do
      do i = 1, 7
        T = tgrid(i)
        react2 = (/ ith_O, ith_H2 /); prod2 = (/ ith_OH, ith_H /)
        kf = rk_O2_O_H2_hydroxyl(T)
        kr = rate_from_detailed_balance(kf, react2, prod2, T)
        kf_back = rate_from_detailed_balance(kr, prod2, react2, T)
        write(*,'(f10.1,a,es22.4)') T, '  O2  O+H2  <-> OH+H   ', abs(kf_back/kf - 1.0d0)
      end do
      do i = 1, 7
        T = tgrid(i)
        react3 = (/ ith_O, ith_H, ith_H /); prod3 = (/ ith_H2O, ith_H2O, ith_H2O /)
        react3(3) = ith_H
        kf = rk_O9_H_OH_assoc(T)
        kr = rate_from_detailed_balance(kf, (/ ith_H, ith_OH /), (/ ith_H2O /), T)
        kf_back = rate_from_detailed_balance(kr, (/ ith_H2O /), (/ ith_H, ith_OH /), T)
        write(*,'(f10.1,a,es22.4)') T, '  O9  H+OH  <-> H2O    ', abs(kf_back/kf - 1.0d0)
      end do
      write(*,*)

! ---------------------------------------------------------------------
      write(*,'(a)') '=== 4. Physical detailed-balance test on O10 ==='
      write(*,'(a)') '    Photochem transcribes OH + OH -> H2O + O (Baulch et al. 1992,'
      write(*,'(a)') '    Lifshitz & Michael 1991).  VULCAN transcribes the OPPOSITE'
      write(*,'(a)') '    direction, O + H2O -> OH + OH (id 5, 250-2400 K).  Reversing the'
      write(*,'(a)') '    first with the Shomate data must reproduce the second.'
      write(*,'(a)') '    Both forward transcriptions are reversed, because the test'
      write(*,'(a)') '    discriminates between them.'
      write(*,'(a)') '       T[K]         K_c      k_rev(2005 fwd) k_rev(1992 fwd)  k_published   r05   r92'
      do i = 1, 7
        T = tgrid(i)
        kf = rk_O10_OH_OH_water(T)
        kc = equilibrium_constant_conc((/ ith_OH, ith_OH /), (/ ith_H2O, ith_O /), T)
        kr = rate_from_detailed_balance(kf, (/ ith_OH, ith_OH /), (/ ith_H2O, ith_O /), T)
        kf_back = rate_from_detailed_balance(2.5d-15*T**1.14d0*exp(-50.0d0/T),&
                     (/ ith_OH, ith_OH /), (/ ith_H2O, ith_O /), T)
        kv = 8.20d-14*T**0.95d0*exp(-8570.0d0/T)
        write(*,'(f10.1,4es15.5,2f7.3)') T, kc, kr, kf_back, kv, kr/kv,      &
             kf_back/kv
      end do
      write(*,*)
      write(*,'(a)') '    The 1992 pairing agrees far better than the 2005 one, which'
      write(*,'(a)') '    suggests VULCAN''s id 5 was itself obtained by reversing the 1992'
      write(*,'(a)') '    OH + OH recommendation rather than transcribed independently'
      write(*,'(a)') '    (the network file gives it no reference).  This test is therefore'
      write(*,'(a)') '    weaker than it looks; the O1 test below is the real one.'
      write(*,*)
      write(*,'(a)') '    Second published test, and the stronger one: Baulch et al.'
      write(*,'(a)') '    (1992) Table 1 p. 418 recommends the reverse of O1 as a rate in'
      write(*,'(a)') '    its own right, H + H2O -> OH + H2 = 7.5e-16 T^1.6 exp(-9270/T),'
      write(*,'(a)') '    300-2500 K, dlog k = +-0.2, data sheet p. 504.  Reversing the'
      write(*,'(a)') '    forward rate with the Shomate data must reproduce it.'
      write(*,'(a)') '       T[K]   k_rev(thermo,2005 fwd)  k_rev(thermo,1992 fwd)   k_published   ratios'
      do i = 1, 7
        T = tgrid(i)
        kr = rate_from_detailed_balance(rk_O1_OH_H2_water(T),                &
                 (/ ith_OH, ith_H2 /), (/ ith_H2O, ith_H /), T)
        kc = rate_from_detailed_balance(1.7d-16*T**1.6d0*exp(-1660.0d0/T),   &
                 (/ ith_OH, ith_H2 /), (/ ith_H2O, ith_H /), T)
        kv = 7.5d-16*T**1.6d0*exp(-9270.0d0/T)
        write(*,'(f10.1,3es20.5,2f9.3)') T, kr, kc, kv, kr/kv, kc/kv
      end do
      write(*,*)
      write(*,'(a)') '    Third published test: the O1 data sheet (Baulch et al. 1992,'
      write(*,'(a)') '    p. 552) prints its own thermodynamic data for OH + H2 -> H + H2O,'
      write(*,'(a)') '    dH(298) = -62.9 kJ/mol and dS(298) = -10.9 J/K/mol, taken there'
      write(*,'(a)') '    from the Sandia Chemkin compilation.  The Shomate table gives:'
      write(*,'(a,f10.4,a)') '       dH(298.15) = ',                             &
           enthalpy_shomate(ith_H2O, 298.15d0) + enthalpy_shomate(ith_H, 298.15d0)&
         - enthalpy_shomate(ith_OH, 298.15d0) - enthalpy_shomate(ith_H2, 298.15d0),&
           ' kJ/mol   (published -62.9)'
      write(*,'(a,f10.4,a)') '       dS(298.15) = ',                             &
           entropy_shomate(ith_H2O, 298.15d0) + entropy_shomate(ith_H, 298.15d0)  &
         - entropy_shomate(ith_OH, 298.15d0) - entropy_shomate(ith_H2, 298.15d0), &
           ' J/K/mol  (published -10.9)'
      write(*,'(a)') '    and the same data sheet fits the equilibrium constant as'
      write(*,'(a)') '    Kp = 0.113 T^0.0639 exp(7680/T) (dn = 0, so Kp = Kc):'
      write(*,'(a)') '       T[K]        K_c(Shomate)     K_p(published)      ratio'
      do i = 1, 7
        T = tgrid(i)
        kc = equilibrium_constant_conc((/ ith_OH, ith_H2 /),                 &
                                       (/ ith_H2O, ith_H /), T)
        kv = 0.113d0*T**0.0639d0*exp(7680.0d0/T)
        write(*,'(f10.1,2es20.5,f12.4)') T, kc, kv, kc/kv
      end do
      write(*,*)
      write(*,'(a)') '    The same reversal applied to O1 and O2, for which no published'
      write(*,'(a)') '    reverse exists in either network (both codes generate it):'
      write(*,'(a)') '       T[K]   k_rev(H2O+H->OH+H2)  k_rev(OH+H->O+H2)'
      do i = 1, 7
        T = tgrid(i)
        kr = rate_from_detailed_balance(rk_O1_OH_H2_water(T),                &
                 (/ ith_OH, ith_H2 /), (/ ith_H2O, ith_H /), T)
        kf = rate_from_detailed_balance(rk_O2_O_H2_hydroxyl(T),              &
                 (/ ith_O, ith_H2 /), (/ ith_OH, ith_H /), T)
        write(*,'(f10.1,2es22.5)') T, kr, kf
      end do
      write(*,*)

! ---------------------------------------------------------------------
      write(*,'(a)') '=== 5. Measured size of each reaction at the HD 189733 b base ==='
      write(*,'(a)') '    P = 0.9386 dyn cm^-2, T = 863.911 K, n = 7.869e12 cm^-3;'
      write(*,'(a)') '    abundances from the P1 Photochem/NCHO run (level 121).'
      write(*,'(a)') '    Rates are gross forward rates [cm^-3 s^-1]; the column "/O1" is'
      write(*,'(a)') '    what the include/exclude verdict of the audit rests on.'
      write(*,*)
      T = T_base
      r1 = rk_O1_OH_H2_water(T)*n_OH*n_H2
      write(*,'(a)') '    channel                                   rate            /O1'
      call report('O1   OH + H2 -> H2O + H         ', r1, r1)
      call report('O1r  H2O + H -> OH + H2  (thermo)',                        &
           rate_from_detailed_balance(rk_O1_OH_H2_water(T),                   &
             (/ ith_OH, ith_H2 /), (/ ith_H2O, ith_H /), T)*n_H2O*n_H, r1)
      call report('O2   O + H2 -> OH + H           ',                         &
           rk_O2_O_H2_hydroxyl(T)*n_O*n_H2, r1)
      call report('O2r  OH + H -> O + H2    (thermo)',                        &
           rate_from_detailed_balance(rk_O2_O_H2_hydroxyl(T),                 &
             (/ ith_O, ith_H2 /), (/ ith_OH, ith_H /), T)*n_OH*n_H, r1)
      call report('O6   O(1D) + H2 -> OH + H       ',                         &
           rk_O6_O1D_H2_hydroxyl()*n_O1D*n_H2, r1)
      call report('O8   O + H + M -> OH + M        ',                         &
           lindemann_rate(rk_O8_O_H_assoc(T), 1.0d-11, n_base)*n_O*n_H, r1)
      call report('O9   H + OH + M -> H2O + M      ',                         &
           lindemann_rate(rk_O9_H_OH_assoc(T), 2.7d-10*exp(-75.0d0/T),        &
             n_base)*n_H*n_OH, r1)
      call report('O10  OH + OH -> H2O + O         ',                         &
           rk_O10_OH_OH_water(T)*n_OH*n_OH, r1)
      call report('O12  O(1D) + H2O -> OH + OH     ',                         &
           rk_O12_O1D_H2O_hydroxyl()*n_O1D*n_H2O, r1)
      write(*,*)
      write(*,'(a,es12.4)') '    Reduced pressure Pr of O8 at this level: ',   &
           rk_O8_O_H_assoc(T)*n_base/1.0d-11
      write(*,'(a,es12.4)') '    Reduced pressure Pr of O9 at this level: ',   &
           rk_O9_H_OH_assoc(T)*n_base/(2.7d-10*exp(-75.0d0/T))
      write(*,'(a,es12.4,a)') '    O(1D) lifetime against O6:  ',              &
           1.0d0/(rk_O6_O1D_H2_hydroxyl()*n_H2), ' s'
      write(*,'(a,es12.4,a)') '    O(1D) lifetime against O12: ',              &
           1.0d0/(rk_O12_O1D_H2O_hydroxyl()*n_H2O), ' s'
      write(*,*)

! ---------------------------------------------------------------------
      write(*,'(a)') '=== 6. Photolysis bands and thresholds ==='
      write(*,'(a)') '    band  interval [A]        qy(OH+H)  qy(H2+O1D)  qy(O+H+H)  qy(OH->O+H)'
      do i = 1, n_fuv_band
        write(*,'(4x,a2,2x,f7.1,a,f7.1,4f12.2)') fuv_band_name(i),            &
             fuv_band_lo_A(i), ' -', fuv_band_hi_A(i),                        &
             qy_H2O_OH_H(i), qy_H2O_H2_O1D(i), qy_H2O_O_H_H(i), qy_OH_O_H(i)
      end do
      write(*,*)
      write(*,'(a)') '    channel threshold energies (dH at 298.15 K from the Shomate table)'
      write(*,'(a,es13.5,a,f8.4,a)') '      H2O -> OH + H    ',                &
           photolysis_threshold_erg(ich_H2O_OH_H), ' erg = ',                  &
           photolysis_threshold_erg(ich_H2O_OH_H)/1.602176634d-12, ' eV'
      write(*,'(a,es13.5,a,f8.4,a)') '      H2O -> H2 + O    ',                &
           photolysis_threshold_erg(ich_H2O_H2_O), ' erg = ',                  &
           photolysis_threshold_erg(ich_H2O_H2_O)/1.602176634d-12, ' eV'
      write(*,'(a,es13.5,a,f8.4,a)') '      H2O -> O + H + H ',                &
           photolysis_threshold_erg(ich_H2O_O_H_H), ' erg = ',                 &
           photolysis_threshold_erg(ich_H2O_O_H_H)/1.602176634d-12, ' eV'
      write(*,'(a,es13.5,a,f8.4,a)') '      OH  -> O + H     ',                &
           photolysis_threshold_erg(ich_OH_O_H), ' erg = ',                    &
           photolysis_threshold_erg(ich_OH_O_H)/1.602176634d-12, ' eV'
      write(*,'(a,es13.5,a,f8.4,a)') '      O(1D) excitation ', e_excite_O1D_erg,&
           ' erg = ', e_excite_O1D_erg/1.602176634d-12, ' eV'
      write(*,*)
      write(*,'(a)') '# end'

      contains

      double precision function tgrid(i) result(T)
      integer, intent(in) :: i
      real(dp), parameter :: tt(7) = (/ 300.0d0, 500.0d0, 864.0d0,        &
                                       1000.0d0, 1500.0d0, 2000.0d0,      &
                                       2500.0d0 /)
      T = tt(i)
      end function

      subroutine conservation_report()
      ! Species of the A2 network with their (H, O, C, charge) counts.  M is
      ! a third body and appears on both sides, so it cancels and is not
      ! listed.  hv carries no nuclei and no charge.
      integer, parameter :: nsp = 8
      character(len=5), parameter :: nm(nsp) = (/ 'H    ', 'H2   ',       &
           'O    ', 'O(1D)', 'OH   ', 'H2O  ', 'CO   ', 'hv   ' /)
      integer, parameter :: cnt(4, nsp) = reshape( (/                     &
      !     H   O   C   q
            1,  0,  0,  0,                                                &
            2,  0,  0,  0,                                                &
            0,  1,  0,  0,                                                &
            0,  1,  0,  0,                                                &
            1,  1,  0,  0,                                                &
            2,  1,  0,  0,                                                &
            0,  1,  1,  0,                                                &
            0,  0,  0,  0 /), (/ 4, nsp /) )
      integer, parameter :: nrx = 10
      character(len=34), parameter :: rxname(nrx) = (/                    &
           'O1   OH + H2   -> H2O + H         ',                          &
           'O2   O  + H2   -> OH  + H         ',                          &
           'O3   H2O + hv  -> OH  + H         ',                          &
           'O4   H2O + hv  -> H2  + O(1D)     ',                          &
           'O5   H2O + hv  -> O   + H + H     ',                          &
           'O6   O(1D)+ H2 -> OH  + H         ',                          &
           'O7   OH  + hv  -> O   + H         ',                          &
           'O8   O + H + M -> OH  + M         ',                          &
           'O9   H + OH+ M -> H2O + M         ',                          &
           'O10  OH + OH   -> H2O + O         ' /)
      ! reactant and product species indices, 0 = unused
      integer, parameter :: irc(3, nrx) = reshape( (/                     &
           5, 2, 0,   3, 2, 0,   6, 8, 0,   6, 8, 0,   6, 8, 0,           &
           4, 2, 0,   5, 8, 0,   3, 1, 0,   1, 5, 0,   5, 5, 0 /),        &
           (/ 3, nrx /) )
      integer, parameter :: ipr(3, nrx) = reshape( (/                     &
           6, 1, 0,   5, 1, 0,   5, 1, 0,   2, 4, 0,   3, 1, 1,           &
           5, 1, 0,   3, 1, 0,   5, 0, 0,   6, 0, 0,   6, 3, 0 /),        &
           (/ 3, nrx /) )
      integer :: ir, ie, k, bal(4)
      logical :: ok_all
      character(len=8), parameter :: elname(4) =                          &
           (/ 'H nuclei', 'O nuclei', 'C nuclei', 'charge  ' /)

      write(*,'(a)') '    species    H   O   C   q'
      do k = 1, nsp
        write(*,'(4x,a8,4i4)') nm(k), cnt(1,k), cnt(2,k), cnt(3,k), cnt(4,k)
      end do
      write(*,*)
      write(*,'(a)') '    reaction                            dH   dO   dC   dq'
      ok_all = .true.
      do ir = 1, nrx
        bal = 0
        do k = 1, 3
          if (ipr(k, ir) > 0) bal = bal + cnt(:, ipr(k, ir))
          if (irc(k, ir) > 0) bal = bal - cnt(:, irc(k, ir))
        end do
        write(*,'(4x,a34,4i5)') rxname(ir), (bal(ie), ie = 1, 4)
        if (any(bal /= 0)) ok_all = .false.
      end do
      write(*,*)
      if (ok_all) then
        write(*,'(a)') '    CONSERVATION: PASS -- every reaction balances H, O, C and charge.'
      else
        write(*,'(a)') '    CONSERVATION: FAIL'
      end if
      write(*,'(a)') '    CO appears in no reaction: decision D4 carries it as an'
      write(*,'(a)') '    unreactive oxygen reservoir, so its C and O move only with the'
      write(*,'(a)') '    gas.  O12 is excluded from the set and is balanced separately:'
      bal = cnt(:, 4) + cnt(:, 6) - cnt(:, 5) - cnt(:, 5)
      write(*,'(a,4i5)') '    O12  O(1D) + H2O -> OH + OH        ',       &
           (bal(ie), ie = 1, 4)
      end subroutine

      subroutine report(label, rate, ref)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: rate, ref
      write(*,'(4x,a,es16.5,es14.4)') label, rate, rate/ref
      end subroutine

      end program a2_m1_rate_check
