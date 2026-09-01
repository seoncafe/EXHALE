      program a2_m2_kinetics_check
      ! Standalone checks of the quantities milestone M2 of
      ! docs/a2_oxygen_option_design.md adds on top of the M1 audit. Nothing
      ! here is part of any build; to regenerate, from src/tests/a2_m2/:
      !
      !   gfortran -O2 -o a2_m2_kinetics_check.x \
      !       ../../modules/lower_atmosphere/oxygen_rates.f90 \
      !       ../../modules/lower_atmosphere/water_photolysis.f90 \
      !       a2_m2_kinetics_check.f90
      !   ./a2_m2_kinetics_check.x > a2_m2_kinetics_check.out
      !
      ! Sections:
      !   0  the atomic-carbon thermodynamic row added at M2, against the
      !      NIST-JANAF gas-phase entry it was transcribed from
      !   1  the CO <-> C + O reservoir over temperature: where it turns over
      !   2  the O / OH / H2O chemical-equilibrium partition, and the
      !      statement gate G3 rests on -- that the two-row oxygen system at
      !      zero photolysis has exactly that partition as its root
      !   3  the FUV band constants and the energy ledger they have to obey

      use oxygen_rates
      use water_photolysis
      implicit none

      integer, parameter :: dp = kind(1.0d0)
      integer  :: i, ib
      real(dp) :: T, nH, nC, nO, nCO, f_oh, f_h2o, f_o
      real(dp) :: nh2, nhi, kc1, kc2, ratio
      real(dp) :: k1, k1r, k2, k2r, res_oh, res_h2o, scal
      real(dp) :: n_oh, n_h2o, n_o0, nOfam
      real(dp) :: T_list(8)
      data T_list / 300.0d0, 500.0d0, 864.0d0, 1140.0d0, 1500.0d0,        &
                    2000.0d0, 2500.0d0, 3000.0d0 /

      call water_photolysis_init

      write(*,'(a)') '=== A2 M2 checks ==='
      write(*,'(a)') ''

      ! ---------------------------------------------------------------
      write(*,'(a)') '--- 0. atomic carbon, the row added at M2 ---'
      write(*,'(a)') 'NIST-JANAF gas-phase C: dfH(298.15) = 716.68'//     &
                     ' kJ/mol, S(298.15) = 158.100 J/mol/K'
      write(*,'(a,f14.6,a,f14.6)') 'this table:  dfH = ',                 &
           enthalpy_shomate(ith_C, 298.15d0), '   S = ',                  &
           entropy_shomate(ith_C, 298.15d0)
      write(*,'(a)') 'NIST-JANAF gas-phase CO: dfH(298.15) ='//           &
                     ' -110.53 kJ/mol'
      write(*,'(a,f14.6)') 'this table:  dfH = ',                         &
           enthalpy_shomate(ith_CO, 298.15d0)
      write(*,'(a,f14.4,a)') 'C=O bond, dH(298.15) of CO -> C + O = ',    &
           enthalpy_shomate(ith_C, 298.15d0)                              &
           + enthalpy_shomate(ith_O, 298.15d0)                            &
           - enthalpy_shomate(ith_CO, 298.15d0), ' kJ/mol'
      write(*,'(a)') ''

      ! ---------------------------------------------------------------
      write(*,'(a)') '--- 1. the CO reservoir, CO <-> C + O ---'
      write(*,'(a)') 'Solar C/H = 2.69e-4, O/H = 4.90e-4 (the metals.inp'//&
                     ' of the regression cases).'
      write(*,'(a)') 'f_CO is the fraction of the CARBON locked in CO;'// &
                     ' 1 - n_CO/n_O is the oxygen left'
      write(*,'(a)') 'to the water family.'
      write(*,'(a)') '     T[K]        n_H        K_c[cm^-3]      n_CO'// &
                     '        f_CO     O left'
      do i = 1, 16
        T  = 500.0d0 + 250.0d0*(i-1)
        nH = 1.0d13*exp(-(T - 900.0d0)/900.0d0)
        if (nH .lt. 1.0d4) nH = 1.0d4
        nC = 2.69d-4*nH
        nO = 4.90d-4*nH
        nCO = co_equilibrium_density(nC, nO, T)
        write(*,'(f9.1,5es12.4)') T, nH,                                  &
             equilibrium_constant_conc((/ ith_CO /), (/ ith_C, ith_O /), T),&
             nCO, nCO/nC, 1.0d0 - nCO/nO
      end do
      write(*,'(a)') ''

      ! ---------------------------------------------------------------
      write(*,'(a)') '--- 2. O / OH / H2O chemical equilibrium (gate G3) ---'
      write(*,'(a)') 'Partition of the free oxygen family at n_H2/n_HI'// &
                     ' = 5, from the two audited'
      write(*,'(a)') 'equilibria O + H2 <-> OH + H and OH + H2 <->'//     &
                     ' H2O + H.'
      write(*,'(a)') '     T[K]     Kc(O2)      Kc(O1)        f_O'//      &
                     '        f_OH       f_H2O'
      nh2 = 5.0d12
      nhi = 1.0d12
      do i = 1, 8
        T = T_list(i)
        kc2 = equilibrium_constant_conc((/ ith_O, ith_H2 /),              &
                                        (/ ith_OH, ith_H /), T)
        kc1 = equilibrium_constant_conc((/ ith_OH, ith_H2 /),             &
                                        (/ ith_H2O, ith_H /), T)
        call oxygen_chemical_equilibrium_fractions(T, nh2, nhi,           &
                                                   f_oh, f_h2o)
        f_o = 1.0d0 - f_oh - f_h2o
        write(*,'(f9.1,5es12.4)') T, kc2, kc1, f_o, f_oh, f_h2o
      end do
      write(*,'(a)') ''
      write(*,'(a)') 'G3: the two oxygen-carrier rows of the coupled'//   &
                     ' system, with every photolysis'
      write(*,'(a)') 'rate at zero, must vanish on exactly that'//        &
                     ' partition. The rows are'
      write(*,'(a)') '  OH : k1r nH2O nHI + k2 nO nH2 - k1 nOH nH2'//     &
                     ' - k2r nOH nHI'
      write(*,'(a)') '  H2O: k1 nOH nH2 - k1r nH2O nHI'
      write(*,'(a)') 'with k1r and k2r from detailed balance. Residuals'//&
                     ' below are scaled by the'
      write(*,'(a)') 'largest term of each row, which is what the'//      &
                     ' solver itself works on.'
      write(*,'(a)') '     T[K]    res_OH/scale  res_H2O/scale'
      nOfam = 1.0d9
      do i = 1, 8
        T   = T_list(i)
        k1  = rk_O1_OH_H2_water(T)
        k1r = rate_from_detailed_balance(k1, (/ ith_OH, ith_H2 /),        &
                                             (/ ith_H2O, ith_H /), T)
        k2  = rk_O2_O_H2_hydroxyl(T)
        k2r = rate_from_detailed_balance(k2, (/ ith_O, ith_H2 /),         &
                                             (/ ith_OH, ith_H /), T)
        call oxygen_chemical_equilibrium_fractions(T, nh2, nhi,           &
                                                   f_oh, f_h2o)
        n_oh  = f_oh *nOfam
        n_h2o = f_h2o*nOfam
        n_o0  = (1.0d0 - f_oh - f_h2o)*nOfam
        res_oh  = k1r*n_h2o*nhi + k2*n_o0*nh2                             &
                - k1*n_oh*nh2 - k2r*n_oh*nhi
        res_h2o = k1*n_oh*nh2 - k1r*n_h2o*nhi
        scal = max(abs(k1r*n_h2o*nhi), abs(k2*n_o0*nh2))
        scal = max(scal, abs(k1*n_oh*nh2))
        write(*,'(f9.1,2es16.6)') T, res_oh/scal, res_h2o/scal
      end do
      write(*,'(a)') ''

      ! ---------------------------------------------------------------
      write(*,'(a)') '--- 3. FUV bands ---'
      write(*,'(a)') 'band   lo[A]   hi[A]   sigma_H2O    sigma_OH'//     &
                     '     <hv>[eV]   <E>H2O[eV]  <E>OH[eV]'
      do ib = 1, n_fuv_band
        write(*,'(2x,a2,2f8.1,2es13.5,3f12.4)') fuv_band_name(ib),        &
             fuv_band_lo_A(ib), fuv_band_hi_A(ib),                        &
             sigma_H2O_band(ib), sigma_OH_band(ib),                       &
             e_photon_flat_band(ib)/1.602176634d-12,                      &
             e_photon_H2O_band(ib)/1.602176634d-12,                       &
             e_photon_OH_band(ib)/1.602176634d-12
      end do
      write(*,'(a)') ''
      write(*,'(a)') 'Heat per dissociation [eV], and the threshold it'// &
                     ' is measured against.'
      write(*,'(a)') 'The heating uses <hv>, not <E>: the band-average'// &
                     ' transmission absorbs every'
      write(*,'(a)') 'photon of the band in the saturated limit, so'//    &
                     ' charging each of them <E> would'
      write(*,'(a)') 'deposit more energy than the band carries (see'//   &
                     ' water_photolysis.f90 sec. 2).'
      write(*,'(a)') 'The last column is what that costs in the'//        &
                     ' optically thin limit.'
      write(*,'(a)') 'band  heat_H2O[eV]  heat_OH[eV]   <E>/<hv> H2O'//   &
                     '   <E>/<hv> OH'
      do ib = 1, n_fuv_band
        write(*,'(2x,a2,2f13.4,2f15.4)') fuv_band_name(ib),               &
             heat_per_water_dissociation(ib)/1.602176634d-12,             &
             heat_per_hydroxyl_dissociation(ib)/1.602176634d-12,          &
             e_photon_H2O_band(ib)/e_photon_flat_band(ib),                &
             e_photon_OH_band(ib)/e_photon_flat_band(ib)
      end do
      write(*,'(a)') ''
      write(*,'(a)') 'Quantum yields and the yield-weighted threshold'//  &
                     ' [eV] they build:'
      write(*,'(a)') 'band   OH+H   H2+O(1D)  O+H+H   Eth_H2O[eV]'
      do ib = 1, n_fuv_band
        write(*,'(2x,a2,3f9.2,f13.4)') fuv_band_name(ib),                 &
             qy_H2O_OH_H(ib), qy_H2O_H2_O1D(ib), qy_H2O_O_H_H(ib),        &
             (e_photon_flat_band(ib)                                      &
              - heat_per_water_dissociation(ib))/1.602176634d-12
      end do
      write(*,'(a,f10.4,a)') 'OH threshold (single merged channel): ',    &
           photolysis_threshold_erg(ich_OH_O_H)/1.602176634d-12, ' eV'
      write(*,'(a)') ''
      write(*,'(a)') 'Beam identity the band ledger of'//                 &
                     ' output/FUV_bands.txt tests: for j = s N exp(-tau)'
      write(*,'(a)') 'and dtau/dr = -s n, the photons absorbed over a'//  &
                     ' column are exactly'
      write(*,'(a)') 'N (1 - exp(-tau)), whatever the density profile.'// &
                     ' A run measures the'
      write(*,'(a)') 'discretization of that identity; it is not'//       &
                     ' reproduced here.'

      end program a2_m2_kinetics_check
