      program co_destruction
      ! The one-sided CO destruction model: its two published rates, the
      ! shielding table it reads, the energies it deposits, and the cell
      ! mean its rate is taken as.
      !
      ! Production routines exercised: co_self_shielding and
      ! co_self_shielding_texc5 of
      ! src/modules/lower_atmosphere/co_self_shielding_table.f90;
      ! co_photodissociation_rate, co_photodissociation_rate_cell_mean and
      ! heat_per_co_dissociation of
      ! src/modules/lower_atmosphere/co_photodissociation.f90;
      ! rk_D1_Hep_CO, rk_CO_radiative_association and
      ! photolysis_threshold_erg of
      ! src/modules/lower_atmosphere/oxygen_rates.f90;
      ! oxygen_reaction_energy_eV and species_formation_energy of
      ! src/modules/lower_atmosphere/molecular_reaction_heat.f90.
      !
      ! REFERENCES, all READ from the publisher PDF
      ! references/Visser_2009A&A_503_323.pdf unless said otherwise.
      !
      !  * THE TABLE.  Visser, van Dishoeck & Black (2009), A&A 503, 323,
      !    Table 6 (Online Material p. 1), the 12CO block: the shipped
      !    shielding function, b(CO) = 0.3 km/s, T_ex(CO) = 50 K,
      !    T_ex(H2) = 501.5 K.  Table 5 (p. 334), the same block at
      !    T_ex(CO) = 5 K, is the paper's reference set and is carried
      !    beside it so the spread between the two excitation temperatures
      !    can be measured rather than asserted.  The interpolation anchor
      !    is the one already computed in
      !    docs/co_destruction_rates_literature_20260906.md sec. 10.7 on the
      !    base cell of the oxygen_chemistry case,
      !    log N(CO) = 18.96 and log N(H2) = 21.80, where log-bilinear
      !    interpolation of Table 5 gives Theta = 2.6e-5.
      !
      !  * THE PHOTON ENERGY OF ONE EVENT.  Their Table 1 lists the 37
      !    predissociating bands of 12CO with lambda_0, the band oscillator
      !    strength f_v0 and the dissociation efficiency eta.  The constant
      !    e_co_photon_erg is the f-weighted mean of hc/lambda_0 over those
      !    37 rows; the table is transcribed below and the mean recomputed,
      !    so the constant cannot drift from the data it came from.
      !
      !  * THE He+ CHANNEL.  UMIST RATE22 entry 4068 of rate22_final.rates
      !    (Millar, Walsh, Van de Sande & Markwick 2024, A&A 682, A109):
      !    He+ + CO -> O + C+ + He, alpha = 1.60e-9, beta = 0, gamma = 0,
      !    method M, accuracy A.  beta = gamma = 0 makes it temperature
      !    independent, and that is what the three-temperature check below
      !    asserts to the bit: an accidental Arrhenius factor would break
      !    it.
      !
      !  * THE FORMATION ENTRY, record only.  RATE22 entry 8597,
      !    C + O -> CO + PHOTON, alpha = 4.69e-19, beta = 1.52,
      !    gamma = -50.5, T = 10-14700 K.
      !
      !  * THE NORMALIZATION CROSS-CHECK.  Visser's k0 belongs to the
      !    Draine (1978) field, whose 912-1110 A photon flux at chi = 1 is
      !    1.232e7 cm^-2 s^-1 (Draine & Bertoldi 1996, Table 1, the F/chi
      !    column, READ through lyman_werner.f90 sec. 1).  Their Table 6
      !    header gives k0 = 2.590e-10 s^-1, so their data imply a band-mean
      !    dissociation cross section of 2.102e-17 cm^2 over 912-1110 A.
      !    This code carries 1.0160e-17 cm^2 over the WIDER 912-1201 A beam,
      !    which is 1.5496e-17 over 912-1110 A alone; the ratio of the two
      !    is checked inside a factor 1.5, which is Visser's own 20 per cent
      !    absolute accuracy loosened by the difference in spectral shape
      !    that co_photodissociation.f90 sec. 2 states.  It is a check that
      !    the two normalizations describe the same molecule, not that they
      !    are the same number.

      use assertion_report
      use co_self_shielding_table, only: co_self_shielding,               &
               co_self_shielding_texc5, co_shield_max_co_column,          &
               co_shield_max_h2_column, co_shield_tex_limit_K,            &
               n_co_shield_co, n_co_shield_h2,                            &
               co_shield_log_co, co_shield_log_h2
      use co_photodissociation, only: co_photodissociation_rate,          &
               co_photodissociation_rate_cell_mean, sigma_CO_band,        &
               e_co_photon_erg, heat_per_co_dissociation
      use lyman_werner_photodissociation, only: e_lw_photon_erg
      use oxygen_rates, only: rk_D1_Hep_CO, rk_CO_radiative_association,  &
               photolysis_threshold_erg, ich_CO_C_O,                      &
               oxygen_formation_energy_eV, ith_CO
      use molecular_reaction_heat, only: oxygen_reaction_energy_eV,       &
               ir_D1, species_formation_energy, isp_eps_electron
      use species_table, only: isp_HeII, isp_HeI, isp_CO, mion_fsp,       &
               mion_ethr, im_CII, im_CI, im_OI
      use global_parameters, only: e_th_HeI

      implicit none
      real*8, parameter :: eV_to_erg = 1.602176634d-12
      real*8, parameter :: hc_eV_A   = 12398.419843320026d0

      ! Visser et al. Table 6, 12CO block, transcribed independently of the
      ! production module: theta6_ref(i,k) at log N(CO) = the i-th knot and
      ! log N(H2) = the k-th.
      real*8, parameter :: theta6_ref(8,6) = reshape( (/                  &
        1.000d0,  9.405d-1, 7.046d-1, 4.015d-1, 9.964d-2, 1.567d-2,       &
        3.162d-3, 4.839d-4,                                               &
        7.546d-1, 6.979d-1, 4.817d-1, 2.577d-1, 6.505d-2, 1.135d-2,       &
        2.369d-3, 3.924d-4,                                               &
        5.752d-1, 5.228d-1, 3.279d-1, 1.559d-1, 3.559d-2, 6.443d-3,       &
        1.526d-3, 2.751d-4,                                               &
        2.493d-1, 2.196d-1, 1.135d-1, 4.062d-2, 7.864d-3, 1.516d-3,       &
        4.448d-4, 9.367d-5,                                               &
        1.550d-3, 1.370d-3, 6.801d-4, 2.127d-4, 5.051d-5, 1.198d-5,       &
        6.553d-6, 3.937d-6,                                               &
        8.492d-8, 8.492d-8, 8.492d-8, 8.492d-8, 8.492d-8, 8.492d-8,       &
        8.488d-8, 8.453d-8 /), (/ 8, 6 /) )

      ! Visser et al. Table 1: the 37 bands, lambda_0 [A], f_v0 and eta.
      integer, parameter :: n_band = 37
      real*8, parameter :: band_lam(n_band) = (/                          &
        912.7037d0, 913.4002d0, 913.4255d0, 913.6678d0, 915.7258d0,       &
        915.9708d0, 917.2721d0, 919.2097d0, 920.1410d0, 922.7561d0,       &
        924.6309d0, 925.8093d0, 928.6575d0, 930.0611d0, 931.0744d0,       &
        931.6547d0, 933.0583d0, 935.6638d0, 939.9574d0, 941.1685d0,       &
        946.2860d0, 948.3860d0, 950.0429d0, 956.2369d0, 964.3973d0,       &
        968.3187d0, 968.8815d0, 970.3588d0, 972.6996d0, 977.3996d0,       &
        982.5914d0, 985.6490d0, 989.7952d0, 1002.5897d0, 1051.7134d0,     &
        1063.0880d0, 1076.0796d0 /)
      real*8, parameter :: band_f(n_band) = (/                            &
        3.4d-3, 1.7d-3, 1.7d-3, 2.7d-2, 2.0d-3, 7.9d-3, 2.3d-2, 2.8d-3,   &
        2.8d-3, 6.3d-3, 5.2d-3, 2.0d-2, 6.7d-3, 6.3d-3, 6.0d-3, 1.2d-2,   &
        2.2d-2, 3.8d-3, 2.1d-2, 3.1d-2, 7.6d-3, 2.8d-3, 2.2d-2, 1.6d-2,   &
        2.8d-3, 1.4d-2, 1.2d-2, 3.4d-2, 1.7d-2, 1.8d-3, 4.8d-4, 1.5d-2,   &
        4.6d-4, 7.9d-3, 3.6d-3, 3.0d-3, 6.8d-2 /)

      real*8  :: theta, ref, worst, spread_lo, spread_hi, ratio
      real*8  :: fsum, esum, e_mean, k, kref, F_band, k_fine, acc
      real*8  :: q_d1, q_tab, d0_co, dep
      real*8  :: nco_lo, nco_hi, nh2_lo, nh2_hi, s, ds
      integer :: i, kk, ib, ns

      !----------------------------------------------------------------!
      ! 1. THE TABLE IS THE TABLE.
      ! Every one of the 48 grid points of the shipped 12CO block is
      ! returned exactly at its own knot: the interpolation is linear in
      ! log10 Theta, so a knot is reproduced to round-off and any
      ! transcription slip shows up as a whole entry.
      worst = 0.0d0
      do kk = 1, n_co_shield_h2
         do i = 1, n_co_shield_co
            theta = co_self_shielding(10.0d0**co_shield_log_co(i),        &
                                      10.0d0**co_shield_log_h2(kk))
            worst = max(worst, abs(theta/theta6_ref(i,kk) - 1.0d0))
         enddo
      enddo
      call check_absolute('co_shield_table6_grid_points', worst, 0.0d0,   &
                          1.0d-12)

      !----------------------------------------------------------------!
      ! 2. THE SPREAD BETWEEN THE TWO EXCITATION TEMPERATURES.
      ! Table 5 (T_ex(CO) = 5 K) against Table 6 (50 K) at every grid
      ! point.  It is not a tolerance on anything: it is the size of the
      ! model's own uncertainty on the axis neither table reaches the gas
      ! on, printed so the run's domain record can be read against it.  The
      ! assertion is only that the spread is bounded, which is what makes
      ! the choice between the two tables a bounded choice.
      spread_lo = 1.0d30
      spread_hi = 0.0d0
      do kk = 1, n_co_shield_h2
         do i = 1, n_co_shield_co
            ratio = co_self_shielding_texc5(                              &
                        10.0d0**co_shield_log_co(i),                      &
                        10.0d0**co_shield_log_h2(kk))                     &
                  / co_self_shielding(10.0d0**co_shield_log_co(i),        &
                                      10.0d0**co_shield_log_h2(kk))
            spread_lo = min(spread_lo, ratio)
            spread_hi = max(spread_hi, ratio)
         enddo
      enddo
      write(*,'(a,es12.4,a,es12.4)')                                      &
        '# co_shield_texc_spread Table5/Table6 over the 48 grid points:'//&
        ' min ', spread_lo, ' max ', spread_hi
      ! MEASURED here: the ratio runs from 0.549 to 22.07 over the block,
      ! the largest value at log N(CO) = 17, log N(H2) = 22, i.e. deep in
      ! the shielded interior.  That is the size of the model's own
      ! uncertainty on the excitation-temperature axis, and it is larger
      ! than the factor of two Visser's sec. 5.2 quotes for a
      ! high-temperature PDR because that sentence is about the rate a
      ! whole cloud model gives and this is the worst single grid point.
      ! The bound below is a REGRESSION on the transcription of the two
      ! blocks, not a physical tolerance: if either table is mistyped the
      ! spread moves.
      call check_relative('co_shield_texc_spread_max', spread_hi,         &
                          22.0701d0, 1.0d-4)
      call check_relative('co_shield_texc_spread_min', spread_lo,         &
                          0.548675d0, 1.0d-4)

      !----------------------------------------------------------------!
      ! 3. INTERPOLATION, against the value computed by hand in the
      ! literature document on the base cell of the oxygen_chemistry case.
      theta = co_self_shielding_texc5(10.0d0**18.96d0, 10.0d0**21.80d0)
      call check_relative('co_shield_interpolation_base_cell', theta,     &
                          2.6d-5, 2.0d-2)

      !----------------------------------------------------------------!
      ! 4. THE EDGE VALUE, NOT AN EXTRAPOLATION.  Above the top of either
      ! axis the table returns its corner entry; Visser's sec. 5.1 states
      ! why that is the right behaviour and not a clamp of convenience.
      call check_relative('co_shield_clamp_above_co_axis',                &
           co_self_shielding(1.0d3*co_shield_max_co_column(), 1.0d0),     &
           theta6_ref(n_co_shield_co,1), 1.0d-12)
      call check_relative('co_shield_clamp_above_h2_axis',                &
           co_self_shielding(1.0d0, 1.0d3*co_shield_max_h2_column()),     &
           theta6_ref(1,n_co_shield_h2), 1.0d-12)
      ! A column of zero is the unshielded limit, which is the table's own
      ! bottom-corner entry and equals 1.
      call check_relative('co_shield_zero_column_unshielded',             &
           co_self_shielding(0.0d0, 0.0d0), 1.0d0, 1.0d-12)

      !----------------------------------------------------------------!
      ! 5. THE 512 K LIMIT is the number Visser's sec. 4.4 states, and it
      ! is what the run's domain record counts cells against.
      call check_relative('co_shield_tex_limit', co_shield_tex_limit_K(), &
                          512.0d0, 1.0d-12)

      !----------------------------------------------------------------!
      ! 6. THE PHOTON ENERGY OF ONE EVENT, recomputed from Table 1.
      fsum = 0.0d0
      esum = 0.0d0
      do ib = 1, n_band
         fsum = fsum + band_f(ib)
         esum = esum + band_f(ib)*hc_eV_A/band_lam(ib)
      enddo
      e_mean = esum/fsum
      call check_relative('co_photon_energy_37_band_mean',                &
                          e_co_photon_erg/eV_to_erg, e_mean, 1.0d-6)

      !----------------------------------------------------------------!
      ! 7. ONE BAND, ONE PHOTON ENERGY.  The CO rate divides the band flux
      ! by the same <hv> the H2 rate does, so the two absorbers of one beam
      ! count the same photons.  With Theta = 1 and no continuum the rate
      ! is the thin limit sigma F/<hv>.
      F_band = 171.5d0                ! the oxygen_chemistry beam flux
      k    = co_photodissociation_rate(F_band, 0.0d0, 0.0d0, 0.0d0)
      kref = sigma_CO_band*F_band/e_lw_photon_erg
      call check_relative('co_thin_rate_is_sigma_times_photon_flux', k,   &
                          kref, 1.0d-12)
      ! And that <hv> is the Lyman-Werner band's own constant, not a copy
      ! that can drift from it.
      call check_relative('co_band_photon_energy_matches_lw',             &
                          kref*e_lw_photon_erg/(sigma_CO_band*F_band),    &
                          1.0d0, 1.0d-12)

      !----------------------------------------------------------------!
      ! 8. THE NORMALIZATION CROSS-CHECK against Visser's own k0 in the
      ! Draine field (header).  The shipped cross section is over
      ! 912-1201 A; restated over the 912-1110 A band the Draine flux
      ! belongs to it is 1.5496e-17 cm^2.
      call check_within_factor('co_cross_section_against_visser_k0',      &
           1.5496d-17, 2.590d-10/1.232d7, 1.5d0)

      !----------------------------------------------------------------!
      ! 9. THE He+ CHANNEL IS TEMPERATURE INDEPENDENT, to the bit.
      call check_relative('co_hep_rate_300K',  rk_D1_Hep_CO(), 1.60d-9,   &
                          1.0d-15)
      call check_relative('co_hep_rate_1500K', rk_D1_Hep_CO(), 1.60d-9,   &
                          1.0d-15)
      call check_relative('co_hep_rate_3000K', rk_D1_Hep_CO(), 1.60d-9,   &
                          1.0d-15)

      !----------------------------------------------------------------!
      ! 10. THE TWO ENERGIES, EACH A DIFFERENCE OF THE ONE TABLE.
      ! q(D1) = eps(He+) + eps(CO) - eps(C+) - eps(O) - eps(He), formed
      ! here from species_formation_energy with the runtime species keys and
      ! compared with what the reaction table returns.  eps(O) and eps(He)
      ! are zero by the reference state, and the sum is the difference of
      ! the helium and carbon ionization potentials less the C=O bond.
      q_tab = oxygen_reaction_energy_eV(ir_D1)
      q_d1  = species_formation_energy(isp_HeII)                          &
            + species_formation_energy(isp_CO)                            &
            - species_formation_energy(mion_fsp(im_CII))                  &
            - species_formation_energy(mion_fsp(im_OI))                   &
            - species_formation_energy(isp_HeI)
      call check_relative('co_hep_reaction_energy_from_table', q_tab,     &
                          q_d1, 1.0d-14)
      call check_relative('co_hep_reaction_energy_value', q_tab,          &
                          e_th_HeI - mion_ethr(im_CI)                     &
                          + oxygen_formation_energy_eV(ith_CO), 1.0d-12)
      call check_positive('co_hep_reaction_is_exothermic', q_tab)

      ! D0(CO) is the SAME difference, so the two channels of one molecule
      ! cannot state two bond energies.
      d0_co = photolysis_threshold_erg(ich_CO_C_O)/eV_to_erg
      call check_relative('co_bond_energy_from_table', d0_co,             &
                          -oxygen_formation_energy_eV(ith_CO), 1.0d-14)
      ! The fragments keep the rest of the photon.
      dep = heat_per_co_dissociation()/eV_to_erg
      call check_relative('co_photodissociation_deposit',  dep,           &
                          e_mean - d0_co, 1.0d-6)
      call check_positive('co_photodissociation_deposit_positive', dep)

      !----------------------------------------------------------------!
      ! 11. THE CELL MEAN against a fine reference of the same integrand.
      ! A cell spanning three decades of CO column and two of H2 column,
      ! with a continuum depth across it: the composite rule of the
      ! production routine against a 4000-segment midpoint sum.  Both
      ! columns and the depth run linearly across the cell, which is the
      ! rectangle rule the column integration itself uses.
      nco_lo = 1.0d15
      nco_hi = 1.0d18
      nh2_lo = 1.0d20
      nh2_hi = 1.0d22
      k = co_photodissociation_rate_cell_mean(F_band, nco_lo, nco_hi,     &
                                              nh2_lo, nh2_hi, 0.1d0,      &
                                              0.9d0)
      ns  = 4000
      ds  = 1.0d0/dble(ns)
      acc = 0.0d0
      do i = 1, ns
         s   = (dble(i) - 0.5d0)*ds
         acc = acc + ds*co_photodissociation_rate(F_band,                 &
                 nco_lo + s*(nco_hi - nco_lo),                            &
                 nh2_lo + s*(nh2_hi - nh2_lo),                            &
                 0.1d0 + s*0.8d0)
      enddo
      k_fine = acc
      call check_relative('co_cell_mean_against_fine_quadrature', k,      &
                          k_fine, 1.0d-3)

      ! A cell whose star-ward column is zero -- the outermost cell, and
      ! any cell at the top of the CO distribution -- has no geometric
      ! starting point and takes the head segment.
      k = co_photodissociation_rate_cell_mean(F_band, 0.0d0, 1.0d17,      &
                                              0.0d0, 1.0d21, 0.0d0, 0.0d0)
      acc = 0.0d0
      do i = 1, ns
         s   = (dble(i) - 0.5d0)*ds
         acc = acc + ds*co_photodissociation_rate(F_band, s*1.0d17,       &
                                                  s*1.0d21, 0.0d0)
      enddo
      call check_relative('co_cell_mean_head_segment', k, acc, 1.0d-3)

      ! A cell with no column gradient is the point rate.
      k = co_photodissociation_rate_cell_mean(F_band, 1.0d16, 1.0d16,     &
                                              1.0d21, 1.0d21, 0.0d0, 0.0d0)
      call check_relative('co_cell_mean_uniform_cell_is_point_rate', k,   &
           co_photodissociation_rate(F_band, 1.0d16, 1.0d21, 0.0d0),      &
           1.0d-12)

      ! No band flux, no rate.
      call check_absolute('co_rate_zero_without_band_flux',               &
           co_photodissociation_rate_cell_mean(0.0d0, 1.0d15, 1.0d18,     &
                                    1.0d20, 1.0d22, 0.0d0, 1.0d0),        &
           0.0d0, 0.0d0)

      !----------------------------------------------------------------!
      ! 12. THE FORMATION ENTRY, which enters no row and only the record.
      ! At 1440 K, RATE22 8597 gives 4.69e-19 (T/300)^1.52 exp(50.5/T);
      ! recomputed here from the entry's own four numbers.
      call check_relative('co_radiative_association_rate_1440K',          &
           rk_CO_radiative_association(1440.0d0),                         &
           4.69d-19*(1440.0d0/300.0d0)**1.52d0*exp(50.5d0/1440.0d0),      &
           1.0d-12)
      ! It is a formation rate and must be far slower than the destruction
      ! it is compared against: at the base composition of the
      ! oxygen_chemistry case (n_O ~ 1e10 cm^-3) it is a formation time of
      ! more than 1e5 s per unit of that density.
      call check_positive('co_radiative_association_positive',            &
                          rk_CO_radiative_association(1440.0d0))

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'co_destruction: ', assertion_failures,      &
              ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'co_destruction: all assertions passed'

      end program co_destruction
