      program atomic_mass_and_radius_constants
      ! The mass of a helium atom in the units the composition table counts
      ! in, and the astronomical constants the readers of `Planet radius
      ! [R_J]` and `Orbital distance [AU]` multiply by: one constant, one
      ! definition, and every module that needs it reads that definition.
      !
      ! Production data exercised: bsp_mass (src/modules/init/species_table.f90),
      ! the astronomical constants of module global_parameters
      ! (src/modules/init/parameters.f90), the copies the Wind-AE front end
      ! holds (src/modules/wind_ae/wae_exhale_input.f90), and the helium mass
      ! of binary_element_diffusion, wae_ic_writer
      ! and lower_column.  The fourth copy, the background collision partners
      ! of the carrier transport, is asserted in src/tests/species_masses/,
      ! whose driver needs the linear-algebra library that module's closure
      ! pulls in.
      !
      ! REFERENCE 1, the helium mass.  bsp_mass is stated in units of the
      ! density normalization mu = 1.67353284e-24 g, which is the hydrogen
      ! ATOM (parameters.f90: m_p + m_e - 13.6 eV/c^2, CODATA 2018), not the
      ! atomic mass unit u.  The helium-4 atom weighs 6.6464790722e-24 g
      ! (4.002603254 u, AME2020, times the CODATA 2018 atomic mass unit), so
      ! the number the table must carry for every helium stage is
      !     m_He/m_H = 6.6464790722e-24 / 1.67353284e-24 = 3.9715259,
      ! and HeH+ carries that plus one hydrogen atom.  The reference is
      ! formed HERE from the two atom masses, so the assertion compares the
      ! table with the measurement and not with a ratio someone rounded: the
      ! tolerance is therefore 1e-9 relative, round-off for a quotient of
      ! two doubles, and a table carrying a rounded 3.9715 (6.5e-6 away) or
      ! a helium of 4.0 m_H (0.72 per cent away) both fail it.  The mass
      ! every helium stage carries is the mass of the same nucleus with its
      ! electrons, so the four stages and the metastable agree to round-off,
      ! and HeH+ closes on the sum of its two nuclei to round-off.
      !
      ! REFERENCE 1a, THE ATOMIC MASS UNIT AGAINST THE CODE'S MASS UNIT.
      ! Standard atomic weights are tabulated in u (species_table
      ! melem_A_u: 15.999 for oxygen), while every mass the code adds into
      ! rho is counted in hydrogen atoms.  One u is amu/mu = 0.99223573
      ! hydrogen atoms, so melem_A, the array calc_rho, calc_mmw,
      ! comp_mass_per_H, the element census closure and the element
      ! diffusion weigh a metal nucleus with, must be melem_A_u times that
      ! factor: a weight left in u overstates every metal mass by 0.78 per
      ! cent.  The factor is again formed here from the two masses, so the
      ! tolerance is 1e-9.  The oxygen-chemistry rows of bsp_mass are built
      ! from the same converted weights, which is what makes one oxygen
      ! nucleus weigh the same whether it is counted through a metal ion
      ! column or through a molecule; that identity is asserted at 1e-15.
      !
      ! REFERENCE 1b, THE COPIES OF THE HELIUM MASS.  Three further modules
      ! weigh a helium atom of their own: the binary element diffusion
      ! (whose two-component closure m_1 n_H + m_He n_He = rho is exact only
      ! while its helium mass is the table's), the Wind-AE IC writer (which
      ! converts a windsoln density to nuclei and
      ! back) and the analytic Koskinen 2022 lower column (whose mean
      ! molecular weight is (1 + (m_He/m_H) f_He)/n_part).  Each must weigh
      ! the same atom as bsp_mass, so these are identities and the tolerance
      ! is zero, with one exception: the Wind-AE tree is compiled into the
      ! standalone wind_ae_ic.x and cannot use the species table, so it
      ! forms the ratio from the two atom masses it hands the solver.  Those
      ! are the same two measured masses, so the assertion is an identity as
      ! well, made at 1e-15 relative to allow for the two quotients being
      ! folded independently.
      !
      ! REFERENCE 2, the Jupiter radius.  One input key names one planet:
      ! `Planet radius [R_J]` is read by input_read (which multiplies by
      ! global_parameters RJ) and by the Wind-AE front end (which multiplies
      ! by its own copy), so the two constants must be one number or a run
      ! carries two planets.  The adopted value is the IAU 2015 nominal
      ! equatorial radius R_J^N(eq) = 7.1492e9 cm (Prsa et al. 2016, AJ 152,
      ! 41, Table 1), against which transiting-planet radii are quoted; the
      ! masses are the IAU 2015 nominal values built from the nominal
      ! GM (1.8982e30 g for Jupiter, 1.98842e33 g for the Sun).  These are
      ! exact identities between two literals, so the tolerance is zero.
      !
      ! REFERENCE 3, the astronomical unit.  `Orbital distance [AU]` is read
      ! by input_read and by the Wind-AE front end, so, like the planet
      ! radius, the two constants must be one number.  The IAU 2012
      ! definition 1 au = 1.495978707e11 m is exact, so this is an identity
      ! between two literals and the tolerance is zero.
      use species_table,    only: bsp_mass, melem_A, melem_A_u,           &
                                  iel_C, iel_O
      use global_parameters, only: RJ, MJ, Msun, AU, mu, m_He_atom, amu,  &
                                   m_He_over_m_H, amu_over_m_H
      use wae_exhale_input, only: RJ_bridge => RJ, MJ_bridge => MJ,       &
                                  MSUN_bridge => MSUN, AU_bridge => AU
      use wae_ic_writer,    only: ICW_mHe_over_mH
      use binary_element_diffusion, only: bed_m_He => m_He_amu,           &
                                          bed_m_H  => m_H_amu
      use lower_column,     only: column_m_He_over_m_H => m_He_over_m_H
      use assertion_report
      implicit none

      ! bsp positions of the base species (species_table, bsp_fsp order):
      ! 1 H I, 2 H II, 3 He I, 4 He II, 5 He III, 6 He 2^3S, 10 HeH+.
      ! and, for the oxygen chemistry, 11 OH, 12 H2O, 13 CO.
      integer, parameter :: ib_HI = 1, ib_HeI = 3, ib_HeII = 4,           &
                            ib_HeIII = 5, ib_HeTR = 6, ib_HeHp = 10,      &
                            ib_OH = 11, ib_H2O = 12, ib_CO = 13
      real*8, parameter :: m_He_g = 6.6464790722d-24 ! helium-4 atom [g]
      real*8, parameter :: m_H_g  = 1.67353284d-24   ! hydrogen atom [g]
      real*8, parameter :: m_u_g  = 1.66053906660d-24! atomic mass unit [g]
      real*8, parameter :: mass_ratio_He_H = m_He_g/m_H_g
      real*8, parameter :: mass_ratio_u_H  = m_u_g/m_H_g

      call check_relative('helium_atom_mass_in_hydrogen_atom_units',      &
           bsp_mass(ib_HeI)/bsp_mass(ib_HI), mass_ratio_He_H, 1.0d-9)
      call check_relative('helium_atom_mass_of_global_parameters',        &
           m_He_atom/mu, mass_ratio_He_H, 1.0d-9)
      call check_relative('helium_mass_ratio_of_global_parameters',       &
           m_He_over_m_H, mass_ratio_He_H, 1.0d-9)

      call check_relative('atomic_mass_unit_in_hydrogen_atom_units',      &
           amu/mu, mass_ratio_u_H, 1.0d-9)
      call check_relative('atomic_mass_unit_ratio_of_global_parameters',  &
           amu_over_m_H, mass_ratio_u_H, 1.0d-9)
      call check_relative('oxygen_weight_converted_from_u_to_hydrogen',   &
           melem_A(iel_O)/melem_A_u(iel_O), mass_ratio_u_H, 1.0d-9)
      call check_relative('iron_weight_converted_from_u_to_hydrogen',     &
           melem_A(10)/melem_A_u(10), mass_ratio_u_H, 1.0d-9)
      call check_relative('hydroxyl_mass_closes_on_its_two_nuclei',       &
           bsp_mass(ib_OH), bsp_mass(ib_HI) + melem_A(iel_O), 1.0d-15)
      call check_relative('water_mass_closes_on_its_three_nuclei',        &
           bsp_mass(ib_H2O), 2.0d0*bsp_mass(ib_HI) + melem_A(iel_O),      &
           1.0d-15)
      call check_relative('carbon_monoxide_mass_closes_on_its_nuclei',    &
           bsp_mass(ib_CO), melem_A(iel_C) + melem_A(iel_O), 1.0d-15)
      call check_relative('helium_ii_mass_equals_helium_i',               &
           bsp_mass(ib_HeII), bsp_mass(ib_HeI), 1.0d-12)
      call check_relative('helium_iii_mass_equals_helium_i',              &
           bsp_mass(ib_HeIII), bsp_mass(ib_HeI), 1.0d-12)
      call check_relative('helium_metastable_mass_equals_helium_i',       &
           bsp_mass(ib_HeTR), bsp_mass(ib_HeI), 1.0d-12)
      call check_relative('heh_plus_mass_closes_on_its_two_nuclei',       &
           bsp_mass(ib_HeHp), bsp_mass(ib_HeI) + bsp_mass(ib_HI), 1.0d-15)

      call check_relative('jupiter_radius_iau_2015_nominal_equatorial',   &
           RJ, 7.1492d9, 0.0d0)
      call check_relative('jupiter_radius_of_wind_ae_front_end',          &
           RJ_bridge, RJ, 0.0d0)
      call check_relative('jupiter_mass_of_wind_ae_front_end',            &
           MJ_bridge, MJ, 0.0d0)
      call check_relative('solar_mass_of_wind_ae_front_end',              &
           MSUN_bridge, Msun, 0.0d0)
      call check_relative('astronomical_unit_of_wind_ae_front_end',       &
           AU_bridge, AU, 0.0d0)

      call check_relative('helium_mass_of_binary_element_diffusion',      &
           bed_m_He, bsp_mass(ib_HeI), 0.0d0)
      call check_relative('hydrogen_mass_of_binary_element_diffusion',    &
           bed_m_H, bsp_mass(ib_HI), 0.0d0)

      call check_relative('helium_mass_of_wind_ae_ic_writer',             &
           ICW_mHe_over_mH, bsp_mass(ib_HeI)/bsp_mass(ib_HI), 1.0d-15)
      call check_relative('helium_mass_of_analytic_lower_column',         &
           column_m_He_over_m_H, bsp_mass(ib_HeI)/bsp_mass(ib_HI), 0.0d0)
      call check_relative('lower_column_helium_mass_is_the_global',       &
           column_m_He_over_m_H, m_He_over_m_H, 0.0d0)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'atomic_mass_and_radius_constants: ',        &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)')                                                     &
           'atomic_mass_and_radius_constants: all assertions passed'

      end program atomic_mass_and_radius_constants
