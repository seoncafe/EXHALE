      program carrier_background_masses
      ! The mass the carrier transport gives each of its background
      ! collision partners.
      !
      ! Production data exercised: carrier_background_mass of
      ! src/modules/lower_atmosphere/diffusive_photochemistry.f90 and
      ! bsp_mass of src/modules/init/species_table.f90.
      !
      ! The molecular carrier transport moves a carrier through a background
      ! of nine species by Blanc's law, and the coefficient of every pair is
      ! the rigid-sphere coefficient of Banks & Kockarts (1973),
      !     D = 1.52e18 (1/A_s + 1/A_t)^(1/2) sqrt(T)/n,
      ! whose A_s and A_t are the two masses in units of the hydrogen atom.
      ! The same masses set the gravitational settling drift.  They must be
      ! the masses calc_rho weighs those species with, i.e. the species
      ! table's, or the friction and the density describe different gases;
      ! the module holding its own list of them is how a helium of 4 m_H
      ! survived the table's move to 3.9715 (decision 14 of
      ! docs/development_plan_20260905_rev3.md section 10.5).  These are
      ! identities between one number and itself, so the tolerance is zero.
      !
      ! Verdict lines follow the convention of src/tests/physics_probe:
      !     PASS|FAIL <name> measured= reference= tol=
      use species_table, only: bsp_mass
      use diffusive_photochemistry, only: carrier_background_mass
      use assertion_report
      implicit none

      ! bsp positions of the base species (species_table, bsp_fsp order):
      ! 1 H I, 2 H II, 3 He I, 4 He II, 5 He III, 6 He 2^3S, 7 H2, 8 H2+,
      ! 9 H3+, 10 HeH+.
      integer, parameter :: ib_HI = 1, ib_HII = 2, ib_HeI = 3,            &
                            ib_HeII = 4, ib_HeIII = 5, ib_H2 = 7,         &
                            ib_H2p = 8, ib_H3p = 9, ib_HeHp = 10
      ! Blanc's-law background partners, in the order the operator lists
      ! them: 1 H I, 2 H II, 3 H2, 4 H2+, 5 H3+, 6 He I, 7 He II,
      ! 8 He III, 9 HeH+.
      integer, parameter :: n_partner = 9
      integer, parameter :: partner_bsp(n_partner) =                      &
           [ ib_HI, ib_HII, ib_H2, ib_H2p, ib_H3p,                        &
             ib_HeI, ib_HeII, ib_HeIII, ib_HeHp ]
      character(len=32), parameter :: partner_name(n_partner) =           &
           [ 'hydrogen_i_mass_of_background    ',                         &
             'hydrogen_ii_mass_of_background   ',                         &
             'h2_mass_of_background            ',                         &
             'h2_plus_mass_of_background       ',                         &
             'h3_plus_mass_of_background       ',                         &
             'helium_i_mass_of_background      ',                         &
             'helium_ii_mass_of_background     ',                         &
             'helium_iii_mass_of_background    ',                         &
             'heh_plus_mass_of_background      ' ]
      integer :: is

      do is = 1, n_partner
         call check_relative(trim(partner_name(is)),                      &
              carrier_background_mass(is), bsp_mass(partner_bsp(is)),     &
              0.0d0)
      enddo

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'carrier_background_masses: ',               &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'carrier_background_masses: all assertions passed'

      end program carrier_background_masses
