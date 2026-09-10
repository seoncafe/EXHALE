      program h3p_recombination_branching
      ! The two product channels of H3+ dissociative recombination, against
      ! the published total and branching they are built from.
      !
      ! Production routines exercised: rk_R6_H3p_dr_H2 and rk_R7_H3p_dr_3H
      ! of src/modules/lower_atmosphere/mol_rates.f90.  Together they are
      ! the H3+ sink of the molecular base, and so they set how much H3+
      ! survives to radiate through h3p_cooling.f90.
      !
      ! REFERENCE.
      !   Larsson, McCall & Orel (2008), Chem. Phys. Lett. 462, 145,
      !   p. 149, give the thermal rate constant of the Kokoouline and
      !   Greene calculation,
      !
      !       alpha(300 K) = (7.2 +- 1.1)e-8 cm^3 s^-1,
      !
      !   "in good agreement with the new storage ring results", and the
      !   three-body branching ratio 0.70 +- 0.07, "in very good agreement
      !   with the CRYRING storage ring results".  Neither of the two
      !   constants in the module is printed in that paper; each is the
      !   total split by the branching,
      !
      !       k(H2 + H) = 0.30 alpha(300 K),  k(3H) = 0.70 alpha(300 K),
      !
      !   carried in temperature by the index 0.65 of the earlier
      !   storage-ring fits, which Yelle (2004) Icarus 170, 167, Table 1
      !   prints as R16a, 2.9e-8 (300/Te)^0.65 (Sundstrom et al. 1994), and
      !   R16b, 8.6e-8 (300/Te)^0.65 (Datz et al. 1995).
      !
      ! WHAT IS ASSERTED.
      !   That the module returns that construction and nothing else: the
      !   300 K total is the published 7.2e-8, the branching at every
      !   temperature is the published 0.70, and each channel follows
      !   (300/Te)^0.65.  These are transcriptions, so the tolerance is
      !   round-off, 1e-12 relative.
      !
      !   The note lines carry the ratio to the older Sundstrom and Datz
      !   pair that Yelle's table and Frelikh & Murray-Clay (2026) Table 1
      !   both use.  Larsson's own reading of that pair, p. 149, is that
      !   "the early results obtained at CRYRING [23,24] and ASTRID [27],
      !   which gave results just above or at 1e-7 cm^3 s^-1, were slightly
      !   too high because of rotational excitations", his Ref. 24 being
      !   Sundstrom et al., so the ratio is a superseded value over the one
      !   that supersedes it, not a disagreement between equals.

      use mol_rates,         only: rk_R6_H3p_dr_H2, rk_R7_H3p_dr_3H
      use assertion_report

      implicit none

      ! Larsson et al. (2008) p. 149.
      real*8, parameter :: alpha_300K = 7.2d-8, branch_3H = 0.70d0
      ! The temperature index of the Sundstrom and Datz fits, as Yelle
      ! (2004) Table 1 prints it, and their two constants.
      real*8, parameter :: index_T = 0.65d0
      real*8, parameter :: k_sundstrom_300K = 2.9d-8, k_datz_300K = 8.6d-8

      ! The base, the H3+ layer above it, and the top of the molecular
      ! region: the whole range over which the extrapolation is used.
      real*8, parameter :: T_K(3) = (/ 300.0d0, 1000.0d0, 3000.0d0 /)

      real*8 :: k_H2_H, k_3H, total, shape, reference
      character(len=64) :: label
      integer :: j

      do j = 1, 3
         shape = (300.0d0/T_K(j))**index_T
         k_H2_H = rk_R6_H3p_dr_H2(T_K(j))
         k_3H   = rk_R7_H3p_dr_3H(T_K(j))
         total  = k_H2_H + k_3H

         reference = alpha_300K*shape
         write(label,'(a,i0,a)') 'h3p_dissociative_recombination_total_T', &
              int(T_K(j)), 'K'
         call check_relative(label, total, reference, 1.0d-12)

         reference = (1.0d0 - branch_3H)*alpha_300K*shape
         write(label,'(a,i0,a)') 'h3p_recombination_to_h2_and_h_T',        &
              int(T_K(j)), 'K'
         call check_relative(label, k_H2_H, reference, 1.0d-12)

         reference = branch_3H*alpha_300K*shape
         write(label,'(a,i0,a)') 'h3p_recombination_to_three_atoms_T',     &
              int(T_K(j)), 'K'
         call check_relative(label, k_3H, reference, 1.0d-12)

         write(label,'(a,i0,a)') 'h3p_three_body_branching_fraction_T',    &
              int(T_K(j)), 'K'
         call check_relative(label, k_3H/total, branch_3H, 1.0d-12)

         write(*,'(a,i0,a,f8.5)')                                          &
              'note: h3p_recombination_ratio_to_sundstrom_datz_total_T',   &
              int(T_K(j)), 'K = ',                                         &
              total/((k_sundstrom_300K + k_datz_300K)*shape)
      enddo

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'h3p_recombination_branching: ',              &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)') 'h3p_recombination_branching: all assertions passed'

      end program h3p_recombination_branching
