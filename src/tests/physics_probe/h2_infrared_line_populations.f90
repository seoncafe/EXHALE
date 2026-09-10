      program h2_infrared_line_populations
      ! The LTE level populations that weight the H2 quadrupole lines of the
      ! infrared coolant are the populations of the caloric equation of
      ! state's level ladder, drawn from the one partition function of the
      ! code.
      !
      ! Production routines exercised:
      !   molecular_infrared_init and h2_line_emission_lte of
      !     src/modules/lower_atmosphere/molecular_infrared_cooling.f90 --
      !     the optically thin LTE emission per H2 molecule, eq. (4) of that
      !     module, f_u = g_u exp(-T_u/T)/Q summed over the line list;
      !   h2_partition_function of src/modules/states/caloric_eos.f90 -- the
      !     one Boltzmann sum over the one H2 level ladder, which the
      !     equilibrium constants of mol_rates and the internal energy of the
      !     caloric equation of state are also built from.
      !
      ! REFERENCE.  The denominator of f_u must be that Q and nothing else:
      ! a coolant whose populations came from a second sum would radiate the
      ! energy of a molecule the equation of state does not carry.  The
      ! emission is therefore reconstructed here, line by line, from
      ! h2_partition_function and the line list of molecular_infrared_data,
      ! and compared with what the module returns.  Tolerance 1e-12: the two
      ! are the same arithmetic, and all that separates them is the
      ! logarithm and exponential of the module's own log-log table lookup,
      ! which is exact at a node to rounding.
      !
      ! The temperatures are nodes of that table, t_tab(k) = 30 (8000/30)^
      ! ((k-1)/400) K, so the lookup returns the tabulated value itself and
      ! the assertion tests the populations rather than an interpolation.
      ! The five chosen nodes span 100 to 8000 K, the range over which the
      ! molecular layer radiates.
      !
      ! DIAGNOSTICS, not assertions:
      !   * the retired local sum.  Until B3b-IR the module recomputed the
      !     same ladder with an exp(-min(E/T, 700)) underflow guard.  Its
      !     ratio to h2_partition_function is printed at every node: the
      !     guard binds only where some level has E/T > 700, i.e. below
      !     51966/700 = 74.2 K, and it then replaces a term that underflows
      !     to zero by one of order 1e-304, which cannot move a sum whose
      !     ground term is 1.  This is the transition record of what the
      !     de-duplication changed.
      !   * the table ceiling.  h2_line_emission_lte is clamped above the
      !     8000 K end of its grid; the direct sum at 10000 K is printed
      !     beside the clamped value so the size of that clamp is on record.
      use molecular_infrared_data
      use molecular_infrared_cooling, only: molecular_infrared_init,      &
                                            h2_line_emission_lte
      use caloric_eos, only: h2_partition_function
      use global_parameters, only: kb_erg
      use assertion_report
      implicit none

      ! Nodes of the module's log-spaced evaluation grid, k = 1..401 over
      ! 30 to 8000 K.  These five are the nodes nearest 100, 300, 1000,
      ! 3000 and 8000 K.
      integer, parameter :: nodes(5) = (/ 87, 166, 253, 331, 401 /)
      real*8, parameter :: t_lo = 3.0d1, t_hi = 8.0d3
      integer, parameter :: n_tab = 401
      real*8  :: tg, e_code, e_ref, q_shared, q_retired
      character(len=48) :: label
      integer :: i

      call molecular_infrared_init(1.0d3)

      do i = 1, size(nodes)
         tg = t_lo*(t_hi/t_lo)**(dble(nodes(i) - 1)/dble(n_tab - 1))
         e_code   = h2_line_emission_lte(tg)
         e_ref    = h2_line_emission_from_partition_function(tg)
         write(label,'(a,i0)') 'h2_line_emission_from_shared_Q_at_',      &
                               nint(tg)
         call check_relative(label, e_code, e_ref, 1.0d-12)

         q_shared  = h2_partition_function(tg)
         q_retired = retired_guarded_level_sum(tg)
         write(*,'(a,f8.1,a,es22.15,a,es9.2)')                            &
              '  DIAGNOSTIC T [K] = ', tg,                                &
              '  Q = ', q_shared,                                         &
              '  retired/Q - 1 = ', q_retired/q_shared - 1.0d0
      enddo

      write(*,'(a,es12.5,a,es12.5)')                                      &
           '  DIAGNOSTIC emission per molecule at 10000 K: clamped = ',   &
           h2_line_emission_lte(1.0d4),                                   &
           '  direct sum = ',                                             &
           h2_line_emission_from_partition_function(1.0d4)

      if (assertion_failures .gt. 0) then
         write(*,'(a,i0,a)') 'h2_infrared_line_populations: ',            &
              assertion_failures, ' assertion(s) failed'
         flush(6)
         stop 1
      endif
      write(*,'(a)')                                                      &
           'h2_infrared_line_populations: all assertions passed'

      contains

      real*8 function h2_line_emission_from_partition_function(tg)        &
           result(e)
      ! Equation (4) of molecular_infrared_cooling, rebuilt here from the
      ! shared partition function: the optically thin LTE power radiated by
      ! one H2 molecule, Sum_l f_u A_l dE_l k_B, with f_u the population of
      ! the upper level of line l in the ladder of Roueff et al. (2019).
      real*8, intent(in) :: tg
      real*8  :: z, fu
      integer :: l

      z = h2_partition_function(tg)
      e = 0.0d0
      do l = 1, n_h2_line
         if (h2_line_Tu(l)/tg .gt. 7.0d2) cycle
         fu = h2_line_gu(l)*exp(-h2_line_Tu(l)/tg)/z
         e  = e + fu*h2_line_A(l)*h2_line_dE(l)*kb_erg
      enddo

      end function h2_line_emission_from_partition_function

      real*8 function retired_guarded_level_sum(tg) result(z)
      ! The level sum molecular_infrared_cooling carried before B3b-IR, kept
      ! here so that no production path can reach it.  Same ladder and same
      ! weights as h2_partition_function; the only difference is the clamp
      ! of the exponent at 700, which returns exp(-700) where the shared sum
      ! underflows to zero.
      real*8, intent(in) :: tg
      integer :: l

      z = 0.0d0
      do l = 1, n_h2_lev
         z = z + h2_lev_g(l)*exp(-min(h2_lev_T(l)/tg, 7.0d2))
      enddo

      end function retired_guarded_level_sum

      end program h2_infrared_line_populations
