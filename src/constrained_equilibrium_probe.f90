      program constrained_equilibrium_probe
      ! Reevaluate a saved constrained-equilibrium cell without hydrodynamics.
      ! EXHALE_CCE_DUMP=<file> writes the first failed full-field continuation;
      ! this program regenerates rate coefficients at the saved state and
      ! reports residual terms, finite-difference Jacobians, and conditioning.
      ! Build with make cce_probe; run ./cce_probe.x <dump file>.

      use constrained_chemical_equilibrium, only: cce_probe_from_dump
      use mol_rates, only: h2_thermochemistry_init

      implicit none
      character(len=512) :: fname
      integer :: nargs

      nargs = command_argument_count()
      if (nargs .lt. 1) then
         write(*,'(A)') ' usage: cce_probe.x <dump file>'
         write(*,'(A)') '   the dump is written by a run with'//          &
                        ' EXHALE_CCE_DUMP set to that file name'
         stop
      endif
      call get_command_argument(1, fname)

      ! The H2 equilibrium-constant table the molecular rates read is built
      ! once, serially, before any of them is evaluated, as the main program
      ! does before its first sweep.
      call h2_thermochemistry_init
      call cce_probe_from_dump(trim(fname))

      end program constrained_equilibrium_probe
