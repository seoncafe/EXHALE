      program constrained_equilibrium_probe
      ! Standalone re-evaluation of one saved constrained-equilibrium state,
      ! without the hydrodynamics: phase A/B of
      ! docs/charge_exchange_cancellation_limit.md.
      !
      ! The state is produced by a normal run with the environment variable
      ! EXHALE_CCE_DUMP set to a file name; the first continuation that fails
      ! to reach the cell's full radiation field writes itself there. This
      ! program reads that file, rebuilds the cell by calling the SAME
      ! coefficient routines the equilibrium sweep calls (so the rate
      ! coefficients are regenerated from temperature and density rather than
      ! transcribed), and reports the residual decomposition, the Jacobian
      ! under several difference rules, and the conditioning.
      !
      ! Not part of the default build. Build it with
      !     make cce_probe
      ! and run it as
      !     ./cce_probe.x <dump file>
      !
      !   EXHALE_CCE_DUMP=/tmp/cell387.dump ./EXHALE.x
      !   ./cce_probe.x /tmp/cell387.dump

      use constrained_chemical_equilibrium, only: cce_probe_from_dump

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

      call cce_probe_from_dump(trim(fname))

      end program constrained_equilibrium_probe
