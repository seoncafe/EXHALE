      program base_bcs_test
      ! C-2 Phase-1 validation: load a seed windsoln, compute base_bcs, print.
      ! Compare against relax_wrapper.base_bcs (Python oracle).
      ! Usage: base_bcs_test <seed.csv>
      use wae_params,       only: par => wae_par
      use wae_continuation, only: load_seed, base_bcs
      implicit none
      character(len=4096) :: seed
      real*8 :: R, Rmax, rho, T
      call get_command_argument(1, seed)
      call load_seed(trim(seed))
      write(*,'(A,4ES16.8)') ' params Lstar/a/Mp/Rp = ',                   &
         par%Lstar, par%semimajor, par%Mp, par%Rp
      call base_bcs(R, Rmax, rho, T)
      write(*,'(A,ES24.16)') ' Rmin     = ', R
      write(*,'(A,ES24.16)') ' rho_rmin = ', rho
      write(*,'(A,ES24.16)') ' T_rmin   = ', T
      write(*,'(A,ES24.16)') ' Rmax     = ', Rmax
      end program base_bcs_test
