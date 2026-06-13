      program ncol_test
      ! C-2 Phase-2 validation: load a converged seed, compute the
      ! self-consistent sonic-point column density goal, print it.
      ! Compare against relax_wrapper.self_consistent_Ncol (Python oracle).
      ! Usage: ncol_test <seed.csv> <spectrum.inp>
      use wae_config,      only: nsp => wae_nspecies
      use wae_params,      only: par => wae_par
      use wae_spectrum,    only: wae_load_spectrum
      use wae_rate_coeffs, only: wae_rate_coeffs_init
      use wae_rnew,        only: rnew_init
      use wae_rrec,        only: rrec_init
      use wae_fe,          only: fe_init
      use wae_ion_pots,    only: ion_pots_init
      use wae_lc_cii,      only: cii_init
      use wae_lc_ciii,     only: ciii_init
      use wae_lc_oii,      only: oii_init
      use wae_lc_oiii,     only: oiii_init
      use wae_continuation,only: load_seed, setup_indices_scales,           &
                                 self_consistent_Ncol
      implicit none
      character(len=4096) :: seed, specf
      real*8 :: goal(nsp)
      integer :: j
      call get_command_argument(1, seed)
      call get_command_argument(2, specf)
      call wae_rate_coeffs_init()
      call rnew_init(); call rrec_init(); call fe_init(); call ion_pots_init()
      call cii_init(); call ciii_init(); call oii_init(); call oiii_init()
      call wae_load_spectrum(trim(specf))
      call load_seed(trim(seed))
      call setup_indices_scales()
      call self_consistent_Ncol(goal)
      do j = 1, nsp
         write(*,'(A,I0,A,ES24.16,A,ES16.8)') ' Ncol_sp goal(', j, ') = ',  &
            goal(j), '   seed BC = ', par%Ncol_sp(j)
      end do
      end program ncol_test
