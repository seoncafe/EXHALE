      program ramp2_test
      ! self-consistent-BC ramp test: load a seed windsoln, ramp the system params to a target with
      ! static_bcs=.false. so the base BCs and sonic-point column density
      ! re-converge as the ramp moves far from the seed (the case that stalls
      ! the static-BC ramp). Reports success + final base BCs + a sanity check on the solution.
      ! Usage: ramp2_test <seed.csv> <spectrum.inp> <rhoscale> \
      !                   <Ftot> <Mp> <Rp> <Mstar> <a> <Lstar>
      use wae_config,      only: nsp => wae_nspecies, m => wae_m, ne => wae_ne
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
      use wae_continuation,only: load_seed, setup_indices_scales, ramp_to,  &
                                 cy, rhoscale
      implicit none
      character(len=4096) :: seed, specf, a(7)
      real*8  :: tgt(6)
      integer :: i, r, nbad
      logical :: okfin

      call get_command_argument(1, seed)
      call get_command_argument(2, specf)
      do i = 1, 7
         call get_command_argument(2+i, a(i))
      end do
      read(a(1),*) rhoscale
      do i = 1, 6
         read(a(1+i),*) tgt(i)
      end do

      call wae_rate_coeffs_init()
      call rnew_init(); call rrec_init(); call fe_init(); call ion_pots_init()
      call cii_init(); call ciii_init(); call oii_init(); call oiii_init()
      call wae_load_spectrum(trim(specf))
      call load_seed(trim(seed))
      call setup_indices_scales()

      write(*,'(A,6ES12.5)') ' seed Ftot/Mp/Rp/Ms/a/Ls = ',                 &
         par%Ftot, par%Mp, par%Rp, par%Mstar, par%semimajor, par%Lstar
      write(*,'(A,6ES12.5)') ' target                  = ', tgt
      write(*,'(A,3ES12.5)') ' seed base Rmin/rho/T     = ',                &
         par%Rmin, par%rho_rmin, par%T_rmin

      r = ramp_to(tgt(1), tgt(2), tgt(3), tgt(4), tgt(5), tgt(6),           &
                  static_bcs=.false.)

      write(*,'(A,I0)')      ' ramp_to result code = ', r
      write(*,'(A,6ES12.5)') ' final Ftot/Mp/Rp/Ms/a/Ls = ',               &
         par%Ftot, par%Mp, par%Rp, par%Mstar, par%semimajor, par%Lstar
      write(*,'(A,3ES12.5)') ' final base Rmin/rho/T     = ',              &
         par%Rmin, par%rho_rmin, par%T_rmin
      write(*,'(A,2ES12.5)') ' final Ncol_sp            = ', par%Ncol_sp
      ! sanity: NaN/sign check on the relaxed solution
      nbad = 0; okfin = .true.
      do i = 1, m
         if (cy(1,i) .le. 0.0d0) nbad = nbad + 1          ! v>0
         if (cy(3,i) .ne. cy(3,i)) okfin = .false.        ! rho NaN
         if (cy(4,i) .le. 0.0d0) okfin = .false.          ! T>0
      end do
      write(*,'(A,I0,A,L1)') ' relaxed soln: nonpos-v pts = ', nbad,        &
         ', finite/positive = ', okfin
      if (r .eq. 0 .and. okfin) then
         write(*,'(A)') ' RESULT: SELF-CONSISTENT-BC RAMP CONVERGED'
      else
         write(*,'(A)') ' RESULT: SELF-CONSISTENT-BC RAMP FAILED'
      end if
      end program ramp2_test
