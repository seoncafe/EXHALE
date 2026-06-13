      program soe_trace_f
      ! Fortran counterpart of wind-ae src/soe_trace.c: the P1-e leaf gate
      ! for wae_get_spQ and wae_get_dYsdr. Loads spectrum + params, then
      ! evaluates both at the sampled I_EQNVARS states and writes the same
      ! CSV layout (spQ, dYsdr[0..nsp-1]) as the C oracle.
      !
      ! Usage: soe_trace_f <samples> <out.csv> <spectrum.inp> <inputs_dir/>
      use wae_config,     only: nsp => wae_nspecies
      use wae_types,      only: wae_i_eqnvars
      use wae_spectrum,   only: wae_load_spectrum, wae_load_params
      use wae_params,     only: wae_load_all_params
      use wae_rate_coeffs,only: wae_rate_coeffs_init
      use wae_rnew,       only: wae_rnew_init => rnew_init
      use wae_rrec,       only: wae_rrec_init => rrec_init
      use wae_fe,         only: wae_fe_init => fe_init
      use wae_ion_pots,   only: wae_ion_pots_init => ion_pots_init
      use wae_lc_cii,     only: cii_init
      use wae_lc_ciii,    only: ciii_init
      use wae_lc_oii,     only: oii_init
      use wae_lc_oiii,    only: oiii_init
      use wae_soe,        only: wae_get_spQ, wae_get_dYsdr
      implicit none

      character(len=4096) :: samples, outcsv, specf, indir, line
      integer :: u, uo, ios, j
      real*8  :: tmp(4 + 2*2), spQ
      real*8, allocatable :: dYsdr(:)
      type(wae_i_eqnvars) :: v

      if (command_argument_count() .lt. 4) then
         write(*,*) 'usage: soe_trace_f <samples> <out.csv> <spectrum.inp> <inputs_dir/>'
         stop 2
      end if
      call get_command_argument(1, samples)
      call get_command_argument(2, outcsv)
      call get_command_argument(3, specf)
      call get_command_argument(4, indir)

      ! data tables + spectrum + params
      call wae_rate_coeffs_init()
      call wae_rnew_init(); call wae_rrec_init(); call wae_fe_init()
      call wae_ion_pots_init()
      call cii_init(); call ciii_init(); call oii_init(); call oiii_init()
      call wae_load_spectrum(trim(specf))
      ! glq_rates uses wae_spectrum's HX/atomic_mass/Ftot copies: fill them
      ! from guess.inp (same values as inputs/*.inp).
      call wae_load_params(trim(indir)//'guess.inp')
      ! soe uses wae_par (full paramlist): read all inputs/*.inp.
      call wae_load_all_params(trim(indir))

      allocate(dYsdr(nsp))
      open(newunit=u,  file=trim(samples), status='old', action='read')
      open(newunit=uo, file=trim(outcsv),  status='replace', action='write')
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (len_trim(line) .eq. 0 .or. line(1:1) .eq. '#') cycle
         read(line,*,iostat=ios) tmp(1:4+2*nsp)
         if (ios .ne. 0) cycle
         v%q = tmp(1); v%v = tmp(2); v%z = tmp(3); v%rho = tmp(4); v%T = tmp(5)
         do j = 1, nsp
            v%Ys(j)   = tmp(5+j)
            v%Ncol(j) = tmp(5+nsp+j)
         end do
         call wae_get_spQ(spQ, v, 1.0d4)
         call wae_get_dYsdr(dYsdr, v, 1.0d4)
         write(uo,'(ES24.17)',advance='no') spQ
         do j = 1, nsp
            write(uo,'(A,ES24.17)',advance='no') ',', dYsdr(j)
         end do
         write(uo,'(A)') ''
      end do
      close(u); close(uo)
      end program soe_trace_f
