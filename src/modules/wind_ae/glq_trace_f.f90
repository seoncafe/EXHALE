      program glq_trace_f
      ! Fortran counterpart of wind-ae src/glq_trace.c (the P1 unit gate).
      ! Loads inputs/spectrum.inp + inputs/guess.inp from the cwd, then
      ! evaluates wae_glq_rates_eval at the sampled (N,Ys) states in the
      ! samples file and writes the same CSV layout as the C oracle.
      !
      ! Usage:  glq_trace_f <samples.txt> <out.csv> <spectrum.inp> <guess.inp>
      use wae_config,   only: nsp => wae_nspecies
      use wae_spectrum, only: wae_load_spectrum, wae_load_params
      use wae_rate_coeffs, only: wae_rate_coeffs_init
      use wae_glq_rates,   only: wae_glq_rates_eval
      implicit none

      character(len=4096) :: samples, outcsv, specf, guessf, line
      integer :: u, uo, ios, j, nion
      real*8  :: N(2), Ys(2)
      real*8, allocatable :: ion(:), heat(:)

      if (command_argument_count() .lt. 4) then
         write(*,*) 'usage: glq_trace_f <samples> <out.csv> <spectrum.inp> <guess.inp>'
         stop 2
      end if
      call get_command_argument(1, samples)
      call get_command_argument(2, outcsv)
      call get_command_argument(3, specf)
      call get_command_argument(4, guessf)

      call wae_rate_coeffs_init()
      call wae_load_spectrum(trim(specf))
      call wae_load_params(trim(guessf))

      nion = nsp*(nsp+1)
      allocate(ion(nion), heat(nsp))

      open(newunit=u,  file=trim(samples), status='old', action='read')
      open(newunit=uo, file=trim(outcsv),  status='replace', action='write')
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (len_trim(line) .eq. 0 .or. line(1:1) .eq. '#') cycle
         read(line,*,iostat=ios) N(1), N(2), Ys(1), Ys(2)
         if (ios .ne. 0) cycle

         call wae_glq_rates_eval(N, Ys, ion, heat)

         do j = 1, nsp
            write(uo,'(ES24.17,A)',advance='no') N(j), ','
         end do
         do j = 1, nsp
            write(uo,'(ES24.17,A)',advance='no') Ys(j), ','
         end do
         do j = 1, nion
            write(uo,'(ES24.17,A)',advance='no') ion(j), ','
         end do
         do j = 1, nsp
            if (j .lt. nsp) then
               write(uo,'(ES24.17,A)',advance='no') heat(j), ','
            else
               write(uo,'(ES24.17)') heat(j)
            end if
         end do
      end do
      close(u); close(uo)
      end program glq_trace_f
