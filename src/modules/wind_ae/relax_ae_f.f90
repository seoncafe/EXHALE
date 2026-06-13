      program relax_ae_f
      ! P1 end-to-end relaxation driver (the relax() of wind-ae relax.c,
      ! minus the outward ODE integration). Loads inputs/guess.inp as the
      ! initial guess, runs the ported NR relaxation (wae_solvde over
      ! wae_difeq), and writes the relaxed solution (first M points) so it
      ! can be diffed against the windsoln_ref first M rows.
      !
      ! Usage: relax_ae_f <guess.inp> <out.csv> <spectrum.inp> <inputs_dir/> [rhoscale]
      use wae_config,     only: nsp => wae_nspecies, m => wae_m,          &
                                ne => wae_ne, nb => wae_nb,               &
                                itmax => wae_itmax, conv => wae_conv,     &
                                slowc => wae_slowc, T0 => wae_T0,         &
                                RHO0 => wae_RHO0,                         &
                                VSCALE=>wae_vscale, ZSCALE=>wae_zscale,   &
                                TEMPSCALE=>wae_tempscale,                 &
                                FPSCALE=>wae_fpscale, NCOLSCALE=>wae_ncolscale
      use wae_params,     only: par => wae_par, wae_load_all_params
      use wae_spectrum,   only: wae_load_spectrum, wae_load_params
      use wae_rate_coeffs,only: wae_rate_coeffs_init
      use wae_rnew,       only: rnew_init
      use wae_rrec,       only: rrec_init
      use wae_fe,         only: fe_init
      use wae_ion_pots,   only: ion_pots_init
      use wae_lc_cii,     only: cii_init
      use wae_lc_ciii,    only: ciii_init
      use wae_lc_oii,     only: oii_init
      use wae_lc_oiii,    only: oiii_init
      use wae_grid,       only: wae_x
      use wae_nr,         only: wae_solvde
      implicit none

      character(len=4096) :: gfile, outcsv, specf, indir, line, rhoarg
      integer :: u, uo, ios, row, j, i, c0
      real*8  :: y(ne, m), s(ne, 2*ne+1)
      real*8, allocatable :: c(:,:,:)
      integer :: indexv(ne)
      real*8  :: scalv(ne), cval(10), r, rad

      if (command_argument_count() .lt. 4) then
         write(*,*) 'usage: relax_ae_f <guess> <out.csv> <spectrum.inp> <inputs_dir/> [rhoscale]'
         stop 2
      end if
      call get_command_argument(1, gfile)
      call get_command_argument(2, outcsv)
      call get_command_argument(3, specf)
      call get_command_argument(4, indir)

      ! data tables + spectrum + params
      call wae_rate_coeffs_init()
      call rnew_init(); call rrec_init(); call fe_init(); call ion_pots_init()
      call cii_init(); call ciii_init(); call oii_init(); call oiii_init()
      call wae_load_spectrum(trim(specf))
      call wae_load_params(trim(indir)//'guess.inp')
      call wae_load_all_params(trim(indir))

      ! --- load guess.inp first M rows into y/x (set_initial_guess) ---
      open(newunit=u, file=trim(gfile), status='old', action='read')
      row = 0
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#' .or. len_trim(line) .eq. 0) cycle
         if (row .ge. m) exit
         call parse_csv(line, cval, 10)
         row = row + 1
         ! guess cols: r,rho,v,T,Ys_HI,Ys_HeI,Ncol_HI,Ncol_HeI,q,z
         wae_x(row) = cval(9)
         y(1,row) = cval(3); y(2,row) = cval(10); y(3,row) = cval(2)
         y(4,row) = cval(4); y(5,row) = cval(5); y(6,row) = cval(6)
         y(7,row) = cval(7); y(8,row) = cval(8)
      end do
      close(u)
      if (row .lt. m) then
         write(*,*) 'guess too short:', row, 'of', m; stop 4
      end if

      ! --- setup_indices (relax.c) ---
      do i = 1, ne
         indexv(i) = i
      end do
      indexv(1) = 3 + 2*nsp
      indexv(2) = 4 + 2*nsp
      indexv(3) = 1
      indexv(4) = 2 + nsp
      do j = 0, nsp-1
         indexv(5+j)     = 2 + j
         indexv(5+j+nsp) = 3 + nsp + j
      end do

      ! --- setup_scales (relax.c); RHOSCALE optionally overridden ---
      scalv(1) = VSCALE; scalv(2) = ZSCALE; scalv(4) = TEMPSCALE
      scalv(3) = 100.0d0
      if (command_argument_count() .ge. 5) then
         call get_command_argument(5, rhoarg); read(rhoarg,*) scalv(3)
      end if
      do j = 0, nsp-1
         scalv(5+j)     = FPSCALE
         scalv(5+nsp+j) = NCOLSCALE
      end do

      ! --- relax (init_relax zeros s; solvde) ---
      allocate(c(ne, ne-nb+1, m+1))
      c = 0.0d0
      s = 0.0d0
      write(*,*) '(relax_ae_f) Starting relaxation...'
      call wae_solvde(itmax, conv, slowc, scalv, indexv, ne, nb, m, y, c, s)
      write(*,*) '(relax_ae_f) Relaxation done.'

      ! --- write relaxed first-M solution (set_relax_soln mapping) ---
      open(newunit=uo, file=trim(outcsv), status='replace', action='write')
      write(uo,'(A)') '#vars: r,rho,v,T,Ys_HI,Ys_HeI,Ncol_HI,Ncol_HeI,q,z'
      do row = 1, m
         r = wae_x(row)*y(2,row) + par%Rmin
         write(uo,'(9(ES24.17,A),ES24.17)')                              &
            r, ',', y(3,row), ',', y(1,row), ',', y(4,row), ',',         &
            y(5,row), ',', y(6,row), ',', y(7,row), ',', y(8,row), ',',  &
            wae_x(row), ',', y(2,row)
      end do
      close(uo)

      contains
      subroutine parse_csv(s_in, vals, n)
      character(len=*), intent(in) :: s_in
      integer, intent(in) :: n
      real*8, intent(out) :: vals(n)
      integer :: p, jj, b, L
      character(len=64) :: tok
      p = 1; L = len_trim(s_in)
      do jj = 1, n
         b = index(s_in(p:L), ',')
         if (b .eq. 0) then
            tok = adjustl(s_in(p:L)); p = L + 1
         else
            tok = adjustl(s_in(p:p+b-2)); p = p + b
         end if
         read(tok,*) vals(jj)
      end do
      end subroutine parse_csv
      end program relax_ae_f
