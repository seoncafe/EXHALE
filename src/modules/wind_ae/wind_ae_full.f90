      program wind_ae_full
      ! Full P1 driver: relaxation (base->sonic) + outward Bulirsch-Stoer
      ! integration (sonic->Rmax), writing the complete windsoln
      ! (M+ADDPTS rows) for diffing against the C windsoln_ref.
      !
      ! Usage: wind_ae_full <guess.inp> <out.csv> <spectrum.inp> <inputs/> [rhoscale]
      use wae_config,     only: nsp => wae_nspecies, m => wae_m,          &
                                ne => wae_ne, nb => wae_nb,               &
                                itmax => wae_itmax, conv => wae_conv,     &
                                slowc => wae_slowc, addpts => wae_addpts, &
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
      use wae_intode,     only: wae_integrate_ode
      implicit none

      character(len=4096) :: gfile, outcsv, specf, indir, line, rhoarg
      integer :: u, uo, ios, row, j, i, ntot, tp
      real*8  :: y(ne, m), s(ne, 2*ne+1)
      real*8, allocatable :: c(:,:,:)
      integer :: indexv(ne)
      real*8  :: scalv(ne), cval(10)
      real*8, allocatable :: sr(:), srho(:), sv(:), sT(:), sq(:), sz(:)
      real*8, allocatable :: sYs(:,:), sNcol(:,:)

      if (command_argument_count() .lt. 4) then
         write(*,*) 'usage: wind_ae_full <guess> <out.csv> <spectrum.inp> <inputs/> [rhoscale]'
         stop 2
      end if
      call get_command_argument(1, gfile)
      call get_command_argument(2, outcsv)
      call get_command_argument(3, specf)
      call get_command_argument(4, indir)

      call wae_rate_coeffs_init()
      call rnew_init(); call rrec_init(); call fe_init(); call ion_pots_init()
      call cii_init(); call ciii_init(); call oii_init(); call oiii_init()
      call wae_load_spectrum(trim(specf))
      call wae_load_params(trim(indir)//'guess.inp')
      call wae_load_all_params(trim(indir))

      ! --- load guess first M rows ---
      open(newunit=u, file=trim(gfile), status='old', action='read')
      row = 0
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#' .or. len_trim(line) .eq. 0) cycle
         if (row .ge. m) exit
         call parse_csv(line, cval, 10)
         row = row + 1
         wae_x(row) = cval(9)
         y(1,row)=cval(3); y(2,row)=cval(10); y(3,row)=cval(2)
         y(4,row)=cval(4); y(5,row)=cval(5); y(6,row)=cval(6)
         y(7,row)=cval(7); y(8,row)=cval(8)
      end do
      close(u)
      if (row .lt. m) stop 'guess too short'

      ! --- setup_indices / setup_scales ---
      do i = 1, ne
         indexv(i) = i
      end do
      indexv(1) = 3 + 2*nsp; indexv(2) = 4 + 2*nsp
      indexv(3) = 1;         indexv(4) = 2 + nsp
      do j = 0, nsp-1
         indexv(5+j) = 2 + j; indexv(5+j+nsp) = 3 + nsp + j
      end do
      scalv(1)=VSCALE; scalv(2)=ZSCALE; scalv(4)=TEMPSCALE; scalv(3)=100.0d0
      if (command_argument_count() .ge. 5) then
         call get_command_argument(5, rhoarg); read(rhoarg,*) scalv(3)
      end if
      do j = 0, nsp-1
         scalv(5+j)=FPSCALE; scalv(5+nsp+j)=NCOLSCALE
      end do

      ! --- relax ---
      allocate(c(ne, ne-nb+1, m+1)); c = 0.0d0; s = 0.0d0
      write(*,*) '(wind_ae_full) relaxing...'
      call wae_solvde(itmax, conv, slowc, scalv, indexv, ne, nb, m, y, c, s)

      ! --- map relaxed solution into full-grid storage ---
      tp = m + addpts
      allocate(sr(tp), srho(tp), sv(tp), sT(tp), sq(tp), sz(tp))
      allocate(sYs(tp,nsp), sNcol(tp,nsp))
      do row = 1, m
         sq(row)   = wae_x(row)
         sz(row)   = y(2,row)
         sr(row)   = wae_x(row)*y(2,row) + par%Rmin
         srho(row) = y(3,row)
         sv(row)   = y(1,row)
         sT(row)   = y(4,row)
         do j = 1, nsp
            sYs(row,j)   = y(4+j,row)
            sNcol(row,j) = y(4+nsp+j,row)
         end do
      end do

      ! --- outward integration ---
      write(*,*) '(wind_ae_full) integrating outward...'
      call wae_integrate_ode(sr, srho, sv, sT, sYs, sNcol, sq, sz, ntot)
      write(*,'(A,I0,A)') ' (wind_ae_full) done, ', ntot, ' points.'

      ! --- write full windsoln ---
      open(newunit=uo, file=trim(outcsv), status='replace', action='write')
      write(uo,'(A)') '#vars: r,rho,v,T,Ys_HI,Ys_HeI,Ncol_HI,Ncol_HeI,q,z'
      do row = 1, ntot
         write(uo,'(9(ES24.17,A),ES24.17)')                              &
            sr(row),',',srho(row),',',sv(row),',',sT(row),',',           &
            sYs(row,1),',',sYs(row,2),',',sNcol(row,1),',',sNcol(row,2), &
            ',',sq(row),',',sz(row)
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
      end program wind_ae_full
