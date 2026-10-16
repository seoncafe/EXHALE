      program ic_writer_test
      ! Standalone driver for wae_ic_writer: read a windsoln.csv
      ! (dimensionalized via its #scales line) + an EXHALE IC_dump grid,
      ! write the EXHALE IC files. Mirrors windae_to_exhale_ic.py args
      ! (plus the three Roche-geometry values wae_write_ic gained for its
      ! base blend, which wind_ae_ic passes as Mstar/Mp, a/Rp and the
      ! tidal-force switch) so the output can be diffed against the
      ! Python converter.
      !
      ! Usage: ic_writer_test <windsoln.csv> <IC_dump.txt> <outdir> \
      !                       <lognbase> <T0> <Rp[cm]> <HeH> \
      !                       <Mstar/Mp> <a/Rp> <tidalforce>
      use wae_ic_writer, only: wae_write_ic
      implicit none
      character(len=4096) :: wfile, gfile, outdir, a4, a5, a6, a7, line
      character(len=4096) :: a8, a9, a10
      integer :: u, ios, nw, ng, i
      real*8  :: lognbase, T0, Rp, HeH, Mrapp, atilde, tidalf
      real*8  :: scales(10)
      real*8, allocatable :: rw(:), rhow(:), vw(:), Tw(:), YsHIw(:), YsHeIw(:)
      real*8, allocatable :: rgrid(:), gtmp(:)
      real*8  :: cval(10)

      call get_command_argument(1, wfile)
      call get_command_argument(2, gfile)
      call get_command_argument(3, outdir)
      call get_command_argument(4, a4); read(a4,*) lognbase
      call get_command_argument(5, a5); read(a5,*) T0
      call get_command_argument(6, a6); read(a6,*) Rp
      call get_command_argument(7, a7); read(a7,*) HeH
      call get_command_argument(8, a8);  read(a8,*)  Mrapp
      call get_command_argument(9, a9);  read(a9,*)  atilde
      call get_command_argument(10,a10); read(a10,*) tidalf

      ! --- read windsoln: #scales then data ---
      scales = 1.0d0
      nw = 0
      open(newunit=u, file=trim(wfile), status='old', action='read')
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (index(line,'#scales:') .gt. 0) then
            call parse_csv_after_colon(line, scales, 8)
         else if (line(1:1) .ne. '#' .and. len_trim(line) .gt. 0) then
            nw = nw + 1
         end if
      end do
      allocate(rw(nw), rhow(nw), vw(nw), Tw(nw), YsHIw(nw), YsHeIw(nw))
      rewind(u)
      i = 0
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#' .or. len_trim(line) .eq. 0) cycle
         call parse_csv(line, cval, 8)
         i = i + 1
         ! cols r,rho,v,T,Ys_HI,Ys_HeI,Ncol_HI,Ncol_HeI ; dimensionalize
         rw(i)   = cval(1)*scales(1)
         rhow(i) = cval(2)*scales(2)
         vw(i)   = cval(3)*scales(3)
         Tw(i)   = cval(4)*scales(4)
         YsHIw(i)  = cval(5)*scales(5)
         YsHeIw(i) = cval(6)*scales(6)
      end do
      close(u)

      ! --- read grid (IC_dump: r in col 1) ---
      ng = 0
      open(newunit=u, file=trim(gfile), status='old', action='read')
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .ne. '#' .and. len_trim(line) .gt. 0) ng = ng + 1
      end do
      allocate(rgrid(ng), gtmp(5))
      rewind(u)
      i = 0
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#' .or. len_trim(line) .eq. 0) cycle
         i = i + 1
         read(line,*) rgrid(i)
      end do
      close(u)

      call wae_write_ic(rw, rhow, vw, Tw, YsHIw, YsHeIw, nw,             &
                        rgrid, ng, lognbase, T0, Rp, HeH,                &
                        Mrapp, atilde, tidalf, trim(outdir))
      write(*,'(A,I0,A,I0)') 'windsoln pts=', nw, ' grid pts=', ng

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
      subroutine parse_csv_after_colon(s_in, vals, n)
      character(len=*), intent(in) :: s_in
      integer, intent(in) :: n
      real*8, intent(out) :: vals(n)
      integer :: a
      a = index(s_in,':')+1
      call parse_csv(adjustl(s_in(a:)), vals, n)
      end subroutine parse_csv_after_colon
      end program ic_writer_test
