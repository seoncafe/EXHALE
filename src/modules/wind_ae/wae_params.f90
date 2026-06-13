      module wae_params
      ! Reads wind-ae inputs/{planet_params,phys_params,bcs,term_ind,
      ! tech_params}.inp into a wae_paramlist, mirroring set_parameters()
      ! in io.c. H0 is derived as in io.c: H0 = CS0^2 Rp^2 /(G Mp).
      ! Faithful to the C field meanings; only the fields the ported
      ! physics needs are populated (enough for get_spQ/get_dYsdr and the
      ! residual equations).
      use wae_config, only: nsp => wae_nspecies, T0 => wae_T0
      use wae_types,  only: wae_paramlist
      implicit none

      ! constants (defs.h)
      real*8, parameter :: wae_G  = 6.67259d-8
      real*8, parameter :: wae_K  = 1.380658d-16
      real*8, parameter :: wae_MH = 1.6733d-24
      real*8, parameter :: wae_CS0 = 1.0d0   ! set in init (sqrt(K*T0/MH))

      type(wae_paramlist) :: wae_par
      real*8 :: wae_cs0_val = 0.0d0

      contains

      subroutine wae_load_all_params(dir)
      ! dir = path to the inputs/ directory (with trailing '/').
      character(len=*), intent(in) :: dir
      wae_cs0_val = sqrt(wae_K*T0/wae_MH)
      call read_planet(trim(dir)//'planet_params.inp')
      call read_phys(  trim(dir)//'phys_params.inp')
      call read_bcs(   trim(dir)//'bcs.inp')
      call read_term(  trim(dir)//'term_ind.inp')
      call read_tech(  trim(dir)//'tech_params.inp')
      ! derived scale height parameter H0 (io.c set_planet_params)
      wae_par%H0 = wae_cs0_val**2 * wae_par%Rp**2 / wae_G / wae_par%Mp
      end subroutine wae_load_all_params

      !--- helpers --------------------------------------------------!

      subroutine val_after_colon(line, s)
      ! text after ':' with any trailing '#...' comment stripped
      character(len=*), intent(in)  :: line
      character(len=*), intent(out) :: s
      integer :: a, h
      a = index(line, ':') + 1
      s = adjustl(line(a:))
      h = index(s, '#')
      if (h .gt. 0) s = s(1:h-1)
      end subroutine val_after_colon

      subroutine read_n_reals(s, vals, n)
      ! parse up to n comma-separated reals from s
      character(len=*), intent(in)  :: s
      integer,          intent(in)  :: n
      real*8,           intent(out) :: vals(n)
      integer :: p, j, b, L
      character(len=64) :: tok
      p = 1; L = len_trim(s)
      do j = 1, n
         b = index(s(p:L), ',')
         if (b .eq. 0) then
            tok = adjustl(s(p:L)); p = L + 1
         else
            tok = adjustl(s(p:p+b-2)); p = p + b
         end if
         read(tok,*) vals(j)
      end do
      end subroutine read_n_reals

      subroutine read_n_ints(s, vals, n)
      character(len=*), intent(in)  :: s
      integer,          intent(in)  :: n
      integer,          intent(out) :: vals(n)
      integer :: p, j, b, L
      character(len=64) :: tok
      p = 1; L = len_trim(s)
      do j = 1, n
         b = index(s(p:L), ',')
         if (b .eq. 0) then
            tok = adjustl(s(p:L)); p = L + 1
         else
            tok = adjustl(s(p:p+b-2)); p = p + b
         end if
         read(tok,*) vals(j)
      end do
      end subroutine read_n_ints

      !--- file readers ---------------------------------------------!

      subroutine read_planet(fname)
      character(len=*), intent(in) :: fname
      character(len=512) :: line, s
      integer :: u, ios
      open(newunit=u, file=fname, status='old', action='read')
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#' .or. len_trim(line) .eq. 0) cycle
         call val_after_colon(line, s)
         if      (index(line,'Mp:')        .eq. 1) then
            read(s,*) wae_par%Mp
         else if (index(line,'Ruv:')       .eq. 1) then
            read(s,*) wae_par%Rp
         else if (index(line,'Mstar:')     .eq. 1) then
            read(s,*) wae_par%Mstar
         else if (index(line,'semimajor:') .eq. 1) then
            read(s,*) wae_par%semimajor
         else if (index(line,'Ftot:')      .eq. 1) then
            read(s,*) wae_par%Ftot
         else if (index(line,'Lstar:')     .eq. 1) then
            read(s,*) wae_par%Lstar
         end if
      end do
      close(u)
      end subroutine read_planet

      subroutine read_phys(fname)
      character(len=*), intent(in) :: fname
      character(len=512) :: line, s
      integer :: u, ios, p, b, L, j
      character(len=32) :: tok
      open(newunit=u, file=fname, status='old', action='read')
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#' .or. len_trim(line) .eq. 0) cycle
         call val_after_colon(line, s)
         if      (index(line,'atomic_mass:') .eq. 1) then
            call read_n_reals(s, wae_par%atomic_mass, nsp)
         else if (index(line,'HX:')          .eq. 1) then
            call read_n_reals(s, wae_par%HX, nsp)
         else if (index(line,'Z:')           .eq. 1) then
            call read_n_ints(s, wae_par%Z, nsp)
         else if (index(line,'Ne:')          .eq. 1) then
            call read_n_ints(s, wae_par%N_e, nsp)
         else if (index(line,'Molec_adjust:').eq. 1) then
            read(s,*) wae_par%molec_adjust
         else if (index(line,'species_name:').eq. 1) then
            p = 1; L = len_trim(s)
            do j = 1, nsp
               b = index(s(p:L), ',')
               if (b .eq. 0) then
                  tok = adjustl(s(p:L)); p = L + 1
               else
                  tok = adjustl(s(p:p+b-2)); p = p + b
               end if
               wae_par%species(j) = trim(tok)
            end do
         end if
      end do
      close(u)
      end subroutine read_phys

      subroutine read_bcs(fname)
      character(len=*), intent(in) :: fname
      character(len=512) :: line, s
      integer :: u, ios
      open(newunit=u, file=fname, status='old', action='read')
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#' .or. len_trim(line) .eq. 0) cycle
         call val_after_colon(line, s)
         if      (index(line,'Rmin:')     .eq. 1) then
            read(s,*) wae_par%Rmin
         else if (index(line,'Rmax:')     .eq. 1) then
            read(s,*) wae_par%Rmax
         else if (index(line,'rho_rmin:') .eq. 1) then
            read(s,*) wae_par%rho_rmin
         else if (index(line,'T_rmin:')   .eq. 1) then
            read(s,*) wae_par%T_rmin
         else if (index(line,'Ys_rmin:')  .eq. 1) then
            call read_n_reals(s, wae_par%Ys_rmin, nsp)
         else if (index(line,'Ncol_sp:')  .eq. 1) then
            call read_n_reals(s, wae_par%Ncol_sp, nsp)
         else if (index(line,'erf_drop:') .eq. 1) then
            call read_n_reals(s, wae_par%erf_drop, 2)
         end if
      end do
      close(u)
      end subroutine read_bcs

      subroutine read_term(fname)
      character(len=*), intent(in) :: fname
      character(len=512) :: line, s
      integer :: u, ios
      open(newunit=u, file=fname, status='old', action='read')
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#' .or. len_trim(line) .eq. 0) cycle
         call val_after_colon(line, s)
         if      (index(line,'lyacool:')        .eq. 1) then
            read(s,*) wae_par%lyacool
         else if (index(line,'tidalforce:')     .eq. 1) then
            read(s,*) wae_par%tidalforce
         else if (index(line,'bolo_heat_cool:') .eq. 1) then
            read(s,*) wae_par%bolo_heat_cool
         else if (index(line,'integrate_outward:') .eq. 1) then
            read(s,*) wae_par%integrate_outward
         end if
      end do
      close(u)
      end subroutine read_term

      subroutine read_tech(fname)
      character(len=*), intent(in) :: fname
      character(len=512) :: line, s
      integer :: u, ios
      open(newunit=u, file=fname, status='old', action='read')
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#' .or. len_trim(line) .eq. 0) cycle
         call val_after_colon(line, s)
         if      (index(line,'breezeparam:') .eq. 1) then
            read(s,*) wae_par%breezeparam
         else if (index(line,'rapidity:')    .eq. 1) then
            read(s,*) wae_par%rapidity
         else if (index(line,'erfn:')        .eq. 1) then
            read(s,*) wae_par%erfn
         else if (index(line,'mach_limit:')  .eq. 1) then
            read(s,*) wae_par%mach_limit
         end if
      end do
      close(u)
      end subroutine read_tech

      end module wae_params
