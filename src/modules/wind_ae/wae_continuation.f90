      module wae_continuation
      ! C-1 continuation driver (static BCs): ramp the system parameters
      ! (Ftot, gravity, star) from a shipped seed windsoln to a target,
      ! re-relaxing at each adaptive step. Mirrors the Python wrapper's
      ! ramp_var/ramp_to with static_bcs=True (the base BCs stay at the
      ! seed values). In-memory guess hand-off replaces the C
      ! windsoln<->guess.inp file round-trip.
      use wae_config,  only: nsp => wae_nspecies, m => wae_m, ne => wae_ne, &
                             nb => wae_nb, itmax => wae_itmax,             &
                             conv => wae_conv, slowc => wae_slowc,         &
                             VSCALE=>wae_vscale, ZSCALE=>wae_zscale,       &
                             TEMPSCALE=>wae_tempscale, FPSCALE=>wae_fpscale, &
                             NCOLSCALE=>wae_ncolscale
      use wae_params,  only: par => wae_par, wae_G, wae_cs0_val, wae_K, wae_MH
      use wae_config,  only: T0 => wae_T0
      use wae_spectrum,only: wae_HX, wae_atomic_mass, wae_Ftot
      use wae_grid,    only: wae_x
      use wae_nr,      only: wae_solvde
      use wae_status,  only: wae_err
      implicit none

      ! relaxation work arrays + indices/scales (set once)
      real*8  :: cy(ne, m)            ! current relaxation solution (the guess)
      integer :: cv(ne)               ! indexv
      real*8  :: sv(ne)               ! scalv
      real*8  :: rhoscale = 100.0d0
      contains

      subroutine setup_indices_scales()
      integer :: i, j
      do i = 1, ne
         cv(i) = i
      end do
      cv(1) = 3 + 2*nsp; cv(2) = 4 + 2*nsp; cv(3) = 1; cv(4) = 2 + nsp
      do j = 0, nsp-1
         cv(5+j) = 2 + j; cv(5+j+nsp) = 3 + nsp + j
      end do
      sv(1) = VSCALE; sv(2) = ZSCALE; sv(3) = rhoscale; sv(4) = TEMPSCALE
      do j = 0, nsp-1
         sv(5+j) = FPSCALE; sv(5+nsp+j) = NCOLSCALE
      end do
      end subroutine setup_indices_scales

      integer function solve()
      ! relax cy in place from itself as the guess; return wae_err.
      real*8, allocatable :: c(:,:,:)
      real*8 :: s(ne, 2*ne+1)
      allocate(c(ne, ne-nb+1, m+1)); c = 0.0d0; s = 0.0d0
      wae_err = 0
      call wae_solvde(itmax, conv, slowc, sv, cv, ne, nb, m, cy, c, s)
      deallocate(c)
      solve = wae_err
      end function solve

      subroutine apply_var(which, val)
      ! set parameter `which` to val and refresh derived quantities
      integer, intent(in) :: which
      real*8,  intent(in) :: val
      select case (which)
      case (1); par%Ftot = val; wae_Ftot = val      ! glq uses wae_Ftot copy
      case (2); par%Mp = val
      case (3); par%Rp = val
      case (4); par%Mstar = val
      case (5); par%semimajor = val
      case (6); par%Lstar = val
      end select
      ! H0 = CS0^2 Rp^2 /(G Mp) depends on Mp, Rp
      par%H0 = wae_cs0_val**2 * par%Rp**2 / wae_G / par%Mp
      end subroutine apply_var

      real*8 function get_var(which)
      integer, intent(in) :: which
      select case (which)
      case (1); get_var = par%Ftot
      case (2); get_var = par%Mp
      case (3); get_var = par%Rp
      case (4); get_var = par%Mstar
      case (5); get_var = par%semimajor
      case (6); get_var = par%Lstar
      end select
      end function get_var

      integer function ramp_var(which, target, label)
      ! adaptive multiplicative ramp of one system parameter (static BCs).
      ! returns 0 on success, 101 on failure.
      integer, intent(in) :: which
      real*8,  intent(in) :: target
      character(len=*), intent(in) :: label
      real*8  :: cur, trial, delta, flip
      real*8  :: ysave(ne, m)
      integer :: failed, ierr
      cur = get_var(which)
      ramp_var = 0
      if (target .eq. 0.0d0 .or. abs(cur-target)/abs(target) .lt. 1.0d-10) then
         write(*,'(A,A,A)') '  ', trim(label), ' already at target.'
         return
      end if
      flip = 1.0d0
      if (cur .gt. target) flip = -1.0d0
      delta = 0.02d0*flip
      failed = 0
      do while (abs(cur-target)/abs(target) .gt. 1.0d-10)
         trial = cur*(1.0d0 + delta)
         if (flip*trial .gt. flip*target) trial = target
         ysave = cy
         call apply_var(which, trial)
         ierr = solve()
         if (ierr .ne. 0) then
            cy = ysave                 ! restore guess
            call apply_var(which, cur) ! restore param
            delta = delta/2.0d0
            failed = failed + 1
            if (failed .gt. 60) then
               write(*,'(A,A,A,ES12.5)') '  RAMP FAIL ', trim(label),    &
                  ' stuck at ', cur
               ramp_var = 101; return
            end if
         else
            cur = trial
            if (.not. (flip .lt. 0.0d0 .and. delta .le. -0.5d0)) delta = delta*2.0d0
            failed = 0
         end if
      end do
      write(*,'(A,A,A,ES12.5)') '  ramped ', trim(label), ' -> ', cur
      end function ramp_var

      integer function ramp_to(Ftot_t, Mp_t, Rp_t, Mstar_t, a_t, Lstar_t)
      ! ramp order matches the wrapper: Ftot -> gravity (Rp,Mp) ->
      ! star (Mstar, semimajor, Lstar).
      real*8, intent(in) :: Ftot_t, Mp_t, Rp_t, Mstar_t, a_t, Lstar_t
      integer :: r
      ramp_to = 0
      r = ramp_var(1, Ftot_t,  'Ftot');      if (r .ne. 0) then; ramp_to=r; return; end if
      r = ramp_var(3, Rp_t,    'Rp');        if (r .ne. 0) then; ramp_to=r; return; end if
      r = ramp_var(2, Mp_t,    'Mp');        if (r .ne. 0) then; ramp_to=r; return; end if
      r = ramp_var(4, Mstar_t, 'Mstar');     if (r .ne. 0) then; ramp_to=r; return; end if
      r = ramp_var(5, a_t,     'semimajor'); if (r .ne. 0) then; ramp_to=r; return; end if
      if (Lstar_t .gt. 0.0d0) then
         r = ramp_var(6, Lstar_t, 'Lstar');  if (r .ne. 0) then; ramp_to=r; return; end if
      end if
      end function ramp_to

      subroutine load_seed(wsfile)
      ! parse a windsoln header (#plnt_prms/#bcs/#phys_prms/#tech/#flags)
      ! into wae_par + wae_spectrum copies, and load the first M data rows
      ! into cy / wae_x. BCs stay fixed at these seed values (static_bcs).
      character(len=*), intent(in) :: wsfile
      character(len=8192) :: line
      integer :: u, ios, row, j
      real*8  :: pp(6), bb(4+2*nsp+2), tc(4), fl(4), cval(10)
      wae_cs0_val = sqrt(wae_K*T0/wae_MH)   ! CS0 (needed for H0)
      open(newunit=u, file=wsfile, status='old', action='read')
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (index(line,'#plnt_prms:') .gt. 0) then
            call csv_after(line, pp, 6)
            par%Mp=pp(1); par%Rp=pp(2); par%Mstar=pp(3)
            par%semimajor=pp(4); par%Ftot=pp(5); par%Lstar=pp(6)
         else if (index(line,'#bcs:') .gt. 0) then
            call csv_after(line, bb, 4+2*nsp+2)
            par%Rmin=bb(1); par%Rmax=bb(2); par%rho_rmin=bb(3); par%T_rmin=bb(4)
            do j = 1, nsp
               par%Ys_rmin(j) = bb(4+j); par%Ncol_sp(j) = bb(4+nsp+j)
            end do
            par%erf_drop(1)=bb(4+2*nsp+1); par%erf_drop(2)=bb(4+2*nsp+2)
         else if (index(line,'#phys_prms:') .gt. 0) then
            call phys_after(line)
         else if (index(line,'#tech:') .gt. 0) then
            call csv_after(line, tc, 4)
            par%breezeparam=tc(1); par%rapidity=tc(2); par%erfn=tc(3); par%mach_limit=tc(4)
         else if (index(line,'#flags:') .gt. 0) then
            call csv_after(line, fl, 4)
            par%lyacool=nint(fl(1)); par%tidalforce=fl(2)
            par%bolo_heat_cool=fl(3); par%integrate_outward=nint(fl(4))
         end if
      end do
      ! derived
      par%H0 = wae_cs0_val**2 * par%Rp**2 / wae_G / par%Mp
      wae_Ftot = par%Ftot
      wae_HX = par%HX
      wae_atomic_mass = par%atomic_mass
      ! data rows -> cy, wae_x (first M)
      rewind(u); row = 0
      do
         read(u,'(A)',iostat=ios) line
         if (ios .ne. 0) exit
         if (line(1:1) .eq. '#' .or. len_trim(line) .eq. 0) cycle
         if (row .ge. m) exit
         call csv_after_pos(line, cval, 10, 0)
         row = row + 1
         wae_x(row)=cval(9)
         cy(1,row)=cval(3); cy(2,row)=cval(10); cy(3,row)=cval(2); cy(4,row)=cval(4)
         cy(5,row)=cval(5); cy(6,row)=cval(6); cy(7,row)=cval(7); cy(8,row)=cval(8)
      end do
      close(u)
      end subroutine load_seed

      ! ---- parsing helpers ----
      subroutine csv_after(line, vals, n)
      character(len=*), intent(in) :: line
      integer, intent(in) :: n
      real*8, intent(out) :: vals(n)
      call csv_after_pos(line, vals, n, index(line,':'))
      end subroutine csv_after

      subroutine csv_after_pos(line, vals, n, after)
      character(len=*), intent(in) :: line
      integer, intent(in) :: n, after
      real*8, intent(out) :: vals(n)
      integer :: p, jj, b, L
      character(len=64) :: tok
      p = after + 1; L = len_trim(line)
      do jj = 1, n
         do while (p .le. L .and. line(p:p) .eq. ' ')
            p = p + 1
         end do
         b = index(line(p:L), ',')
         if (b .eq. 0) then
            tok = adjustl(line(p:L)); p = L + 1
         else
            tok = adjustl(line(p:p+b-2)); p = p + b
         end if
         read(tok,*) vals(jj)
      end do
      end subroutine csv_after_pos

      subroutine phys_after(line)
      ! #phys_prms: HX(1..nsp), name(1..nsp), amass(1..nsp), molec_adjust
      character(len=*), intent(in) :: line
      integer :: p, jj, b, L
      character(len=64) :: tok
      p = index(line,':') + 1; L = len_trim(line)
      do jj = 1, nsp
         call next_tok(line, p, L, tok); read(tok,*) par%HX(jj)
      end do
      do jj = 1, nsp
         call next_tok(line, p, L, tok); par%species(jj) = trim(tok)
      end do
      do jj = 1, nsp
         call next_tok(line, p, L, tok); read(tok,*) par%atomic_mass(jj)
      end do
      call next_tok(line, p, L, tok); read(tok,*) par%molec_adjust
      ! Z/N_e are not in the windsoln header; assume H/He
      par%Z(1)=1; par%Z(2)=2; par%N_e(1)=1; par%N_e(2)=2
      end subroutine phys_after

      subroutine next_tok(line, p, L, tok)
      character(len=*), intent(in) :: line
      integer, intent(inout) :: p
      integer, intent(in) :: L
      character(len=*), intent(out) :: tok
      integer :: b
      do while (p .le. L .and. line(p:p) .eq. ' ')
         p = p + 1
      end do
      b = index(line(p:L), ',')
      if (b .eq. 0) then
         tok = adjustl(line(p:L)); p = L + 1
      else
         tok = adjustl(line(p:p+b-2)); p = p + b
      end if
      end subroutine next_tok

      end module wae_continuation
