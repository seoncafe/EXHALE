      module wae_intode
      ! Outward (sonic point -> Rmax) IVP integration, ported from
      ! wind-ae intode.c + the NR Bulirsch-Stoer stack (odeint, bsstep,
      ! mmid, pzextr). The ODE system reuses the already-ported
      ! derivative functions. y-vector layout here: y(1)=rho, y(2)=v,
      ! y(3)=T, y(4..3+nsp)=Ys, y(4+nsp..3+2nsp)=Ncol  (NUM_EQNS=3+2nsp).
      use wae_config, only: nsp => wae_nspecies, m => wae_m,             &
                            addpts => wae_addpts
      use wae_types,  only: wae_i_eqnvars
      use wae_params, only: par => wae_par
      use wae_soe,    only: wae_get_dvdr, wae_get_drhodr, wae_get_dTdr,   &
                            wae_get_dYsdr, wae_get_dNcoldr
      implicit none

      integer, parameter :: NUMEQ = 3 + 2*nsp
      real*8,  parameter :: ACCURACY  = 1.0d-13
      real*8,  parameter :: STEP_GUESS = 1.0d-1
      real*8,  parameter :: STEP_MIN  = 0.0d0

      ! intode globals
      real*8 :: thez = 0.0d0

      ! pzextr / bsstep shared extrapolation tableau
      real*8, allocatable :: d_tab(:,:), x_tab(:)

      ! bsstep persistent state
      integer, parameter :: KMAXX = 8, IMAXX = KMAXX+1
      logical :: bs_first = .true.
      integer :: bs_kmax = 0, bs_kopt = 0
      real*8  :: bs_epsold = -1.0d0, bs_xnew = 0.0d0
      real*8  :: bs_a(0:IMAXX)
      real*8  :: bs_alf(0:KMAXX, 0:KMAXX)
      integer :: bs_nseq(0:IMAXX) = [0,2,4,6,8,10,12,14,16,18]
      contains

      !---------------------------------------------------------------!
      subroutine wae_integrate_ode(sol_r, sol_rho, sol_v, sol_T,        &
                                   sol_Ys, sol_Ncol, sol_q, sol_z, npts_out)
      ! sol_* are full-grid arrays (length >= M+ADDPTS); rows 1..M hold
      ! the relaxed solution, rows M+1..M+ADDPTS are filled here.
      real*8, intent(inout) :: sol_r(:), sol_rho(:), sol_v(:), sol_T(:)
      real*8, intent(inout) :: sol_Ys(:,:), sol_Ncol(:,:)
      real*8, intent(inout) :: sol_q(:), sol_z(:)
      integer, intent(out)  :: npts_out
      real*8  :: ystart(NUMEQ), rsp, start, end, delta
      integer :: i, j, k, nok, nbad

      ! init_odeint: critical point = row M
      ystart(1) = sol_rho(m)
      ystart(2) = sol_v(m)
      ystart(3) = sol_T(m)
      do j = 1, nsp
         ystart(3+j)     = sol_Ys(m,j)
         ystart(3+nsp+j) = sol_Ncol(m,j)
      end do
      thez = sol_z(m)
      bs_first = .true.; bs_epsold = -1.0d0   ! fresh BS state per run

      if (.not. allocated(d_tab)) allocate(d_tab(NUMEQ, KMAXX), x_tab(KMAXX))

      rsp = sol_r(m)
      if (par%Rmax - rsp .lt. 0.1d0*rsp) then
         delta = 0.1d0*rsp/addpts
      else
         delta = (par%Rmax - rsp)/addpts
      end if
      end = rsp
      do i = 0, addpts-1
         start = end
         end = start + delta
         call odeint(ystart, NUMEQ, start, end, ACCURACY, STEP_GUESS, STEP_MIN, nok, nbad)
         ! fix unphysical Ys
         do j = 1, nsp
            k = 3 + j
            if (ystart(k) .gt. 1.0d0) ystart(k) = 1.0d0
            if (ystart(k) .lt. 0.0d0) ystart(k) = 0.0d0
         end do
         ! set_odeint_soln at row M+1+i
         call store_soln(m+1+i, end, ystart, sol_r, sol_rho, sol_v,     &
                         sol_T, sol_Ys, sol_Ncol, sol_q, sol_z)
      end do
      npts_out = m + addpts
      end subroutine wae_integrate_ode

      subroutine store_soln(row, end, ystart, sol_r, sol_rho, sol_v,    &
                            sol_T, sol_Ys, sol_Ncol, sol_q, sol_z)
      integer, intent(in) :: row
      real*8, intent(in)  :: end, ystart(NUMEQ)
      real*8, intent(inout) :: sol_r(:), sol_rho(:), sol_v(:), sol_T(:)
      real*8, intent(inout) :: sol_Ys(:,:), sol_Ncol(:,:), sol_q(:), sol_z(:)
      integer :: j
      sol_rho(row) = ystart(1)
      sol_v(row)   = ystart(2)
      sol_T(row)   = ystart(3)
      do j = 1, nsp
         sol_Ys(row,j)   = ystart(3+j)
         sol_Ncol(row,j) = ystart(3+nsp+j)
      end do
      sol_r(row) = end
      sol_z(row) = sol_z(1)
      sol_q(row) = (sol_r(row) - sol_r(1))/sol_z(row)
      end subroutine store_soln

      !---------------------------------------------------------------!
      subroutine ode_system(xr, y, dydx)
      real*8, intent(in)  :: xr, y(NUMEQ)
      real*8, intent(out) :: dydx(NUMEQ)
      type(wae_i_eqnvars) :: vars
      real*8 :: temp(nsp)
      integer :: j
      vars%r = xr
      vars%rho = y(1); vars%v = y(2); vars%T = y(3)
      do j = 1, nsp
         vars%Ys(j)   = y(3+j)
         vars%Ncol(j) = y(3+nsp+j)
      end do
      vars%q = (xr - par%Rmin)/thez
      vars%z = thez
      call wae_get_dNcoldr(temp, vars)
      do j = 1, nsp
         dydx(3+nsp+j) = temp(j)
      end do
      call wae_get_dYsdr(temp, vars, 1.0d4)     ! must come second
      do j = 1, nsp
         dydx(3+j) = temp(j)
      end do
      call wae_get_dvdr(dydx(2), vars)
      call wae_get_drhodr(dydx(1), vars, dydx(2))
      call wae_get_dTdr(dydx(3), vars, dydx(1), temp, 1.0d4)
      end subroutine ode_system

      !---------------------------------------------------------------!
      subroutine odeint(ystart, nvar, x1, x2, eps, h1, hmin, nok, nbad)
      integer, intent(in)    :: nvar
      real*8,  intent(inout) :: ystart(nvar)
      real*8,  intent(in)    :: x1, x2, eps, h1, hmin
      integer, intent(out)   :: nok, nbad
      integer, parameter :: MAXSTP = 1000000
      real*8,  parameter :: TINY = 1.0d-30
      integer :: nstp, i
      real*8  :: xx, hnext, hdid, h, yscal(nvar), yy(nvar), dydx(nvar)

      xx = x1
      h = sign(h1, x2-x1)
      nok = 0; nbad = 0
      yy = ystart
      do nstp = 1, MAXSTP
         call ode_system(xx, yy, dydx)
         do i = 1, nvar
            yscal(i) = abs(yy(i)) + abs(dydx(i)*h) + TINY
         end do
         if ((xx+h-x2)*(xx+h-x1) .gt. 0.0d0) h = x2 - xx
         call bsstep(yy, dydx, nvar, xx, h, eps, yscal, hdid, hnext)
         if (hdid .eq. h) then
            nok = nok + 1
         else
            nbad = nbad + 1
         end if
         if ((xx-x2)*(x2-x1) .ge. 0.0d0) then
            ystart = yy
            return
         end if
         if (abs(hnext) .le. hmin) stop '(odeint) Step size too small'
         h = hnext
      end do
      stop '(odeint) Too many steps'
      end subroutine odeint

      !---------------------------------------------------------------!
      subroutine bsstep(y, dydx, nv, xx, htry, eps, yscal, hdid, hnext)
      integer, intent(in)    :: nv
      real*8,  intent(inout) :: y(nv), dydx(nv), xx
      real*8,  intent(in)    :: htry, eps, yscal(nv)
      real*8,  intent(out)   :: hdid, hnext
      real*8,  parameter :: SAFE1=0.25d0, SAFE2=0.7d0, REDMAX=1.0d-5
      real*8,  parameter :: REDMIN=0.7d0, TINY=1.0d-30, SCALMX=0.1d0
      integer :: i, iq, k, kk, km
      real*8  :: eps1, errmax, fact, h, red, scale, work, wrkmin, xest
      real*8  :: err(KMAXX), yerr(nv), ysav(nv), yseq(nv)
      logical :: reduct, exitflag
      red = 1.0d0; km = 0; scale = 1.0d0; k = 0

      if (eps .ne. bs_epsold) then
         hnext = -1.0d29; bs_xnew = -1.0d29
         eps1 = SAFE1*eps
         bs_a(1) = bs_nseq(1) + 1.0d0
         do k = 1, KMAXX
            bs_a(k+1) = bs_a(k) + bs_nseq(k+1)
         end do
         do iq = 2, KMAXX
            do k = 1, iq-1
               bs_alf(k,iq) = eps1**((bs_a(k+1)-bs_a(iq+1)) /            &
                              ((bs_a(iq+1)-bs_a(1)+1.0d0)*(2*k+1)))
            end do
         end do
         bs_epsold = eps
         do bs_kopt = 2, KMAXX-1
            if (bs_a(bs_kopt+1) .gt. bs_a(bs_kopt)*bs_alf(bs_kopt-1,bs_kopt)) exit
         end do
         bs_kmax = bs_kopt
      end if
      h = htry
      ysav = y
      if (xx .ne. bs_xnew .or. h .ne. hnext) then
         bs_first = .true.
         bs_kopt = bs_kmax
      end if
      reduct = .false.
      do
         exitflag = .false.
         do k = 1, bs_kmax
            bs_xnew = xx + h
            if (bs_xnew .eq. xx) stop '(bsstep) step size underflow'
            call mmid(ysav, dydx, nv, xx, h, bs_nseq(k), yseq)
            xest = (h/bs_nseq(k))**2
            call pzextr(k, xest, yseq, y, yerr, nv)
            if (k .ne. 1) then
               errmax = TINY
               do i = 1, nv
                  errmax = max(errmax, abs(yerr(i)/yscal(i)))
               end do
               errmax = errmax/eps
               km = k - 1
               err(km) = (errmax/SAFE1)**(1.0d0/(2*km+1))
            end if
            if (k .ne. 1 .and. (k .ge. bs_kopt-1 .or. bs_first)) then
               if (errmax .lt. 1.0d0) then
                  exitflag = .true.; exit
               end if
               if (k .eq. bs_kmax .or. k .eq. bs_kopt+1) then
                  red = SAFE2/err(km); exit
               else if (k .eq. bs_kopt .and. bs_alf(bs_kopt-1,bs_kopt) .lt. err(km)) then
                  red = 1.0d0/err(km); exit
               else if (bs_kopt .eq. bs_kmax .and. bs_alf(km,bs_kmax-1) .lt. err(km)) then
                  red = bs_alf(km,bs_kmax-1)*SAFE2/err(km); exit
               else if (bs_alf(km,bs_kopt) .lt. err(km)) then
                  red = bs_alf(km,bs_kopt-1)/err(km); exit
               end if
            end if
         end do
         if (exitflag) exit
         red = min(red, REDMIN)
         red = max(red, REDMAX)
         h = h*red
         reduct = .true.
      end do
      xx = bs_xnew
      hdid = h
      bs_first = .false.
      wrkmin = 1.0d35
      do kk = 1, km
         fact = max(err(kk), SCALMX)
         work = fact*bs_a(kk+1)
         if (work .lt. wrkmin) then
            scale = fact; wrkmin = work; bs_kopt = kk+1
         end if
      end do
      hnext = h/scale
      if (bs_kopt .ge. k .and. bs_kopt .ne. bs_kmax .and. .not. reduct) then
         fact = max(scale/bs_alf(bs_kopt-1,bs_kopt), SCALMX)
         if (bs_a(bs_kopt+1)*fact .le. wrkmin) then
            hnext = h/fact
            bs_kopt = bs_kopt + 1
         end if
      end if
      end subroutine bsstep

      !---------------------------------------------------------------!
      subroutine mmid(y, dydx, nvar, xs, htot, nstep, yout)
      integer, intent(in)  :: nvar, nstep
      real*8,  intent(in)  :: y(nvar), dydx(nvar), xs, htot
      real*8,  intent(out) :: yout(nvar)
      integer :: n, i
      real*8  :: xr, swap, h2, h, ym(nvar), yn(nvar)
      h = htot/nstep
      do i = 1, nvar
         ym(i) = y(i)
         yn(i) = y(i) + h*dydx(i)
      end do
      xr = xs + h
      call ode_system(xr, yn, yout)
      h2 = 2.0d0*h
      do n = 2, nstep
         do i = 1, nvar
            swap = ym(i) + h2*yout(i)
            ym(i) = yn(i)
            yn(i) = swap
         end do
         xr = xr + h
         call ode_system(xr, yn, yout)
      end do
      do i = 1, nvar
         yout(i) = 0.5d0*(ym(i) + yn(i) + h*yout(i))
      end do
      end subroutine mmid

      !---------------------------------------------------------------!
      subroutine pzextr(iest, xest, yest, yz, dy, nv)
      integer, intent(in)    :: iest, nv
      real*8,  intent(in)    :: xest, yest(nv)
      real*8,  intent(inout) :: yz(nv), dy(nv)
      integer :: k1, j
      real*8  :: q, f2, f1, delta, cc(nv)
      x_tab(iest) = xest
      do j = 1, nv
         dy(j) = yest(j); yz(j) = yest(j)
      end do
      if (iest .eq. 1) then
         do j = 1, nv
            d_tab(j,1) = yest(j)
         end do
      else
         do j = 1, nv
            cc(j) = yest(j)
         end do
         do k1 = 1, iest-1
            delta = 1.0d0/(x_tab(iest-k1) - xest)
            f1 = xest*delta
            f2 = x_tab(iest-k1)*delta
            do j = 1, nv
               q = d_tab(j,k1)
               d_tab(j,k1) = dy(j)
               delta = cc(j) - q
               dy(j) = f1*delta
               cc(j) = f2*delta
               yz(j) = yz(j) + dy(j)
            end do
         end do
         do j = 1, nv
            d_tab(j,iest) = dy(j)
         end do
      end if
      end subroutine pzextr

      end module wae_intode
