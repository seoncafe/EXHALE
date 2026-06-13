      module wae_nr
      ! Numerical Recipes relaxation stack (solvde / pinvs / red / bksub),
      ! ported verbatim from wind-ae ext/nrecipes/{solvde,pinvs,red,bksub}.c.
      ! Index conventions, array shapes (y(ne,m), s(ne,2*ne+1),
      ! c(ne,ne-nb+1,m+1)) and the difeq call order are preserved exactly;
      ! Fortran's native 1-based indexing replaces the C 1-based emulation.
      ! No algorithmic changes (per the port plan).
      use wae_difeq, only: wae_difeq_eval
      use wae_status, only: wae_err
      implicit none
      contains

      !---------------------------------------------------------------!
      subroutine wae_solvde(itmax, conv, slowc, scalv, indexv,         &
                            ne, nb, m, y, c, s)
      integer, intent(in)    :: itmax, ne, nb, m
      real*8,  intent(in)    :: conv, slowc
      real*8,  intent(in)    :: scalv(ne)
      integer, intent(in)    :: indexv(ne)
      real*8,  intent(inout) :: y(ne, m)
      real*8,  intent(inout) :: c(ne, ne-nb+1, m+1)
      real*8,  intent(inout) :: s(ne, 2*ne+1)

      integer :: ic1, ic2, ic3, ic4, it, j, j1, j2, j3, j4, j5, j6
      integer :: j7, j8, j9, jc1, jcf, jv, k, k1, k2, km, kp, nvars
      integer :: kmax(ne)
      real*8  :: ermax(ne)
      real*8  :: err, errj, fac, vmax, vz

      k1 = 1
      k2 = m
      nvars = ne*m
      j1 = 1
      j2 = nb
      j3 = nb + 1
      j4 = ne
      j5 = j4 + j1
      j6 = j4 + j2
      j7 = j4 + j3
      j8 = j4 + j4
      j9 = j8 + j1
      ic1 = 1
      ic2 = ne - nb
      ic3 = ic2 + 1
      ic4 = ne
      jc1 = 1
      jcf = ic3
      wae_err = 0
      do it = 1, itmax
         k = k1
         call wae_difeq_eval(k, k1, k2, j9, ic3, ic4, indexv, ne, s, y)
         if (wae_err .ne. 0) return
         call wae_pinvs(ic3, ic4, j5, j9, jc1, k1, c, s, ne, nb, m)
         if (wae_err .ne. 0) return
         do k = k1+1, k2
            kp = k - 1
            call wae_difeq_eval(k, k1, k2, j9, ic1, ic4, indexv, ne, s, y)
            if (wae_err .ne. 0) return
            call wae_red(ic1, ic4, j1, j2, j3, j4, j9, ic3, jc1, jcf,  &
                         kp, c, s, ne, nb, m)
            call wae_pinvs(ic1, ic4, j3, j9, jc1, k, c, s, ne, nb, m)
            if (wae_err .ne. 0) return
         end do
         k = k2 + 1
         call wae_difeq_eval(k, k1, k2, j9, ic1, ic2, indexv, ne, s, y)
         if (wae_err .ne. 0) return
         call wae_red(ic1, ic2, j5, j6, j7, j8, j9, ic3, jc1, jcf,     &
                      k2, c, s, ne, nb, m)
         call wae_pinvs(ic1, ic2, j7, j9, jcf, k2+1, c, s, ne, nb, m)
         if (wae_err .ne. 0) return
         call wae_bksub(ne, nb, jcf, k1, k2, c, m)
         err = 0.0d0
         do j = 1, ne
            jv = indexv(j)
            errj = 0.0d0
            vmax = 0.0d0
            km = 0
            do k = k1, k2
               vz = abs(c(jv, 1, k))
               if (vz .gt. vmax) then
                  vmax = vz
                  km = k
               end if
               errj = errj + vz
            end do
            err = err + errj/scalv(j)
            ermax(j) = c(jv, 1, km)/scalv(j)
            kmax(j) = km
         end do
         err = err/nvars
         if (err .gt. slowc) then
            fac = slowc/err
         else
            fac = 1.0d0
         end if
         do j = 1, ne
            jv = indexv(j)
            do k = k1, k2
               y(j, k) = y(j, k) - fac*c(jv, 1, k)
            end do
         end do
         write(*,'(/,A7,2X,A9,2X,A9)') 'Iter.', 'Error', 'FAC'
         write(*,'(I7,2X,ES9.3,2X,ES9.3)') it, err, fac
         if (err .lt. conv) return
      end do
      wae_err = 4   ! too many iterations -- graceful fail for continuation
      end subroutine wae_solvde

      !---------------------------------------------------------------!
      subroutine wae_pinvs(ie1, ie2, je1, jsf, jc1, k, c, s, ne, nb, m)
      integer, intent(in)    :: ie1, ie2, je1, jsf, jc1, k, ne, nb, m
      real*8,  intent(inout) :: c(ne, ne-nb+1, m+1)
      real*8,  intent(inout) :: s(ne, 2*ne+1)

      integer :: js1, jpiv, jp, je2, jcoff, j, irow, ipiv, id, icoff, i
      integer, allocatable :: indxr(:)
      real*8,  allocatable :: pscl(:)
      real*8  :: pivinv, piv, dum, big

      allocate(indxr(ie1:ie2), pscl(ie1:ie2))
      jpiv = 0
      ipiv = 0
      je2 = je1 + ie2 - ie1
      js1 = je2 + 1
      do i = ie1, ie2
         big = 0.0d0
         do j = je1, je2
            if (abs(s(i,j)) .gt. big) big = abs(s(i,j))
         end do
         if (big .eq. 0.0d0) then
            wae_err = 1; deallocate(indxr, pscl); return
         end if
         pscl(i) = 1.0d0/big
         indxr(i) = 0
      end do
      do id = ie1, ie2
         piv = 0.0d0
         do i = ie1, ie2
            if (indxr(i) .eq. 0) then
               big = 0.0d0
               do j = je1, je2
                  if (abs(s(i,j)) .gt. big) then
                     jp = j
                     big = abs(s(i,j))
                  end if
               end do
               if (big*pscl(i) .gt. piv) then
                  ipiv = i
                  jpiv = jp
                  piv = big*pscl(i)
               end if
            end if
         end do
         if (s(ipiv,jpiv) .eq. 0.0d0) then
            wae_err = 1; deallocate(indxr, pscl); return
         end if
         indxr(ipiv) = jpiv
         pivinv = 1.0d0/s(ipiv,jpiv)
         do j = je1, jsf
            s(ipiv,j) = s(ipiv,j)*pivinv
         end do
         s(ipiv,jpiv) = 1.0d0
         do i = ie1, ie2
            if (indxr(i) .ne. jpiv) then
               if (s(i,jpiv) .ne. 0.0d0) then
                  dum = s(i,jpiv)
                  do j = je1, jsf
                     s(i,j) = s(i,j) - dum*s(ipiv,j)
                  end do
                  s(i,jpiv) = 0.0d0
               end if
            end if
         end do
      end do

      jcoff = jc1 - js1
      icoff = ie1 - je1
      do i = ie1, ie2
         irow = indxr(i) + icoff
         if (irow .lt. 1) stop '(wae_pinvs) irow < 1'
         do j = js1, jsf
            c(irow, j+jcoff, k) = s(i,j)
         end do
      end do
      deallocate(indxr, pscl)
      end subroutine wae_pinvs

      !---------------------------------------------------------------!
      subroutine wae_red(iz1, iz2, jz1, jz2, jm1, jm2, jmf, ic1,       &
                         jc1, jcf, kc, c, s, ne, nb, m)
      integer, intent(in)    :: iz1, iz2, jz1, jz2, jm1, jm2, jmf
      integer, intent(in)    :: ic1, jc1, jcf, kc, ne, nb, m
      real*8,  intent(inout) :: c(ne, ne-nb+1, m+1)
      real*8,  intent(inout) :: s(ne, 2*ne+1)

      integer :: loff, l, j, ic, i
      real*8  :: vx

      loff = jc1 - jm1
      ic = ic1
      do j = jz1, jz2
         do l = jm1, jm2
            vx = c(ic, l+loff, kc)
            do i = iz1, iz2
               s(i,l) = s(i,l) - s(i,j)*vx
            end do
         end do
         vx = c(ic, jcf, kc)
         do i = iz1, iz2
            s(i,jmf) = s(i,jmf) - s(i,j)*vx
         end do
         ic = ic + 1
      end do
      end subroutine wae_red

      !---------------------------------------------------------------!
      subroutine wae_bksub(ne, nb, jf, k1, k2, c, m)
      integer, intent(in)    :: ne, nb, jf, k1, k2, m
      real*8,  intent(inout) :: c(ne, ne-nb+1, m+1)

      integer :: nbf, im, kp, k, j, i
      real*8  :: xx

      nbf = ne - nb
      im = 1
      do k = k2, k1, -1
         if (k .eq. k1) im = nbf + 1
         kp = k + 1
         do j = 1, nbf
            xx = c(j, jf, kp)
            do i = im, ne
               c(i, jf, k) = c(i, jf, k) - c(i, j, k)*xx
            end do
         end do
      end do
      do k = k1, k2
         kp = k + 1
         do i = 1, nb
            c(i, 1, k) = c(i+nbf, jf, k)
         end do
         do i = 1, nbf
            c(i+nb, 1, k) = c(i, jf, kp)
         end do
      end do
      end subroutine wae_bksub

      end module wae_nr
