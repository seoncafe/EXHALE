      module wae_ic_writer
      ! Convert a Wind-AE windsoln onto the EXHALE grid and write
      ! output/{Hydro_ioniz,Ion_species}_IC.txt for EXHALE's `Load IC`.
      ! Faithful Fortran port of src/utils/windae_to_exhale_ic.py
      ! (validated 2026-06-13). np_interp() matches numpy.interp (linear,
      ! endpoint-clamped). Below Wind-AE's Rmin: hydrostatic blend to the
      ! EXHALE base anchor; above Rmax: flux-conserving extrapolation.
      implicit none
      real*8, parameter :: ICW_mH = 1.6726d-24
      real*8, parameter :: ICW_kB = 1.380649d-16
      contains

      real*8 function np_interp(xq, xp, fp, n)
      ! linear interpolation with endpoint clamping (numpy.interp), xp asc
      integer, intent(in) :: n
      real*8,  intent(in) :: xq, xp(n), fp(n)
      integer :: lo, hi, mid
      real*8  :: t
      if (xq .le. xp(1)) then
         np_interp = fp(1); return
      else if (xq .ge. xp(n)) then
         np_interp = fp(n); return
      end if
      lo = 1; hi = n
      do while (hi - lo .gt. 1)
         mid = (lo+hi)/2
         if (xp(mid) .le. xq) then
            lo = mid
         else
            hi = mid
         end if
      end do
      t = (xq - xp(lo))/(xp(hi) - xp(lo))
      np_interp = fp(lo) + t*(fp(hi) - fp(lo))
      end function np_interp

      subroutine wae_write_ic(rw, rhow, vw, Tw, YsHIw, YsHeIw, nw,        &
                              rgrid, ng, lognbase, T0, Rp, HeH, outdir)
      ! rw[cm]/rhow[g/cc]/vw[cm/s]/Tw[K]/Ys* : windsoln (length nw, r asc)
      ! rgrid : EXHALE grid radius [Rp] (length ng, incl. ghosts)
      integer, intent(in) :: nw, ng
      real*8,  intent(in) :: rw(nw), rhow(nw), vw(nw), Tw(nw)
      real*8,  intent(in) :: YsHIw(nw), YsHeIw(nw), rgrid(ng)
      real*8,  intent(in) :: lognbase, T0, Rp, HeH
      character(len=*), intent(in) :: outdir
      real*8  :: rwR(nw), nnuc_w(nw), lognnuc_w(nw), fHIw(nw), fHeIw(nw)
      real*8  :: lognnuc(ng), T(ng), v(ng), fHI(ng), fHeI(ng), r(ng)
      real*8  :: nbase, rmin_w, rmax_w, b0_eff, ln10, sblend, logw_at, hyd_at
      real*8  :: nn_rmin, v_rmin, T_rmin_w, fHI_rmin, fHeI_rmin
      real*8  :: nH, nHe, nHI, nHII, nHeI, nHeII, ne, ntot, p, nmass
      integer :: i, i1, u
      ln10 = log(10.0d0)
      nbase = 10.0d0**lognbase
      r = rgrid
      do i = 1, nw
         rwR(i) = rw(i)/Rp
         nnuc_w(i) = (rhow(i)/(ICW_mH*(1.0d0+4.0d0*HeH)))*(1.0d0+HeH)
         lognnuc_w(i) = log10(nnuc_w(i))
         fHIw(i)  = min(max(YsHIw(i),  0.0d0), 1.0d0)
         fHeIw(i) = min(max(YsHeIw(i), 0.0d0), 1.0d0)
      end do
      rmin_w = rwR(1); rmax_w = rwR(nw)

      ! mid-region interpolation
      do i = 1, ng
         lognnuc(i) = np_interp(r(i), rwR, lognnuc_w, nw)
         T(i)    = np_interp(r(i), rwR, Tw, nw)
         v(i)    = np_interp(r(i), rwR, vw, nw)
         fHI(i)  = np_interp(r(i), rwR, fHIw, nw)
         fHeI(i) = np_interp(r(i), rwR, fHeIw, nw)
      end do

      ! b0_eff from the windsoln scale height just above Rmin
      i1 = 1
      do i = 1, nw
         if (rwR(i) .ge. rmin_w + 0.1d0) then
            i1 = i; exit
         end if
      end do
      b0_eff = (log(nnuc_w(1)) - log(nnuc_w(i1))) /                      &
               (1.0d0/rwR(1) - 1.0d0/rwR(i1))
      logw_at = np_interp(rmin_w, rwR, lognnuc_w, nw)
      hyd_at  = lognbase + b0_eff*(1.0d0/rmin_w - 1.0d0)/ln10
      nn_rmin = 10.0d0**np_interp(rmin_w, rwR, lognnuc_w, nw)
      v_rmin  = np_interp(rmin_w, rwR, vw, nw)
      T_rmin_w  = np_interp(rmin_w, rwR, Tw, nw)
      fHI_rmin  = np_interp(rmin_w, rwR, fHIw, nw)
      fHeI_rmin = np_interp(rmin_w, rwR, fHeIw, nw)

      do i = 1, ng
         if (r(i) .le. rmin_w) then
            ! below Rmin: hydrostatic blend to EXHALE base anchor
            sblend = (1.0d0 - 1.0d0/r(i)) / (1.0d0 - 1.0d0/rmin_w)
            sblend = min(max(sblend, 0.0d0), 1.0d0)
            lognnuc(i) = (lognbase + b0_eff*(1.0d0/r(i) - 1.0d0)/ln10)   &
                         + (logw_at - hyd_at)*sblend
            T(i)    = T0 + (T_rmin_w - T0)*sblend
            fHI(i)  = 1.0d0 + (fHI_rmin  - 1.0d0)*sblend
            fHeI(i) = 1.0d0 + (fHeI_rmin - 1.0d0)*sblend
            v(i) = v_rmin*(nn_rmin*rmin_w**2)/(10.0d0**lognnuc(i)*r(i)**2)
         else if (r(i) .ge. rmax_w) then
            ! above Rmax: frozen T,v; flux-conserving density
            T(i) = Tw(nw); v(i) = vw(nw)
            lognnuc(i) = lognnuc_w(nw) + 2.0d0*log10(rmax_w) - 2.0d0*log10(r(i))
            fHI(i) = fHIw(nw); fHeI(i) = fHeIw(nw)
         end if
      end do

      open(newunit=u, file=trim(outdir)//'/Hydro_ioniz_IC.txt',         &
           status='replace', action='write')
      write(u,'(A)') '# Wind-AE -> EXHALE IC (converted)'
      write(u,'(A)') '# columns r[Rp] n[cm-3] v[cm/s] p[cgs] T[K] '//    &
                     'heat[erg/cm3/s] cool[erg/cm3/s]'
      block
        integer :: jj
        real*8 :: fh, fhe, nn
        do jj = 1, ng
           fh  = min(max(fHI(jj),  0.0d0), 1.0d0)
           fhe = min(max(fHeI(jj), 0.0d0), 1.0d0)
           nn  = 10.0d0**lognnuc(jj)
           nH  = nn/(1.0d0+HeH); nHe = nH*HeH
           nHI = fh*nH; nHII = (1.0d0-fh)*nH
           nHeI = fhe*nHe; nHeII = (1.0d0-fhe)*nHe
           ne = nHII + nHeII
           ntot = nH + nHe
           p = (ntot + ne)*ICW_kB*T(jj)
           nmass = nHI + nHII + 4.0d0*(nHeI + nHeII)
           write(u,'(1X,7(ES16.10,1X))') r(jj), nmass, v(jj), p, T(jj), 0.0d0, 0.0d0
        end do
      end block
      close(u)

      open(newunit=u, file=trim(outdir)//'/Ion_species_IC.txt',         &
           status='replace', action='write')
      write(u,'(A)') '# Wind-AE -> EXHALE IC (converted)'
      write(u,'(A)') '# columns r[Rp] HI HII HeI HeII HeIII HeITR'
      block
        integer :: jj
        real*8 :: fh, fhe, nn
        do jj = 1, ng
           fh  = min(max(fHI(jj),  0.0d0), 1.0d0)
           fhe = min(max(fHeI(jj), 0.0d0), 1.0d0)
           nn  = 10.0d0**lognnuc(jj)
           nH  = nn/(1.0d0+HeH); nHe = nH*HeH
           write(u,'(1X,7(ES16.10,1X))') r(jj), fh*nH, (1.0d0-fh)*nH,     &
                 fhe*nHe, (1.0d0-fhe)*nHe, 0.0d0, 0.0d0
        end do
      end block
      close(u)
      end subroutine wae_write_ic

      end module wae_ic_writer
