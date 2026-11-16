      module wae_ic_writer
      ! Write {Hydro_ioniz,Ion_species}_IC.txt in the supplied outdir.
      ! Fortran port of src/utils/windae_to_exhale_ic.py.
      ! np_interp is linear and endpoint-clamped, matching numpy.interp.
      ! Below Rmin use the hydrostatic base blend; above Rmax conserve flux.
      !
      ! Differences from src/utils/windae_to_exhale_ic.py (2026-10-27,
      ! md/windae_defect_fixes_20261027.md):
      ! - Scale-height point (D14): the script takes the first wind point at
      !   least 0.1 Rp above Rmin and divides 0/0 when none exists.  Here
      !   the farthest wind point is used then, with a warning; a wind
      !   solution with no point above Rmin stops the writer.
      ! - Base blend (D15): the weight (1 - 1/r)/(1 - 1/Rmin) is 0/0 when
      !   the wind base is at Rmin = 1 Rp.  That case has no blend region:
      !   points below Rmin take the hydrostatic baseline (weight 0), the
      !   point at Rmin the wind value (weight 1), the limits of the
      !   formula at r <= 1 and r = Rmin.
      implicit none
      ! H atom mass, the density unit the IC is written in (matches mu in
      ! parameters.f90 and wae_MH).
      real*8, parameter :: ICW_mH = 1.67353284d-24
      ! Helium atom mass [g] and the helium/hydrogen mass ratio the density
      ! conversions below weigh helium with.  These are the SAME two masses
      ! wae_exhale_input.f90 hands the Wind-AE solver as atomic_mass(1) and
      ! atomic_mass(2), so the windsoln that comes back and the IC written
      ! from it count the same gas; the Wind-AE tree is compiled into the
      ! standalone wind_ae_ic.x and therefore cannot use global_parameters or
      ! species_table, where the same ratio is the definition
      ! (bsp_mass(He I), which global_parameters forms from the same two
      ! atom masses, so the two ratios are the same number, 3.9715259).
      real*8, parameter :: ICW_mHe = 6.6464790722d-24
      real*8, parameter :: ICW_mHe_over_mH = ICW_mHe/ICW_mH
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

      real*8 function phiR(r, Mrapp, atilde, tidalf)
      ! Roche potential / b0, matching grav_field.f90::phi and wae_soe d1_phi.
      ! tidalf=0 gives planetary -1/r; tidalf=1 includes stellar and
      ! centrifugal terms. The base density uses this potential before
      ! blending to the wind; the blend need not be hydrostatic.
      real*8, intent(in) :: r, Mrapp, atilde, tidalf
      phiR = -1.0d0/r                                                     &
             - tidalf*( Mrapp/(atilde - r)                               &
               + (1.0d0 + Mrapp)/(2.0d0*atilde**3)                       &
                 *(atilde*Mrapp/(1.0d0 + Mrapp) - r)**2 )
      end function phiR

      subroutine wae_write_ic(rw, rhow, vw, Tw, YsHIw, YsHeIw, nw,        &
                              rgrid, ng, lognbase, T0, Rp, HeH,           &
                              Mrapp, atilde, tidalf, outdir)
      ! rw[cm]/rhow[g/cc]/vw[cm/s]/Tw[K]/Ys* : windsoln (length nw, r asc)
      ! rgrid : EXHALE grid radius [Rp] (length ng, incl. ghosts)
      ! Mrapp=Mstar/Mp, atilde=a/Rp, tidalf=tidalforce: set the Roche base blend
      integer, intent(in) :: nw, ng
      real*8,  intent(in) :: rw(nw), rhow(nw), vw(nw), Tw(nw)
      real*8,  intent(in) :: YsHIw(nw), YsHeIw(nw), rgrid(ng)
      real*8,  intent(in) :: lognbase, T0, Rp, HeH, Mrapp, atilde, tidalf
      character(len=*), intent(in) :: outdir
      real*8  :: rwR(nw), nnuc_w(nw), lognnuc_w(nw), fHIw(nw), fHeIw(nw)
      real*8  :: lognnuc(ng), T(ng), v(ng), fHI(ng), fHeI(ng), r(ng)
      real*8  :: nbase, rmin_w, rmax_w, b0_eff, ln10, sblend, logw_at, hyd_at
      real*8, parameter :: dr_scale_height = 0.1d0   ! [Rp]
      real*8  :: nn_rmin, v_rmin, T_rmin_w, fHI_rmin, fHeI_rmin
      real*8  :: nH, nHe, nHI, nHII, nHeI, nHeII, ne, ntot, p, nmass
      integer :: i, i1, u
      ln10 = log(10.0d0)
      nbase = 10.0d0**lognbase
      r = rgrid
      do i = 1, nw
         rwR(i) = rw(i)/Rp
         nnuc_w(i) = (rhow(i)/(ICW_mH*(1.0d0+ICW_mHe_over_mH*HeH)))     &
                     *(1.0d0+HeH)
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

      ! b0_eff from the windsoln scale height just above Rmin: the first
      ! wind point at least dr_scale_height above Rmin, else the farthest
      ! wind point (D14).
      i1 = nw
      do i = 1, nw
         if (rwR(i) .ge. rmin_w + dr_scale_height) then
            i1 = i; exit
         end if
      end do
      if (.not. (rwR(i1) .gt. rwR(1))) then
         write(*,'(A)') ' (wae_write_ic) the wind solution has no point'//  &
                        ' above its base; no scale height can be formed.'
         error stop '(wae_write_ic) degenerate wind solution'
      end if
      if (rwR(i1) .lt. rmin_w + dr_scale_height) then
         write(*,'(A,F8.5,A)') ' (wae_write_ic) WARNING: no wind point 0.1 Rp'// &
            ' above Rmin; scale height taken from the farthest point,',      &
            rwR(i1) - rmin_w, ' Rp above Rmin.'
      end if
      b0_eff = (log(nnuc_w(1)) - log(nnuc_w(i1))) /                      &
               (1.0d0/rwR(1) - 1.0d0/rwR(i1))
      logw_at = np_interp(rmin_w, rwR, lognnuc_w, nw)
      ! Roche-consistent hydrostatic baseline: ln n = ln n_base
      !   - (phiR(r) - phiR(1)) * b0_eff,  reducing to the legacy spherical
      ! b0_eff*(1/r - 1) when tidalf = 0.
      hyd_at  = lognbase - (phiR(rmin_w, Mrapp, atilde, tidalf)            &
                - phiR(1.0d0, Mrapp, atilde, tidalf))*b0_eff/ln10
      nn_rmin = 10.0d0**np_interp(rmin_w, rwR, lognnuc_w, nw)
      v_rmin  = np_interp(rmin_w, rwR, vw, nw)
      T_rmin_w  = np_interp(rmin_w, rwR, Tw, nw)
      fHI_rmin  = np_interp(rmin_w, rwR, fHIw, nw)
      fHeI_rmin = np_interp(rmin_w, rwR, fHeIw, nw)

      do i = 1, ng
         if (r(i) .le. rmin_w) then
            ! Below Rmin, blend the hydrostatic baseline (weight 0 at
            ! r = 1) to the wind (weight 1 at Rmin).  With Rmin <= 1 there
            ! is no blend region: baseline below Rmin, wind at Rmin (D15).
            if (rmin_w .gt. 1.0d0) then
               sblend = (1.0d0 - 1.0d0/r(i)) / (1.0d0 - 1.0d0/rmin_w)
               sblend = min(max(sblend, 0.0d0), 1.0d0)
            else if (r(i) .ge. rmin_w) then
               sblend = 1.0d0
            else
               sblend = 0.0d0
            end if
            lognnuc(i) = (lognbase - (phiR(r(i), Mrapp, atilde, tidalf)  &
                         - phiR(1.0d0, Mrapp, atilde, tidalf))*b0_eff/ln10) &
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
           nmass = nHI + nHII + ICW_mHe_over_mH*(nHeI + nHeII)
           write(u,'(1X,7(ES17.10,1X))') r(jj), nmass, v(jj), p, T(jj), 0.0d0, 0.0d0
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
           write(u,'(1X,7(ES17.10,1X))') r(jj), fh*nH, (1.0d0-fh)*nH,     &
                 fhe*nHe, (1.0d0-fhe)*nHe, 0.0d0, 0.0d0
        end do
      end block
      close(u)
      end subroutine wae_write_ic

      end module wae_ic_writer
