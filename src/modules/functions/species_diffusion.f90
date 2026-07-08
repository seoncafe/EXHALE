      module species_diffusion
      ! He/H diffusive separation for the single-fluid EXHALE wind.
      !
      ! EXHALE's hydro evolves only (rho, momentum, energy); the composition
      ! f_sp is re-solved each step by LOCAL ionization equilibrium, which
      ! conserves the element ratio -> He/H is otherwise frozen at the input
      ! HeH.  This module adds a genuine transport equation for the He element
      ! ratio  fHe = n_He/n_H, advected at the bulk velocity v plus a
      ! molecular-diffusion drift of He relative to H:
      !
      !   d fHe/dt + v d fHe/dr
      !        = (1/n_H) (1/r^2) d/dr[ r^2 n_H D ( d fHe/dr + fHe * G ) ]
      !
      !   G = (m_He - m_H) g / (k T)          [gravitational settling, 1/cm]
      !   D = 1.52e18 (1/m_H+1/m_He)^1/2 T^1/2 / n_tot   [cm^2/s]  (Banks &
      !                                        Kockarts 1973 binary He-in-H)
      !
      ! Design notes: docs/design_hehe_diffusion.md.  Current scope:
      !   * He element only; metals stay frozen to H (scaled with H).
      !   * neutral-atom mass difference Dm = 3 amu (ambipolar deferred; for
      !     He++ in H+ the ambipolar-corrected settling force ~ 3 m_H g, so the
      !     neutral value is a good first approximation).
      !   * thermal diffusion alpha_T = 0.
      ! Gated on he_diffusion (default .false.): when off this module is never
      ! entered, so metals-off / flag-off runs are byte-identical to before.
      !
      ! The operator works entirely on the mass-normalized fractions f_sp
      ! (sum_s m_s f_sp_s = 1), so no dimensional density bookkeeping is needed
      ! for the composition rescale; only the drift RHS uses cgs units.

      use global_parameters
      use grav_func,     only: Dphi
      use species_table, only: isp_HI, isp_HII, isp_HeI, isp_HeII,        &
                               isp_HeIII, isp_HeTR, n_mion, mion_fsp,      &
                               n_melem, melem_i0, melem_top, melem_A

      implicit none
      private
      public :: he_diffusion_step

      ! Neutral He-H mass difference and binary-diffusion prefactor.
      real*8, parameter :: m_H_amu  = 1.0d0
      real*8, parameter :: m_He_amu = 4.0d0
      real*8, parameter :: dm_amu   = m_He_amu - m_H_amu          ! = 3
      real*8, parameter :: m_amu_g  = 1.6726d-24                  ! amu in g
      ! D = Dpref * T^1/2 / n_tot,  Dpref = 1.52e18*(1/mH+1/mHe)^1/2
      real*8, parameter :: Dpref    = 1.52d18*1.118033989d0       ! = 1.699e18


      contains

      ! ------------------------------------------------------------------ !

      subroutine he_diffusion_step(rho, v, Tcode, f_sp, dt_code)
      ! Advance the He/H element ratio one relaxation step (sub-cycled for
      ! diffusion stability) and write the new element split back into f_sp.
      ! rho, v, Tcode are the current adimensional primitives; dt_code is the
      ! adimensional relaxation timestep.  f_sp is modified in place.

      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: rho, v, Tcode
      real*8, dimension(1-Ng:N+Ng,n_species), intent(inout) :: f_sp
      real*8, dimension(1-Ng:N+Ng),           intent(in)    :: dt_code

      real*8, dimension(1-Ng:N+Ng) :: sumH, sumHe, cmet, fHe
      real*8, dimension(1-Ng:N+Ng) :: nH_phys, ntot_phys, TK, Dco, Gco
      real*8, dimension(1-Ng:N+Ng) :: zbHe, zbH, dmeff  ! mean charges; eff. Dm
      real*8, dimension(1-Ng:N+Ng) :: rp, rep, dt_phys
      real*8, dimension(1-Ng:N+Ng) :: nHe, nHl        ! He number dens.; n_H
      real*8, dimension(1-Ng:N+Ng) :: nX, nXold, DcoX, GcoX, zbX ! metal work
      ! NB: local time scale is named tscale (NOT t0) -- a local `t0` would
      ! collide (Fortran is case-insensitive) with the GLOBAL temperature
      ! normalization T0 used in TK = Tcode*T0, silently zeroing TK.
      real*8 :: tscale, MHnew, rH, rHe, mX, DprefX, fXbase, rX
      integer :: j, im, i0m, top, k

      if (.not. he_diffusion) return
      if (.not. thereis_He)   return

      ! --- current element sums and ratio from f_sp (mass sum_s m_s f_sp = 1)
      sumH  = f_sp(:,isp_HI)  + f_sp(:,isp_HII)
      sumHe = f_sp(:,isp_HeI) + f_sp(:,isp_HeII) + f_sp(:,isp_HeIII)
      if (thereis_HeITR) sumHe = sumHe + f_sp(:,isp_HeTR)
      ! metal mass per unit rho = 1 - (H mass) - (He mass); m_H=1, m_He=4 amu
      cmet = 1.0d0 - sumH - m_He_amu*sumHe          ! metal mass fraction
      where (sumH .lt. 1.0d-30) sumH = 1.0d-30
      fHe = sumHe/sumH                              ! n_He/n_H (m_H=1)

      ! --- dimensional fields for the drift RHS (cgs)
      TK        = Tcode*T0
      nH_phys   = sumH*rho*n0                        ! n_H  [cm^-3]
      ntot_phys = (sumH + sumHe)*rho*n0              ! background N [cm^-3]
      where (ntot_phys .lt. 1.0d0) ntot_phys = 1.0d0
      where (TK        .lt. 1.0d0) TK        = 1.0d0
      Dco = Dpref*sqrt(TK)/ntot_phys                 ! D [cm^2/s]

      ! P2b -- ambipolar-corrected effective settling mass difference.  In the
      ! ionized wind the polarization field eE = -(1/n_e) dp_e/dr (Koskinen+2013)
      ! lifts ions; for an H+ background it supports ~half a proton weight per
      ! unit charge, so a species of charge Z has effective mass m - Z*m_H/2.
      ! The He-vs-H relative settling mass is then
      !   Dm_eff = (m_He - Zbar_He*m_H/2) - (m_H - Zbar_H*m_H/2)
      !          = 3 - 0.5*(Zbar_He - Zbar_H)   [amu],
      ! -> 3 at the neutral base, 2.5 in the fully ionized (He++/H+) wind.
      ! Zbar are the mean charges from the local ionization state (f_sp).
      ! Set he_ambipolar=.false. to recover the constant neutral value (3).
      if (he_ambipolar) then
         do j = 1-Ng, N+Ng
            zbHe(j) = (f_sp(j,isp_HeII) + 2.0d0*f_sp(j,isp_HeIII))        &
                      / max(sumHe(j),1.0d-30)
            zbH(j)  =  f_sp(j,isp_HII) / max(sumH(j),1.0d-30)
         enddo
         dmeff = dm_amu - 0.5d0*(zbHe - zbH)
      else
         dmeff = dm_amu
      endif
      do j = 1-Ng, N+Ng
         ! settling coefficient G = Dm_eff m_amu g /(k T)   [1/cm]
         Gco(j) = dmeff(j)*m_amu_g*(Dphi(r(j))*v0*v0/R0)/(kb_erg*TK(j))
      enddo
      ! P2c -- thermal diffusion: add alpha_T d(lnT)/dr to the settling
      ! coefficient (drift ~ D[df/dr + f(Dm g/kT + alpha_T dlnT/dr)]).  Default
      ! he_alphaT=0 (no-op).  Central difference on ln T; endpoints one-sided.
      if (he_alphaT .ne. 0.0d0) then
         do j = 2-Ng, N+Ng-1
            Gco(j) = Gco(j) + he_alphaT*(log(TK(j+1))-log(TK(j-1)))       &
                     / max((r(j+1)-r(j-1))*R0, 1.0d0)
         enddo
      endif

      rp  = r*R0
      rep = r_edg*R0

      ! --- CONSERVATIVE implicit (backward-Euler, tridiagonal) update of the
      ! He number density n_He.  Both advection (bulk v) and the diffusive
      ! flux enter ONE conservative divergence of the total face flux
      !   J_{j+1/2} = n_He v  -  D n_H ( d f/dr + f G ),   f = n_He/n_H,
      ! so the scheme is flux-consistent (the earlier non-conservative
      ! material-advection form was not, which drained He).  n_H is lagged
      ! from the current f_sp; the settling scale height is tiny at the cold
      ! base so the solve is implicit (unconditionally stable).  Cell-by-cell
      ! dt_phys (local-timestepping): the fixed point dn_He/dt=0 is
      ! dt-independent.
      tscale  = R0/v0
      dt_phys = dt_code*tscale

      nHl = nH_phys
      where (nHl .lt. 1.0d-30) nHl = 1.0d-30
      nHe = sumHe*rho*n0                                ! current He density
      where (nHe .lt. 0.0d0)   nHe = 0.0d0

      ! Conservative implicit advection-diffusion-settling solve for n_He,
      ! relative to the background n_H, via the shared element kernel.
      call solve_1elem(nHe, nHl, Dco, Gco, HeH, dt_phys, rp, rep, v)

      ! new He/H ratio; cap at HeH.  A near-base pile-up (settling concentrates
      ! He against the fixed reservoir base; the lagged n_H cannot feed back to
      ! limit f=n_He/n_H) is held by this physical bound -- He/H cannot exceed
      ! the reservoir value in a settling+escape column.  A self-consistent n_H
      ! Picard iteration was tried and does NOT help: the concentration is real,
      ! so n_H is driven to 0 and f still diverges.  The aloft separation (the
      ! observable) sits below HeH and is unaffected by the cap.
      fHe = nHe/nHl
      where (fHe .gt. HeH) fHe = HeH
      fHe(1-Ng:1)   = HeH                               ! base + inner ghosts
      fHe(N+1:N+Ng) = fHe(N)                            ! outer ghosts

      ! --- write the new element split back into f_sp, conserving mass
      !     (sum_s m_s f_sp = 1) with metals frozen to H (metal/H ratio fixed,
      !     so metals scale by the same rH as H).  With MH = sumH (m_H=1),
      !     He mass = m_He*fHe*MHnew and metal mass = (cmet/sumH)*MHnew, mass
      !     balance  MHnew (1 + m_He*fHe + cmet/sumH) = 1  gives MHnew below.
      !     When fHe is unchanged this returns MHnew=sumH (rH=rHe=1), i.e.
      !     f_sp is untouched -> byte-identical no-op.
      do j = 1-Ng, N+Ng
         MHnew = 1.0d0/(1.0d0 + m_He_amu*fHe(j) + cmet(j)/sumH(j))
         if (MHnew .lt. 1.0d-30) MHnew = 1.0d-30
         rH  = MHnew/sumH(j)
         if (sumHe(j) .gt. 1.0d-30) then
            rHe = fHe(j)*MHnew/sumHe(j)
         else
            rHe = 1.0d0
         endif
         f_sp(j,isp_HI)    = f_sp(j,isp_HI)   *rH
         f_sp(j,isp_HII)   = f_sp(j,isp_HII)  *rH
         f_sp(j,isp_HeI)   = f_sp(j,isp_HeI)  *rHe
         f_sp(j,isp_HeII)  = f_sp(j,isp_HeII) *rHe
         f_sp(j,isp_HeIII) = f_sp(j,isp_HeIII)*rHe
         if (thereis_HeITR) f_sp(j,isp_HeTR) = f_sp(j,isp_HeTR)*rHe
         do im = 1, n_mion                            ! metals moved with H
            f_sp(j,mion_fsp(im)) = f_sp(j,mion_fsp(im))*rH
         enddo
      enddo

      ! --- P2d: diffuse each trace metal element relative to H.  Metals are
      ! trace (~1e-3 of the mass), so their diffusion does not feed back on
      ! n_H; each element is diffused independently against the (post-rescale)
      ! background n_H with its own mass melem_A and mean charge (ambipolar).
      ! The He rescale above left metals frozen to H (metal/H = base value);
      ! this loop adds their diffusive separation on top.  Default OFF.
      if (he_metal_diffusion .and. thereis_metals) then
         nHl = (f_sp(:,isp_HI)+f_sp(:,isp_HII))*rho*n0     ! updated n_H
         where (nHl .lt. 1.0d-30) nHl = 1.0d-30
         do im = 1, n_melem
            i0m = melem_i0(im)
            top = melem_top(im)
            mX  = melem_A(im)
            ! element total density and mean charge (stage k has charge k)
            nX  = 0.0d0
            zbX = 0.0d0
            do k = 0, top
               nX  = nX  + f_sp(:,mion_fsp(i0m+k))*rho*n0
               zbX = zbX + dble(k)*f_sp(:,mion_fsp(i0m+k))*rho*n0
            enddo
            where (nX .gt. 1.0d-30)
               zbX = zbX/nX
            elsewhere
               zbX = 0.0d0
            end where
            nXold = nX
            ! element diffusion coefficient and settling (ambipolar + thermal)
            DprefX = 1.52d18*sqrt(1.0d0 + 1.0d0/mX)   ! Banks&Kockarts, X-in-H
            DcoX = DprefX*sqrt(TK)/ntot_phys
            if (he_ambipolar) then
               dmeff = (mX - m_H_amu) - 0.5d0*(zbX - zbH)
            else
               dmeff = mX - m_H_amu
            endif
            do j = 1-Ng, N+Ng
               GcoX(j) = dmeff(j)*m_amu_g*(Dphi(r(j))*v0*v0/R0)           &
                         /(kb_erg*TK(j))
            enddo
            if (he_alphaT .ne. 0.0d0) then
               do j = 2-Ng, N+Ng-1
                  GcoX(j) = GcoX(j) + he_alphaT*(log(TK(j+1))-log(TK(j-1)))&
                            / max((r(j+1)-r(j-1))*R0, 1.0d0)
               enddo
            endif
            fXbase = nXold(1)/nHl(1)                       ! reservoir metal/H
            call solve_1elem(nX, nHl, DcoX, GcoX, fXbase, dt_phys, rp, rep, v)
            ! Rescale this element's stages to the diffused total.  The cap
            ! min(nX, fXbase*nHl) enforces metal/H <= reservoir (same physics
            ! as the He cap); rX itself is NOT clamped to <=1 -- diffusion may
            ! legitimately replenish a previously depleted cell, and clamping
            ! would make depletion a one-way ratchet.  When the old element
            ! density is negligible (stage ratios meaningless to scale), seed
            ! the solved amount into the neutral stage instead; the ionization
            ! equilibrium re-partitions the stages on the next solve.
            do j = 1-Ng, N+Ng
               rX = min(nX(j), fXbase*nHl(j))              ! target density
               if (rX .lt. 0.0d0) rX = 0.0d0
               if (nXold(j) .gt. 1.0d-25*nHl(j)) then
                  rX = rX/nXold(j)                         ! scale factor
                  do k = 0, top
                     f_sp(j,mion_fsp(i0m+k)) = f_sp(j,mion_fsp(i0m+k))*rX
                  enddo
               else if (rX .gt. 1.0d-25*nHl(j)) then
                  ! element returned to an exhausted cell: re-seed (neutral)
                  f_sp(j,mion_fsp(i0m)) = rX/max(rho(j)*n0, 1.0d-30)
                  do k = 1, top
                     f_sp(j,mion_fsp(i0m+k)) = 0.0d0
                  enddo
               endif
            enddo
         enddo
      endif

      end subroutine he_diffusion_step

      ! ------------------------------------------------------------------ !

      ! Shared kernel: one conservative implicit (backward-Euler, tridiagonal)
      ! advection-diffusion-settling solve for an element number density nX
      ! diffusing relative to a background n_H (nHl), with element diffusion
      ! coefficient Dco, settling coefficient Gco, and fixed reservoir base
      ! ratio fbase (nX/nHl at the base).  Used for He and each trace metal.
      ! nX is intent(inout): supply the current density, receive the solved one.
      subroutine solve_1elem(nX, nHl, Dco, Gco, fbase, dt_phys, rp, rep, v)
      real*8, dimension(1-Ng:N+Ng), intent(inout) :: nX
      real*8, dimension(1-Ng:N+Ng), intent(in)    :: nHl, Dco, Gco, dt_phys
      real*8, dimension(1-Ng:N+Ng), intent(in)    :: rp, rep, v
      real*8,                       intent(in)    :: fbase

      real*8, dimension(0:N)       :: Lf, Rf
      real*8, dimension(1-Ng:N+Ng) :: aa, bb, cc, dd, cp, dp
      real*8 :: dr_f, nHf, Df, Gf, DK, Kj, vf, mden, PL, PR, AL, AR, nX_base
      integer :: j

      ! face flux coefficients  J(j) = Lf(j) nX(j) + Rf(j) nX(j+1) at r_{j+1/2}
      Lf = 0.0d0
      Rf = 0.0d0
      do j = 1, N-1
         dr_f = max(rp(j+1)-rp(j), 1.0d0)
         nHf  = 0.5d0*(nHl(j)+nHl(j+1))
         Df   = 0.5d0*(Dco(j)+Dco(j+1))
         Gf   = 0.5d0*(Gco(j)+Gco(j+1))
         DK   = Df + he_kzz
         ! diffusion (central) + settling flux -nHf*Df*Gf*f.  Peclet-based
         ! hybrid on the settling drift: use central differencing where it is
         ! well-resolved and stable (settling Peclet |Df*Gf|*dr < 2*DK, e.g.
         ! light He -- more accurate/less numerical diffusion), and first-order
         ! upwind where not (heavy metals, large Gf -- keeps the tridiagonal
         ! diagonally dominant / NaN-free).  Gf>0 settles inward -> upwind from
         ! the cell above (j+1); Gf<0 -> from below (j).
         if (abs(Df*Gf)*dr_f .lt. 2.0d0*DK) then
            PL =  nHf/nHl(j)  *(DK/dr_f - 0.5d0*Df*Gf)
            PR = -nHf/nHl(j+1)*(DK/dr_f + 0.5d0*Df*Gf)
         else if (Gf .ge. 0.0d0) then
            PL =  nHf/nHl(j)  *(DK/dr_f)
            PR = -nHf/nHl(j+1)*(DK/dr_f + Df*Gf)
         else
            PL =  nHf/nHl(j)  *(DK/dr_f - Df*Gf)
            PR = -nHf/nHl(j+1)*(DK/dr_f)
         endif
         vf   = 0.5d0*(v(j)+v(j+1))*v0
         if (vf .ge. 0.0d0) then
            AL = vf ; AR = 0.0d0
         else
            AL = 0.0d0 ; AR = vf
         endif
         Lf(j) = AL + PL
         Rf(j) = AR + PR
      enddo
      vf    = max(v(N)*v0, 0.0d0)
      Lf(N) = vf
      Rf(N) = 0.0d0

      do j = 2, N
         Kj    = 1.0d0/(rp(j)**2*max(rep(j)-rep(j-1),1.0d0))
         aa(j) = -Kj*rep(j-1)**2*Lf(j-1)
         bb(j) =  1.0d0/dt_phys(j)                                        &
                + Kj*(rep(j)**2*Lf(j) - rep(j-1)**2*Rf(j-1))
         cc(j) =  Kj*rep(j)**2*Rf(j)
         dd(j) =  nX(j)/dt_phys(j)
      enddo
      cc(N) = 0.0d0

      nX_base = fbase*nHl(1)                              ! base Dirichlet
      nX(1-Ng:1) = nX_base
      dd(2) = dd(2) - aa(2)*nX_base
      aa(2) = 0.0d0

      cp(2) = cc(2)/bb(2)
      dp(2) = dd(2)/bb(2)
      do j = 3, N
         mden  = bb(j) - aa(j)*cp(j-1)
         cp(j) = cc(j)/mden
         dp(j) = (dd(j) - aa(j)*dp(j-1))/mden
      enddo
      nX(N) = dp(N)
      do j = N-1, 2, -1
         nX(j) = dp(j) - cp(j)*nX(j+1)
      enddo
      do j = 2, N
         if (nX(j) .lt. 0.0d0) nX(j) = 0.0d0
      enddo

      end subroutine solve_1elem

      ! End of module
      end module species_diffusion
