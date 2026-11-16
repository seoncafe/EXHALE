      module wae_soe
      ! Wind-AE soe.c: mean molecular weight, recombination, heating/cooling,
      ! and the coupled wind derivatives, including sonic-point linearization.
      ! Residual equations and their dispatcher are in wae_eqns.
      !
      ! Differences from upstream Wind-AE (soe.c), decided 2026-10-27 to
      ! follow the physics rather than the upstream code
      ! (md/windae_defect_fixes_20261027.md):
      ! - Electron density in get_spQ: upstream sums n_e only inside
      !   `if (lyacool)`, so recombination cooling sees n_e = 0 when Ly-alpha
      !   cooling (and with it the line cooling) is off.  Here n_e is summed from the composition
      !   unconditionally (D16).
      ! - Recombination cooling temperature: upstream evaluates the H fit
      !   and the 3/2 k_B T term of the other ions with the scaled
      !   T/T0 (T0 = 1e4 K).  Here both take T in K (D17); the H fit is
      !   Black (1981, MNRAS 197, 553, Table 3), defined for T in K.
      ! - Mean molecular weight: upstream counts one free electron at most
      !   for each nucleus (2 - Ys); here the electrons already removed
      !   from the lowest state of the species (Z - N_e, as in the n_e sum)
      !   are counted too, 2 + (Z - N_e) - Ys.  Identical for H and He.
      ! - dmu/dr (D18, known, on hold by user decision 2026-10-27): upstream
      !   uses -mu^2 * sum_j HX_j dYs_j/dr, which drops the m_H/m_j weights of
      !   get_mu and the dependence of the molecular blend erfc(v) on v(r).
      !   The exact derivative is known but not used: without the blend's
      !   implicit heat source the HD 209458 b tidal ramp stalls.  This port
      !   keeps the upstream form.
      use wae_config,   only: nsp => wae_nspecies, T0 => wae_T0,          &
                              RHO0 => wae_RHO0, gamma_atomic=>wae_gamma_atomic
      use wae_types,    only: wae_i_eqnvars
      use wae_params,   only: par => wae_par, CS0 => wae_cs0_val
      use wae_glq_rates,only: wae_glq_rates_eval
      use wae_status,   only: wae_err
      use wae_rnew,     only: rnew
      use wae_rrec,     only: rrec
      use wae_fe,       only: fe
      use wae_lc_cii,   only: CII
      use wae_lc_ciii,  only: CIII
      use wae_lc_oii,   only: OII
      use wae_lc_oiii,  only: OIII
      implicit none

      ! physical constants. CODATA values, matching wae_params (see the note
      ! there); defs.h originally carried K = 1.380658d-16.
      real*8, parameter :: wae_K   = 1.380649d-16
      real*8, parameter :: wae_PI  = 3.141592653589793d0
      real*8, parameter :: wae_SIGSB = 5.6705d-5
      real*8, parameter :: LYACOOL_COEFF = -7.5d-19
      ! Ly-alpha collisional excitation cooling, Black (1981, MNRAS 197,
      ! 553, Table 3): 7.5e-19 exp(-118348/T) n_e n(H) erg cm^-3 s^-1;
      ! the threshold is stored in units of T0 because v%T is T/T0.
      real*8, parameter :: LYACOOL_TEMP  = 118348.0d0/T0

      ! soe.c globals
      real*8 :: wae_erf_norm = 1.0d0
      real*8 :: wae_dvdr_slope = 0.0d0
      real*8 :: wae_dvdr_last  = 0.0d0
      real*8 :: wae_q_last     = 0.0d0
      integer, parameter :: QCOMP = 100
      contains

      subroutine wae_get_gamma(g)
      real*8, intent(out) :: g
      g = gamma_atomic
      end subroutine wae_get_gamma

      subroutine wae_get_rad(rad, v)
      real*8, intent(out) :: rad
      type(wae_i_eqnvars), intent(in) :: v
      rad = par%Rmin + v%q*v%z
      end subroutine wae_get_rad

      real*8 function wae_atomic_particle_number(kv)
      ! Free particles (nuclei plus electrons) in the atomic gas, in units
      ! of 1/m_H per gram: sum_j (m_H/m_j) HX_j (2 + (Z_j - N_e,j) - Ys_j).
      ! Each species carries Z - N_e electrons already removed in its lowest
      ! state and one more in its ionized fraction, as in the n_e sums.
      type(wae_i_eqnvars), intent(in) :: kv
      integer :: j
      real*8  :: g
      g = 0.0d0
      do j = 1, nsp
         g = g + (par%atomic_mass(1)/par%atomic_mass(j))*par%HX(j)*       &
                 (2.0d0 + dble(par%Z(j) - par%N_e(j)) - kv%Ys(j))
      end do
      wae_atomic_particle_number = g
      end function wae_atomic_particle_number

      real*8 function wae_molecular_fraction(kv)
      ! Weight of the molecular layer, bolo_heat_cool*erfc((v-v_drop)/w)
      ! normalized by wae_erf_norm; it blends mu and the bolometric terms.
      type(wae_i_eqnvars), intent(in) :: kv
      wae_molecular_fraction = par%bolo_heat_cool*                        &
         (1.0d0 - erf((kv%v - par%erf_drop(1))/par%erf_drop(2)))/wae_erf_norm
      end function wae_molecular_fraction

      subroutine wae_get_mu(mu, kv)
      ! Mean molecular weight in units of m_H (species 1 must be HI):
      ! mu = (1 - s)/g + s*molec_adjust, with g the atomic particle number
      ! and s the molecular-layer weight.
      real*8, intent(out) :: mu
      type(wae_i_eqnvars), intent(in) :: kv
      real*8  :: serf
      serf = wae_molecular_fraction(kv)
      mu = (1.0d0 - serf)/wae_atomic_particle_number(kv)                  &
           + serf*par%molec_adjust
      end subroutine wae_get_mu

      subroutine wae_get_alpharec(alpharec, kv, calc_lower_ion_state)
      ! CLOUDY rrfit recombination coefficients (rnew/rrec/fe tables).
      real*8, intent(out) :: alpharec(nsp)
      type(wae_i_eqnvars), intent(in) :: kv
      integer, intent(in) :: calc_lower_ion_state
      integer :: j, i, iz, in
      real*8  :: tt, T, r
      T = kv%T*T0
      do j = 1, nsp
         iz = par%Z(j)
         in = par%N_e(j)
         if (calc_lower_ion_state .eq. 0) then
            if (in .lt. iz) then
               in = in + 1
            else
               alpharec(j) = 0.0d0
               cycle
            end if
         end if
         r = 0.0d0
         if (in .lt. 3 .or. in .eq. 11 .or. (iz .gt. 5 .and. iz .lt. 9)   &
             .or. iz .eq. 10) then
            do i = 1, 126
               if (int(rnew(i,1)) .eq. iz .and. int(rnew(i,2)) .eq. in) then
                  tt = sqrt(T/rnew(i,5))
                  r = rnew(i,3) / ( tt*(tt+1.0d0)**(1.0d0-rnew(i,4)) *     &
                      (1.0d0+sqrt(T/rnew(i,6)))**(1.0d0+rnew(i,4)) )
               end if
            end do
         else
            tt = T*1.0d-4
            if (iz .eq. 26 .and. in .lt. 13) then
               do i = 1, 10
                  if (int(fe(i,1)) .eq. in) then
                     r = fe(i,2)/tt**(fe(i,3)+fe(i,4)*log10(tt))
                  end if
               end do
            else
               do i = 1, 349
                  if (int(rrec(i,1)) .eq. iz .and. int(rrec(i,2)) .eq. in) then
                     r = rrec(i,3)/tt**rrec(i,4)
                  end if
               end do
            end if
         end if
         alpharec(j) = r
      end do
      end subroutine wae_get_alpharec

      subroutine wae_get_dNcoldr(dNdr, v)
      real*8, intent(out) :: dNdr(nsp)
      type(wae_i_eqnvars), intent(in) :: v
      integer :: j
      real*8  :: nH
      do j = 1, nsp
         nH = (RHO0/par%atomic_mass(j))*par%HX(j)*v%rho
         dNdr(j) = -v%Ys(j)*nH
         dNdr(j) = dNdr(j)*par%Rp/1.0d17     ! NCOL0
      end do
      end subroutine wae_get_dNcoldr

      subroutine wae_get_drhodr(drhodr, v, dvdr)
      real*8, intent(out) :: drhodr
      type(wae_i_eqnvars), intent(in) :: v
      real*8, intent(in)  :: dvdr
      real*8 :: rad
      call wae_get_rad(rad, v)
      drhodr = -v%rho*(2.0d0/rad + dvdr/v%v)
      end subroutine wae_get_drhodr

      subroutine wae_get_dYsdr(dYsdr, v, kk)
      ! Neutral-fraction gradient: recombination minus primary/secondary
      ! ionization, divided by the flow speed; radius is in Rp units.
      real*8, intent(out) :: dYsdr(nsp)
      type(wae_i_eqnvars), intent(in) :: v
      real*8, intent(in)  :: kk
      integer :: q, qq, eps
      real*8  :: ion_rate, rec_rate, ne, secondary
      real*8  :: alpharec(nsp), n_tot_q, n0_q, n_ion_q, n0_qq
      real*8  :: ionarr(nsp*(nsp+1)), heatarr(nsp)

      ne = 0.0d0
      do q = 1, nsp
         n_tot_q = RHO0*v%rho*par%HX(q)/par%atomic_mass(q)
         n0_q = v%Ys(q)*n_tot_q
         n_ion_q = (1.0d0 - v%Ys(q))*n_tot_q
         eps = par%Z(q) - par%N_e(q)
         ne = ne + n0_q*dble(eps) + n_ion_q*dble(eps+1)
      end do

      call wae_get_alpharec(alpharec, v, 1)
      call wae_glq_rates_eval(v%Ncol, v%Ys, ionarr, heatarr)

      do q = 1, nsp
         n_tot_q = par%HX(q)/par%atomic_mass(q)        ! over rho
         n0_q = v%Ys(q)*n_tot_q
         n_ion_q = (1.0d0 - v%Ys(q))*n_tot_q
         ion_rate = n0_q*ionarr((q-1)*(nsp+1) + 1)     ! primary
         secondary = 0.0d0
         do qq = 1, nsp
            n0_qq = v%Ys(qq)*par%HX(qq)/par%atomic_mass(qq)
            secondary = secondary + n0_qq*ionarr((qq-1)*(nsp+1) + 1 + q)
            ion_rate  = ion_rate  + n0_qq*ionarr((qq-1)*(nsp+1) + 1 + q)
         end do
         rec_rate = alpharec(q)*ne*n_ion_q
         dYsdr(q) = (rec_rate - ion_rate)/(v%v*n_tot_q) * par%Rp/CS0
      end do
      end subroutine wae_get_dYsdr

      subroutine wae_get_spQ(spQ, v, kk)
      ! Net specific heating/cooling in CS0^3/Rp units: photoheating,
      ! bolometric exchange, Ly-alpha, line, and recombination cooling.
      ! v%T is T/T0; every rate fit below is evaluated at T = v%T*T0 in K.
      real*8, intent(out) :: spQ
      type(wae_i_eqnvars), intent(in) :: v
      real*8, intent(in)  :: kk
      integer :: j, sp, line, eps
      real*8  :: n0overrho, nIONoverrho, nIONjoverrho
      real*8  :: ne, n_ion_j, n0_j, n_tot_j
      real*8  :: spQ_photoheat, spQ_lyacool, spQ_reccool, spQ_linecool
      real*8  :: spQ_boloheat, spQ_bolocool, serf
      real*8  :: kappa_opt, kappa_IR, nHIoverrho
      real*8  :: alpha_hi(nsp), alpha_lo(nsp), T_K
      real*8  :: ionarr(nsp*(nsp+1)), heatarr(nsp)

      kappa_opt = 4.0d-3
      kappa_IR  = 1.0d-2
      spQ = 0.0d0
      T_K = v%T*T0

      ! Electron density from the composition: Z - N_e electrons from the
      ! lowest state of each species and one more from its ionized
      ! fraction (the sum of wae_get_dYsdr).  Needed by the recombination
      ! cooling whatever the Ly-alpha switch (D16).
      ne = 0.0d0
      do j = 1, nsp
         n_tot_j = RHO0*v%rho*par%HX(j)/par%atomic_mass(j)
         n0_j = v%Ys(j)*n_tot_j
         n_ion_j = (1.0d0 - v%Ys(j))*n_tot_j
         eps = par%Z(j) - par%N_e(j)
         ne = ne + n0_j*dble(eps) + n_ion_j*dble(eps+1)
      end do

      ! photo-ionization heating
      call wae_glq_rates_eval(v%Ncol, v%Ys, ionarr, heatarr)
      spQ_photoheat = 0.0d0
      do j = 1, nsp
         n0overrho = v%Ys(j)*par%HX(j)/par%atomic_mass(j)
         spQ_photoheat = spQ_photoheat + n0overrho*heatarr(j)
      end do
      spQ_photoheat = spQ_photoheat*(par%Rp/CS0**3)
      spQ = spQ + spQ_photoheat

      if (par%bolo_heat_cool .gt. 0.0d0) then
         serf = 1.0d0 - erf((v%v - par%erf_drop(1))/par%erf_drop(2))
         serf = serf/wae_erf_norm
         spQ_boloheat = par%Lstar/(4.0d0*wae_PI*par%semimajor**2) *       &
                        (kappa_opt*serf + 0.25d0*kappa_IR*serf)
         spQ = spQ + par%bolo_heat_cool*spQ_boloheat*(par%Rp/CS0**3)
         spQ_bolocool = -2.0d0*wae_SIGSB*(v%T*T0)**4*kappa_IR*serf
         spQ = spQ + par%bolo_heat_cool*spQ_bolocool*(par%Rp/CS0**3)
      end if

      if (par%lyacool .eq. 1) then
         nHIoverrho = v%Ys(1)*par%HX(1)/par%atomic_mass(1)
         spQ_lyacool = LYACOOL_COEFF*exp(-LYACOOL_TEMP/v%T)*nHIoverrho*ne
         spQ_lyacool = spQ_lyacool*(par%Rp/CS0**3)
         spQ = spQ + spQ_lyacool

         ! line cooling (only triggers for CI/CII/OI/OII species)
         do j = 1, nsp
            if (trim(par%species(j)) .eq. 'CI') then
               sp = j
               nIONjoverrho = (1.0d0-v%Ys(sp))*par%HX(sp)/par%atomic_mass(sp)
               do line = 1, 3
                  spQ_linecool = nIONjoverrho*ne*CII(line,2)*             &
                     exp(-CII(line,3)/(v%T*T0)) / (ne*(1.0d0+CII(line,4)/ne))
                  spQ_linecool = spQ_linecool*(par%Rp/CS0**3)/1.014d0
                  spQ = spQ + spQ_linecool
               end do
            end if
            if (trim(par%species(j)) .eq. 'CII') then
               sp = j
               nIONjoverrho = (1.0d0-v%Ys(sp))*par%HX(sp)/par%atomic_mass(sp)
               do line = 1, 2
                  spQ_linecool = nIONjoverrho*ne*CIII(line,2)*            &
                     exp(-CIII(line,3)/(v%T*T0)) / (ne*(1.0d0+CIII(line,4)/ne))
                  spQ_linecool = spQ_linecool*(par%Rp/CS0**3)
                  spQ = spQ + spQ_linecool
               end do
            end if
            if (trim(par%species(j)) .eq. 'OI') then
               sp = j
               nIONjoverrho = (1.0d0-v%Ys(sp))*par%HX(sp)/par%atomic_mass(sp)
               do line = 1, 4
                  spQ_linecool = nIONjoverrho*ne*OII(line,2)*             &
                     exp(-OII(line,3)/(v%T*T0)) / (ne*(1.0d0+OII(line,4)/ne))
                  spQ_linecool = spQ_linecool*(par%Rp/CS0**3)
                  spQ = spQ + spQ_linecool
               end do
            end if
            if (trim(par%species(j)) .eq. 'OII') then
               sp = j
               nIONjoverrho = (1.0d0-v%Ys(sp))*par%HX(sp)/par%atomic_mass(sp)
               do line = 1, 4
                  spQ_linecool = nIONjoverrho*ne*OIII(line,2)*            &
                     exp(-OIII(line,3)/(v%T*T0)) / (ne*(1.0d0+OIII(line,4)/ne))
                  spQ_linecool = spQ_linecool*(par%Rp/CS0**3)
                  spQ = spQ + spQ_linecool
               end do
            end if
         end do
      end if

      ! Recombination cooling (D17: T in K).  H: case A fit of Black (1981,
      ! MNRAS 197, 553, Table 3), 2.85e-27 T^1/2 (5.914 - 0.5 ln T
      ! + 0.01184 T^1/3) n_e n(H+) erg cm^-3 s^-1, stated valid for
      ! 5e3 < T < 5e5 K; the bracket stays positive at all T (its minimum,
      ! 0.14 at T = 2e6 K), so it is used as is outside that range.
      ! Other ions: alpha_rec n_e n_ion (3/2) k_B T, the mean thermal energy
      ! carried by a recombining electron (Broome 2026 PhD thesis, sec.
      ! "Recombination Cooling").
      call wae_get_alpharec(alpha_hi, v, 1)
      call wae_get_alpharec(alpha_lo, v, 0)
      spQ_reccool = 0.0d0
      do j = 1, nsp
         nIONoverrho = (1.0d0-v%Ys(j))*par%HX(j)/par%atomic_mass(j)
         n0overrho   = v%Ys(j)*par%HX(j)/par%atomic_mass(j)
         if (par%Z(j) .eq. 1) then
            spQ_reccool = spQ_reccool - 2.85d-27*ne*nIONoverrho*sqrt(T_K) &
               *(5.914d0 - 0.5d0*log(T_K) + 0.01184d0*T_K**(1.0d0/3.0d0))
         else
            spQ_reccool = spQ_reccool - alpha_hi(j)*nIONoverrho*ne*       &
                          (1.5d0*wae_K*T_K)
         end if
      end do
      spQ_reccool = spQ_reccool*(par%Rp/CS0**3)
      spQ = spQ + spQ_reccool
      end subroutine wae_get_spQ

      subroutine wae_get_dvdr(dvdr, v)
      ! velocity derivative (analytic + erf-weighted linearization near
      ! the sonic point); uses the dvdr_slope/last/q_last globals.
      real*8, intent(out) :: dvdr
      type(wae_i_eqnvars), intent(in) :: v
      real*8 :: mu, rad, gamma, Mach, spQ, d1_phi, norm_a, qinv
      real*8 :: weight, r, s, n

      dvdr = 0.0d0
      call wae_get_gamma(gamma)
      call wae_get_mu(mu, v)
      Mach = v%v/sqrt(gamma*v%T/mu)
      if (Mach .ne. Mach) Mach = 0.0d0
      n = par%erfn
      if (Mach .le. 1.0d0) then
         r = par%rapidity
         s = r/par%mach_limit - r + n
      else
         r = 80.0d0
         s = r/0.999d0 - r + n
      end if
      if (Mach .lt. r/(r+s+n) .or. Mach .gt. (r+2.0d0*(s+n))/(r+s+n)) then
         weight = 1.0d0
      else if (Mach .le. 1.0d0) then
         if (Mach .ge. par%mach_limit) then
            weight = 0.0d0
         else
            weight = 0.5d0*erf(r*(1.0d0/Mach-1.0d0)-s) + 0.5d0*erf(n)
         end if
      else
         if (Mach .le. (r+2.0d0*(s-n))/(r+s-n)) then
            weight = 0.0d0
         else
            weight = 0.5d0*erf(r*(1.0d0/(2.0d0-Mach)-1.0d0)-s) + 0.5d0*erf(n)
         end if
      end if
      if (weight .ne. weight .or. weight .lt. 0.0d0 .or. weight .gt. 1.0d0) then
         wae_err = 2; dvdr = 0.0d0; return   ! graceful fail for continuation
      end if
      if (weight .ne. 0.0d0) then
         call wae_get_rad(rad, v)
         call wae_get_spQ(spQ, v, 1.0d4)
         norm_a = par%semimajor/par%Rp
         d1_phi = 1.0d0/(rad*rad)
         qinv = par%Mstar/par%Mp
         d1_phi = d1_phi + ( -qinv/(norm_a-rad)**2                       &
                  + (qinv*(norm_a-rad) - rad)/norm_a**3 )*par%tidalforce
         d1_phi = d1_phi*(par%Rp/par%H0)
         dvdr = 2.0d0*gamma*v%T/(mu*rad) - (gamma-1.0d0)*spQ/v%v - d1_phi
         dvdr = dvdr*v%v/(v%v*v%v - gamma*v%T/mu)
         dvdr = dvdr*weight
      end if
      if (weight .ne. 1.0d0) then
         dvdr = dvdr + (1.0d0-weight)*(wae_dvdr_slope*(v%q-wae_q_last)    &
                + wae_dvdr_last)
      end if
      end subroutine wae_get_dvdr

      subroutine wae_get_dTdr(dTdr, v, drhodr, dYsdr, kk)
      ! dT/dr = (gamma-1)(Q mu/v + (T/rho) drho/dr) + (T/mu) dmu/dr.
      real*8, intent(out) :: dTdr
      type(wae_i_eqnvars), intent(in) :: v
      real*8, intent(in)  :: drhodr, dYsdr(nsp), kk
      real*8 :: mu, gamma, dmudr, spQ, dgdr
      integer :: j
      call wae_get_mu(mu, v)
      ! dmu/dr as upstream: -mu^2 sum_j HX_j dYs_j/dr.  The exact derivative
      ! of wae_get_mu (m_H/m_j weights and the erfc blend term) is known; it
      ! is not used because, without the blend's implicit heat source, the
      ! tidal ramp of HD 209458 b stalls (md/windae_defect_fixes_20261027.md).
      ! To Be Checked: D18 on hold (user decision 2026-10-27), see md/windae_defect_fixes_20261027.md
      dgdr = 0.0d0
      do j = 1, nsp
         dgdr = dgdr - par%HX(j)*dYsdr(j)
      end do
      dmudr = -mu*mu*dgdr
      call wae_get_gamma(gamma)
      call wae_get_spQ(spQ, v, kk)
      dTdr = (gamma-1.0d0)*(spQ/v%v*mu + v%T/v%rho*drhodr) + v%T/mu*dmudr
      end subroutine wae_get_dTdr

      real*8 function wae_get_component(v, compnum)
      type(wae_i_eqnvars), intent(in) :: v
      integer, intent(in) :: compnum
      real*8 :: comp
      if (compnum .eq. QCOMP) then
         call wae_get_spQ(comp, v, 1.0d4)
         comp = comp*v%rho
      else
         write(*,*) '(wae_get_component) unknown component'; stop 704
      end if
      wae_get_component = comp
      end function wae_get_component

      subroutine wae_get_component_derivs(derivs, v, compnum)
      use wae_types, only: wae_varlist
      type(wae_varlist), intent(out) :: derivs
      type(wae_i_eqnvars), intent(in) :: v
      integer, intent(in) :: compnum
      real*8 :: dvar, cpl, cmi
      type(wae_i_eqnvars) :: vm
      integer :: j
      real*8, parameter :: DERIVDIV = 1.0d7

      dvar = v%rho/DERIVDIV; vm = v
      vm%rho = v%rho + dvar/2.0d0; cpl = wae_get_component(vm, compnum)
      vm%rho = v%rho - dvar/2.0d0; cmi = wae_get_component(vm, compnum)
      derivs%rho = (cpl-cmi)/dvar

      dvar = v%v/DERIVDIV; vm = v
      vm%v = v%v + dvar/2.0d0; cpl = wae_get_component(vm, compnum)
      vm%v = v%v - dvar/2.0d0; cmi = wae_get_component(vm, compnum)
      derivs%v = (cpl-cmi)/dvar

      dvar = v%T/DERIVDIV; vm = v
      vm%T = v%T + dvar/2.0d0; cpl = wae_get_component(vm, compnum)
      vm%T = v%T - dvar/2.0d0; cmi = wae_get_component(vm, compnum)
      derivs%T = (cpl-cmi)/dvar

      do j = 1, nsp
         dvar = v%Ys(j)/DERIVDIV; vm = v
         vm%Ys(j) = v%Ys(j) + dvar/2.0d0; cpl = wae_get_component(vm, compnum)
         vm%Ys(j) = v%Ys(j) - dvar/2.0d0; cmi = wae_get_component(vm, compnum)
         derivs%Ys(j) = (cpl-cmi)/dvar
      end do
      do j = 1, nsp
         dvar = v%Ncol(j)/DERIVDIV; vm = v
         vm%Ncol(j) = v%Ncol(j) + dvar/2.0d0; cpl = wae_get_component(vm, compnum)
         vm%Ncol(j) = v%Ncol(j) - dvar/2.0d0; cmi = wae_get_component(vm, compnum)
         derivs%Ncol(j) = (cpl-cmi)/dvar
      end do
      dvar = v%z/DERIVDIV; vm = v
      vm%z = v%z + dvar/2.0d0; cpl = wae_get_component(vm, compnum)
      vm%z = v%z - dvar/2.0d0; cmi = wae_get_component(vm, compnum)
      derivs%z = (cpl-cmi)/dvar
      end subroutine wae_get_component_derivs

      subroutine wae_linearize_dvdr_crit(x, y, m, ne)
      ! Setup once per iteration of dvdr_slope/last/q_last (L'Hopital at
      ! the critical point). Ported from soe.c:linearize_dvdr_crit.
      use wae_types, only: wae_varlist
      integer, intent(in) :: m, ne
      real*8,  intent(in) :: x(m), y(ne, m)
      integer :: k, j
      real*8  :: Mach, rad, gamma, mu, last_Mach
      real*8  :: dYsdr(nsp), dNdr(nsp), dmudr, dgdr
      real*8  :: d1_phi, d2_phi, norm_a, qinv, spQ, spQ1, spQ2
      real*8  :: a, b, c, sqrtarg, dvdr_crit
      type(wae_varlist)   :: dQ
      type(wae_i_eqnvars) :: v

      call wae_get_gamma(gamma)
      last_Mach = par%mach_limit/(1.0d0 + 2.0d0*(par%erfn*par%mach_limit  &
                  /par%rapidity))

      k = m - 1
      do
         v%q = x(k); v%v = y(1,k); v%z = y(2,k); v%rho = y(3,k); v%T = y(4,k)
         do j = 1, nsp
            v%Ys(j) = y(j+4,k); v%Ncol(j) = y(j+4+nsp,k)
         end do
         call wae_get_mu(mu, v)
         Mach = v%v/sqrt(gamma*v%T/mu)
         if (Mach .lt. last_Mach) then
            v%q = x(k+1); v%v = y(1,k+1); v%z = y(2,k+1)
            v%rho = y(3,k+1); v%T = y(4,k+1)
            do j = 1, nsp
               v%Ys(j) = y(j+4,k+1); v%Ncol(j) = y(j+4+nsp,k+1)
            end do
            exit
         end if
         k = k - 1
         if (k .le. 0) exit
      end do

      call wae_get_rad(rad, v)
      call wae_get_spQ(spQ, v, 1.0d4)
      norm_a = par%semimajor/par%Rp
      d1_phi = 1.0d0/(rad*rad)
      qinv = par%Mstar/par%Mp
      d1_phi = d1_phi + ( -qinv/(norm_a-rad)**2                          &
               + (qinv*(norm_a-rad) - rad)/norm_a**3 )*par%tidalforce
      d1_phi = d1_phi*(par%Rp/par%H0)
      wae_dvdr_last = 2.0d0*gamma*v%T/(mu*rad) - (gamma-1.0d0)*spQ/v%v - d1_phi
      wae_dvdr_last = wae_dvdr_last*v%v/(v%v*v%v - gamma*v%T/mu)
      wae_q_last = v%q

      ! critical point (k = M)
      v%q = x(m); v%v = y(1,m); v%z = y(2,m); v%rho = y(3,m); v%T = y(4,m)
      do j = 1, nsp
         v%Ys(j) = y(j+4,m); v%Ncol(j) = y(j+4+nsp,m)
      end do
      call wae_get_rad(rad, v)
      call wae_get_mu(mu, v)
      call wae_get_spQ(spQ, v, 1.0d4)
      call wae_get_component_derivs(dQ, v, QCOMP)
      call wae_get_dYsdr(dYsdr, v, 1.0d4)
      call wae_get_dNcoldr(dNdr, v)
      ! dmu/dr as upstream; the exact derivative is on hold (see wae_get_dTdr).
      ! To Be Checked: D18 on hold (user decision 2026-10-27), see md/windae_defect_fixes_20261027.md
      dgdr = 0.0d0
      do j = 1, nsp
         dgdr = dgdr - par%HX(j)*dYsdr(j)
      end do
      dmudr = -mu*mu*dgdr

      norm_a = par%semimajor/par%Rp
      d2_phi = -2.0d0/rad**3
      qinv = par%Mstar/par%Mp
      d2_phi = d2_phi + (-2.0d0*qinv/(norm_a-rad)**3                      &
               - (1.0d0+qinv)/norm_a**3)*par%tidalforce
      d2_phi = d2_phi*(par%Rp/par%H0)

      spQ1 = dQ%v + dQ%rho*(-v%rho/v%v) + dQ%T*(-(gamma-1.0d0)*v%T/v%v)
      spQ1 = spQ1/v%rho
      spQ2 = dQ%rho*(-2.0d0*v%rho/rad)                                    &
             + dQ%T*(v%T/mu*dmudr + (gamma-1.0d0)*spQ*mu/v%v              &
                     - 2.0d0*(gamma-1.0d0)*v%T/rad)
      do j = 1, nsp
         spQ2 = spQ2 + dQ%Ys(j)*min(dYsdr(j), 0.0d0) + dQ%Ncol(j)*dNdr(j)
         if (dYsdr(j) .gt. 0.0d0) &
            write(*,*) 'WARNING: dYsdr_', trim(par%species(j)), ' =', dYsdr(j), ' > 0'
      end do
      spQ2 = spQ2/v%rho

      a = 2.0d0*v%v + gamma*(gamma-1.0d0)*v%T/(mu*v%v)
      b = (gamma-1.0d0)*(4.0d0*gamma*v%T/(mu*rad) - gamma*spQ/v%v + spQ1)
      c = (gamma-1.0d0)**2*(v%v*d2_phi/(gamma-1.0d0)**2                   &
          + spQ2/(gamma-1.0d0)                                           &
          + 2.0d0*gamma*(2.0d0*gamma-1.0d0)*v%T*v%v/(mu*((gamma-1.0d0)*rad)**2) &
          - 2.0d0*spQ/rad)
      sqrtarg = b*b - 4.0d0*a*c
      if (sqrtarg .lt. 0.0d0) then
         wae_err = 3; return   ! graceful fail for continuation
      end if
      dvdr_crit = (-b + sqrt(sqrtarg))/(2.0d0*a)
      wae_dvdr_slope = (dvdr_crit - wae_dvdr_last)/(v%q - wae_q_last)
      end subroutine wae_linearize_dvdr_crit

      end module wae_soe
