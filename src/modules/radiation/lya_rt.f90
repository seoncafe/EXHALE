   module lya_rt
   ! Phase 3b: Ly-alpha radiative transfer by the escape-probability method,
   ! computed IN-LINE (every timestep) from the current state.
   !
   ! Because we use the Neufeld(1990)/Harrington analytic escape probability rather
   ! than a Monte Carlo, the field is cheap and is evaluated directly inside the
   ! run (like the Phase-3a H(n=2) feedback) instead of an offline pass. The result
   ! is the Voigt-profile-averaged Ly-alpha mean intensity J_lya(r) (Huang et al.
   ! 2023 Eq. 11), which excited_hydrogen uses to pump H(n=2) when jlya_mode = 2.
   ! The separate option of reading an externally-computed (e.g. real Monte Carlo)
   ! J_lya(r) profile is kept as jlya_mode = 1 (load_jlya_rt in excited_hydrogen).
   !
   ! Why escape probability, not Monte Carlo: the WASP-121b column is extremely
   ! thick at line center (tau0 ~ 5e7, a.tau0 ~ 1e4), so an acceleration-free MC is
   ! intractable and core-skipping is excluded. In the a.tau0 >> 1 regime the
   ! Neufeld/Harrington wing-escape solution is accurate.
   !
   ! METHOD (plane-parallel, substellar column):
   !   Dnu_D(r) = nu_lya sqrt(2 kT/mu)/c                       Doppler width
   !   a(r)     = (A_2p1s/4pi)/Dnu_D                           Voigt parameter
   !   tau(r)   = top-down trapezoid of nhi*C_lya/Dnu_D        line-center depth
   !   beta(r)  = pi^(-1/4) sqrt(a/tau)  (capped at 1)         wing escape prob.
   !              (photon escapes once the wing optical depth at x* ~ (a.tau)^1/2
   !               drops to unity; Neufeld 1990, Harrington 1973)
   !   P(r)     = alpha_B ne nhii + C_1s2p ne nhi              Ly-alpha production
   !              (recombination cascade, 1 per case-B recomb; electron-impact
   !               1s->2p -- the two internal sources used by Huang)
   !   J_int(r) = (2 h nu^3/c^2)(g1s/g2p) P / (A_2p1s beta nhi)
   !              (trapped source-function buildup: n2p = P/(A beta), Eq. 11)
   !   J_star(r)= xi F_Lya_star/(4 pi Dnu_D) exp(-tau)         stellar beam (dayside
   !              dilution xi as for ground-state EUV; attenuated to the photosphere)
   !   J_lya = J_int + J_star.

   use global_parameters

   implicit none

   ! ----- Ly-alpha atomic data (cgs; same values as excited_hydrogen.f90) ----- !
   real*8, parameter :: lA_lya  = 1215.6701d-8       ! Ly-alpha wavelength [cm]
   real*8, parameter :: nu_lya  = c_light/lA_lya     ! Ly-alpha frequency [s^-1]
   real*8, parameter :: A_2p1s  = 6.3d8              ! A(2p->1s) [s^-1]
   real*8, parameter :: g1s_l   = 2.0d0, g2p_l = 6.0d0
   real*8, parameter :: f_lya   = 0.4162d0           ! Ly-alpha oscillator strength
   real*8, parameter :: C_lya   = 1.49736d-2*f_lya   ! sqrt(pi) e^2/(m_e c) f [cm^2 Hz]
   real*8, parameter :: beta_c  = 0.7511255d0        ! pi^(-1/4), wing-escape prefactor

   ! Diagnostic split of J_lya into internal (recomb+collisional) and stellar
   ! contributions, filled each call for output/Excited_H.txt.
   real*8, dimension(1-Ng:N+Ng) :: jint_arr = 0.0d0, jstar_arr = 0.0d0

   contains

   ! --------------------------------------------------------------- !

   subroutine jlya_escape_prob(T_K, nhi, nhii, ne, v_in, Jlya, tau)
   ! Voigt-averaged Ly-alpha mean intensity J_lya(r) [erg s^-1 cm^-2 Hz^-1 sr^-1]
   ! by the escape-probability method, and the line-center optical depth tau(r).
   ! T_K/nhi/nhii/ne are physical (cgs) per-cell arrays; v_in is the dimensionless
   ! radial velocity (v_in*v0 = cm/s). Called in-line from excited_H_update.
   real*8, dimension(1-Ng:N+Ng), intent(in)  :: T_K, nhi, nhii, ne, v_in
   real*8, dimension(1-Ng:N+Ng), intent(out) :: Jlya, tau

   integer :: j, jm, jp
   real*8, dimension(1-Ng:N+Ng) :: DnuD
   real*8 :: avoigt, beta_esc, Tl, t4, aB, C1s2p, Prec, Pcol, Jint, Jstar, Jpref, xi
   real*8 :: vth, Xs, x1, Tstar, Dnu_star
   real*8 :: C_sob, dvdr, tau_sob, beta_sob, beta_tot

   ! Doppler width (needed for both tau and the source function).
   do j = 1-Ng, N+Ng
      DnuD(j) = nu_lya*sqrt(2.0d0*kb_erg*max(T_K(j),1.0d0)/mu)/c_light
   enddo

   ! Top-down line-center optical depth to the outer surface.
   tau(N+Ng) = 0.0d0
   do j = N+Ng-1, 1-Ng, -1
      tau(j) = tau(j+1)                                                     &
             + 0.5d0*( nhi(j)  *C_lya/max(DnuD(j),  1.0d-30)                &
                     + nhi(j+1)*C_lya/max(DnuD(j+1),1.0d-30) )             &
             *(r(j+1) - r(j))*R0
   enddo

   ! Dayside dilution for the stellar beam (same factor as ground-state EUV).
   if      (index(appx_mth,'Rate/4') .gt. 0) then
      xi = 0.25d0
   else if (index(appx_mth,'Rate/2') .gt. 0) then
      xi = 0.5d0
   else
      xi = 1.0d0
   endif

   Jpref = 2.0d0*hp_erg*nu_lya**3.0d0/c_light**2.0d0     ! 2 h nu^3 / c^2
   ! Stellar Ly-alpha line width [Hz]: the BROAD stellar profile sets the
   ! profile-averaged mean intensity (Jbar = J_nu(nu0) for a stellar field flat
   ! across the narrow planetary line), NOT the planetary Doppler width Dnu_D.
   ! Dividing by Dnu_D would make J_star rise spuriously as the wind cools.
   Dnu_star = nu_lya*(dv_star_lya*1.0d5)/c_light
   ! Sobolev (velocity-gradient) escape coefficient: tau_S = C_sob n1s/|dv/dr|,
   ! C_sob = (pi e^2/m_e c) f_lya lambda = sqrt(pi) C_lya lambda.
   C_sob = sqrt(pi)*C_lya*lA_lya
   do j = 1-Ng, N+Ng
      avoigt   = (A_2p1s/(4.0d0*pi))/max(DnuD(j),1.0d-30)
      beta_esc = min(1.0d0, beta_c*sqrt(avoigt/max(tau(j),1.0d-30)))

      ! Sobolev escape via the wind velocity gradient -- a parallel channel to the
      ! static wing escape: the photon Doppler-shifts out of resonance over the
      ! Sobolev length. |dv/dr| in cgs = (v0/R0)|d v_in/d r|.
      jm = max(j-1, 1-Ng); jp = min(j+1, N+Ng)
      dvdr     = (v0/R0)*abs(v_in(jp) - v_in(jm))/max(abs(r(jp) - r(jm)),1.0d-30)
      tau_sob  = C_sob*max(nhi(j),0.0d0)/max(dvdr,1.0d-30)
      beta_sob = min(1.0d0,(1.0d0 - exp(-min(tau_sob,200.0d0)))/max(tau_sob,1.0d-30))
      ! Total escape probability: escape via the static OR the Sobolev channel.
      beta_tot = 1.0d0 - (1.0d0 - beta_esc)*(1.0d0 - beta_sob)

      Tl    = max(T_K(j), 1.0d0)
      t4    = Tl/1.0d4
      aB    = 2.54d-13*t4**(-0.8163d0 - 0.0208d0*log(t4))      ! case-B [cm^3/s]
      C1s2p = 1.71d-8*(1.0d0/t4)**0.077d0*exp(-118400.0d0/Tl)  ! 1s->2p [cm^3/s]
      Prec  = aB*ne(j)*nhii(j)
      Pcol  = C1s2p*ne(j)*nhi(j)

      ! Trapped internal field, escape-probability closure Jbar = S(1-beta):
      ! S = (2hv3/c2)(g1s/g2p) n2p/n1s with the trapped pile-up n2p = P/(A beta).
      ! The (1-beta) factor makes Jbar -> S in the thick limit and -> 0 when thin
      ! (beta -> 1), instead of the source function S diverging as n1s -> 0 in the
      ! ionized outer wind (which spuriously raised Jbar outward).
      Jint  = Jpref*(g1s_l/g2p_l)*(Prec + Pcol)*(1.0d0 - beta_tot)         &
            /(A_2p1s*max(beta_tot,1.0d-30)*max(nhi(j),1.0d-30))

      ! Stellar beam: resonantly SCATTERED, not destroyed -- so it is not killed by
      ! exp(-tau); instead the broad stellar line penetrates through its wings and
      ! FILLS the region it reaches with the free-streaming mean intensity
      ! F/(4 pi Dnu_D). A stellar photon at Doppler offset x reaches depth r if
      ! |x| > x1 = sqrt(a tau/sqrt(pi)) (planetary wing optical depth 1); the
      ! fraction of a Gaussian stellar profile of Doppler width Xs = dv_star/v_th
      ! beyond x1 is T_star = erfc(x1/(sqrt(2) Xs)) -- smooth, so the penetration
      ! photosphere is not a sharp edge (a flat profile gives a kink). (No 1/beta
      ! buildup: that applies to a VOLUMETRIC source, not an incident beam.)
      ! The stellar term still dominates the internal sources where it penetrates.
      vth   = sqrt(2.0d0*kb_erg*Tl/mu)                     ! thermal speed [cm/s]
      Xs    = (dv_star_lya*1.0d5)/max(vth,1.0d-30)         ! stellar width [Doppler]
      x1    = sqrt(avoigt*tau(j)/sqrt(pi))                 ! wing-escape frequency
      Tstar = erfc(x1/(sqrt(2.0d0)*max(Xs,1.0d-30)))
      ! Bounded trapping buildup of the penetrating stellar beam: the scattered
      ! photons raise the mean intensity above free-streaming by
      !   E = 1 + (lya_star_boost - 1)*(1 - beta)*T_star,
      ! which is ~lya_star_boost where the beam fully penetrates and is trapped
      ! (T_star~1, beta<<1, the wind) and -> 1 both in the thin outer wind
      ! (beta -> 1) and at the deep penetration edge (T_star -> 0), so it does not
      ! over-count the base. (The full 1/beta blows up; a flat cap over-enhances
      ! the deep edge.) lya_star_boost tuned to Huang+2023 Fig. 11.
      Jstar = xi*F_Lya_star/(4.0d0*pi*Dnu_star)*Tstar                        &
            * (1.0d0 + (lya_star_boost - 1.0d0)*(1.0d0 - beta_tot)*Tstar)

      jint_arr(j)  = Jint
      jstar_arr(j) = Jstar
      Jlya(j) = Jint + Jstar
   enddo

   end subroutine jlya_escape_prob

   ! --------------------------------------------------------------- !

   end module lya_rt
